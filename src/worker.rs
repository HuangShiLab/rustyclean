use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::Duration;

use anyhow::Result;
use console::style;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use tokio::sync::Semaphore;
use tracing::{error, info, warn};

use crate::checkpoint::CheckpointManager;
use crate::config::Config;
use crate::pipeline::{self, Checkpoint, PipelineStage};
use crate::sample::Sample;

// ============================================================================
// Progress tracking
// ============================================================================

pub struct ProgressTracker {
    pub total: AtomicUsize,
    pub completed: AtomicUsize,
    pub failed: AtomicUsize,
    pub skipped: AtomicUsize,
    pub multi_progress: MultiProgress,
}

impl ProgressTracker {
    pub fn new(total: usize) -> Self {
        Self {
            total: AtomicUsize::new(total),
            completed: AtomicUsize::new(0),
            failed: AtomicUsize::new(0),
            skipped: AtomicUsize::new(0),
            multi_progress: MultiProgress::new(),
        }
    }

    pub fn create_bar(&self, sample_id: &str) -> ProgressBar {
        let sty = ProgressStyle::default_spinner()
            .template("{spinner:.green} {prefix:20} {msg} {elapsed_precise}")
            .unwrap();

        let bar = self.multi_progress.add(ProgressBar::new_spinner());
        bar.set_style(sty);
        bar.set_prefix(sample_id.to_string());
        bar.enable_steady_tick(Duration::from_millis(100));
        bar
    }

    pub fn report_final(&self) {
        let total = self.total.load(Ordering::Relaxed);
        let completed = self.completed.load(Ordering::Relaxed);
        let failed = self.failed.load(Ordering::Relaxed);
        let skipped = self.skipped.load(Ordering::Relaxed);

        println!("\n{}", style("=".repeat(60)).cyan());
        println!("{}", style("Pipeline Complete").bold().cyan());
        println!("  Total:     {}", style(total).bold());
        println!("  Success:   {}", style(completed).green().bold());
        println!("  Skipped:   {}", style(skipped).yellow().bold());
        println!("  Failed:    {}", style(failed).red().bold());
        println!("{}\n", style("=".repeat(60)).cyan());
    }
}

// ============================================================================
// Pipeline summary
// ============================================================================

#[derive(Default, Debug)]
pub struct PipelineSummary {
    pub successful: Vec<String>,
    pub failed: Vec<(String, Option<String>)>,
    pub skipped: Vec<String>,
}

// ============================================================================
// Worker pool
// ============================================================================

pub async fn run_pipeline(
    samples: Vec<Sample>,
    config: Config,
    checkpoint_mgr: Arc<CheckpointManager>,
    work_base: PathBuf,
) -> Result<PipelineSummary> {
    let total = samples.len();
    let progress = Arc::new(ProgressTracker::new(total));
    let semaphore = Arc::new(Semaphore::new(config.execution.workers));
    let config = Arc::new(config);

    let cancel = tokio_util::sync::CancellationToken::new();

    // Install Ctrl-C handler
    let ctrl_c_handle = {
        let cancel = cancel.clone();
        tokio::spawn(async move {
            let _ = tokio::signal::ctrl_c().await;
            warn!("Ctrl-C received, cancelling remaining work...");
            cancel.cancel();
        })
    };

    let mut handles = Vec::with_capacity(total);

    for sample in samples {
        let permit = semaphore.clone().acquire_owned().await?;
        let config = config.clone();
        let checkpoint_mgr = checkpoint_mgr.clone();
        let progress = progress.clone();
        let cancel = cancel.clone();
        let work_base = work_base.clone();

        let handle = tokio::spawn(async move {
            let _permit = permit; // held until task completes

            if cancel.is_cancelled() {
                return;
            }

            // Initialize checkpoint
            let input_hash = match sample.compute_hash().await {
                Ok(h) => h,
                Err(e) => {
                    error!(sample = %sample.id, "Failed to compute input hash: {}", e);
                    progress.failed.fetch_add(1, Ordering::Relaxed);
                    return;
                }
            };

            let mut checkpoint = match checkpoint_mgr.load(&sample.id).await {
                Ok(Some(cp)) if cp.input_hash == input_hash && cp.is_complete() => {
                    info!(sample = %sample.id, "Already completed, skipping");
                    progress.skipped.fetch_add(1, Ordering::Relaxed);
                    return;
                }
                Ok(Some(cp)) if cp.input_hash == input_hash && cp.can_resume() => {
                    info!(sample = %sample.id, stage = ?cp.stage, "Resuming from checkpoint");
                    cp
                }
                _ => Checkpoint::new(sample.id.clone(), input_hash),
            };

            let bar = progress.create_bar(&sample.id);
            let work_dir = work_base.join(&sample.id);
            if let Err(e) = tokio::fs::create_dir_all(&work_dir).await {
                error!(sample = %sample.id, "Failed to create work dir: {}", e);
                progress.failed.fetch_add(1, Ordering::Relaxed);
                bar.finish_with_message(style("x Failed").red().to_string());
                return;
            }

            let max_attempts = config.execution.retry_attempts + 1;

            loop {
                checkpoint.increment_attempt();
                let _ = checkpoint_mgr.save(&checkpoint).await;

                bar.set_message(format!("attempt {}/{}", checkpoint.attempts, max_attempts));

                match pipeline::execute_pipeline(&sample, &mut checkpoint, &config, &work_dir)
                    .await
                {
                    Ok(()) => {
                        let _ = checkpoint_mgr.save(&checkpoint).await;
                        bar.finish_with_message(style("OK").green().to_string());
                        progress.completed.fetch_add(1, Ordering::Relaxed);
                        return;
                    }
                    Err(e) if checkpoint.attempts < max_attempts && !cancel.is_cancelled() => {
                        warn!(
                            sample = %sample.id,
                            attempt = checkpoint.attempts,
                            "Failed: {}, retrying...", e
                        );
                        // Reset stage to allow retry from the failed point
                        tokio::time::sleep(Duration::from_secs(3)).await;
                    }
                    Err(e) => {
                        checkpoint.record_error(e.to_string());
                        let _ = checkpoint_mgr.save(&checkpoint).await;
                        error!(sample = %sample.id, "Failed permanently: {}", e);
                        bar.finish_with_message(style("x Failed").red().to_string());
                        progress.failed.fetch_add(1, Ordering::Relaxed);
                        return;
                    }
                }
            }
        });

        handles.push(handle);
    }

    // Wait for all tasks
    for handle in handles {
        let _ = handle.await;
    }

    // The Ctrl-C listener task would otherwise keep the tokio runtime alive
    // after all pipeline work has finished, causing the process to hang.
    ctrl_c_handle.abort();
    let _ = ctrl_c_handle.await;

    progress.report_final();

    // Build summary
    let mut summary = PipelineSummary::default();
    let all_checkpoints = checkpoint_mgr.load_all().await.unwrap_or_default();
    for (id, cp) in &all_checkpoints {
        match cp.stage {
            PipelineStage::Completed => summary.successful.push(id.clone()),
            PipelineStage::Failed => summary.failed.push((id.clone(), cp.last_error.clone())),
            _ => {}
        }
    }

    Ok(summary)
}
