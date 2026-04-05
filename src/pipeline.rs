use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use tokio::fs;
use tokio::process::Command;
use tracing::info;

use crate::config::Config;
use crate::error::RustycleanError;
use crate::sample::Sample;

// ============================================================================
// Pipeline stages (ordered for comparison)
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum PipelineStage {
    Pending,
    FastpRunning,
    FastpComplete,
    Kraken2Running,
    Kraken2Complete,
    Validating,
    Completed,
    Failed,
}

impl std::fmt::Display for PipelineStage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

// ============================================================================
// Metrics
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FastpMetrics {
    pub input_reads: u64,
    pub output_reads: u64,
    pub q20_rate: f64,
    pub q30_rate: f64,
    pub gc_content: f64,
    pub adapter_trimmed: u64,
    pub too_short_reads: u64,
    pub low_quality_reads: u64,
    pub output_r1: PathBuf,
    pub output_r2: Option<PathBuf>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Kraken2Metrics {
    pub classified_reads: u64,
    pub unclassified_reads: u64,
    pub human_reads: u64,
    pub contamination_percent: f64,
    pub output_r1: PathBuf,
    pub output_r2: Option<PathBuf>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationResult {
    pub passed: bool,
    pub file_size_bytes: u64,
    pub validation_timestamp: DateTime<Utc>,
    pub errors: Vec<String>,
}

// ============================================================================
// Checkpoint
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Checkpoint {
    pub version: u32,
    pub sample_id: String,
    pub stage: PipelineStage,
    pub input_hash: u64,
    pub created_at: DateTime<Utc>,
    pub updated_at: DateTime<Utc>,
    pub attempts: u32,
    pub fastp_metrics: Option<FastpMetrics>,
    pub kraken2_metrics: Option<Kraken2Metrics>,
    pub validation_result: Option<ValidationResult>,
    pub last_error: Option<String>,
}

impl Checkpoint {
    pub const CURRENT_VERSION: u32 = 1;

    pub fn new(sample_id: String, input_hash: u64) -> Self {
        let now = Utc::now();
        Self {
            version: Self::CURRENT_VERSION,
            sample_id,
            stage: PipelineStage::Pending,
            input_hash,
            created_at: now,
            updated_at: now,
            attempts: 0,
            fastp_metrics: None,
            kraken2_metrics: None,
            validation_result: None,
            last_error: None,
        }
    }

    pub fn transition(&mut self, new_stage: PipelineStage) {
        self.stage = new_stage;
        self.updated_at = Utc::now();
    }

    pub fn increment_attempt(&mut self) {
        self.attempts += 1;
        self.updated_at = Utc::now();
    }

    pub fn record_error(&mut self, error: String) {
        self.last_error = Some(error);
        self.transition(PipelineStage::Failed);
    }

    pub fn can_resume(&self) -> bool {
        matches!(
            self.stage,
            PipelineStage::Pending
                | PipelineStage::Failed
                | PipelineStage::FastpRunning
                | PipelineStage::Kraken2Running
        )
    }

    pub fn is_complete(&self) -> bool {
        self.stage == PipelineStage::Completed
    }
}

// ============================================================================
// Pipeline execution
// ============================================================================

/// Run the full pipeline for a single sample.
pub async fn execute_pipeline(
    sample: &Sample,
    checkpoint: &mut Checkpoint,
    config: &Config,
    work_dir: &Path,
) -> Result<()> {
    // Stage 1: fastp
    if checkpoint.stage < PipelineStage::FastpComplete {
        checkpoint.transition(PipelineStage::FastpRunning);
        run_fastp(sample, checkpoint, config, work_dir).await?;
        checkpoint.transition(PipelineStage::FastpComplete);
    }

    // Stage 2: kraken2
    if checkpoint.stage < PipelineStage::Kraken2Complete {
        checkpoint.transition(PipelineStage::Kraken2Running);
        run_kraken2(sample, checkpoint, config, work_dir).await?;
        checkpoint.transition(PipelineStage::Kraken2Complete);
    }

    // Stage 3: validation & move output
    checkpoint.transition(PipelineStage::Validating);
    validate_and_finalize(sample, checkpoint, config).await?;
    checkpoint.transition(PipelineStage::Completed);

    Ok(())
}

// ============================================================================
// fastp
// ============================================================================

async fn run_fastp(
    sample: &Sample,
    checkpoint: &mut Checkpoint,
    config: &Config,
    work_dir: &Path,
) -> Result<()> {
    let fastp_cfg = &config.tools.fastp;
    let trimmed_r1 = work_dir.join("trimmed_R1.fastq.gz");
    let trimmed_r2 = if sample.is_paired() {
        Some(work_dir.join("trimmed_R2.fastq.gz"))
    } else {
        None
    };
    let json_report = work_dir.join("fastp.json");

    let mut cmd = Command::new("fastp");
    cmd.arg("--in1").arg(&sample.r1)
        .arg("--out1").arg(&trimmed_r1)
        .arg("--json").arg(&json_report)
        .arg("--thread").arg(fastp_cfg.threads.to_string())
        .arg("--compression").arg(fastp_cfg.compression_level.to_string());

    if let (Some(r2), Some(out2)) = (&sample.r2, &trimmed_r2) {
        cmd.arg("--in2").arg(r2).arg("--out2").arg(out2);
        if fastp_cfg.detect_adapters {
            cmd.arg("--detect_adapter_for_pe");
        }
    }

    if fastp_cfg.cut_front {
        cmd.arg("--cut_front");
    }
    if fastp_cfg.cut_tail {
        cmd.arg("--cut_tail");
    }

    cmd.arg("--qualified_quality_phred")
        .arg(fastp_cfg.qualified_quality_phred.to_string())
        .arg("--length_required")
        .arg(fastp_cfg.length_required.to_string());

    let output = cmd
        .output()
        .await
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to run fastp: {}", e)))?;

    if !output.status.success() {
        return Err(
            RustycleanError::ToolExecution(String::from_utf8_lossy(&output.stderr).to_string()).into(),
        );
    }

    // Parse fastp JSON
    let metrics = parse_fastp_json(&json_report, &trimmed_r1, trimmed_r2.as_deref()).await?;
    checkpoint.fastp_metrics = Some(metrics);

    Ok(())
}

async fn parse_fastp_json(
    json_path: &Path,
    out_r1: &Path,
    out_r2: Option<&Path>,
) -> Result<FastpMetrics> {
    let content = fs::read_to_string(json_path).await?;
    let json: serde_json::Value = serde_json::from_str(&content)?;

    let summary = json
        .get("summary")
        .context("fastp JSON missing 'summary'")?;
    let before = summary
        .get("before_filtering")
        .context("fastp JSON missing 'before_filtering'")?;
    let after = summary
        .get("after_filtering")
        .context("fastp JSON missing 'after_filtering'")?;

    Ok(FastpMetrics {
        input_reads: before
            .get("total_reads")
            .and_then(|v| v.as_u64())
            .unwrap_or(0),
        output_reads: after
            .get("total_reads")
            .and_then(|v| v.as_u64())
            .unwrap_or(0),
        q20_rate: after
            .get("q20_rate")
            .and_then(|v| v.as_f64())
            .unwrap_or(0.0),
        q30_rate: after
            .get("q30_rate")
            .and_then(|v| v.as_f64())
            .unwrap_or(0.0),
        gc_content: after
            .get("gc_content")
            .and_then(|v| v.as_f64())
            .unwrap_or(0.0),
        adapter_trimmed: json
            .get("adapter_cutting")
            .and_then(|v| v.get("adapter_trimmed_reads"))
            .and_then(|v| v.as_u64())
            .unwrap_or(0),
        too_short_reads: json
            .get("filtering_result")
            .and_then(|v| v.get("too_short_reads"))
            .and_then(|v| v.as_u64())
            .unwrap_or(0),
        low_quality_reads: json
            .get("filtering_result")
            .and_then(|v| v.get("low_quality_reads"))
            .and_then(|v| v.as_u64())
            .unwrap_or(0),
        output_r1: out_r1.to_path_buf(),
        output_r2: out_r2.map(|p| p.to_path_buf()),
    })
}

// ============================================================================
// kraken2
// ============================================================================

async fn run_kraken2(
    sample: &Sample,
    checkpoint: &mut Checkpoint,
    config: &Config,
    work_dir: &Path,
) -> Result<()> {
    let kraken_cfg = &config.tools.kraken2;
    let fastp = checkpoint
        .fastp_metrics
        .as_ref()
        .context("fastp metrics missing before kraken2 stage")?;

    let kraken_report = work_dir.join("kraken2.report");

    // kraken2 --unclassified-out with paired needs a special #-based naming
    let mut cmd = Command::new("kraken2");
    cmd.arg("--db")
        .arg(&kraken_cfg.db_path)
        .arg("--threads")
        .arg(kraken_cfg.threads.to_string())
        .arg("--confidence")
        .arg(kraken_cfg.confidence_threshold.to_string())
        .arg("--minimum-hit-groups")
        .arg(kraken_cfg.minimum_hit_groups.to_string())
        .arg("--report")
        .arg(&kraken_report)
        .arg("--gzip-compressed");

    let (output_r1, output_r2);

    if sample.is_paired() {
        // For paired reads, kraken2 uses # in --unclassified-out to split R1/R2
        let unclass_base = work_dir.join("clean#.fastq.gz");
        output_r1 = work_dir.join("clean_1.fastq.gz");
        output_r2 = Some(work_dir.join("clean_2.fastq.gz"));

        cmd.arg("--paired")
            .arg("--unclassified-out")
            .arg(&unclass_base)
            .arg(&fastp.output_r1)
            .arg(fastp.output_r2.as_ref().unwrap());
    } else {
        output_r1 = work_dir.join("clean.fastq.gz");
        output_r2 = None;

        cmd.arg("--unclassified-out")
            .arg(&output_r1)
            .arg(&fastp.output_r1);
    }

    let output = cmd
        .output()
        .await
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to run kraken2: {}", e)))?;

    let stderr = String::from_utf8_lossy(&output.stderr);
    if stderr.contains("Error:") || stderr.contains("FATAL") {
        return Err(RustycleanError::ToolExecution(stderr.to_string()).into());
    }

    // Parse report
    let metrics = parse_kraken_report(&kraken_report, output_r1, output_r2).await?;
    checkpoint.kraken2_metrics = Some(metrics);

    Ok(())
}

async fn parse_kraken_report(
    report_path: &Path,
    output_r1: PathBuf,
    output_r2: Option<PathBuf>,
) -> Result<Kraken2Metrics> {
    let content = fs::read_to_string(report_path).await?;
    let mut classified = 0u64;
    let mut unclassified = 0u64;
    let mut human = 0u64;

    for line in content.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }

        let count: u64 = fields[1].trim().parse().unwrap_or(0);
        let rank = fields[3].trim();
        let name = fields[5].trim();

        if name == "unclassified" {
            unclassified = count;
        } else if name == "root" {
            classified = count;
        } else if rank == "S" && name == "Homo sapiens" {
            human = count;
        }
    }

    let total = classified + unclassified;
    let contamination = if total > 0 {
        (human as f64 / total as f64) * 100.0
    } else {
        0.0
    };

    Ok(Kraken2Metrics {
        classified_reads: classified,
        unclassified_reads: unclassified,
        human_reads: human,
        contamination_percent: contamination,
        output_r1,
        output_r2,
    })
}

// ============================================================================
// Validation & finalize
// ============================================================================

async fn validate_and_finalize(
    sample: &Sample,
    checkpoint: &mut Checkpoint,
    config: &Config,
) -> Result<()> {
    let kraken = checkpoint
        .kraken2_metrics
        .as_ref()
        .context("kraken2 metrics missing at validation stage")?;

    let mut errors = Vec::new();

    // Check output file size
    let metadata = fs::metadata(&kraken.output_r1)
        .await
        .map_err(|e| RustycleanError::ValidationFailed(sample.id.clone(), format!("R1 output missing: {}", e)))?;

    if metadata.len() < config.validation.min_output_size_bytes {
        errors.push(format!(
            "R1 output too small: {} bytes",
            metadata.len()
        ));
    }

    if let Some(r2) = &kraken.output_r2 {
        let m2 = fs::metadata(r2)
            .await
            .map_err(|e| RustycleanError::ValidationFailed(sample.id.clone(), format!("R2 output missing: {}", e)))?;
        if m2.len() < config.validation.min_output_size_bytes {
            errors.push(format!("R2 output too small: {} bytes", m2.len()));
        }
    }

    // Check contamination
    if kraken.contamination_percent > config.validation.max_contamination_percent {
        errors.push(format!(
            "Contamination too high: {:.2}% (max {}%)",
            kraken.contamination_percent, config.validation.max_contamination_percent
        ));
    }

    let passed = errors.is_empty();

    checkpoint.validation_result = Some(ValidationResult {
        passed,
        file_size_bytes: metadata.len(),
        validation_timestamp: Utc::now(),
        errors: errors.clone(),
    });

    if !passed {
        return Err(
            RustycleanError::ValidationFailed(sample.id.clone(), errors.join("; ")).into(),
        );
    }

    // Move final outputs to destination
    fs::create_dir_all(&sample.output_dir).await?;

    let final_r1 = sample.output_dir.join(format!("{}_clean_R1.fastq.gz", sample.id));
    fs::rename(&kraken.output_r1, &final_r1).await?;

    if let Some(r2_src) = &kraken.output_r2 {
        let final_r2 = sample.output_dir.join(format!("{}_clean_R2.fastq.gz", sample.id));
        fs::rename(r2_src, &final_r2).await?;
    }

    info!(
        sample = %sample.id,
        unclassified_reads = kraken.unclassified_reads,
        human_reads = kraken.human_reads,
        contamination = format!("{:.2}%", kraken.contamination_percent),
        "Sample validated and finalized"
    );

    Ok(())
}
