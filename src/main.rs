mod checkpoint;
mod cli;
mod config;
mod error;
mod pipeline;
mod sample;
mod worker;

use std::sync::Arc;

use anyhow::{bail, Result};
use clap::Parser;
use tracing::info;
use tracing_subscriber::EnvFilter;

use crate::checkpoint::CheckpointManager;
use crate::cli::{Cli, HostRemovalModeCli};
use crate::config::{Config, ExecutionConfig, FastpConfig, HostRemovalConfig, ToolConfig, ValidationConfig};
use crate::sample::{parse_sample_list, sample_from_paths, validate_inputs};

#[tokio::main]
async fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info")),
        )
        .init();

    let cli = Cli::parse();

    // Load or build config
    let mut config = if let Some(config_path) = &cli.config {
        let content = tokio::fs::read_to_string(config_path).await?;
        toml::from_str::<Config>(&content)?
    } else {
        Config {
            execution: ExecutionConfig::default(),
            tools: ToolConfig {
                fastp: FastpConfig::default(),
                host_removal: HostRemovalConfig::default(),
            },
            validation: ValidationConfig::default(),
        }
    };

    // Apply CLI overrides
    if let Some(w) = cli.workers {
        config.execution.workers = w;
    }
    if let Some(t) = cli.threads {
        config.tools.fastp.threads = t;
        config.tools.host_removal = set_host_removal_threads(config.tools.host_removal, t);
    }
    config.tools.fastp.skip_qc = cli.skip_qc;

    // Set host-removal backend and paths
    let host_index = cli.host_index.as_deref();
    config.tools.host_removal = match cli.host_removal_mode {
        HostRemovalModeCli::Kraken2 => {
            let db_path = cli.kraken2_db.unwrap_or_else(|| {
                match &config.tools.host_removal {
                    HostRemovalConfig::Kraken2 { db_path, .. } => db_path.clone(),
                    _ => std::path::PathBuf::from("/db/minikraken2_v2_8GB"),
                }
            });
            HostRemovalConfig::Kraken2 {
                db_path,
                threads: config.tools.host_removal.threads(),
                confidence_threshold: 0.0,
                minimum_hit_groups: 2,
                memory_mapping: cli.kraken2_memory_mapping,
                bowtie2_recheck: cli.bowtie2_recheck,
                bowtie2_index_prefix: cli.host_index.clone(),
            }
        }
        HostRemovalModeCli::Minimap2 => {
            let index_path = host_index
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| std::path::PathBuf::from("/db/human_t2t_hla.mmi"));
            HostRemovalConfig::Minimap2 {
                index_path,
                threads: config.tools.host_removal.threads(),
            }
        }
        HostRemovalModeCli::Bowtie2 => {
            let index_prefix = host_index
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| std::path::PathBuf::from("/db/human_t2t_hla"));
            HostRemovalConfig::Bowtie2 {
                index_prefix,
                threads: config.tools.host_removal.threads(),
            }
        }
        HostRemovalModeCli::Sylph => {
            let db_path = host_index
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| std::path::PathBuf::from("/db/human_t2t_hla.sylph.syldb"));
            HostRemovalConfig::Sylph {
                db_path,
                threads: config.tools.host_removal.threads(),
            }
        }
        HostRemovalModeCli::Centrifuge => {
            let db_path = host_index
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| std::path::PathBuf::from("/db/human_t2t_hla_cf"));
            HostRemovalConfig::Centrifuge {
                db_path,
                threads: config.tools.host_removal.threads(),
            }
        }
        HostRemovalModeCli::Auto => {
            let kraken2_db_path = cli.kraken2_db.unwrap_or_else(|| {
                match &config.tools.host_removal {
                    HostRemovalConfig::Kraken2 { db_path, .. } => db_path.clone(),
                    HostRemovalConfig::Auto { kraken2_db_path, .. } => kraken2_db_path.clone(),
                    _ => std::path::PathBuf::from("/db/minikraken2_v2_8GB"),
                }
            });
            let bowtie2_index_prefix = host_index
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| {
                    match &config.tools.host_removal {
                        HostRemovalConfig::Bowtie2 { index_prefix, .. } => index_prefix.clone(),
                        HostRemovalConfig::Auto { bowtie2_index_prefix, .. } => bowtie2_index_prefix.clone(),
                        _ => std::path::PathBuf::from("/db/human_t2t_hla"),
                    }
                });

            // Validate that both required databases are available for auto mode.
            if !kraken2_db_path.exists() {
                bail!("--kraken2-db path does not exist: {}", kraken2_db_path.display());
            }
            // bowtie2 index prefix: check at least the .1.bt2 file exists
            let bt2_test = bowtie2_index_prefix.with_extension("1.bt2");
            if !bt2_test.exists() {
                bail!("bowtie2 index not found at prefix: {}", bowtie2_index_prefix.display());
            }

            HostRemovalConfig::Auto {
                kraken2_db_path,
                bowtie2_index_prefix,
                threads: config.tools.host_removal.threads(),
                host_pct_low_threshold: cli.auto_low_threshold,
                host_pct_high_threshold: cli.auto_high_threshold,
                reads_high_threshold: cli.auto_reads_threshold,
                user_host_pct: cli.host_pct,
                survey: cli.auto_survey,
                survey_n_reads: cli.auto_survey_nreads,
                survey_threads: cli.auto_survey_threads,
                bowtie2_recheck: cli.bowtie2_recheck,
            }
        }
    };

    if cli.resume {
        config.execution.resume = true;
    }

    config.validation.max_contamination_percent = cli.max_contamination;

    // Build sample list
    let samples = if let Some(list_path) = &cli.samples {
        parse_sample_list(list_path, &cli.output).await?
    } else {
        let r1 = cli.r1.clone().unwrap();
        let r2 = cli.r2.clone();
        vec![sample_from_paths(r1, r2, &cli.output)?]
    };

    if samples.is_empty() {
        bail!("No samples to process");
    }

    info!(
        "Loaded {} sample(s), host-removal mode: {}",
        samples.len(),
        config.tools.host_removal.mode()
    );

    // Validate inputs
    validate_inputs(&samples).await?;

    // Dry run
    if cli.dry_run {
        println!("Dry run - {} sample(s):", samples.len());
        for s in &samples {
            let mode = if s.is_paired() { "PE" } else { "SE" };
            println!(
                "  [{}] {}: {} {}",
                mode,
                s.id,
                s.r1.display(),
                s.r2.as_ref().map_or(String::new(), |p| p.display().to_string()),
            );
        }
        return Ok(());
    }

    // Check external tools
    check_tools(&config.tools.host_removal, config.tools.fastp.skip_qc)?;

    // Checkpoint manager
    let checkpoint_mgr = Arc::new(CheckpointManager::new(cli.checkpoint_dir.clone()).await?);

    // Work directory for intermediate files
    let work_base = cli.checkpoint_dir.join("work");
    tokio::fs::create_dir_all(&work_base).await?;

    // Run pipeline
    let summary = worker::run_pipeline(samples, config, checkpoint_mgr.clone(), work_base).await?;

    // Clean checkpoints if requested
    if cli.clean {
        let cleaned = checkpoint_mgr.clean_completed(&summary.successful).await?;
        info!("Cleaned {} completed checkpoint(s)", cleaned);
    }

    if !summary.failed.is_empty() {
        eprintln!("\nFailed samples:");
        for (id, err) in &summary.failed {
            eprintln!("  {}: {}", id, err.as_deref().unwrap_or("unknown error"));
        }
        bail!("{} sample(s) failed", summary.failed.len());
    }

    Ok(())
}

fn set_host_removal_threads(cfg: HostRemovalConfig, threads: usize) -> HostRemovalConfig {
    match cfg {
        HostRemovalConfig::Kraken2 { db_path, confidence_threshold, minimum_hit_groups, memory_mapping, bowtie2_recheck, bowtie2_index_prefix, .. } => {
            HostRemovalConfig::Kraken2 { db_path, threads, confidence_threshold, minimum_hit_groups, memory_mapping, bowtie2_recheck, bowtie2_index_prefix }
        }
        HostRemovalConfig::Minimap2 { index_path, .. } => {
            HostRemovalConfig::Minimap2 { index_path, threads }
        }
        HostRemovalConfig::Bowtie2 { index_prefix, .. } => {
            HostRemovalConfig::Bowtie2 { index_prefix, threads }
        }
        HostRemovalConfig::Sylph { db_path, .. } => {
            HostRemovalConfig::Sylph { db_path, threads }
        }
        HostRemovalConfig::Centrifuge { db_path, .. } => {
            HostRemovalConfig::Centrifuge { db_path, threads }
        }
        HostRemovalConfig::Auto {
            kraken2_db_path,
            bowtie2_index_prefix,
            host_pct_low_threshold,
            host_pct_high_threshold,
            reads_high_threshold,
            user_host_pct,
            survey,
            survey_n_reads,
            survey_threads,
            bowtie2_recheck,
            ..
        } => {
            HostRemovalConfig::Auto {
                kraken2_db_path,
                bowtie2_index_prefix,
                threads,
                host_pct_low_threshold,
                host_pct_high_threshold,
                reads_high_threshold,
                user_host_pct,
                survey,
                survey_n_reads,
                survey_threads,
                bowtie2_recheck,
            }
        }
    }
}

fn check_tools(host_removal: &HostRemovalConfig, skip_qc: bool) -> Result<()> {
    if !skip_qc {
        for tool in &["fastp"] {
            which::which(tool).map_err(|_| {
                anyhow::anyhow!(
                    "'{}' not found in PATH. Please install it first.",
                    tool
                )
            })?;
        }
    }

    // For auto mode, both bowtie2 and kraken2 must be available.
    if matches!(host_removal, HostRemovalConfig::Auto { .. }) {
        for tool in &["bowtie2", "kraken2"] {
            which::which(tool).map_err(|_| {
                anyhow::anyhow!(
                    "'{}' not found in PATH. Auto mode requires both bowtie2 and kraken2.",
                    tool
                )
            })?;
        }
        return Ok(());
    }

    let backend_tools: Vec<&str> = match host_removal {
        HostRemovalConfig::Kraken2 { bowtie2_recheck, .. } => {
            let mut tools = vec!["kraken2"];
            if *bowtie2_recheck {
                tools.extend(&["bowtie2", "samtools"]);
            }
            tools
        }
        HostRemovalConfig::Minimap2 { .. } => vec!["minimap2"],
        HostRemovalConfig::Bowtie2 { .. } => vec!["bowtie2", "samtools"],
        HostRemovalConfig::Sylph { .. } => vec!["sylph"],
        HostRemovalConfig::Centrifuge { .. } => vec!["centrifuge"],
        HostRemovalConfig::Auto { .. } => unreachable!(),
    };
    for tool in backend_tools {
        which::which(tool).map_err(|_| {
            anyhow::anyhow!(
                "'{}' not found in PATH. Please install it first.",
                tool
            )
        })?;
    }
    Ok(())
}
