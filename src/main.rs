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
use crate::cli::Cli;
use crate::config::{Config, ExecutionConfig, FastpConfig, Kraken2Config, ToolConfig, ValidationConfig};
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
                kraken2: Kraken2Config::default(),
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
        config.tools.kraken2.threads = t;
    }
    if let Some(db) = &cli.kraken2_db {
        config.tools.kraken2.db_path = db.clone();
    }
    if cli.resume {
        config.execution.resume = true;
    }

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

    info!("Loaded {} sample(s)", samples.len());

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
    check_tools()?;

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

fn check_tools() -> Result<()> {
    for tool in &["fastp", "kraken2"] {
        which::which(tool).map_err(|_| {
            anyhow::anyhow!(
                "'{}' not found in PATH. Please install it first.",
                tool
            )
        })?;
    }
    Ok(())
}
