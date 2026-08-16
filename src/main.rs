mod checkpoint;
mod cli;
mod config;
mod error;
mod pipeline;
mod sample;
mod worker;

use std::path::PathBuf;
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
    let workers_from_cli = cli.workers.is_some();
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
            let db_path = cli.sylph_db
                .unwrap_or_else(|| std::path::PathBuf::from("/db/human_t2t.syldb"));
            let bowtie2_index_prefix = host_index
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| std::path::PathBuf::from("/db/human_t2t_hla"));

            if !db_path.exists() {
                bail!("--sylph-db path does not exist: {}", db_path.display());
            }
            let bt2_test = bowtie2_index_prefix.with_extension("1.bt2");
            if !bt2_test.exists() {
                bail!("bowtie2 index not found at prefix for sylph backend: {}", bowtie2_index_prefix.display());
            }

            HostRemovalConfig::Sylph {
                db_path,
                bowtie2_index_prefix,
                threads: config.tools.host_removal.threads(),
                min_ani: cli.sylph_min_ani,
                min_eff_cov: cli.sylph_min_cov,
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
            let sylph_db_path = cli.sylph_db
                .unwrap_or_else(|| {
                    match &config.tools.host_removal {
                        HostRemovalConfig::Sylph { db_path, .. } => db_path.clone(),
                        HostRemovalConfig::Auto { sylph_db_path, .. } => sylph_db_path.clone(),
                        _ => std::path::PathBuf::from("/db/human_t2t.syldb"),
                    }
                });
            let bowtie2_index_prefix = host_index
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| {
                    match &config.tools.host_removal {
                        HostRemovalConfig::Bowtie2 { index_prefix, .. } => index_prefix.clone(),
                        HostRemovalConfig::Sylph { bowtie2_index_prefix, .. } => bowtie2_index_prefix.clone(),
                        HostRemovalConfig::Auto { bowtie2_index_prefix, .. } => bowtie2_index_prefix.clone(),
                        _ => std::path::PathBuf::from("/db/human_t2t_hla"),
                    }
                });
            // Optional Kraken2 database for explicit fallback. If not provided,
            // auto mode simply cannot fall back to kraken2.
            let kraken2_db_path = cli.kraken2_db.or_else(|| {
                match &config.tools.host_removal {
                    HostRemovalConfig::Kraken2 { db_path, .. } => Some(db_path.clone()),
                    HostRemovalConfig::Auto { kraken2_db_path, .. } => kraken2_db_path.clone(),
                    _ => None,
                }
            });

            // Validate required databases for auto mode.
            if !sylph_db_path.exists() {
                bail!("--sylph-db path does not exist: {}", sylph_db_path.display());
            }
            let bt2_test = bowtie2_index_prefix.with_extension("1.bt2");
            if !bt2_test.exists() {
                bail!("bowtie2 index not found at prefix: {}", bowtie2_index_prefix.display());
            }
            if let Some(ref kdb) = kraken2_db_path {
                if !kdb.exists() {
                    bail!("--kraken2-db path does not exist: {}", kdb.display());
                }
            }

            HostRemovalConfig::Auto {
                sylph_db_path,
                bowtie2_index_prefix,
                kraken2_db_path,
                threads: config.tools.host_removal.threads(),
                host_pct_low_threshold: cli.auto_low_threshold,
                host_pct_high_threshold: cli.auto_high_threshold,
                reads_high_threshold: cli.auto_reads_threshold,
                user_host_pct: cli.host_pct,
                survey: cli.auto_survey,
                survey_n_reads: cli.auto_survey_nreads,
                survey_threads: cli.auto_survey_threads,
                sylph_min_ani: cli.sylph_min_ani,
                sylph_min_eff_cov: cli.sylph_min_cov,
                memory_mapping: cli.kraken2_memory_mapping,
                bowtie2_recheck: cli.bowtie2_recheck,
            }
        }
    };

    // Cap worker count by available memory when the user did not explicitly set it.
    // Each worker loads its own database copy, so concurrent memory demand scales
    // with worker count. This prevents OOM on memory-constrained nodes.
    if !workers_from_cli {
        let cpu_workers = config.execution.workers;
        let mem_capped = memory_cap_workers(&config.tools.host_removal, cpu_workers);
        if mem_capped != cpu_workers {
            info!(
                "Capping parallel workers by available memory: {} -> {}",
                cpu_workers, mem_capped
            );
            config.execution.workers = mem_capped;
        }
    }

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

/// Return available memory in kB, preferring cgroup limits over /proc/meminfo.
fn available_memory_kb() -> Option<u64> {
    // cgroup v2
    if let Ok(s) = std::fs::read_to_string("/sys/fs/cgroup/memory.max") {
        if let Ok(bytes) = s.trim().parse::<u64>() {
            return Some(bytes / 1024);
        }
    }
    // cgroup v1
    if let Ok(s) = std::fs::read_to_string("/sys/fs/cgroup/memory/memory.limit_in_bytes") {
        if let Ok(bytes) = s.trim().parse::<u64>() {
            return Some(bytes / 1024);
        }
    }
    // Fallback to system memory
    if let Ok(s) = std::fs::read_to_string("/proc/meminfo") {
        for line in s.lines() {
            if line.starts_with("MemAvailable:") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let Ok(kb) = parts[1].parse::<u64>() {
                        return Some(kb);
                    }
                }
            }
        }
    }
    None
}

/// Estimate the resident database size in kB for the chosen backend.
fn estimate_db_size_kb(host_removal: &HostRemovalConfig) -> Option<u64> {
    let paths: Vec<PathBuf> = match host_removal {
        HostRemovalConfig::Kraken2 { db_path, .. } => {
            vec![db_path.join("hash.k2d")]
        }
        HostRemovalConfig::Auto { sylph_db_path, bowtie2_index_prefix, .. } => {
            // Auto mode's memory peak is dominated by the sylph+bowtie2 branch.
            let mut paths = vec![sylph_db_path.clone()];
            paths.extend([
                bowtie2_index_prefix.with_extension("1.bt2"),
                bowtie2_index_prefix.with_extension("2.bt2"),
                bowtie2_index_prefix.with_extension("3.bt2"),
                bowtie2_index_prefix.with_extension("4.bt2"),
                bowtie2_index_prefix.with_extension("rev.1.bt2"),
                bowtie2_index_prefix.with_extension("rev.2.bt2"),
            ]);
            paths
        }
        HostRemovalConfig::Bowtie2 { index_prefix, .. } => {
            vec![
                index_prefix.with_extension("1.bt2"),
                index_prefix.with_extension("2.bt2"),
                index_prefix.with_extension("3.bt2"),
                index_prefix.with_extension("4.bt2"),
                index_prefix.with_extension("rev.1.bt2"),
                index_prefix.with_extension("rev.2.bt2"),
            ]
        }
        HostRemovalConfig::Minimap2 { index_path, .. } => vec![index_path.clone()],
        HostRemovalConfig::Centrifuge { db_path, .. } => {
            vec![
                db_path.with_extension("1.cf"),
                db_path.with_extension("2.cf"),
                db_path.with_extension("3.cf"),
            ]
        }
        HostRemovalConfig::Sylph { db_path, bowtie2_index_prefix, .. } => {
            // Worst-case memory for the sylph+bowtie2 pipeline: both databases
            // may need to be resident when sylph signals host presence.
            let mut paths = vec![db_path.clone()];
            paths.extend([
                bowtie2_index_prefix.with_extension("1.bt2"),
                bowtie2_index_prefix.with_extension("2.bt2"),
                bowtie2_index_prefix.with_extension("3.bt2"),
                bowtie2_index_prefix.with_extension("4.bt2"),
                bowtie2_index_prefix.with_extension("rev.1.bt2"),
                bowtie2_index_prefix.with_extension("rev.2.bt2"),
            ]);
            paths
        }
    };

    let mut total_bytes: u64 = 0;
    for p in paths {
        if let Ok(m) = std::fs::metadata(&p) {
            total_bytes += m.len();
        }
    }

    if total_bytes == 0 {
        None
    } else {
        Some(total_bytes / 1024)
    }
}

/// Cap worker count so that concurrent database copies fit in available memory.
fn memory_cap_workers(host_removal: &HostRemovalConfig, cpu_workers: usize) -> usize {
    let Some(db_kb) = estimate_db_size_kb(host_removal) else {
        return cpu_workers;
    };
    let Some(avail_kb) = available_memory_kb() else {
        return cpu_workers;
    };

    // Reserve 20% headroom for the OS, tool overhead, and intermediate files.
    let usable_kb = (avail_kb as f64 * 0.8).max(1.0);
    let max_by_mem = (usable_kb / db_kb as f64).max(1.0) as usize;

    cpu_workers.min(max_by_mem)
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
        HostRemovalConfig::Sylph { db_path, bowtie2_index_prefix, min_ani, min_eff_cov, .. } => {
            HostRemovalConfig::Sylph { db_path, bowtie2_index_prefix, threads, min_ani, min_eff_cov }
        }
        HostRemovalConfig::Centrifuge { db_path, .. } => {
            HostRemovalConfig::Centrifuge { db_path, threads }
        }
        HostRemovalConfig::Auto {
            sylph_db_path,
            bowtie2_index_prefix,
            kraken2_db_path,
            host_pct_low_threshold,
            host_pct_high_threshold,
            reads_high_threshold,
            user_host_pct,
            survey,
            survey_n_reads,
            survey_threads,
            sylph_min_ani,
            sylph_min_eff_cov,
            memory_mapping,
            bowtie2_recheck,
            ..
        } => {
            HostRemovalConfig::Auto {
                sylph_db_path,
                bowtie2_index_prefix,
                kraken2_db_path,
                threads,
                host_pct_low_threshold,
                host_pct_high_threshold,
                reads_high_threshold,
                user_host_pct,
                survey,
                survey_n_reads,
                survey_threads,
                sylph_min_ani,
                sylph_min_eff_cov,
                memory_mapping,
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

    // For auto mode, sylph and bowtie2 are required; kraken2 is only needed
    // when an explicit fallback is configured.
    if matches!(host_removal, HostRemovalConfig::Auto { .. }) {
        for tool in &["sylph", "bowtie2", "samtools"] {
            which::which(tool).map_err(|_| {
                anyhow::anyhow!(
                    "'{}' not found in PATH. Auto mode requires sylph, bowtie2 and samtools.",
                    tool
                )
            })?;
        }
        if let HostRemovalConfig::Auto { kraken2_db_path: Some(_), .. } = host_removal {
            which::which("kraken2").map_err(|_| {
                anyhow::anyhow!("'kraken2' not found in PATH but --kraken2-db was provided for auto mode fallback.")
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
        HostRemovalConfig::Sylph { .. } => vec!["sylph", "bowtie2", "samtools"],
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
