use std::collections::HashSet;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use chrono::{DateTime, Utc};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};
use tokio::fs;
use tokio::process::Command;
use tracing::info;

use crate::config::{Config, HostRemovalConfig};
use crate::error::RustycleanError;
use crate::sample::Sample;

// ============================================================================
// Pipeline stages (ordered for comparison)
// Stages below are named historically after Kraken2 but are used generically
// for all host-removal backends (kraken2, minimap2, bowtie2).
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

/// Host-removal metrics.  Field names retain the Kraken2 heritage for
/// backwards compatibility of checkpoint files, but the values are filled in
/// uniformly for every backend:
///   - human_reads: reads identified as host
///   - unclassified_reads: reads kept after host removal (microbial + unclassified)
///   - classified_reads: same as human_reads (host reads are the "classified" target)
///   - contamination_percent: human / (human + kept) * 100
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
    /// Resolved backend when using auto mode (for verification/logging).
    #[serde(default)]
    pub auto_backend: Option<String>,
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
            auto_backend: None,
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
    // Stage 1: fastp (skipped when --skip-qc is requested)
    if checkpoint.stage < PipelineStage::FastpComplete {
        checkpoint.transition(PipelineStage::FastpRunning);
        if config.tools.fastp.skip_qc {
            info!(sample = %sample.id, "QC skipped by --skip-qc; using raw input reads");
            let metrics = build_skip_qc_metrics(sample)
                .await
                .map_err(|e| RustycleanError::ToolExecution(format!("failed to build skip-QC metrics: {}", e)))?;
            checkpoint.fastp_metrics = Some(metrics);
        } else {
            run_fastp(sample, checkpoint, config, work_dir).await?;
        }
        checkpoint.transition(PipelineStage::FastpComplete);
    }

    // Stage 2: host removal (backend selected by config.tools.host_removal)
    if checkpoint.stage < PipelineStage::Kraken2Complete {
        checkpoint.transition(PipelineStage::Kraken2Running);
        run_host_removal(sample, checkpoint, config, work_dir).await?;
        checkpoint.transition(PipelineStage::Kraken2Complete);
    }

    // Stage 3: validation & move output
    checkpoint.transition(PipelineStage::Validating);
    validate_and_finalize(sample, checkpoint, config).await?;
    checkpoint.transition(PipelineStage::Completed);

    Ok(())
}

async fn run_host_removal(
    sample: &Sample,
    checkpoint: &mut Checkpoint,
    config: &Config,
    work_dir: &Path,
) -> Result<()> {
    match &config.tools.host_removal {
        HostRemovalConfig::Kraken2 { .. } => run_kraken2(sample, checkpoint, config, work_dir).await,
        HostRemovalConfig::Minimap2 { .. } => run_minimap2(sample, checkpoint, config, work_dir).await,
        HostRemovalConfig::Bowtie2 { .. } => run_bowtie2(sample, checkpoint, config, work_dir).await,
        HostRemovalConfig::Sylph { .. } => run_sylph(sample, checkpoint, config, work_dir).await,
        HostRemovalConfig::Centrifuge { .. } => run_centrifuge(sample, checkpoint, config, work_dir).await,
        HostRemovalConfig::Auto { .. } => {
            let resolved = resolve_auto_config(sample, checkpoint, config, work_dir).await?;
            checkpoint.auto_backend = Some(resolved.mode().to_string());
            let mut resolved_config = config.clone();
            resolved_config.tools.host_removal = resolved;
            match &resolved_config.tools.host_removal {
                HostRemovalConfig::Kraken2 { .. } => run_kraken2(sample, checkpoint, &resolved_config, work_dir).await,
                HostRemovalConfig::Bowtie2 { .. } => run_bowtie2(sample, checkpoint, &resolved_config, work_dir).await,
                _ => bail!("internal error: auto mode resolved to unsupported backend: {}", resolved_config.tools.host_removal.mode()),
            }
        }
    }
}

// ============================================================================
// Auto backend resolution
// ============================================================================

/// Resolve the `auto` backend into a concrete kraken2 or bowtie2 config.
///
/// Decision order:
/// 1. If `user_host_pct` is provided, use it directly.
/// 2. Else if `survey` is enabled, sample N reads from the fastp-trimmed R1,
///    map with bowtie2, and estimate host percentage.
/// 3. Else fall back to a size-based heuristic using fastp input reads.
async fn resolve_auto_config(
    sample: &Sample,
    checkpoint: &Checkpoint,
    config: &Config,
    work_dir: &Path,
) -> Result<HostRemovalConfig> {
    let auto_cfg = match &config.tools.host_removal {
        HostRemovalConfig::Auto { .. } => config.tools.host_removal.clone(),
        _ => bail!("internal error: resolve_auto_config called with non-auto config"),
    };
    let (
        kraken2_db_path,
        bowtie2_index_prefix,
        threads,
        low_thr,
        high_thr,
        reads_thr,
        user_host_pct,
        survey,
        survey_n_reads,
        survey_threads,
        bowtie2_recheck,
    ) = match &auto_cfg {
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
        } => (
            kraken2_db_path.clone(),
            bowtie2_index_prefix.clone(),
            *threads,
            *host_pct_low_threshold,
            *host_pct_high_threshold,
            *reads_high_threshold,
            *user_host_pct,
            *survey,
            *survey_n_reads,
            *survey_threads,
            *bowtie2_recheck,
        ),
        _ => unreachable!(),
    };

    let fastp = checkpoint
        .fastp_metrics
        .as_ref()
        .context("fastp metrics missing before auto backend resolution")?;
    let input_reads = fastp.input_reads;

    let host_pct = if let Some(pct) = user_host_pct {
        info!(
            sample = %sample.id,
            host_pct = pct,
            "auto mode: using user-provided host percentage"
        );
        pct
    } else if survey {
        let survey_path = work_dir.join("auto_survey_R1.fastq.gz");
        sample_fastq_reads(&fastp.output_r1, &survey_path, survey_n_reads)?;
        let mapped = run_bowtie2_survey(&survey_path, &bowtie2_index_prefix, survey_threads, work_dir).await?;
        let surveyed = survey_n_reads.min(input_reads).max(1);
        let pct = (mapped as f64 / surveyed as f64) * 100.0;
        info!(
            sample = %sample.id,
            surveyed_reads = surveyed,
            mapped_reads = mapped,
            estimated_host_pct = format!("{:.2}", pct),
            "auto mode: completed survey"
        );
        pct
    } else {
        // Fallback: assume low-host for small samples so bowtie2 is used.
        let pct = if input_reads < reads_thr { low_thr / 2.0 } else { high_thr + 10.0 };
        info!(
            sample = %sample.id,
            input_reads = input_reads,
            estimated_host_pct = format!("{:.2}", pct),
            "auto mode: using size-based fallback"
        );
        pct
    };

    let chosen = choose_auto_backend(host_pct, input_reads, low_thr, high_thr, reads_thr);
    info!(
        sample = %sample.id,
        host_pct = format!("{:.2}", host_pct),
        input_reads = input_reads,
        chosen_backend = chosen,
        "auto mode: selected backend"
    );

    let resolved = match chosen {
        "kraken2" => HostRemovalConfig::Kraken2 {
            db_path: kraken2_db_path,
            threads,
            confidence_threshold: 0.0,
            minimum_hit_groups: 2,
            memory_mapping: false,
            bowtie2_recheck,
            bowtie2_index_prefix: Some(bowtie2_index_prefix.clone()),
        },
        _ => HostRemovalConfig::Bowtie2 {
            index_prefix: bowtie2_index_prefix,
            threads,
        },
    };

    Ok(resolved)
}

fn choose_auto_backend(
    host_pct: f64,
    input_reads: u64,
    low_threshold: f64,
    high_threshold: f64,
    reads_threshold: u64,
) -> &'static str {
    if host_pct < low_threshold {
        "bowtie2"
    } else if host_pct > high_threshold && input_reads > reads_threshold {
        "kraken2"
    } else {
        // Conservative default: bowtie2 unless clearly a large high-host sample.
        "bowtie2"
    }
}

/// Randomly sample `n_reads` from a (possibly gzipped) FASTQ file into `dst`.
/// Prefers `seqtk sample` when available for unbiased reservoir sampling;
/// falls back to an in-memory reservoir sample otherwise.
fn sample_fastq_reads(src: &Path, dst: &Path, n_reads: u64) -> Result<()> {
    if which::which("seqtk").is_ok() {
        let output = std::process::Command::new("seqtk")
            .arg("sample")
            .arg("-s")
            .arg("42")
            .arg(src)
            .arg(n_reads.to_string())
            .stdout(std::process::Stdio::piped())
            .output()
            .map_err(|e| RustycleanError::ToolExecution(format!("failed to run seqtk sample: {}", e)))?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(RustycleanError::ToolExecution(format!("seqtk sample failed: {}", stderr)).into());
        }

        let file = std::fs::File::create(dst)
            .map_err(|e| RustycleanError::ToolExecution(format!("failed to create survey FASTQ {}: {}", dst.display(), e)))?;
        let mut writer = GzEncoder::new(BufWriter::new(file), Compression::fast());
        writer.write_all(&output.stdout)
            .map_err(|e| RustycleanError::ToolExecution(format!("failed to write seqtk survey output: {}", e)))?;
        writer.finish()
            .map_err(|e| RustycleanError::ToolExecution(format!("failed to finalize survey gzip: {}", e)))?;
        return Ok(());
    }

    // Fallback: reservoir sampling (one pass, unbiased).
    let mut reader = open_fastq_reader(src)?;
    let file = std::fs::File::create(dst)
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to create survey FASTQ {}: {}", dst.display(), e)))?;
    let mut writer = GzEncoder::new(BufWriter::new(file), Compression::fast());

    let mut reservoir: Vec<Vec<u8>> = Vec::with_capacity(n_reads as usize);
    let mut record = Vec::with_capacity(1024);
    let mut line_count = 0u64;
    let mut read_index = 0u64;
    let mut line = Vec::new();
    let mut rng = StdRng::seed_from_u64(42);

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)
            .map_err(|e| RustycleanError::ToolExecution(format!("survey read error: {}", e)))?;
        if n == 0 {
            break;
        }
        record.extend_from_slice(&line);
        line_count += 1;

        if line_count == 4 {
            if read_index < n_reads {
                reservoir.push(record.clone());
            } else {
                let j = rng.gen_range(0..=read_index);
                if j < n_reads {
                    reservoir[j as usize] = record.clone();
                }
            }
            record.clear();
            line_count = 0;
            read_index += 1;
        }
    }

    for rec in reservoir {
        writer.write_all(&rec)
            .map_err(|e| RustycleanError::ToolExecution(format!("survey write error: {}", e)))?;
    }

    writer.finish()
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to finalize survey gzip: {}", e)))?;

    Ok(())
}

/// Run a quick bowtie2 survey on a small FASTQ and return the number of mapped reads.
async fn run_bowtie2_survey(
    reads: &Path,
    index_prefix: &Path,
    threads: usize,
    work_dir: &Path,
) -> Result<u64> {
    let sam_output = work_dir.join("auto_survey.sam");

    let output = Command::new("bowtie2")
        .arg("-x").arg(index_prefix)
        .arg("--very-fast-local")
        .arg("-p").arg(threads.to_string())
        .arg("-U").arg(reads)
        .arg("-S").arg(&sam_output)
        .output()
        .await
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to run bowtie2 survey: {}", e)))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(RustycleanError::ToolExecution(format!("bowtie2 survey failed: {}", stderr)).into());
    }

    let (mapped, _) = tokio::task::spawn_blocking({
        let sam_output = sam_output.clone();
        move || count_sam_mapped(&sam_output)
    })
    .await
    .map_err(|e| RustycleanError::ToolExecution(format!("failed to parse survey SAM: {}", e)))??;

    Ok(mapped)
}

fn count_sam_mapped(path: &Path) -> Result<(u64, u64)> {
    let file = std::fs::File::open(path)
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to open SAM {}: {}", path.display(), e)))?;
    let reader = BufReader::new(file);
    let mut mapped = 0u64;
    let mut unmapped = 0u64;

    for line in reader.lines() {
        let line = line.map_err(|e| RustycleanError::ToolExecution(format!("failed to read SAM {}: {}", path.display(), e)))?;
        if line.is_empty() || line.starts_with('@') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }
        let flag: u16 = fields[1].parse().unwrap_or(0);
        if flag & 0x4 == 0 {
            mapped += 1;
        } else {
            unmapped += 1;
        }
    }

    Ok((mapped, unmapped))
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

/// Build synthetic fastp metrics when QC is skipped.  The input read count is
/// required for AUTO mode backend selection and for consistent reporting.
async fn build_skip_qc_metrics(sample: &Sample) -> Result<FastpMetrics> {
    let r1 = sample.r1.clone();
    let input_reads = tokio::task::spawn_blocking(move || count_fastq_records(&r1))
        .await
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to count input reads: {}", e)))??;

    Ok(FastpMetrics {
        input_reads,
        output_reads: input_reads,
        // Quality metrics are unknown when fastp is skipped.
        q20_rate: 0.0,
        q30_rate: 0.0,
        gc_content: 0.0,
        adapter_trimmed: 0,
        too_short_reads: 0,
        low_quality_reads: 0,
        output_r1: sample.r1.clone(),
        output_r2: sample.r2.clone(),
    })
}

/// Count FASTQ records in a (possibly gzipped) file.  Uses seqtk when available
/// for speed, otherwise falls back to a streaming line counter.
fn count_fastq_records(path: &Path) -> Result<u64> {
    if let Ok(seqtk) = which::which("seqtk") {
        let output = std::process::Command::new(seqtk)
            .arg("comp")
            .arg(path)
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .output()
            .map_err(|e| RustycleanError::ToolExecution(format!("failed to run seqtk comp: {}", e)))?;
        if output.status.success() {
            let stdout = String::from_utf8_lossy(&output.stdout);
            let records = stdout.lines().filter(|l| !l.is_empty()).count() as u64;
            return Ok(records);
        }
    }

    // Fallback: count lines / 4
    let reader = open_fastq_reader(path)?;
    let mut lines = 0u64;
    for line in reader.lines() {
        line.map_err(|e| RustycleanError::ToolExecution(format!("failed to read {}: {}", path.display(), e)))?;
        lines += 1;
    }
    Ok(lines / 4)
}

// ============================================================================
// Read ID normalization
// ============================================================================

fn normalize_read_id(id: &str) -> String {
    let mut id = id.split('#').next().unwrap_or(id).to_string();
    if let Some((base, suffix)) = id.rsplit_once('/') {
        if suffix == "1" || suffix == "2" {
            id = base.to_string();
        }
    }
    id
}

// ============================================================================
// Host-removal backends
// ============================================================================

/// Shared helper: build final output paths, filter fastp outputs by removing
/// the read IDs in `human_ids`, and store metrics on the checkpoint.
async fn finalize_host_removal(
    sample: &Sample,
    checkpoint: &mut Checkpoint,
    human_ids: HashSet<String>,
    classified_reads: u64,
    unclassified_reads: u64,
) -> Result<()> {
    let fastp = checkpoint
        .fastp_metrics
        .as_ref()
        .context("fastp metrics missing before host-removal finalize")?;

    fs::create_dir_all(&sample.output_dir).await?;
    let final_r1 = sample.output_dir.join(format!("{}_clean_R1.fastq.gz", sample.id));
    let final_r2 = if sample.is_paired() {
        Some(sample.output_dir.join(format!("{}_clean_R2.fastq.gz", sample.id)))
    } else {
        None
    };

    let tmp_r1 = sample.output_dir.join(format!(".{}_clean_R1.fastq.gz.tmp", sample.id));
    let tmp_r2 = if sample.is_paired() {
        Some(sample.output_dir.join(format!(".{}_clean_R2.fastq.gz.tmp", sample.id)))
    } else {
        None
    };

    let mut filter_tasks: Vec<(PathBuf, PathBuf, HashSet<String>)> = vec![
        (fastp.output_r1.clone(), tmp_r1.clone(), human_ids.clone()),
    ];

    if let (Some(r2_src), Some(tmp_r2)) = (&fastp.output_r2, &tmp_r2) {
        filter_tasks.push((r2_src.clone(), tmp_r2.clone(), human_ids.clone()));
    }

    let mut filter_handles = Vec::new();
    for (src, dst, ids) in filter_tasks {
        let handle = tokio::task::spawn_blocking(move || filter_fastq_file(&src, &dst, &ids));
        filter_handles.push(handle);
    }

    for handle in filter_handles {
        handle
            .await
            .map_err(|e| RustycleanError::ToolExecution(format!("FASTQ filter task failed: {}", e)))??;
    }

    fs::rename(&tmp_r1, &final_r1).await?;
    if let (Some(tmp_r2), Some(final_r2)) = (&tmp_r2, &final_r2) {
        fs::rename(tmp_r2, final_r2).await?;
    }

    let total = classified_reads + unclassified_reads;
    let contamination = if total > 0 {
        (classified_reads as f64 / total as f64) * 100.0
    } else {
        0.0
    };

    checkpoint.kraken2_metrics = Some(Kraken2Metrics {
        classified_reads,
        unclassified_reads,
        human_reads: classified_reads,
        contamination_percent: contamination,
        output_r1: final_r1,
        output_r2: final_r2,
    });

    Ok(())
}

// ----------------------------------------------------------------------------
// kraken2
// ----------------------------------------------------------------------------

async fn run_kraken2(
    sample: &Sample,
    checkpoint: &mut Checkpoint,
    config: &Config,
    work_dir: &Path,
) -> Result<()> {
    let kraken_cfg = match &config.tools.host_removal {
        HostRemovalConfig::Kraken2 { db_path, threads, confidence_threshold, minimum_hit_groups, memory_mapping, bowtie2_recheck, bowtie2_index_prefix } => {
            (db_path.clone(), *threads, *confidence_threshold, *minimum_hit_groups, *memory_mapping, *bowtie2_recheck, bowtie2_index_prefix.clone())
        }
        _ => bail!("internal error: run_kraken2 called with non-kraken2 config"),
    };
    let (db_path, threads, confidence_threshold, minimum_hit_groups, memory_mapping, bowtie2_recheck, bowtie2_index_prefix) = kraken_cfg;

    let fastp = checkpoint
        .fastp_metrics
        .as_ref()
        .context("fastp metrics missing before kraken2 stage")?;

    let kraken_report = work_dir.join("kraken2.report");
    let kraken_output = work_dir.join("kraken2.output.txt");

    let mut cmd = Command::new("kraken2");
    cmd.arg("--db")
        .arg(&db_path)
        .arg("--threads")
        .arg(threads.to_string())
        .arg("--confidence")
        .arg(confidence_threshold.to_string())
        .arg("--minimum-hit-groups")
        .arg(minimum_hit_groups.to_string())
        .arg("--report")
        .arg(&kraken_report)
        .arg("--output")
        .arg(&kraken_output)
        .arg("--gzip-compressed");

    if memory_mapping {
        cmd.arg("--memory-mapping");
    }

    if sample.is_paired() {
        cmd.arg("--paired")
            .arg(&fastp.output_r1)
            .arg(fastp.output_r2.as_ref().unwrap());
    } else {
        cmd.arg(&fastp.output_r1);
    }

    let output = cmd
        .output()
        .await
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to run kraken2: {}", e)))?;

    let stderr = String::from_utf8_lossy(&output.stderr);
    if stderr.contains("Error:") || stderr.contains("FATAL") {
        return Err(RustycleanError::ToolExecution(stderr.to_string()).into());
    }

    // Parse Kraken2 per-read output to identify human and unclassified reads
    let classification = tokio::task::spawn_blocking({
        let kraken_output = kraken_output.clone();
        move || parse_kraken_output(&kraken_output)
    })
    .await
    .map_err(|e| RustycleanError::ToolExecution(format!("failed to parse kraken2 output: {}", e)))??;

    let mut human_ids = classification.human_ids;

    // Optional Bowtie2 re-check of Kraken2-unclassified reads against the host index.
    if bowtie2_recheck {
        let index_prefix = bowtie2_index_prefix
            .as_ref()
            .context("--bowtie2-recheck requires --host-index (bowtie2 index prefix)")?;
        info!(
            sample = %sample.id,
            unclassified_reads = classification.unclassified_ids.len(),
            "running Bowtie2 re-check on Kraken2-unclassified reads"
        );
        let recheck_human_ids = run_bowtie2_recheck(
            sample,
            &fastp.output_r1,
            fastp.output_r2.as_deref(),
            &classification.unclassified_ids,
            index_prefix,
            threads,
            work_dir,
        ).await?;
        let additional_host = recheck_human_ids.len() as u64;
        human_ids.extend(recheck_human_ids);
        info!(
            sample = %sample.id,
            additional_host_reads = additional_host,
            "Bowtie2 re-check complete"
        );
    }

    // Derive counts from the report for metrics.
    // We report consistently with alignment backends:
    //   human_reads     = reads classified as Homo sapiens (host)
    //   unclassified    = reads kept after filtering (microbial + unclassified)
    let (_classified_total, _unclassified_total, human) = tokio::task::spawn_blocking({
        let kraken_report = kraken_report.clone();
        move || parse_kraken_report_counts(&kraken_report)
    })
    .await
    .map_err(|e| RustycleanError::ToolExecution(format!("failed to parse kraken2 report: {}", e)))??;

    // If recheck added host reads, adjust the reported human count and unclassified count.
    let reported_human = if bowtie2_recheck {
        human_ids.len() as u64
    } else {
        human
    };
    let kept = fastp.output_reads.saturating_sub(reported_human);
    finalize_host_removal(sample, checkpoint, human_ids, reported_human, kept).await?;

    Ok(())
}

/// Re-align Kraken2-unclassified reads with Bowtie2 against the host index and
/// return the read IDs that map to the host.
async fn run_bowtie2_recheck(
    _sample: &Sample,
    r1: &Path,
    r2: Option<&Path>,
    unclassified_ids: &HashSet<String>,
    index_prefix: &Path,
    threads: usize,
    work_dir: &Path,
) -> Result<HashSet<String>> {
    if unclassified_ids.is_empty() {
        return Ok(HashSet::new());
    }

    // Extract unclassified reads from the QC-filtered input.
    let recheck_r1 = work_dir.join("recheck_R1.fastq.gz");
    let recheck_r2 = r2.map(|_| work_dir.join("recheck_R2.fastq.gz"));

    let ids = unclassified_ids.clone();
    let r1_src = r1.to_path_buf();
    let r1_dst = recheck_r1.clone();
    tokio::task::spawn_blocking(move || extract_fastq_reads(&r1_src, &r1_dst, &ids))
        .await
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to extract recheck R1 reads: {}", e)))??;

    if let (Some(r2_src), Some(r2_dst)) = (r2, recheck_r2.as_ref()) {
        let ids = unclassified_ids.clone();
        let r2_src = r2_src.to_path_buf();
        let r2_dst = r2_dst.clone();
        tokio::task::spawn_blocking(move || extract_fastq_reads(&r2_src, &r2_dst, &ids))
            .await
            .map_err(|e| RustycleanError::ToolExecution(format!("failed to extract recheck R2 reads: {}", e)))??;
    }

    // Run Bowtie2 on the extracted reads.
    let sam_output = work_dir.join("bowtie2_recheck.sam");
    let mut cmd = Command::new("bowtie2");
    cmd.arg("-x").arg(index_prefix)
        .arg("--very-fast-local")
        .arg("-p").arg(threads.to_string());

    if let Some(r2_path) = &recheck_r2 {
        cmd.arg("-1").arg(&recheck_r1)
            .arg("-2").arg(r2_path);
    } else {
        cmd.arg("-U").arg(&recheck_r1);
    }

    cmd.arg("-S").arg(&sam_output);

    let output = cmd
        .output()
        .await
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to run bowtie2 recheck: {}", e)))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(RustycleanError::ToolExecution(format!("bowtie2 recheck failed: {}", stderr)).into());
    }

    // Parse SAM to obtain mapped IDs.
    let (mapped_ids, _, _) = tokio::task::spawn_blocking({
        let sam_output = sam_output.clone();
        move || parse_sam_mapped_ids(&sam_output)
    })
    .await
    .map_err(|e| RustycleanError::ToolExecution(format!("failed to parse bowtie2 recheck SAM: {}", e)))??;

    Ok(mapped_ids)
}

/// Parsed per-read Kraken2 output.
#[derive(Debug, Default)]
struct Kraken2ReadClassification {
    /// Read IDs classified as Homo sapiens (taxid 9606).
    human_ids: HashSet<String>,
    /// Read IDs classified as unclassified (status "U").
    unclassified_ids: HashSet<String>,
}

fn parse_kraken_output(path: &Path) -> Result<Kraken2ReadClassification> {
    let file = std::fs::File::open(path)
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to open kraken2 output: {}", e)))?;
    let reader = BufReader::new(file);
    let mut result = Kraken2ReadClassification::default();

    for line in reader.lines() {
        let line = line.map_err(|e| RustycleanError::ToolExecution(format!("failed to read kraken2 output: {}", e)))?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        let status = fields[0];
        let read_id = fields[1];
        let taxid: i64 = fields[2].parse().unwrap_or(0);
        let normalized = normalize_read_id(read_id);

        if status == "C" && taxid == 9606 {
            result.human_ids.insert(normalized);
        } else if status == "U" {
            result.unclassified_ids.insert(normalized);
        }
    }

    Ok(result)
}

fn parse_kraken_report_counts(path: &Path) -> Result<(u64, u64, u64)> {
    let file = std::fs::File::open(path)
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to open kraken2 report: {}", e)))?;
    let reader = BufReader::new(file);
    let mut classified = 0u64;
    let mut unclassified = 0u64;
    let mut human = 0u64;

    for line in reader.lines() {
        let line = line.map_err(|e| RustycleanError::ToolExecution(format!("failed to read kraken2 report: {}", e)))?;
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

    Ok((classified, unclassified, human))
}

// ----------------------------------------------------------------------------
// minimap2
// ----------------------------------------------------------------------------

async fn run_minimap2(
    sample: &Sample,
    checkpoint: &mut Checkpoint,
    config: &Config,
    work_dir: &Path,
) -> Result<()> {
    let (index_path, threads) = match &config.tools.host_removal {
        HostRemovalConfig::Minimap2 { index_path, threads } => (index_path.clone(), *threads),
        _ => bail!("internal error: run_minimap2 called with non-minimap2 config"),
    };

    let fastp = checkpoint
        .fastp_metrics
        .as_ref()
        .context("fastp metrics missing before minimap2 stage")?;

    let sam_output = work_dir.join("minimap2.sam");

    let mut cmd = Command::new("minimap2");
    cmd.arg("-x").arg("sr")
        .arg("-t").arg(threads.to_string())
        .arg("-a") // output SAM
        .arg("--secondary=no")
        .arg("-o").arg(&sam_output)
        .arg(&index_path)
        .arg(&fastp.output_r1);

    if sample.is_paired() {
        cmd.arg(fastp.output_r2.as_ref().unwrap());
    }

    let output = cmd
        .output()
        .await
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to run minimap2: {}", e)))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(RustycleanError::ToolExecution(format!("minimap2 failed: {}", stderr)).into());
    }

    let (mapped_ids, mapped, unmapped) = tokio::task::spawn_blocking({
        let sam_output = sam_output.clone();
        move || parse_sam_mapped_ids(&sam_output)
    })
    .await
    .map_err(|e| RustycleanError::ToolExecution(format!("failed to parse minimap2 SAM: {}", e)))??;

    finalize_host_removal(sample, checkpoint, mapped_ids, mapped, unmapped).await?;

    Ok(())
}

// ----------------------------------------------------------------------------
// bowtie2
// ----------------------------------------------------------------------------

async fn run_bowtie2(
    sample: &Sample,
    checkpoint: &mut Checkpoint,
    config: &Config,
    work_dir: &Path,
) -> Result<()> {
    let (index_prefix, threads) = match &config.tools.host_removal {
        HostRemovalConfig::Bowtie2 { index_prefix, threads } => (index_prefix.clone(), *threads),
        _ => bail!("internal error: run_bowtie2 called with non-bowtie2 config"),
    };

    let fastp = checkpoint
        .fastp_metrics
        .as_ref()
        .context("fastp metrics missing before bowtie2 stage")?;

    fs::create_dir_all(&sample.output_dir).await?;
    let final_r1 = sample.output_dir.join(format!("{}_clean_R1.fastq.gz", sample.id));
    let final_r2 = if sample.is_paired() {
        Some(sample.output_dir.join(format!("{}_clean_R2.fastq.gz", sample.id)))
    } else {
        None
    };

    let count_before = work_dir.join("bowtie2_count_before.txt");
    let count_after = work_dir.join("bowtie2_count_after.txt");

    let index = index_prefix.display().to_string();
    let r1 = fastp.output_r1.display().to_string();
    let t = threads.to_string();
    let c_before = count_before.display().to_string();
    let c_after = count_after.display().to_string();

    let pipeline = if sample.is_paired() {
        let r2 = fastp.output_r2.as_ref().unwrap().display().to_string();
        let out1 = final_r1.display().to_string();
        let out2 = final_r2.as_ref().unwrap().display().to_string();
        format!(
            concat!(
                "set -euo pipefail; ",
                "bowtie2 -x \"{}\" -1 \"{}\" -2 \"{}\" --very-fast-local -p {} -k 1 --mm ",
                "  | tee >(samtools view -F 2304 -c - > \"{}\") ",
                "  | samtools view -hf 12 - ",
                "  | tee >(samtools view -F 2304 -c - > \"{}\") ",
                "  | samtools fastq --threads {} -c 6 -1 \"{}\" -2 \"{}\" -0 /dev/null -s /dev/null"
            ),
            index, r1, r2, t, c_before, c_after, t, out1, out2
        )
    } else {
        let out1 = final_r1.display().to_string();
        format!(
            concat!(
                "set -euo pipefail; ",
                "bowtie2 -x \"{}\" -U \"{}\" --very-fast-local -p {} -k 1 --mm ",
                "  | tee >(samtools view -F 2304 -c - > \"{}\") ",
                "  | samtools view -hf 4 - ",
                "  | tee >(samtools view -F 2304 -c - > \"{}\") ",
                "  | samtools fastq --threads {} -c 6 -0 \"{}\""
            ),
            index, r1, t, c_before, c_after, t, out1
        )
    };
    let output = Command::new("bash")
        .arg("-c")
        .arg(&pipeline)
        .output()
        .await
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to run bowtie2 streaming pipeline: {}", e)))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(RustycleanError::ToolExecution(format!("bowtie2 streaming pipeline failed: {}", stderr)).into());
    }

    let before = read_count_file(&count_before).await?;
    let after = read_count_file(&count_after).await?;
    let classified_reads = before.saturating_sub(after);
    let unclassified_reads = after;
    let total = before;
    let contamination = if total > 0 {
        (classified_reads as f64 / total as f64) * 100.0
    } else {
        0.0
    };

    checkpoint.kraken2_metrics = Some(Kraken2Metrics {
        classified_reads,
        unclassified_reads,
        human_reads: classified_reads,
        contamination_percent: contamination,
        output_r1: final_r1,
        output_r2: final_r2,
    });

    validate_and_finalize(sample, checkpoint, config).await?;

    Ok(())
}

async fn read_count_file(path: &Path) -> Result<u64> {
    let content = fs::read_to_string(path)
        .await
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to read count file {}: {}", path.display(), e)))?;
    content
        .trim()
        .parse::<u64>()
        .map_err(|e| RustycleanError::ToolExecution(format!("invalid count in {}: {}", path.display(), e)).into())
}


// ----------------------------------------------------------------------------
// centrifuge
// ----------------------------------------------------------------------------

const HUMAN_TAXID: u64 = 9606;

async fn run_centrifuge(
    sample: &Sample,
    checkpoint: &mut Checkpoint,
    config: &Config,
    work_dir: &Path,
) -> Result<()> {
    let (db_path, threads) = match &config.tools.host_removal {
        HostRemovalConfig::Centrifuge { db_path, threads } => (db_path.clone(), *threads),
        _ => bail!("internal error: run_centrifuge called with non-centrifuge config"),
    };

    let fastp = checkpoint
        .fastp_metrics
        .as_ref()
        .context("fastp metrics missing before centrifuge stage")?;

    let class_output = work_dir.join("centrifuge_classifications.tsv");
    let report_output = work_dir.join("centrifuge_report.tsv");

    let mut cmd = Command::new("centrifuge");
    cmd.arg("-x").arg(&db_path)
        .arg("-S").arg(&class_output)
        .arg("--report-file").arg(&report_output)
        .arg("--threads").arg(threads.to_string());

    if sample.is_paired() {
        cmd.arg("-1").arg(&fastp.output_r1)
            .arg("-2").arg(fastp.output_r2.as_ref().unwrap());
    } else {
        cmd.arg("-U").arg(&fastp.output_r1);
    }

    let output = cmd
        .output()
        .await
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to run centrifuge: {}", e)))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(RustycleanError::ToolExecution(format!("centrifuge failed: {}", stderr)).into());
    }

    let human_ids = tokio::task::spawn_blocking({
        let class_output = class_output.clone();
        move || parse_centrifuge_classifications(&class_output)
    })
    .await
    .map_err(|e| RustycleanError::ToolExecution(format!("failed to parse centrifuge output: {}", e)))??;

    let classified_reads = human_ids.len() as u64;
    let total_reads = fastp.output_reads;
    let unclassified_reads = total_reads.saturating_sub(classified_reads);

    finalize_host_removal(sample, checkpoint, human_ids, classified_reads, unclassified_reads).await?;

    Ok(())
}

fn parse_centrifuge_classifications(path: &Path) -> Result<HashSet<String>> {
    let file = std::fs::File::open(path)
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to open centrifuge classifications {}: {}", path.display(), e)))?;
    let reader = BufReader::new(file);
    let mut human_ids = HashSet::new();

    for (idx, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| RustycleanError::ToolExecution(format!("read error in centrifuge classifications: {}", e)))?;
        if idx == 0 || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        if let Ok(taxid) = fields[2].parse::<u64>() {
            if taxid == HUMAN_TAXID {
                human_ids.insert(normalize_read_id(fields[0]));
            }
        }
    }

    Ok(human_ids)
}

// ----------------------------------------------------------------------------
// sylph
// ----------------------------------------------------------------------------

async fn run_sylph(
    _sample: &Sample,
    _checkpoint: &mut Checkpoint,
    _config: &Config,
    _work_dir: &Path,
) -> Result<()> {
    // sylph 0.9.x does not provide per-read classification.  It can only
    // estimate sample-level containment / abundance.  Read-level host removal
    // would require a different sylph version or another tool.
    bail!(
        "sylph backend is not supported for read-level host removal: installed sylph lacks a per-read classify command. Use kraken2, minimap2, or bowtie2 instead."
    )
}

// ----------------------------------------------------------------------------
// SAM parsing shared by minimap2 / bowtie2
// ----------------------------------------------------------------------------

fn parse_sam_mapped_ids(path: &Path) -> Result<(HashSet<String>, u64, u64)> {
    let file = std::fs::File::open(path)
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to open SAM {}: {}", path.display(), e)))?;
    let reader = BufReader::new(file);
    let mut mapped_ids = HashSet::new();
    let mut mapped = 0u64;
    let mut unmapped = 0u64;

    for line in reader.lines() {
        let line = line.map_err(|e| RustycleanError::ToolExecution(format!("failed to read SAM {}: {}", path.display(), e)))?;
        if line.is_empty() || line.starts_with('@') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }
        let read_id = fields[0];
        let flag: u16 = fields[1].parse().unwrap_or(0);

        // 0x4 = READ_UNMAPPED
        if flag & 0x4 == 0 {
            mapped_ids.insert(normalize_read_id(read_id));
            mapped += 1;
        } else {
            unmapped += 1;
        }
    }

    Ok((mapped_ids, mapped, unmapped))
}

// ============================================================================
// FASTQ filtering
// ============================================================================

fn open_fastq_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = std::fs::File::open(path)
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to open FASTQ {}: {}", path.display(), e)))?;
    if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn filter_fastq_file(src: &Path, dst: &Path, human_ids: &HashSet<String>) -> Result<()> {
    let mut reader = open_fastq_reader(src)?;
    let file = std::fs::File::create(dst)
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to create FASTQ {}: {}", dst.display(), e)))?;
    let mut writer = GzEncoder::new(BufWriter::new(file), Compression::fast());

    let mut record = Vec::with_capacity(1024);
    let mut line_count = 0;
    let mut line = Vec::new();

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)
            .map_err(|e| RustycleanError::ToolExecution(format!("read error: {}", e)))?;
        if n == 0 {
            if line_count > 0 {
                return Err(RustycleanError::ToolExecution("unexpected EOF in FASTQ".to_string()).into());
            }
            break;
        }
        record.extend_from_slice(&line);
        line_count += 1;

        if line_count == 4 {
            // Parse the first line to obtain the read ID
            let first_end = record.iter().position(|&b| b == b'\n').unwrap_or(record.len());
            let first = std::str::from_utf8(&record[..first_end])
                .map_err(|e| RustycleanError::ToolExecution(format!("invalid UTF-8 in FASTQ name: {}", e)))?;
            let read_id = first.strip_prefix('@').unwrap_or(first);
            let read_id = read_id.split_whitespace().next().unwrap_or(read_id);
            let normalized = normalize_read_id(read_id);

            if !human_ids.contains(&normalized) {
                writer.write_all(&record)
                    .map_err(|e| RustycleanError::ToolExecution(format!("write error: {}", e)))?;
            }

            record.clear();
            line_count = 0;
        }
    }

    writer.finish()
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to finalize gzip: {}", e)))?;

    Ok(())
}

/// Extract reads whose normalized ID is in `keep_ids` from `src` to `dst`.
fn extract_fastq_reads(src: &Path, dst: &Path, keep_ids: &HashSet<String>) -> Result<()> {
    let mut reader = open_fastq_reader(src)?;
    let file = std::fs::File::create(dst)
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to create FASTQ {}: {}", dst.display(), e)))?;
    let mut writer = GzEncoder::new(BufWriter::new(file), Compression::fast());

    let mut record = Vec::with_capacity(1024);
    let mut line_count = 0;
    let mut line = Vec::new();

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)
            .map_err(|e| RustycleanError::ToolExecution(format!("read error: {}", e)))?;
        if n == 0 {
            if line_count > 0 {
                return Err(RustycleanError::ToolExecution("unexpected EOF in FASTQ".to_string()).into());
            }
            break;
        }
        record.extend_from_slice(&line);
        line_count += 1;

        if line_count == 4 {
            let first_end = record.iter().position(|&b| b == b'\n').unwrap_or(record.len());
            let first = std::str::from_utf8(&record[..first_end])
                .map_err(|e| RustycleanError::ToolExecution(format!("invalid UTF-8 in FASTQ name: {}", e)))?;
            let read_id = first.strip_prefix('@').unwrap_or(first);
            let read_id = read_id.split_whitespace().next().unwrap_or(read_id);
            let normalized = normalize_read_id(read_id);

            if keep_ids.contains(&normalized) {
                writer.write_all(&record)
                    .map_err(|e| RustycleanError::ToolExecution(format!("write error: {}", e)))?;
            }

            record.clear();
            line_count = 0;
        }
    }

    writer.finish()
        .map_err(|e| RustycleanError::ToolExecution(format!("failed to finalize gzip: {}", e)))?;

    Ok(())
}

// Kept for backwards compatibility of any external callers/tests.
#[allow(dead_code)]
async fn parse_kraken_report(
    report_path: &Path,
    output_r1: PathBuf,
    output_r2: Option<PathBuf>,
) -> Result<Kraken2Metrics> {
    let (classified, unclassified, human) = parse_kraken_report_counts(report_path)?;
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
        .context("host-removal metrics missing at validation stage")?;

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

    // Final outputs are already written to sample.output_dir by the backend.
    // Validate that they exist and are non-empty.
    if !fs::metadata(&kraken.output_r1).await.map(|m| m.len() > 0).unwrap_or(false) {
        return Err(RustycleanError::ValidationFailed(sample.id.clone(), "R1 output is missing or empty".to_string()).into());
    }

    if let Some(r2) = &kraken.output_r2 {
        if !fs::metadata(r2).await.map(|m| m.len() > 0).unwrap_or(false) {
            return Err(RustycleanError::ValidationFailed(sample.id.clone(), "R2 output is missing or empty".to_string()).into());
        }
    }

    info!(
        sample = %sample.id,
        unclassified_reads = kraken.unclassified_reads,
        human_reads = kraken.human_reads,
        contamination = format!("{:.2}%", kraken.contamination_percent),
        auto_backend = ?checkpoint.auto_backend,
        "Sample validated and finalized"
    );

    Ok(())
}
