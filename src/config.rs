use std::path::PathBuf;

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Config {
    pub execution: ExecutionConfig,
    pub tools: ToolConfig,
    pub validation: ValidationConfig,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ExecutionConfig {
    pub workers: usize,
    pub resume: bool,
    pub retry_attempts: u32,
    pub sample_timeout_minutes: u64,
}

impl Default for ExecutionConfig {
    fn default() -> Self {
        Self {
            workers: (num_cpus::get() / 2).max(1),
            resume: true,
            retry_attempts: 2,
            sample_timeout_minutes: 120,
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ToolConfig {
    pub fastp: FastpConfig,
    pub host_removal: HostRemovalConfig,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct FastpConfig {
    pub threads: usize,
    pub detect_adapters: bool,
    pub cut_front: bool,
    pub cut_tail: bool,
    pub qualified_quality_phred: u8,
    pub length_required: u32,
    pub compression_level: u8,
    /// Skip fastp QC and feed raw reads directly to the host-removal backend.
    #[serde(default)]
    pub skip_qc: bool,
}

impl Default for FastpConfig {
    fn default() -> Self {
        Self {
            threads: 4,
            detect_adapters: true,
            cut_front: true,
            cut_tail: true,
            qualified_quality_phred: 20,
            length_required: 50,
            compression_level: 6,
            skip_qc: false,
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(tag = "mode")]
pub enum HostRemovalConfig {
    #[serde(rename = "kraken2")]
    Kraken2 {
        db_path: PathBuf,
        threads: usize,
        confidence_threshold: f64,
        minimum_hit_groups: u32,
        #[serde(default = "default_memory_mapping")]
        memory_mapping: bool,
        /// Re-align Kraken2-unclassified reads with Bowtie2 against the host index.
        #[serde(default)]
        bowtie2_recheck: bool,
        /// Bowtie2 index prefix used for the optional recheck step.
        bowtie2_index_prefix: Option<PathBuf>,
    },
    #[serde(rename = "minimap2")]
    Minimap2 {
        index_path: PathBuf,
        threads: usize,
    },
    #[serde(rename = "bowtie2")]
    Bowtie2 {
        index_prefix: PathBuf,
        threads: usize,
    },
    #[serde(rename = "sylph")]
    Sylph {
        db_path: PathBuf,
        threads: usize,
    },
    #[serde(rename = "centrifuge")]
    Centrifuge {
        db_path: PathBuf,
        threads: usize,
    },
    #[serde(rename = "auto")]
    Auto {
        kraken2_db_path: PathBuf,
        bowtie2_index_prefix: PathBuf,
        threads: usize,
        // Thresholds for backend selection
        host_pct_low_threshold: f64,
        host_pct_high_threshold: f64,
        reads_high_threshold: u64,
        // User-provided host percentage (0-100). When Some, survey is skipped.
        user_host_pct: Option<f64>,
        // Survey configuration
        survey: bool,
        survey_n_reads: u64,
        survey_threads: usize,
        /// Use Kraken2 --memory-mapping when the auto backend resolves to kraken2.
        #[serde(default = "default_memory_mapping")]
        memory_mapping: bool,
        /// Re-align Kraken2-unclassified reads with Bowtie2 against the host index
        /// when the auto backend resolves to kraken2.
        #[serde(default)]
        bowtie2_recheck: bool,
    },
}

impl HostRemovalConfig {
    pub fn mode(&self) -> &'static str {
        match self {
            HostRemovalConfig::Kraken2 { .. } => "kraken2",
            HostRemovalConfig::Minimap2 { .. } => "minimap2",
            HostRemovalConfig::Bowtie2 { .. } => "bowtie2",
            HostRemovalConfig::Sylph { .. } => "sylph",
            HostRemovalConfig::Centrifuge { .. } => "centrifuge",
            HostRemovalConfig::Auto { .. } => "auto",
        }
    }

    pub fn resolved_mode(&self) -> &'static str {
        match self {
            HostRemovalConfig::Kraken2 { .. } => "kraken2",
            HostRemovalConfig::Minimap2 { .. } => "minimap2",
            HostRemovalConfig::Bowtie2 { .. } => "bowtie2",
            HostRemovalConfig::Sylph { .. } => "sylph",
            HostRemovalConfig::Centrifuge { .. } => "centrifuge",
            HostRemovalConfig::Auto { .. } => "auto-unresolved",
        }
    }

    pub fn threads(&self) -> usize {
        match self {
            HostRemovalConfig::Kraken2 { threads, .. } => *threads,
            HostRemovalConfig::Minimap2 { threads, .. } => *threads,
            HostRemovalConfig::Bowtie2 { threads, .. } => *threads,
            HostRemovalConfig::Sylph { threads, .. } => *threads,
            HostRemovalConfig::Centrifuge { threads, .. } => *threads,
            HostRemovalConfig::Auto { threads, .. } => *threads,
        }
    }
}

fn default_memory_mapping() -> bool {
    true
}

impl Default for HostRemovalConfig {
    fn default() -> Self {
        HostRemovalConfig::Kraken2 {
            db_path: PathBuf::from("/db/minikraken2_v2_8GB"),
            threads: 4,
            confidence_threshold: 0.0,
            minimum_hit_groups: 2,
            memory_mapping: false,
            bowtie2_recheck: false,
            bowtie2_index_prefix: None,
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ValidationConfig {
    pub min_output_size_bytes: u64,
    pub max_contamination_percent: f64,
}

impl Default for ValidationConfig {
    fn default() -> Self {
        Self {
            min_output_size_bytes: 1024,
            max_contamination_percent: 100.0,
        }
    }
}
