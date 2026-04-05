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
    pub kraken2: Kraken2Config,
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
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Kraken2Config {
    pub db_path: PathBuf,
    pub threads: usize,
    pub confidence_threshold: f64,
    pub minimum_hit_groups: u32,
}

impl Default for Kraken2Config {
    fn default() -> Self {
        Self {
            db_path: PathBuf::from("/db/minikraken2_v2_8GB"),
            threads: 4,
            confidence_threshold: 0.1,
            minimum_hit_groups: 2,
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
            max_contamination_percent: 5.0,
        }
    }
}
