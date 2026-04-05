use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "rustyclean")]
#[command(about = "High-performance metagenome QC and host removal pipeline using fastp + kraken2")]
#[command(version, arg_required_else_help = true)]
pub struct Cli {
    /// Forward reads (R1) fastq(.gz) file
    #[arg(long, group = "input")]
    pub r1: Option<PathBuf>,

    /// Reverse reads (R2) fastq(.gz) file (optional, for paired-end)
    #[arg(long, requires = "r1")]
    pub r2: Option<PathBuf>,

    /// Sample list file (TSV: sample_id, r1_path[, r2_path])
    #[arg(short, long, group = "input")]
    pub samples: Option<PathBuf>,

    /// Output directory
    #[arg(short, long, default_value = "rustyclean_output")]
    pub output: PathBuf,

    /// Kraken2 database path
    #[arg(long)]
    pub kraken2_db: Option<PathBuf>,

    /// Configuration file (TOML)
    #[arg(short, long)]
    pub config: Option<PathBuf>,

    /// Checkpoint directory for resume capability
    #[arg(long, default_value = ".rustyclean_checkpoints")]
    pub checkpoint_dir: PathBuf,

    /// Number of parallel workers
    #[arg(short, long)]
    pub workers: Option<usize>,

    /// Number of threads per tool (fastp/kraken2)
    #[arg(short, long)]
    pub threads: Option<usize>,

    /// Resume from previous checkpoints
    #[arg(long)]
    pub resume: bool,

    /// Clean completed checkpoints after run
    #[arg(long)]
    pub clean: bool,

    /// Dry run: validate inputs without processing
    #[arg(long)]
    pub dry_run: bool,
}
