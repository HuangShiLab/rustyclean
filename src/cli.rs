use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "rustyclean")]
#[command(about = "High-performance metagenome QC and host removal pipeline using fastp + kraken2/minimap2/bowtie2/kraken2+bowtie2-recheck/centrifuge/auto")]
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

    /// Host-removal backend. `auto` selects bowtie2 for low-host samples and
    /// kraken2 with bowtie2 recheck for high-host samples based on user-provided host
    /// percentage or a light-weight survey.
    #[arg(long, value_enum, visible_alias = "mode", default_value = "kraken2")]
    pub host_removal_mode: HostRemovalModeCli,

    /// Expected host contamination percentage (0-100). When provided with
    /// `--host-removal-mode auto`, it bypasses the survey and directly
    /// selects the backend.
    #[arg(long, value_name = "PCT")]
    pub host_pct: Option<f64>,

    /// Enable light-weight survey for `--host-removal-mode auto`.
    /// The first `--auto-survey-nreads` reads are mapped with bowtie2 to
    /// estimate host percentage when `--host-pct` is not given.
    #[arg(long)]
    pub auto_survey: bool,

    /// Number of reads to survey in auto mode (default: 100k).
    #[arg(long, default_value = "100000")]
    pub auto_survey_nreads: u64,

    /// Threads used for the auto survey (default: 2).
    #[arg(long, default_value = "2")]
    pub auto_survey_threads: usize,

    /// Low-host threshold (%). Below this, auto mode prefers bowtie2.
    #[arg(long, default_value = "10.0")]
    pub auto_low_threshold: f64,

    /// High-host threshold (%). Above this with a large sample, auto mode
    /// prefers kraken2.
    #[arg(long, default_value = "30.0")]
    pub auto_high_threshold: f64,

    /// Large-sample read-count threshold. Combined with high-host threshold
    /// to switch to kraken2.
    #[arg(long, default_value = "20000000")]
    pub auto_reads_threshold: u64,

    /// Kraken2 database path (for kraken2 / auto mode)
    #[arg(long)]
    pub kraken2_db: Option<PathBuf>,

    /// Use Kraken2 --memory-mapping to share the database via mmap
    /// instead of loading it into RAM (default: disabled).
    #[arg(long, default_value_t = false)]
    pub kraken2_memory_mapping: bool,

    /// After Kraken2 classification, re-align the reads Kraken2 called host
    /// against the index given here, and keep the ones Bowtie2 cannot place on
    /// the host. The pass exists to lower the number of microbial reads
    /// discarded by mistake, at the cost of additional runtime.
    ///
    /// Supplying a Bowtie2 index prefix enables the pass; omitting the flag
    /// disables it. The index may differ from --host-index, so the survey and
    /// the verification pass can use different references.
    #[arg(long, value_name = "BOWTIE2_INDEX_PREFIX")]
    pub bowtie2_recheck: Option<PathBuf>,

    /// Path to human reference index:
    /// - minimap2: .mmi file
    /// - bowtie2: index prefix
    /// - sylph: bowtie2 index prefix used for the actual removal step
    /// - centrifuge: index prefix (minus trailing .X.cf)
    /// - auto: bowtie2 index prefix used for the survey and low-host branch
    #[arg(long)]
    pub host_index: Option<PathBuf>,

    /// Sylph database path (.syldb) for the sylph backend and for the
    /// high-host branch of auto mode.
    #[arg(long)]
    pub sylph_db: Option<PathBuf>,

    /// Minimum Adjusted_ANI (%) reported by sylph for a sample to be considered
    /// host-positive and passed to Bowtie2 removal (default: 95.0).
    #[arg(long, default_value_t = 95.0)]
    pub sylph_min_ani: f64,

    /// Minimum effective coverage reported by sylph for a sample to be considered
    /// host-positive and passed to Bowtie2 removal (default: 0.0005).
    #[arg(long, default_value_t = 0.0005)]
    pub sylph_min_cov: f64,

    /// Maximum allowed host contamination percent in output (default: 100.0)
    #[arg(long, default_value_t = 100.0)]
    pub max_contamination: f64,

    /// Configuration file (TOML)
    #[arg(short, long)]
    pub config: Option<PathBuf>,

    /// Checkpoint directory for resume capability
    #[arg(long, default_value = ".rustyclean_checkpoints")]
    pub checkpoint_dir: PathBuf,

    /// Number of parallel workers
    #[arg(short, long)]
    pub workers: Option<usize>,

    /// Number of threads per tool (fastp/host-removal)
    #[arg(short, long)]
    pub threads: Option<usize>,

    /// Skip fastp QC and feed raw reads directly to the host-removal backend.
    /// Useful when reads have already been quality-controlled, or when comparing
    /// only the host-removal step against tools that do not perform QC.
    #[arg(long, default_value_t = false)]
    pub skip_qc: bool,

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

#[derive(Debug, Clone, Copy, Default, clap::ValueEnum)]
pub enum HostRemovalModeCli {
    #[default]
    Kraken2,
    Minimap2,
    Bowtie2,
    Sylph,
    Centrifuge,
    Auto,
}



#[cfg(test)]
mod recheck_flag_tests {
    use super::*;
    use clap::Parser;

    fn cli(extra: &[&str]) -> Cli {
        let mut argv = vec!["rustyclean", "--r1", "x.fq.gz", "--kraken2-db", "/db"];
        argv.extend_from_slice(extra);
        Cli::try_parse_from(argv).expect("args should parse")
    }

    #[test]
    fn absent_flag_disables_the_pass() {
        assert!(cli(&[]).bowtie2_recheck.is_none());
    }

    #[test]
    fn supplying_an_index_enables_the_pass_and_carries_it() {
        let c = cli(&["--bowtie2-recheck", "/db/human_t2t_hla"]);
        assert_eq!(c.bowtie2_recheck.as_deref(), Some(std::path::Path::new("/db/human_t2t_hla")));
    }

    #[test]
    fn the_flag_requires_a_value() {
        assert!(Cli::try_parse_from(
            ["rustyclean", "--r1", "x.fq.gz", "--kraken2-db", "/db", "--bowtie2-recheck"]
        ).is_err());
    }

    #[test]
    fn recheck_index_is_independent_of_host_index() {
        let c = cli(&["--host-index", "/idx/survey", "--bowtie2-recheck", "/idx/verify"]);
        assert_eq!(c.host_index.as_deref(), Some(std::path::Path::new("/idx/survey")));
        assert_eq!(c.bowtie2_recheck.as_deref(), Some(std::path::Path::new("/idx/verify")));
    }
}
