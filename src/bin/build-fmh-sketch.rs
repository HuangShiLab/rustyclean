//! Build a human-specific FracMinHash k-mer sketch for the `fmh` host-removal
//! backend of RustyClean.
//!
//! Two-pass, memory-efficient construction:
//!   1. Pass 1 streams the host genome(s) and builds a FracMinHash set of
//!      canonical k-mers (only hashes below `u64::MAX / scale` are kept, so
//!      the set holds ~genome/scale entries).
//!   2. Pass 2 streams the microbial reference FASTA(s); sampled microbial
//!      k-mers are checked against the host set and any shared hash is
//!      removed. No full microbial k-mer set is ever built, so large
//!      pangenomes stay cheap.

use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::Parser;

use rustyclean::fmh::{passes_threshold, stream_fasta_lines, FmhSketch, KmerRoller};

#[derive(Parser, Debug)]
#[command(
    name = "rustyclean-build-fmh-sketch",
    about = "Build a human-specific FracMinHash k-mer sketch for RustyClean's fmh host-removal backend",
    version
)]
struct Args {
    /// Host genome FASTA file(s), possibly gzipped (.fa[.gz], multi-FASTA ok).
    /// May be given multiple times.
    #[arg(long, required = true)]
    host: Vec<PathBuf>,

    /// Microbial reference FASTA file(s) or directories of .fa/.fasta(.gz)
    /// files used to subtract shared (non-host-specific) k-mers. May be given
    /// multiple times.
    #[arg(long, required = true)]
    microbes: Vec<PathBuf>,

    /// Output sketch path (.fmh).
    #[arg(long, short)]
    output: PathBuf,

    /// K-mer size (1-32).
    #[arg(long, default_value_t = 31)]
    kmer: u32,

    /// FracMinHash scale denominator: keep ~1/scale of all k-mers.
    #[arg(long, default_value_t = 40)]
    scale: u32,
}

fn main() -> Result<()> {
    let args = Args::parse();

    if !(1..=32).contains(&args.kmer) {
        bail!("--kmer must be in 1..=32, got {}", args.kmer);
    }
    if args.scale == 0 {
        bail!("--scale must be >= 1");
    }

    for host in &args.host {
        if !host.is_file() {
            bail!("--host file not found: {}", host.display());
        }
    }
    let microbe_files = expand_microbe_inputs(&args.microbes)?;
    if microbe_files.is_empty() {
        bail!("no microbial FASTA files found in --microbes inputs");
    }

    // Pass 1: host FMH set.
    let mut sketch = FmhSketch::new(args.kmer, args.scale);
    let mut total_host_kmers = 0u64;
    for host in &args.host {
        let mut roller = KmerRoller::new(args.kmer);
        stream_fasta_lines(host, |seq, record_start| {
            if record_start {
                roller.reset();
            }
            for &b in seq {
                if let Some(h) = roller.push(b) {
                    total_host_kmers += 1;
                    sketch.insert_mixed(h);
                }
            }
        })
        .with_context(|| format!("failed to stream host FASTA {}", host.display()))?;
    }

    let host_specific = sketch.hashes.len() as u64;

    // Pass 2: subtract sampled microbial k-mers shared with the host set.
    let mut removed = 0u64;
    for microbe in &microbe_files {
        let mut roller = KmerRoller::new(args.kmer);
        stream_fasta_lines(microbe, |seq, record_start| {
            if record_start {
                roller.reset();
            }
            for &b in seq {
                if let Some(h) = roller.push(b) {
                    if passes_threshold(h, args.scale) && sketch.hashes.remove(&h) {
                        removed += 1;
                    }
                }
            }
        })
        .with_context(|| format!("failed to stream microbial FASTA {}", microbe.display()))?;
    }

    sketch
        .save(&args.output)
        .with_context(|| format!("failed to write sketch {}", args.output.display()))?;

    println!(
        "Built FMH sketch: {} (k={}, scale={}, host_kmers={}, host_specific={}, shared_removed={})",
        args.output.display(),
        args.kmer,
        args.scale,
        total_host_kmers,
        host_specific,
        removed
    );

    Ok(())
}

/// Expand `--microbes` inputs (files and/or directories) into a sorted list
/// of FASTA files.
fn expand_microbe_inputs(inputs: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for input in inputs {
        if input.is_dir() {
            let mut dir_files: Vec<PathBuf> = Vec::new();
            for entry in std::fs::read_dir(input)
                .with_context(|| format!("failed to read directory {}", input.display()))?
            {
                let entry = entry?;
                let path = entry.path();
                if path.is_file() && is_fasta_name(&path) {
                    dir_files.push(path);
                }
            }
            dir_files.sort();
            files.extend(dir_files);
        } else if input.is_file() {
            files.push(input.clone());
        } else {
            bail!("--microbes path not found: {}", input.display());
        }
    }
    files.sort();
    files.dedup();
    Ok(files)
}

fn is_fasta_name(path: &Path) -> bool {
    let name = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let name = name.strip_suffix(".gz").unwrap_or(&name).to_string();
    name.ends_with(".fa") || name.ends_with(".fasta") || name.ends_with(".fna")
}
