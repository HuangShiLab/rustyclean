# RustyClean

A high-performance metagenome QC and host removal pipeline written in Rust. RustyClean chains [fastp](https://github.com/OpenGene/fastp) (quality control & adapter trimming) and [Kraken2](https://github.com/DerrickWood/kraken2) (taxonomic classification & host removal) into a streamlined workflow for metagenomic sequencing data.

## Features

- **Dual input mode** -- direct FASTQ(.gz) file input or batch processing via sample list
- **Single-end & paired-end** -- automatically adapts the pipeline based on input
- **Parallel processing** -- concurrent sample processing with configurable worker count
- **Checkpoint & resume** -- saves progress automatically; interrupted runs can be resumed
- **Retry on failure** -- configurable retry attempts per sample
- **Validation** -- checks output file size and host contamination rate after processing
- **Graceful shutdown** -- Ctrl-C cleanly cancels in-progress work

## Prerequisites

Install the following tools and ensure they are available in `$PATH`:

- [fastp](https://github.com/OpenGene/fastp) (>= 0.23)
- [Kraken2](https://github.com/DerrickWood/kraken2) (>= 2.1) with a pre-built database

## Installation

```bash
git clone https://github.com/HuangShiLab/rustyclean.git
cd rustyclean
cargo build --release

# The binary is at target/release/rustyclean
```

## Quick Start

```bash
# Single sample, paired-end
rustyclean --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
           --kraken2-db /path/to/kraken2_db \
           -o output/

# Single sample, single-end
rustyclean --r1 sample.fastq.gz \
           --kraken2-db /path/to/kraken2_db \
           -o output/

# Batch mode with sample list
rustyclean --samples list.txt \
           --kraken2-db /path/to/kraken2_db \
           -o output/ -w 4 -t 8
```

## Sample List Format

A tab-separated file with 2 or 3 columns (no header required). Lines starting with `#` are treated as comments.

```
# sample_id    R1_path                      R2_path (optional)
sampleA        /data/sampleA_R1.fastq.gz    /data/sampleA_R2.fastq.gz
sampleB        /data/sampleB_R1.fastq.gz    /data/sampleB_R2.fastq.gz
sampleC        /data/sampleC.fastq.gz
```

- **2 columns** (sample_id + R1): single-end mode
- **3 columns** (sample_id + R1 + R2): paired-end mode

## CLI Options

```
Usage: rustyclean [OPTIONS]

Options:
      --r1 <R1>                  Forward reads (R1) fastq(.gz) file
      --r2 <R2>                  Reverse reads (R2) fastq(.gz) file (paired-end)
  -s, --samples <SAMPLES>        Sample list file (TSV)
  -o, --output <OUTPUT>          Output directory [default: rustyclean_output]
      --kraken2-db <KRAKEN2_DB>  Kraken2 database path
  -c, --config <CONFIG>          Configuration file (TOML)
      --checkpoint-dir <DIR>     Checkpoint directory [default: .rustyclean_checkpoints]
  -w, --workers <WORKERS>        Number of parallel workers
  -t, --threads <THREADS>        Number of threads per tool (fastp/kraken2)
      --resume                   Resume from previous checkpoints
      --clean                    Clean completed checkpoints after run
      --dry-run                  Validate inputs without processing
  -h, --help                     Print help
  -V, --version                  Print version
```
## Pipeline

```
Input FASTQ(.gz)
      |
      v
  +---------+
  |  fastp   |  Quality control, adapter trimming, read filtering
  +---------+
      |
      v
  +---------+
  | kraken2  |  Taxonomic classification, host (human) read removal
  +---------+
      |
      v
  Validation    Check output size & contamination rate
      |
      v
  Clean FASTQ   {sample_id}_clean_R1.fastq.gz [+ _R2.fastq.gz]
```

## Output

For each sample, the final output is written to the output directory:

```
output/
  sampleA/
    sampleA_clean_R1.fastq.gz
    sampleA_clean_R2.fastq.gz   # paired-end only
  sampleB/
    sampleB_clean_R1.fastq.gz
```

## License

MIT
