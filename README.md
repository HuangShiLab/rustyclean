# RustyClean

A high-performance metagenome QC and host removal pipeline written in Rust. RustyClean chains [fastp](https://github.com/OpenGene/fastp) (quality control & adapter trimming) with a choice of host-removal backends into a streamlined workflow for metagenomic sequencing data.

## Features

- **Multiple host-removal backends** -- Kraken2, minimap2, Bowtie2, sylph, Centrifuge, deacon, or adaptive `auto` mode
- **AUTO mode** -- with `--deacon-index`, every sample runs the deacon minimizer backend as Tier-1 (per-read cost independent of host fraction) with an automatic Bowtie2 recheck for high-host samples; without a deacon index it surveys a small subset of reads and selects Bowtie2 (low-host) or sylph+bowtie2 (high-host samples)
- **`--skip-qc` mode** -- bypass fastp and feed raw reads directly to the host-removal backend, useful for already-QC'd data or fair benchmarking of host removal only
- **Dual input mode** -- direct FASTQ(.gz) file input or batch processing via sample list
- **Single-end & paired-end** -- automatically adapts the pipeline based on input
- **Parallel processing** -- concurrent sample processing with configurable worker count; when workers are not specified, RustyClean caps concurrency by available memory to avoid loading more database copies than fit in RAM
- **Checkpoint & resume** -- saves progress automatically; interrupted runs can be resumed
- **Retry on failure** -- configurable retry attempts per sample
- **Validation** -- checks output file size and host contamination rate after processing
- **Graceful shutdown** -- Ctrl-C cleanly cancels in-progress work

## Host-removal backends

| Backend | Description | Best for |
|---------|-------------|----------|
| `kraken2` | k-mer based taxonomic classification against a host-specific Kraken2 database; removes reads classified as *Homo sapiens* (taxid 9606). Default human database is T2T-only. | Large, high-host samples; when microbial context is also useful |
| `bowtie2` | Short-read alignment against a host reference index | Low-host samples; fastest when host fraction is small |
| `sylph` | Fast k-mer sketch prefilter (`sylph query`) followed by Bowtie2 read-level removal for host-positive samples | Very fast screening of large cohorts; only runs full alignment when host signal is detected |
| `minimap2` | Long- or short-read alignment (`-x sr`) | Long reads or when a minimap2 index is preferred |
| `centrifuge` | Compressed FM-index taxonomic classification | Alternative k-mer classifier; removes human taxid 9606 reads |
| `deacon` | Minimizer-based depletion with the external `deacon filter -d` tool against an index built by `deacon index build`; reads meeting the `--deacon-abs-threshold` / `--deacon-rel-threshold` minimizer-hit thresholds are discarded | Very fast, memory-light host removal when a deacon index is available; default Tier-1 backend in `auto` mode when `--deacon-index` is given |
| `auto` | With `--deacon-index`: deacon as Tier-1 for every sample, plus a Bowtie2 recheck of deacon-retained reads when the removed proportion reaches `--recheck-threshold` (default 0.3). Without `--deacon-index`: surveys reads with Bowtie2, estimates host %, then picks `bowtie2` or `sylph` | General use; deacon Tier-1 when an index is available, otherwise legacy routing |

## Databases

RustyClean is not restricted to a single host reference. You can point any backend to a custom database or index built from the host genome of interest. Common use cases include human (e.g. GRCh38, T2T-CHM13), mouse, rat, pig, rice, monkey, and other plant or animal host genomes.

**Default human database.** For human metagenomes, the recommended default is a **T2T-only Kraken2 database** built from the T2T-CHM13v2.0 assembly. This database is smaller and faster than mixed multi-host libraries while retaining high accuracy for human host removal. A mixed multi-host Kraken2 database ("Kraken16", containing human plus other common hosts) is available as an optional taxonomy-aware mode when you expect cross-species contamination or want taxonomic context. GRCh38.p14-based Kraken2 databases are also supported.

**Deacon human index (panhuman-1).** For deacon-based removal (explicit `--host-removal-mode deacon`, or `auto` mode with `--deacon-index`), the recommended human index is the prebuilt **panhuman-1** (k31w15 minimizers): a union of human pangenome references minus FDA-ARGOS bacteria and RefSeq viral sequences, so microbial reads are not depleted along with host reads. The index is user-provided via `--deacon-index`; for non-human hosts, build a custom index with `deacon index build`. Alongside it, the Kraken2 databases described above (T2T-only, GRCh38, or mixed multi-host) remain supported for the kraken2 backend and the legacy auto-routing fallback.

| Backend | Database / index type | How to specify |
|---------|----------------------|----------------|
| `kraken2` | Pre-built Kraken2 database containing the host taxon | `--kraken2-db /path/to/kraken2_db` |
| `bowtie2` | Bowtie2 index prefix (files `{prefix}.1.bt2`, `{prefix}.2.bt2`, ...) | `--host-index /path/to/bowtie2_index_prefix` |
| `sylph` | sylph sketch database (`.syldb`) plus a Bowtie2 index prefix for full removal | `--sylph-db /path/to/human.syldb` and `--host-index /path/to/bowtie2_index_prefix` |
| `minimap2` | Minimap2 index file (`.mmi`) | `--host-index /path/to/index.mmi` |
| `centrifuge` | Centrifuge index prefix (files `{prefix}.1.cf`, `{prefix}.2.cf`, ...) | `--host-index /path/to/centrifuge_index_prefix` |
| `deacon` | Deacon minimizer index built by `deacon index build` (e.g. the prebuilt human panhuman-1 index, k31w15) | `--deacon-index /path/to/deacon_index` |
| `auto` | Deacon index for Tier-1 (with a Bowtie2 index prefix for the recheck), or a sylph database + Bowtie2 index prefix for legacy routing | `--deacon-index ...` and `--host-index ...`, or `--sylph-db ...` and `--host-index ...` |

For `kraken2` and `centrifuge`, the database only needs to contain the host lineage (e.g. *Homo sapiens*, taxid 9606). For `bowtie2` and `minimap2`, build the index directly from the host reference FASTA. This makes RustyClean applicable across diverse host species and metagenome types (saliva, vaginal, gut, plant root, etc.).

### Reference genomes used in benchmarks

The following NCBI RefSeq assemblies were used to build the host indices in our benchmark studies. You can use the same references or substitute your own host genome of interest.

**Note:** For human host removal, the default reference is **T2T-CHM13v2.0**. GRCh38.p14 is provided as an alternative standard reference. They are not combined; choose one according to your study.

| Host | Assembly | NCBI accessions / download links |
|------|----------|----------------------------------|
| Human (GRCh38) | GRCh38.p14 | `GCF_000001405.40_GRCh38.p14_genomic.fna.gz` |
| Human (T2T+HLA) | T2T-CHM13v2.0 + HLA | Hostile-prepared `human-t2t-hla.fa.gz` (T2T-CHM13v2.0 plus HLA contigs) |
| Mouse | GRCm39 | `GCF_000001635.27_GRCm39_genomic.fna.gz` |
| Rat | mRatBN7.2 | `GCF_015227675.2_mRatBN7.2_genomic.fna.gz` |
| Pig | Sscrofa11.1 | `GCF_000003025.6_Sscrofa11.1_genomic.fna.gz` |
| Rice | IRGSP-1.0 | `GCF_001433935.1_IRGSP-1.0_genomic.fna.gz` |
| Monkey (rhesus) | Mmul_10 | `GCF_003339765.1_Mmul_10_genomic.fna.gz` |

#### Pre-built indices on HKU HPC2021

If you have access to the HKU HPC2021 cluster, the following pre-built indices are available under `/lustre1/g/aos_shihuang/databases/`:

| Index | Path | Backend | Note |
|-------|------|---------|------|
| Human GRCh38.p14 FASTA | `/lustre1/g/aos_shihuang/databases/human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz` | Bowtie2 / minimap2 | Build your own index with `bowtie2-build` or `minimap2 -d` |
| Human GRCh38.p14 minimap2 | `/lustre1/g/aos_shihuang/databases/human/GRCh38.p14.mmi` | minimap2 | Ready to use |
| Human hg39 Bowtie2 | `/lustre1/g/aos_shihuang/databases/kneaddata/hg_39` | Bowtie2 | KneadData-compatible human index |
| Human T2T+HLA | `/home/shihuang/.local/share/hostile/human-t2t-hla` | Bowtie2 / minimap2 | T2T-CHM13v2.0 + HLA sequences prepared by Hostile |
| Human T2T+HLA (copy, incomplete) | `/lustre1/g/aos_shihuang/databases/rustyclean_alt/human_t2t_hla` | Bowtie2 / minimap2 | Index build incomplete; use the Hostile path above |
| Cross-species multi-host Bowtie2 | `/lustre1/g/aos_shihuang/databases/host_genomes_cross/multi_host_bt2` | Bowtie2 | Human + mouse + rat + pig + rice + monkey combined |
| Human T2T-only Kraken2 | `/lustre1/g/aos_shihuang/databases/rustyclean_human_t2t_only/kraken2/t2t_only` | Kraken2 | **Default human database** (T2T-CHM13v2.0 only) |
| Human T2T-only sylph | `/lustre1/g/aos_shihuang/databases/rustyclean_human_t2t_only/sylph/human_t2t.syldb` | sylph | Fast k-mer sketch of T2T-CHM13v2.0; used with a Bowtie2 index for full read-level removal |
| Mixed multi-host Kraken2 ("Kraken16") | `/lustre1/g/aos_shihuang/databases/kraken2/kraken16` | Kraken2 | Optional taxonomy-aware mode; contains human lineage (taxid 9606) plus microbial genomes |

For Kraken2/Centrifuge, only the host lineage (e.g. taxid 9606 for human) needs to be present in the database. For Bowtie2 and minimap2, build the index directly from the reference FASTA:

```bash
# Bowtie2 index
bowtie2-build host.fa host_index_prefix

# Minimap2 index
minimap2 -x sr -d host_index.mmi host.fa
```

#### Building the default human T2T-only Kraken2 database

The default human Kraken2 database is built from the T2T-CHM13v2.0 assembly only. Example:

```bash
DB_DIR="/path/to/rustyclean_human_t2t_only"
mkdir -p "${DB_DIR}/kraken2/t2t_only"

# Download T2T-CHM13v2.0
wget -P "${DB_DIR}" \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz

# Build Kraken2 database
kraken2-build --download-taxonomy --db "${DB_DIR}/kraken2/t2t_only" --threads 8
kraken2-build --add-to-library \
  "${DB_DIR}/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" \
  --db "${DB_DIR}/kraken2/t2t_only" --threads 8
kraken2-build --build --db "${DB_DIR}/kraken2/t2t_only" --threads 8

# Use it
rustyclean --r1 sample.fastq.gz --kraken2-db "${DB_DIR}/kraken2/t2t_only" -o output/
```

## Prerequisites

Install the following tools and ensure they are available in `$PATH`. Only the tools required by your chosen backend need to be installed.

- [fastp](https://github.com/OpenGene/fastp) (>= 0.23) -- required unless `--skip-qc` is used
- [Kraken2](https://github.com/DerrickWood/kraken2) (>= 2.1) with a pre-built database -- for `--host-removal-mode kraken2` or `auto`
- [Bowtie2](https://github.com/BenLangmead/bowtie2) (>= 2.4) and [samtools](http://www.htslib.org/) -- for `--host-removal-mode bowtie2`, `--host-removal-mode sylph`, or `auto`
- [minimap2](https://github.com/lh3/minimap2) -- for `--host-removal-mode minimap2`
- [Centrifuge](https://ccb.jhu.edu/software/centrifuge/) -- for `--host-removal-mode centrifuge`
- [sylph](https://github.com/bluenote-1577/sylph) (>= 0.9) -- for `--host-removal-mode sylph`
- [deacon](https://github.com/dnbaker/deacon) (>= 0.17; `cargo install deacon`) -- for `--host-removal-mode deacon`, or as the default Tier-1 backend in `auto` mode when `--deacon-index` is given. The prebuilt **panhuman-1** index (k31w15) is the recommended human index; see [Databases](#databases)

## Installation

```bash
git clone https://github.com/HuangShiLab/rustyclean.git
cd rustyclean
cargo build --release

# The binary is at target/release/rustyclean
```

## Quick Start

```bash
# Single sample, paired-end, kraken2 mode
rustyclean --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
           --kraken2-db /path/to/kraken2_db \
           -o output/

# Single sample, single-end, kraken2 mode
rustyclean --r1 sample.fastq.gz \
           --kraken2-db /path/to/kraken2_db \
           -o output/

# AUTO mode: adaptively choose bowtie2/kraken2 based on a 100k-read survey
rustyclean --r1 sample.fastq.gz \
           --host-removal-mode auto \
           --kraken2-db /path/to/kraken2_db \
           --host-index /path/to/bowtie2_index_prefix \
           --auto-survey \
           -o output/ -t 8

# sylph prefilter + Bowtie2 removal (very fast for cohort screening)
rustyclean --r1 sample.fastq.gz \
           --host-removal-mode sylph \
           --sylph-db /path/to/human_t2t.syldb \
           --host-index /path/to/bowtie2_index_prefix \
           -o output/ -t 8

# Skip QC and run host removal only (fair comparison with tools like Hostile)
rustyclean --r1 sample.fastq.gz \
           --host-removal-mode bowtie2 \
           --host-index /path/to/bowtie2_index_prefix \
           --skip-qc \
           -o output/ -t 8

# Batch mode with sample list
rustyclean --samples list.txt \
           --host-removal-mode auto \
           --kraken2-db /path/to/kraken2_db \
           --host-index /path/to/bowtie2_index_prefix \
           --auto-survey \
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
      --r1 <R1>                        Forward reads (R1) fastq(.gz) file
      --r2 <R2>                        Reverse reads (R2) fastq(.gz) file (paired-end)
  -s, --samples <SAMPLES>              Sample list file (TSV)
  -o, --output <OUTPUT>                Output directory [default: rustyclean_output]
      --host-removal-mode <MODE>       Host-removal backend: kraken2, minimap2, bowtie2,
                                       centrifuge, auto [default: auto]
      --host-pct <PCT>                 Expected host contamination % (0-100); used by auto mode
      --auto-survey                    Enable lightweight survey for auto mode
      --auto-survey-nreads <N>         Reads to survey [default: 100000]
      --auto-survey-threads <N>        Threads for auto survey [default: 2]
      --auto-low-threshold <PCT>       Low-host threshold for auto mode [default: 10.0]
      --auto-high-threshold <PCT>      High-host threshold for auto mode [default: 30.0]
      --auto-reads-threshold <N>       Large-sample read threshold for auto mode [default: 20000000]
      --kraken2-db <KRAKEN2_DB>        Kraken2 database path
      --kraken2-memory-mapping         Use Kraken2 --memory-mapping
      --sylph-db <PATH>                sylph sketch database (.syldb) for sylph backend
      --sylph-min-ani <PCT>            Minimum Adjusted_ANI (%) for sylph host-positive call [default: 95.0]
      --sylph-min-cov <FLOAT>          Minimum effective coverage for sylph host-positive call [default: 0.0005]
      --host-index <PATH>              Host index path (minimap2 .mmi, bowtie2 prefix,
                                       centrifuge prefix, sylph .syldb, or auto survey index)
      --deacon-index <PATH>            Deacon minimizer index; required for `--host-removal-mode
                                       deacon`, and enables deacon as Tier-1 backend in auto mode
      --recheck-threshold <FLOAT>      Removed-proportion (0-1) triggering the auto-mode
                                       bowtie2 recheck of deacon-retained reads [default: 0.3]
      --max-contamination <PCT>        Max allowed host contamination in output [default: 100.0]
      --skip-qc                        Skip fastp QC and use raw reads for host removal
  -c, --config <CONFIG>                Configuration file (TOML)
      --checkpoint-dir <DIR>           Checkpoint directory [default: .rustyclean_checkpoints]
  -w, --workers <WORKERS>              Number of parallel workers
  -t, --threads <THREADS>              Number of threads per tool
      --resume                         Resume from previous checkpoints
      --clean                          Clean completed checkpoints after run
      --dry-run                        Validate inputs without processing
  -h, --help                           Print help
  -V, --version                        Print version
```

### Memory-aware worker cap

If `-w/--workers` is not set, RustyClean estimates the resident database size
(for example, `hash.k2d` for Kraken2, the Bowtie2 index files, or the minimap2
`.mmi`) and the available memory (cgroup limit first, then `/proc/meminfo`
`MemAvailable`). It then caps the default CPU-based worker count so that
concurrent workers do not collectively exceed ~80% of available RAM. This
prevents out-of-memory failures when many samples are processed in parallel on
shared-memory nodes. You can override the cap by explicitly setting `-w`.

## Pipeline

### Default (QC + host removal)

```
Input FASTQ(.gz)
      |
      v
  +---------+
  |  fastp  |  Quality control, adapter trimming, read filtering
  +---------+
      |
      v
  +------------------+
  | host-removal     |  kraken2 / bowtie2 / sylph / minimap2 / centrifuge / auto
  | (selected backend)
  +------------------+
      |
      v
  Validation           Check output size & contamination rate
      |
      v
  Clean FASTQ          {sample_id}_clean_R1.fastq.gz [+ _R2.fastq.gz]
```

### `--skip-qc` (host removal only)

```
Input FASTQ(.gz)
      |
      v
  +------------------+
  | host-removal     |  kraken2 / bowtie2 / sylph / minimap2 / centrifuge / auto
  | (selected backend)
  +------------------+
      |
      v
  Validation
      |
      v
  Clean FASTQ
```

## Output

### Final clean FASTQ

For each sample, the final output is written to the output directory:

```
output/
  sampleA/
    sampleA_clean_R1.fastq.gz
    sampleA_clean_R2.fastq.gz   # paired-end only
  sampleB/
    sampleB_clean_R1.fastq.gz
```

### Checkpoint directory

By default RustyClean writes per-sample checkpoints and intermediate files under `.rustyclean_checkpoints/`:

```
.rustyclean_checkpoints/
  {sample_id}.json           # Resume checkpoint (pipeline stage, metrics, input hash)
  work/
    {sample_id}/
      fastp.json             # fastp QC report (absent when --skip-qc is used)
      trimmed_R1.fastq.gz    # fastp-trimmed reads (absent when --skip-qc is used)
      trimmed_R2.fastq.gz    # paired-end only
      kraken2.report         # Kraken2 report (kraken2/auto-kraken2 only)
      kraken2.output.txt     # Per-read Kraken2 classifications
      ...                    # Backend-specific intermediate files
```

**Checkpoint JSON** (`{sample_id}.json`) contains:
- Current pipeline stage (`Pending`, `FastpRunning`, `FastpComplete`, `Kraken2Running`, ...)
- Input file hash (for resume consistency)
- `fastp_metrics`: read counts, Q20/Q30, GC content, adapter trimming stats (when QC ran)
- `kraken2_metrics`: host reads kept, unclassified/kept reads, contamination %, output paths
- `auto_backend`: backend chosen by `auto` mode (e.g. `bowtie2` or `kraken2`)
- `validation_result`: pass/fail status, output file size, errors

### Log files

RustyClean prints structured logs to stderr. Redirecting stderr to a file gives a log like:

```
2024-01-15T08:30:12Z  INFO rustyclean: Loaded 2 sample(s), host-removal mode: auto
2024-01-15T08:30:13Z  INFO rustyclean::pipeline: auto mode: selected backend sample=sampleA host_pct="45.20" input_reads=30000000 chosen_backend="kraken2"
2024-01-15T08:35:45Z  INFO rustyclean::pipeline: Sample validated and finalized sample=sampleA unclassified_reads=16500000 human_reads=13500000 contamination="45.00%" auto_backend=Some("kraken2")
```

Key log fields:
- `host-removal mode`: the backend requested on the CLI
- `auto mode: selected backend`: backend chosen by `auto` mode, with estimated host % and input read count
- `Sample validated and finalized`: final kept reads, removed host reads, and contamination rate
- `auto_backend`: confirms which backend was actually used

### fastp JSON

When QC is enabled, `fastp.json` inside the checkpoint work directory is the standard fastp report. It contains pre/post-filtering read counts, quality metrics, adapter trimming statistics, and filtering results.

## License

MIT
