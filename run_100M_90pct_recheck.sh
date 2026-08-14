#!/bin/bash
#SBATCH --job-name=rc_100M90_recheck
#SBATCH --partition=amd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=06:00:00
#SBATCH --output=/home/shihuang/rc_100M90_recheck_%j.out
#SBATCH --error=/home/shihuang/rc_100M90_recheck_%j.err

set -e
source ~/.cargo/env 2>/dev/null || true
source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate /lustre1/g/aos_shihuang/rustyclean-paper/.conda_envs/rustyclean-benchmark
export PATH="/group/aos_shihuang/conda/envs/fastp/bin:/group/aos_shihuang/conda/envs/bowtie2/bin:/group/aos_shihuang/conda/envs/kraken2/bin:/lustre1/g/aos_shihuang/tools/samtools/samtools-1.21:${PATH}"

RC=/lustre1/g/aos_shihuang/rustyclean/target/release/rustyclean
DATA=/scr/u/shihuang/rustyclean-paper/data/enhanced
OUT=/scr/u/shihuang/rustyclean-paper/bowtie2_recheck_results/100M_90pct_high_lognormal_SE/rep_1
KRAKEN_DB=/lustre1/g/aos_shihuang/databases/kraken2/kraken16
HOST_INDEX=/lustre1/g/aos_shihuang/databases/kneaddata/hg_39

rm -rf "$OUT"
mkdir -p "$OUT"
log="$OUT/run.log"
time_log="$OUT/time.log"

/usr/bin/time -v "$RC" \
    --r1 "$DATA/100M_90pct_high_lognormal_SE/reads.fastq.gz" \
    --host-removal-mode auto \
    --kraken2-db "$KRAKEN_DB" \
    --host-index "$HOST_INDEX" \
    --auto-survey \
    --bowtie2-recheck \
    --kraken2-memory-mapping \
    --checkpoint-dir "$OUT/.checkpoints" \
    -o "$OUT" \
    -t 8 \
    --clean \
    > "$log" 2> "$time_log"

echo "Done"
