#!/bin/bash
#SBATCH --job-name=rc_bt2recheck
#SBATCH --partition=amd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/home/shihuang/rc_bt2recheck_%j.out
#SBATCH --error=/home/shihuang/rc_bt2recheck_%j.err

set -e
source ~/.cargo/env 2>/dev/null || true
source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate /lustre1/g/aos_shihuang/rustyclean-paper/.conda_envs/rustyclean-benchmark

# Add tool-specific conda envs to PATH (the benchmark env does not include them)
export PATH="/group/aos_shihuang/conda/envs/fastp/bin:/group/aos_shihuang/conda/envs/bowtie2/bin:/group/aos_shihuang/conda/envs/kraken2/bin:/lustre1/g/aos_shihuang/tools/samtools/samtools-1.21:${PATH}"

RC=/lustre1/g/aos_shihuang/rustyclean/target/release/rustyclean
PROJECT=/scr/u/shihuang/rustyclean-paper
DATA=$PROJECT/data/enhanced
OUT=$PROJECT/bowtie2_recheck_results
METRICS=$PROJECT/results_v2/metrics/performance_bowtie2_recheck.csv
KRAKEN_DB=/lustre1/g/aos_shihuang/databases/kraken2/kraken16
HOST_INDEX=/lustre1/g/aos_shihuang/databases/kneaddata/hg_39

mkdir -p "$OUT"

DATASETS=(
    "30M_50pct_high_skewed_SE"
    "60M_90pct_high_lognormal_SE"
    "100M_50pct_high_lognormal_SE"
    "100M_90pct_high_lognormal_SE"
)

# Convert h:mm:ss or m:ss to seconds
parse_time() {
    local t="$1"
    local s=0
    if [[ "$t" == *:* ]]; then
        local n
        n=$(echo "$t" | awk -F: '{print NF}')
        if [ "$n" -eq 2 ]; then
            s=$(echo "$t" | awk -F: '{print ($1*60)+$2}')
        elif [ "$n" -eq 3 ]; then
            s=$(echo "$t" | awk -F: '{print ($1*3600)+($2*60)+$3}')
        fi
    else
        s="$t"
    fi
    printf "%.2f" "$s"
}

# Write header
printf "tool,dataset,rep,runtime_seconds,max_memory_kb,timestamp\n" > "$METRICS"

for dataset in "${DATASETS[@]}"; do
    R1="$DATA/$dataset/reads.fastq.gz"
    if [ ! -f "$R1" ]; then
        echo "WARNING: $R1 not found, skipping $dataset"
        continue
    fi
    for rep in 1 2 3; do
        ds_out="$OUT/$dataset/rep_$rep"
        rm -rf "$ds_out"
        mkdir -p "$ds_out"
        log="$ds_out/run.log"
        time_log="$ds_out/time.log"
        echo "=== $dataset rep $rep ===" | tee -a "$log"
        /usr/bin/time -v "$RC" \
            --r1 "$R1" \
            --host-removal-mode auto \
            --kraken2-db "$KRAKEN_DB" \
            --host-index "$HOST_INDEX" \
            --auto-survey \
            --bowtie2-recheck \
            --kraken2-memory-mapping \
            --checkpoint-dir "$ds_out/.checkpoints" \
            -o "$ds_out" \
            -t 8 \
            --clean \
            > "$log" 2> "$time_log"
        rt=$(grep "Elapsed (wall clock) time:" "$time_log" | sed 's/.*: //' || echo "0:0")
        mem=$(grep "Maximum resident set size (kbytes):" "$time_log" | awk '{print $NF}' || echo "0")
        rt_sec=$(parse_time "$rt")
        ts=$(date -Iseconds)
        printf "rustyclean_bt2recheck,%s,%s,%s,%s,%s\n" "$dataset" "$rep" "$rt_sec" "$mem" "$ts" >> "$METRICS"
        echo "  runtime=${rt_sec}s memory=${mem}kb" | tee -a "$log"
    done
done

echo "Benchmark complete. Metrics: $METRICS"
