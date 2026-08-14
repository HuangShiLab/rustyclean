#!/usr/bin/env python3
"""Compare bowtie2-recheck results with original rustyclean results."""
import csv
import statistics
from pathlib import Path

DATASETS = [
    "30M_50pct_high_skewed_SE",
    "60M_90pct_high_lognormal_SE",
    "100M_50pct_high_lognormal_SE",
    "100M_90pct_high_lognormal_SE",
]


def parse_runtime(s):
    s = s.strip()
    if s == "unknown":
        return float("nan")
    parts = s.split(":")
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    elif len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    else:
        return float(parts[0])


def summarize(path, tool_name, runtime_col="runtime_seconds"):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r["tool"] == tool_name and r["dataset"] in DATASETS:
                rows.append(r)
    summary = {}
    for d in DATASETS:
        vals = [parse_runtime(r[runtime_col]) for r in rows if r["dataset"] == d]
        mems = [float(r["max_memory_kb"]) for r in rows if r["dataset"] == d]
        if vals:
            summary[d] = {
                "mean_s": statistics.mean(vals),
                "sd_s": statistics.stdev(vals) if len(vals) > 1 else 0,
                "mean_mem_kb": statistics.mean(mems),
            }
    return summary


def main():
    orig = summarize(
        "/scr/u/shihuang/rustyclean-paper/results_v2/metrics/performance.csv",
        "rustyclean",
    )
    bt2 = summarize(
        "/lustre1/g/aos_shihuang/rustyclean-paper/results_v2/metrics/performance_bowtie2_recheck.csv",
        "rustyclean_bt2recheck",
    )

    print("Dataset | Original RC (min) | BT2-recheck (min) | Slowdown | Original Mem (GB) | BT2 Mem (GB)")
    print("---|---|---|---|---|---")
    for d in DATASETS:
        o = orig.get(d, {})
        b = bt2.get(d, {})
        if not o or not b:
            print(f"{d} | missing | missing | - | - | -")
            continue
        o_min = o["mean_s"] / 60
        b_min = b["mean_s"] / 60
        slowdown = b_min / o_min if o_min else float("nan")
        print(
            f"{d} | {o_min:.1f} ± {o['sd_s']/60:.1f} | {b_min:.1f} ± {b['sd_s']/60:.1f} | "
            f"{slowdown:.2f}x | {o['mean_mem_kb']/1e6:.1f} | {b['mean_mem_kb']/1e6:.1f}"
        )


if __name__ == "__main__":
    main()
