#!/usr/bin/env python3
"""Re-parse rustyclean bowtie2-recheck time.log files into a metrics CSV."""
import csv
import os
from datetime import datetime
from pathlib import Path

BASE = Path("/scr/u/shihuang/rustyclean-paper/bowtie2_recheck_results")
OUT = Path("/lustre1/g/aos_shihuang/rustyclean-paper/results_v2/metrics/performance_bowtie2_recheck.csv")

OUT.parent.mkdir(parents=True, exist_ok=True)


def parse_time(t: str) -> float:
    t = t.strip()
    parts = t.split(":")
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    elif len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    else:
        return float(parts[0])


rows = []
for d in sorted(BASE.iterdir()):
    if not d.is_dir():
        continue
    dataset = d.name
    for r in sorted(d.iterdir()):
        if not r.is_dir() or not r.name.startswith("rep_"):
            continue
        rep = r.name.replace("rep_", "")
        tl = r / "time.log"
        if not tl.exists() or tl.stat().st_size == 0:
            continue
        text = tl.read_text()
        rt = "0:0"
        mem = "0"
        for line in text.splitlines():
            if "Elapsed (wall clock) time" in line:
                rt = line.split(": ", 1)[1].strip()
            if "Maximum resident set size (kbytes)" in line:
                mem = line.split(":", 1)[1].strip()
                mem = mem.split()[0].strip()
        try:
            rt_sec = f"{parse_time(rt):.2f}"
        except Exception:
            rt_sec = rt
        ts = datetime.fromtimestamp(tl.stat().st_mtime).astimezone().isoformat()
        rows.append({
            "tool": "rustyclean_bt2recheck",
            "dataset": dataset,
            "rep": rep,
            "runtime_seconds": rt_sec,
            "max_memory_kb": mem,
            "timestamp": ts,
        })

with open(OUT, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=["tool", "dataset", "rep", "runtime_seconds", "max_memory_kb", "timestamp"])
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} rows to {OUT}")
for row in rows:
    print(row)
