#!/usr/bin/env python3
"""
01_prepare_data.py

Summarize Zephyr respiratory FASTA files.

Outputs:
    results/pool_summary.csv
"""

from pathlib import Path
import gzip

import pandas as pd
from Bio import SeqIO

DATA_DIR = Path("data")
RESULT_DIR = Path("results")

RESULT_DIR.mkdir(exist_ok=True)

summary = []


def open_fasta(path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


files = sorted(DATA_DIR.glob("*.respiratory.fasta.gz"))

if not files:
    raise FileNotFoundError("No Zephyr FASTA files found in data/")

print(f"Found {len(files)} files\n")

for fasta in files:

    print(f"Processing {fasta.name}")

    read_lengths = []
    total_gc = 0
    total_bp = 0

    with open_fasta(fasta) as handle:

        for record in SeqIO.parse(handle, "fasta"):

            seq = str(record.seq).upper()

            L = len(seq)

            read_lengths.append(L)

            total_bp += L
            total_gc += seq.count("G")
            total_gc += seq.count("C")

    gc = (100 * total_gc / total_bp) if total_bp else 0

    summary.append({

        "pool": fasta.name.replace(".respiratory.fasta.gz", ""),

        "reads": len(read_lengths),

        "bases": total_bp,

        "mean_length": round(sum(read_lengths) / len(read_lengths), 2),

        "median_length": round(pd.Series(read_lengths).median(), 2),

        "min_length": min(read_lengths),

        "max_length": max(read_lengths),

        "gc_percent": round(gc, 2)

    })

summary = pd.DataFrame(summary)

summary.sort_values("reads", ascending=False, inplace=True)

summary.to_csv(

    RESULT_DIR / "pool_summary.csv",

    index=False

)

print("\nSaved results/pool_summary.csv")
print(summary)