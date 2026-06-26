#!/usr/bin/env python3
"""
03_genome_coverage.py

Compute genome coverage statistics from the aligned BAM file.

Inputs:
    results/aligned.bam

Outputs:
    results/genome_coverage.csv
    plots/genome_coverage.png
"""

from pathlib import Path
import subprocess

import pandas as pd
import matplotlib.pyplot as plt

RESULT_DIR = Path("results")
PLOT_DIR = Path("plots")

BAM_FILE = RESULT_DIR / "aligned.bam"

RESULT_DIR.mkdir(exist_ok=True)
PLOT_DIR.mkdir(exist_ok=True)

REF_NAMES = {
    "NC_001617.1": "Rhinovirus-A",
    "NC_009996.1": "Rhinovirus-B",
    "NC_004148.2": "Metapneumovirus",
    "MN908947.3": "SARS-CoV-2",
    "NC_038235.1": "RSV-A",
    "NC_001803.1": "RSV-B",
    "NC_007373.1": "Influenza-A",
    "NC_002023.1": "Influenza-B",
    "NC_003461.1": "Parainfluenza-1",
    "NC_003462.2": "Parainfluenza-2",
    "NC_001796.2": "Parainfluenza-3",
    "NC_001798.2": "Parainfluenza-4",
    "NC_001405.1": "Adenovirus-C",
}

if not BAM_FILE.exists():
    raise FileNotFoundError(f"{BAM_FILE} not found. Run 02_taxonomic_classification.py first.")

print("Reading BAM header...")

header = subprocess.run(
    ["samtools", "view", "-H", str(BAM_FILE)],
    capture_output=True,
    text=True,
    check=True,
).stdout

reference_lengths = {}

for line in header.splitlines():
    if line.startswith("@SQ"):
        fields = dict(x.split(":", 1) for x in line.split("\t")[1:])
        reference_lengths[fields["SN"]] = int(fields["LN"])

rows = []

print(f"Computing coverage for {len(reference_lengths)} reference genomes...\n")

for ref, genome_length in reference_lengths.items():

    result = subprocess.run(
        [
            "samtools",
            "depth",
            "-a",
            "-r",
            ref,
            str(BAM_FILE),
        ],
        capture_output=True,
        text=True,
        check=True,
    )

    depths = [
        int(line.split("\t")[2])
        for line in result.stdout.splitlines()
        if line
    ]

    if not depths:
        continue

    covered = sum(d > 0 for d in depths)

    rows.append(
        {
            "virus": REF_NAMES.get(ref, ref),
            "reference": ref,
            "genome_length": genome_length,
            "covered_bases": covered,
            "coverage_percent": round(100 * covered / genome_length, 2),
            "mean_depth": round(sum(depths) / genome_length, 2),
            "max_depth": max(depths),
        }
    )

coverage = pd.DataFrame(rows)

coverage.sort_values(
    "coverage_percent",
    ascending=False,
    inplace=True,
)

coverage.to_csv(
    RESULT_DIR / "genome_coverage.csv",
    index=False,
)

print("Coverage summary:\n")
print(coverage)

# -------------------------
# Plot
# -------------------------

detected = coverage[coverage["coverage_percent"] > 0].copy()

plt.figure(figsize=(9,6))

plt.barh(
    detected["virus"],
    detected["coverage_percent"],
)

plt.xlabel("Genome Coverage (%)")
plt.ylabel("Virus")
plt.title("Respiratory Virus Genome Coverage")

plt.xlim(0, 100)

plt.tight_layout()

plt.savefig(
    PLOT_DIR / "genome_coverage.png",
    dpi=150,
)

plt.close()

print("\nSaved:")
print("  results/genome_coverage.csv")
print("  plots/genome_coverage.png")