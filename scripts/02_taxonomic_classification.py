#!/usr/bin/env python3
"""
02_taxonomic_classification.py

Taxonomic classification of Zephyr respiratory reads via minimap2 alignment
to reference respiratory virus genomes.

Inputs:
    data/*.respiratory.fasta.gz
    references/respiratory_viruses.fasta

Outputs:
    results/aligned.bam
    results/taxonomy_summary.csv
    results/coverage_summary.csv
    plots/coverage_summary.png
"""

import os
import shutil
import subprocess
import time
import argparse
from pathlib import Path
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt
from Bio import Entrez

# ── paths ──────────────────────────────────────────────────────────────────
DATA_DIR   = Path("data")
REF_DIR    = Path("references")
RESULT_DIR = Path("results")
PLOT_DIR   = Path("plots")
LOG_FILE   = RESULT_DIR / "minimap2.log"

for d in (REF_DIR, RESULT_DIR, PLOT_DIR):
    d.mkdir(exist_ok=True)

REF_FASTA = REF_DIR / "respiratory_viruses.fasta"
SAM_FILE  = RESULT_DIR / "aligned.sam"
BAM_FILE  = RESULT_DIR / "aligned.bam"

# ── reference accessions ───────────────────────────────────────────────────
ACCESSIONS = [
    "NC_001617", "NC_009996",   # Rhinovirus A/B
    "NC_004148",                 # Metapneumovirus
    "MN908947",                  # SARS-CoV-2
    "NC_038235", "NC_001803",   # RSV A/B
    "NC_007373", "NC_002023",   # Influenza A/B
    "NC_003461", "NC_003462",   # Parainfluenza 1/2
    "NC_001796", "NC_001798",   # Parainfluenza 3/4
    "NC_001405",                 # Adenovirus C
]

REF_NAMES = {
    "NC_001617.1": "Rhinovirus-A",
    "NC_009996.1": "Rhinovirus-B",
    "NC_004148.2": "Metapneumovirus",
    "MN908947.3":  "SARS-CoV-2",
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

# ── tool check ─────────────────────────────────────────────────────────────
def check_tools() -> None:
    for tool in ["minimap2", "samtools"]:
        if shutil.which(tool) is None:
            raise RuntimeError(f"{tool} is not installed or not in PATH.")

# ── step 1: fetch references ───────────────────────────────────────────────
def fetch_references() -> None:
    if REF_FASTA.exists() and REF_FASTA.stat().st_size > 0:
        print("References already present — skipping fetch.")
        return

    # fix 1: email from environment, not hardcoded
    Entrez.email = os.environ.get("ENTREZ_EMAIL")
    if not Entrez.email:
        raise RuntimeError("Please set the ENTREZ_EMAIL environment variable.")

    print(f"Fetching {len(ACCESSIONS)} reference genomes from NCBI...")
    with REF_FASTA.open("w") as f:
        for acc in ACCESSIONS:
            print(f"  {acc}")
            # fix 2: use context manager for Entrez handle
            with Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text") as handle:
                f.write(handle.read())
            time.sleep(0.4)

    print(f"Saved {REF_FASTA}\n")

# ── step 2: align with minimap2 ────────────────────────────────────────────
def align() -> None:
    fasta_files = sorted(DATA_DIR.glob("*.respiratory.fasta.gz"))
    if not fasta_files:
        raise FileNotFoundError("No FASTA files found in data/")

    print(f"Aligning {len(fasta_files)} FASTA files with minimap2...")

    # fix 3: use context manager for log file
    with SAM_FILE.open("w") as sam_out, LOG_FILE.open("w") as log:
        try:
            subprocess.run(
                ["minimap2", "-ax", "map-ont", str(REF_FASTA)]
                + [str(f) for f in fasta_files],
                stdout=sam_out,
                stderr=log,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            # fix 4: catch subprocess failures
            print(f"minimap2 failed: {e}")
            raise

    print("Sorting and indexing BAM...")

    # fix 5: no shell=True — use Popen pipe instead
    try:
        view = subprocess.Popen(
            ["samtools", "view", "-bS", str(SAM_FILE)],
            stdout=subprocess.PIPE,
        )
        subprocess.run(
            ["samtools", "sort", "-o", str(BAM_FILE)],
            stdin=view.stdout,
            check=True,
        )
        view.stdout.close()
        view.wait()

        subprocess.run(["samtools", "index", str(BAM_FILE)], check=True)
    except subprocess.CalledProcessError as e:
        print(f"samtools failed: {e}")
        raise

    print(f"BAM ready: {BAM_FILE}\n")

# ── step 3: taxonomy counts ────────────────────────────────────────────────
def taxonomy_counts() -> pd.DataFrame:
    print("Counting mapped reads per pool per virus...")
    counts: dict = defaultdict(lambda: defaultdict(int))

    with SAM_FILE.open() as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 6:
                continue
            flag = int(fields[1])
            ref  = fields[2]
            if flag & 4 or ref == "*":
                continue

            # fix 6: document read ID format assumption explicitly
            # Zephyr read IDs are expected as: <pool-name>_<read-suffix>
            read_id = fields[0]
            pool = read_id.split("_")[0] if "_" in read_id else read_id

            virus = REF_NAMES.get(ref, ref)
            counts[pool][virus] += 1

    rows = [
        {"pool": pool, "virus": virus, "mapped_reads": n}
        for pool, refs in counts.items()
        for virus, n in refs.items()
    ]

    # fix 7: handle empty DataFrame
    df = pd.DataFrame(rows)
    if df.empty:
        print("No mapped reads found.")
        return df

    df = df.sort_values(["pool", "mapped_reads"], ascending=[True, False])
    df.to_csv(RESULT_DIR / "taxonomy_summary.csv", index=False)
    print(f"Saved results/taxonomy_summary.csv\n")
    return df

# ── step 4: coverage breadth ───────────────────────────────────────────────
def coverage_breadth() -> pd.DataFrame:
    print("Computing genome coverage breadth...")

    header = subprocess.run(
        ["samtools", "view", "-H", str(BAM_FILE)],
        capture_output=True, text=True
    ).stdout

    ref_lengths: dict = {}
    for line in header.splitlines():
        if line.startswith("@SQ"):
            parts = dict(f.split(":", 1) for f in line.split("\t")[1:])
            ref_lengths[parts["SN"]] = int(parts["LN"])

    rows = []
    for ref, length in ref_lengths.items():
        result = subprocess.run(
            ["samtools", "depth", "-a", "-r", ref, str(BAM_FILE)],
            capture_output=True, text=True
        )
        depths = [int(l.split("\t")[2]) for l in result.stdout.splitlines() if l]
        if not depths:
            continue
        covered = sum(1 for d in depths if d > 0)
        rows.append({
            "virus":            REF_NAMES.get(ref, ref),
            "reference":        ref,
            "genome_length":    length,
            "covered_bases":    covered,
            "coverage_breadth": round(100 * covered / length, 2),
            "mean_depth":       round(sum(depths) / length, 3),
        })

    # fix 8: handle empty DataFrame
    df = pd.DataFrame(rows)
    if df.empty:
        print("No coverage computed.")
        return df

    df = df.sort_values("coverage_breadth", ascending=False)
    df.to_csv(RESULT_DIR / "coverage_summary.csv", index=False)
    print(df[["virus", "coverage_breadth", "mean_depth"]].to_string())
    print(f"\nSaved results/coverage_summary.csv\n")
    return df

# ── step 5: plot ───────────────────────────────────────────────────────────
def plot_coverage(df: pd.DataFrame) -> None:
    detected = df[df["coverage_breadth"] > 0].copy()
    if detected.empty:
        print("No detected viruses to plot.")
        return

    # fix 9: sort before plotting
    detected = detected.sort_values("coverage_breadth", ascending=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Respiratory Virus Detection — Zephyr Pools", fontsize=13, fontweight="bold")

    ax = axes[0]
    colors = ["#2ecc71" if v > 90 else "#f39c12" if v > 50 else "#e74c3c"
              for v in detected["coverage_breadth"]]
    ax.barh(detected["virus"], detected["coverage_breadth"], color=colors)
    ax.axvline(90, color="green",  linestyle="--", alpha=0.5, label="90% threshold")
    ax.axvline(50, color="orange", linestyle="--", alpha=0.5, label="50% threshold")
    ax.set_xlabel("Coverage Breadth (%)")
    ax.set_title("Genome Coverage Breadth")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 105)

    ax = axes[1]
    ax.barh(detected["virus"], detected["mean_depth"], color="#3498db")
    ax.set_xlabel("Mean Depth (x)")
    ax.set_title("Mean Sequencing Depth")
    ax.set_xscale("log")

    plt.tight_layout()
    out = PLOT_DIR / "coverage_summary.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")
    plt.close()

# ── CLI ────────────────────────────────────────────────────────────────────
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Zephyr taxonomic classification")
    parser.add_argument("--fetch",  action="store_true", help="Fetch reference genomes")
    parser.add_argument("--align",  action="store_true", help="Run minimap2 alignment")
    parser.add_argument("--count",  action="store_true", help="Count mapped reads")
    parser.add_argument("--cover",  action="store_true", help="Compute coverage breadth")
    parser.add_argument("--plot",   action="store_true", help="Generate plots")
    parser.add_argument("--all",    action="store_true", help="Run all steps (default)")
    return parser.parse_args()

# ── main ───────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    args = parse_args()
    run_all = args.all or not any([args.fetch, args.align, args.count, args.cover, args.plot])

    check_tools()

    if run_all or args.fetch:
        fetch_references()
    if run_all or args.align:
        align()

    tax_df = pd.DataFrame()
    cov_df = pd.DataFrame()

    if run_all or args.count:
        tax_df = taxonomy_counts()
    if run_all or args.cover:
        cov_df = coverage_breadth()
    if run_all or args.plot:
        if cov_df.empty:
            cov_df = pd.read_csv(RESULT_DIR / "coverage_summary.csv")
        plot_coverage(cov_df)

    print("\nDone — all results in results/ and plots/")