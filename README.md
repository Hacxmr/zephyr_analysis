# Zephyr Metagenomic Biosurveillance Analysis

**Candidate:** Mitali Raj

**Time Taken:** ~4 hours

---

## Overview

This repository contains my implementation for the coding components of the **Frontier Biodefense Fellowship – Metagenomic Biosurveillance Work Test**.

The pipeline analyzes pooled Oxford Nanopore Technologies (ONT) respiratory metatranscriptomic sequencing data from the Zephyr program to:

- summarize sequencing pools,
- classify respiratory viruses,
- estimate genome coverage, and
- explore sequence similarity using embedding-based clustering.

The implementation prioritizes simplicity, reproducibility, and modularity while using established bioinformatics tools together with lightweight machine learning techniques.

---

## Repository Structure

```text
zephyr_analysis/
├── data/
│   └── README.md
├── references/
│   └── respiratory_viruses.fasta
├── results/
│   ├── pool_summary.csv
│   ├── taxonomy_summary.csv
│   ├── genome_coverage.csv
│   ├── read_embeddings.csv
│   └── cluster_summary.csv
├── plots/
│   ├── coverage_summary.png
│   ├── genome_coverage.png
│   └── read_clusters.png
├── scripts/
│   ├── 01_prepare_data.py
│   ├── 02_taxonomic_classification.py
│   ├── 03_genome_coverage.py
│   └── 04_embedding_clustering.py
├── pyproject.toml
├── uv.lock
├── README.md
└── .gitignore
```

---

## Pipeline

### Step 1 – Data Preparation

```bash
python scripts/01_prepare_data.py
```

Summarizes each sequencing pool by computing:

- Number of reads
- Total sequenced bases
- Mean read length
- Median read length
- Minimum read length
- Maximum read length
- GC percentage

**Output**

```text
results/pool_summary.csv
```

---

### Step 2 – Taxonomic Classification

```bash
python scripts/02_taxonomic_classification.py
```

This script:

- downloads respiratory virus reference genomes from NCBI (if absent),
- aligns ONT reads using **minimap2**,
- sorts and indexes alignments using **samtools**,
- counts mapped reads for each respiratory virus, and
- computes coverage statistics.

**Outputs**

```text
results/taxonomy_summary.csv
results/coverage_summary.csv
plots/coverage_summary.png
```

---

### Step 3 – Genome Coverage

```bash
python scripts/03_genome_coverage.py
```

Computes genome-level statistics from the aligned BAM file, including:

- covered bases,
- genome coverage percentage,
- mean sequencing depth, and
- maximum sequencing depth.

**Outputs**

```text
results/genome_coverage.csv
plots/genome_coverage.png
```

---

### Step 4 – Embedding-Based Clustering

```bash
python scripts/04_embedding_clustering.py
```

Each sequencing read is represented using normalized **4-mer frequency embeddings**.

The workflow consists of:

1. Computing normalized 4-mer frequencies.
2. Reducing dimensionality using PCA.
3. Clustering reads using K-Means.
4. Visualizing clusters in two-dimensional embedding space.

**Outputs**

```text
results/read_embeddings.csv
results/cluster_summary.csv
plots/read_clusters.png
```

---

## Results

### Detected Respiratory Viruses

| Virus | Genome Coverage |
|------------------|----------------:|
| Rhinovirus-A | 100.00% |
| Rhinovirus-B | 99.06% |
| SARS-CoV-2 | 92.02% |
| Metapneumovirus | 70.06% |

No meaningful genome coverage was observed for:

- RSV A
- RSV B
- Influenza A
- Influenza B
- Parainfluenza 1–4
- Adenovirus C

---

### Coverage Analysis

Genome coverage breadth and sequencing depth were jointly used to assess detection confidence.

- **Rhinovirus-A:** Complete genome recovery (100%)
- **Rhinovirus-B:** Nearly complete genome recovery (99.06%)
- **SARS-CoV-2:** High-confidence detection (92.02%)
- **Metapneumovirus:** Moderate-confidence detection with partial genome recovery (70.06%)

Near-complete genome coverage across independent genomic regions provides substantially stronger evidence than isolated read alignments.

---

### Embedding-Based Clustering

Sequence embeddings generated five clusters after PCA projection and K-Means clustering.

The resulting clusters indicate that nucleotide composition alone captures meaningful biological structure. Cluster assignments were broadly consistent with the taxonomic classification, suggesting that similar viral genomes occupy nearby regions in embedding space.

Embedding-based clustering provides an orthogonal validation of reference-based taxonomic assignment and offers a framework for exploring potentially divergent or previously uncharacterized viral sequences.

---

## Software Requirements

### External Tools

- minimap2
- samtools

### Python Environment

Using **uv**:

```bash
uv sync
```

Or install manually:

```bash
pip install pandas biopython matplotlib scikit-learn
```

---

## Running the Pipeline

```bash
python scripts/01_prepare_data.py

python scripts/02_taxonomic_classification.py

python scripts/03_genome_coverage.py

python scripts/04_embedding_clustering.py
```

---

## Data

The Zephyr respiratory FASTA files are **not included** in this repository.

Place the downloaded `*.respiratory.fasta.gz` files inside the `data/` directory before running the pipeline.

---

## Summary

This repository implements a modular metagenomic biosurveillance workflow combining:

- sequencing pool quality assessment,
- reference-based taxonomic classification,
- genome coverage analysis for confidence estimation, and
- embedding-based clustering for exploratory sequence analysis.

The pipeline uses widely adopted bioinformatics tools (**minimap2** and **samtools**) together with lightweight machine learning methods to produce reproducible and interpretable analyses suitable for pathogen surveillance.