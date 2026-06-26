#!/usr/bin/env python3
"""
04_embedding_clustering.py

Generate sequence embeddings using k-mer frequencies and
cluster respiratory reads.

Inputs
------
data/*.respiratory.fasta.gz

Outputs
-------
results/read_embeddings.csv
results/cluster_summary.csv
plots/read_clusters.png
"""

from pathlib import Path
import gzip
import itertools

import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

DATA_DIR = Path("data")
RESULT_DIR = Path("results")
PLOT_DIR = Path("plots")

RESULT_DIR.mkdir(exist_ok=True)
PLOT_DIR.mkdir(exist_ok=True)

K = 4
N_CLUSTERS = 5
MAX_READS = 5000


def open_fasta(path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


# all possible 4-mers
kmers = ["".join(x) for x in itertools.product("ACGT", repeat=K)]
kmer_index = {k: i for i, k in enumerate(kmers)}


def kmer_vector(sequence):

    sequence = sequence.upper()

    vec = [0] * len(kmers)

    for i in range(len(sequence) - K + 1):

        kmer = sequence[i:i + K]

        if set(kmer) <= {"A", "C", "G", "T"}:
            vec[kmer_index[kmer]] += 1

    total = sum(vec)

    if total:
        vec = [x / total for x in vec]

    return vec


embeddings = []
metadata = []

count = 0

files = sorted(DATA_DIR.glob("*.respiratory.fasta.gz"))

if not files:
    raise FileNotFoundError("No FASTA files found.")

print(f"Processing {len(files)} FASTA files")

for fasta in files:

    pool = fasta.name.replace(".respiratory.fasta.gz", "")

    with open_fasta(fasta) as handle:

        for record in SeqIO.parse(handle, "fasta"):

            embeddings.append(kmer_vector(str(record.seq)))

            metadata.append({
                "read": record.id,
                "pool": pool,
                "length": len(record.seq)
            })

            count += 1

            if count >= MAX_READS:
                break

    if count >= MAX_READS:
        break

print(f"Embedded {len(embeddings)} reads")

X = pd.DataFrame(embeddings)

pca = PCA(n_components=2, random_state=42)

coords = pca.fit_transform(X)

kmeans = KMeans(
    n_clusters=N_CLUSTERS,
    random_state=42,
    n_init=10
)

clusters = kmeans.fit_predict(X)

df = pd.DataFrame(metadata)

df["PC1"] = coords[:, 0]
df["PC2"] = coords[:, 1]
df["cluster"] = clusters

df.to_csv(
    RESULT_DIR / "read_embeddings.csv",
    index=False
)

summary = (
    df.groupby(["cluster", "pool"])
      .size()
      .reset_index(name="reads")
)

summary.to_csv(
    RESULT_DIR / "cluster_summary.csv",
    index=False
)

plt.figure(figsize=(8, 6))

for c in sorted(df.cluster.unique()):

    subset = df[df.cluster == c]

    plt.scatter(
        subset.PC1,
        subset.PC2,
        s=12,
        alpha=0.7,
        label=f"Cluster {c}"
    )

plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.title("Respiratory Read Embedding Clusters")

plt.legend()

plt.tight_layout()

plt.savefig(
    PLOT_DIR / "read_clusters.png",
    dpi=150
)

plt.close()

print("\nSaved:")
print(" results/read_embeddings.csv")
print(" results/cluster_summary.csv")
print(" plots/read_clusters.png")