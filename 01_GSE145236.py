# -----------------------------
# Imports
# -----------------------------
from __future__ import annotations
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.io import mmread
import os
import time
import matplotlib.pyplot as plt

# -----------------------------
# Setup
# -----------------------------
SEED = 42
np.random.seed(SEED)

OUTPUT_DIR = "scanpy_01"
os.makedirs(OUTPUT_DIR, exist_ok=True)

sc.settings.figdir = OUTPUT_DIR
sc.set_figure_params(dpi=80, facecolor="white")

runtime_log = []

def timed(name, func):
    start = time.time()
    result = func()
    end = time.time()
    runtime_log.append((name, end - start))
    print(f"{name}: {end - start:.2f} sec")
    return result

# -----------------------------
# Load Data
# -----------------------------
matrix = mmread("GSM4319249_injured/matrix.mtx.gz").T.tocsr()
genes = pd.read_csv("GSM4319249_injured/features.tsv.gz", header=None, sep="\t")
barcodes = pd.read_csv("GSM4319249_injured/barcodes.tsv.gz", header=None)

adata = sc.AnnData(matrix)
adata.var_names = genes[1]
adata.obs_names = barcodes[0]
adata.var_names_make_unique()

print(adata)

# -----------------------------
# QC Metrics
# -----------------------------
adata.var["mt"] = adata.var_names.str.startswith("mt-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

timed("QC metrics", lambda: sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
))

# QC plots
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save="_qc_violin_01.png",
    show=False
)

sc.pl.scatter(
    adata,
    "total_counts",
    "n_genes_by_counts",
    color="pct_counts_mt",
    save="_qc_scatter_01.png",
    show=False
)

# -----------------------------
# Filtering
# -----------------------------
timed("Filter cells", lambda: sc.pp.filter_cells(adata, min_genes=100))
timed("Filter genes", lambda: sc.pp.filter_genes(adata, min_cells=3))

adata = adata[adata.obs["pct_counts_mt"] < 10, :].copy()
print("After filtering:", adata.shape)

# -----------------------------
# Normalization
# -----------------------------
adata.layers["counts"] = adata.X.copy()

timed("Normalize", lambda: sc.pp.normalize_total(adata))
timed("Log1p", lambda: sc.pp.log1p(adata))

# -----------------------------
# HVG
# -----------------------------
timed("HVG", lambda: sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat"
))

adata = adata[:, adata.var.highly_variable].copy()

# -----------------------------
# Scaling
# -----------------------------
timed("Scaling", lambda: sc.pp.scale(adata))

# -----------------------------
# PCA (ARPACK)
# -----------------------------
timed("PCA", lambda: sc.tl.pca(
    adata,
    svd_solver="arpack",
    random_state=SEED
))

# Save PCA outputs
pd.DataFrame(adata.obsm["X_pca"][:, :2]).to_csv(f"{OUTPUT_DIR}/pca_embeddings_01.csv", index=False)
pd.DataFrame(adata.uns["pca"]["variance_ratio"]).to_csv(f"{OUTPUT_DIR}/variance_explained_01.csv", index=False)

# -----------------------------
# Neighbors + Clustering
# -----------------------------
timed("Neighbors", lambda: sc.pp.neighbors(adata, n_pcs=10))
timed("Leiden", lambda: sc.tl.leiden(adata, random_state=SEED))

# -----------------------------
# UMAP
# -----------------------------
timed("UMAP", lambda: sc.tl.umap(adata, random_state=SEED))

# =============================
# 🔥 PLOTS (with legends + axis control)
# =============================

# -----------------------------
# PCA Plot
# -----------------------------
fig = sc.pl.pca(
    adata,
    color="leiden",
    legend_loc="right margin",
    title="PCA (PC1 vs PC2)",
    show=False,
    return_fig=True
)

ax = fig.axes[0]
ax.invert_xaxis()
ax.tick_params(labelsize=10)
ax.locator_params(nbins=5)

plt.savefig(f"{OUTPUT_DIR}/pca_pc1_pc2_01.png", bbox_inches="tight")
plt.close()

# -----------------------------
# UMAP Plot
# -----------------------------
fig = sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="right margin",
    title="UMAP (clusters)",
    show=False,
    return_fig=True
)

ax = fig.axes[0]
ax.invert_xaxis()
ax.tick_params(labelsize=10)
ax.locator_params(nbins=5)

plt.savefig(f"{OUTPUT_DIR}/umap_clusters_01.png", bbox_inches="tight")
plt.close()

# -----------------------------
# Marker Genes
# -----------------------------
timed("Markers", lambda: sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon"
))

result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names

markers = []
for g in groups:
    for i in range(10):
        markers.append({
            "cluster": g,
            "gene": result["names"][g][i],
            "logFC": result["logfoldchanges"][g][i]
        })

pd.DataFrame(markers).to_csv(f"{OUTPUT_DIR}/top_markers_01.csv", index=False)

# -----------------------------
# Top marker genes
# -----------------------------
top_genes = [result["names"][g][0] for g in groups]
print("Top genes:", top_genes)

# -----------------------------
# UMAP marker overlays
# -----------------------------
fig = sc.pl.umap(
    adata,
    color=top_genes,
    ncols=3,
    show=False,
    return_fig=True
)

for ax in fig.axes:
    ax.invert_xaxis()
    ax.tick_params(labelsize=8)

plt.savefig(f"{OUTPUT_DIR}/umap_top_markers_01.png", bbox_inches="tight")
plt.close()

# -----------------------------
# PCA marker overlays
# -----------------------------
fig = sc.pl.pca(
    adata,
    color=top_genes,
    ncols=3,
    show=False,
    return_fig=True
)

for ax in fig.axes:
    ax.invert_xaxis()
    ax.tick_params(labelsize=8)

plt.savefig(f"{OUTPUT_DIR}/pca_top_markers_01.png", bbox_inches="tight")
plt.close()

# -----------------------------
# Runtime Summary
# -----------------------------
runtime_df = pd.DataFrame(runtime_log, columns=["Step", "Time_sec"])
print(runtime_df)

runtime_df.to_csv(f"{OUTPUT_DIR}/runtime_summary_scanpy_01.csv", index=False)