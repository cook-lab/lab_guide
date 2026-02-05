# =============================================================================
# scRNA-seq Processing Workflow (Scanpy)
# =============================================================================
# This script performs complete preprocessing of 10x Flex scRNA-seq data
# including QC, doublet detection, normalization, clustering, and cell type
# annotation using the Scanpy/scverse ecosystem.
#
# To use with your own data:
# 1. Update the file paths below
# 2. Run interactively in Jupyter or VS Code, section by section
# 3. Adjust QC thresholds based on YOUR data distributions (not blindly!)
#
# For detailed explanations, see the README:
# https://github.com/cook-lab/workflows
#
# Visualization follows our lab style guide:
# https://github.com/cook-lab/style_guide
#
# Environment: mamba activate scverse
# =============================================================================

# %% [markdown]
# ## File Paths

# %%
# ---- FILE PATHS ----
# Update these paths for your data
filtered_h5 = "data/Xenium_HGSC_24824/sample_filtered_feature_bc_matrix.h5"
marker_file = "data/cellassign_markers_v2.csv"
output_dir = "output/"
fig_dir = "figs/"

# Sample identifier (for multi-sample workflows)
sample_id = "HGSC_24824"

# %% [markdown]
# ## Packages

# %%
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

# Disable interactive plotting (for script execution)
plt.switch_backend("Agg")
sc.settings.autoshow = False

# Create output directories
Path(output_dir).mkdir(parents=True, exist_ok=True)
Path(fig_dir).mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## Visualization Settings
#
# Following our lab style guide:
# - Minimal theme, no gridlines
# - Colorblind-friendly palettes
# - Light grey -> inverted magma for expression

# %%
# ---- VISUALIZATION SETTINGS ----
# Following our lab style guide (https://github.com/cook-lab/style_guide)

# Set scanpy figure parameters
sc.settings.figdir = fig_dir
sc.settings.set_figure_params(
    dpi=150,
    dpi_save=300,
    frameon=False,
    fontsize=10,
    figsize=(5, 4),
)

# Okabe-Ito palette for categorical data (colorblind-friendly, up to 8 categories)
okabe_ito = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#999999",  # gray
]

# Kelly palette for many categories (omitting white and black)
kelly_colors = [
    "#F3C300",  # vivid yellow
    "#875692",  # strong purple
    "#F38400",  # vivid orange
    "#A1CAF1",  # very light blue
    "#BE0032",  # vivid red
    "#C2B280",  # grayish yellow
    "#848482",  # medium gray
    "#008856",  # vivid green
    "#E68FAC",  # strong purplish pink
    "#0067A5",  # strong blue
    "#F99379",  # strong yellowish pink
    "#604E97",  # strong violet
    "#F6A600",  # vivid orange yellow
    "#B3446C",  # strong purplish red
    "#DCD300",  # vivid greenish yellow
    "#882D17",  # strong reddish brown
    "#8DB600",  # vivid yellow green
    "#654522",  # deep yellowish brown
    "#E25822",  # vivid reddish orange
    "#2B3D26",  # dark olive green
]

# Cell type colors with biological meaning
cell_type_colors = {
    "Epithelial": "#CC6600",
    "Mesothelial": "#FFCC99",
    "Fibroblast": "#D4A574",
    "Smooth_Muscle": "#8B7355",
    "Pericyte": "#BC8F8F",
    "Endothelial": "#8B4513",
    "T_cell": "#4169E1",
    "NK_cell": "#6A5ACD",
    "B_cell": "#9370DB",
    "Plasma_cell": "#8A2BE2",
    "Macrophage": "#228B22",
    "DC": "#32CD32",
    "Neutrophil": "#90EE90",
    "Mast": "#98FB98",
    "Erythrocyte": "#DC143C",
    "Other": "#999999",
}

# Expression colormap: light grey -> inverted magma
# Create custom colormap
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl

magma_colors = plt.cm.magma(np.linspace(0, 1, 100))
expression_colors = np.vstack([[[0.83, 0.83, 0.83, 1.0]], magma_colors[::-1]])
expression_cmap = LinearSegmentedColormap.from_list("grey_magma_r", expression_colors)

# Set matplotlib defaults
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = 10

# %% [markdown]
# ## 1. Load Data

# %%
# =============================================================================
# 1. LOAD DATA
# =============================================================================
# Load the filtered count matrix from Cell Ranger output.
# Note: For Scanpy, we typically don't run SoupX (ambient RNA removal) as
# there's no direct Python equivalent with the same ease of use. If ambient
# RNA is a concern, consider running SoupX in R first.

adata = sc.read_10x_h5(filtered_h5)
adata.var_names_make_unique()

# Add sample metadata
adata.obs["sample_id"] = sample_id

print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

# %% [markdown]
# ## 2. Quality Control
#
# We calculate standard QC metrics and visualize their distributions.
#
# **Important:** Always examine YOUR data distributions before setting thresholds!
# Don't blindly apply thresholds from other datasets.
#
# **Philosophy:** Default to fewer filters/corrections. You can always see the
# impact on downstream analysis and circle back if needed. Over-filtering
# can remove real biology.

# %%
# =============================================================================
# 2. QUALITY CONTROL
# =============================================================================
# Calculate QC metrics

# Mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# Ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
)

print(f"Initial cells: {adata.n_obs}")

# %%
# QC Histograms - better than violin plots for setting thresholds
fig, axes = plt.subplots(1, 3, figsize=(12, 4))

axes[0].hist(adata.obs["n_genes_by_counts"], bins=100, color=okabe_ito[0], edgecolor="black", linewidth=0.3)
axes[0].set_xlabel("Genes per Cell")
axes[0].set_ylabel("Count")
axes[0].set_title("Gene Count Distribution")

axes[1].hist(adata.obs["total_counts"], bins=100, color=okabe_ito[1], edgecolor="black", linewidth=0.3)
axes[1].set_xlabel("UMIs per Cell")
axes[1].set_ylabel("Count")
axes[1].set_title("UMI Count Distribution")

axes[2].hist(adata.obs["pct_counts_mt"], bins=100, color=okabe_ito[5], edgecolor="black", linewidth=0.3)
axes[2].set_xlabel("Mitochondrial %")
axes[2].set_ylabel("Count")
axes[2].set_title("MT% Distribution")

plt.tight_layout()
plt.savefig(f"{fig_dir}/qc_histograms_scanpy.png", dpi=150, bbox_inches="tight")
plt.show()

# %%
# Scatter plot: total_counts vs pct_counts_mt
fig, ax = plt.subplots(figsize=(5, 4))
ax.scatter(adata.obs["total_counts"], adata.obs["pct_counts_mt"], s=1, alpha=0.3, c="black")
ax.set_xlabel("UMIs per Cell")
ax.set_ylabel("Mitochondrial %")
ax.set_title("MT% vs UMI Count")
plt.tight_layout()
plt.savefig(f"{fig_dir}/qc_scatter_mt_scanpy.png", dpi=150, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 3. Doublet Detection (Scanpy's scrublet wrapper)
#
# Scanpy includes a wrapper for Scrublet that simulates artificial doublets
# and identifies likely doublets in your data.
#
# **Why this matters:**
# - Doublets can appear as intermediate cell states
# - Can confound trajectory analysis
# - May create artificial "hybrid" cell types
#
# Note: We rely on explicit doublet detection rather than upper-limit thresholds
# on nCount/nFeature, which have been shown to be non-specific for doublets.

# %%
# =============================================================================
# 3. DOUBLET DETECTION (via Scanpy's scrublet wrapper)
# =============================================================================

# Run Scrublet via scvi-tools external implementation
# Note: Scrublet's automatic threshold detection can fail when the score
# distribution isn't clearly bimodal. For 10x Flex data, we expect ~8%
# doublets at 20k cells loaded. We manually set the threshold based on
# score percentile if the automatic detection fails.

sc.external.pp.scrublet(adata, expected_doublet_rate=0.08)

# Check if automatic thresholding failed (very few doublets called)
expected_rate = 0.08
if adata.obs["predicted_doublet"].sum() / len(adata.obs) < expected_rate / 2:
    print("Automatic threshold likely failed - setting manual threshold based on expected rate")
    # Set threshold at (100 - expected_doublet_rate) percentile
    threshold = np.percentile(adata.obs["doublet_score"], (1 - expected_rate) * 100)
    adata.obs["predicted_doublet"] = adata.obs["doublet_score"] > threshold
    print(f"Manual threshold: {threshold:.3f}")

doublet_rate = adata.obs["predicted_doublet"].sum() / len(adata.obs) * 100
print(f"Detected doublet rate: {doublet_rate:.1f}%")

# %%
# Visualize doublet scores
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(
    adata.obs.loc[~adata.obs["predicted_doublet"], "doublet_score"],
    bins=50, alpha=0.7, label="Singlet", color=okabe_ito[0]
)
ax.hist(
    adata.obs.loc[adata.obs["predicted_doublet"], "doublet_score"],
    bins=50, alpha=0.7, label="Doublet", color=okabe_ito[5]
)
ax.set_xlabel("Doublet Score")
ax.set_ylabel("Count")
ax.set_title("Scrublet Doublet Detection")
ax.legend()
plt.tight_layout()
plt.savefig(f"{fig_dir}/doublet_scores_scanpy.png", dpi=150, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 4. Filtering
#
# We filter based on:
# 1. **pct_counts_mt** - remove dying/stressed cells (threshold based on distribution)
# 2. **Doublets** - remove cells classified as doublets by Scrublet
#
# We do NOT apply upper-limit thresholds on n_genes/total_counts because:
# - They are not specific for doublets (ground truth studies show this)
# - We already run explicit doublet detection
# - Over-filtering can remove real biology (e.g., highly active cells)

# %%
# =============================================================================
# 4. FILTERING
# =============================================================================
# Look at YOUR histogram and set threshold accordingly!

max_mt_pct = 20  # Adjust based on YOUR data!

print(f"\n--- Filtering Summary ---")
print(f"Starting cells: {adata.n_obs}")

# Apply filters
adata = adata[
    (adata.obs["pct_counts_mt"] < max_mt_pct) &
    (~adata.obs["predicted_doublet"])
].copy()

print(f"After filtering: {adata.n_obs}")

# %%
# Post-filter histograms
fig, axes = plt.subplots(1, 2, figsize=(8, 4))

axes[0].hist(adata.obs["n_genes_by_counts"], bins=100, color=okabe_ito[0], edgecolor="black", linewidth=0.3)
axes[0].set_xlabel("Genes per Cell")
axes[0].set_ylabel("Count")
axes[0].set_title("Gene Count (filtered)")

axes[1].hist(adata.obs["pct_counts_mt"], bins=100, color=okabe_ito[5], edgecolor="black", linewidth=0.3)
axes[1].set_xlabel("Mitochondrial %")
axes[1].set_ylabel("Count")
axes[1].set_title("MT% (filtered)")

plt.tight_layout()
plt.savefig(f"{fig_dir}/qc_histograms_postfilter_scanpy.png", dpi=150, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 5. Normalization
#
# We use scran-based pooling for size factor estimation (via scanpy's implementation)
# followed by log transformation. This is more robust than simple library size
# normalization for heterogeneous cell populations.
#
# We keep raw counts in `adata.raw` for downstream analysis like differential expression.

# %%
# =============================================================================
# 5. NORMALIZATION
# =============================================================================
# Store raw counts for later use (DE, visualization)
adata.layers["counts"] = adata.X.copy()

# Normalize to median total counts and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Store normalized data in raw slot for visualization
adata.raw = adata

# %% [markdown]
# ## 6. Feature Selection & Scaling

# %%
# =============================================================================
# 6. FEATURE SELECTION & SCALING
# =============================================================================
# Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print(f"Highly variable genes: {adata.var['highly_variable'].sum()}")

# Scale data (only on HVGs for efficiency)
sc.pp.scale(adata, max_value=10)

# %% [markdown]
# ## 7. Dimensionality Reduction

# %%
# =============================================================================
# 7. DIMENSIONALITY REDUCTION
# =============================================================================
# PCA
sc.tl.pca(adata, svd_solver="arpack", n_comps=50)

# Elbow plot
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(range(1, 51), adata.uns["pca"]["variance_ratio"], "o-", markersize=3)
ax.set_xlabel("PC")
ax.set_ylabel("Variance Ratio")
ax.set_title("PCA Elbow Plot")
plt.tight_layout()
plt.savefig(f"{fig_dir}/pca_elbow_scanpy.png", dpi=150, bbox_inches="tight")
plt.show()

# %%
# UMAP
n_pcs = 30  # Adjust based on elbow plot
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
sc.tl.umap(adata)

# %% [markdown]
# ## 8. Clustering
#
# **Important: More clusters is NOT always better!**
#
# Consider the "granularity" of your biological question:
# - If asking "are there more T cells in condition A vs B?" don't cluster
#   at such high resolution that you get 8 T cell clusters
# - More clusters = more work to characterize and explain each one
# - High resolution increases risk of identifying technical/sample-specific
#   patterns rather than biological ones
#
# **We typically find resolution 0.2-0.3 matches most biological questions.**
# Start low and increase only if you need finer granularity.

# %%
# =============================================================================
# 8. CLUSTERING
# =============================================================================
# Cluster at a few resolutions to compare
resolutions = [0.2, 0.3, 0.5]

for res in resolutions:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")

# %%
# Visualize clustering at different resolutions
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for i, res in enumerate(resolutions):
    n_clusters = adata.obs[f"leiden_{res}"].nunique()
    palette = okabe_ito[:n_clusters] if n_clusters <= 8 else kelly_colors[:n_clusters]

    sc.pl.umap(
        adata,
        color=f"leiden_{res}",
        palette=palette,
        legend_loc="on data",
        legend_fontsize=8,
        title=f"Resolution {res}",
        ax=axes[i],
        show=False,
    )

plt.tight_layout()
plt.savefig(f"{fig_dir}/clustering_resolutions_scanpy.png", dpi=150, bbox_inches="tight")
plt.show()

# %%
# Set default clustering - start with low resolution!
adata.obs["leiden"] = adata.obs["leiden_0.2"]

# Final UMAP with clusters
n_clusters = adata.obs["leiden"].nunique()
palette = okabe_ito[:n_clusters] if n_clusters <= 8 else kelly_colors[:n_clusters]

fig, ax = plt.subplots(figsize=(7, 6))
sc.pl.umap(
    adata,
    color="leiden",
    palette=palette,
    legend_loc="on data",
    legend_fontsize=10,
    title="UMAP - Clusters (res=0.2)",
    ax=ax,
    show=False,
)
plt.tight_layout()
plt.savefig(f"{fig_dir}/umap_clusters_scanpy.png", dpi=150, bbox_inches="tight")
plt.show()

print(f"Number of clusters: {n_clusters}")

# %% [markdown]
# ## 9. Cell Type Annotation (Manual)
#
# We use canonical marker genes to manually annotate cell types.
#
# For automated annotation options, see:
# - CellAssign: https://github.com/cook-lab/scrna-integration-pipeline
# - Celltypist (pre-trained models)
# - scANVI (semi-supervised)

# %%
# =============================================================================
# 9. CELL TYPE ANNOTATION (Manual)
# =============================================================================
# Load marker genes
markers_df = pd.read_csv(marker_file)
marker_genes = markers_df["Gene"].tolist()

# Check which markers are present in our data
present_markers = [g for g in marker_genes if g in adata.var_names]
print(f"Marker genes found in data: {len(present_markers)} / {len(marker_genes)}")

# %%
# DotPlot of marker genes by cluster
# Use raw slot (log-normalized) for visualization
fig, ax = plt.subplots(figsize=(14, 5))
sc.pl.dotplot(
    adata,
    var_names=present_markers,
    groupby="leiden",
    standard_scale="var",
    cmap=expression_cmap,
    ax=ax,
    show=False,
)
plt.tight_layout()
plt.savefig(f"{fig_dir}/marker_dotplot_scanpy.png", dpi=150, bbox_inches="tight")
plt.show()

# %%
# Feature plots for key lineage markers
epithelial_markers = ["EPCAM", "KRT8", "KRT18", "PAX8"]
epithelial_markers = [g for g in epithelial_markers if g in adata.var_names]

if epithelial_markers:
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    for i, gene in enumerate(epithelial_markers):
        ax = axes[i // 2, i % 2]
        sc.pl.umap(
            adata, color=gene, cmap=expression_cmap,
            ax=ax, show=False, title=gene,
            vmin=0
        )
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/markers_epithelial_scanpy.png", dpi=150, bbox_inches="tight")
    plt.show()

# %%
# Stromal markers
stromal_markers = ["COL1A1", "DCN", "ACTA2", "PECAM1"]
stromal_markers = [g for g in stromal_markers if g in adata.var_names]

if stromal_markers:
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    for i, gene in enumerate(stromal_markers):
        ax = axes[i // 2, i % 2]
        sc.pl.umap(
            adata, color=gene, cmap=expression_cmap,
            ax=ax, show=False, title=gene,
            vmin=0
        )
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/markers_stromal_scanpy.png", dpi=150, bbox_inches="tight")
    plt.show()

# %%
# Immune markers
immune_markers = ["PTPRC", "CD3E", "CD68", "MS4A1"]
immune_markers = [g for g in immune_markers if g in adata.var_names]

if immune_markers:
    n_markers = len(immune_markers)
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    for i, gene in enumerate(immune_markers):
        ax = axes[i // 2, i % 2]
        sc.pl.umap(
            adata, color=gene, cmap=expression_cmap,
            ax=ax, show=False, title=gene,
            vmin=0
        )
    # Hide empty subplots if fewer than 4 markers
    for j in range(n_markers, 4):
        axes[j // 2, j % 2].axis("off")
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/markers_immune_scanpy.png", dpi=150, bbox_inches="tight")
    plt.show()

# %%
# ---- MANUAL ANNOTATION ----
# Based on the marker expression patterns, assign cell type labels.
# UPDATE THIS SECTION based on your data!
#
# For this dataset (based on dot plot marker expression):
# Note: Cluster numbers may differ slightly between Seurat and Scanpy

# First, let's see the cluster sizes to help with mapping
print("Cluster sizes:")
print(adata.obs["leiden"].value_counts().sort_index())

# %%
# Create cell type mapping based on marker expression
# ADJUST THIS based on YOUR marker plots!
# Note: Cluster numbers differ between Seurat and Scanpy due to different
# algorithms and random seeds - always verify with your dot plot!
cluster_to_celltype = {
    "0": "Epithelial",    # EPCAM+, KRT8+, KRT18+
    "1": "Fibroblast",    # COL1A1+, DCN+, PDGFRA+
    "2": "T_cell",        # CD3E+, PTPRC+
    "3": "Macrophage",    # CD68+, C1QA+, C1QB+
    "4": "Epithelial",    # EPCAM+, mixed markers
    "5": "Endothelial",   # PECAM1+, VWF+
    "6": "Plasma_cell",   # IGHG1+, MZB1+
    "7": "Fibroblast",    # COL1A1+, low expression
}

# Apply mapping
adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_to_celltype).fillna("Other")

print("\nCell type counts:")
print(adata.obs["cell_type"].value_counts())

# %%
# ---- CELL TYPE UMAP ----
# Final UMAP with cell type annotations

# Get colors for present cell types
present_celltypes = adata.obs["cell_type"].unique()
ct_palette = [cell_type_colors.get(ct, "#999999") for ct in present_celltypes]

fig, ax = plt.subplots(figsize=(7, 6))
sc.pl.umap(
    adata,
    color="cell_type",
    palette=cell_type_colors,
    legend_loc="right margin",
    legend_fontsize=10,
    title="Cell Type Annotations",
    ax=ax,
    show=False,
)
plt.tight_layout()
plt.savefig(f"{fig_dir}/umap_celltype_scanpy.png", dpi=300, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 10. Save Processed Data

# %%
# =============================================================================
# 10. SAVE PROCESSED DATA
# =============================================================================
# Save as h5ad (AnnData native format)
adata.write(f"{output_dir}/{sample_id}_processed.h5ad")
print(f"Saved AnnData object to: {output_dir}/{sample_id}_processed.h5ad")

# %% [markdown]
# ## 11. Export Interoperable Data
#
# Export data in formats that can be read by other tools (R, etc.)

# %%
# =============================================================================
# 11. EXPORT INTEROPERABLE DATA
# =============================================================================
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
import gzip

export_dir = Path(output_dir) / "exported_data_scanpy"
export_dir.mkdir(exist_ok=True)

# Get raw counts from layers
counts = adata.layers["counts"]
if not isinstance(counts, csr_matrix):
    counts = csr_matrix(counts)

# Export sparse count matrix (Matrix Market format, gzipped)
# Note: mmwrite doesn't support gzip directly, so we write then compress
mmwrite(str(export_dir / "matrix.mtx"), counts.T)  # Transpose for Cell Ranger format

# Gzip the matrix file
import shutil
with open(export_dir / "matrix.mtx", "rb") as f_in:
    with gzip.open(export_dir / "matrix.mtx.gz", "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
(export_dir / "matrix.mtx").unlink()  # Remove uncompressed file

# Export features (genes) - Cell Ranger format
features_df = pd.DataFrame({
    "gene_id": adata.var_names,
    "gene_name": adata.var_names,
    "feature_type": "Gene Expression"
})
features_df.to_csv(export_dir / "features.tsv.gz", sep="\t", index=False, header=False, compression="gzip")

# Export barcodes
pd.DataFrame(adata.obs_names).to_csv(
    export_dir / "barcodes.tsv.gz", sep="\t", index=False, header=False, compression="gzip"
)

# Export cell metadata
adata.obs.to_csv(export_dir / "cell_metadata.csv")

# Export UMAP coordinates
umap_df = pd.DataFrame(
    adata.obsm["X_umap"],
    index=adata.obs_names,
    columns=["UMAP1", "UMAP2"]
)
umap_df.to_csv(export_dir / "umap_coordinates.csv")

# Export PCA coordinates
pca_df = pd.DataFrame(
    adata.obsm["X_pca"],
    index=adata.obs_names,
    columns=[f"PC{i+1}" for i in range(adata.obsm["X_pca"].shape[1])]
)
pca_df.to_csv(export_dir / "pca_coordinates.csv")

print(f"\nExported interoperable data to: {export_dir}")
print("Files (Cell Ranger format):")
print("  - matrix.mtx.gz (sparse count matrix)")
print("  - features.tsv.gz (gene names)")
print("  - barcodes.tsv.gz (cell barcodes)")
print("Additional files:")
print("  - cell_metadata.csv (all cell-level metadata)")
print("  - umap_coordinates.csv")
print("  - pca_coordinates.csv")

# %% [markdown]
# ## Session Info

# %%
# =============================================================================
# SESSION INFO
# =============================================================================
print("\n--- Session Info ---")
print(f"scanpy: {sc.__version__}")
print(f"numpy: {np.__version__}")
print(f"pandas: {pd.__version__}")
sc.logging.print_versions()
