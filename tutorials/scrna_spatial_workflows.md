# Analysis Workflows for Single-Cell and Spatial Genomics

> Reproducible workflows for single-cell and spatial transcriptomics analysis. A resource for lab members and AI agents to ensure quality and consistency.

## Table of Contents

- [Overview](#overview)
- [Setup](#setup)
- [Single-Cell RNA-seq](#single-cell-rna-seq)
  - [Processing in R (Seurat)](#processing-in-r-seurat)
    - [File Paths](#file-paths)
    - [Load Data](#load-data)
    - [Ambient RNA Removal (SoupX)](#ambient-rna-removal-soupx)
    - [Quality Control](#quality-control)
    - [Doublet Detection (scDblFinder)](#doublet-detection-scdblfinder)
    - [Filtering](#filtering)
    - [Normalization](#normalization)
    - [Dimensionality Reduction](#dimensionality-reduction)
    - [Clustering](#clustering)
    - [Cell Type Annotation](#cell-type-annotation)
    - [Save and Export Data](#save-and-export-data)
  - [Processing in Python (Scanpy)](#processing-in-python-scanpy)
- [Spatial Transcriptomics](#spatial-transcriptomics)
  - [Processing Xenium Data (Voyager)](#processing-xenium-data-voyager)
    - [Why SpatialFeatureExperiment?](#why-spatialfeatureexperiment)
    - [Understanding SFE Geometry](#understanding-sfe-geometry)
    - [Spatial Visualization](#spatial-visualization---centroids-vs-polygons)
  - [Annotation with scRNA-seq Reference (SingleR)](#annotation-with-scrna-seq-reference-singler)
    - [Why Use scRNA-seq as a Reference?](#why-use-scrna-seq-as-a-reference)
    - [Examine Annotation Quality](#examine-annotation-quality)
- [Data Interoperability](#data-interoperability)
- [Multi-Sample Considerations](#multi-sample-considerations)
- [Resources](#resources)

---

## Overview

This repository provides complete, self-contained workflow scripts for common single-cell and spatial transcriptomics analyses. Each script:

- **Runs top-to-bottom** without modification for the demo data
- **Is heavily documented** with explanations of the "why" not just the "what"
- **Follows our [lab style guide](../guides/visualization_style_guide.md)** for visualization
- **Exports interoperable formats** for R/Python cross-compatibility

### Philosophy

**Default to fewer filters and corrections.** You can always see the impact on downstream analysis and circle back if needed. Over-filtering can remove real biology.

### Demo Data

The repository includes demo data from a high-grade serous ovarian carcinoma (HGSC) sample:
- **10x Flex scRNA-seq** (`data/Xenium_HGSC_24824/`)
- **10x Xenium spatial** (`data/Flex_HGSC_24824/`)
- **Cell type markers** (`data/cellassign_markers_v2.csv`)

---

## Setup

### R Packages

```r
# Core packages
install.packages(c("Seurat", "ggplot2", "patchwork", "dplyr", "Matrix"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "SingleCellExperiment",
  "scDblFinder",
  "SingleR",
  "SpatialFeatureExperiment",
  "Voyager"
))

# CRAN
install.packages("SoupX")
```

### Python Environment

For Scanpy workflows, use the `scverse` conda environment:

```bash
mamba activate scverse
```

---

## Single-Cell RNA-seq

### Processing in R (Seurat)

**Complete script:** [`scripts/scrna_processing_seurat.R`](scripts/scrna_processing_seurat.R)

This workflow processes 10x Flex scRNA-seq data through QC, normalization, clustering, and annotation.

---

#### File Paths

Update these paths for your own data:

```r
# ---- FILE PATHS ----
filtered_h5 <- "data/Xenium_HGSC_24824/sample_filtered_feature_bc_matrix.h5"
raw_h5 <- "data/Xenium_HGSC_24824/sample_raw_feature_bc_matrix.h5"
marker_file <- "data/cellassign_markers_v2.csv"
output_dir <- "output/"
fig_dir <- "figs/"

sample_id <- "HGSC_24824"
```

---

#### Packages and Theme

We load all required packages and set up a consistent visualization theme following our [style guide](../guides/visualization_style_guide.md):

```r
library(Seurat)
library(SoupX)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)

# Okabe-Ito palette for categorical data (colorblind-friendly, up to 8 categories)
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

# Kelly palette for many categories (omitting white and black)
kelly_colors <- c(
  "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032",
  "#C2B280", "#848482", "#008856", "#E68FAC", "#0067A5",
  "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
  "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26"
)

# Expression color scale: light grey (zero) -> inverted magma (high)
expression_colors <- c("lightgrey", rev(viridisLite::magma(100)))

# Lab ggplot2 theme (minimal, no gridlines)
theme_lab <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(color = "black"),
      legend.background = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      strip.background = element_blank()
    )
}

theme_set(theme_lab())
```

---

#### Load Data

We load both filtered and raw count matrices. The raw matrix is needed for SoupX ambient RNA correction.

```r
filtered_counts <- Read10X_h5(filtered_h5)
raw_counts <- Read10X_h5(raw_h5)

cat("Filtered matrix:", nrow(filtered_counts), "genes x", ncol(filtered_counts), "cells\n")
cat("Raw matrix:", nrow(raw_counts), "genes x", ncol(raw_counts), "droplets\n")
```

```
Filtered matrix: 18129 genes x 20315 cells
Raw matrix: 39186 genes x 613580 droplets
```

---

#### Ambient RNA Removal (SoupX)

**Why this matters:** Ambient RNA from lysed cells contaminates all droplets. This can cause false positive gene expression, particularly for lowly-expressed markers, and can confound cell type identification.

```r
# Ensure same genes in both matrices
common_genes <- intersect(rownames(filtered_counts), rownames(raw_counts))
filtered_counts <- filtered_counts[common_genes, ]
raw_counts <- raw_counts[common_genes, ]

# Create SoupChannel and estimate contamination
soup_channel <- SoupChannel(tod = raw_counts, toc = filtered_counts)

# Quick clustering for SoupX (doesn't need to be perfect)
seu_temp <- CreateSeuratObject(filtered_counts)
seu_temp <- NormalizeData(seu_temp, verbose = FALSE)
seu_temp <- FindVariableFeatures(seu_temp, verbose = FALSE)
seu_temp <- ScaleData(seu_temp, verbose = FALSE)
seu_temp <- RunPCA(seu_temp, verbose = FALSE)
seu_temp <- FindNeighbors(seu_temp, dims = 1:20, verbose = FALSE)
seu_temp <- FindClusters(seu_temp, resolution = 0.5, verbose = FALSE)

soup_channel <- setClusters(soup_channel, setNames(seu_temp$seurat_clusters, colnames(seu_temp)))
soup_channel <- autoEstCont(soup_channel, verbose = FALSE)

cat("Estimated contamination fraction:", round(soup_channel$fit$rhoEst, 3), "\n")

adjusted_counts <- adjustCounts(soup_channel)
```

```
Estimated contamination fraction: 0.129
```

In this dataset, SoupX estimates ~13% ambient RNA contamination.

---

#### Quality Control

We calculate standard QC metrics and visualize their distributions to inform filtering.

**Important:** Always examine YOUR data distributions before setting thresholds! Don't blindly apply thresholds from other datasets.

```r
seu <- CreateSeuratObject(
  counts = adjusted_counts,
  project = sample_id,
  min.cells = 3,
  min.features = 200
)

seu$sample_id <- sample_id
seu[["percent_mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent_ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
```

**Histograms for QC metrics:**

Histograms are better than violin plots for setting thresholds because you can see the exact distribution shape.

```r
p_hist_ngene <- ggplot(seu@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100, fill = okabe_ito[1], color = "black", linewidth = 0.2) +
  labs(x = "Genes per Cell", y = "Count", title = "Gene Count Distribution") +
  theme_lab()

p_hist_mt <- ggplot(seu@meta.data, aes(x = percent_mt)) +
  geom_histogram(bins = 100, fill = okabe_ito[6], color = "black", linewidth = 0.2) +
  labs(x = "Mitochondrial %", y = "Count", title = "MT% Distribution") +
  theme_lab()
```

![QC Histograms](figs/qc_histograms.png)

The histograms show:
- Most cells have 500-3000 genes
- MT% has a long tail of high-MT cells that should be filtered

---

#### Doublet Detection (scDblFinder)

**Why this matters:** Doublets (droplets containing 2+ cells) can appear as intermediate cell states, confound trajectory analysis, and create artificial "hybrid" cell types.

We rely on explicit doublet detection rather than upper-limit thresholds on nCount/nFeature, which have been shown to be non-specific for doublets in ground truth studies.

```r
sce <- as.SingleCellExperiment(seu)
set.seed(42)
sce <- scDblFinder(sce)

seu$scDblFinder_score <- sce$scDblFinder.score
seu$scDblFinder_class <- sce$scDblFinder.class

doublet_rate <- mean(seu$scDblFinder_class == "doublet")
cat("Detected doublet rate:", round(doublet_rate * 100, 1), "%\n")
```

```
Detected doublet rate: 14.9%
```

![Doublet Scores](figs/doublet_scores.png)

---

#### Filtering

We filter based on:
1. **percent_mt** - remove dying/stressed cells (threshold based on distribution)
2. **Doublets** - remove cells classified as doublets by scDblFinder

We do NOT apply upper-limit thresholds on nCount/nFeature because:
- They are not specific for doublets (ground truth studies show this)
- We already run explicit doublet detection
- Over-filtering can remove real biology (e.g., highly active cells)

```r
# Look at YOUR histogram and set threshold accordingly!
max_mt_pct <- 20

seu <- subset(
  seu,
  subset = percent_mt < max_mt_pct &
    scDblFinder_class == "singlet"
)
```

```
--- Filtering Summary ---
Starting cells: 20071
After filtering: 16015
Cells removed: 4056 (20.2%)
```

![QC Histograms Post-filter](figs/qc_histograms_postfilter.png)

---

#### Normalization

SCTransform v2 performs variance stabilization and regresses out sequencing depth. It's recommended for dimensionality reduction, clustering, and integration.

**Important:** While SCTransform is great for embedding/clustering, it's still common practice to use log-normalized counts for visualization (FeaturePlots, violin plots) and differential expression analysis. We run both normalizations.

```r
# SCTransform for embedding/clustering
seu <- SCTransform(
  seu,
  method = "glmGamPoi",
  vst.flavor = "v2",
  verbose = FALSE
)

# Also compute standard log-normalization for visualization/DE
seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)
seu <- ScaleData(seu, assay = "RNA", verbose = FALSE)
```

---

#### Dimensionality Reduction

```r
seu <- RunPCA(seu, verbose = FALSE)
```

![PCA Elbow Plot](figs/pca_elbow.png)

The elbow plot shows variance explained by each PC. Choose the point where the curve starts to flatten—here, ~30 PCs.

```r
n_pcs <- 30
seu <- RunUMAP(seu, dims = 1:n_pcs, verbose = FALSE)
```

---

#### Clustering

**Important: More clusters is NOT always better!**

Consider the "granularity" of your biological question:
- If asking "are there more T cells in condition A vs B?" don't cluster at such high resolution that you get 8 T cell clusters
- More clusters = more work to characterize and explain each one
- High resolution increases risk of identifying technical/sample-specific patterns rather than biological ones

**We typically find resolution 0.2-0.3 matches most biological questions.** Start low and increase only if you need finer granularity.

```r
seu <- FindNeighbors(seu, dims = 1:n_pcs, verbose = FALSE)

# Compare a few resolutions
resolutions <- c(0.2, 0.3, 0.5)
for (res in resolutions) {
  seu <- FindClusters(seu, resolution = res, verbose = FALSE)
}
```

![Clustering Resolutions](figs/clustering_resolutions.png)

We default to resolution 0.2 (11 clusters) for this dataset:

```r
Idents(seu) <- "SCT_snn_res.0.2"
seu$seurat_clusters <- seu$SCT_snn_res.0.2
```

![UMAP Clusters](figs/umap_clusters.png)

---

#### Cell Type Annotation

We use canonical marker genes to manually annotate cell types. The marker list (`data/cellassign_markers_v2.csv`) contains genes for major cell types.

```r
markers_df <- read.csv(marker_file)

# DotPlot - use RNA assay (log-normalized) for visualization
p_dotplot <- DotPlot(
  seu,
  features = unique(unlist(marker_list)),
  assay = "RNA",
  cluster.idents = TRUE
) +
  scale_color_gradientn(colors = expression_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

![Marker DotPlot](figs/marker_dotplot.png)

**Feature plots** with our style guide color scale (light grey → inverted magma):

![Epithelial Markers](figs/markers_epithelial.png)
![Stromal Markers](figs/markers_stromal.png)
![Immune Markers](figs/markers_immune.png)

Based on marker expression, assign cell type labels:

```r
seu$cell_type <- case_when(
  seu$seurat_clusters %in% c(0, 2, 4, 9) ~ "Epithelial",
  seu$seurat_clusters %in% c(1, 10) ~ "Fibroblast",
  seu$seurat_clusters == 3 ~ "Macrophage",
  seu$seurat_clusters %in% c(5, 6) ~ "Plasma_cell",
  seu$seurat_clusters == 7 ~ "T_cell",
  seu$seurat_clusters == 8 ~ "Endothelial",
  TRUE ~ "Other"
)
```

**Cell type UMAP with shadowtext labels:**

![Cell Type UMAP](figs/umap_celltype.png)

**Automated annotation options:**
- [CellAssign](https://github.com/cook-lab/scrna-integration-pipeline) - Probabilistic annotation using marker gene matrix
- [SingleR](https://bioconductor.org/packages/SingleR/) - Reference-based annotation
- [Celltypist](https://www.celltypist.org/) - Pre-trained machine learning models

---

#### Save and Export Data

Save the processed Seurat object:

```r
saveRDS(seu, file.path(output_dir, paste0(sample_id, "_processed.rds")))
```

Export in Cell Ranger-compatible format (gzipped) for use with other tools:

```r
# Sparse count matrix
counts_matrix <- LayerData(seu, assay = "RNA", layer = "counts")
mtx_file <- file.path(export_dir, "matrix.mtx")
writeMM(counts_matrix, mtx_file)
system(paste("gzip -f", mtx_file))

# Features (genes) - Cell Ranger format
features_df <- data.frame(
  gene_id = rownames(seu),
  gene_name = rownames(seu),
  feature_type = "Gene Expression"
)
write.table(features_df, file.path(export_dir, "features.tsv"), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
system(paste("gzip -f", file.path(export_dir, "features.tsv")))

# Barcodes
write.table(colnames(seu), file.path(export_dir, "barcodes.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
system(paste("gzip -f", file.path(export_dir, "barcodes.tsv")))
```

Exported files:
- `matrix.mtx.gz` - Sparse count matrix
- `features.tsv.gz` - Gene names (Cell Ranger format)
- `barcodes.tsv.gz` - Cell barcodes
- `cell_metadata.csv` - All cell-level metadata
- `umap_coordinates.csv`, `pca_coordinates.csv` - Embeddings

---

### Processing in Python (Scanpy)

**Complete script:** [`scripts/scrna_processing_scanpy.py`](scripts/scrna_processing_scanpy.py)

**Environment:** `mamba activate scverse`

This workflow mirrors the Seurat workflow using the Scanpy/scverse ecosystem.

---

#### Key Differences from Seurat

| Step | Seurat | Scanpy |
|------|--------|--------|
| Ambient RNA | SoupX | Not included (run SoupX in R first if needed) |
| Doublet detection | scDblFinder | Scrublet (via `sc.external.pp.scrublet`) |
| Normalization | SCTransform | `normalize_total` + `log1p` |
| HVG selection | SCTransform | `highly_variable_genes` |
| Clustering | Leiden | Leiden |

---

#### Doublet Detection

Scrublet's automatic threshold can fail when the score distribution isn't clearly bimodal. We fall back to a percentile-based threshold:

```python
sc.external.pp.scrublet(adata, expected_doublet_rate=0.08)

# If automatic threshold fails, use percentile-based approach
expected_rate = 0.08
if adata.obs["predicted_doublet"].sum() / len(adata.obs) < expected_rate / 2:
    threshold = np.percentile(adata.obs["doublet_score"], (1 - expected_rate) * 100)
    adata.obs["predicted_doublet"] = adata.obs["doublet_score"] > threshold
```

![Doublet Scores (Scanpy)](figs/doublet_scores_scanpy.png)

---

#### QC and Filtering

Same philosophy as Seurat—examine histograms and set thresholds based on YOUR data:

![QC Histograms (Scanpy)](figs/qc_histograms_scanpy.png)

---

#### Clustering

![Clustering Resolutions (Scanpy)](figs/clustering_resolutions_scanpy.png)

Note: Cluster numbers differ between Seurat and Scanpy due to different random seeds and algorithm implementations. Always verify annotations with marker plots.

---

#### Cell Type Annotation

![Marker DotPlot (Scanpy)](figs/marker_dotplot_scanpy.png)

![Cell Type UMAP (Scanpy)](figs/umap_celltype_scanpy.png)

---

#### Save and Export

```python
# Save as h5ad
adata.write(f"{output_dir}/{sample_id}_processed.h5ad")

# Export Cell Ranger format for R interoperability
from scipy.io import mmwrite
mmwrite(export_dir / "matrix.mtx", counts.T)
# Then gzip the files
```

---

## Spatial Transcriptomics

### Processing Xenium Data (Voyager)

**Complete script:** [`scripts/spatial_processing_voyager.R`](scripts/spatial_processing_voyager.R)

This workflow processes 10x Xenium spatial transcriptomics data using the SpatialFeatureExperiment (SFE) and Voyager packages.

#### Why SpatialFeatureExperiment?

SpatialFeatureExperiment is an extension of SingleCellExperiment designed for spatial transcriptomics. Key advantages:

1. **Stores both centroids AND segmentation polygons** - You can choose the right representation for your visualization
2. **sf integration** - Full access to spatial statistics and geometry operations
3. **Voyager integration** - Rich visualization and spatial analysis functions

---

#### File Paths

```r
xenium_dir <- "data/Flex_HGSC_24824/"
output_dir <- "output/"
fig_dir <- "figs/"

sample_id <- "HGSC_24824_spatial"
```

---

#### Load and Preprocess Xenium Data

```r
library(SpatialFeatureExperiment)
library(Voyager)
library(BiocNeighbors)
library(scico)

sfe <- readXenium(xenium_dir, add_molecules = FALSE)

# IMPORTANT: Use gene symbols instead of Ensembl IDs
rownames(sfe) <- rowData(sfe)$Symbol

cat("Loaded SFE object:\n")
cat("  Cells:", ncol(sfe), "\n")
cat("  Genes (before filtering):", nrow(sfe), "\n")
```

```
Loaded SFE object:
  Cells: 167780
  Genes (before filtering): 487
```

**Remove control probes** - Xenium includes negative controls that should be removed:

```r
ctrl_probes <- c("NegControlProbe", "NegControlCodeword", "UnassignedCodeword")
to_drop <- grepl(paste(ctrl_probes, collapse = "|"), rownames(sfe))
sfe <- sfe[!to_drop, ]
cat("Genes after removing controls:", nrow(sfe), "\n")
```

```
Genes after removing controls: 480
```

**Filter low-quality cells and debris:**

```r
# Remove cells with very low counts
sfe <- sfe[, sfe$transcript_counts > 10]

# Remove isolated cells (likely debris) using k-nearest neighbor distances
g <- findKNN(spatialCoords(sfe)[, 1:2], k = 5, BNPARAM = AnnoyParam())
max_dist <- matrixStats::rowMaxs(g$distance)
min_dist <- matrixStats::rowMins(g$distance)
sfe$main_tissue <- !(max_dist > 100 | min_dist > 60)
sfe <- sfe[, sfe$main_tissue]
```

**Normalize by cell area** (not library size!):

```r
# For spatial data, larger cells capture more transcripts naturally
sfe <- logNormCounts(sfe, size.factor = sfe$cell_area)
```

---

#### Understanding SFE Geometry

**This is the key concept:** SpatialFeatureExperiment contains multiple geometry types that serve different purposes:

| Geometry | Type | Use Case |
|----------|------|----------|
| **Centroids** | POINT | Fast plotting, whole-tissue overviews |
| **Cell segmentation** | POLYGON | Detailed ROI views, showing cell boundaries |
| **Nucleus segmentation** | POLYGON | When you need nucleus-specific analysis |

```r
# See available geometries
colGeometryNames(sfe)
```

```
[1] "centroids" "cellSeg"   "nucSeg"
```

```r
# Cell centroids are stored in spatialCoords
head(spatialCoords(sfe))
```

```
              x_centroid y_centroid
AAAAACCCGGTATTTC  4691.655   2440.437
AAAAAGCAGTTTCGTA  2939.019   4355.684
AAAAAGCCGGTGCGCT  4746.851   2396.016
...
```

**When to use each:**
- **Centroids:** Whole-tissue visualization with many cells (fast, low memory)
- **Cell polygons:** Detailed ROI visualization (shows actual cell shapes)

---

#### Quality Control

```r
sfe <- scuttle::addPerCellQCMetrics(sfe)
```

![Spatial QC Histograms](figs/spatial_qc_histograms.png)

---

#### Spatial Visualization - Centroids vs Polygons

Following the style guide, we use:
- **Whole tissue:** `ggrastr::geom_point_rast()` with rasterized points (shape=16, size=0.1, alpha=0.8)
- **ROI:** Voyager's `plotSpatialFeature()` with segmentation polygons (linewidth=0.2, color="grey15")
- **Expression colors:** Inverted scico "lapaz" palette for histology-like appearance

**Whole tissue using centroids** (fast):

```r
library(ggrastr)
library(scico)

df <- as.data.frame(spatialCoords(sfe))
df$sum <- sfe$sum

p_counts_tissue <- ggplot(df, aes(x = x_centroid, y = y_centroid)) +
  geom_point_rast(aes(color = sum),
                  shape = 16, size = 0.1, alpha = 0.8,
                  raster.dpi = 300) +
  scale_color_gradientn(
    colours = scico(30, palette = "lapaz", direction = -1),
    name = "Total\nCounts"
  ) +
  coord_fixed() +
  theme_void()
```

![Spatial Counts Centroids](figs/spatial_counts_centroids.png)

**ROI with cell segmentation polygons** (detailed):

```r
# Define bounding box for ROI (500 microns)
coords <- spatialCoords(sfe)
center_x <- median(coords[, 1])
center_y <- median(coords[, 2])

bbox_roi <- c(
  xmin = center_x - 250,
  xmax = center_x + 250,
  ymin = center_y - 250,
  ymax = center_y + 250
)

# Plot with cell boundaries using Voyager
p_roi_counts <- plotSpatialFeature(
  sfe,
  features = "sum",
  colGeometryName = "cellSeg",
  bbox = bbox_roi,
  linewidth = 0.2,
  color = "grey15",
  alpha = 0.8
) +
  scale_fill_gradientn(
    colours = scico(30, palette = "lapaz", direction = -1),
    name = "Total\nCounts"
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )
```

![Spatial ROI Polygons](figs/spatial_roi_polygons.png)

---

#### Gene Expression Visualization

```r
epithelial_markers <- c("EPCAM", "KRT8", "KRT18", "PAX8")
epithelial_markers <- epithelial_markers[epithelial_markers %in% rownames(sfe)]

p_epi_spatial <- plotSpatialFeature(
  sfe,
  features = epithelial_markers[1],
  colGeometryName = "centroids",
  size = 0.1,
  alpha = 0.8
) +
  scale_color_gradientn(colors = expression_colors, name = epithelial_markers[1]) +
  theme_void()
```

---

#### Custom Spatial Visualization with ggplot2

For maximum control, extract geometries and use ggplot2 directly:

```r
# Subset SFE to ROI
in_roi <- coords[, 1] >= bbox_roi["xmin"] & coords[, 1] <= bbox_roi["xmax"] &
  coords[, 2] >= bbox_roi["ymin"] & coords[, 2] <= bbox_roi["ymax"]
sfe_roi <- sfe[, in_roi]

# Extract cell boundaries as sf object
cell_geom <- colGeometry(sfe_roi, "cellSeg")
cell_sf <- cell_geom
cell_sf$sum <- sfe_roi$sum

# Custom plot with style guide parameters
ggplot(cell_sf) +
  geom_sf(aes(fill = sum), color = "grey15", linewidth = 0.2, alpha = 0.8) +
  scale_fill_gradientn(
    colours = scico(30, palette = "lapaz", direction = -1),
    name = "Total\nCounts"
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )
```

![Custom Spatial ROI](figs/spatial_custom_roi.png)

---

### Annotation with scRNA-seq Reference (SingleR)

**Complete script:** [`scripts/spatial_annotation_singler.R`](scripts/spatial_annotation_singler.R)

**Prerequisites:**
1. Run `scripts/scrna_processing_seurat.R` first to generate the annotated scRNA-seq reference
2. Run `scripts/spatial_processing_voyager.R` to generate the processed spatial data

#### Why Use scRNA-seq as a Reference?

- **Higher gene coverage** - scRNA-seq typically captures more genes than spatial platforms
- **Well-characterized cell types** - Careful manual annotation from clustering analysis
- **Same biological context** - Using your own scRNA-seq ensures biological relevance
- **Transfers knowledge** - Leverages the mature scRNA-seq annotation ecosystem for spatial data

---

#### Load Reference and Query Data

```r
library(SingleR)
library(SpatialFeatureExperiment)
library(Seurat)

# Load annotated scRNA-seq reference
seu <- readRDS("output/HGSC_24824_processed.rds")

# Load processed spatial data
sfe <- readRDS("output/HGSC_24824_spatial_processed.rds")
```

---

#### Prepare for SingleR

SingleR requires log-normalized expression matrices and cell type labels:

```r
# Convert Seurat to SingleCellExperiment
sce_ref <- as.SingleCellExperiment(seu, assay = "RNA")

# Reference labels from manual annotation
ref_labels <- seu$cell_type

# Find common genes
# Note: Xenium panel has ~480 genes; scRNA-seq has ~18,000
# Common genes are limited to what's in the Xenium panel
common_genes <- intersect(rownames(sce_ref), rownames(sfe))
cat("Common genes:", length(common_genes), "\n")
```

```
Common genes: ~450
```

The number of common genes depends on your Xenium panel and how many genes from the panel are also in your scRNA-seq reference (e.g., mitochondrial genes may be excluded).

---

#### Run SingleR Annotation

```r
singler_results <- SingleR(
  test = sfe_query,
  ref = sce_ref,
  labels = ref_labels,
  de.method = "classic",
  assay.type.ref = "logcounts",
  assay.type.test = "logcounts"
)
```

---

#### Examine Annotation Quality

**Score heatmap** - shows correlation scores for each cell type:

![SingleR Score Heatmap](figs/singler_score_heatmap.png)

**Delta distribution** - confidence of predictions (higher = more confident):

```r
# Delta = difference between best and second-best score
# Low delta cells may be ambiguous or transitional
```

![SingleR Delta Distribution](figs/singler_delta_distribution.png)

---

#### Spatial Visualization of Cell Types

**Whole tissue - centroids:**

![Spatial SingleR Cell Types](figs/spatial_singler_celltype.png)

**Annotation confidence map:**

![Spatial SingleR Confidence](figs/spatial_singler_confidence.png)

**ROI with cell boundaries:**

![Spatial SingleR ROI Cell Types](figs/spatial_singler_roi_celltype.png)

**Per-cell-type spatial distribution:**

![Spatial SingleR Cell Type Facet](figs/spatial_singler_celltype_facet.png)

---

#### Compare Cell Type Proportions

```r
# Compare scRNA-seq reference vs spatial annotations
comparison_df <- data.frame(
  cell_type = all_types,
  scRNAseq = as.numeric(ref_props[all_types]),
  Spatial = as.numeric(spatial_props[all_types])
)
```

![SingleR Proportion Comparison](figs/singler_proportion_comparison.png)

Note: Proportions may differ due to:
- Spatial sampling bias (different tissue regions)
- Platform-specific detection sensitivity
- Reference annotation quality

---

## Data Interoperability

### Loading Seurat exports in Scanpy

The exported files follow Cell Ranger conventions and can be loaded with `scanpy.read_10x_mtx`:

```python
import scanpy as sc
import pandas as pd

# Load the Cell Ranger format files
adata = sc.read_10x_mtx("output/exported_data/")

# Add metadata
metadata = pd.read_csv("output/exported_data/cell_metadata.csv", index_col=0)
adata.obs = metadata

# Add embeddings
umap = pd.read_csv("output/exported_data/umap_coordinates.csv", index_col="barcode")
adata.obsm["X_umap"] = umap.values

pca = pd.read_csv("output/exported_data/pca_coordinates.csv", index_col="barcode")
adata.obsm["X_pca"] = pca.values
```

---

## Multi-Sample Considerations

The workflows in this repository focus on single-sample processing. When working with **multiple samples**, consider:

1. **Batch effects:** Technical variation between samples can obscure biological signal. Integration methods help align samples while preserving biological differences.

2. **Integration options:**
   - [Harmony](https://github.com/immunogenomics/harmony) - Fast, works well for most cases
   - [scVI/scANVI](https://scvi-tools.org/) - Deep learning approach, handles complex batch effects
   - Seurat's CCA/RPCA integration

3. **Our integration pipeline:** For multi-sample experiments, we recommend using our [scRNA-seq integration pipeline](https://github.com/cook-lab/scrna-integration-pipeline) which implements multiple integration methods with benchmarking.

4. **Sample metadata:** Always include a `sample_id` column in your metadata to track which cells came from which sample.

---

## Resources

### Lab Resources
- [Lab Style Guide](../guides/visualization_style_guide.md) - Visualization standards
- [scRNA-seq Integration Pipeline](https://github.com/cook-lab/scrna-integration-pipeline) - Multi-sample integration

### External Resources
- [Seurat Documentation](https://satijalab.org/seurat/)
- [Scanpy Documentation](https://scanpy.readthedocs.io/)
- [Voyager Documentation](https://pachterlab.github.io/voyager/)
- [scverse Tutorials](https://scverse.org/)

### Marker Gene Resources
- The `data/cellassign_markers_v2.csv` file contains markers for major cell types in ovarian/peritoneal tissues
- [CellMarker Database](http://xteam.xbio.top/CellMarker/)
- [PanglaoDB](https://panglaodb.se/)
