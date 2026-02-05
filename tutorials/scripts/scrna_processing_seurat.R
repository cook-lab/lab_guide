# ============================================================================
# scRNA-seq Processing Workflow (Seurat)
# ============================================================================
# This script performs complete preprocessing of 10x Flex scRNA-seq data
# including ambient RNA removal, QC, doublet detection, normalization,
# clustering, and cell type annotation.
#
# To use with your own data:
# 1. Update the file paths below
# 2. Run interactively in RStudio, section by section
# 3. Adjust QC thresholds based on YOUR data distributions (not blindly!)
#
# For detailed explanations, see the README:
# https://github.com/cook-lab/workflows
#
# Visualization follows our lab style guide:
# https://github.com/cook-lab/style_guide
# ============================================================================

# ---- FILE PATHS ----
# Update these paths for your data
filtered_h5 <- "data/Xenium_HGSC_24824/sample_filtered_feature_bc_matrix.h5"
raw_h5 <- "data/Xenium_HGSC_24824/sample_raw_feature_bc_matrix.h5"
marker_file <- "data/cellassign_markers_v2.csv"
output_dir <- "output/"
fig_dir <- "figs/"

# Sample identifier (for multi-sample workflows)
sample_id <- "HGSC_24824"

# ---- PACKAGES ----
library(Seurat)
library(SoupX)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggrastr)
library(shadowtext)

# ---- VISUALIZATION THEME ----
# Following our lab style guide (https://github.com/cook-lab/style_guide)
# - Minimal theme, no gridlines
# - Sans-serif fonts (Arial/Helvetica)
# - Colorblind-friendly palettes

# Okabe-Ito palette for categorical data (colorblind-friendly, up to 8 categories)
okabe_ito <- c(
"#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999"  # gray
)

# Kelly palette for many categories (omitting colors 1-2: white and black)
# Use when you have more than 8 categories
kelly_colors <- c(
  "#F3C300", # vivid yellow
  "#875692", # strong purple
  "#F38400", # vivid orange
  "#A1CAF1", # very light blue
  "#BE0032", # vivid red
  "#C2B280", # grayish yellow
  "#848482", # medium gray
  "#008856", # vivid green
  "#E68FAC", # strong purplish pink
  "#0067A5", # strong blue
  "#F99379", # strong yellowish pink
  "#604E97", # strong violet
  "#F6A600", # vivid orange yellow
  "#B3446C", # strong purplish red
  "#DCD300", # vivid greenish yellow
  "#882D17", # strong reddish brown
  "#8DB600", # vivid yellow green
  "#654522", # deep yellowish brown
  "#E25822", # vivid reddish orange
  "#2B3D26"  # dark olive green
)

# Cell type palette - organized by lineage with strategic visual hierarchy
# From lab style guide: fibroblasts recede as background, epithelial stands out
# - Epithelial: warm orange (prominent, often cell of interest)
# - Mesenchymal: neutral nude/cream (background stromal)
# - Lymphoid: blue-to-purple gradient (cool contrast)
# - Myeloid: green family (distinct from lymphoid)
# - Vascular: muted maroon/rose
cell_type_colors <- c(
  # Epithelial - warm, prominent
  "Epithelial"    = "#E6A141",
  "Mesothelial"   = "#D4A574",
  # Mesenchymal - neutral background
  "Fibroblast"    = "#DDD5CA",
  "Mesenchymal"   = "#DDD5CA",
  "Smooth_Muscle" = "#D14E6C",
  # Vascular
  "Pericyte"      = "#B87A7A",
  "Endothelial"   = "#7D4E4E",
  # Lymphoid - blue-purple gradient
  "T_cell"        = "#87CEFA",
  "NK_cell"       = "#56AFC4",
  "B_cell"        = "#5665B6",
  "Plasma_cell"   = "#8A5DAF",
  # Myeloid - green family
  "Macrophage"    = "#8FBC8F",
  "DC"            = "#2E8B57",
  "Neutrophil"    = "#6B8E23",
  "Mast"          = "#8B9B6B",
  # Other
  "Erythrocyte"   = "#CD5C5C",
  "Other"         = "#A0A0A0"
)

# Expression color scale: light grey (zero) -> inverted magma (high)
# This follows our style guide for continuous expression data
expression_colors <- c("lightgrey", rev(viridisLite::magma(100)))

# Lab ggplot2 theme
theme_lab <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold")
    )
}

theme_set(theme_lab())

# Create output directories if they don't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. LOAD DATA
# ============================================================================
# We load both filtered and raw count matrices. The raw matrix is needed for
# SoupX ambient RNA correction - it contains droplets that Cell Ranger
# identified as empty, which helps estimate the ambient RNA profile.

# Load filtered counts (cells identified by Cell Ranger)
filtered_counts <- Read10X_h5(filtered_h5)

# Load raw counts (all droplets, including empty ones)
raw_counts <- Read10X_h5(raw_h5)

# Quick look at the data dimensions
cat("Filtered matrix:", nrow(filtered_counts), "genes x", ncol(filtered_counts), "cells\n")
cat("Raw matrix:", nrow(raw_counts), "genes x", ncol(raw_counts), "droplets\n")

# ============================================================================
# 2. AMBIENT RNA REMOVAL (SoupX)
# ============================================================================
# SoupX estimates the ambient RNA profile from empty droplets and removes this
# contamination from cell-containing droplets. This is particularly important
# for samples with high ambient RNA (e.g., tissues with lysed cells).
#
# Why this matters:
# - Ambient RNA can cause false positive gene expression
# - Particularly problematic for lowly-expressed markers
# - Can confound cell type identification

# Ensure both matrices have the same genes (required for SoupX)
# This can happen with 10x Flex data where probe sets may differ
common_genes <- intersect(rownames(filtered_counts), rownames(raw_counts))
cat("Common genes between filtered and raw:", length(common_genes), "\n")

filtered_counts <- filtered_counts[common_genes, ]
raw_counts <- raw_counts[common_genes, ]

# Create SoupChannel object
soup_channel <- SoupChannel(tod = raw_counts, toc = filtered_counts)

# SoupX needs clustering to estimate contamination.
# Quick preliminary clustering (doesn't need to be perfect).
seu_temp <- CreateSeuratObject(filtered_counts)
seu_temp <- NormalizeData(seu_temp, verbose = FALSE)
seu_temp <- FindVariableFeatures(seu_temp, verbose = FALSE)
seu_temp <- ScaleData(seu_temp, verbose = FALSE)
seu_temp <- RunPCA(seu_temp, verbose = FALSE)
seu_temp <- FindNeighbors(seu_temp, dims = 1:20, verbose = FALSE)
seu_temp <- FindClusters(seu_temp, resolution = 0.5, verbose = FALSE)

# Add clustering to SoupChannel and estimate contamination
soup_channel <- setClusters(soup_channel, setNames(seu_temp$seurat_clusters, colnames(seu_temp)))
soup_channel <- autoEstCont(soup_channel, verbose = FALSE)

cat("Estimated contamination fraction:", round(soup_channel$fit$rhoEst, 3), "\n")

# Adjust counts (remove ambient RNA)
# Round to nearest integer for consistency with count-based tools
adjusted_counts <- round(adjustCounts(soup_channel))

# Clean up
rm(seu_temp, soup_channel, raw_counts, filtered_counts)
gc()

# ============================================================================
# 3. CREATE SEURAT OBJECT & CALCULATE QC METRICS
# ============================================================================
# Now we create our main Seurat object with the SoupX-corrected counts
# and calculate standard QC metrics.

seu <- CreateSeuratObject(
  counts = adjusted_counts,
  project = sample_id,
  min.cells = 3,       # Only keep genes expressed in >= 3 cells
  min.features = 200   # Only keep cells with >= 200 genes
)

seu$sample_id <- sample_id

# Calculate mitochondrial percentage (high MT% = stressed/dying cells)
seu[["percent_mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# Calculate ribosomal percentage (can be useful for identifying certain cell states)
seu[["percent_ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

cat("Initial cells after SoupX:", ncol(seu), "\n")

# ============================================================================
# 4. QC VISUALIZATION
# ============================================================================
# Visualize QC metrics to inform filtering decisions.
# IMPORTANT: Always examine YOUR data distributions before setting thresholds!
# Don't blindly apply thresholds from other datasets.
#
# Philosophy: Default to FEWER filters/corrections. You can always see the
# impact on downstream analysis and circle back if needed. Over-filtering
# can remove real biology.

# Histograms are better than violin plots for setting thresholds
# because you can see the exact distribution shape

p_hist_ngene <- ggplot(seu@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100, fill = okabe_ito[1], color = "black", linewidth = 0.2) +
  labs(x = "Genes per Cell", y = "Count", title = "Gene Count Distribution") +
  theme_lab()

p_hist_ncount <- ggplot(seu@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100, fill = okabe_ito[2], color = "black", linewidth = 0.2) +
  labs(x = "UMIs per Cell", y = "Count", title = "UMI Count Distribution") +
  theme_lab()

p_hist_mt <- ggplot(seu@meta.data, aes(x = percent_mt)) +
  geom_histogram(bins = 100, fill = okabe_ito[6], color = "black", linewidth = 0.2) +
  labs(x = "Mitochondrial %", y = "Count", title = "MT% Distribution") +
  theme_lab()

p_qc_histograms <- p_hist_ngene + p_hist_ncount + p_hist_mt
ggsave(file.path(fig_dir, "qc_histograms.png"), p_qc_histograms,
       width = 12, height = 4, dpi = 150)

# Scatter plot: nCount vs percent_mt
# High MT% cells often have lower total counts (dying cells)
p_scatter_mt <- ggplot(seu@meta.data, aes(x = nCount_RNA, y = percent_mt)) +
  geom_point(alpha = 0.3, size = 0.5) +
  labs(x = "UMIs per Cell", y = "Mitochondrial %", title = "MT% vs UMI Count") +
  theme_lab()

ggsave(file.path(fig_dir, "qc_scatter_mt.png"), p_scatter_mt,
       width = 5, height = 4, dpi = 150)

# ============================================================================
# 5. DOUBLET DETECTION (scDblFinder)
# ============================================================================
# scDblFinder identifies likely doublets (droplets containing 2+ cells).
# It simulates artificial doublets and trains a classifier to identify them.
#
# Why this matters:
# - Doublets can appear as intermediate cell states
# - Can confound trajectory analysis
# - May create artificial "hybrid" cell types
#
# Note: We rely on explicit doublet detection rather than upper-limit thresholds
# on nCount/nFeature, which have been shown to be non-specific for doublets.

sce <- as.SingleCellExperiment(seu)
set.seed(42)
sce <- scDblFinder(sce)

seu$scDblFinder_score <- sce$scDblFinder.score
seu$scDblFinder_class <- sce$scDblFinder.class

doublet_rate <- mean(seu$scDblFinder_class == "doublet")
cat("Detected doublet rate:", round(doublet_rate * 100, 1), "%\n")

# Visualize doublet scores
p_doublet <- ggplot(seu@meta.data, aes(x = scDblFinder_score, fill = scDblFinder_class)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("singlet" = okabe_ito[1], "doublet" = okabe_ito[6])) +
  labs(x = "Doublet Score", y = "Count", fill = "Classification",
       title = "scDblFinder Doublet Detection") +
  theme_lab()

ggsave(file.path(fig_dir, "doublet_scores.png"), p_doublet,
       width = 6, height = 4, dpi = 150)

rm(sce)
gc()

# ============================================================================
# 6. FILTERING
# ============================================================================
# Apply QC filters based on the visualizations above.
#
# IMPORTANT: Look at YOUR histograms and set thresholds accordingly!
# The values below are based on THIS dataset - yours may differ.
#
# We filter on:
# 1. percent_mt - remove dying/stressed cells (set based on distribution)
# 2. Doublets - remove cells classified as doublets by scDblFinder
#
# We do NOT apply upper-limit thresholds on nCount/nFeature because:
# - They are not specific for doublets (ground truth studies show this)
# - We already run explicit doublet detection
# - Over-filtering can remove real biology (e.g., highly active cells)

# ---- SET YOUR THRESHOLD HERE ----
# Look at the MT% histogram and set a threshold that removes the long tail
# of high-MT cells while keeping the main population
max_mt_pct <- 20  # Adjust based on YOUR data!

cat("\n--- Filtering Summary ---\n")
cat("Starting cells:", ncol(seu), "\n")

n_prefilter <- ncol(seu)

seu <- subset(
  seu,
  subset = percent_mt < max_mt_pct &
    scDblFinder_class == "singlet"
)

cat("After filtering:", ncol(seu), "\n")
cat("Cells removed:", n_prefilter - ncol(seu),
    "(", round((n_prefilter - ncol(seu)) / n_prefilter * 100, 1), "%)\n")

# Post-filter histograms
p_hist_post_ngene <- ggplot(seu@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100, fill = okabe_ito[1], color = "black", linewidth = 0.2) +
  labs(x = "Genes per Cell", y = "Count", title = "Gene Count (filtered)") +
  theme_lab()

p_hist_post_mt <- ggplot(seu@meta.data, aes(x = percent_mt)) +
  geom_histogram(bins = 100, fill = okabe_ito[6], color = "black", linewidth = 0.2) +
  labs(x = "Mitochondrial %", y = "Count", title = "MT% (filtered)") +
  theme_lab()

p_qc_post <- p_hist_post_ngene + p_hist_post_mt
ggsave(file.path(fig_dir, "qc_histograms_postfilter.png"), p_qc_post,
       width = 8, height = 4, dpi = 150)

# ============================================================================
# 7. NORMALIZATION
# ============================================================================
# SCTransform v2 performs variance stabilization and regresses out sequencing
# depth. It's recommended for dimensionality reduction, clustering, and
# integration.
#
# IMPORTANT: While SCTransform is great for embedding/clustering, it's still
# common practice to use log-normalized counts for:
# - Visualization (FeaturePlots, violin plots)
# - Differential expression analysis
#
# We run both normalizations and set the default assay appropriately.

# Increase memory limit for parallel processing (needed for large datasets)
options(future.globals.maxSize = 2000 * 1024^2)  # 2GB

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

# ============================================================================
# 8. DIMENSIONALITY REDUCTION
# ============================================================================
# PCA for initial dimensionality reduction, then UMAP for visualization.

seu <- RunPCA(seu, verbose = FALSE)

# Elbow plot to determine number of PCs
p_elbow <- ElbowPlot(seu, ndims = 50) +
  theme_lab() +
  labs(title = "PCA Elbow Plot",
       subtitle = "Choose dims where curve flattens")

ggsave(file.path(fig_dir, "pca_elbow.png"), p_elbow,
       width = 6, height = 4, dpi = 150)

n_pcs <- 30  # Adjust based on elbow plot

seu <- RunUMAP(seu, dims = 1:n_pcs, verbose = FALSE)

# ============================================================================
# 9. CLUSTERING
# ============================================================================
# Leiden clustering groups cells with similar expression profiles.
#
# IMPORTANT: More clusters is NOT always better!
# Consider the "granularity" of your biological question:
# - If asking "are there more T cells in condition A vs B?" don't cluster
#   at such high resolution that you get 8 T cell clusters
# - More clusters = more work to characterize and explain each one
# - High resolution increases risk of identifying technical/sample-specific
#   patterns rather than biological ones
#
# We typically find resolution 0.2-0.3 matches most biological questions.
# Start low and increase only if you need finer granularity.

seu <- FindNeighbors(seu, dims = 1:n_pcs, verbose = FALSE)

# Cluster at a few resolutions to compare
resolutions <- c(0.2, 0.3, 0.5)
for (res in resolutions) {
  seu <- FindClusters(seu, resolution = res, verbose = FALSE)
}

# Visualize clustering at different resolutions
# Use Kelly palette for many clusters
cluster_cols <- paste0("SCT_snn_res.", resolutions)
p_res_list <- lapply(cluster_cols, function(col) {
  n_clusters <- length(unique(seu@meta.data[[col]]))
  # Choose appropriate palette based on number of clusters
  if (n_clusters <= 8) {
    pal <- okabe_ito[1:n_clusters]
  } else {
    pal <- kelly_colors[1:n_clusters]
  }

  DimPlot(seu, group.by = col, label = TRUE, label.size = 3) +
    scale_color_manual(values = pal) +
    theme_lab() +
    theme(legend.position = "none") +
    labs(title = paste0("Resolution ", gsub("SCT_snn_res.", "", col)))
})

p_resolutions <- wrap_plots(p_res_list, ncol = 3)
ggsave(file.path(fig_dir, "clustering_resolutions.png"), p_resolutions,
       width = 12, height = 4, dpi = 150)

# Set default clustering - start with low resolution!
Idents(seu) <- "SCT_snn_res.0.2"
seu$seurat_clusters <- seu$SCT_snn_res.0.2

# Final UMAP with clusters
n_clusters <- length(unique(seu$seurat_clusters))
if (n_clusters <= 8) {
  cluster_pal <- okabe_ito[1:n_clusters]
} else {
  cluster_pal <- kelly_colors[1:n_clusters]
}

p_umap_clusters <- DimPlot(seu, label = TRUE, label.size = 4) +
  scale_color_manual(values = cluster_pal) +
  theme_lab() +
  labs(title = "UMAP - Clusters (res=0.2)")

ggsave(file.path(fig_dir, "umap_clusters.png"), p_umap_clusters,
       width = 7, height = 6, dpi = 150)

# ============================================================================
# 10. CELL TYPE ANNOTATION (Manual)
# ============================================================================
# Manual annotation using canonical marker genes.
# We use the lab's marker gene list for major cell types.
#
# For automated annotation options, see:
# - CellAssign: https://github.com/cook-lab/scrna-integration-pipeline
# - SingleR (reference-based)
# - Celltypist (pre-trained models)

markers_df <- read.csv(marker_file)
marker_genes <- markers_df$Gene

present_markers <- marker_genes[marker_genes %in% rownames(seu)]
cat("Marker genes found in data:", length(present_markers), "/", length(marker_genes), "\n")

# Create marker list for each cell type
cell_types <- colnames(markers_df)[-1]
marker_list <- lapply(cell_types, function(ct) {
  markers_df$Gene[markers_df[[ct]] == 1]
})
names(marker_list) <- cell_types
marker_list <- lapply(marker_list, function(x) x[x %in% rownames(seu)])
marker_list <- marker_list[sapply(marker_list, length) > 0]

# DotPlot of marker genes by cluster
# Use RNA assay (log-normalized) for visualization
p_dotplot <- DotPlot(
  seu,
  features = unique(unlist(marker_list)),
  assay = "RNA",
  cluster.idents = TRUE
) +
  theme_lab() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colors = expression_colors) +
  labs(title = "Marker Gene Expression by Cluster")

ggsave(file.path(fig_dir, "marker_dotplot.png"), p_dotplot,
       width = 14, height = 5, dpi = 150)

# Feature plots for key lineage markers
# Use RNA assay for visualization with light grey -> inverted magma scale

# Helper function for feature plots with correct color scale
plot_markers <- function(seu, markers, title_prefix) {
  markers <- markers[markers %in% rownames(seu)]
  if (length(markers) == 0) return(NULL)

  FeaturePlot(seu, features = markers, ncol = 2, order = TRUE) &
    scale_color_gradientn(colors = expression_colors) &
    theme_lab()
}

# Epithelial markers
epithelial_markers <- c("EPCAM", "KRT8", "KRT18", "PAX8")
p_epi <- plot_markers(seu, epithelial_markers, "Epithelial")
if (!is.null(p_epi)) {
  ggsave(file.path(fig_dir, "markers_epithelial.png"), p_epi,
         width = 8, height = 8, dpi = 150)
}

# Stromal markers
stromal_markers <- c("COL1A1", "DCN", "ACTA2", "PECAM1")
p_stromal <- plot_markers(seu, stromal_markers, "Stromal")
if (!is.null(p_stromal)) {
  ggsave(file.path(fig_dir, "markers_stromal.png"), p_stromal,
         width = 8, height = 8, dpi = 150)
}

# Immune markers
immune_markers <- c("PTPRC", "CD3E", "CD68", "MS4A1")
p_immune <- plot_markers(seu, immune_markers, "Immune")
if (!is.null(p_immune)) {
  ggsave(file.path(fig_dir, "markers_immune.png"), p_immune,
         width = 8, height = 8, dpi = 150)
}

# ---- MANUAL ANNOTATION ----
# Based on the marker expression patterns, assign cell type labels.
# UPDATE THIS SECTION based on your data!
#
# For this dataset (based on dot plot marker expression):
# - Cluster 1: COL1A1+, DCN+ -> Fibroblast
# - Cluster 3: CD68+, C1QA+ -> Macrophage
# - Cluster 7: CD3E+, PTPRC+ -> T_cell
# - Clusters 0, 2, 4, 9: EPCAM+, KRT8+, KRT18+ -> Epithelial
# - Cluster 8: PECAM1+, VWF+ -> Endothelial
# - Clusters 5, 6: IGHG1+, MZB1+ -> Plasma_cell
# - Cluster 10: Low marker expression, some COL1A1 -> Fibroblast (uncertain)

seu$cell_type <- case_when(
  seu$seurat_clusters %in% c(0, 2, 4, 9) ~ "Epithelial",
  seu$seurat_clusters %in% c(1, 10) ~ "Fibroblast",
  seu$seurat_clusters == 3 ~ "Macrophage",
  seu$seurat_clusters %in% c(5, 6) ~ "Plasma_cell",
  seu$seurat_clusters == 7 ~ "T_cell",
  seu$seurat_clusters == 8 ~ "Endothelial",
  TRUE ~ "Other"
)

# ---- CELL TYPE UMAP WITH SHADOWTEXT LABELS ----
# Following lab style guide: direct labels with shadowtext for legibility

# Calculate centroids for each cell type
umap_df <- data.frame(
  UMAP1 = Embeddings(seu, "umap")[, 1],
  UMAP2 = Embeddings(seu, "umap")[, 2],
  cell_type = seu$cell_type
)

label_positions <- umap_df %>%
  group_by(cell_type) %>%
  summarize(
    UMAP1 = mean(UMAP1),
    UMAP2 = mean(UMAP2),
    .groups = "drop"
  )

# UMAP with cell type colors and shadowtext labels
p_celltype <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point_rast(
    aes(color = cell_type),
    size = 0.3,
    alpha = 0.8,
    raster.dpi = 300
  ) +
  geom_shadowtext(
    data = label_positions,
    aes(label = cell_type),
    size = 3.5,
    color = "black",
    bg.color = "white",
    fontface = "bold"
  ) +
  scale_color_manual(values = cell_type_colors) +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Cell Type Annotations")

ggsave(file.path(fig_dir, "umap_celltype.png"), p_celltype,
       width = 6, height = 5, dpi = 300)

# ============================================================================
# 11. SAVE PROCESSED DATA
# ============================================================================

saveRDS(seu, file.path(output_dir, paste0(sample_id, "_processed.rds")))
cat("Saved Seurat object to:", file.path(output_dir, paste0(sample_id, "_processed.rds")), "\n")

# ============================================================================
# 12. EXPORT INTEROPERABLE DATA
# ============================================================================
# Export data in Cell Ranger-compatible format for use with other tools.
# Files are gzipped following Cell Ranger conventions.

export_dir <- file.path(output_dir, "exported_data")
dir.create(export_dir, showWarnings = FALSE)

# Export sparse count matrix (Matrix Market format, gzipped)
counts_matrix <- LayerData(seu, assay = "RNA", layer = "counts")

# Write matrix.mtx.gz
mtx_file <- file.path(export_dir, "matrix.mtx")
writeMM(counts_matrix, mtx_file)
system(paste("gzip -f", mtx_file))

# Export features (genes) - Cell Ranger format: gene_id, gene_name, feature_type
features_df <- data.frame(
  gene_id = rownames(seu),
  gene_name = rownames(seu),
  feature_type = "Gene Expression"
)
features_file <- file.path(export_dir, "features.tsv")
write.table(features_df, features_file, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
system(paste("gzip -f", features_file))

# Export barcodes
barcodes_file <- file.path(export_dir, "barcodes.tsv")
write.table(colnames(seu), barcodes_file,
            row.names = FALSE, col.names = FALSE, quote = FALSE)
system(paste("gzip -f", barcodes_file))

# Export cell metadata (not gzipped - this is our addition)
write.csv(seu@meta.data, file.path(export_dir, "cell_metadata.csv"), row.names = TRUE)

# Export UMAP coordinates
umap_coords <- as.data.frame(Embeddings(seu, "umap"))
umap_coords$barcode <- rownames(umap_coords)
write.csv(umap_coords, file.path(export_dir, "umap_coordinates.csv"), row.names = FALSE)

# Export PCA coordinates
pca_coords <- as.data.frame(Embeddings(seu, "pca"))
pca_coords$barcode <- rownames(pca_coords)
write.csv(pca_coords, file.path(export_dir, "pca_coordinates.csv"), row.names = FALSE)

cat("\nExported interoperable data to:", export_dir, "\n")
cat("Files (Cell Ranger format):\n")
cat("  - matrix.mtx.gz (sparse count matrix)\n")
cat("  - features.tsv.gz (gene names)\n")
cat("  - barcodes.tsv.gz (cell barcodes)\n")
cat("Additional files:\n")
cat("  - cell_metadata.csv (all cell-level metadata)\n")
cat("  - umap_coordinates.csv\n")
cat("  - pca_coordinates.csv\n")

# ============================================================================
# SESSION INFO
# ============================================================================
cat("\n--- Session Info ---\n")
sessionInfo()
