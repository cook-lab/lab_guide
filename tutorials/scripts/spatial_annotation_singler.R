# ============================================================================
# Spatial Data Annotation with SingleR
# ============================================================================
# This script annotates spatial transcriptomics data (SpatialFeatureExperiment)
# using SingleR with a scRNA-seq reference dataset as the reference.
#
# SingleR is a reference-based annotation tool that:
# 1. Correlates each cell's expression profile with reference cell types
# 2. Uses iterative refinement to identify the best match
# 3. Provides confidence scores for each annotation
#
# Why use scRNA-seq as a reference for spatial data?
# - scRNA-seq typically has higher gene coverage and sensitivity
# - Well-characterized cell types from careful manual annotation
# - Same tissue/sample context ensures biological relevance
# - Transfers annotations from a data type where clustering/annotation
#   is more established to spatial data
#
# Prerequisites:
# - Run scrna_processing_seurat.R first to generate the annotated scRNA-seq
# - Run spatial_processing_voyager.R to generate the processed spatial data
#
# For detailed explanations, see the README:
# https://github.com/cook-lab/workflows
#
# Visualization follows our lab style guide:
# https://github.com/cook-lab/style_guide
# ============================================================================

# ---- FILE PATHS ----
scrna_rds <- "output/HGSC_24824_processed.rds"  # Annotated scRNA-seq Seurat object
spatial_rds <- "output/HGSC_24824_spatial_processed.rds"  # Processed SFE object
output_dir <- "output/"
fig_dir <- "figs/"

sample_id <- "HGSC_24824_spatial"

# ---- PACKAGES ----
library(SingleR)
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(Voyager)
library(SingleCellExperiment)
library(Seurat)
library(scuttle)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrastr)
library(sf)
library(scico)
library(viridis)

# ---- VISUALIZATION SETTINGS ----
# Following our lab style guide (https://github.com/cook-lab/style_guide)

# Okabe-Ito palette for categorical data (≤8 categories)
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

# Cell type palette - organized by lineage with strategic visual hierarchy
# From lab style guide: fibroblasts recede as background, epithelial stands out
celltype_palette <- c(
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

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. LOAD DATA
# ============================================================================
# Load both the scRNA-seq reference (Seurat) and spatial query (SFE)

cat("Loading scRNA-seq reference...\n")
seu <- readRDS(scrna_rds)
cat("  Cells:", ncol(seu), "\n")
cat("  Genes:", nrow(seu), "\n")
cat("  Cell types:", paste(unique(seu$cell_type), collapse = ", "), "\n")

cat("\nLoading spatial data...\n")
sfe <- readRDS(spatial_rds)
cat("  Cells:", ncol(sfe), "\n")
cat("  Genes:", nrow(sfe), "\n")

# ============================================================================
# 2. PREPARE REFERENCE
# ============================================================================
# SingleR expects SingleCellExperiment objects with:
# - Log-normalized expression matrix
# - Cell type labels
#
# We convert the Seurat object and use our manual annotations as labels.

cat("\nPreparing scRNA-seq reference for SingleR...\n")

# Convert Seurat to SingleCellExperiment
sce_ref <- as.SingleCellExperiment(seu, assay = "RNA")

# Ensure we have log-normalized data
# The 'logcounts' assay should be present from NormalizeData()
if (!"logcounts" %in% assayNames(sce_ref)) {
  # If not present, check if there's a 'data' layer we can use
  if ("data" %in% assayNames(sce_ref)) {
    assay(sce_ref, "logcounts") <- assay(sce_ref, "data")
  } else {
    stop("No log-normalized counts found in Seurat object")
  }
}

# Use our manual cell type annotations as reference labels
ref_labels <- seu$cell_type
cat("Reference cell types:\n")
print(table(ref_labels))

# ============================================================================
# 3. PREPARE QUERY (SPATIAL DATA)
# ============================================================================
# The spatial data (SFE) is already a SingleCellExperiment subclass.
# We need to ensure it has log-normalized counts for SingleR.

cat("\nPreparing spatial query data...\n")

# SFE already has logcounts from our Voyager workflow
if (!"logcounts" %in% assayNames(sfe)) {
  stop("No log-normalized counts in SFE. Run logNormCounts() first.")
}

# Find common genes between reference and query
# Note: scRNA-seq has more genes than Xenium panel
common_genes <- intersect(rownames(sce_ref), rownames(sfe))
cat("Common genes between reference and spatial:", length(common_genes), "\n")
cat("  Reference total genes:", nrow(sce_ref), "\n")
cat("  Spatial panel genes:", nrow(sfe), "\n")

# Subset to common genes
sce_ref <- sce_ref[common_genes, ]
sfe_query <- sfe[common_genes, ]

# ============================================================================
# 4. RUN SINGLER ANNOTATION
# ============================================================================
# SingleR performs reference-based annotation using correlation of expression
# profiles. It iteratively refines predictions using only marker genes for
# the most likely cell types.
#
# Key parameters:
# - labels: Reference cell type labels
# - de.method: Method to identify marker genes ("classic" = fast, "wilcox" = robust)
# - assay.type.ref / assay.type.test: Which assay to use (logcounts)

cat("\nRunning SingleR annotation...\n")
cat("This may take a few minutes depending on cell numbers...\n")

# Run SingleR
singler_results <- SingleR(
  test = sfe_query,
  ref = sce_ref,
  labels = ref_labels,
  de.method = "classic",  # Fast differential expression
  assay.type.ref = "logcounts",
  assay.type.test = "logcounts"
)

cat("SingleR annotation complete.\n")

# ============================================================================
# 5. EXAMINE RESULTS
# ============================================================================
# SingleR returns:
# - labels: Best annotation for each cell
# - scores: Correlation scores for each cell type
# - delta.next: Difference between best and second-best score (confidence)
# - pruned.labels: Labels with low-confidence predictions set to NA

cat("\n--- SingleR Results Summary ---\n")

# Predicted cell type distribution
cat("\nPredicted cell types:\n")
print(table(singler_results$labels))

# Check pruned labels (low confidence = NA)
n_pruned <- sum(is.na(singler_results$pruned.labels))
cat("\nLow-confidence cells (pruned):", n_pruned,
    "(", round(n_pruned / ncol(sfe_query) * 100, 1), "%)\n")

# Add annotations to SFE object
sfe$singler_labels <- singler_results$labels
sfe$singler_pruned <- singler_results$pruned.labels
sfe$singler_delta <- singler_results$delta.next

# ============================================================================
# 6. VISUALIZE ANNOTATION QUALITY
# ============================================================================

# ---- 6a. Score Heatmap ----
# Shows the correlation scores for each cell type across all cells
# Good annotations show high scores for one cell type and low for others

p_scoreheat <- plotScoreHeatmap(singler_results)
ggsave(file.path(fig_dir, "singler_score_heatmap.png"), p_scoreheat,
       width = 8, height = 6, dpi = 150)

# ---- 6b. Delta Distribution ----
# Delta = difference between best and second-best score
# Higher delta = more confident annotation
# Low delta cells may be ambiguous or transitional states

delta_df <- data.frame(
  delta = singler_results$delta.next,
  cell_type = singler_results$labels
)

p_delta <- ggplot(delta_df, aes(x = cell_type, y = delta, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = celltype_palette) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = "Predicted Cell Type",
    y = "Delta (Confidence)",
    title = "SingleR Annotation Confidence",
    subtitle = "Higher delta = more confident prediction"
  )

ggsave(file.path(fig_dir, "singler_delta_distribution.png"), p_delta,
       width = 8, height = 5, dpi = 150)

# ---- 6c. Delta Histogram ----
p_delta_hist <- ggplot(delta_df, aes(x = delta)) +
  geom_histogram(bins = 50, fill = okabe_ito[1], color = "black", linewidth = 0.2) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(
    x = "Delta Score",
    y = "Count",
    title = "Annotation Confidence Distribution",
    subtitle = "Red line: suggested threshold for low-confidence cells"
  )

ggsave(file.path(fig_dir, "singler_delta_histogram.png"), p_delta_hist,
       width = 6, height = 4, dpi = 150)

# ============================================================================
# 7. SPATIAL VISUALIZATION OF CELL TYPES
# ============================================================================
# Following style guide:
# - Whole tissue: rasterized points with celltype_palette
# - ROI: cell segmentation polygons with celltype_palette
# - shape=16, size=0.1, alpha=0.8, coord_fixed(), theme_void()

# ---- 7a. Whole Tissue - Centroids ----
coords <- spatialCoords(sfe)
spatial_df <- data.frame(
  x = coords[, 1],
  y = coords[, 2],
  cell_type = sfe$singler_labels,
  delta = sfe$singler_delta
)

# Ensure palette covers all cell types
present_types <- unique(spatial_df$cell_type)
type_pal <- celltype_palette[names(celltype_palette) %in% present_types]

# Add any missing types with gray
missing_types <- setdiff(present_types, names(type_pal))
if (length(missing_types) > 0) {
  extra_colors <- rep("#A0A0A0", length(missing_types))
  names(extra_colors) <- missing_types
  type_pal <- c(type_pal, extra_colors)
}

p_spatial_celltype <- ggplot(spatial_df, aes(x = x, y = y)) +
  geom_point_rast(aes(color = cell_type),
                  shape = 16, size = 0.1, alpha = 0.8,
                  raster.dpi = 300) +
  scale_color_manual(values = type_pal, name = "Cell Type") +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 2)) +
  coord_fixed() +
  theme_void() +
  labs(title = "SingleR Cell Type Annotations")

ggsave(file.path(fig_dir, "spatial_singler_celltype.png"), p_spatial_celltype,
       width = 7, height = 5, dpi = 300)

# ---- 7b. Confidence Map ----
# Use inverted lapaz for quantitative confidence
p_spatial_confidence <- ggplot(spatial_df, aes(x = x, y = y)) +
  geom_point_rast(aes(color = delta),
                  shape = 16, size = 0.1, alpha = 0.8,
                  raster.dpi = 300) +
  scale_color_gradientn(
    colours = scico(30, palette = "lapaz", direction = -1),
    name = "Confidence\n(Delta)"
  ) +
  coord_fixed() +
  theme_void() +
  labs(title = "Annotation Confidence Map")

ggsave(file.path(fig_dir, "spatial_singler_confidence.png"), p_spatial_confidence,
       width = 6, height = 5, dpi = 300)

# ---- 7c. ROI with Cell Segmentation Polygons ----
# For detailed view, use cell segmentation polygons

# Define ROI (same as in Voyager workflow)
center_x <- median(coords[, 1])
center_y <- median(coords[, 2])

bbox_roi <- c(
  xmin = center_x - 250,
  xmax = center_x + 250,
  ymin = center_y - 250,
  ymax = center_y + 250
)

in_roi <- coords[, 1] >= bbox_roi["xmin"] & coords[, 1] <= bbox_roi["xmax"] &
  coords[, 2] >= bbox_roi["ymin"] & coords[, 2] <= bbox_roi["ymax"]

cat("\nROI contains", sum(in_roi), "cells\n")

# Plot ROI with cell type colors using Voyager
p_roi_celltype <- plotSpatialFeature(
  sfe,
  features = "singler_labels",
  colGeometryName = "cellSeg",
  bbox = bbox_roi,
  linewidth = 0.2,
  color = "grey15",
  alpha = 0.8
) +
  scale_fill_manual(values = type_pal, name = "Cell Type") +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(title = "ROI Cell Types (SingleR)")

ggsave(file.path(fig_dir, "spatial_singler_roi_celltype.png"), p_roi_celltype,
       width = 5, height = 4, dpi = 300, bg = "white")

# ---- 7d. Per-Cell-Type Spatial Maps ----
# Faceted plot showing each cell type separately

p_celltype_facet <- ggplot(spatial_df, aes(x = x, y = y)) +
  geom_point_rast(
    data = spatial_df %>% select(-cell_type),
    color = "lightgrey",
    shape = 16,
    size = 0.05,
    alpha = 0.3,
    raster.dpi = 100
  ) +
  geom_point_rast(
    aes(color = cell_type),
    shape = 16,
    size = 0.1,
    alpha = 0.8,
    raster.dpi = 100
  ) +
  scale_color_manual(values = type_pal) +
  facet_wrap(~cell_type, ncol = 3) +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold")
  ) +
  labs(title = "Cell Type Spatial Distribution")

ggsave(file.path(fig_dir, "spatial_singler_celltype_facet.png"), p_celltype_facet,
       width = 12, height = 10, dpi = 150)

# ============================================================================
# 8. COMPARE WITH REFERENCE
# ============================================================================
# Summary statistics comparing spatial and scRNA-seq cell type proportions

ref_props <- table(seu$cell_type) / ncol(seu)
spatial_props <- table(sfe$singler_labels) / ncol(sfe)

# Create comparison data frame
all_types <- union(names(ref_props), names(spatial_props))
comparison_df <- data.frame(
  cell_type = all_types,
  scRNAseq = as.numeric(ref_props[all_types]),
  Spatial = as.numeric(spatial_props[all_types])
)
comparison_df[is.na(comparison_df)] <- 0

# Reshape for plotting
comparison_long <- tidyr::pivot_longer(
  comparison_df,
  cols = c("scRNAseq", "Spatial"),
  names_to = "dataset",
  values_to = "proportion"
)

p_comparison <- ggplot(comparison_long, aes(x = cell_type, y = proportion, fill = dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("scRNAseq" = okabe_ito[1], "Spatial" = okabe_ito[2])) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Cell Type",
    y = "Proportion",
    fill = "Dataset",
    title = "Cell Type Proportions: scRNA-seq vs Spatial",
    subtitle = "Based on SingleR annotations for spatial data"
  )

ggsave(file.path(fig_dir, "singler_proportion_comparison.png"), p_comparison,
       width = 10, height = 5, dpi = 150)

cat("\n--- Cell Type Proportion Comparison ---\n")
print(comparison_df)

# ============================================================================
# 9. SAVE ANNOTATED DATA
# ============================================================================

saveRDS(sfe, file.path(output_dir, paste0(sample_id, "_annotated.rds")))
cat("\nSaved annotated SFE object to:",
    file.path(output_dir, paste0(sample_id, "_annotated.rds")), "\n")

# Also save the SingleR results for reference
saveRDS(singler_results, file.path(output_dir, paste0(sample_id, "_singler_results.rds")))
cat("Saved SingleR results to:",
    file.path(output_dir, paste0(sample_id, "_singler_results.rds")), "\n")

# ============================================================================
# SESSION INFO
# ============================================================================
cat("\n--- Session Info ---\n")
sessionInfo()
