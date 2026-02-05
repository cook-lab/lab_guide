# ============================================================================
# Spatial Transcriptomics Processing Workflow (Voyager)
# ============================================================================
# This script processes 10x Xenium spatial transcriptomics data using the
# SpatialFeatureExperiment and Voyager packages for visualization and
# spatial statistics.
#
# To use with your own data:
# 1. Update the file paths below
# 2. Run interactively in RStudio, section by section
#
# For detailed explanations, see the README:
# https://github.com/cook-lab/workflows
#
# Visualization follows our lab style guide:
# https://github.com/cook-lab/style_guide
# ============================================================================

# ---- FILE PATHS ----
xenium_dir <- "data/Flex_HGSC_24824/"
output_dir <- "output/"
fig_dir <- "figs/"

sample_id <- "HGSC_24824_spatial"

# ---- PACKAGES ----
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(Voyager)
library(SingleCellExperiment)
library(scuttle)
library(BiocNeighbors)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrastr)
library(sf)
library(scico)
library(viridis)

# ---- VISUALIZATION THEME ----
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
# 1. LOAD XENIUM DATA
# ============================================================================
# SpatialFeatureExperiment (SFE) is an extension of SingleCellExperiment
# designed for spatial transcriptomics data. It stores:
# - Count matrix (genes x cells)
# - Cell metadata (colData)
# - Gene metadata (rowData)
# - Spatial coordinates (spatialCoords)
# - Cell/nucleus segmentation polygons (colGeometries)
#
# Key advantage: SFE stores BOTH cell centroids AND segmentation polygons,
# allowing flexible visualization choices.

sfe <- readXenium(xenium_dir, add_molecules = FALSE)

# IMPORTANT: Use gene symbols instead of Ensembl IDs for interpretability
rownames(sfe) <- rowData(sfe)$Symbol

cat("Loaded SFE object:\n")
cat("  Cells:", ncol(sfe), "\n")
cat("  Genes (before filtering):", nrow(sfe), "\n")

# ============================================================================
# 2. REMOVE CONTROL PROBES
# ============================================================================
# Xenium includes control probes that should be removed before analysis:
# - NegControlProbe: Negative control probes
# - NegControlCodeword: Negative control codewords
# - UnassignedCodeword: Probes that couldn't be decoded
#
# These are NOT biological genes and will confuse downstream analysis.

ctrl_probes <- c("NegControlProbe", "NegControlCodeword", "UnassignedCodeword")
to_drop <- grepl(paste(ctrl_probes, collapse = "|"),
                 rownames(sfe), ignore.case = FALSE)

cat("\nRemoving control probes:", sum(to_drop), "\n")
sfe <- sfe[!to_drop, ]
cat("Genes after removing controls:", nrow(sfe), "\n")

# ============================================================================
# 3. UNDERSTANDING SFE GEOMETRY
# ============================================================================
# SpatialFeatureExperiment contains multiple geometry types:
#
# 1. Cell centroids (POINT geometry) - stored in spatialCoords()
#    - Fast to plot, good for whole-tissue overviews
#    - Use when you have many cells or need quick visualization
#
# 2. Cell segmentation polygons (POLYGON geometry) - stored in colGeometries()
#    - Shows actual cell boundaries from segmentation
#    - Better for ROI visualization, seeing cell shapes
#    - More memory intensive to plot
#
# Let's explore what geometries are available:

cat("\nAvailable column geometries:\n")
print(colGeometryNames(sfe))

# Cell centroids are the default spatialCoords
cat("\nSpatial coordinates (centroids) - first 5 cells:\n")
print(head(spatialCoords(sfe)))

# Cell boundaries are stored as polygons
cat("\nCell boundary geometry type:\n")
print(class(colGeometry(sfe, "cellSeg")))

# ============================================================================
# 4. QUALITY CONTROL
# ============================================================================
# Standard QC metrics for spatial data

# Calculate QC metrics
sfe <- addPerCellQCMetrics(sfe)

# Look at QC metrics
cat("\nQC metrics available:\n")
print(names(colData(sfe)))

# QC Histograms
p_hist_counts <- ggplot(as.data.frame(colData(sfe)), aes(x = sum)) +
  geom_histogram(bins = 100, fill = okabe_ito[1], color = "black", linewidth = 0.2) +
  labs(x = "Total Counts per Cell", y = "Count", title = "Total Counts Distribution") +
  theme_classic()

p_hist_genes <- ggplot(as.data.frame(colData(sfe)), aes(x = detected)) +
  geom_histogram(bins = 100, fill = okabe_ito[2], color = "black", linewidth = 0.2) +
  labs(x = "Genes Detected per Cell", y = "Count", title = "Genes Detected Distribution") +
  theme_classic()

p_qc_hist <- p_hist_counts + p_hist_genes
ggsave(file.path(fig_dir, "spatial_qc_histograms.png"), p_qc_hist,
       width = 10, height = 4, dpi = 150)

# ============================================================================
# 5. FILTERING
# ============================================================================
# Filter cells with very low counts (likely debris or segmentation errors)

min_counts <- 10
n_before <- ncol(sfe)
sfe <- sfe[, sfe$transcript_counts > min_counts]
cat("\nCells after filtering (>", min_counts, "counts):", ncol(sfe),
    "(removed", n_before - ncol(sfe), ")\n")

# ============================================================================
# 6. REMOVE ISOLATED CELLS (DEBRIS)
# ============================================================================
# Some cells may be isolated debris far from the main tissue.
# We identify these using k-nearest neighbor distances.
# Cells very far from their neighbors are likely broken debris.

cat("\nIdentifying isolated cells (likely debris)...\n")
g <- findKNN(spatialCoords(sfe)[, 1:2], k = 5, BNPARAM = AnnoyParam())
max_dist <- matrixStats::rowMaxs(g$distance)
min_dist <- matrixStats::rowMins(g$distance)

# Cells with very large nearest-neighbor distances are isolated
sfe$main_tissue <- !(max_dist > 100 | min_dist > 60)

n_isolated <- sum(!sfe$main_tissue)
cat("Isolated cells identified:", n_isolated, "\n")

sfe <- sfe[, sfe$main_tissue]
cat("Cells after removing isolated debris:", ncol(sfe), "\n")

# ============================================================================
# 7. NORMALIZATION
# ============================================================================
# IMPORTANT: For spatial data, normalize by CELL AREA, not library size!
# Larger cells naturally capture more transcripts; this is not a technical
# artifact but a biological reality we need to account for.

sfe <- logNormCounts(sfe, size.factor = sfe$cell_area)
cat("\nNormalized counts by cell area\n")

# ============================================================================
# 8. SPATIAL VISUALIZATION - WHOLE TISSUE (CENTROIDS)
# ============================================================================
# For whole-tissue visualization with many cells, use rasterized points
# with cell centroids. This is fast and produces reasonable file sizes.
#
# Following style guide:
# - Use ggrastr::geom_point_rast() for large datasets
# - shape=16 (filled circle), size=0.1, alpha=0.8
# - coord_fixed() to preserve tissue geometry
# - theme_void() for clean spatial plots

df <- as.data.frame(spatialCoords(sfe))
df$sum <- sfe$sum
df$detected <- sfe$detected

# Total counts - whole tissue
# Use inverted scico "lapaz" for quantitative spatial data (histology-like)
p_counts_tissue <- ggplot(df, aes(x = x_centroid, y = y_centroid)) +
  geom_point_rast(aes(color = sum),
                  shape = 16, size = 0.1, alpha = 0.8,
                  raster.dpi = 300) +
  scale_color_gradientn(
    colours = scico(30, palette = "lapaz", direction = -1),
    name = "Total\nCounts"
  ) +
  coord_fixed() +
  theme_void() +
  labs(title = "Total Counts (Centroids)")

ggsave(file.path(fig_dir, "spatial_counts_centroids.png"), p_counts_tissue,
       width = 6, height = 5, dpi = 300)

# Genes detected - whole tissue
p_genes_tissue <- ggplot(df, aes(x = x_centroid, y = y_centroid)) +
  geom_point_rast(aes(color = detected),
                  shape = 16, size = 0.1, alpha = 0.8,
                  raster.dpi = 300) +
  scale_color_gradientn(
    colours = scico(30, palette = "lapaz", direction = -1),
    name = "Genes\nDetected"
  ) +
  coord_fixed() +
  theme_void() +
  labs(title = "Genes Detected (Centroids)")

ggsave(file.path(fig_dir, "spatial_genes_centroids.png"), p_genes_tissue,
       width = 6, height = 5, dpi = 300)

# ============================================================================
# 9. SPATIAL VISUALIZATION - ROI WITH SEGMENTATION (POLYGONS)
# ============================================================================
# For detailed views, use cell segmentation polygons with Voyager's
# plotSpatialFeature(). This shows actual cell boundaries.
#
# Following style guide:
# - Segmentation line color: grey15
# - Line width: 0.2
# - Alpha: 0.8
# - White background

# Define a region of interest (ROI) - adjust based on your tissue
# Here we pick a central region
coords <- spatialCoords(sfe)
center_x <- median(coords[, 1])
center_y <- median(coords[, 2])

# Define bounding box for ROI (500 microns)
bbox_roi <- c(
  xmin = center_x - 250,
  xmax = center_x + 250,
  ymin = center_y - 250,
  ymax = center_y + 250
)

# Count cells in ROI
in_roi <- coords[, 1] >= bbox_roi["xmin"] & coords[, 1] <= bbox_roi["xmax"] &
  coords[, 2] >= bbox_roi["ymin"] & coords[, 2] <= bbox_roi["ymax"]
cat("\nROI contains", sum(in_roi), "cells\n")

# Plot ROI with cell segmentation polygons - QC metric
# Use inverted lapaz for quantitative expression
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
  ) +
  labs(title = "ROI with Cell Segmentation")

ggsave(file.path(fig_dir, "spatial_roi_polygons.png"), p_roi_counts,
       width = 4, height = 3.5, dpi = 300, bg = "white")

# ============================================================================
# 10. GENE EXPRESSION VISUALIZATION
# ============================================================================
# Visualize expression of marker genes spatially
# Use inverted lapaz palette for expression on tissue (histology-like)

# Check which markers are present
epithelial_markers <- c("EPCAM", "KRT8", "KRT18", "PAX8")
epithelial_markers <- epithelial_markers[epithelial_markers %in% rownames(sfe)]

if (length(epithelial_markers) > 0) {
  marker <- epithelial_markers[1]

  # Whole tissue - centroids
  df$expr <- as.numeric(logcounts(sfe)[marker, ])

  p_expr_tissue <- ggplot(df, aes(x = x_centroid, y = y_centroid)) +
    geom_point_rast(aes(color = expr),
                    shape = 16, size = 0.1, alpha = 0.8,
                    raster.dpi = 300) +
    scale_color_gradientn(
      colours = scico(30, palette = "lapaz", direction = -1),
      name = marker
    ) +
    coord_fixed() +
    theme_void() +
    labs(title = paste(marker, "Expression"))

  ggsave(file.path(fig_dir, paste0("spatial_", marker, ".png")),
         p_expr_tissue, width = 6, height = 5, dpi = 300)

  # ROI with segmentation
  p_expr_roi <- plotSpatialFeature(
    sfe,
    features = marker,
    colGeometryName = "cellSeg",
    bbox = bbox_roi,
    linewidth = 0.2,
    color = "grey15"
  ) +
    scale_fill_gradientn(
      colours = scico(30, palette = "lapaz", direction = -1),
      name = marker
    ) +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", color = NA)
    )

  ggsave(file.path(fig_dir, paste0("spatial_roi_", marker, ".png")),
         p_expr_roi, width = 4, height = 3.5, dpi = 300, bg = "white")
}

# Stromal markers
stromal_markers <- c("COL1A1", "DCN", "ACTA2")
stromal_markers <- stromal_markers[stromal_markers %in% rownames(sfe)]

if (length(stromal_markers) > 0) {
  marker <- stromal_markers[1]
  df$expr <- as.numeric(logcounts(sfe)[marker, ])

  p_stromal_tissue <- ggplot(df, aes(x = x_centroid, y = y_centroid)) +
    geom_point_rast(aes(color = expr),
                    shape = 16, size = 0.1, alpha = 0.8,
                    raster.dpi = 300) +
    scale_color_gradientn(
      colours = scico(30, palette = "lapaz", direction = -1),
      name = marker
    ) +
    coord_fixed() +
    theme_void() +
    labs(title = paste(marker, "Expression"))

  ggsave(file.path(fig_dir, paste0("spatial_", marker, ".png")),
         p_stromal_tissue, width = 6, height = 5, dpi = 300)
}

# Immune markers
immune_markers <- c("PTPRC", "CD3E", "CD68")
immune_markers <- immune_markers[immune_markers %in% rownames(sfe)]

if (length(immune_markers) > 0) {
  marker <- immune_markers[1]
  df$expr <- as.numeric(logcounts(sfe)[marker, ])

  p_immune_tissue <- ggplot(df, aes(x = x_centroid, y = y_centroid)) +
    geom_point_rast(aes(color = expr),
                    shape = 16, size = 0.1, alpha = 0.8,
                    raster.dpi = 300) +
    scale_color_gradientn(
      colours = scico(30, palette = "lapaz", direction = -1),
      name = marker
    ) +
    coord_fixed() +
    theme_void() +
    labs(title = paste(marker, "Expression"))

  ggsave(file.path(fig_dir, paste0("spatial_", marker, ".png")),
         p_immune_tissue, width = 6, height = 5, dpi = 300)
}

# ============================================================================
# 11. CUSTOM SPATIAL VISUALIZATION WITH GGPLOT2
# ============================================================================
# For maximum control, extract geometries and use ggplot2 directly
# This is useful when Voyager's built-in functions don't offer enough control

# Subset SFE to ROI for custom plotting
sfe_roi <- sfe[, in_roi]

# Extract cell boundaries as sf object
cell_geom <- colGeometry(sfe_roi, "cellSeg")
cell_data <- as.data.frame(colData(sfe_roi))

# Add data to geometry
cell_sf <- cell_geom
cell_sf$sum <- cell_data$sum
cell_sf$detected <- cell_data$detected

# Custom plot with our style
p_custom_roi <- ggplot(cell_sf) +
  geom_sf(aes(fill = sum), color = "grey15", linewidth = 0.2, alpha = 0.8) +
  scale_fill_gradientn(
    colours = scico(30, palette = "lapaz", direction = -1),
    name = "Total\nCounts"
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(title = "Custom ROI Plot")

ggsave(file.path(fig_dir, "spatial_custom_roi.png"), p_custom_roi,
       width = 4, height = 3.5, dpi = 300, bg = "white")

# ============================================================================
# 12. SAVE PROCESSED DATA
# ============================================================================

saveRDS(sfe, file.path(output_dir, paste0(sample_id, "_processed.rds")))
cat("\nSaved SFE object to:", file.path(output_dir, paste0(sample_id, "_processed.rds")), "\n")

# ============================================================================
# SESSION INFO
# ============================================================================
cat("\n--- Session Info ---\n")
sessionInfo()
