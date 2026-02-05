# =============================================================================
# SPATIAL VISUALIZATION DEMO
# =============================================================================
# Demonstrates lab-standard aesthetics for spatial transcriptomics data
# using Xenium HGSC data processed with the Voyager/SpatialFeatureExperiment
# workflow.
# =============================================================================

# --- Libraries ---
library(Voyager)
library(SFEData)
library(SingleCellExperiment)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(ggplot2)
library(dplyr)
library(scuttle)
library(BiocNeighbors)
library(InSituType)
library(Seurat)
library(ggrastr)
library(shadowtext)
library(scico)
library(viridis)

theme_set(theme_bw())
options(timeout = Inf)

# =============================================================================
# DATA LOADING & PREPROCESSING
# =============================================================================
# Based on standard processing pipeline outlined in Voyager vignette

sfe <- readXenium("~/Projects/style_guide/data/HGSC_SP24_24824/",
                  add_molecules = FALSE)
rownames(sfe) <- rowData(sfe)$Symbol  # Use gene symbols instead of Ensembl IDs

# Remove control probes from matrix
ctrl_probes <- c("NegControlProbe", "NegControlCodeword", "UnassignedCodeword")
to_drop <- grepl(paste(ctrl_probes, collapse = "|"),
                 rownames(sfe), ignore.case = FALSE)
sfe <- sfe[!to_drop, ]

# Quality control: filter low-count cells
sfe <- sfe[, sfe$transcript_counts > 10]

# Remove cells isolated away from tissue (broken debris)
g <- findKNN(spatialCoords(sfe)[, 1:2], k = 5, BNPARAM = AnnoyParam())
max_dist <- matrixStats::rowMaxs(g$distance)
min_dist <- matrixStats::rowMins(g$distance)
sfe$main_tissue <- !(max_dist > 100 | min_dist > 60)
sfe <- sfe[, sfe$main_tissue]

# Normalize based on cell area (NOT library size for spatial data)
sfe <- logNormCounts(sfe, size.factor = sfe$cell_area)

# =============================================================================
# CELL TYPE ANNOTATION
# =============================================================================
# Using InSituType with scRNA-seq reference for automated annotation

ref <- readRDS("~/Projects/style_guide/data/seurat_ovarian_atlas_500ds.rds")

getPseudobulk <- function(seu_obj, celltype_meta) {
  celltype <- factor(seu_obj@meta.data[, celltype_meta])
  names(celltype) <- colnames(seu_obj)
  mat <- seu_obj[["RNA"]]$counts

  mat.summary <- do.call(cbind, lapply(levels(celltype), function(s) {
    cells <- names(celltype)[celltype == s]
    pseudobulk <- rowMeans(mat[, cells])
    return(pseudobulk)
  }))
  colnames(mat.summary) <- levels(celltype)
  return(mat.summary)
}

annotateData <- function(sfe, ref, celltype_meta, genes_remove) {
  message("Getting pseudobulk for reference")
  ref_mat <- getPseudobulk(ref, celltype_meta)
  query_mat <- counts(sfe)
  query_mat <- as(query_mat, "dgCMatrix")
  query_mat <- query_mat[!rownames(query_mat) %in% genes_remove, ]

  message("Annotating spatial data")
  neg <- rep(0, ncol(sfe))
  names(neg) <- colnames(sfe)

  insitutype_res <- insitutypeML(
    x = t(query_mat),
    neg = neg,
    reference_profiles = ref_mat
  )

  sfe$cell_type <- insitutype_res$clust
  return(sfe)
}

sfe <- annotateData(sfe, ref, "celltype_level1", c("CAPS"))


# =============================================================================
# COLOR PALETTE
# =============================================================================
# Modified from lab standards to match this dataset's cell type labels

celltype_palette <- c(
  # Epithelial
  "Epithelial"    = "#E6A141",
  "Mesothelial"   = "#D4A574",
  # Mesenchymal
  "Mesenchymal"   = "#DDD5CA",
  "Smooth_Muscle" = "#D14E6C",
  # Vascular
  "Pericyte_SM"   = "#B87A7A",
  "Endothelial"   = "#7D4E4E",
  # Lymphoid
  "T cell"        = "#87CEFA",
  "NK cell"       = "#56AFC4",
  "B cell"        = "#5665B6",
  "Plasma cell"   = "#8A5DAF",
  # Myeloid
  "Macrophage"    = "#8FBC8F",
  "DC"            = "#2E8B57",
  "Neutrophil"    = "#6B8E23",
  "Mast"          = "#8B9B6B",
  # Other
  "Erythrocyte"   = "#CD5C5C",
  "Other"         = "#A0A0A0"
)

# =============================================================================
# WHOLE TISSUE VISUALIZATION
# =============================================================================
# For standard dot plots, ggplot2 + ggrastr is more flexible than Voyager's
# built-in functions. We can control rasterization, point sizes, etc.

df <- as.data.frame(spatialCoords(sfe))
df$cell_type <- sfe$cell_type

p_tissue <- ggplot(df, aes(x = x_centroid, y = y_centroid)) +
  geom_point_rast(aes(color = cell_type),
                  shape = 16, size = 0.1, alpha = 0.8,
                  raster.dpi = 300) +
  scale_color_manual(values = celltype_palette) +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 3)) +
  coord_fixed() +
  theme_void()
p_tissue

ggsave("figures/spatial_tissue_celltype.png", p_tissue,
       width = 6, height = 5, dpi = 300)


# =============================================================================
# SEGMENTED CELL VISUALIZATION (ROI)
# =============================================================================
# SpatialFeatureExperiment/Voyager handles cell polygons cleanly

# Define a region of interest
bbox_roi <- c(xmin = 7594, xmax = 8307, ymin = -8702, ymax = -8057)

# Cell type with segmentation outlines
p_seg_celltype <- plotSpatialFeature(
  sfe,
  features = "cell_type",
  colGeometryName = "cellSeg",
  bbox = bbox_roi,
  linewidth = 0.2,
  color = "grey15",
  alpha = 0.8
) +
  scale_fill_manual(values = celltype_palette) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )
p_seg_celltype

ggsave("figures/spatial_segmented_celltype.png", p_seg_celltype,
       width = 4, height = 3.5, dpi = 300, bg = "white")

# Expression overlay with scico "lapaz" (inverted) - resembles histology staining
p_seg_expr <- plotSpatialFeature(
  sfe,
  features = "IFIT1",
  colGeometryName = "cellSeg",
  bbox = bbox_roi,
  linewidth = 0.2,
  color = "grey15"
) +
  scale_fill_gradientn(colours = scico(30, palette = "lapaz", direction = -1)) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )
p_seg_expr

ggsave("figures/spatial_segmented_expression.png", p_seg_expr,
       width = 4, height = 3.5, dpi = 300, bg = "white")


# =============================================================================
# UMAP VISUALIZATION
# =============================================================================
# Build Seurat object from SFE data for dimensionality reduction

obj <- CreateSeuratObject(
  counts = counts(sfe),
  meta.data = as.data.frame(colData(sfe))
)
obj[["RNA"]]$data <- assay(sfe, "logcounts")
VariableFeatures(obj) <- rownames(obj)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:15, verbose = FALSE)

# --- UMAP: Cell Type (with legend) ---
df_umap <- data.frame(
  UMAP1 = Embeddings(obj, "umap")[, 1],
  UMAP2 = Embeddings(obj, "umap")[, 2],
  cell_type = obj$cell_type
)

p_umap_celltype <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
  geom_point_rast(aes(color = cell_type),
                  shape = 16, size = 0.1, alpha = 0.8,
                  raster.dpi = 300) +
  scale_color_manual(values = celltype_palette) +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 3)) +
  theme_void()
p_umap_celltype

ggsave("figures/umap_celltype.png", p_umap_celltype,
       width = 5, height = 3.5, dpi = 300)

# --- UMAP: Cell Type (with shadowtext labels, no legend) ---
# shadowtext provides outlined text that remains legible over any background
label_loc <- df_umap %>%
  group_by(cell_type) %>%
  summarize(
    UMAP1 = mean(UMAP1),
    UMAP2 = mean(UMAP2),
    .groups = "drop"
  )

p_umap_labeled <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) +
  geom_point_rast(aes(color = cell_type),
                  shape = 16, size = 0.1, alpha = 0.8,
                  raster.dpi = 300) +
  geom_shadowtext(
    data = label_loc,
    aes(label = cell_type),
    size = 3, color = "black", bg.color = "white"
  ) +
  scale_color_manual(values = celltype_palette) +
  theme_void() +
  theme(legend.position = "none")
p_umap_labeled

ggsave("figures/umap_celltype_labeled.png", p_umap_labeled,
       width = 4.5, height = 4, dpi = 300)

# --- UMAP: Gene Expression ---
# Light grey for zero expression, inverted magma for non-zero values
plotGene <- function(gene, obj) {
  df <- obj@meta.data
  df$UMAP1 <- Embeddings(obj, "umap")[, 1]
  df$UMAP2 <- Embeddings(obj, "umap")[, 2]
  df$Gene <- obj[["RNA"]]$data[gene, ]

  ggplot(df, aes(x = UMAP1, y = UMAP2)) +
    geom_point_rast(aes(color = Gene),
                    shape = 16, size = 0.1, alpha = 0.75,
                    raster.dpi = 300) +
    scale_color_gradientn(
      colours = c("lightgrey", rev(magma(100))),
      name = "log2(TP10k)",
      guide = guide_colorbar(
        ticks.colour = "black",
        frame.colour = "black",
        barwidth = 0.8
      )
    ) +
    ggtitle(gene) +
    theme_void() +
    theme(
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.title = element_text(size = 12, face = "italic")
    )
}

p_umap_expr <- plotGene("PDGFRA", obj)
p_umap_expr

ggsave("figures/umap_expression.png", p_umap_expr,
       width = 4, height = 3.5, dpi = 300)

# =============================================================================
# EXPORT NOTES
# =============================================================================
# Keep dimensions compact - these visualizations work well in publications
# and presentations at small sizes. Avoid large whitespace margins.
