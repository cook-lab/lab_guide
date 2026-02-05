# scRNA-seq and Spatial Transcriptomics Skill

Guide single-cell and spatial transcriptomics analysis following Cook Lab best practices.

## When to Use

Use this skill when:
- Processing scRNA-seq data (Seurat or Scanpy)
- Processing spatial transcriptomics data (Xenium, Visium, etc.)
- Performing QC, normalization, clustering, or annotation
- Transferring annotations between modalities
- The user asks about single-cell or spatial analysis workflows

## Reference

Full tutorial: `tutorials/scrna_spatial_workflows.md`
Scripts: `tutorials/scripts/`

---

## Core Philosophy

### Default to Fewer Filters
You can always see the impact on downstream analysis and circle back if needed. **Over-filtering removes real biology.**

### Examine YOUR Data
Never blindly apply thresholds from other datasets. Always visualize distributions and set thresholds based on the specific data.

### Readability Over Efficiency
Code should be understandable and reproducible. Structure for interactive execution in IDE (RStudio, Jupyter).

---

## scRNA-seq Processing

### Recommended Workflow Order

```
1. Load data (filtered + raw matrices)
2. Ambient RNA removal (SoupX)
3. Create Seurat/AnnData object
4. Calculate QC metrics
5. Visualize QC distributions (histograms!)
6. Doublet detection (scDblFinder/Scrublet)
7. Filter cells
8. Normalize (SCTransform + log-norm)
9. Dimensionality reduction (PCA → UMAP)
10. Clustering (start low resolution)
11. Cell type annotation
12. Save and export
```

### QC Metrics

Always calculate:
```r
seu[["percent_mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent_ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
```

**Visualize with histograms** (better than violins for threshold setting):
```r
ggplot(seu@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100) +
  theme_lab()
```

### Doublet Detection

**Use explicit doublet detection, NOT upper nFeature/nCount thresholds.**

Upper-limit thresholds are not specific for doublets (ground truth studies show this). We rely on:

**Seurat/R:**
```r
library(scDblFinder)
sce <- as.SingleCellExperiment(seu)
sce <- scDblFinder(sce)
seu$scDblFinder_class <- sce$scDblFinder.class
```

**Scanpy/Python:**
```python
sc.external.pp.scrublet(adata, expected_doublet_rate=0.08)

# Scrublet's auto-threshold can fail; use percentile fallback
if adata.obs["predicted_doublet"].sum() / len(adata.obs) < 0.04:
    threshold = np.percentile(adata.obs["doublet_score"], 92)
    adata.obs["predicted_doublet"] = adata.obs["doublet_score"] > threshold
```

### Filtering

Filter based on:
1. **percent_mt** — remove dying/stressed cells (threshold based on YOUR distribution)
2. **Doublets** — remove cells classified as doublets

```r
seu <- subset(seu,
  subset = percent_mt < 20 &
    scDblFinder_class == "singlet"
)
```

**Do NOT apply upper-limit thresholds on nCount/nFeature** — we already run explicit doublet detection.

### Normalization

**Run both normalizations:**

```r
# SCTransform for embedding/clustering
seu <- SCTransform(seu, method = "glmGamPoi", vst.flavor = "v2")

# Log-normalization for visualization and DE
seu <- NormalizeData(seu, assay = "RNA")
```

- **SCTransform:** Use for PCA, UMAP, clustering, integration
- **Log-normalized:** Use for FeaturePlots, violin plots, differential expression

### Clustering Resolution

**More clusters is NOT always better.**

- High resolution increases risk of technical/sample-specific patterns
- More clusters = more work to characterize each one
- **Start at resolution 0.2-0.3** and increase only if you need finer granularity

```r
# Compare a few resolutions
for (res in c(0.2, 0.3, 0.5)) {
  seu <- FindClusters(seu, resolution = res)
}
# Visualize and choose the one matching your biological question
```

### Cell Type Annotation

**Options:**
1. **Manual** — DotPlots + FeaturePlots with known markers
2. **CellAssign** — Probabilistic annotation with marker matrix
3. **SingleR** — Reference-based annotation
4. **Celltypist** — Pre-trained ML models

```r
# DotPlot for marker visualization (use RNA assay!)
DotPlot(seu, features = markers, assay = "RNA") +
  scale_color_gradientn(colors = c("lightgrey", rev(viridis::magma(100))))
```

---

## Spatial Transcriptomics

### Recommended Stack

**R:** SpatialFeatureExperiment + Voyager
- Stores centroids AND segmentation polygons
- Rich spatial statistics via sf integration
- Best visualization support

### Loading Xenium Data

```r
library(SpatialFeatureExperiment)
library(Voyager)

sfe <- readXenium(xenium_dir, add_molecules = FALSE)

# Use gene symbols (not Ensembl IDs)
rownames(sfe) <- rowData(sfe)$Symbol

# Remove control probes
ctrl_probes <- c("NegControlProbe", "NegControlCodeword", "UnassignedCodeword")
sfe <- sfe[!grepl(paste(ctrl_probes, collapse = "|"), rownames(sfe)), ]
```

### Spatial QC and Filtering

```r
# Remove low-count cells
sfe <- sfe[, sfe$transcript_counts > 10]

# Remove isolated debris using k-nearest neighbor distances
library(BiocNeighbors)
g <- findKNN(spatialCoords(sfe), k = 5, BNPARAM = AnnoyParam())
max_dist <- matrixStats::rowMaxs(g$distance)
min_dist <- matrixStats::rowMins(g$distance)
sfe$main_tissue <- !(max_dist > 100 | min_dist > 60)
sfe <- sfe[, sfe$main_tissue]
```

### Normalization

**Normalize by cell area, not library size:**
```r
sfe <- logNormCounts(sfe, size.factor = sfe$cell_area)
```

Larger cells naturally capture more transcripts.

### Geometry Types

| Geometry | Type | Use Case |
|----------|------|----------|
| Centroids | POINT | Fast whole-tissue plots |
| Cell segmentation | POLYGON | Detailed ROI views |
| Nucleus segmentation | POLYGON | Nucleus-specific analysis |

```r
# Check available geometries
colGeometryNames(sfe)  # "centroids", "cellSeg", "nucSeg"
```

### Visualization

**Whole tissue (centroids):**
```r
library(ggrastr)

df <- as.data.frame(spatialCoords(sfe))
df$value <- sfe$sum

ggplot(df, aes(x = x_centroid, y = y_centroid)) +
  geom_point_rast(aes(color = value),
                  shape = 16, size = 0.1, alpha = 0.8, raster.dpi = 300) +
  scale_color_scico(palette = "lapaz", direction = -1) +
  coord_fixed() +
  theme_void()
```

**ROI with cell boundaries:**
```r
bbox_roi <- c(xmin = x - 250, xmax = x + 250, ymin = y - 250, ymax = y + 250)

plotSpatialFeature(sfe, features = "sum",
                   colGeometryName = "cellSeg",
                   bbox = bbox_roi,
                   linewidth = 0.2, color = "grey15") +
  scale_fill_scico(palette = "lapaz", direction = -1)
```

---

## Reference-Based Annotation (SingleR)

Transfer annotations from scRNA-seq to spatial data:

```r
library(SingleR)

# Prepare reference (annotated scRNA-seq)
sce_ref <- as.SingleCellExperiment(seu, assay = "RNA")
ref_labels <- seu$cell_type

# Find common genes (limited by spatial panel)
common_genes <- intersect(rownames(sce_ref), rownames(sfe))

# Run SingleR
singler_results <- SingleR(
  test = sfe[common_genes, ],
  ref = sce_ref[common_genes, ],
  labels = ref_labels,
  de.method = "classic"
)

sfe$singler_label <- singler_results$labels
sfe$singler_delta <- singler_results$delta.next
```

**Check annotation quality:**
- Score heatmap: shows correlation with each cell type
- Delta distribution: higher = more confident (low delta = ambiguous)

---

## Data Export for Interoperability

### Seurat → Cell Ranger Format (for Scanpy)

```r
counts_matrix <- LayerData(seu, assay = "RNA", layer = "counts")

# Matrix
writeMM(counts_matrix, "matrix.mtx")
system("gzip -f matrix.mtx")

# Features
features_df <- data.frame(
  gene_id = rownames(seu),
  gene_name = rownames(seu),
  feature_type = "Gene Expression"
)
write.table(features_df, "features.tsv", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
system("gzip -f features.tsv")

# Barcodes
write.table(colnames(seu), "barcodes.tsv",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
system("gzip -f barcodes.tsv")

# Metadata
write.csv(seu@meta.data, "cell_metadata.csv")
```

### Load in Scanpy

```python
adata = sc.read_10x_mtx("exported_data/")
metadata = pd.read_csv("exported_data/cell_metadata.csv", index_col=0)
adata.obs = metadata
```

---

## Quick Reference

### Seurat Workflow
```r
# After loading data
seu <- CreateSeuratObject(counts, min.cells = 3, min.features = 200)
seu[["percent_mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
# ... doublet detection ...
seu <- subset(seu, subset = percent_mt < 20 & scDblFinder_class == "singlet")
seu <- SCTransform(seu, method = "glmGamPoi", vst.flavor = "v2")
seu <- NormalizeData(seu, assay = "RNA")
seu <- RunPCA(seu)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.2)
seu <- RunUMAP(seu, dims = 1:30)
```

### Scanpy Workflow
```python
# After loading data
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
sc.external.pp.scrublet(adata)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
adata = adata[~adata.obs["predicted_doublet"]].copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, resolution=0.2)
sc.tl.umap(adata)
```

### Key Thresholds (Starting Points)

| Parameter | Typical Range | Notes |
|-----------|---------------|-------|
| min.features | 200-500 | Lower for spatial (fewer genes) |
| percent_mt | 10-25% | Depends on tissue type |
| Clustering resolution | 0.2-0.5 | Start low |
| PCA dimensions | 20-50 | Check elbow plot |
| Doublet rate | 5-15% | Depends on loading density |

---

## Environment

**R:** Local installation with Seurat, SpatialFeatureExperiment, Voyager

**Python:**
- `mamba activate scverse` for standard scanpy workflows
- `mamba activate claude_code` if new tools must be installed
