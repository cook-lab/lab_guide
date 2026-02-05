# Data Visualization Skill

Generate publication-ready visualizations following Cook Lab style standards.

## When to Use

Use this skill when:
- Creating any data visualization in R/ggplot2
- Plotting single-cell or spatial transcriptomics data
- Generating figures for publications or presentations
- The user asks for plots, figures, UMAPs, heatmaps, etc.

## Critical Rules

**These rules are non-negotiable. Violating them produces incorrect output.**

### Color Palettes

| Data Type | Palette | Code |
|-----------|---------|------|
| Continuous (expression, scores) | `magma` or `viridis` | `scale_color_viridis_c(option = "magma")` |
| Diverging (log2FC, z-scores) | `scico::vik` | `scale_fill_scico(palette = "vik", midpoint = 0)` |
| Categorical ≤8 | Okabe-Ito | `scale_color_manual(values = okabe_ito)` |
| Categorical >8 | Kelly or celltype_palette | `scale_color_manual(values = kelly(n))` |
| Cell types | `celltype_palette` | `scale_color_manual(values = celltype_palette)` |
| Spatial expression | `scico::lapaz` inverted | `scale_fill_scico(palette = "lapaz", direction = -1)` |
| UMAP expression | grey → magma | `scale_color_gradientn(colours = c("lightgrey", rev(magma(100))))` |

**NEVER use:**
- Rainbow/jet palettes
- Red-green combinations
- Diverging palettes for non-centered data (expression levels are NOT centered)

### Plot Type Selection

| Purpose | Use | NEVER Use |
|---------|-----|-----------|
| Distribution comparison | Strip plot, box plot, violin | Bar plot with error bars |
| Proportions/composition | Stacked bar | Pie chart |
| Small n (<50) | Strip plot with points | Violin, histogram |
| Large n (>200) | Violin, density | Individual points alone |

### Theme

**Always use `theme_lab()` or equivalent clean theme:**

```r
theme_lab <- function(base_family = "") {
  theme_classic(base_family = base_family) %+replace%
    theme(
      # Explicit font sizes for readable, consistent appearance
      text = element_text(color = "black"),
      plot.title = element_text(size = 14, margin = margin(b = 8)),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      strip.text = element_text(size = 11),
      # Clean appearance
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    )
}
```

**Key points:**
- Explicit sizes: title=14, axis.title=12, axis.text=10, legend=9-10
- No gridlines, no panel border (only axis lines)
- For categorical x-axis: add `theme(axis.text.x = element_text(angle = 45, hjust = 1))`

---

## Palettes Reference

### Okabe-Ito (categorical ≤8)

```r
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)
```

### Lab Cell Type Palette

```r
celltype_palette <- c(
  # Epithelial
  "Epi_Secretory"       = "#E6A141",
  "Epi_Secretory_Polar" = "#9A7D55",
  "Epi_Ciliated"        = "#E07850",
  "Mesothelial"         = "#D4A574",
  # Mesenchymal
  "Fibroblast"          = "#DDD5CA",
  "Smooth_Muscle"       = "#D14E6C",
  # Vascular
  "Pericyte"            = "#B87A7A",
  "Endothelial"         = "#7D4E4E",
  # Lymphoid
  "T_cell"              = "#87CEFA",
  "NK_cell"             = "#56AFC4",
  "B_cell"              = "#5665B6",
  "Plasma_cell"         = "#8A5DAF",
  # Myeloid
  "Macrophage"          = "#8FBC8F",
  "DC"                  = "#2E8B57",
  "Neutrophil"          = "#6B8E23",
  "Mast"                = "#8B9B6B",
  # Other
  "Erythrocyte"         = "#CD5C5C",
  "Other"               = "#A0A0A0"
)
```

---

## Code Templates

### UMAP (Cell Types)

```r
library(ggplot2)
library(ggrastr)

df <- data.frame(
  UMAP1 = Embeddings(obj, "umap")[, 1],
  UMAP2 = Embeddings(obj, "umap")[, 2],
  cell_type = obj$cell_type
)

p <- ggplot(df, aes(x = UMAP1, y = UMAP2)) +
  geom_point_rast(aes(color = cell_type),
                  shape = 16, size = 0.1, alpha = 0.8,
                  raster.dpi = 300) +
  scale_color_manual(values = celltype_palette) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_void()

ggsave("umap_celltype.png", p, width = 5, height = 3.5, dpi = 300)
```

**Required parameters for UMAP/spatial:**
- `shape = 16` (filled circle)
- `size = 0.1` (very small to prevent overplotting)
- `alpha = 0.8`
- `raster.dpi = 300`
- Legend: `override.aes = list(size = 2)`

### UMAP with Direct Labels

```r
library(shadowtext)
library(dplyr)

label_loc <- df %>%
  group_by(cell_type) %>%
  summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2), .groups = "drop")

p <- ggplot(df, aes(x = UMAP1, y = UMAP2)) +
  geom_point_rast(aes(color = cell_type),
                  shape = 16, size = 0.1, alpha = 0.8, raster.dpi = 300) +
  geom_shadowtext(data = label_loc, aes(label = cell_type),
                  size = 3, color = "black", bg.color = "white") +
  scale_color_manual(values = celltype_palette) +
  theme_void() +
  theme(legend.position = "none")
```

### Feature Plot (Gene Expression on UMAP)

```r
library(viridis)

df$expression <- obj[["RNA"]]$data["GENE", ]

p <- ggplot(df, aes(x = UMAP1, y = UMAP2)) +
  geom_point_rast(aes(color = expression),
                  shape = 16, size = 0.1, alpha = 0.75, raster.dpi = 300) +
  scale_color_gradientn(
    colours = c("lightgrey", rev(magma(100))),
    name = "Expression"
  ) +
  ggtitle("GENE") +
  theme_void() +
  theme(plot.title = element_text(face = "italic"))
```

**Critical:** Use `c("lightgrey", rev(magma(100)))` so non-expressing cells are grey.

### Spatial Tissue Plot

```r
df <- as.data.frame(spatialCoords(sfe))
df$cell_type <- sfe$cell_type

p <- ggplot(df, aes(x = x_centroid, y = y_centroid)) +
  geom_point_rast(aes(color = cell_type),
                  shape = 16, size = 0.1, alpha = 0.8, raster.dpi = 300) +
  scale_color_manual(values = celltype_palette) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  coord_fixed() +
  theme_void()
```

**Critical:** Always use `coord_fixed()` for spatial data.

### Spatial Expression (Segmented Cells)

```r
library(Voyager)
library(scico)

plotSpatialFeature(sfe, features = "GENE", colGeometryName = "cellSeg",
                   bbox = bbox_roi, linewidth = 0.2, color = "grey15") +
  scale_fill_scico(palette = "lapaz", direction = -1) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", color = NA))
```

### Violin + Box Plot

```r
ggplot(data, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = okabe_ito) +
  theme_lab() +
  theme(legend.position = "none")
```

### Strip Plot (Small n)

```r
ggplot(data, aes(x = group, y = value, color = group)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5,
               color = "black", linewidth = 0.4) +
  scale_color_manual(values = okabe_ito) +
  theme_lab() +
  theme(legend.position = "none")
```

### Stacked Bar (Proportions)

```r
ggplot(prop_data, aes(x = condition, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = celltype_palette) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion") +
  theme_lab() +
  theme(panel.grid.major.x = element_blank())
```

### Dot Plot (Markers)

```r
ggplot(dot_data, aes(x = gene, y = cell_type)) +
  geom_point(aes(size = pct_exp, fill = avg_exp), shape = 21, color = "black") +
  scale_size_continuous(range = c(0.5, 5), name = "% Expressing") +
  scale_fill_viridis_c(option = "magma", name = "Expression") +
  theme_lab() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

### Heatmap (ComplexHeatmap)

```r
library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(seq(-2, 2, length.out = 100),
                       colorRampPalette(brewer.pal(9, "RdBu"))(100))

Heatmap(mat, name = "Z-score", col = col_fun,
        cluster_rows = TRUE, cluster_columns = TRUE,
        clustering_method_rows = "ward.D2",
        column_names_rot = 45,
        rect_gp = gpar(col = "black", lwd = 1),
        width = ncol(mat) * unit(5, "mm"),
        height = nrow(mat) * unit(5, "mm"),
        use_raster = FALSE)
```

**Critical for heatmaps:**
- Always set explicit `width` and `height` for square cells
- Always cluster or order meaningfully
- Cap color range to handle outliers
- Use `rect_gp = gpar(col = "black", lwd = 1)` for cell borders

---

## Export Dimensions

### Critical: Optimize Whitespace

**The #1 visualization mistake: exporting at dimensions that are too large.**

When a 10" × 8" figure gets scaled to fit a 3.5" column or PowerPoint slide, text becomes unreadable. **Export at dimensions close to your final display size.**

| Target | Export Size |
|--------|-------------|
| Publication single column | 4-5" wide |
| Half a PowerPoint slide | 5-6" wide |
| Full slide figure | 9-10" wide |

**Test:** Open exported file at 100% zoom. All text should be readable.

**If scaling down >50% to fit destination, you're exporting too large.**

| Plot Type | Width | Height |
|-----------|-------|--------|
| Whole tissue | 6-7" | 5-6" |
| ROI segmented | 5" | 4" |
| UMAP with legend | 6" | 4" |
| UMAP labeled | 5" | 4.5" |
| Feature plot | 5" | 4" |
| Standard working size | 7" | 5" |

```r
# Working/exploration
ggsave("fig_draft.png", p, width = 7, height = 5, dpi = 150)

# Publication (export larger, will scale to column width)
ggsave("fig.pdf", p, width = 7, height = 5)

# Presentation
ggsave("fig_slides.png", p, width = 10, height = 5.6, dpi = 300)
```

---

## Checklist Before Generating a Plot

1. [ ] **Export size appropriate?** Close to final display size—not too large
2. [ ] **Palette appropriate?** Sequential for expression, diverging for centered data, qualitative for categories
3. [ ] **No forbidden palettes?** No rainbow, no red-green
3. [ ] **Plot type appropriate?** No bar plots for distributions, no pie charts
4. [ ] **Theme clean?** No gridlines, no bold, no panel border
5. [ ] **UMAP/spatial params?** shape=16, size=0.1, alpha=0.8, raster.dpi=300
6. [ ] **Legend readable?** override.aes for small points
7. [ ] **Spatial aspect ratio?** coord_fixed()
8. [ ] **Heatmap square?** Explicit width/height with unit()

---

## Full Reference

See `guides/visualization_style_guide.md` for complete documentation including:
- NPG figure specifications
- All palette hex codes with visual examples
- Detailed anti-pattern explanations
- Additional plot type examples
