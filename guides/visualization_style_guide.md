# Data Visualization Style Guide
## Cook Lab — Single-Cell & Spatial Genomics

*Version 0.6*

---

## Table of Contents

1. [Philosophy & Principles](#philosophy--principles)
2. [Nature Publishing Group Guidelines](#nature-publishing-group-guidelines)
3. [Color Selection](#color-selection)
4. [Typography & Legibility](#typography--legibility)
5. [ggplot2 Theme & Defaults](#ggplot2-theme--defaults)
6. [Plot Types & When to Use Them](#plot-types--when-to-use-them)
7. [Common Anti-Patterns to Avoid](#common-anti-patterns-to-avoid)
8. [Domain-Specific Visualizations](#domain-specific-visualizations)

---

## Philosophy & Principles

This guide draws heavily from [Martin Krzywinski's "Points of View"](https://mk.bcgsc.ca/pointsofview/) column in *Nature Methods* and the principles articulated in [Friends Don't Let Friends Make Bad Graphs](https://github.com/cxli233/FriendsDontLetFriends).

### Core Tenets

1. **Clarity over decoration.** Every visual element should serve a purpose. Remove chartjunk, excessive gridlines, and redundant labels.

2. **Show the data.** Prefer plot types that reveal the underlying distribution rather than obscuring it behind summary statistics.

3. **Accessibility first.** Use colorblind-friendly palettes. Design for grayscale printing when possible.

4. **Publication & presentation ready.** Text must be legible at final reproduction size. Optimize whitespace—neither cramped nor wasteful.

5. **Aesthetics matter.** Beautiful figures communicate professionalism and care. Strive for visual harmony.

6. **Consistency across the lab.** Shared palettes and themes build a recognizable visual identity and reduce cognitive load for collaborators and reviewers.

---

## Nature Publishing Group Guidelines

When preparing figures for Nature journals, follow these specifications from the [Nature Research Figure Guide](https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/).

### Font Requirements

| Element | Specification |
|---------|---------------|
| **Font type** | Sans-serif only (Arial or Helvetica) |
| **Body text** | 5–7 pt (minimum to maximum) |
| **Panel labels** | 8 pt, bold, lowercase (a, b, c...), upright (not italic) |
| **Greek characters** | Symbol font |
| **Critical rule** | Never outline text — keep editable |

### Figure Dimensions

| Width | Measurement | Use |
|-------|-------------|-----|
| Single column | 88 mm (3.46") | Default for simple figures |
| 1.5 column | 120 mm (4.72") | Medium complexity |
| Double column | 180 mm (7.09") | Complex multi-panel figures |

### Color & Image Quality

- **Color space:** RGB (not CMYK)
- **Resolution:** Minimum 300 dpi for photographs; 450 dpi optimal
- **Accessibility:** Must use colorblind-friendly palettes
- **Text color:** Avoid colored text; use keylines and legends instead

### Export Format

- **Preferred:** PDF or EPS with vector artwork
- **All text, scale bars, boxes:** Vector (not rasterized)
- **Embedded images:** Minimum 450 dpi

### Design Rules

| Do | Don't |
|----|-------|
| Include axis lines and tick marks | Use gridlines or drop shadows |
| Label axes with units in parentheses | Use patterned fills |
| Keep scale bars on separate layers | Overlap text with data |
| Use consistent fonts throughout | Mix font families |

---

## Color Selection

Color is arguably the most important and most frequently misused element in scientific visualization.

### General Rules

| Rule | Rationale |
|------|-----------|
| **Colorblind-friendly** | ~8% of males have color vision deficiency. Red-green combinations are particularly problematic. |
| **Perceptually uniform** | Equal steps in data should produce equal steps in perceived color. Rainbow palettes fail this. |
| **Print-safe** | Many publications still print in grayscale. Luminance should encode key information. |
| **Meaningful defaults** | Diverging palettes for data with a meaningful center; sequential for unidirectional data. |

### Continuous Data

For continuous variables (gene expression, scores, gradients), use **perceptually uniform sequential palettes**.

#### Recommended Palettes

**Viridis family** (available in base R 4.0+ and `viridisLite`):
- `viridis` — General purpose, excellent contrast
- `magma` — Warm tones, good for heatmaps
- `plasma` — High contrast, vivid
- `inferno` — Similar to magma with more yellow
- `cividis` — Optimized for deuteranopia

```r
library(ggplot2)
library(viridisLite)

ggplot(data, aes(x, y, fill = expression)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma")
```

**scico palettes** (via the `scico` package):
- `batlow` — Perceptually uniform rainbow alternative
- `roma` — Diverging, centered at white
- `vik` — Diverging blue-white-red
- `berlin` — Diverging blue-white-orange
- `lajolla` — Sequential warm tones
- `lapaz` — Sequential, good for spatial expression overlays

```r
library(scico)

# Sequential
ggplot(data, aes(x, y, fill = value)) +
  geom_tile() +
  scale_fill_scico(palette = "batlow")

# Diverging (e.g., log2 fold change centered at 0)
ggplot(data, aes(x, y, fill = log2fc)) +
  geom_tile() +
  scale_fill_scico(palette = "vik", midpoint = 0)
```

#### When to Use Diverging vs Sequential

| Data Type | Palette Type | Example Palettes |
|-----------|--------------|------------------|
| Expression levels, counts, scores | Sequential | `viridis`, `magma`, `batlow` |
| Log2 fold change, z-scores, centered metrics | Diverging | `vik`, `roma`, `berlin` |
| Correlation coefficients | Diverging (centered at 0) | `roma`, `vik` |

**Anti-pattern:** Using a diverging palette (red-white-blue) for expression data that has no meaningful center point assigns false significance to arbitrary values.

### Discrete/Categorical Data

For categorical variables (cell types, clusters, conditions), use **qualitative palettes** designed for maximum distinguishability.

#### Palette Hierarchy (in order of preference)

1. **Okabe-Ito** (8 colors) — *Lab default for small categorical sets*
   - Designed specifically for colorblind accessibility
   - Clean, professional appearance

```r
# Base R (4.0+)
palette.colors(palette = "Okabe-Ito")

# Via ggokabeito package
library(ggokabeito)
ggplot(data, aes(x, y, color = condition)) +
  geom_point() +
  scale_color_okabe_ito()

# Manual specification
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)
```

2. **RColorBrewer qualitative palettes** (≤12 colors)
   - `Set1`, `Set2`, `Dark2`, `Paired` are safe choices
   - Avoid `Set3` for serious publications (too pastel)

```r
library(RColorBrewer)
ggplot(data, aes(x, y, fill = cluster)) +
  geom_bar() +
  scale_fill_brewer(palette = "Set2")
```

3. **Kelly's 22 Colors** (large categorical sets) — *Lab default for >9 categories*
   - Designed for maximum contrast between sequential colors
   - First 9 optimized for colorblind viewers

```r
library(pals)
kelly_colors <- kelly(22)

ggplot(data, aes(x, y, color = cell_type)) +
  geom_point() +
  scale_color_manual(values = kelly(n_distinct(data$cell_type)))
```

#### Number of Categories Decision Tree

```
How many categories?
├── ≤8  -> Okabe-Ito
├── 9-12 -> RColorBrewer (Set2, Paired)
└── >12 -> Kelly's 22 OR lab cell type palette
```

---

### Lab Cell Type Palette

For spatial transcriptomics and single-cell visualizations, we use a custom palette optimized for:

1. **Visual hierarchy** — Background cells (fibroblasts, mesenchymal) are neutral; cells of interest pop
2. **Lineage coherence** — Related cell types share color families
3. **Tissue context** — Colors that look natural when mapped to tissue coordinates

#### Design Principles

| Cell Category | Color Strategy | Rationale |
|---------------|----------------|-----------|
| **Fibroblast/Mesenchymal** | Warm nude/beige | Should recede as the "canvas"; warm tones avoid clinical feel |
| **Epithelial** | Orange-brown family | Stand out from stroma; polarized secretory is darker brown |
| **Lymphoid (T, NK, B, Plasma)** | Blue-purple gradient | Family resemblance; cool tones contrast with warm epithelium |
| **Myeloid (Mac, DC, Neutrophil)** | Green family | Distinct from lymphoid; natural tissue association |
| **Mast cells** | Sage/olive | Myeloid-adjacent but muted |
| **Endothelial** | Muted maroon | Blood vessel association |
| **Pericyte** | Muted rose | Vascular support; contrasts with maroon endothelial |
| **Erythrocyte** | Red/coral | Blood association |
| **Other/Unknown** | Gray | Neutral, doesn't draw attention |

#### The Palette (v1.2)

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

scale_color_celltype <- function(...) scale_color_manual(values = celltype_palette, ...)
scale_fill_celltype <- function(...) scale_fill_manual(values = celltype_palette, ...)
```

---

## Typography & Legibility

**The cardinal rule: Text must be readable at final reproduction size.**

### Minimum Sizes

| Context | Minimum Font Size |
|---------|-------------------|
| Nature journals body text | 5–7 pt |
| Nature panel labels | 8 pt bold |
| General figure panels | 6–8 pt after scaling |
| Axis titles | 7–9 pt |
| Poster figures | Scale up 2–3x |
| Slides | 14 pt minimum for any text |

### Font Recommendations

- **Sans-serif for figures:** Arial, Helvetica
- **Consistency:** Use the same font family throughout a manuscript
- **Avoid:** Times New Roman in figures (serifs become muddy at small sizes)

---

## ggplot2 Theme & Defaults

### Lab Standard Theme

```r
theme_lab <- function(base_family = "") {
  theme_classic(base_family = base_family) %+replace%
    theme(
      # Text elements - explicit sizes for consistent appearance
      text = element_text(color = "black"),
      plot.title = element_text(size = 14, margin = margin(b = 8)),
      plot.subtitle = element_text(size = 10, color = "gray40", margin = margin(b = 8)),
      axis.title = element_text(size = 12),
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6), angle = 90),
      axis.text = element_text(size = 10, color = "black"),

      # Axis lines (x and y only, no border)
      axis.line = element_line(color = "black", linewidth = 0.5),

      # Legend
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.position = "right",

      # Panel (no grid, no border)
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),

      # Strip (facets)
      strip.background = element_blank(),
      strip.text = element_text(size = 11, margin = margin(b = 4)),

      # Margins
      plot.margin = margin(8, 8, 8, 8)
    )
}

# Set as session default
theme_set(theme_lab())
```

**Key design decisions:**
- Based on `theme_classic()` — no gridlines, no panel border
- Explicit font sizes: title (14), axis titles (12), axis text (10), legend (9-10)
- Design at readable sizes, export at dimensions that scale appropriately for publication
- For categorical x-axis labels: add `theme(axis.text.x = element_text(angle = 45, hjust = 1))`

### Standard Figure Dimensions

Design at readable sizes, then export at dimensions appropriate for the context:

```r
# Working/exploration - comfortable viewing size
ggsave("fig_draft.png", p, width = 7, height = 5, dpi = 150)

# Publication - export larger, will be scaled down to column width
# Single column target: 88mm (~3.5"), but export at 2x for quality
ggsave("fig_publication.pdf", p, width = 7, height = 5)

# Presentations (16:9 aspect)
ggsave("fig_slides.png", p, width = 10, height = 5.6, dpi = 300)
```

**Note:** NPG column widths are 88mm (single), 120mm (1.5), 180mm (double). Design at comfortable sizes with readable fonts; the figure will scale appropriately when placed in the manuscript.

### Whitespace Optimization

**The most common visualization mistake: exporting at dimensions that are too large.**

When a 10" × 8" figure gets scaled down to fit a 3.5" column or a PowerPoint slide, the text becomes unreadable and the data gets lost in excessive margins. The solution is to **export at dimensions close to your final display size**.

#### The Problem

```r
# BAD: Exporting too large
ggsave("fig.png", p, width = 12, height = 10)
# When this gets shrunk to fit a slide, 10pt text becomes 4pt
```

#### The Solution

1. **Know your target size.** Where will this figure appear?
   - Publication single column: ~3.5" wide
   - PowerPoint slide (with text): ~5-6" wide
   - Full slide figure: ~10" wide

2. **Export at or near target size.** Don't rely on scaling down.
   ```r
   # For a figure that will occupy half a slide
   ggsave("fig.png", p, width = 5, height = 4, dpi = 300)

   # For a publication figure
   ggsave("fig.pdf", p, width = 5, height = 4)
   ```

3. **Test the output.** Open the exported file and view at actual size (100% zoom). Can you read all text? Is the data visible?

#### Whitespace Checklist

- [ ] **Plot margins:** Are they minimal? Large margins waste space when scaled.
- [ ] **Legend position:** Does it fit naturally or create awkward whitespace?
- [ ] **Aspect ratio:** Does it match the target space? A 1:1 figure in a wide slot wastes half the area.
- [ ] **Text legibility:** At final display size, is all text ≥6pt?
- [ ] **Data density:** Is the actual data using most of the plot area?

#### Common Fixes

```r
# Reduce plot margins
theme(plot.margin = margin(4, 4, 4, 4))

# Move legend to reduce width
theme(legend.position = "bottom")
guides(color = guide_legend(nrow = 2))

# Remove unnecessary legend
theme(legend.position = "none")

# Adjust aspect ratio to fit target
coord_fixed(ratio = 0.8)  # or let it be flexible
```

**Rule of thumb:** If you're consistently scaling figures down by more than 50% to fit their destination, you're exporting too large.

---

## Plot Types & When to Use Them

### Decision Framework

```
What are you showing?
│
├── Distribution of a continuous variable
│   ├── n < 50 -> Strip plot / jittered points
│   ├── n = 50-200 -> Box plot + points OR violin + points
│   └── n > 200 -> Violin plot OR density ridges
│
├── Comparison of means/medians
│   └── NEVER bar plot -> Use strip/box/violin
│
├── Relationship between two continuous variables
│   ├── Few points -> Scatter plot
│   └── Many points -> Hex bin OR contour density
│
├── Composition / Proportions
│   ├── Single group -> Stacked bar (one bar)
│   ├── Multiple groups -> Stacked bar OR faceted
│   └── NEVER pie chart
│
├── Expression matrix / Heatmap
│   ├── Always cluster or order meaningfully
│   ├── Check for outliers before setting color scale
│   └── Consider capping at 95th/99th percentile
│
└── Spatial data
    └── See Domain-Specific section
```

### Distribution Comparison

**Never use bar plots with error bars** — they hide the underlying distribution.

```r
# Strip plot with mean crossbar
ggplot(data, aes(x = group, y = value, color = group)) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5,
               color = "black", linewidth = 0.4) +
  scale_color_okabe_ito() +
  theme_lab()
```

### Violin + Box Plot (Large n)

```r
ggplot(data, aes(x = cell_type, y = expression, fill = cell_type)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = celltype_palette) +
  theme_lab() +
  theme(legend.position = "none")
```

### Cell Type Composition (Stacked Bar)

```r
ggplot(prop_data, aes(x = condition, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = celltype_palette, name = "Cell Type") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion") +
  theme_lab() +
  theme(panel.grid.major.x = element_blank())
```

### Marker Gene Dot Plots

```r
ggplot(dot_data, aes(x = gene, y = cell_type)) +
  geom_point(aes(size = pct_exp, fill = avg_exp), shape = 21, color = "black") +
  scale_size_continuous(range = c(0.5, 5), name = "% Expressing",
                        labels = scales::percent) +
  scale_fill_viridis_c(option = "magma", name = "Avg Expression") +
  theme_lab() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

### Heatmaps (ComplexHeatmap)

```r
library(ComplexHeatmap)
library(circlize)

# Color scale with capped range
col_fun <- colorRamp2(
  seq(-2, 2, length.out = 100),
  colorRampPalette(brewer.pal(9, "RdBu"))(100)
)

Heatmap(
  mat,
  name = "Z-score",
  col = col_fun,

  # Clustering
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "ward.D2",

  # Labels
  show_row_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 9),

  # Cell appearance
  rect_gp = gpar(col = "black", lwd = 1),

  # Square cells (critical)
  width = ncol(mat) * unit(5, "mm"),
  height = nrow(mat) * unit(5, "mm"),

  use_raster = FALSE
)
```

---

## Common Anti-Patterns to Avoid

### 1. Bar Plots for Distributions
**Problem:** Hides distribution shape.
**Solution:** Strip plots, box plots, or violins.

### 2. Violin Plots with n < 50
**Problem:** Kernel density estimation unreliable.
**Solution:** Just plot the points.

### 3. Bidirectional Color Scales for Unidirectional Data
**Problem:** Using red-white-blue for expression falsely implies the midpoint is meaningful.
**Solution:** Sequential palette with one hue.

### 4. Pie Charts (Ever)
**Problem:** Humans read angles poorly.
**Solution:** Stacked or grouped bar charts.

### 5. Unordered Heatmaps
**Problem:** Random order obscures patterns.
**Solution:** Always cluster or sort meaningfully.

### 6. Rainbow Color Scales
**Problem:** Not perceptually uniform; colorblind-unfriendly.
**Solution:** Viridis, magma, or scico palettes.

### 7. Histograms with n < 100
**Problem:** Bin choice dramatically affects appearance.
**Solution:** Strip plots or dot plots.

---

## Domain-Specific Visualizations

### UMAP Styling

```r
library(ggrastr)

df <- data.frame(
  UMAP1 = Embeddings(obj, "umap")[, 1],
  UMAP2 = Embeddings(obj, "umap")[, 2],
  cell_type = obj$cell_type
)

ggplot(df, aes(x = UMAP1, y = UMAP2)) +
  geom_point_rast(aes(color = cell_type),
                  shape = 16, size = 0.1, alpha = 0.8,
                  raster.dpi = 300) +
  scale_color_manual(values = celltype_palette) +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 3)) +
  theme_void()
```

### UMAP with Direct Labels (shadowtext)

```r
library(shadowtext)
library(dplyr)

label_loc <- df %>%
  group_by(cell_type) %>%
  summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2), .groups = "drop")

ggplot(df, aes(x = UMAP1, y = UMAP2)) +
  geom_point_rast(aes(color = cell_type),
                  shape = 16, size = 0.1, alpha = 0.8, raster.dpi = 300) +
  geom_shadowtext(data = label_loc, aes(label = cell_type),
                  size = 3, color = "black", bg.color = "white") +
  scale_color_manual(values = celltype_palette) +
  theme_void() +
  theme(legend.position = "none")
```

### Feature Plots (Gene Expression on UMAP)

Use **light grey for zero** transitioning to **inverted magma** for non-zero values:

```r
library(viridis)

plotGene <- function(gene, obj) {
  df <- obj@meta.data
  df$UMAP1 <- Embeddings(obj, "umap")[, 1]
  df$UMAP2 <- Embeddings(obj, "umap")[, 2]
  df$Gene <- obj[["RNA"]]$data[gene, ]

  ggplot(df, aes(x = UMAP1, y = UMAP2)) +
    geom_point_rast(aes(color = Gene),
                    shape = 16, size = 0.1, alpha = 0.75, raster.dpi = 300) +
    scale_color_gradientn(
      colours = c("lightgrey", rev(magma(100))),
      name = "Expression",
      guide = guide_colorbar(ticks.colour = "black", frame.colour = "black", barwidth = 0.8)
    ) +
    ggtitle(gene) +
    theme_void() +
    theme(plot.title = element_text(size = 12, face = "italic"))
}
```

### Spatial Tissue Visualization

```r
library(ggrastr)

df <- as.data.frame(spatialCoords(sfe))
df$cell_type <- sfe$cell_type

ggplot(df, aes(x = x_centroid, y = y_centroid)) +
  geom_point_rast(aes(color = cell_type),
                  shape = 16, size = 0.1, alpha = 0.8, raster.dpi = 300) +
  scale_color_manual(values = celltype_palette) +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 3)) +
  coord_fixed() +
  theme_void()
```

### Segmented Cell Expression (ROI)

For spatial expression on segmented cells, use **scico `lapaz` inverted**:

```r
library(Voyager)
library(scico)

plotSpatialFeature(
  sfe,
  features = "GENE",
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
```

---

## References

- [Nature Research Figure Guide](https://research-figure-guide.nature.com/)
- [Martin Krzywinski's Points of View](https://mk.bcgsc.ca/pointsofview/)
- [Friends Don't Let Friends Make Bad Graphs](https://github.com/cxli233/FriendsDontLetFriends)
- [Fabio Crameri's Scientific Colour Maps](https://www.fabiocrameri.ch/colourmaps/)
- [Okabe & Ito Color Universal Design](https://jfly.uni-koeln.de/color/)
- [Data Visualization by Claus Wilke](https://clauswilke.com/dataviz/)
