---
title: "[Phase] Report: [Project Name]"
author: "[Your Name]"
date: today
---

# Summary

[2-3 sentence overview of this analysis phase. What was done and what are the key findings?]

**Scripts covered:** `01_xxx.R` through `03_xxx.R`

---

# Methods

[Brief description of the analytical approach. Include key parameters and software versions.]

---

# Results

## [Section 1: e.g., Quality Control]

[Describe findings with embedded statistics. Be specific with numbers.]

![Description of figure](../results/figures/figure_name.png){width=80%}

*Figure 1: [Informative caption including key statistics, e.g., "UMAP visualization of 10,234 cells colored by cluster identity (n=12 clusters, resolution=0.5)"]*

| Metric | Value |
|--------|-------|
| Total cells | 12,000 |
| Cells passing QC | 10,234 (85.3%) |
| Median genes/cell | 2,456 |
| Median UMIs/cell | 8,234 |

## [Section 2]

[Continue with additional results sections as needed.]

---

# Key Decisions

- **[Decision 1]:** [Rationale supported by data, e.g., "Used 500 gene minimum threshold based on distribution showing clear bimodal separation at this cutoff"]

- **[Decision 2]:** [Rationale]

---

# Issues & Limitations

- [Any problems encountered, caveats, or limitations to note]

---

# Next Steps

- [ ] [Immediate next task]
- [ ] [Following task]

---

# Session Info

```
R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)

Packages:
- Seurat 5.0.1
- tidyverse 2.0.0
- [etc.]
```
