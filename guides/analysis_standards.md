# Cook Lab Analysis Standards

Guidelines for reproducible, well-documented computational analysis—optimized for both human collaboration and AI-assisted workflows.

---

## Core Principles

1. **Plans are anchors.** Document the analysis plan before writing code. This prevents goal drift and gives both humans and AI a reference point.

2. **Logs are memory.** Maintain a living log of what's been done. This counters volatile context and makes it easy to resume work.

3. **Code is sequential.** Number scripts to make execution order explicit. Anyone should be able to reproduce the analysis by running scripts in order.

4. **Reports are evidence.** Every claim needs embedded statistics and figures. Visual interpretations require quantitative support.

5. **Failures are tracked.** Log what didn't work, not just what did. Silent failures compound into hidden errors.

---

## Project Structure

```
project/
├── README.md                 # Short what/why/how-to-reproduce
├── PROJECT_SPEC.md           # Research questions, goals, approach
├── ANALYSIS_PLAN.md          # Detailed analysis plan (the anchor)
├── ANALYSIS_LOG.md           # Living log of completed work
├── data/
│   ├── imaging/              # Organized by assay type
│   ├── spatial/
│   ├── single_cell/
│   └── external/             # Downloaded/public datasets (versioned, logged)
├── metadata/
│   └── samples.csv           # Single source of truth for sample info
├── scripts/
│   ├── 00_setup.R            # Environment, packages, paths
│   ├── 01_load_data.R        # Data import and initial QC
│   ├── 02_preprocessing.R    # Filtering, normalization, etc.
│   └── sandbox/              # Exploratory work (not part of main pipeline)
├── shellscripts/             # HPC job scripts (Slurm/PBS)
├── output/                   # Reproducible intermediates (gitignored)
├── figs/                     # Exploratory figure exports
├── reports/                  # Phase reports (markdown + PDF)
└── docs/
    └── manuscript/
        └── figures/          # Final paper figures + generation scripts
```

### Directory Descriptions

| Directory | Purpose |
|-----------|---------|
| `data/` | All project data, organized by assay type (imaging, spatial, single_cell). Subdirectories are flexible per project. `external/` holds downloaded/public datasets separate from lab-generated data. |
| `metadata/` | Single source of truth for sample information. `samples.csv` maps sample IDs to subjects, conditions, batches, tissues, library prep, file paths, and covariates. |
| `scripts/` | Numbered analysis scripts (00_, 01_, etc.) for reproducible pipeline. `sandbox/` nested inside for exploratory work. |
| `shellscripts/` | HPC job scripts (Slurm/PBS submission files). Keep parameterized so they regenerate `output/` from `data/`. |
| `output/` | Fully reproducible outputs—tables, matrices, processed objects, QC summaries. Nothing hand-edited; safe to delete and regenerate from code. |
| `figs/` | Working area for exploratory figures—quick plots, diagnostics, drafts. Intentionally fluid with many iterations. |
| `reports/` | Phase reports documenting analysis results with embedded statistics and figures. |
| `docs/manuscript/figures/` | Final paper figures with generation scripts. Versionable and reproducible, separate from exploratory graphics. |

**Principle:** Sub-directories within each top-level folder can be flexible and project-specific, but these top-level directories keep things in the right area.

---

## The Analysis Plan (ANALYSIS_PLAN.md)

The analysis plan is your anchor document. Write it before coding and update it as the plan evolves.

### Template

```markdown
# Analysis Plan: [Project Name]

**Last Updated:** [YYYY-MM-DD]
**Status:** [Planning | In Progress | Complete]

## Objective
[What question are we answering with this analysis?]

## Input Data
- [Dataset 1]: [description, location, n samples/cells]
- [Dataset 2]: [description, location]

## Analysis Phases

### Phase 1: [Name, e.g., "Quality Control"]
**Goal:** [What this phase accomplishes]
**Approach:**
1. [Step 1]
2. [Step 2]
**Success criteria:** [How we know this phase is complete]
**Output:** [Expected deliverables]

### Phase 2: [Name]
...

## Key Parameters & Decisions
| Parameter | Value | Rationale |
|-----------|-------|-----------|
| [e.g., min_genes] | [500] | [Standard QC threshold for 10x data] |

## Known Limitations
- [Limitation 1]
- [Limitation 2]
```

---

## The Analysis Log (ANALYSIS_LOG.md)

A living document that tracks everything done during the analysis. Update after every discrete step (script completion, key decision, failed attempt).

### Template

```markdown
# Analysis Log: [Project Name]

## Status
**Current phase:** [e.g., Phase 2 - Clustering]
**Current task:** [What's in progress or next]
**Last updated:** [YYYY-MM-DD HH:MM]

## Quick Summary
[2-3 sentences on current state. What has been done? What's the immediate next step?]

---

## Log

### [YYYY-MM-DD HH:MM] - [Script or task name]
**Type:** [Script | Decision | Failed attempt | Exploration]
**Phase:** [Phase name]
**Script:** [e.g., `02_preprocessing.R` or `sandbox/explore_clustering.R`]
**Status:** [Complete | Staged | Abandoned]

**What was done:**
- [Specific action with numbers]

**Key outputs:**
- [File created: path/to/file]
- [Result: key numbers or findings]

**Decisions:**
- [Decision]: [Rationale]

**Issues:**
- [Issue]: [Resolution or "UNRESOLVED"]

---

### [YYYY-MM-DD HH:MM] - [Previous entry]
...
```

### Log Guidelines

- **Update after every discrete step**: script completion, key decision, or failed attempt
- **Be specific with numbers**: "Filtered to 8,234 cells (removed 1,766 low-quality)" not "Filtered cells"
- **Log failures explicitly**: What was tried, why it failed, what was done instead
- **Track file outputs**: Note paths to generated data, figures, tables
- **Distinguish script types**: Note if work is in sandbox vs staged in scripts/
- **Timestamp precisely**: Include time when multiple updates per day

---

## Script Conventions

### Numbering
- Use two-digit prefixes: `00_`, `01_`, `02_`, ...
- `00_setup.R` is always environment setup
- Scripts should be runnable in order to reproduce the analysis

### Structure
Every script should include:

```r
# =============================================================================
# Script: 01_load_data.R
# Project: [Project name]
# Purpose: [What this script does]
# Author: [Name]
# Date: [YYYY-MM-DD]
# =============================================================================

# Dependencies ----------------------------------------------------------------
source("scripts/00_setup.R")

# [Main code organized in sections] -------------------------------------------

# Save outputs ----------------------------------------------------------------
saveRDS(object, "output/object_name.rds")

# Session info ----------------------------------------------------------------
sessionInfo()
```

### Naming
- Use snake_case for scripts and variables
- Be descriptive: `02_filter_and_normalize.R` not `02_process.R`
- Sandbox scripts use descriptive names without numbers (e.g., `explore_clustering.R`)

---

## Sandbox Workflow

The sandbox (`scripts/sandbox/`) is for dynamic exploration—testing parameters, trying approaches, iterating on visualizations. Once an approach is finalized, stage it as a numbered script in the parent `scripts/` directory.

### Sandbox vs Staged Scripts

| Aspect | Sandbox (`scripts/sandbox/`) | Staged (`scripts/`) |
|--------|------------------------------|---------------------|
| Purpose | Exploration, iteration | Reproducible pipeline |
| Naming | Descriptive (e.g., `explore_clustering.R`) | Numbered (e.g., `03_clustering.R`) |
| State | Work in progress | Finalized |
| Dependencies | May have loose dependencies | Strict sequential order |
| Outputs | Temporary, may overwrite | Versioned, documented |

### Working in the Sandbox

1. **Create exploration scripts** with descriptive names:
   ```
   scripts/sandbox/
   ├── explore_clustering_resolution.R
   ├── compare_normalization_methods.R
   └── test_marker_genes.R
   ```

2. **Iterate freely**: Try different parameters, visualizations, approaches

3. **Log exploration** in ANALYSIS_LOG.md with `Type: Exploration`:
   ```markdown
   ### 2025-02-05 14:30 - Exploring clustering resolution
   **Type:** Exploration
   **Script:** `sandbox/explore_clustering_resolution.R`
   **Status:** In progress

   **What was done:**
   - Tested resolutions 0.3, 0.5, 0.8, 1.0
   - Resolution 0.5 gives 12 clusters with good separation

   **Decisions:**
   - Proceeding with resolution 0.5: balances granularity vs interpretability
   ```

4. **Save intermediate outputs** to the exploratory figures directory:
   ```r
   # In sandbox scripts, save to figs/ (exploratory area)
   ggsave("figs/clustering_comparison.png", p)
   ```

### Staging a Script

When exploration is complete and you've settled on an approach:

1. **Create the staged script** in `scripts/` with the next number:
   ```r
   # scripts/03_clustering.R - staged from sandbox/explore_clustering_resolution.R
   ```

2. **Clean up the code**:
   - Remove dead ends and commented alternatives
   - Add proper header with documentation
   - Ensure it sources `00_setup.R`
   - Save outputs to `output/` and final figures to `docs/manuscript/figures/`

3. **Update ANALYSIS_LOG.md** to mark the staging:
   ```markdown
   ### 2025-02-05 16:00 - Staged clustering script
   **Type:** Script
   **Script:** `scripts/03_clustering.R`
   **Status:** Staged
   **Staged from:** `scripts/sandbox/explore_clustering_resolution.R`

   **What was done:**
   - Finalized clustering at resolution 0.5
   - 12 clusters identified

   **Key outputs:**
   - `output/seurat_clustered.rds`
   - `docs/manuscript/figures/umap_clusters.png`
   ```

4. **Archive or delete sandbox script** (optional):
   - Keep if it contains useful reference code
   - Delete if fully captured in staged script

### Sandbox Guidelines

- **Don't skip the sandbox** for non-trivial decisions—exploration is valuable
- **Log sandbox work** even if it doesn't lead anywhere (failed attempts matter)
- **Don't stage prematurely**—make sure the approach is settled before cleaning up
- **Keep sandbox outputs separate** from main results to avoid confusion

---

## Reports

Reports document the results of each analysis phase with embedded evidence.

### When to Write Reports
- After completing each major phase (QC, clustering, differential expression, etc.)
- When making key decisions that affect downstream analysis
- At project milestones

### Report Structure

```markdown
# [Phase] Report: [Project Name]

**Date:** [YYYY-MM-DD]
**Analyst:** [Name]
**Scripts:** [List of scripts covered]

## Summary
[2-3 sentence overview of this phase]

## Methods
[Brief description of what was done]

## Results

### [Section 1]
[Text describing findings with embedded statistics]

![Figure description](../results/figures/figure_name.png)
*Figure 1: [Caption with key statistics]*

| Metric | Value |
|--------|-------|
| [e.g., Cells passing QC] | [10,234] |

### [Section 2]
...

## Key Decisions
- [Decision]: [Rationale supported by data]

## Issues & Limitations
- [Any problems encountered or caveats]

## Next Steps
- [What follows from this phase]
```

### Report Guidelines

1. **Every claim needs evidence**
   - Good: "Cluster 5 shows elevated MKI67 expression (log2FC = 2.3, adj.p < 0.001), suggesting proliferative cells"
   - Bad: "Cluster 5 appears to be proliferating based on the UMAP"

2. **Quantify visual observations**
   - If referencing a UMAP or heatmap, include supporting statistics
   - Report cluster sizes, expression values, statistical tests

3. **Embed figures inline**
   - Use relative paths: `![](../figs/fig.png)` for exploratory, `![](../docs/manuscript/figures/fig.png)` for final
   - Include informative captions with key numbers

4. **Be conservative with interpretation**
   - Distinguish between "shows" (data) and "suggests" (interpretation)
   - Flag when conclusions depend on visual assessment

### Rendering to PDF

Reports are written in markdown and rendered to PDF using Quarto (preferred) or Pandoc.

#### Setup (one-time per project)

Copy the rendering config to your project's `reports/` directory:

```bash
# From lab_guide templates
cp templates/reports/_quarto.yml your_project/reports/
cp templates/reports/Makefile your_project/reports/
```

#### Rendering

```bash
# Render a single report
cd reports/
make 01_qc_report.pdf

# Render all reports
make

# Or use the standalone script from anywhere
./render_report.sh reports/01_qc_report.md --open
```

#### Fonts

Reports use system fonts (Helvetica Neue on macOS, no installation required). If on Windows or Linux, edit `_quarto.yml`:
- Windows: Arial + Consolas
- Linux: DejaVu Sans + DejaVu Sans Mono

#### Troubleshooting

If rendering fails:
1. Check Quarto is installed: `quarto --version`
2. Check XeLaTeX is installed: `xelatex --version` (comes with MacTeX or TinyTeX)
3. Install TinyTeX if needed: `quarto install tinytex`

---

## Working with AI Assistants

Guidelines for effective AI-assisted analysis:

### Before Starting
- Ensure PROJECT_SPEC.md and ANALYSIS_PLAN.md exist
- Point AI to these documents at the start of each session
- Share the current ANALYSIS_LOG.md for context

### During Analysis
- Ask AI to update ANALYSIS_LOG.md after completing work
- Request explicit logging of any errors or failed approaches
- For visual interpretation tasks, ask for quantitative support

### Handoff Protocol
When ending a session or switching contexts:
1. Update ANALYSIS_LOG.md with current state
2. Note any in-progress work
3. List immediate next steps
4. Save all modified scripts

### Red Flags
Watch for these AI analysis pitfalls:
- **Over-interpretation of visual patterns**: Always request statistical backing
- **Silent failures**: Ask explicitly "Did anything fail or produce warnings?"
- **Goal drift**: Check claims against ANALYSIS_PLAN.md
- **Hallucinated results**: Verify key numbers in the actual outputs

---

## Checklist

### Starting a new analysis
- [ ] PROJECT_SPEC.md exists with clear research questions
- [ ] ANALYSIS_PLAN.md drafted with phases and success criteria
- [ ] ANALYSIS_LOG.md initialized
- [ ] Directory structure created
- [ ] `00_setup.R` created with environment details

### Completing a phase
- [ ] All scripts numbered and documented
- [ ] ANALYSIS_LOG.md updated with completed work
- [ ] Report written with embedded statistics and figures
- [ ] Report rendered to PDF
- [ ] Key decisions documented with rationale

### Ending a work session
- [ ] ANALYSIS_LOG.md updated
- [ ] Scripts saved and documented
- [ ] Next steps noted
