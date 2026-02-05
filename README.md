# Cook Lab Guide

Standards, practices, and Claude Code skills for the Cook Lab at OHRI/University of Ottawa.

## What's Here

```
lab_guide/
├── guides/
│   ├── analysis_standards.md      # Project structure, logging, reports
│   └── visualization_style_guide.md  # Palettes, themes, domain-specific plots
├── templates/
│   ├── project_spec.md            # Research project template
│   └── reports/                   # Quarto config, Makefile for PDF reports
└── .claude/skills/                # Claude Code skills
    ├── analysis-assistant.md      # Analysis workflow guidance
    ├── visualization.md           # Plotting standards
    ├── scrna-spatial.md           # scRNA-seq and spatial processing
    └── project-spec.md            # Project specification creation
```

## Quick Setup

```bash
# Clone the lab guide
git clone https://github.com/cook-lab/lab_guide.git ~/.cook-lab-guide

# Install Claude Code skills
mkdir -p ~/.claude/skills
cp ~/.cook-lab-guide/.claude/skills/* ~/.claude/skills/

# Verify
ls ~/.claude/skills/
```

## Skills Overview

| Skill | Purpose |
|-------|---------|
| **analysis-assistant** | Project structure, logging, sandbox workflow, reports |
| **visualization** | Lab theme, palettes, UMAP/spatial/heatmap templates |
| **scrna-spatial** | scRNA-seq and spatial transcriptomics best practices |
| **project-spec** | Creating research project specifications |

## Using in Your Projects

### Reference standards in your CLAUDE.md

```markdown
## Lab Standards

Follow Cook Lab standards from ~/.cook-lab-guide/:
- Analysis: guides/analysis_standards.md
- Visualization: guides/visualization_style_guide.md
```

### Copy templates

```bash
# Report setup
cp -r ~/.cook-lab-guide/templates/reports/ your_project/reports/

# Project spec
cp ~/.cook-lab-guide/templates/project_spec.md your_project/PROJECT_SPEC.md
```

## Updating

```bash
cd ~/.cook-lab-guide && git pull
cp ~/.cook-lab-guide/.claude/skills/* ~/.claude/skills/
```

## Related Repositories

- [workflows](https://github.com/cook-lab/workflows) — scRNA-seq and spatial tutorials with runnable scripts
- [scrna-integration-pipeline](https://github.com/cook-lab/scrna-integration-pipeline) — Multi-sample integration for HPC
- [tumour_segmentation](https://github.com/cook-lab/tumour_segmentation) — Tumour region segmentation for spatial data

## About the Lab

The Cook Lab combines single-cell and spatial genomics with patient-derived models to decode gynecologic disease. Based at the Ottawa Hospital Research Institute and University of Ottawa.

- **Website:** [cooklab.ca](https://www.cooklab.ca)
- **PI:** Dr. David Cook
