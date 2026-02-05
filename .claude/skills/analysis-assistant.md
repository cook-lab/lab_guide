# Analysis Assistant

Guide reproducible computational analysis for Cook Lab research projects following lab standards.

## When to Use

Use this skill when:
- Starting a new analysis project
- Resuming work on an existing analysis
- Writing or updating analysis plans, logs, or reports
- The user asks for help with bioinformatics/computational analysis
- Exploring approaches in the sandbox before staging final scripts

## Critical: Anchoring & Memory

**AI analysis work is prone to goal drift, volatile memory, and hidden failures.** Counter this by:

1. **Always check for anchor documents first:**
   - `PROJECT_SPEC.md` - What are we trying to answer?
   - `ANALYSIS_PLAN.md` - What's the plan?
   - `ANALYSIS_LOG.md` - What's been done?

2. **Update the log after every discrete step** - Not just sessions, but each script/task/decision

3. **Explicitly track failures** - Log what didn't work, not just successes

4. **Reference the plan** - Before making decisions, check if they align with stated goals

## Reference Document

Full standards: `guides/analysis_standards.md`

---

## Log Format (Updated)

Log at the granularity of each discrete step. Use this format:

```markdown
### [YYYY-MM-DD HH:MM] - [Script or task name]
**Type:** [Script | Decision | Failed attempt | Exploration]
**Phase:** [Phase name]
**Script:** [e.g., `02_preprocessing.R` or `sandbox/explore_clustering.R`]
**Status:** [Complete | Staged | Abandoned | In progress]

**What was done:**
- [Specific action with numbers]

**Key outputs:**
- [File created: path/to/file]
- [Result: key numbers or findings]

**Decisions:**
- [Decision]: [Rationale]

**Issues:**
- [Issue]: [Resolution or "UNRESOLVED"]
```

### Log Entry Types

- **Script**: Completed work on a numbered script in `scripts/`
- **Exploration**: Work in `sandbox/` testing approaches
- **Decision**: Key parameter or approach decision (even without code changes)
- **Failed attempt**: Something that didn't work (important to track!)

---

## Sandbox Workflow

Use the sandbox for dynamic exploration. Stage to `scripts/` when finalized.

### When to Use Sandbox

- Testing different parameters (clustering resolution, thresholds)
- Comparing approaches (normalization methods, integration strategies)
- Iterating on visualizations
- Any non-trivial decision that benefits from exploration

### Sandbox → Staged Workflow

1. **Explore in sandbox:**
   ```
   sandbox/explore_clustering_resolution.R
   ```

2. **Log exploration** (Type: Exploration)

3. **When approach is settled, stage:**
   - Create `scripts/03_clustering.R`
   - Clean up code, add proper header
   - Save outputs to `data/processed/` and `results/`

4. **Log staging** (Type: Script, note "Staged from: sandbox/xxx.R")

### Key Rule

**Don't stage prematurely.** Make sure the approach is settled before moving out of sandbox. It's okay to have multiple exploration sessions before staging.

---

## Report Standards

### Evidence Requirements

**Every claim must have quantitative support:**

| Type | Good | Bad |
|------|------|-----|
| Cell identity | "Cluster 5 is enriched for MKI67 (log2FC=2.3, p<1e-10)" | "Cluster 5 looks proliferative" |
| Visual pattern | "UMAP shows separation (silhouette=0.72)" | "Clear separation on UMAP" |
| Comparison | "Method A outperforms B (AUC 0.89 vs 0.76)" | "Method A works better" |

**For visual observations:**
- Compute supporting statistics
- Report effect sizes and p-values
- If purely visual, flag explicitly: *"Visual inspection suggests X (requires quantitative validation)"*

### Rendering Reports

Reports render to PDF with clean sans-serif styling (system fonts, no installation needed). Setup:

```bash
# Copy config to project
cp templates/reports/_quarto.yml project/reports/
cp templates/reports/Makefile project/reports/

# Render
cd project/reports/
make 01_qc_report.pdf

# Or render all
make
```

---

## Workflow by Stage

### Starting a New Analysis

1. **Check for PROJECT_SPEC.md** (create if missing using project-spec skill)

2. **Create ANALYSIS_PLAN.md:**
   ```markdown
   # Analysis Plan: [Project Name]

   **Last Updated:** [YYYY-MM-DD]
   **Status:** Planning

   ## Objective
   [What question are we answering?]

   ## Input Data
   - [Dataset]: [description, location]

   ## Analysis Phases

   ### Phase 1: [Name]
   **Goal:** [What this accomplishes]
   **Approach:**
   1. [Step]
   **Success criteria:** [How we know it's complete]
   **Output:** [Deliverables]

   ## Key Parameters & Decisions
   | Parameter | Value | Rationale |
   |-----------|-------|-----------|

   ## Known Limitations
   ```

3. **Initialize ANALYSIS_LOG.md:**
   ```markdown
   # Analysis Log: [Project Name]

   ## Status
   **Current phase:** Phase 1
   **Current task:** [First task]
   **Last updated:** [YYYY-MM-DD HH:MM]

   ## Quick Summary
   Analysis initialized. Starting with [first task].

   ---

   ## Log

   ### [YYYY-MM-DD HH:MM] - Project initialized
   **Type:** Script
   **Phase:** Setup
   **Script:** `00_setup.R`
   **Status:** Complete

   **What was done:**
   - Created analysis plan
   - Set up project structure
   - Created setup script

   **Next steps:**
   - [First actual task]
   ```

4. **Create directory structure:**
   ```
   scripts/
   data/raw/
   data/processed/
   results/figures/
   results/tables/
   reports/
   sandbox/
   sandbox/outputs/
   ```

5. **Create `scripts/00_setup.R`**

6. **Copy report config:**
   ```bash
   cp templates/reports/_quarto.yml reports/
   cp templates/reports/Makefile reports/
   ```

### Resuming an Analysis

1. **Read anchor documents in order:**
   - PROJECT_SPEC.md
   - ANALYSIS_PLAN.md
   - ANALYSIS_LOG.md (focus on Status section and recent entries)

2. **Confirm state with user:**
   > "Based on the log, you're in [Phase] working on [task]. Last completed: [X]. Should I continue with [next step]?"

### During Analysis

#### After Every Discrete Step

Update ANALYSIS_LOG.md with:
- Timestamp and task name
- Type (Script/Exploration/Decision/Failed attempt)
- What was done (specific, with numbers)
- Key outputs (file paths)
- Any decisions or issues

#### When Something Fails

**Always log it:**
```markdown
### 2025-02-05 15:30 - Failed: Integration with Harmony
**Type:** Failed attempt
**Phase:** Integration
**Script:** `sandbox/test_harmony.R`
**Status:** Abandoned

**What was done:**
- Attempted Harmony integration with default parameters
- Also tried theta=2, theta=4

**Issues:**
- Over-correction: cell types no longer separated after integration
- Decided to try RPCA instead

**Next steps:**
- Test RPCA integration in sandbox
```

#### Before Making Decisions

Check against ANALYSIS_PLAN.md:
- Does this align with the stated objective?
- Is this within the planned approach?
- If deviating, log the decision with rationale

### Writing Reports

After completing each major phase:

1. **Create report** in `reports/` using template
2. **Embed all figures** with informative captions
3. **Include statistics** for every claim
4. **Render to PDF:** `make [report].pdf`
5. **Log the report** in ANALYSIS_LOG.md

### Ending a Session

Before stopping:

1. **Update ANALYSIS_LOG.md Status section:**
   ```markdown
   ## Status
   **Current phase:** Phase 2 - Clustering
   **Current task:** Testing clustering parameters in sandbox
   **Last updated:** 2025-02-05 17:30

   ## Quick Summary
   Completed preprocessing (10,234 cells). Currently exploring clustering
   resolutions in sandbox. Resolution 0.5 looks promising but need to
   verify marker expression before staging.
   ```

2. **Summarize for user:**
   > "Updated the analysis log. Current state: [X]. Next steps: [Y]. The log has full details for resuming."

---

## Common Pitfalls to Avoid

| Pitfall | How to Avoid |
|---------|--------------|
| Over-interpreting visuals | Always compute supporting statistics |
| Silent failures | Log every failed attempt with "Type: Failed attempt" |
| Goal drift | Check decisions against ANALYSIS_PLAN.md |
| Skipping logging | Update after every discrete step, not just sessions |
| Premature staging | Keep exploring in sandbox until approach is settled |
| Vague claims | Be specific: numbers, file paths, parameters |
| Lost context | Keep Status section current for easy resumption |

---

## Quick Reference

### Directory Structure
```
project/
├── PROJECT_SPEC.md
├── ANALYSIS_PLAN.md
├── ANALYSIS_LOG.md
├── scripts/          # Numbered, finalized scripts
├── sandbox/          # Exploration scripts
│   └── outputs/      # Temporary sandbox outputs
├── data/raw/         # Immutable input data
├── data/processed/   # Generated data objects
├── results/figures/  # Final figures
├── results/tables/   # Final tables
└── reports/          # Markdown + PDF reports
```

### Script Naming
- `scripts/00_setup.R` - Always first
- `scripts/01_load_data.R`, `02_preprocessing.R`, etc.
- `sandbox/explore_[topic].R` - Exploration scripts

### Log Entry Types
- **Script**: Completed numbered script
- **Exploration**: Sandbox work
- **Decision**: Key choice made
- **Failed attempt**: Something that didn't work
