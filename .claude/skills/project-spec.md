# Project Spec Creator

Guide users through creating or updating project specification documents for Cook Lab research projects.

## When to Use

Use this skill when the user wants to:
- Create a new project spec for a research project
- Update or refine an existing project spec
- Review a project spec for completeness
- Convert informal project notes into a structured spec

## Template Location

The project spec template is located at: `templates/project_spec.md`

## Workflow

### Creating a New Project Spec

1. **Gather core information** through conversation:
   - Project name and current status
   - Who is leading the project
   - 2-3 sentence summary of what the project is about

2. **Define the research questions:**
   - Ask: "What are the key questions this project aims to answer?"
   - Help refine vague questions into specific, testable hypotheses
   - Prioritize: identify the primary question vs. secondary/exploratory questions

3. **Establish goals and deliverables:**
   - Ask: "What does success look like for this project?"
   - Distinguish between goals (outcomes) and deliverables (tangible outputs)
   - Common deliverables in the lab: datasets, analysis pipelines, figures, manuscripts, presentations

4. **Document the approach:**
   - Experimental strategy: samples, models, assays
   - Computational strategy: analysis methods, pipelines, tools
   - For genomics projects, clarify: sequencing platform, expected cell/sample numbers, analysis frameworks

5. **Catalog resources:**
   - Samples & models (patient samples, cell lines, PDX, organoids)
   - Data types (scRNA-seq, spatial transcriptomics, bulk RNA-seq, etc.)
   - Key software and tools
   - Public datasets being integrated

6. **Identify team and collaborators:**
   - Who is working on this within the lab?
   - External collaborators and their roles

7. **Set timeline and milestones:**
   - Major milestones with target dates
   - Be realistic about timelines for genomics projects (data generation, analysis iterations, manuscript preparation)

### Updating an Existing Project Spec

1. Read the current spec
2. Ask what aspects need updating
3. Make targeted updates while preserving the overall structure
4. Add an entry to the "Notes & Updates" section with the date and summary of changes

## Writing Guidelines

- **Be specific:** "Characterize immune cell populations in HGSC" is better than "Analyze the data"
- **Use active voice:** "Identify markers of treatment resistance" not "Markers will be identified"
- **Keep it scannable:** Use bullet points and tables; avoid long paragraphs
- **Include context:** Brief background helps onboard new team members and AI agents
- **Update regularly:** A stale spec is worse than no spec

## Common Project Types in Cook Lab

### Single-cell RNA-seq Projects
- Typical deliverables: processed Seurat/AnnData object, cell type annotations, differential expression results, figures
- Key resources to document: sample metadata, sequencing depth, clustering parameters

### Spatial Transcriptomics Projects
- Typical deliverables: spatial object, cell segmentation results, niche analysis, spatial statistics
- Key resources to document: platform (Xenium, Visium, etc.), tissue sections, region annotations

### Biomarker Discovery Projects
- Typical deliverables: candidate biomarker list, validation cohort results, diagnostic performance metrics
- Key resources to document: discovery cohort, validation cohort, assay platforms

### Methods Development Projects
- Typical deliverables: software package, benchmarking results, documentation, publication
- Key resources to document: test datasets, comparison methods, compute resources

## Output

Save new project specs to: `projects/[project-name]/PROJECT_SPEC.md`

If a project directory doesn't exist, create it.
