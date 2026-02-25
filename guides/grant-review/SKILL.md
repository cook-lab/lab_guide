---
name: grant-review
description: >
  Multi-agent evaluation of research grant proposals and resubmissions. Use when the user wants
  to review, evaluate, critique, or assess a grant proposal, grant resubmission, or funding
  application. Supports parallel agent review from multiple independent perspectives, synthesis
  of findings into a thematic report, and generation of a comprehensive Word document. Triggers
  include: "review this grant", "evaluate this proposal", "assess this resubmission", "critique
  this application", or any request involving grant/proposal feedback across agencies (CIHR, NIH,
  ERC, NSF, NSERC, etc.). Also triggers when the user provides grant-related documents (proposals,
  reviewer comments, response to reviewers) and asks for feedback or evaluation.
---

# Grant Review

Multi-perspective agent-based evaluation of research grant proposals with synthesized reporting.

## Workflow

1. **Gather documents** from the grant directory
2. **Extract text** from all documents (.docx, .pdf)
3. **Launch 2-4 review agents in parallel** with distinct perspectives
4. **Synthesize** all agent reports into a thematic analysis
5. **Generate Word document** using the docx skill

## Step 1: Gather and Extract Documents

Scan the grant directory for:
- **Research proposal** (.docx or .pdf) -- the core document all agents evaluate
- **Previous reviewer comments** (.pdf) -- needed for resubmission evaluations
- **Response to reviewers** (.docx) -- needed for resubmission evaluations

Extract text from .docx files using `textutil -convert txt -stdout <file>` (macOS) or equivalent.
Extract text from PDFs using the Read tool. Save extracted proposal text to a temp file so agents
can read it by path.

## Step 2: Select and Launch Review Agents

Launch agents in parallel using the Task tool (`subagent_type: "general-purpose"`). Select the
agent panel based on what documents are available.

### For new proposals (no prior reviews):
- **Agent A**: Fresh Science & Narrative reviewer
- **Agent B**: Methods & Feasibility specialist
- **Agent C** (optional): Grant Strategy reviewer

### For resubmissions (prior reviews + response document available):
- **Agent A**: Fresh Science & Narrative reviewer (has NOT seen previous reviews)
- **Agent B**: Previous Reviewer Perspective (has seen reviews + response document)
- **Agent C**: Methods & Feasibility specialist (has NOT seen previous reviews)
- **Agent D**: Grant Strategy & Response Document reviewer (has seen reviews + response)

Agent A and C provide unbiased fresh perspectives. Agent B and D evaluate whether prior concerns
were addressed. This separation is important -- do not give previous reviews to fresh reviewers.

See [references/agent-prompts.md](references/agent-prompts.md) for detailed prompt templates.

### Prompt construction principles

- Give each agent the full proposal text (by file path)
- For resubmission agents (B, D): include full text of previous reviews AND response to reviewers
  directly in the prompt (not just file paths, since agents work more reliably with inline text)
- Define a structured output format (numbered sections with headers)
- Ask for specific, actionable items -- not vague praise or criticism
- Instruct agents to list strengths AND weaknesses with concrete evidence from the proposal

## Step 3: Synthesize Agent Reports

This is the most important step. Do NOT simply concatenate agent reports.

### Organize by theme, not by agent

Catalog every weakness and suggestion from all agents, then group thematically:
- Causality / interpretive challenges
- Model systems and their limitations
- Human cohort / clinical design
- Statistical rigor and reproducibility
- Aim-specific concerns
- New vulnerabilities introduced by revisions (for resubmissions)

### For each item, note:
- **Agent attribution**: Which agents flagged it (e.g., "Agents A, B, C" or "Agent C only")
- **Consensus level**: All agents, majority, or single-agent finding
- **Perspective annotation**: Brief editorial assessment of validity, severity, and ease of fix

### Critical rule: more info is better for weaknesses

Do NOT filter to consensus-only items. A concern flagged by a single agent may be the most
critical finding in the report. Consensus is a metric for *likely reviewer impact*, but
single-agent findings flagged as critical should be prominently featured with a note like
"[Agent C -- critical point not raised by other agents]".

### Strengths: consensus first, then additional

List consensus strengths (flagged by all agents) first, then additional strengths noted by 1-2
agents with attribution.

## Step 4: Generate the Word Document

Use the docx skill to create a professionally formatted report. See
[references/report-structure.md](references/report-structure.md) for the full document structure,
section templates, and formatting guidance.

### Key report sections:
1. **Title page** with grant metadata and scoring context
2. **Part 1: Strengths** (consensus + additional)
3. **Part 2: All Weaknesses** (thematic, with agent attribution + perspective annotations)
4. **Part 3: Response to Reviewers Critique** (for resubmissions)
5. **Part 4: Scorecard** (concern-by-concern tracking table, for resubmissions)
6. **Part 5: All Suggestions** (tiered by impact)
7. **Part 6: Bottom Line** (overall assessment)

### Perspective annotations

For each weakness, include a purple italic "Perspective:" paragraph that provides editorial
judgment: Is this concern valid? How severe? How easy to fix? Will reviewers care? This is the
most valuable part of the report -- it turns a list of concerns into actionable guidance.

Example labels: "CRITICAL FIX", "EASY FIX with meaningful impact", "Valid but standard limitation
in the field", "Low priority", "Worth a brief mention".

### Suggestion tiers

Organize all suggestions into tiers:
- **Tier 1**: Response to reviewers fixes (highest impact, easiest changes)
- **Tier 2**: Proposal text revisions (high priority)
- **Tier 3**: Experimental design improvements
- **Tier 4**: Analytical/computational enhancements

## Adapting to Different Agencies

The agent prompts in the reference file use CIHR as the default. Adapt terminology:
- **CIHR**: Scientific Officer Notes, committee discussion, 0-5 scoring
- **NIH**: Summary Statement, study section critique, 1-9 scoring (1=best)
- **NSF**: Panel Summary, individual reviews, E/V/G/F/P ratings
- **ERC**: Evaluation Summary Report, panel comments

Adjust agent prompts to reference the appropriate review structure and scoring conventions.
