# Report Structure and Formatting

Guide for generating the comprehensive Word document evaluation report.

## Document Generation

Use the `docx` skill (docx-js via Node.js) to generate a professional .docx file. Install with
`npm install -g docx` if needed. Run with `NODE_PATH=/opt/homebrew/lib/node_modules node script.js`
(macOS) or `NODE_PATH=$(npm root -g) node script.js` (cross-platform) if global modules are not
found.

## Color Palette

```javascript
const DARK_BLUE = "1B3A5C";   // H1 headings, table headers
const MED_BLUE = "2C5F8A";    // H2/H3 headings, callout box borders
const LIGHT_BLUE_BG = "E8F0F8"; // Callout box background
const GREEN = "2E7D32";       // "Fully Addressed" ratings, recommendations
const RED = "B71C1C";         // "Core Challenge" callouts
const AMBER = "E65100";       // "Partially Addressed" ratings, warnings
const GRAY = "666666";        // Metadata text, agent attribution
const PURPLE = "4A148C";      // Perspective annotations
```

## Document Sections

### Title Page
- Grant title, PI, institution, application number
- Scoring context callout box (original score, individual scores, target score, committee)
- Agent panel description (which perspectives were used)
- Date

### Part 1: Strengths
- **Consensus Strengths** (flagged by all agents): Numbered list with bold lead sentence +
  supporting detail. Each item should be 2-4 sentences explaining WHY this is a strength.
- **Additional Strengths** (1-2 agents): Bullet list with agent attribution in brackets.

### Part 2: All Identified Weaknesses (thematic)

Organize into themes. Each theme gets an H2 heading. Each item within a theme gets:

1. **H3 heading**: Numbered (e.g., "1.1 Fibrosis is not isolated from other changes")
2. **Body paragraph**: Full description of the concern (2-4 sentences)
3. **Agent tag line**: Gray italic text: `[Flagged by: Agents A, B, C]` or
   `[Agent C -- critical point not raised by other agents]`
4. **Perspective annotation**: Purple italic paragraph with "Perspective: " bold prefix.
   Editorial assessment of validity, severity, ease of fix, and strategic advice.

For the most critical items, use a colored callout box:
- Red box (`RED` border, `"FBE9E7"` background): Core challenges flagged by all agents
- Amber box (`AMBER` border, `"YELLOW_BG"` background): Important warnings

### Part 3: Response to Reviewers Critique (resubmissions only)

Sections:
- **A. Tone and Diplomacy**: Numbered list of problematic passages with specific quotes and
  suggested rewrites
- **B. Structural Gaps**: Missing sections (SO Notes, unaddressed concerns, no summary of
  revisions, no closing statement)
- **C. Quality of Individual Responses**: Which work well, which need improvement
- **D. Recommended Document Structure**: Suggested reordering

### Part 4: Scorecard (resubmissions only)

Table with two columns:
- Column 1: Reviewer concern (e.g., "R2-W1: Effect size / biological importance")
- Column 2: Rating with color-coded background
  - "FULLY ADDRESSED" (green background)
  - "PARTIALLY ADDRESSED" (amber background)
  - "NOT ADDRESSED" (red background)

Summary line below table: "Fully Addressed: X/Y | Partially Addressed: Z/Y"

### Part 5: All Suggestions for Improvement

Organized in 4 tiers, each with a numbered list:

- **Tier 1: Response to Reviewers (Highest Impact)**: Fixes to the response document --
  typically the easiest and highest-value changes
- **Tier 2: Proposal Revisions -- High Priority**: Text and figure changes addressing
  persistent reviewer concerns
- **Tier 3: Proposal Revisions -- Experimental Design**: Methodological improvements
- **Tier 4: Analytical and Computational Enhancements**: Additional approaches

Each item: bold action phrase + normal explanatory text.

### Part 6: Bottom Line

Single callout box (blue) containing:
- Opening bold statement (overall assessment in one sentence)
- 2-3 paragraphs covering: scientific strength, biggest risk to score, deepest scientific
  challenge, most concrete feasibility concern
- Closing bold statement on fundability prospects

## Helper Functions

Standard helper functions for the docx-js generation script:

```javascript
// Heading helpers
function h1(text) { /* HeadingLevel.HEADING_1, DARK_BLUE, size 32 */ }
function h2(text) { /* HeadingLevel.HEADING_2, MED_BLUE, size 26 */ }
function h3(text) { /* bold, MED_BLUE, size 22 */ }

// Text helpers
function p(text) { /* normal paragraph, size 21 */ }
function bold(text) { /* bold TextRun */ }
function normal(text) { /* normal TextRun */ }
function italic(text) { /* italic TextRun */ }

// Agent attribution
function tagLine(agents) {
  // Gray italic: [Flagged by: {agents}]
}

// Perspective annotation
function perspective(text) {
  // Purple italic with bold "Perspective: " prefix, indented left 360
}

// Callout box (colored border + background containing child paragraphs)
function calloutBox(children, borderColor, bgColor) { }

// Scorecard table row
function ratingRow(concern, rating, ratingColor) { }

// Bullets and numbered lists
function bullet(runs) { }
function numItem(runs, numberingRef) { }
```

## Formatting Principles

- **Font**: Arial throughout, size 21 (10.5pt) for body text
- **Page margins**: 1 inch all sides (1440 twips)
- **Headers**: Right-aligned, italic, gray, size 18
- **Footers**: Centered page numbers
- **Page breaks**: Before each major part
- **Spacing**: Paragraphs have `after: 100`, headings have larger before/after values
- **Callout boxes**: Implemented as single-cell tables with colored borders and shading
