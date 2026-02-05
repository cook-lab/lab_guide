#!/bin/bash
# Render a single markdown report to PDF with Cook Lab styling
#
# Usage:
#   ./render_report.sh reports/01_qc_report.md
#   ./render_report.sh reports/01_qc_report.md --open  # Open after rendering
#
# Requirements:
#   - Quarto (preferred) OR Pandoc with XeLaTeX
#   - Inter font installed (or modify font settings below)

set -e

# Check arguments
if [ $# -lt 1 ]; then
    echo "Usage: $0 <report.md> [--open]"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="${INPUT_FILE%.md}.pdf"
OPEN_AFTER=false

if [ "$2" = "--open" ]; then
    OPEN_AFTER=true
fi

# Check input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File not found: $INPUT_FILE"
    exit 1
fi

# Determine which renderer to use
if command -v quarto &> /dev/null; then
    echo "Rendering with Quarto: $INPUT_FILE → $OUTPUT_FILE"

    # Check for _quarto.yml in the report directory
    REPORT_DIR=$(dirname "$INPUT_FILE")
    if [ -f "$REPORT_DIR/_quarto.yml" ]; then
        quarto render "$INPUT_FILE" --to pdf
    else
        # Use inline Quarto options with system fonts
        quarto render "$INPUT_FILE" --to pdf \
            -M mainfont:"Helvetica Neue" \
            -M sansfont:"Helvetica Neue" \
            -M monofont:"Menlo" \
            -M fontsize:11pt \
            -M geometry:margin=1in \
            -M colorlinks:true \
            -M linkcolor:NavyBlue
    fi

elif command -v pandoc &> /dev/null; then
    echo "Rendering with Pandoc: $INPUT_FILE → $OUTPUT_FILE"

    pandoc "$INPUT_FILE" -o "$OUTPUT_FILE" \
        --pdf-engine=xelatex \
        -V geometry:margin=1in \
        -V mainfont="Helvetica Neue" \
        -V sansfont="Helvetica Neue" \
        -V monofont="Menlo" \
        -V fontsize=11pt \
        -V colorlinks=true \
        -V linkcolor=NavyBlue \
        -V urlcolor=NavyBlue \
        --highlight-style=github
else
    echo "Error: Neither Quarto nor Pandoc found. Please install one of them."
    echo "  - Quarto: https://quarto.org/docs/get-started/"
    echo "  - Pandoc: https://pandoc.org/installing.html"
    exit 1
fi

echo "✓ Generated: $OUTPUT_FILE"

# Open if requested
if [ "$OPEN_AFTER" = true ]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then
        open "$OUTPUT_FILE"
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        xdg-open "$OUTPUT_FILE"
    fi
fi
