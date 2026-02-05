# Makefile for rendering Cook Lab analysis reports
# Copy this file to your project's reports/ directory
#
# Usage:
#   make                  - Render all .md files to PDF
#   make 01_qc_report.pdf - Render specific report
#   make clean            - Remove all generated PDFs
#   make list             - List available reports

# Find all markdown files (excluding README and templates)
MD_FILES := $(filter-out README.md _%.md, $(wildcard *.md))
PDF_FILES := $(MD_FILES:.md=.pdf)

# Default target: render all reports
all: $(PDF_FILES)

# Pattern rule: convert .md to .pdf using Quarto
%.pdf: %.md _quarto.yml
	@echo "Rendering $< → $@"
	@quarto render $< --to pdf --quiet
	@echo "✓ Generated $@"

# Alternative: use Pandoc if Quarto not available
# Uncomment below and comment out the Quarto rule above
# %.pdf: %.md
# 	@echo "Rendering $< → $@"
# 	@pandoc $< -o $@ \
# 		--pdf-engine=xelatex \
# 		-V geometry:margin=1in \
# 		-V mainfont="Inter" \
# 		-V sansfont="Inter" \
# 		-V monofont="JetBrains Mono" \
# 		-V fontsize=11pt \
# 		-V colorlinks=true \
# 		-V linkcolor=NavyBlue \
# 		--highlight-style=github
# 	@echo "✓ Generated $@"

# List available reports
list:
	@echo "Available reports:"
	@for f in $(MD_FILES); do echo "  - $$f"; done
	@echo ""
	@echo "Generated PDFs:"
	@for f in $(wildcard *.pdf); do echo "  - $$f"; done

# Clean generated PDFs
clean:
	@echo "Removing generated PDFs..."
	@rm -f $(PDF_FILES)
	@echo "✓ Clean complete"

# Watch for changes and auto-render (requires fswatch or entr)
watch:
	@echo "Watching for changes... (Ctrl+C to stop)"
	@fswatch -o $(MD_FILES) | xargs -n1 -I{} make all

# Help
help:
	@echo "Cook Lab Report Rendering"
	@echo ""
	@echo "Usage:"
	@echo "  make              Render all reports to PDF"
	@echo "  make <name>.pdf   Render specific report"
	@echo "  make list         List available reports"
	@echo "  make clean        Remove all PDFs"
	@echo "  make watch        Auto-render on file changes"
	@echo ""
	@echo "Requirements:"
	@echo "  - Quarto (preferred) or Pandoc with XeLaTeX"
	@echo "  - Inter font (or edit _quarto.yml for alternatives)"

.PHONY: all list clean watch help
