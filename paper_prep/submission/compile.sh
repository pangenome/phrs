#!/bin/bash
# Compile paper.tex to paper.pdf using pdflatex + biber
# Usage: bash compile.sh
set -euo pipefail

cd "$(dirname "$0")"

echo "==> Pass 1: pdflatex"
pdflatex -interaction=nonstopmode -file-line-error paper.tex

echo "==> bibtex"
bibtex paper

echo "==> Pass 2: pdflatex"
pdflatex -interaction=nonstopmode -file-line-error paper.tex

echo "==> Pass 3: pdflatex"
pdflatex -interaction=nonstopmode -file-line-error paper.tex

echo "==> Done: paper.pdf"
