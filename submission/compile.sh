#!/usr/bin/env bash
set -euxo pipefail

PAPER=paper

echo "=== Pass 1: pdflatex ==="
pdflatex -interaction=nonstopmode "$PAPER"

echo "=== bibtex: main bibliography ==="
bibtex "$PAPER"

echo "=== bibtex: Meth bibliography ==="
bibtex Meth

echo "=== bibtex: Supp bibliography ==="
bibtex Supp || true   # Supp may be empty; tolerate failure

echo "=== Pass 2: pdflatex ==="
pdflatex -interaction=nonstopmode "$PAPER"

echo "=== Pass 3: pdflatex ==="
pdflatex -interaction=nonstopmode "$PAPER"

echo "=== Done: $PAPER.pdf ==="
