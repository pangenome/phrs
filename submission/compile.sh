#!/usr/bin/env bash
set -euxo pipefail

: "${PAPER:=paper}"
: "${LATEX:=pdflatex}"
: "${BIBTEX:=bibtex}"
: "${LATEX_FLAGS:=-interaction=nonstopmode -halt-on-error}"

rm -f "$PAPER.bbl" "$PAPER.blg" Meth.bbl Meth.blg Supp.bbl Supp.blg

echo "=== Pass 1: $LATEX ==="
"$LATEX" $LATEX_FLAGS "$PAPER"

echo "=== bibtex: main bibliography ==="
"$BIBTEX" "$PAPER"

echo "=== bibtex: Meth bibliography ==="
"$BIBTEX" Meth

echo "=== bibtex: Supp bibliography ==="
"$BIBTEX" Supp || true   # Supp may be empty; tolerate failure

echo "=== Pass 2: $LATEX ==="
"$LATEX" $LATEX_FLAGS "$PAPER"

echo "=== Pass 3: $LATEX ==="
"$LATEX" $LATEX_FLAGS "$PAPER"

echo "=== Done: $PAPER.pdf ==="
