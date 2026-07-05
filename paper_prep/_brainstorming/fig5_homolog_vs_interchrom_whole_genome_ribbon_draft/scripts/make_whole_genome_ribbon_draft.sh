#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

# Labels match the published Fig5 whole-genome ribbon (PAN011 father h1 / PAN027
# pat child / PAN011 father h2). Inputs default to the repo-vendored data/ copies
# (self-contained); PDF/PNG need rsvg-convert (guix), else SVG only.
python3 paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft/scripts/plot_whole_genome_ribbon_draft.py \
  --query-label "PAN027 pat child" \
  --target-h1-label "PAN011 father h1" \
  --target-h2-label "PAN011 father h2"
