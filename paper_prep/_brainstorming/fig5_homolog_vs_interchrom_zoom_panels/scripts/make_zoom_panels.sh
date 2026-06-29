#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

python3 paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/scripts/extract_zoom_windows.py
/home/erikg/.guix-profile/bin/Rscript \
  paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/scripts/plot_zoom_panels.R

