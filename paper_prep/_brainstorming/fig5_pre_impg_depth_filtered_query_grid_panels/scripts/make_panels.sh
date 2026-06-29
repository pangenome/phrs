#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

python3 paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_query_grid_panels/scripts/extract_pre_impg_query_grid_segments.py
/home/erikg/.guix-profile/bin/Rscript paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_query_grid_panels/scripts/plot_pre_impg_query_grid_panel.R \
  paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_query_grid_panels

