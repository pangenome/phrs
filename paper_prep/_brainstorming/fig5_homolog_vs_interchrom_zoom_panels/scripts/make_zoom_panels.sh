#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

out_dir="paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels"
class_winners="paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs/PAN027pat_vs_PAN011_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz"

if [[ -r "$class_winners" ]]; then
  python3 "$out_dir/scripts/extract_zoom_windows.py"
elif [[ -r "$out_dir/zoom_window_segments.tsv" && -r "$out_dir/zoom_panel_summary.tsv" ]]; then
  echo "Using committed zoom panel TSV snapshots; upstream class-winner gzip not found: $class_winners" >&2
else
  echo "Missing upstream class-winner gzip and committed zoom panel TSV snapshots" >&2
  exit 1
fi

/home/erikg/.guix-profile/bin/Rscript \
  "$out_dir/scripts/plot_zoom_panels.R"
