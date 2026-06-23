#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/moosefs/erikg/phrs}"
PACKAGE_DIR="${PACKAGE_DIR:-$REPO_ROOT/paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_scaffold_jump_sensitivity}"

"$PACKAGE_DIR/scripts/build_filter_tasks.py"
"$PACKAGE_DIR/scripts/summarize_candidate_windows.py"
Rscript "$PACKAGE_DIR/scripts/plot_chr3_heatmap.R" "$PACKAGE_DIR"

printf 'Validation complete:\n'
printf '  tasks: %s\n' "$(tail -n +2 "$PACKAGE_DIR/filter_tasks.tsv" | wc -l)"
printf '  summary rows: %s\n' "$(tail -n +2 "$PACKAGE_DIR/candidate_window_summary.tsv" | wc -l)"
printf '  manifest rows: %s\n' "$(tail -n +2 "$PACKAGE_DIR/filtered_paf_manifest.tsv" | wc -l)"
