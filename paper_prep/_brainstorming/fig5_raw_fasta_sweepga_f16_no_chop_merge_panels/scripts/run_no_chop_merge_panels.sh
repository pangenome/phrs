#!/usr/bin/env bash
set -euo pipefail

panel_dir="${1:-paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_no_chop_merge_panels}"
source_package="${2:-/moosefs/erikg/phrs/.wg-worktrees/agent-2649/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16}"
sweepga_bin="${SWEEPGA:-/home/erikg/.cargo/bin/sweepga}"
scratch_base="${SWEEPGA_DEVSHM_BASE:-${TMPDIR:-/tmp}}"

mkdir -p "$panel_dir/filtered_paf" "$panel_dir/logs" "$panel_dir/work"

python3 "$panel_dir/scripts/run_raw_no_chop_filters.py" \
  --panel-dir "$panel_dir" \
  --source-package "$source_package" \
  --sweepga "$sweepga_bin" \
  --scratch-base "$scratch_base"

python3 "$panel_dir/scripts/extract_raw_no_chop_merge_segments.py" \
  --panel-dir "$panel_dir" \
  --source-package "$source_package"

Rscript "$panel_dir/scripts/plot_raw_no_chop_merge_panel.R" "$panel_dir"
