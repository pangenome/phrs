#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

OUT_DIR="paper_prep/_brainstorming/fig5_extended_maternal_whole_genome_ribbons"
CLASS_WINNER_DIR="/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs"
FAI_DIR="/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs"
PLOTTER="paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft/scripts/plot_whole_genome_ribbon_draft.py"

python3 "$PLOTTER" \
  --comparison-id PAN027mat_vs_PAN010_joint \
  --class-winners "$CLASS_WINNER_DIR/PAN027mat_vs_PAN010_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz" \
  --query-fai "$FAI_DIR/PAN027mat_vs_PAN010_joint.query.fa.fai" \
  --target-fai "$FAI_DIR/PAN027mat_vs_PAN010_joint.target.fa.fai" \
  --output-dir "$OUT_DIR/PAN027mat_vs_PAN010_joint" \
  --query-label "PAN027 maternal" \
  --target-h1-label "PAN010 h1" \
  --target-h2-label "PAN010 h2" \
  --layer-label "mother-child" \
  --panel-label "C"

python3 "$PLOTTER" \
  --comparison-id PAN028mat_vs_PAN027_joint \
  --class-winners "$CLASS_WINNER_DIR/PAN028mat_vs_PAN027_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz" \
  --query-fai "$FAI_DIR/PAN028mat_vs_PAN027_joint.query.fa.fai" \
  --target-fai "$FAI_DIR/PAN028mat_vs_PAN027_joint.target.fa.fai" \
  --output-dir "$OUT_DIR/PAN028mat_vs_PAN027_joint" \
  --query-label "PAN028 maternal" \
  --target-h1-label "PAN027 h1" \
  --target-h2-label "PAN027 h2" \
  --layer-label "mother-child" \
  --panel-label "D"
