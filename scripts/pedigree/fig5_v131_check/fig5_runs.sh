#!/bin/bash
# Fig5 drawn runs from class winners
# args: CLASS_WINNERS_TSV_GZ, CID, FAI_DIR, OUT_DIR
# Runs the unchanged submission/scripts/figures/fig5_ribbon_plotter.py (copy at /tmp/guarracino.fig5_ribbon_plotter.py,
# only ROOT parents[4] -> parents[1] changed). Writes OUT_DIR/CID.whole_genome_ribbon_runs.tsv (+ svg, summary).
set -euo pipefail
CW=$1; CID=$2; FAI=$3; OUT=$4
P=/tmp/guarracino.fig5_ribbon_plotter.py
case $CID in
  PAN027pat_*) L=("PAN027 paternal" "PAN011 h1" "PAN011 h2");;
  PAN027mat_*) L=("PAN027 maternal" "PAN010 h1" "PAN010 h2");;
  PAN028mat_*) L=("PAN028 maternal" "PAN027 h1" "PAN027 h2");;
  *) echo "unknown CID $CID" >&2; exit 1;;
esac
[ -s "$CW" ] || { echo "missing $CW" >&2; exit 1; }
mkdir -p "$OUT"
echo "[$(date -Is)] fig5 runs $CID from $CW"
FIG5_HIDE_LEGEND=1 python3 $P --comparison-id $CID --class-winners "$CW" \
  --query-fai $FAI/$CID.query.fa.fai --target-fai $FAI/$CID.target.fa.fai --output-dir "$OUT" \
  --query-label "${L[0]}" --target-h1-label "${L[1]}" --target-h2-label "${L[2]}" --panel-label "" 2>&1 | tail -5
R=$OUT/$CID.whole_genome_ribbon_runs.tsv
echo "[$(date -Is)] drawn runs: $(($(wc -l < $R) - 1)) in $R"
