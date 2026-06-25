#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"
SCAN_DIR="paper_prep/_brainstorming/fig5_whole_genome_existing_paf_impg_like_scan"

python3 "$SCAN_DIR/scripts/build_paf_input_manifest.py" \
  --bin-size 2000 \
  --out "$SCAN_DIR/manifests/paf_inputs.tsv"

python3 -m py_compile \
  "$SCAN_DIR/scripts/build_paf_input_manifest.py" \
  "$SCAN_DIR/scripts/summarize_existing_paf_impg_like_scan.py" \
  "$SCAN_DIR/scripts/write_report.py"

test -s "$SCAN_DIR/manifests/paf_inputs.tsv"
