#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

python3 "$SCRIPT_DIR/summarize_paf.py"

for cid in PAN027pat_vs_PAN011_joint PAN027mat_vs_PAN010_joint PAN028mat_vs_PAN027_joint; do
    test -s "$PACKAGE_DIR/raw_paf/${cid}.sweepga_many_many_j0.paf.gz"
    test -s "$PACKAGE_DIR/chopped_paf/${cid}.chopped_l${PAF_CHOP_LENGTH:-500000}_o${PAF_CHOP_OVERLAP:-0}.paf.gz"
    for filter_id in many_many_chopped one_one_chopped one_many_chopped two_many_chopped four_many_chopped; do
        test -s "$PACKAGE_DIR/filtered_paf/${cid}.${filter_id}.paf.gz"
    done
done

find "$PACKAGE_DIR/raw_paf" "$PACKAGE_DIR/chopped_paf" "$PACKAGE_DIR/filtered_paf" -name '*.paf.gz' -print0 | xargs -0 -r -n 1 gzip -t
test -s "$PACKAGE_DIR/summaries/input_manifest.tsv"
test -s "$PACKAGE_DIR/summaries/slurm_jobs.tsv"
test -s "$PACKAGE_DIR/summaries/chop_manifest.tsv"
test -s "$PACKAGE_DIR/summaries/filter_manifest.tsv"
test -s "$PACKAGE_DIR/summaries/paf_file_summary.tsv"
echo "Validated whole-genome sweepGA package outputs"
