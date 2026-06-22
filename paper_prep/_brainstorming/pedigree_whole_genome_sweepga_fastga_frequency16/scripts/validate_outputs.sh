#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

test -s "$PACKAGE_DIR/summaries/sweepga_binary.tsv"
test -s "$PACKAGE_DIR/summaries/fastga_binary.tsv"
test -s "$PACKAGE_DIR/summaries/slurm_jobs.tsv"
test -s "$PACKAGE_DIR/summaries/raw_chr3_support.tsv"
test -s "$PACKAGE_DIR/summaries/frequency_sensitivity_summary.tsv"
test -s "$PACKAGE_DIR/summaries/input_manifest.tsv"

grep -R -- '--fastga-frequency 16' "$PACKAGE_DIR/logs" >/dev/null
grep -R -- '--scaffold-jump 0' "$PACKAGE_DIR/logs" >/dev/null
grep -R -- '/dev/shm' "$PACKAGE_DIR/logs" >/dev/null
grep -R -- '/home/erikg/.cargo/bin/sweepga' "$PACKAGE_DIR/logs" >/dev/null
grep -R -- 'FastGA -1:.* -f16 ' "$PACKAGE_DIR/logs" >/dev/null
grep -q $'PAN027_chr9q_chr3q_PHR_candidate\tPAN027pat_vs_PAN011_joint' "$PACKAGE_DIR/summaries/raw_chr3_support.tsv"
grep -q $'PAN028_chr9q_chr3q_PHR_candidate\tPAN028mat_vs_PAN027_joint' "$PACKAGE_DIR/summaries/raw_chr3_support.tsv"
if grep -q 'PATHOLOGICAL_NO_RAW_PAF' "$PACKAGE_DIR/summaries/slurm_jobs.tsv"; then
    test -s "$PACKAGE_DIR/summaries/pathological_runtime.tsv"
else
    test -s "$PACKAGE_DIR/summaries/paf_file_summary.tsv"
    n_raw="$(find "$PACKAGE_DIR/raw_paf" -type f -name '*.paf.gz' | wc -l)"
    if [[ "$n_raw" -ne 3 ]]; then
        echo "expected 3 raw PAFs, found $n_raw" >&2
        exit 1
    fi
fi
if grep -q $'chr3\t[1-9][0-9]*\t' "$PACKAGE_DIR/summaries/raw_chr3_support.tsv"; then
    test -s "$PACKAGE_DIR/summaries/chop_manifest.tsv"
    test -s "$PACKAGE_DIR/summaries/filter_manifest.tsv"
    grep -q 'four_many_chopped' "$PACKAGE_DIR/summaries/filter_manifest.tsv"
else
    if [[ -e "$PACKAGE_DIR/summaries/chop_manifest.tsv" || -e "$PACKAGE_DIR/summaries/filter_manifest.tsv" ]]; then
        echo "chop/filter manifests present despite no raw chr3 support" >&2
        exit 1
    fi
fi
if git -C "$PACKAGE_DIR/../../.." status --short -- submission | grep -q .; then
    echo "submission/ was modified" >&2
    exit 1
fi

echo "Validated frequency16 sweepGA/FastGA package outputs"
