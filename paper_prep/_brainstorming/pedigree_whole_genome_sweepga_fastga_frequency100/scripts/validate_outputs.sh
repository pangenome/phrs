#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

test -s "$PACKAGE_DIR/summaries/sweepga_binary.tsv"
test -s "$PACKAGE_DIR/summaries/fastga_binary.tsv"
test -s "$PACKAGE_DIR/summaries/slurm_jobs.tsv"
test -s "$PACKAGE_DIR/summaries/raw_chr3_support.tsv"
test -s "$PACKAGE_DIR/summaries/frequency_sensitivity_summary.tsv"
test -s "$PACKAGE_DIR/summaries/pathological_runtime.tsv"

grep -R -- '--fastga-frequency 100' "$PACKAGE_DIR/logs" >/dev/null
grep -R -- '--scaffold-jump 0' "$PACKAGE_DIR/logs" >/dev/null
grep -R -- '/dev/shm' "$PACKAGE_DIR/logs" >/dev/null
grep -R -- '/home/erikg/.cargo/bin/sweepga' "$PACKAGE_DIR/logs" >/dev/null
grep -R -- 'FastGA -1:.* -f100 ' "$PACKAGE_DIR/logs" >/dev/null
grep -q 'PATHOLOGICAL_NO_RAW_PAF' "$PACKAGE_DIR/summaries/slurm_jobs.tsv"
grep -q 'PATHOLOGICAL_NO_RAW_PAF' "$PACKAGE_DIR/summaries/raw_chr3_support.tsv"
grep -q 'PATHOLOGICAL_NO_RAW_PAF' "$PACKAGE_DIR/summaries/frequency_sensitivity_summary.tsv"
grep -q $'PAN027_chr9q_chr3q_PHR_candidate\tPAN027pat_vs_PAN011_joint' "$PACKAGE_DIR/summaries/raw_chr3_support.tsv"
grep -q $'PAN028_chr9q_chr3q_PHR_candidate\tPAN028mat_vs_PAN027_joint' "$PACKAGE_DIR/summaries/raw_chr3_support.tsv"
if find "$PACKAGE_DIR/raw_paf" -type f -name '*.paf.gz' | grep -q .; then
    echo "unexpected raw PAF present after pathological frequency100 run" >&2
    exit 1
fi

echo "Validated frequency100 sweepGA/FastGA package outputs"
