#!/usr/bin/env bash
set -euo pipefail

if [[ -n "${PACKAGE_DIR:-}" ]]; then
    SCRIPT_DIR="$PACKAGE_DIR/scripts"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
fi
mkdir -p "$PACKAGE_DIR/chopped_paf" "$PACKAGE_DIR/summaries" "$PACKAGE_DIR/logs"

{
    date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
    echo "host=$(hostname)"
    echo "job_id=${SLURM_JOB_ID:-manual}"
    echo "chop_length=${PAF_CHOP_LENGTH:-500000}"
    echo "overlap=${PAF_CHOP_OVERLAP:-0}"
    python3 "$SCRIPT_DIR/chop_paf.py" --chop-length "${PAF_CHOP_LENGTH:-500000}" --overlap "${PAF_CHOP_OVERLAP:-0}"
    find "$PACKAGE_DIR/chopped_paf" -name '*.paf.gz' -print0 | xargs -0 -r -n 1 gzip -t
    date -u '+finished_utc=%Y-%m-%dT%H:%M:%SZ'
} 2>&1 | tee "$PACKAGE_DIR/logs/chop_matrix.${SLURM_JOB_ID:-manual}.log"
