#!/usr/bin/env bash
set -euo pipefail

if [[ -n "${PACKAGE_DIR:-}" ]]; then
    SCRIPT_DIR="$PACKAGE_DIR/scripts"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
fi

mkdir -p "$PACKAGE_DIR/logs" "$PACKAGE_DIR/summaries" "$PACKAGE_DIR/inputs"
{
    date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
    echo "host=$(hostname)"
    echo "job_id=${SLURM_JOB_ID:-manual}"
    echo "package_dir=$PACKAGE_DIR"
    echo "scope=full whole-genome source FASTAs; no 500kb window FASTA"
    python3 "$SCRIPT_DIR/prepare_inputs.py"
    date -u '+finished_utc=%Y-%m-%dT%H:%M:%SZ'
} 2>&1 | tee "$PACKAGE_DIR/logs/prepare_inputs.${SLURM_JOB_ID:-manual}.log"
