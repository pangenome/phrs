#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
mkdir -p "$PACKAGE_DIR/logs" "$PACKAGE_DIR/summaries"

stdout="$PACKAGE_DIR/logs/filter_matrix.%j.out"
DEVSHM_BASE="${SWEEPGA_DEVSHM_BASE:-/dev/shm}"
submit_out="$(sbatch --parsable \
    --export=ALL,PACKAGE_DIR="$PACKAGE_DIR",PAF_CHOP_LENGTH="${PAF_CHOP_LENGTH:-500000}",PAF_CHOP_OVERLAP="${PAF_CHOP_OVERLAP:-0}",SWEEPGA_DEVSHM_BASE="$DEVSHM_BASE" \
    --job-name="wgsg_filter" \
    --cpus-per-task="${FILTER_CPUS:-4}" \
    --mem="${FILTER_MEM:-32G}" \
    --time="${FILTER_TIME:-12:00:00}" \
    --output="$stdout" \
    "$SCRIPT_DIR/run_filter_matrix.sh")"
job_id="${submit_out%%;*}"
printf "stage\tjob_id\tsweepga_devshm_base\tstdout\nfilter\t%s\t%s\t%s\n" "$job_id" "$DEVSHM_BASE" "$stdout" > "$PACKAGE_DIR/summaries/filter_slurm_jobs.tsv"
echo "$job_id"
