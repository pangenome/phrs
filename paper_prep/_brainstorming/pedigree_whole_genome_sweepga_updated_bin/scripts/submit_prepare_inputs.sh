#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
mkdir -p "$PACKAGE_DIR/logs" "$PACKAGE_DIR/summaries"

stdout="$PACKAGE_DIR/logs/prepare_inputs.%j.out"
submit_out="$(sbatch --parsable \
    --export=ALL,PACKAGE_DIR="$PACKAGE_DIR" \
    --job-name="wgsg_prepare" \
    --cpus-per-task="${PREPARE_CPUS:-4}" \
    --mem="${PREPARE_MEM:-32G}" \
    --time="${PREPARE_TIME:-08:00:00}" \
    --output="$stdout" \
    "$SCRIPT_DIR/run_prepare_inputs.sh")"
job_id="${submit_out%%;*}"
printf "stage\tjob_id\tstdout\nprepare_inputs\t%s\t%s\n" "$job_id" "$stdout" > "$PACKAGE_DIR/summaries/prepare_slurm_jobs.tsv"
echo "$job_id"
