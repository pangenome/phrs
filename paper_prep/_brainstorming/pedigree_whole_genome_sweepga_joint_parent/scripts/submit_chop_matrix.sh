#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
mkdir -p "$PACKAGE_DIR/logs" "$PACKAGE_DIR/summaries"

stdout="$PACKAGE_DIR/logs/chop_matrix.%j.out"
submit_out="$(sbatch --parsable \
    --export=ALL,PACKAGE_DIR="$PACKAGE_DIR",PAF_CHOP_LENGTH="${PAF_CHOP_LENGTH:-500000}",PAF_CHOP_OVERLAP="${PAF_CHOP_OVERLAP:-0}" \
    --job-name="wgsg_chop" \
    --cpus-per-task="${CHOP_CPUS:-2}" \
    --mem="${CHOP_MEM:-16G}" \
    --time="${CHOP_TIME:-04:00:00}" \
    --output="$stdout" \
    "$SCRIPT_DIR/run_chop_matrix.sh")"
job_id="${submit_out%%;*}"
printf "stage\tjob_id\tstdout\tchop_length_bp\toverlap_bp\nchop\t%s\t%s\t%s\t%s\n" "$job_id" "$stdout" "${PAF_CHOP_LENGTH:-500000}" "${PAF_CHOP_OVERLAP:-0}" > "$PACKAGE_DIR/summaries/chop_slurm_jobs.tsv"
echo "$job_id"
