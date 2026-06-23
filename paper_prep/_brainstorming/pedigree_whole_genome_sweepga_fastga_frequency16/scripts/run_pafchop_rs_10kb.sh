#!/usr/bin/env bash
set -euo pipefail

if [[ -n "${PACKAGE_DIR:-}" ]]; then
    SCRIPT_DIR="$PACKAGE_DIR/scripts"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
fi
REPO_ROOT="$(cd "$PACKAGE_DIR/../../.." && pwd)"
PAFCHOP_DIR="${PAFCHOP_DIR:-$REPO_ROOT/paper_prep/_brainstorming/pafchop-rs}"
PAFCHOP_BIN="${PAFCHOP_BIN:-$REPO_ROOT/target/release/pafchop}"
CHOP_LENGTH="${PAF_CHOP_LENGTH:-10000}"
OVERLAP="${PAF_CHOP_OVERLAP:-0}"
CHOP_MODE="${PAF_CHOP_MODE:-query-grid}"
JOBS="${PAFCHOP_JOBS:-3}"
THREADS="${PAFCHOP_THREADS_PER_JOB:-8}"
case "$CHOP_MODE" in
    row-start|query-grid) ;;
    *) echo "invalid PAF_CHOP_MODE: $CHOP_MODE" >&2; exit 2 ;;
esac
MODE_SUFFIX=""
if [[ "$CHOP_MODE" == "query-grid" ]]; then
    MODE_SUFFIX="_query_grid"
fi

{
    date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
    echo "host=$(hostname)"
    echo "job_id=${SLURM_JOB_ID:-manual}"
    echo "package_dir=$PACKAGE_DIR"
    echo "pafchop_dir=$PAFCHOP_DIR"
    echo "pafchop_bin=$PAFCHOP_BIN"
    echo "pafchop_sha256=$(sha256sum "$PAFCHOP_BIN" | awk '{print $1}')"
    echo "chop_length=$CHOP_LENGTH"
    echo "overlap=$OVERLAP"
    echo "chunk_mode=$CHOP_MODE"
    echo "jobs=$JOBS"
    echo "threads_per_job=$THREADS"
    echo "command=PAFCHOP_BIN=$PAFCHOP_BIN PAF_CHOP_MODE=$CHOP_MODE PAFCHOP_JOBS=$JOBS PAFCHOP_THREADS_PER_JOB=$THREADS $PAFCHOP_DIR/scripts/chop_package_parallel.sh $PACKAGE_DIR $CHOP_LENGTH $CHOP_MODE"
    PAFCHOP_BIN="$PAFCHOP_BIN" PAFCHOP_JOBS="$JOBS" PAFCHOP_THREADS_PER_JOB="$THREADS" \
        "$PAFCHOP_DIR/scripts/chop_package_parallel.sh" "$PACKAGE_DIR" "$CHOP_LENGTH" "$CHOP_MODE"
    cp "$PACKAGE_DIR/summaries/chop_manifest_l${CHOP_LENGTH}_o${OVERLAP}${MODE_SUFFIX}.tsv" "$PACKAGE_DIR/summaries/chop_manifest${MODE_SUFFIX}.tsv"
    date -u '+finished_utc=%Y-%m-%dT%H:%M:%SZ'
} 2>&1 | tee "$PACKAGE_DIR/logs/pafchop_rs_l${CHOP_LENGTH}_o${OVERLAP}${MODE_SUFFIX}.${SLURM_JOB_ID:-manual}.log"
