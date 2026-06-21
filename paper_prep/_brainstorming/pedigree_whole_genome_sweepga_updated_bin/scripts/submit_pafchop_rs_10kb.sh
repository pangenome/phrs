#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$PACKAGE_DIR/../../.." && pwd)"
PAFCHOP_DIR="$REPO_ROOT/paper_prep/_brainstorming/pafchop-rs"
PAFCHOP_BIN="${PAFCHOP_BIN:-$REPO_ROOT/target/release/pafchop}"
CHOP_LENGTH="${PAF_CHOP_LENGTH:-10000}"
OVERLAP="${PAF_CHOP_OVERLAP:-0}"
THREADS="${PAFCHOP_THREADS_PER_JOB:-8}"
JOBS="${PAFCHOP_JOBS:-3}"

if [[ "$CHOP_LENGTH" != "10000" ]]; then
    echo "refusing to submit non-10kb chop length: $CHOP_LENGTH" >&2
    exit 2
fi
if [[ "$OVERLAP" != "0" ]]; then
    echo "refusing to submit nonzero overlap: $OVERLAP" >&2
    exit 2
fi

mkdir -p "$PACKAGE_DIR/logs" "$PACKAGE_DIR/summaries"
cargo test --manifest-path "$PAFCHOP_DIR/Cargo.toml"
cargo build --release --manifest-path "$PAFCHOP_DIR/Cargo.toml"
test -x "$PAFCHOP_BIN"

stdout="$PACKAGE_DIR/logs/pafchop_rs_l${CHOP_LENGTH}_o${OVERLAP}.%j.out"
submit_out="$(sbatch --parsable \
    --export=ALL,PACKAGE_DIR="$PACKAGE_DIR",PAFCHOP_DIR="$PAFCHOP_DIR",PAFCHOP_BIN="$PAFCHOP_BIN",PAF_CHOP_LENGTH="$CHOP_LENGTH",PAF_CHOP_OVERLAP="$OVERLAP",PAFCHOP_JOBS="$JOBS",PAFCHOP_THREADS_PER_JOB="$THREADS" \
    --job-name="wgsg_pafchop_10kb" \
    --cpus-per-task="${CHOP_CPUS:-24}" \
    --mem="${CHOP_MEM:-64G}" \
    --time="${CHOP_TIME:-04:00:00}" \
    --output="$stdout" \
    "$SCRIPT_DIR/run_pafchop_rs_10kb.sh")"
job_id="${submit_out%%;*}"
printf "stage\tjob_id\tstdout\tchop_length_bp\toverlap_bp\tpafchop_bin\tpafchop_sha256\tjobs\tthreads_per_job\npafchop_rs_10kb\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$job_id" "$stdout" "$CHOP_LENGTH" "$OVERLAP" "$PAFCHOP_BIN" "$(sha256sum "$PAFCHOP_BIN" | awk '{print $1}')" "$JOBS" "$THREADS" > "$PACKAGE_DIR/summaries/chop_slurm_jobs.tsv"
echo "$job_id"
