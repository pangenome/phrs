#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: $0 PACKAGE_DIR [CHOP_LENGTH=10000]" >&2
    exit 2
fi

PACKAGE_DIR="$(cd "$1" && pwd)"
CHOP_LENGTH="${2:-10000}"
JOBS="${PAFCHOP_JOBS:-3}"
THREADS_PER_JOB="${PAFCHOP_THREADS_PER_JOB:-8}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CRATE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="${PAFCHOP_BIN:-$CRATE_DIR/target/release/pafchop}"

cargo build --release --manifest-path "$CRATE_DIR/Cargo.toml"
mkdir -p "$PACKAGE_DIR/chopped_paf_l${CHOP_LENGTH}_o0" "$PACKAGE_DIR/summaries/pafchop_l${CHOP_LENGTH}_o0"

export PACKAGE_DIR CHOP_LENGTH THREADS_PER_JOB SCRIPT_DIR BIN
tail -n +2 "$PACKAGE_DIR/config/comparisons.tsv" | cut -f1 \
    | xargs -r -n 1 -P "$JOBS" bash -c '
        set -euo pipefail
        cid="$1"
        raw="$PACKAGE_DIR/raw_paf/${cid}.sweepga_many_many_j0.paf.gz"
        out="$PACKAGE_DIR/chopped_paf_l${CHOP_LENGTH}_o0/${cid}.chopped_l${CHOP_LENGTH}_o0.paf.gz"
        summary="$PACKAGE_DIR/summaries/pafchop_l${CHOP_LENGTH}_o0/${cid}.summary.tsv"
        "$SCRIPT_DIR/chop_one.sh" "$raw" "$out" "$summary" "$cid" "$CHOP_LENGTH" "$THREADS_PER_JOB"
    ' _

first_summary="$(find "$PACKAGE_DIR/summaries/pafchop_l${CHOP_LENGTH}_o0" -name '*.summary.tsv' | sort | head -n 1)"
{
    head -n 1 "$first_summary"
    for f in "$PACKAGE_DIR/summaries/pafchop_l${CHOP_LENGTH}_o0/"*.summary.tsv; do
        tail -n +2 "$f"
    done
} > "$PACKAGE_DIR/summaries/chop_manifest_l${CHOP_LENGTH}_o0.tsv"
