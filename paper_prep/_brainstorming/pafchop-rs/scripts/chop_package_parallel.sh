#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: $0 PACKAGE_DIR [CHOP_LENGTH=10000] [CHOP_MODE=row-start|query-grid]" >&2
    exit 2
fi

PACKAGE_DIR="$(cd "$1" && pwd)"
CHOP_LENGTH="${2:-10000}"
CHOP_MODE="${3:-${PAF_CHOP_MODE:-row-start}}"
JOBS="${PAFCHOP_JOBS:-3}"
THREADS_PER_JOB="${PAFCHOP_THREADS_PER_JOB:-8}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CRATE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="${PAFCHOP_BIN:-$CRATE_DIR/target/release/pafchop}"
case "$CHOP_MODE" in
    row-start|query-grid) ;;
    *) echo "invalid CHOP_MODE: $CHOP_MODE" >&2; exit 2 ;;
esac
MODE_SUFFIX=""
if [[ "$CHOP_MODE" == "query-grid" ]]; then
    MODE_SUFFIX="_query_grid"
fi

cargo build --release --manifest-path "$CRATE_DIR/Cargo.toml"
mkdir -p "$PACKAGE_DIR/chopped_paf_l${CHOP_LENGTH}_o0${MODE_SUFFIX}" "$PACKAGE_DIR/summaries/pafchop_l${CHOP_LENGTH}_o0${MODE_SUFFIX}"

export PACKAGE_DIR CHOP_LENGTH CHOP_MODE THREADS_PER_JOB SCRIPT_DIR BIN MODE_SUFFIX
tail -n +2 "$PACKAGE_DIR/config/comparisons.tsv" | cut -f1 \
    | xargs -r -n 1 -P "$JOBS" bash -c '
        set -euo pipefail
        cid="$1"
        raw="$PACKAGE_DIR/raw_paf/${cid}.sweepga_many_many_j0.paf.gz"
        out="$PACKAGE_DIR/chopped_paf_l${CHOP_LENGTH}_o0${MODE_SUFFIX}/${cid}.chopped_l${CHOP_LENGTH}_o0${MODE_SUFFIX}.paf.gz"
        summary="$PACKAGE_DIR/summaries/pafchop_l${CHOP_LENGTH}_o0${MODE_SUFFIX}/${cid}.summary.tsv"
        "$SCRIPT_DIR/chop_one.sh" "$raw" "$out" "$summary" "$cid" "$CHOP_LENGTH" "$THREADS_PER_JOB" "$CHOP_MODE"
    ' _

first_summary="$(find "$PACKAGE_DIR/summaries/pafchop_l${CHOP_LENGTH}_o0${MODE_SUFFIX}" -name '*.summary.tsv' | sort | head -n 1)"
{
    head -n 1 "$first_summary"
    for f in "$PACKAGE_DIR/summaries/pafchop_l${CHOP_LENGTH}_o0${MODE_SUFFIX}/"*.summary.tsv; do
        tail -n +2 "$f"
    done
} > "$PACKAGE_DIR/summaries/chop_manifest_l${CHOP_LENGTH}_o0${MODE_SUFFIX}.tsv"
