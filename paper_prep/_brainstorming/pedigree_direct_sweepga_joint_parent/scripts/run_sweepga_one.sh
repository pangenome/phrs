#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: $0 <comparison_id> [num_mappings] [scaffold_jump]" >&2
    exit 2
fi

COMPARISON_ID="$1"
NUM_MAPPINGS="${2:-many:many}"
SCAFFOLD_JUMP="${3:-0}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="${PACKAGE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
SWEEPGA="${SWEEPGA:-/home/erikg/.cargo/bin/sweepga}"

QUERY="$PACKAGE_DIR/inputs/${COMPARISON_ID}.query.fa"
TARGET="$PACKAGE_DIR/inputs/${COMPARISON_ID}.target.fa"
OUT="$PACKAGE_DIR/raw_paf/${COMPARISON_ID}.sweepga_${NUM_MAPPINGS/:/_}_j${SCAFFOLD_JUMP}.paf.gz"
LOG="$PACKAGE_DIR/logs/${COMPARISON_ID}.sweepga_${NUM_MAPPINGS/:/_}_j${SCAFFOLD_JUMP}.log"

if [[ ! -s "$QUERY" || ! -s "$TARGET" ]]; then
    echo "missing input FASTA for $COMPARISON_ID; run scripts/prepare_inputs.py first" >&2
    exit 1
fi

SCRATCH_BASE="${SLURM_TMPDIR:-/tmp}"
TMPDIR="$(mktemp -d "$SCRATCH_BASE/sweepga.${COMPARISON_ID}.${SLURM_JOB_ID:-manual}.XXXXXX")"
cleanup() {
    rm -rf "$TMPDIR"
}
trap cleanup EXIT
export TMPDIR

{
    date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
    echo "host=$(hostname)"
    echo "job_id=${SLURM_JOB_ID:-manual}"
    echo "comparison_id=$COMPARISON_ID"
    echo "sweepga=$SWEEPGA"
    "$SWEEPGA" --version || true
    TMP_PAF="$TMPDIR/${COMPARISON_ID}.paf"
    echo "coordinate_scope=500kb_telomeric_source_windows; paf_offsets_are_local_0_500kb_window_coordinates"
    echo "command=$SWEEPGA --fastga --num-mappings $NUM_MAPPINGS --scaffold-jump $SCAFFOLD_JUMP --temp-dir $TMPDIR --output-file $TMP_PAF $QUERY $TARGET"
    "$SWEEPGA" --fastga \
        --num-mappings "$NUM_MAPPINGS" \
        --scaffold-jump "$SCAFFOLD_JUMP" \
        --temp-dir "$TMPDIR" \
        --output-file "$TMP_PAF" \
        "$QUERY" "$TARGET"
    gzip -c "$TMP_PAF" > "$OUT"
    date -u '+finished_utc=%Y-%m-%dT%H:%M:%SZ'
    echo "output=$OUT"
    echo "bytes=$(wc -c < "$OUT")"
} 2>&1 | tee "$LOG"
