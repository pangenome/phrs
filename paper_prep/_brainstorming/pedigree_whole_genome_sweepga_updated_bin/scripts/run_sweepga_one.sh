#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: $0 <comparison_id> [num_mappings] [scaffold_jump]" >&2
    exit 2
fi

COMPARISON_ID="$1"
NUM_MAPPINGS="${2:-many:many}"
SCAFFOLD_JUMP="${3:-0}"
if [[ -n "${PACKAGE_DIR:-}" ]]; then
    SCRIPT_DIR="$PACKAGE_DIR/scripts"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
fi
SWEEPGA="${SWEEPGA:-/home/erikg/.cargo/bin/sweepga}"
SWEEPGA_WHICH="$(command -v sweepga || true)"
SWEEPGA_REALPATH="$(realpath "$SWEEPGA")"
SWEEPGA_SHA256="$(sha256sum "$SWEEPGA" | awk '{print $1}')"

QUERY="$PACKAGE_DIR/inputs/${COMPARISON_ID}.query.fa"
TARGET="$PACKAGE_DIR/inputs/${COMPARISON_ID}.target.fa"
OUT="$PACKAGE_DIR/raw_paf/${COMPARISON_ID}.sweepga_${NUM_MAPPINGS/:/_}_j${SCAFFOLD_JUMP}.paf.gz"
LOG="$PACKAGE_DIR/logs/${COMPARISON_ID}.sweepga_${NUM_MAPPINGS/:/_}_j${SCAFFOLD_JUMP}.${SLURM_JOB_ID:-manual}.log"

if [[ ! -s "$QUERY" || ! -s "$TARGET" ]]; then
    echo "missing full whole-genome input FASTA for $COMPARISON_ID; run scripts/prepare_inputs.py via Slurm first" >&2
    exit 1
fi

SCRATCH_BASE="${SWEEPGA_DEVSHM_BASE:-/dev/shm}"
if [[ ! -d "$SCRATCH_BASE" || ! -w "$SCRATCH_BASE" ]]; then
    echo "required sweepGA/FastGA scratch base is not writable: $SCRATCH_BASE" >&2
    exit 1
fi
TMPDIR="$(mktemp -d "$SCRATCH_BASE/sg.${SLURM_JOB_ID:-manual}.XXXXXX")"
cleanup() {
    rm -rf "$TMPDIR"
}
trap cleanup EXIT
export TMPDIR
STAGED_QUERY="$TMPDIR/q.fa"
STAGED_TARGET="$TMPDIR/t.fa"
cp "$QUERY" "$STAGED_QUERY"
cp "$TARGET" "$STAGED_TARGET"

{
    date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
    echo "host=$(hostname)"
    echo "job_id=${SLURM_JOB_ID:-manual}"
    echo "comparison_id=$COMPARISON_ID"
    echo "coordinate_scope=full whole-genome assembly haplotype FASTA records"
    echo "scratch=$TMPDIR"
    echo "scratch_policy=sweepGA/FastGA graph/database/source-adjacent temporaries are explicitly placed under /dev/shm; SLURM_TMPDIR is not used for sweepGA scratch"
    echo "scratch_fastga_sources=$STAGED_QUERY,$STAGED_TARGET"
    echo "sweepga=$SWEEPGA"
    echo "sweepga_which=$SWEEPGA_WHICH"
    echo "sweepga_realpath=$SWEEPGA_REALPATH"
    echo "sweepga_sha256=$SWEEPGA_SHA256"
    "$SWEEPGA" --version || true
    "$SWEEPGA" --help | sed -n '1,80p' || true
    TMP_PAF="$TMPDIR/${COMPARISON_ID}.paf"
    echo "command=$SWEEPGA --fastga --num-mappings $NUM_MAPPINGS --scaffold-jump $SCAFFOLD_JUMP --temp-dir $TMPDIR --output-file $TMP_PAF $STAGED_QUERY $STAGED_TARGET"
    "$SWEEPGA" --fastga \
        --num-mappings "$NUM_MAPPINGS" \
        --scaffold-jump "$SCAFFOLD_JUMP" \
        --temp-dir "$TMPDIR" \
        --output-file "$TMP_PAF" \
        "$STAGED_QUERY" "$STAGED_TARGET"
    gzip -c "$TMP_PAF" > "$OUT"
    gzip -t "$OUT"
    date -u '+finished_utc=%Y-%m-%dT%H:%M:%SZ'
    echo "output=$OUT"
    echo "bytes=$(wc -c < "$OUT")"
} 2>&1 | tee "$LOG"
