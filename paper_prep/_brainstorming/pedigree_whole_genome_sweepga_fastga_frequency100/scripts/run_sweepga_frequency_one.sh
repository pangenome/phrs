#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: $0 <comparison_id> [fastga_frequency] [num_mappings] [scaffold_jump]" >&2
    exit 2
fi

COMPARISON_ID="$1"
FASTGA_FREQUENCY="${2:-100}"
NUM_MAPPINGS="${3:-many:many}"
SCAFFOLD_JUMP="${4:-0}"

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
LABEL="frequency${FASTGA_FREQUENCY}_many_many_j${SCAFFOLD_JUMP}"
OUT="$PACKAGE_DIR/raw_paf/${COMPARISON_ID}.sweepga_${LABEL}.paf.gz"
LOG="$PACKAGE_DIR/logs/${COMPARISON_ID}.sweepga_${LABEL}.${SLURM_JOB_ID:-manual}.command.log"

if [[ ! -s "$QUERY" || ! -s "$TARGET" ]]; then
    echo "missing full whole-genome input FASTA for $COMPARISON_ID; expected symlinks under $PACKAGE_DIR/inputs" >&2
    exit 1
fi

mkdir -p "$PACKAGE_DIR/raw_paf" "$PACKAGE_DIR/logs"

SCRATCH_BASE="${SWEEPGA_DEVSHM_BASE:-/dev/shm}"
if [[ ! -d "$SCRATCH_BASE" || ! -w "$SCRATCH_BASE" ]]; then
    echo "required sweepGA/FastGA scratch base is not writable: $SCRATCH_BASE" >&2
    exit 1
fi

TMPDIR="$(mktemp -d "$SCRATCH_BASE/sg.freq${FASTGA_FREQUENCY}.${SLURM_JOB_ID:-manual}.XXXXXX")"
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
    echo "parameter_set=frequency${FASTGA_FREQUENCY}"
    echo "coordinate_scope=full whole-genome assembly haplotype FASTA records"
    echo "scratch=$TMPDIR"
    echo "scratch_base=$SCRATCH_BASE"
    echo "scratch_policy=sweepGA/FastGA scratch is explicitly under /dev/shm; SLURM_TMPDIR is not used for sweepGA scratch"
    echo "scratch_fastga_sources=$STAGED_QUERY,$STAGED_TARGET"
    echo "sweepga=$SWEEPGA"
    echo "sweepga_which=$SWEEPGA_WHICH"
    echo "sweepga_realpath=$SWEEPGA_REALPATH"
    echo "sweepga_sha256=$SWEEPGA_SHA256"
    "$SWEEPGA" --version || true
    "$SWEEPGA" --help | sed -n '1,140p' || true
    echo "check_fastga_begin"
    "$SWEEPGA" --check-fastga || true
    echo "check_fastga_end"
    TMP_PAF="$TMPDIR/${COMPARISON_ID}.frequency${FASTGA_FREQUENCY}.paf"
    echo "command=$SWEEPGA --fastga --fastga-frequency $FASTGA_FREQUENCY --num-mappings $NUM_MAPPINGS --scaffold-jump $SCAFFOLD_JUMP --temp-dir $TMPDIR --output-file $TMP_PAF $STAGED_QUERY $STAGED_TARGET"
    "$SWEEPGA" --fastga \
        --fastga-frequency "$FASTGA_FREQUENCY" \
        --num-mappings "$NUM_MAPPINGS" \
        --scaffold-jump "$SCAFFOLD_JUMP" \
        --temp-dir "$TMPDIR" \
        --output-file "$TMP_PAF" \
        "$STAGED_QUERY" "$STAGED_TARGET"
    gzip -c "$TMP_PAF" > "$OUT"
    gzip -t "$OUT"
    sha256sum "$OUT" > "$OUT.sha256"
    date -u '+finished_utc=%Y-%m-%dT%H:%M:%SZ'
    echo "output=$OUT"
    echo "bytes=$(wc -c < "$OUT")"
    echo "sha256=$(awk '{print $1}' "$OUT.sha256")"
} 2>&1 | tee "$LOG"
