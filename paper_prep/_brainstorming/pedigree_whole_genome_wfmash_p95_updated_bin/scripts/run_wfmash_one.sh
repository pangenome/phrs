#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: run_wfmash_one.sh COMPARISON_ID PARAMETER_SET" >&2
    exit 2
fi

COMPARISON_ID="$1"
PARAMETER_SET="$2"
PACKAGE_DIR="${PACKAGE_DIR:?PACKAGE_DIR must be exported by submit_wfmash_matrix.sh}"
THREADS="${SLURM_CPUS_PER_TASK:-${WFMASH_THREADS:-32}}"
SCRATCH_BASE="${WFMASH_SCRATCH_BASE:-/dev/shm}"
WFMASH_BIN="${WFMASH_BIN:-/home/erikg/bin/wfmash}"
WFMASH_RUN_LABEL="${WFMASH_RUN_LABEL:-updated_bin_v0.24.2-12-ge040aa10}"
WFMASH_SUFFIX="${WFMASH_SUFFIX:-wfmash-v0.24.2-12-ge040aa10}"

QUERY="$PACKAGE_DIR/inputs/${COMPARISON_ID}.query.fa"
TARGET="$PACKAGE_DIR/inputs/${COMPARISON_ID}.target.fa"
PARAMETERS="$PACKAGE_DIR/config/wfmash_parameters.tsv"
OUT_DIR="$PACKAGE_DIR/raw_paf/$WFMASH_RUN_LABEL"
LOG_DIR="$PACKAGE_DIR/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

if [[ ! -s "$QUERY" || ! -s "$TARGET" ]]; then
    echo "Missing whole-genome input FASTA(s): query=$QUERY target=$TARGET" >&2
    exit 3
fi
if [[ ! -x "$WFMASH_BIN" ]]; then
    echo "Missing executable wfmash binary: $WFMASH_BIN" >&2
    exit 5
fi

OPTS="$(awk -F'\t' -v p="$PARAMETER_SET" 'NR > 1 && $1 == p {print $3}' "$PARAMETERS")"
if [[ -z "$OPTS" ]]; then
    echo "Parameter set not found: $PARAMETER_SET" >&2
    exit 4
fi

JOB_ID="${SLURM_JOB_ID:-manual}"
SCRATCH="$SCRATCH_BASE/wfmash.${JOB_ID}.${WFMASH_RUN_LABEL}.${COMPARISON_ID}.${PARAMETER_SET}"
TMP_BASE="$SCRATCH/tmp"
LOCAL_TARGET="$SCRATCH/target.fa"
LOCAL_QUERY="$SCRATCH/query.fa"
LOCAL_PAF="$SCRATCH/${COMPARISON_ID}.${PARAMETER_SET}.paf"
LOCAL_GZ="$LOCAL_PAF.gz"
FINAL_GZ="$OUT_DIR/${COMPARISON_ID}.${PARAMETER_SET}.${WFMASH_SUFFIX}.paf.gz"
CMD_LOG="$LOG_DIR/${WFMASH_RUN_LABEL}.${COMPARISON_ID}.${PARAMETER_SET}.${JOB_ID}.command.log"

cleanup() {
    rm -rf "$SCRATCH"
}
trap cleanup EXIT

mkdir -p "$SCRATCH" "$TMP_BASE"
{
    echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "host=$(hostname)"
    echo "job_id=$JOB_ID"
    echo "run_label=$WFMASH_RUN_LABEL"
    echo "comparison_id=$COMPARISON_ID"
    echo "parameter_set=$PARAMETER_SET"
    echo "threads=$THREADS"
    echo "scratch=$SCRATCH"
    echo "tmp_base=$TMP_BASE"
    echo "target=$TARGET"
    echo "query=$QUERY"
    echo "output=$FINAL_GZ"
    echo "wfmash_bin=$WFMASH_BIN"
    echo "wfmash_which=$(command -v wfmash || true)"
    echo "wfmash_realpath=$(readlink -f "$WFMASH_BIN")"
    echo "wfmash_sha256=$(sha256sum "$WFMASH_BIN" | awk '{print $1}')"
    echo "wfmash_version=$("$WFMASH_BIN" --version 2>&1 | tr '\n' ' ')"
    echo "wfmash_help_first_line=$("$WFMASH_BIN" --help 2>&1 | head -1)"
    echo "wfmash_help_relevant=$("$WFMASH_BIN" --help 2>&1 | grep -E '^[[:space:]]+(-p|-w|-l|-n|-f|-M|-B)|threads' | head -20 | tr '\n' ' ')"
} | tee "$CMD_LOG"

cp "$TARGET" "$LOCAL_TARGET"
cp "$QUERY" "$LOCAL_QUERY"
[[ -s "$TARGET.fai" ]] && cp "$TARGET.fai" "$LOCAL_TARGET.fai" || samtools faidx "$LOCAL_TARGET"
[[ -s "$QUERY.fai" ]] && cp "$QUERY.fai" "$LOCAL_QUERY.fai" || samtools faidx "$LOCAL_QUERY"

echo "command=$WFMASH_BIN $OPTS -t $THREADS -B $TMP_BASE $LOCAL_TARGET $LOCAL_QUERY > $LOCAL_PAF" | tee -a "$CMD_LOG"
# shellcheck disable=SC2086
"$WFMASH_BIN" $OPTS -t "$THREADS" -B "$TMP_BASE" "$LOCAL_TARGET" "$LOCAL_QUERY" > "$LOCAL_PAF"

bgzip -@ "$THREADS" -f "$LOCAL_PAF"
cp "$LOCAL_GZ" "$FINAL_GZ"
sha256sum "$FINAL_GZ" | tee "$FINAL_GZ.sha256" | tee -a "$CMD_LOG"
wc -c "$FINAL_GZ" | tee -a "$CMD_LOG"
echo "status=OK" | tee -a "$CMD_LOG"
