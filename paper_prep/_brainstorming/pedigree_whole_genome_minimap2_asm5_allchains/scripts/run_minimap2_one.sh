#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: run_minimap2_one.sh COMPARISON_ID PARAMETER_SET" >&2
    exit 2
fi

COMPARISON_ID="$1"
PARAMETER_SET="$2"
PACKAGE_DIR="${PACKAGE_DIR:?PACKAGE_DIR must be exported by submit_minimap2_matrix.sh}"
THREADS="${SLURM_CPUS_PER_TASK:-${MINIMAP2_THREADS:-32}}"
SCRATCH_BASE="${MINIMAP2_SCRATCH_BASE:-/dev/shm}"
MINIMAP2_BIN="${MINIMAP2_BIN:-/home/erikg/bin/minimap2}"
MINIMAP2_RUN_LABEL="${MINIMAP2_RUN_LABEL:-v2.31-r1302}"
MINIMAP2_SUFFIX="${MINIMAP2_SUFFIX:-minimap2-v2.31-r1302}"

QUERY="$PACKAGE_DIR/inputs/${COMPARISON_ID}.query.fa"
TARGET="$PACKAGE_DIR/inputs/${COMPARISON_ID}.target.fa"
PARAMETERS="$PACKAGE_DIR/config/minimap2_parameters.tsv"
OUT_DIR="$PACKAGE_DIR/raw_paf/$MINIMAP2_RUN_LABEL"
LOG_DIR="$PACKAGE_DIR/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

if [[ ! -s "$QUERY" || ! -s "$TARGET" ]]; then
    echo "Missing whole-genome input FASTA(s): query=$QUERY target=$TARGET" >&2
    exit 3
fi
if [[ ! -x "$MINIMAP2_BIN" ]]; then
    echo "Missing executable minimap2 binary: $MINIMAP2_BIN" >&2
    exit 5
fi

OPTS="$(awk -F'\t' -v p="$PARAMETER_SET" 'NR > 1 && $1 == p {print $3}' "$PARAMETERS")"
if [[ -z "$OPTS" ]]; then
    echo "Parameter set not found: $PARAMETER_SET" >&2
    exit 4
fi

JOB_ID="${SLURM_JOB_ID:-manual}"
SCRATCH="$SCRATCH_BASE/minimap2.${JOB_ID}.${MINIMAP2_RUN_LABEL}.${COMPARISON_ID}.${PARAMETER_SET}"
LOCAL_TARGET="$SCRATCH/target.fa"
LOCAL_QUERY="$SCRATCH/query.fa"
LOCAL_GZ="$SCRATCH/${COMPARISON_ID}.${PARAMETER_SET}.paf.gz"
FINAL_GZ="$OUT_DIR/${COMPARISON_ID}.${PARAMETER_SET}.${MINIMAP2_SUFFIX}.paf.gz"
CMD_LOG="$LOG_DIR/${MINIMAP2_RUN_LABEL}.${COMPARISON_ID}.${PARAMETER_SET}.${JOB_ID}.command.log"

cleanup() {
    rm -rf "$SCRATCH"
}
trap cleanup EXIT

mkdir -p "$SCRATCH"
{
    echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "host=$(hostname)"
    echo "job_id=$JOB_ID"
    echo "run_label=$MINIMAP2_RUN_LABEL"
    echo "comparison_id=$COMPARISON_ID"
    echo "parameter_set=$PARAMETER_SET"
    echo "threads=$THREADS"
    echo "scratch=$SCRATCH"
    echo "target=$TARGET"
    echo "query=$QUERY"
    echo "scope=full whole-genome child-haplotype query FASTA and full whole-genome joint-parent target FASTA"
    echo "output=$FINAL_GZ"
    echo "minimap2_bin=$MINIMAP2_BIN"
    echo "minimap2_which=$(command -v minimap2 || true)"
    echo "minimap2_realpath=$(readlink -f "$MINIMAP2_BIN")"
    echo "minimap2_sha256=$(sha256sum "$MINIMAP2_BIN" | awk '{print $1}')"
    echo "minimap2_version=$("$MINIMAP2_BIN" --version 2>&1 | tr '\n' ' ')"
    echo "pigz_which=$(command -v pigz || true)"
    echo "pigz_version=$(pigz --version 2>&1 | tr '\n' ' ')"
} | tee "$CMD_LOG"

cp "$TARGET" "$LOCAL_TARGET"
cp "$QUERY" "$LOCAL_QUERY"
[[ -s "$TARGET.fai" ]] && cp "$TARGET.fai" "$LOCAL_TARGET.fai" || samtools faidx "$LOCAL_TARGET"
[[ -s "$QUERY.fai" ]] && cp "$QUERY.fai" "$LOCAL_QUERY.fai" || samtools faidx "$LOCAL_QUERY"

echo "command=$MINIMAP2_BIN $OPTS -t $THREADS $LOCAL_TARGET $LOCAL_QUERY | pigz -p $THREADS > $LOCAL_GZ" | tee -a "$CMD_LOG"
# shellcheck disable=SC2086
"$MINIMAP2_BIN" $OPTS -t "$THREADS" "$LOCAL_TARGET" "$LOCAL_QUERY" | pigz -p "$THREADS" > "$LOCAL_GZ"

cp "$LOCAL_GZ" "$FINAL_GZ"
sha256sum "$FINAL_GZ" | tee "$FINAL_GZ.sha256" | tee -a "$CMD_LOG"
wc -c "$FINAL_GZ" | tee -a "$CMD_LOG"
echo "status=OK" | tee -a "$CMD_LOG"
