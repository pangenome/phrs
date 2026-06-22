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
THREADS="${PAFCHOP_THREADS_PER_JOB:-8}"
RUN_LABEL="${MINIMAP2_RUN_LABEL:-v2.31-r1302}"
PARAMETER_SET="${MINIMAP2_PARAMETER_SET:-asm5_allchains}"
MINIMAP2_SUFFIX="${MINIMAP2_SUFFIX:-minimap2-v2.31-r1302}"

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
    echo "threads_per_job=$THREADS"
    cargo build --release --manifest-path "$PAFCHOP_DIR/Cargo.toml"
    mkdir -p "$PACKAGE_DIR/chopped_paf_l${CHOP_LENGTH}_o${OVERLAP}" "$PACKAGE_DIR/summaries/pafchop_l${CHOP_LENGTH}_o${OVERLAP}"
    printf "comparison_id\tinput_paf\toutput_paf\tchop_length\toverlap\tthreads\tstatus\n" > "$PACKAGE_DIR/summaries/chop_manifest_l${CHOP_LENGTH}_o${OVERLAP}.tsv"
    tail -n +2 "$PACKAGE_DIR/config/comparisons.tsv" | cut -f1 | while read -r cid; do
        [[ -n "$cid" ]] || continue
        raw="$PACKAGE_DIR/raw_paf/$RUN_LABEL/${cid}.${PARAMETER_SET}.${MINIMAP2_SUFFIX}.paf.gz"
        out="$PACKAGE_DIR/chopped_paf_l${CHOP_LENGTH}_o${OVERLAP}/${cid}.chopped_l${CHOP_LENGTH}_o${OVERLAP}.paf.gz"
        summary="$PACKAGE_DIR/summaries/pafchop_l${CHOP_LENGTH}_o${OVERLAP}/${cid}.summary.tsv"
        echo "command=$PAFCHOP_DIR/scripts/chop_one.sh $raw $out $summary $cid $CHOP_LENGTH $THREADS"
        "$PAFCHOP_DIR/scripts/chop_one.sh" "$raw" "$out" "$summary" "$cid" "$CHOP_LENGTH" "$THREADS"
        printf "%s\t%s\t%s\t%s\t%s\t%s\tOK\n" "$cid" "$raw" "$out" "$CHOP_LENGTH" "$OVERLAP" "$THREADS" >> "$PACKAGE_DIR/summaries/chop_manifest_l${CHOP_LENGTH}_o${OVERLAP}.tsv"
    done
    cp "$PACKAGE_DIR/summaries/chop_manifest_l${CHOP_LENGTH}_o${OVERLAP}.tsv" "$PACKAGE_DIR/summaries/chop_manifest.tsv"
    date -u '+finished_utc=%Y-%m-%dT%H:%M:%SZ'
} 2>&1 | tee "$PACKAGE_DIR/logs/pafchop_rs_l${CHOP_LENGTH}_o${OVERLAP}.${SLURM_JOB_ID:-manual}.log"
