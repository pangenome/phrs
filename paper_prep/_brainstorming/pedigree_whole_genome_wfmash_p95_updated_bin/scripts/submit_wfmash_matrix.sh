#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
COMPARISONS="$PACKAGE_DIR/config/comparisons.tsv"
PARAMETERS="$PACKAGE_DIR/config/wfmash_parameters.tsv"
STATUS="$PACKAGE_DIR/summaries/wfmash_jobs.tsv"
WFMASH_BIN="${WFMASH_BIN:-/home/erikg/bin/wfmash}"
WFMASH_RUN_LABEL="${WFMASH_RUN_LABEL:-updated_bin_v0.24.2-12-ge040aa10}"
WFMASH_SUFFIX="${WFMASH_SUFFIX:-wfmash-v0.24.2-12-ge040aa10}"

mkdir -p "$PACKAGE_DIR/logs" "$PACKAGE_DIR/summaries" "$PACKAGE_DIR/raw_paf"
SCRATCH_BASE="${WFMASH_SCRATCH_BASE:-/dev/shm}"
CPUS="${WFMASH_CPUS:-32}"
MEM="${WFMASH_MEM:-192G}"
TIME="${WFMASH_TIME:-72:00:00}"

printf "run_label\tcomparison_id\tparameter_set\tjob_id\tstatus\tstdout\tstderr\toutput_paf\toutput_size_bytes\tsha256\twfmash_options\tscratch_base\twfmash_bin\n" > "$STATUS"

PARAMETER_FILTER="${WFMASH_PARAMETER_FILTER:-literal_p95}"

tail -n +2 "$COMPARISONS" | cut -f1 | while read -r cid; do
    [[ -n "$cid" ]] || continue
    tail -n +2 "$PARAMETERS" | cut -f1,3 | while IFS=$'\t' read -r parameter_set opts; do
        [[ -n "$parameter_set" ]] || continue
        [[ "$PARAMETER_FILTER" == "all" || "$PARAMETER_FILTER" == "$parameter_set" ]] || continue
        stdout="$PACKAGE_DIR/logs/${cid}.${parameter_set}.%j.out"
        stderr="$PACKAGE_DIR/logs/${cid}.${parameter_set}.%j.err"
        out_paf="$PACKAGE_DIR/raw_paf/${WFMASH_RUN_LABEL}/${cid}.${parameter_set}.${WFMASH_SUFFIX}.paf.gz"
        submit_out="$(sbatch --parsable \
            --export=ALL,PACKAGE_DIR="$PACKAGE_DIR",WFMASH_SCRATCH_BASE="$SCRATCH_BASE",WFMASH_BIN="$WFMASH_BIN",WFMASH_RUN_LABEL="$WFMASH_RUN_LABEL",WFMASH_SUFFIX="$WFMASH_SUFFIX" \
            --job-name="wfm_${parameter_set}_${cid}" \
            --cpus-per-task="$CPUS" \
            --mem="$MEM" \
            --time="$TIME" \
            --output="$stdout" \
            --error="$stderr" \
            "$SCRIPT_DIR/run_wfmash_one.sh" "$cid" "$parameter_set")"
        job_id="${submit_out%%;*}"
        stdout_actual="${stdout//%j/$job_id}"
        stderr_actual="${stderr//%j/$job_id}"
        printf "%s\t%s\t%s\t%s\tSUBMITTED\t%s\t%s\t%s\t\t\t%s\t%s\t%s\n" \
            "$WFMASH_RUN_LABEL" "$cid" "$parameter_set" "$job_id" "$stdout_actual" "$stderr_actual" "$out_paf" "$opts" "$SCRATCH_BASE" "$WFMASH_BIN" >> "$STATUS"
    done
done

echo "Wrote $STATUS"
