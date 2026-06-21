#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
COMPARISONS="$PACKAGE_DIR/config/comparisons.tsv"
PARAMETERS="$PACKAGE_DIR/config/minimap2_parameters.tsv"
STATUS="$PACKAGE_DIR/summaries/slurm_jobs.tsv"
MINIMAP2_BIN="${MINIMAP2_BIN:-/home/erikg/bin/minimap2}"
MINIMAP2_RUN_LABEL="${MINIMAP2_RUN_LABEL:-v2.31-r1302}"
MINIMAP2_SUFFIX="${MINIMAP2_SUFFIX:-minimap2-v2.31-r1302}"

mkdir -p "$PACKAGE_DIR/logs" "$PACKAGE_DIR/summaries" "$PACKAGE_DIR/raw_paf"
SCRATCH_BASE="${MINIMAP2_SCRATCH_BASE:-/dev/shm}"
CPUS="${MINIMAP2_CPUS:-32}"
MEM="${MINIMAP2_MEM:-192G}"
TIME="${MINIMAP2_TIME:-72:00:00}"

printf "run_label\tcomparison_id\tparameter_set\tjob_id\tstatus\tstdout\tstderr\toutput_paf\toutput_size_bytes\tsha256\tminimap2_options\tscratch_base\tminimap2_bin\n" > "$STATUS"

PARAMETER_FILTER="${MINIMAP2_PARAMETER_FILTER:-asm5_allchains}"

tail -n +2 "$COMPARISONS" | cut -f1 | while read -r cid; do
    [[ -n "$cid" ]] || continue
    tail -n +2 "$PARAMETERS" | cut -f1,3 | while IFS=$'\t' read -r parameter_set opts; do
        [[ -n "$parameter_set" ]] || continue
        [[ "$PARAMETER_FILTER" == "all" || "$PARAMETER_FILTER" == "$parameter_set" ]] || continue
        stdout="$PACKAGE_DIR/logs/${cid}.${parameter_set}.%j.out"
        stderr="$PACKAGE_DIR/logs/${cid}.${parameter_set}.%j.err"
        out_paf="$PACKAGE_DIR/raw_paf/${MINIMAP2_RUN_LABEL}/${cid}.${parameter_set}.${MINIMAP2_SUFFIX}.paf.gz"
        submit_out="$(sbatch --parsable \
            --export=ALL,PACKAGE_DIR="$PACKAGE_DIR",MINIMAP2_SCRATCH_BASE="$SCRATCH_BASE",MINIMAP2_BIN="$MINIMAP2_BIN",MINIMAP2_RUN_LABEL="$MINIMAP2_RUN_LABEL",MINIMAP2_SUFFIX="$MINIMAP2_SUFFIX" \
            --job-name="mm2_${parameter_set}_${cid}" \
            --cpus-per-task="$CPUS" \
            --mem="$MEM" \
            --time="$TIME" \
            --output="$stdout" \
            --error="$stderr" \
            "$SCRIPT_DIR/run_minimap2_one.sh" "$cid" "$parameter_set")"
        job_id="${submit_out%%;*}"
        stdout_actual="${stdout//%j/$job_id}"
        stderr_actual="${stderr//%j/$job_id}"
        printf "%s\t%s\t%s\t%s\tSUBMITTED\t%s\t%s\t%s\t\t\t%s\t%s\t%s\n" \
            "$MINIMAP2_RUN_LABEL" "$cid" "$parameter_set" "$job_id" "$stdout_actual" "$stderr_actual" "$out_paf" "$opts" "$SCRATCH_BASE" "$MINIMAP2_BIN" >> "$STATUS"
    done
done

echo "Wrote $STATUS"
