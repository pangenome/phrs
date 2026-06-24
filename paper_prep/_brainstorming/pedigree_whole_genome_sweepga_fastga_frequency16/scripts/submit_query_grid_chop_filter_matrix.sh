#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$PACKAGE_DIR/../../.." && pwd)"
PAFCHOP_DIR="${PAFCHOP_DIR:-$REPO_ROOT/paper_prep/_brainstorming/pafchop-rs}"
PAFCHOP_BIN="${PAFCHOP_BIN:-$REPO_ROOT/target/release/pafchop}"
RAW_PAF_DIR="${RAW_PAF_DIR:-$PACKAGE_DIR/raw_paf}"
LENGTHS="${QUERY_GRID_CHOP_LENGTHS:-10000 5000 2000 1000}"
CPUS="${QUERY_GRID_CPUS_PER_TASK:-8}"
MEM="${QUERY_GRID_MEM:-48G}"
TIME="${QUERY_GRID_TIME:-24:00:00}"
JOBS_PARALLEL="${QUERY_GRID_ARRAY_PARALLEL:-6}"
TASKS_TSV="$PACKAGE_DIR/summaries/query_grid_chop_filter_tasks.tsv"
SUBMIT_TSV="$PACKAGE_DIR/summaries/query_grid_chop_filter_submission.tsv"

mkdir -p "$PACKAGE_DIR/summaries" "$PACKAGE_DIR/logs/query_grid_chop_filter"

cargo build --release --manifest-path "$PAFCHOP_DIR/Cargo.toml"

printf "comparison_id\tchop_length_bp\n" > "$TASKS_TSV"
while IFS=$'\t' read -r cid _rest; do
    [[ "$cid" == "comparison_id" || -z "$cid" ]] && continue
    raw="$RAW_PAF_DIR/${cid}.sweepga_frequency16_many_many_j0.paf.gz"
    if [[ ! -s "$raw" ]]; then
        echo "missing raw PAF for $cid: $raw" >&2
        exit 1
    fi
    for length in $LENGTHS; do
        printf "%s\t%s\n" "$cid" "$length" >> "$TASKS_TSV"
    done
done < "$PACKAGE_DIR/config/comparisons.tsv"

task_count="$(( $(wc -l < "$TASKS_TSV") - 1 ))"
if [[ "$task_count" -lt 1 ]]; then
    echo "no query-grid chop/filter tasks generated" >&2
    exit 1
fi

job_output="$PACKAGE_DIR/logs/query_grid_chop_filter/slurm-%A_%a.out"
sbatch_output="$(
    sbatch \
        --job-name=fig5_qgrid_f16 \
        --cpus-per-task="$CPUS" \
        --mem="$MEM" \
        --time="$TIME" \
        --array="1-${task_count}%${JOBS_PARALLEL}" \
        --output="$job_output" \
        --export="ALL,PACKAGE_DIR=$PACKAGE_DIR,RAW_PAF_DIR=$RAW_PAF_DIR,PAFCHOP_BIN=$PAFCHOP_BIN,QUERY_GRID_TASKS_TSV=$TASKS_TSV,SWEEPGA_DEVSHM_BASE=/dev/shm" \
        "$SCRIPT_DIR/run_query_grid_chop_filter_one.sh"
)"
job_id="$(awk '{print $NF}' <<< "$sbatch_output")"
{
    printf "submitted_utc\tjob_id\tarray_spec\ttask_count\tlengths\tcpus_per_task\tmem\ttime\tarray_parallel\traw_paf_dir\tpafchop_bin\tpafchop_sha256\tsweepga\tpigz\ttasks_tsv\tcommand\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$job_id" "1-${task_count}%${JOBS_PARALLEL}" "$task_count" "$LENGTHS" "$CPUS" "$MEM" "$TIME" "$JOBS_PARALLEL" \
        "$RAW_PAF_DIR" "$PAFCHOP_BIN" "$(sha256sum "$PAFCHOP_BIN" | awk '{print $1}')" "${SWEEPGA:-/home/erikg/.cargo/bin/sweepga}" "$(command -v pigz)" "$TASKS_TSV" \
        "sbatch --array=1-${task_count}%${JOBS_PARALLEL} --cpus-per-task=$CPUS --mem=$MEM --time=$TIME $SCRIPT_DIR/run_query_grid_chop_filter_one.sh"
} > "$SUBMIT_TSV"

echo "$sbatch_output"
echo "Wrote $TASKS_TSV"
echo "Wrote $SUBMIT_TSV"
