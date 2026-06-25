#!/usr/bin/env bash
set -euo pipefail

if [[ -n "${PACKAGE_DIR:-}" ]]; then
    SCRIPT_DIR="$PACKAGE_DIR/scripts"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
fi
REPO_ROOT="$(cd "$PACKAGE_DIR/../../.." && pwd)"

TASKS_TSV="${QUERY_GRID_TASKS_TSV:-$PACKAGE_DIR/summaries/query_grid_chop_filter_tasks.tsv}"
TASK_INDEX="${QUERY_GRID_TASK_INDEX:-${SLURM_ARRAY_TASK_ID:-}}"
RAW_PAF_DIR="${RAW_PAF_DIR:-$PACKAGE_DIR/raw_paf}"
PAFCHOP_BIN="${PAFCHOP_BIN:-$REPO_ROOT/paper_prep/_brainstorming/pafchop-rs/target/release/pafchop}"
SWEEPGA="${SWEEPGA:-/home/erikg/.cargo/bin/sweepga}"
PIGZ="${PIGZ:-$(command -v pigz)}"
THREADS="${SLURM_CPUS_PER_TASK:-${QUERY_GRID_THREADS:-8}}"
PIGZ_COMPRESSION_LEVEL="${PIGZ_COMPRESSION_LEVEL:--6}"
SCRATCH_BASE="${SWEEPGA_DEVSHM_BASE:-/dev/shm}"
CHOP_MODE="query-grid"
OVERLAP="0"
FILTER_ID="one_to_one_ani_o0"
NUM_MAPPINGS="1:1"
SCAFFOLD_JUMP="0"
SCORING="ani"
FILTER_OVERLAP="0"

if [[ -z "$TASK_INDEX" ]]; then
    echo "QUERY_GRID_TASK_INDEX or SLURM_ARRAY_TASK_ID is required" >&2
    exit 2
fi
if [[ "$TASK_INDEX" -lt 1 ]]; then
    echo "task index is 1-based; got $TASK_INDEX" >&2
    exit 2
fi
if [[ ! -r "$TASKS_TSV" ]]; then
    echo "missing task matrix: $TASKS_TSV" >&2
    exit 2
fi
if [[ ! -x "$PAFCHOP_BIN" ]]; then
    echo "missing executable pafchop binary: $PAFCHOP_BIN" >&2
    exit 2
fi
if [[ ! -x "$SWEEPGA" ]]; then
    echo "missing executable sweepGA binary: $SWEEPGA" >&2
    exit 2
fi
if [[ ! -x "$PIGZ" ]]; then
    echo "missing executable pigz: $PIGZ" >&2
    exit 2
fi
if [[ ! -d "$SCRATCH_BASE" || ! -w "$SCRATCH_BASE" ]]; then
    echo "required /dev/shm scratch base is not writable: $SCRATCH_BASE" >&2
    exit 2
fi

row="$(awk -v n="$TASK_INDEX" 'BEGIN{FS="\t"} NR==n+1 {print $0}' "$TASKS_TSV")"
if [[ -z "$row" ]]; then
    echo "no row $TASK_INDEX in $TASKS_TSV" >&2
    exit 2
fi
IFS=$'\t' read -r comparison_id chop_length <<< "$row"

mkdir -p \
    "$PACKAGE_DIR/chopped_paf_qgrid_l${chop_length}_o0" \
    "$PACKAGE_DIR/filtered_paf_chop_sensitivity_query_grid/l${chop_length}" \
    "$PACKAGE_DIR/summaries/pafchop_qgrid_l${chop_length}_o0" \
    "$PACKAGE_DIR/summaries/query_grid_chop_filter_status" \
    "$PACKAGE_DIR/logs/query_grid_chop_filter"

raw="$RAW_PAF_DIR/${comparison_id}.sweepga_frequency16_many_many_j0.paf.gz"
chopped="$PACKAGE_DIR/chopped_paf_qgrid_l${chop_length}_o0/${comparison_id}.chopped_l${chop_length}_o0_query_grid.paf.gz"
chop_summary="$PACKAGE_DIR/summaries/pafchop_qgrid_l${chop_length}_o0/${comparison_id}.summary.tsv"
filtered="$PACKAGE_DIR/filtered_paf_chop_sensitivity_query_grid/l${chop_length}/${comparison_id}.${FILTER_ID}.chopped_l${chop_length}_o0_query_grid.paf.gz"
status="$PACKAGE_DIR/summaries/query_grid_chop_filter_status/${comparison_id}.l${chop_length}.tsv"
log="$PACKAGE_DIR/logs/query_grid_chop_filter/${comparison_id}.l${chop_length}.${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-$TASK_INDEX}.log"

tmp_chopped="${chopped}.tmp.${SLURM_JOB_ID:-manual}.${TASK_INDEX}"
tmp_summary="${chop_summary}.tmp.${SLURM_JOB_ID:-manual}.${TASK_INDEX}"
tmp_filtered="${filtered}.tmp.${SLURM_JOB_ID:-manual}.${TASK_INDEX}"
scratch_dir="$(mktemp -d "$SCRATCH_BASE/sweepga_qgrid.${SLURM_JOB_ID:-manual}.${TASK_INDEX}.XXXXXX")"
started_utc="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
host="$(hostname)"
pafchop_sha256="$(sha256sum "$PAFCHOP_BIN" | awk '{print $1}')"
sweepga_sha256="$(sha256sum "$SWEEPGA" | awk '{print $1}')"
pigz_sha256="$(sha256sum "$PIGZ" | awk '{print $1}')"
chop_command="$PIGZ -dc $raw | $PAFCHOP_BIN --length $chop_length --overlap 0 --chunk-mode query-grid --threads $THREADS --comparison-id $comparison_id --summary $tmp_summary | $PIGZ $PIGZ_COMPRESSION_LEVEL -p $THREADS > $tmp_chopped"
filter_command="$PIGZ -dc $chopped > $scratch_dir/input.paf && $SWEEPGA --threads $THREADS --num-mappings 1:1 --scaffold-jump 0 --scoring ani --overlap 0 --output-file $scratch_dir/filtered.paf $scratch_dir/input.paf && $PIGZ $PIGZ_COMPRESSION_LEVEL -p $THREADS -c $scratch_dir/filtered.paf > $tmp_filtered"

finish_status="FAILED"
cleanup() {
    rc=$?
    finished_utc="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    rm -rf "$scratch_dir"
    if [[ "$rc" -ne 0 ]]; then
        rm -f "$tmp_chopped" "$tmp_summary" "$tmp_filtered"
    fi
    printf "comparison_id\tchop_length_bp\tfilter_id\tjob_id\tarray_task_id\thost\tstarted_utc\tfinished_utc\tstatus\traw_paf\tchopped_paf\tfiltered_paf\tchop_summary\tchunk_mode\toverlap_bp\tnum_mappings\tscaffold_jump\tscoring\tfilter_overlap\tthreads\tpigz_compression_level\tscratch_dir\tpafchop_bin\tpafchop_sha256\tsweepga_bin\tsweepga_sha256\tpigz_bin\tpigz_sha256\tchop_command\tfilter_command\n" > "$status"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$comparison_id" "$chop_length" "$FILTER_ID" "${SLURM_JOB_ID:-manual}" "${SLURM_ARRAY_TASK_ID:-$TASK_INDEX}" "$host" "$started_utc" "$finished_utc" "$finish_status" \
        "$raw" "$chopped" "$filtered" "$chop_summary" "$CHOP_MODE" "$OVERLAP" "$NUM_MAPPINGS" "$SCAFFOLD_JUMP" "$SCORING" "$FILTER_OVERLAP" "$THREADS" "$PIGZ_COMPRESSION_LEVEL" "$scratch_dir" \
        "$PAFCHOP_BIN" "$pafchop_sha256" "$SWEEPGA" "$sweepga_sha256" "$PIGZ" "$pigz_sha256" "$chop_command" "$filter_command" >> "$status"
    exit "$rc"
}
trap cleanup EXIT

exec > >(tee "$log") 2>&1
{
    echo "started_utc=$started_utc"
    echo "host=$host"
    echo "job_id=${SLURM_JOB_ID:-manual}"
    echo "array_task_id=${SLURM_ARRAY_TASK_ID:-$TASK_INDEX}"
    echo "comparison_id=$comparison_id"
    echo "chop_length=$chop_length"
    echo "chunk_mode=$CHOP_MODE"
    echo "threads=$THREADS"
    echo "pigz_compression_level=$PIGZ_COMPRESSION_LEVEL"
    echo "scratch_dir=$scratch_dir"
    echo "raw_paf=$raw"
    echo "chopped_paf=$chopped"
    echo "filtered_paf=$filtered"
    echo "pafchop_bin=$PAFCHOP_BIN"
    echo "pafchop_sha256=$pafchop_sha256"
    echo "sweepga_bin=$SWEEPGA"
    echo "sweepga_sha256=$sweepga_sha256"
    echo "pigz_bin=$PIGZ"
    echo "pigz_sha256=$pigz_sha256"
    echo "chop_command=$chop_command"
    [[ -s "$raw" ]]
    "$PIGZ" -dc "$raw" \
        | "$PAFCHOP_BIN" --length "$chop_length" --overlap 0 --chunk-mode "$CHOP_MODE" --threads "$THREADS" --comparison-id "$comparison_id" --summary "$tmp_summary" \
        | "$PIGZ" "$PIGZ_COMPRESSION_LEVEL" -p "$THREADS" > "$tmp_chopped"
    mv "$tmp_chopped" "$chopped"
    mv "$tmp_summary" "$chop_summary"
    "$PIGZ" -t "$chopped"
    sha256sum "$chopped" > "${chopped}.sha256"

    echo "filter_command=$filter_command"
    "$PIGZ" -dc "$chopped" > "$scratch_dir/input.paf"
    "$SWEEPGA" \
        --threads "$THREADS" \
        --num-mappings "$NUM_MAPPINGS" \
        --scaffold-jump "$SCAFFOLD_JUMP" \
        --scoring "$SCORING" \
        --overlap "$FILTER_OVERLAP" \
        --output-file "$scratch_dir/filtered.paf" \
        "$scratch_dir/input.paf"
    "$PIGZ" "$PIGZ_COMPRESSION_LEVEL" -p "$THREADS" -c "$scratch_dir/filtered.paf" > "$tmp_filtered"
    mv "$tmp_filtered" "$filtered"
    "$PIGZ" -t "$filtered"
    sha256sum "$filtered" > "${filtered}.sha256"
    finish_status="OK"
    date -u '+finished_utc=%Y-%m-%dT%H:%M:%SZ'
}
