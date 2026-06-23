#!/usr/bin/env bash
#SBATCH --job-name=fig5_f16_filter_modes
#SBATCH --partition=tux
#SBATCH --nodelist=tux05
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=08:00:00

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/moosefs/erikg/phrs}"
PANEL_DIR="${PANEL_DIR:-$REPO_ROOT/paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_chop_filter_sensitivity_panels}"
SOURCE_PACKAGE="${SOURCE_PACKAGE:-$REPO_ROOT/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16}"
SWEEPGA="${SWEEPGA:-/home/erikg/.cargo/bin/sweepga}"
SCRATCH_BASE="${SWEEPGA_DEVSHM_BASE:-/dev/shm}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
JOB_ID="${SLURM_JOB_ID:-manual}"

mkdir -p "$PANEL_DIR/logs" "$PANEL_DIR/filtered_paf"

if [[ ! -x "$SWEEPGA" ]]; then
    echo "missing sweepga: $SWEEPGA" >&2
    exit 1
fi
if [[ ! -d "$SCRATCH_BASE" || ! -w "$SCRATCH_BASE" ]]; then
    echo "scratch base is not writable: $SCRATCH_BASE" >&2
    exit 1
fi

task_file="$PANEL_DIR/filter_tasks.tsv"
if [[ ! -s "$task_file" ]]; then
    echo "missing task file: $task_file" >&2
    exit 1
fi

line="$(awk -F '\t' -v task="$TASK_ID" 'NR == task + 1 { print }' "$task_file")"
if [[ -z "$line" ]]; then
    echo "no task line for array id $TASK_ID" >&2
    exit 1
fi

IFS=$'\t' read -r comparison_id chop_length_bp filter_mode filter_label scaffold_jump scoring num_mappings source_paf output_paf <<< "$line"
log="$PANEL_DIR/logs/${comparison_id}.l${chop_length_bp}.${filter_mode}.${JOB_ID}_${TASK_ID}.log"

{
    date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
    echo "job_id=$JOB_ID"
    echo "array_task_id=$TASK_ID"
    echo "host=$(hostname)"
    echo "comparison_id=$comparison_id"
    echo "chop_length_bp=$chop_length_bp"
    echo "filter_mode=$filter_mode"
    echo "filter_label=$filter_label"
    echo "scaffold_jump=$scaffold_jump"
    echo "scoring=$scoring"
    echo "num_mappings=$num_mappings"
    echo "source_paf=$source_paf"
    echo "output_paf=$output_paf"
    echo "sweepga=$SWEEPGA"
    echo "sweepga_sha256=$(sha256sum "$SWEEPGA" | awk '{print $1}')"
    echo "threads=$THREADS"
    echo "scratch_base=$SCRATCH_BASE"
} > "$log"

if [[ ! -s "$source_paf" ]]; then
    echo "missing source PAF: $source_paf" | tee -a "$log" >&2
    exit 1
fi

if [[ -s "$output_paf" ]] && pigz -t -p "$THREADS" "$output_paf"; then
    echo "output exists and validates; skipping" >> "$log"
    exit 0
fi

work="$(mktemp -d "$SCRATCH_BASE/fig5_filter_modes.${comparison_id}.l${chop_length_bp}.${filter_mode}.${JOB_ID}_${TASK_ID}.XXXXXX")"
trap 'rm -rf "$work"' EXIT

pigz -t -p "$THREADS" "$source_paf"
pigz -dc -p "$THREADS" "$source_paf" > "$work/input.paf"

"$SWEEPGA" \
    --num-mappings "$num_mappings" \
    --scaffold-jump "$scaffold_jump" \
    --scoring "$scoring" \
    --threads "$THREADS" \
    --output-file "$work/filtered.paf" \
    "$work/input.paf" >> "$log" 2>&1

tmp_out="${output_paf}.tmp.${JOB_ID}_${TASK_ID}"
rm -f "$tmp_out"
pigz -p "$THREADS" -c "$work/filtered.paf" > "$tmp_out"
pigz -t -p "$THREADS" "$tmp_out"
mv "$tmp_out" "$output_paf"
sha256sum "$output_paf" > "${output_paf}.sha256"

{
    date -u '+finished_utc=%Y-%m-%dT%H:%M:%SZ'
    sha256sum "$source_paf" "$output_paf"
} >> "$log"
