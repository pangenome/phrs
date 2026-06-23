#!/usr/bin/env bash
#SBATCH --job-name=fig5_f16_sj
#SBATCH --partition=tux
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=10:00:00
#SBATCH --array=1-128%8

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/moosefs/erikg/phrs}"
PACKAGE_DIR="${PACKAGE_DIR:-$REPO_ROOT/paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_scaffold_jump_sensitivity}"
SWEEPGA="${SWEEPGA:-/home/erikg/.cargo/bin/sweepga}"
SCRATCH_BASE="${SWEEPGA_DEVSHM_BASE:-/dev/shm}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
JOB_ID="${SLURM_JOB_ID:-manual}"

mkdir -p "$PACKAGE_DIR/logs" "$PACKAGE_DIR/filtered_paf"

if [[ ! -x "$SWEEPGA" ]]; then
    echo "missing sweepga: $SWEEPGA" >&2
    exit 1
fi
if [[ ! -d "$SCRATCH_BASE" || ! -w "$SCRATCH_BASE" ]]; then
    echo "scratch base is not writable: $SCRATCH_BASE" >&2
    exit 1
fi

task_file="$PACKAGE_DIR/filter_tasks.tsv"
if [[ ! -s "$task_file" ]]; then
    "$PACKAGE_DIR/scripts/build_filter_tasks.py"
fi

line="$(awk -F '\t' -v task="$TASK_ID" 'NR == task + 1 { print }' "$task_file")"
if [[ -z "$line" ]]; then
    echo "no task line for array id $TASK_ID" >&2
    exit 1
fi

IFS=$'\t' read -r task_index comparison_id source_paf output_paf cell_id scaffold_jump min_aln_length scoring num_mappings scaffold_mass overlap command_template <<< "$line"
log="$PACKAGE_DIR/logs/${comparison_id}.${cell_id}.${JOB_ID}_${TASK_ID}.log"

{
    date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
    echo "job_id=$JOB_ID"
    echo "array_task_id=$TASK_ID"
    echo "host=$(hostname)"
    echo "comparison_id=$comparison_id"
    echo "cell_id=$cell_id"
    echo "scaffold_jump=$scaffold_jump"
    echo "min_aln_length=$min_aln_length"
    echo "scoring=$scoring"
    echo "num_mappings=$num_mappings"
    echo "scaffold_mass=$scaffold_mass"
    echo "overlap=$overlap"
    echo "source_paf=$source_paf"
    echo "output_paf=$output_paf"
    echo "sweepga=$SWEEPGA"
    echo "sweepga_version=$("$SWEEPGA" --version)"
    echo "sweepga_sha256=$(sha256sum "$SWEEPGA" | awk '{print $1}')"
    echo "threads=$THREADS"
    echo "scratch_base=$SCRATCH_BASE"
    echo "command_template=$command_template"
} > "$log"

if [[ ! -s "$source_paf" ]]; then
    echo "missing source PAF: $source_paf" | tee -a "$log" >&2
    exit 1
fi

if [[ -s "$output_paf" ]] && pigz -t -p "$THREADS" "$output_paf"; then
    echo "output exists and validates; skipping" >> "$log"
    exit 0
fi

work="$(mktemp -d "$SCRATCH_BASE/fig5_f16_sj.${comparison_id}.${cell_id}.${JOB_ID}_${TASK_ID}.XXXXXX")"
trap 'rm -rf "$work"' EXIT
export THREADS
export WORK="$work"

pigz -t -p "$THREADS" "$source_paf"
pigz -dc -p "$THREADS" "$source_paf" > "$work/input.paf"

cmd=("$SWEEPGA" --num-mappings "$num_mappings" --scaffold-jump "$scaffold_jump" --scaffold-mass "$scaffold_mass" --scoring "$scoring")
if [[ "$min_aln_length" != "default" ]]; then
    cmd+=(--min-aln-length "$min_aln_length")
fi
if [[ "$overlap" != "default" ]]; then
    cmd+=(--overlap "$overlap")
fi
cmd+=(--threads "$THREADS" --output-file "$work/filtered.paf" "$work/input.paf")

printf 'command_executed=' >> "$log"
printf '%q ' "${cmd[@]}" >> "$log"
printf '\n' >> "$log"
"${cmd[@]}" >> "$log" 2>&1

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
