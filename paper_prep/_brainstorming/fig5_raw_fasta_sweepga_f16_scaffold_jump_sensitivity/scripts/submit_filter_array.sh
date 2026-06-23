#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/moosefs/erikg/phrs}"
PACKAGE_DIR="${PACKAGE_DIR:-$REPO_ROOT/paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_scaffold_jump_sensitivity}"

"$PACKAGE_DIR/scripts/build_filter_tasks.py"
task_count="$(tail -n +2 "$PACKAGE_DIR/filter_tasks.tsv" | wc -l)"
if [[ "$task_count" -lt 1 ]]; then
    echo "no filter tasks generated" >&2
    exit 1
fi

job_line="$(sbatch --array="1-${task_count}%8" "$PACKAGE_DIR/scripts/run_filter_array.sbatch")"
job_id="$(awk '{print $NF}' <<< "$job_line")"
{
    printf 'submitted_utc\tjob_id\ttask_count\tarray_limit\tsbatch_line\n'
    printf '%s\t%s\t%s\t%s\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$job_id" "$task_count" "8" "$job_line"
} > "$PACKAGE_DIR/slurm_jobs.tsv"
echo "$job_line"
