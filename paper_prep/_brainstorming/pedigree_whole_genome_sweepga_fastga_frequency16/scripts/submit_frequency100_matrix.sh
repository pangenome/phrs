#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
COMPARISONS="$PACKAGE_DIR/config/comparisons.tsv"
mkdir -p "$PACKAGE_DIR/logs" "$PACKAGE_DIR/summaries"

JOBS="$PACKAGE_DIR/summaries/slurm_jobs.tsv"
printf 'stage\tcomparison_id\tjob_id\tparameter_set\tfastga_frequency\tnum_mappings\tscaffold_jump\tsweepga_devshm_base\tstdout\n' > "$JOBS"

tail -n +2 "$COMPARISONS" | while IFS=$'\t' read -r comparison_id _rest; do
    [[ -n "$comparison_id" ]] || continue
    stdout="$PACKAGE_DIR/logs/${comparison_id}.frequency16_many_many_j0.%j.out"
    job_id="$(
        sbatch --parsable \
            --job-name="sg-f16-${comparison_id}" \
            --cpus-per-task="${SLURM_CPUS_PER_TASK_OVERRIDE:-32}" \
            --mem="${SLURM_MEM_OVERRIDE:-180G}" \
            --time="${SLURM_TIME_OVERRIDE:-24:00:00}" \
            --output="$stdout" \
            --export=ALL,PACKAGE_DIR="$PACKAGE_DIR",SWEEPGA=/home/erikg/.cargo/bin/sweepga,SWEEPGA_DEVSHM_BASE=/dev/shm \
            "$SCRIPT_DIR/run_sweepga_frequency_one.sh" "$comparison_id" 16 many:many 0
    )"
    printf 'raw_frequency16\t%s\t%s\tfrequency16\t16\tmany:many\t0\t/dev/shm\t%s\n' "$comparison_id" "$job_id" "$stdout" >> "$JOBS"
done

cat "$JOBS"
