#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
COMPARISONS="$PACKAGE_DIR/config/comparisons.tsv"
STATUS="$PACKAGE_DIR/summaries/slurm_jobs.tsv"

mkdir -p "$PACKAGE_DIR/logs" "$PACKAGE_DIR/summaries" "$PACKAGE_DIR/raw_paf"
DEVSHM_BASE="${SWEEPGA_DEVSHM_BASE:-/dev/shm}"
printf "stage\tcomparison_id\tjob_id\tnum_mappings\tscaffold_jump\tsweepga_devshm_base\tstdout\n" > "$STATUS"

tail -n +2 "$COMPARISONS" | cut -f1 | while read -r cid; do
    [[ -n "$cid" ]] || continue
    stdout="$PACKAGE_DIR/logs/${cid}.many_many_j0.%j.out"
    submit_out="$(sbatch --parsable \
        --export=ALL,PACKAGE_DIR="$PACKAGE_DIR",SWEEPGA_DEVSHM_BASE="$DEVSHM_BASE" \
        --job-name="wgsg_${cid}" \
        --cpus-per-task="${SWEEPGA_CPUS:-32}" \
        --mem="${SWEEPGA_MEM:-192G}" \
        --time="${SWEEPGA_TIME:-72:00:00}" \
        --output="$stdout" \
        "$SCRIPT_DIR/run_sweepga_one.sh" "$cid" many:many 0)"
    job_id="${submit_out%%;*}"
    printf "raw_many_many\t%s\t%s\tmany:many\t0\t%s\t%s\n" "$cid" "$job_id" "$DEVSHM_BASE" "$stdout" >> "$STATUS"
done

echo "Wrote $STATUS"
