#!/usr/bin/env bash
#SBATCH --job-name=fig5-impg-compact
#SBATCH --partition=workers
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=paper_prep/_brainstorming/fig5_whole_genome_existing_paf_impg_like_scan/logs/slurm-compact-%j.out
#SBATCH --error=paper_prep/_brainstorming/fig5_whole_genome_existing_paf_impg_like_scan/logs/slurm-compact-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(git rev-parse --show-toplevel)}"

SCAN_DIR="paper_prep/_brainstorming/fig5_whole_genome_existing_paf_impg_like_scan"
mkdir -p "$SCAN_DIR/logs" "$SCAN_DIR/summaries"

{
  date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ"
  echo "hostname=$(hostname)"
  echo "SLURM_JOB_ID=${SLURM_JOB_ID:-not_slurm}"
  echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-unset}"
  echo "mode=compact_worker_output_recovery"
  echo "awk=$(command -v awk)"
  echo "sort=$(command -v sort)"
} | tee "$SCAN_DIR/logs/runtime_environment.compact.${SLURM_JOB_ID:-manual}.txt"

bash "$SCAN_DIR/scripts/write_compact_summaries_from_worker_outputs.sh"
python3 "$SCAN_DIR/scripts/write_report.py"

date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" | tee -a "$SCAN_DIR/logs/runtime_environment.compact.${SLURM_JOB_ID:-manual}.txt"
