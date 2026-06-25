#!/usr/bin/env bash
#SBATCH --job-name=fig5-impg-recover
#SBATCH --partition=workers
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --output=paper_prep/_brainstorming/fig5_whole_genome_existing_paf_impg_like_scan/logs/slurm-recover-%j.out
#SBATCH --error=paper_prep/_brainstorming/fig5_whole_genome_existing_paf_impg_like_scan/logs/slurm-recover-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(git rev-parse --show-toplevel)}"

SCAN_DIR="paper_prep/_brainstorming/fig5_whole_genome_existing_paf_impg_like_scan"
mkdir -p "$SCAN_DIR/logs" "$SCAN_DIR/manifests" "$SCAN_DIR/summaries"

{
  date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ"
  echo "hostname=$(hostname)"
  echo "SLURM_JOB_ID=${SLURM_JOB_ID:-not_slurm}"
  echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-unset}"
  echo "mode=reuse_worker_outputs"
  echo "python=$(command -v python3)"
  echo "pigz=$(command -v pigz || true)"
  echo "sort=$(command -v sort)"
} | tee "$SCAN_DIR/logs/runtime_environment.recover.${SLURM_JOB_ID:-manual}.txt"

python3 "$SCAN_DIR/scripts/summarize_existing_paf_impg_like_scan.py" \
  --manifest "$SCAN_DIR/manifests/paf_inputs.tsv" \
  --chrom-sizes data/chm13.chrom.sizes \
  --outdir "$SCAN_DIR/summaries" \
  --tmpdir "$SCAN_DIR/summaries/tmp_worker_bin_support" \
  --aggregate-bp 10000,50000 \
  --reuse-worker-outputs

python3 "$SCAN_DIR/scripts/write_report.py"

date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" | tee -a "$SCAN_DIR/logs/runtime_environment.recover.${SLURM_JOB_ID:-manual}.txt"
