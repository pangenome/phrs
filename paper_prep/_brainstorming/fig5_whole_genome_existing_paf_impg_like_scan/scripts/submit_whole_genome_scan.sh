#!/usr/bin/env bash
#SBATCH --job-name=fig5-impg-like-paf-scan
#SBATCH --partition=workers
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=12:00:00
#SBATCH --output=paper_prep/_brainstorming/fig5_whole_genome_existing_paf_impg_like_scan/logs/slurm-%j.out
#SBATCH --error=paper_prep/_brainstorming/fig5_whole_genome_existing_paf_impg_like_scan/logs/slurm-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(git rev-parse --show-toplevel)}"

SCAN_DIR="paper_prep/_brainstorming/fig5_whole_genome_existing_paf_impg_like_scan"
mkdir -p "$SCAN_DIR/logs" "$SCAN_DIR/manifests" "$SCAN_DIR/summaries"

{
  date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ"
  echo "hostname=$(hostname)"
  echo "SLURM_JOB_ID=${SLURM_JOB_ID:-not_slurm}"
  echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-unset}"
  echo "python=$(command -v python3)"
  echo "pigz=$(command -v pigz || true)"
} | tee "$SCAN_DIR/logs/runtime_environment.${SLURM_JOB_ID:-manual}.txt"

python3 "$SCAN_DIR/scripts/build_paf_input_manifest.py" \
  --bin-size 2000 \
  --out "$SCAN_DIR/manifests/paf_inputs.tsv"

python3 "$SCAN_DIR/scripts/summarize_existing_paf_impg_like_scan.py" \
  --manifest "$SCAN_DIR/manifests/paf_inputs.tsv" \
  --chrom-sizes data/chm13.chrom.sizes \
  --outdir "$SCAN_DIR/summaries" \
  --tmpdir "$SCAN_DIR/summaries/tmp_worker_bin_support" \
  --aggregate-bp 10000,50000

python3 "$SCAN_DIR/scripts/write_report.py"

date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" | tee -a "$SCAN_DIR/logs/runtime_environment.${SLURM_JOB_ID:-manual}.txt"
