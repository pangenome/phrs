#!/usr/bin/env bash
#SBATCH --job-name=chm13-cn-cal-merge
#SBATCH --partition=workers
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_logs/calibration-combine-%j.out
#SBATCH --error=slurm_logs/calibration-combine-%j.err

set -euo pipefail

: "${REPO_ROOT:?submit with --export=ALL,REPO_ROOT=/absolute/repository/path}"
: "${RESULT_DIR:?set the absolute directory containing chunk-0.json ... chunk-3.json}"
CHUNKS="${CHUNKS:-4}"
OUTPUT="${OUTPUT:-${RESULT_DIR}/combined.json}"
ANALYSIS_DIR="${REPO_ROOT}/paper_prep/_brainstorming/chm13_copy_enrichment"
files=()
for (( index=0; index<CHUNKS; index++ )); do
  files+=("${RESULT_DIR}/chunk-${index}.json")
done
cd "${ANALYSIS_DIR}"
echo "job_id=${SLURM_JOB_ID} host=$(hostname) cpus=${SLURM_CPUS_PER_TASK} chunks=${CHUNKS} result_dir=${RESULT_DIR}"
time -p guix shell python python-numpy -- python3 combine_calibration.py \
  "${files[@]}" --output "${OUTPUT}"
echo "combined_calibration=${OUTPUT}"
