#!/usr/bin/env bash
#SBATCH --job-name=chm13-cn-cal
#SBATCH --partition=workers
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_logs/calibration-%A_%a.out
#SBATCH --error=slurm_logs/calibration-%A_%a.err

set -euo pipefail

: "${REPO_ROOT:?submit with --export=ALL,REPO_ROOT=/absolute/repository/path}"
CHUNKS="${CHUNKS:-4}"
PSEUDO_OBSERVATIONS="${PSEUDO_OBSERVATIONS:-250}"
PERMUTATIONS="${PERMUTATIONS:-9999}"
EXACT_DRAWS="${EXACT_DRAWS:-100000}"
MASTER_SEED="${MASTER_SEED:-2026071391}"
RESULT_DIR="${RESULT_DIR:-${SLURM_SUBMIT_DIR}/calibration_results}"
ANALYSIS_DIR="${REPO_ROOT}/paper_prep/_brainstorming/chm13_copy_enrichment"

mkdir -p "${RESULT_DIR}"
cd "${ANALYSIS_DIR}"
echo "job_id=${SLURM_JOB_ID} array_job=${SLURM_ARRAY_JOB_ID} chunk=${SLURM_ARRAY_TASK_ID}/${CHUNKS} host=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem=4G seed=${MASTER_SEED} pseudo=${PSEUDO_OBSERVATIONS} permutations=${PERMUTATIONS} exact_draws=${EXACT_DRAWS}"
time -p guix shell python python-numpy -- python3 calibration_suite.py \
  --chunk-index "${SLURM_ARRAY_TASK_ID}" --chunks "${CHUNKS}" \
  --pseudo-observations "${PSEUDO_OBSERVATIONS}" --permutations "${PERMUTATIONS}" \
  --exact-draws "${EXACT_DRAWS}" --seed "${MASTER_SEED}" \
  --output "${RESULT_DIR}/chunk-${SLURM_ARRAY_TASK_ID}.json"
echo "calibration_chunk=${RESULT_DIR}/chunk-${SLURM_ARRAY_TASK_ID}.json"
