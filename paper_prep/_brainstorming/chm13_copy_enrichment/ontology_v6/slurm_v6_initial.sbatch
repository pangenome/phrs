#!/usr/bin/env bash
#SBATCH --job-name=phr-v6-copy-ontology
#SBATCH --partition=workers
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=results/logs/initial-%j.out
#SBATCH --error=results/logs/initial-%j.err

set -euo pipefail

: "${REPO_ROOT:?submit with --export=ALL,REPO_ROOT=/absolute/worktree}"
ANALYSIS_DIR="${REPO_ROOT}/paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v6"
PERMUTATIONS="${PERMUTATIONS:-99999}"
BATCH_SIZE="${BATCH_SIZE:-1000}"

if [[ "${PERMUTATIONS}" != "99999" ]]; then
  echo "V6 initial screen is frozen at exactly 99,999 permutations" >&2
  exit 2
fi

mkdir -p "${ANALYSIS_DIR}/results/logs"
cd "${ANALYSIS_DIR}"
echo "job_id=${SLURM_JOB_ID} host=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem=${SLURM_MEM_PER_NODE:-unknown} permutations=${PERMUTATIONS} batch_size=${BATCH_SIZE}"
time -p guix shell python python-numpy -- \
  python3 v6_engine.py run-initial \
    --permutations "${PERMUTATIONS}" \
    --batch-size "${BATCH_SIZE}"
echo "completion_status=SUCCESS output=${ANALYSIS_DIR}/results/work/primary"
