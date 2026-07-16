#!/usr/bin/env bash
#SBATCH --job-name=phr-v6-finalize
#SBATCH --partition=workers
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=results/logs/finalize-%j.out
#SBATCH --error=results/logs/finalize-%j.err

set -euo pipefail

: "${REPO_ROOT:?submit with --export=ALL,REPO_ROOT=/absolute/worktree}"
: "${V6_INITIAL_SLURM_JOB_ID:?export the completed initial-screen job ID}"
ANALYSIS_DIR="${REPO_ROOT}/paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v6"

cd "${ANALYSIS_DIR}"
echo "job_id=${SLURM_JOB_ID} host=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem=${SLURM_MEM_PER_NODE:-unknown} initial_job_id=${V6_INITIAL_SLURM_JOB_ID}"
time -p guix shell python python-numpy -- \
  python3 v6_engine.py finalize
echo "completion_status=SUCCESS output=${ANALYSIS_DIR}/results/release"
