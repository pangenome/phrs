#!/usr/bin/env bash
#SBATCH --job-name=phr-v6-targeted
#SBATCH --partition=workers
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=results/logs/targeted-%j.out
#SBATCH --error=results/logs/targeted-%j.err

set -euo pipefail
: "${REPO_ROOT:?submit with --export=ALL,REPO_ROOT=/absolute/worktree}"
ANALYSIS_DIR="${REPO_ROOT}/paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v6"
cd "${ANALYSIS_DIR}"
echo "job_id=${SLURM_JOB_ID} host=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem=${SLURM_MEM_PER_NODE:-unknown} seed=2026071301"
time -p guix shell python python-numpy -- python3 v6_targeted_extension.py
echo "completion_status=SUCCESS output=${ANALYSIS_DIR}/results/release/SELECTIVE_EXTENSION_STATUS.json"
