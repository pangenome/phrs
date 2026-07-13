#!/usr/bin/env bash
#SBATCH --job-name=chm13-cn-bench
#SBATCH --partition=workers
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_logs/benchmark-%j.out
#SBATCH --error=slurm_logs/benchmark-%j.err

set -euo pipefail

: "${REPO_ROOT:?submit with --export=ALL,REPO_ROOT=/absolute/repository/path}"
BENCH_PERMUTATIONS="${BENCH_PERMUTATIONS:-5000}"
BENCH_TERMS="${BENCH_TERMS:-15000}"
ANNOTATIONS_PER_LOCUS="${ANNOTATIONS_PER_LOCUS:-12}"
FINAL_PERMUTATIONS="${FINAL_PERMUTATIONS:-99999}"
MASTER_SEED="${MASTER_SEED:-2026071301}"
BENCH_OUTPUT="${BENCH_OUTPUT:-${SLURM_SUBMIT_DIR}/benchmark_results.slurm-${SLURM_JOB_ID}.json}"
ANALYSIS_DIR="${REPO_ROOT}/paper_prep/_brainstorming/chm13_copy_enrichment"

cd "${ANALYSIS_DIR}"
echo "job_id=${SLURM_JOB_ID} host=$(hostname) cpus=${SLURM_CPUS_PER_TASK} seed=${MASTER_SEED} permutations=${BENCH_PERMUTATIONS} terms=${BENCH_TERMS}"
time -p guix shell python python-numpy -- \
  python3 benchmark_COPY_engine.py \
    --permutations "${BENCH_PERMUTATIONS}" \
    --terms "${BENCH_TERMS}" \
    --annotations-per-locus "${ANNOTATIONS_PER_LOCUS}" \
    --final-permutations "${FINAL_PERMUTATIONS}" \
    --seed "${MASTER_SEED}" \
    --output "${BENCH_OUTPUT}"
echo "benchmark_output=${BENCH_OUTPUT}"
