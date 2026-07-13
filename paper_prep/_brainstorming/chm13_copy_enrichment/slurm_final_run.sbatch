#!/usr/bin/env bash
#SBATCH --job-name=chm13-cn-final
#SBATCH --partition=workers
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_logs/final-%j.out
#SBATCH --error=slurm_logs/final-%j.err

set -euo pipefail

: "${REPO_ROOT:?submit with --export=ALL,REPO_ROOT=/absolute/repository/path}"
: "${MODE:?set MODE to primary, terminal, or adjacent}"
: "${MASTER_SEED:?set the prespecified decimal seed}"
: "${OUTPUT_DIR:?set a new absolute output directory}"
: "${ARGS_FILE:?set an absolute file containing one additional CLI argument per line}"
PERMUTATIONS="${PERMUTATIONS:-99999}"
BATCH_SIZE="${BATCH_SIZE:-1000}"
ANALYSIS_DIR="${REPO_ROOT}/paper_prep/_brainstorming/chm13_copy_enrichment"

case "${MODE}:${MASTER_SEED}" in
  primary:2026071301|terminal:2026071302|adjacent:2026071303) ;;
  *) echo "mode/seed does not match the prespecified seed table" >&2; exit 2 ;;
esac
if (( PERMUTATIONS < 99999 )); then
  echo "final Slurm runs require at least 99,999 permutations" >&2
  exit 2
fi
if [[ ! -s "${ARGS_FILE}" ]]; then
  echo "ARGS_FILE is missing or empty: ${ARGS_FILE}" >&2
  exit 2
fi
mapfile -t EXTRA_ARGS < "${ARGS_FILE}"

cd "${ANALYSIS_DIR}"
echo "job_id=${SLURM_JOB_ID} host=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mode=${MODE} seed=${MASTER_SEED} permutations=${PERMUTATIONS} batch_size=${BATCH_SIZE}"
time -p guix shell python python-numpy -- \
  python3 COPY_engine.py \
    --mode "${MODE}" \
    --seed "${MASTER_SEED}" \
    --permutations "${PERMUTATIONS}" \
    --batch-size "${BATCH_SIZE}" \
    --output "${OUTPUT_DIR}" \
    "${EXTRA_ARGS[@]}"
echo "analysis_output=${OUTPUT_DIR}"
