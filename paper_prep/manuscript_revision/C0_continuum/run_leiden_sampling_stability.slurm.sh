#!/usr/bin/env bash
# Slurm-ready specification for the missing D-1 sampling-stability test.
# This is not submitted by default. It documents the intended cluster run and
# fails early unless the worker supplies a concrete implementation script.

#SBATCH --job-name=phrs-c0c-leiden-stability
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=paper_prep/manuscript_revision/C0_continuum/logs/leiden-stability-%j.out
#SBATCH --error=paper_prep/manuscript_revision/C0_continuum/logs/leiden-stability-%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/moosefs/erikg/phrs}"
OUT_DIR="${OUT_DIR:-${REPO_ROOT}/paper_prep/manuscript_revision/C0_continuum/leiden_sampling_stability}"
SCRIPT="${SCRIPT:-${REPO_ROOT}/paper_prep/manuscript_revision/C0_continuum/run_leiden_sampling_stability.R}"

mkdir -p "${OUT_DIR}" "$(dirname "${SLURM_OUTPUT:-paper_prep/manuscript_revision/C0_continuum/logs/placeholder}")"

cat >&2 <<EOF
[C0c/D1] This wrapper is a run specification, not a completed analysis.

Required implementation:
  ${SCRIPT}

Minimum expected outputs:
  ${OUT_DIR}/arm_leiden_phr_bootstrap_partition_stability.tsv
  ${OUT_DIR}/arm_leiden_phr_bootstrap_pair_coclustering.tsv
  ${OUT_DIR}/sequence_leiden_sample_subsampling_stability.tsv

Design:
  1. Resample PHRs or samples/haplotypes.
  2. Recompute the arm-level distance matrix from the canonical per-PHR or
     per-sequence similarity contributions, not from the fixed 41 x 41 matrix.
  3. Re-run the same Leiden scan/selection rule as the baseline.
  4. Report ARI/NMI versus the baseline assignments, selected resolution,
     selected number of communities, and per-pair co-clustering probability.
  5. For sequence-level stability, subsample samples/haplotypes, rebuild the
     k-NN graph over retained PHRs, and compare projected labels to the
     baseline k=75, resolution=0.8, 50-community partition.
EOF

if [[ ! -s "${SCRIPT}" ]]; then
  echo "[C0c/D1] Missing implementation script: ${SCRIPT}" >&2
  exit 2
fi

cd "${REPO_ROOT}"
Rscript "${SCRIPT}" --out-dir "${OUT_DIR}"
