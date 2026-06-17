#!/bin/bash
#SBATCH --job-name=c0b_seq_continuum
#SBATCH --partition=allnodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=08:00:00
#SBATCH --output=/moosefs/erikg/phrs/.wg-worktrees/agent-2465/paper_prep/manuscript_revision/C0_continuum/logs/c0b_seq_continuum.%j.out
#SBATCH --error=/moosefs/erikg/phrs/.wg-worktrees/agent-2465/paper_prep/manuscript_revision/C0_continuum/logs/c0b_seq_continuum.%j.err

set -euo pipefail

cd /moosefs/erikg/phrs/.wg-worktrees/agent-2465

echo "=== C0b sequence-level continuum job ==="
echo "Started: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Node: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-NA}"
echo "Memory request: 80G"

DIST_RDS="/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.dist_matrix.rds"
ASSIGNMENTS="/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv"
OUTDIR="paper_prep/manuscript_revision/C0_continuum"

echo "Input RDS:"
ls -lh "${DIST_RDS}"
echo "Assignment table:"
ls -lh "${ASSIGNMENTS}"

Rscript "${OUTDIR}/scripts/analyze_sequence_continuum.R" \
  --dist-rds "${DIST_RDS}" \
  --assignments "${ASSIGNMENTS}" \
  --outdir "${OUTDIR}" \
  --bin-width 0.005

echo "Outputs:"
find "${OUTDIR}/results" -maxdepth 1 -type f -printf '%s\t%p\n' | sort -nr
ls -lh "${OUTDIR}/C0b_sequence_level_report.md"
echo "Completed: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
