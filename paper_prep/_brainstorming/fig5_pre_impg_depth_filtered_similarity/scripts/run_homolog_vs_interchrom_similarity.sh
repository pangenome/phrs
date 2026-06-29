#!/usr/bin/env bash
#SBATCH --job-name=fig5-impg-hom-int
#SBATCH --partition=tux
#SBATCH --nodes=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=700G
#SBATCH --time=12:00:00
#SBATCH --output=paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/logs/homolog_vs_interchrom-%j.out
#SBATCH --error=paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/logs/homolog_vs_interchrom-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-/moosefs/erikg/phrs}"

OUT="paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity"
COMPARISON_ID="${COMPARISON_ID:-PAN027pat_vs_PAN011_joint}"
BASES="${BASES:-1:1 4:4 10:10}"
WINDOW_SIZE="${WINDOW_SIZE:-2000}"
MAX_WINDOW_DEPTH="${MAX_WINDOW_DEPTH:-100}"
MAX_POST_IMPG_CANDIDATES="${MAX_POST_IMPG_CANDIDATES:-5000}"
THREADS="${SLURM_CPUS_PER_TASK:-96}"
SCRATCH_BASE="${SCRATCH_BASE:-/dev/shm}"

IMPG="${IMPG:-/home/erikg/.cargo/bin/impg}"
PIGZ="${PIGZ:-$(command -v pigz)}"

QUERY_FASTA="/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/${COMPARISON_ID}.query.fa"
TARGET_FASTA="/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/${COMPARISON_ID}.target.fa"

mkdir -p "$OUT"/{logs,metadata,outputs,summaries}

if [[ ! -x "$IMPG" ]]; then echo "missing impg: $IMPG" >&2; exit 1; fi
if [[ ! -x "$PIGZ" ]]; then echo "missing pigz: $PIGZ" >&2; exit 1; fi
for required in "$QUERY_FASTA" "$TARGET_FASTA"; do
  if [[ ! -s "$required" ]]; then echo "missing required input: $required" >&2; exit 1; fi
done

export LC_ALL=C
export RAYON_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"

SCRATCH="$(mktemp -d "$SCRATCH_BASE/fig5_impg_class.${COMPARISON_ID}.${SLURM_JOB_ID:-manual}.XXXXXX")"
trap 'rm -rf "$SCRATCH"' EXIT

{
  date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ"
  echo "hostname=$(hostname)"
  echo "slurm_job_id=${SLURM_JOB_ID:-manual}"
  echo "slurm_cpus_per_task=$THREADS"
  echo "comparison_id=$COMPARISON_ID"
  echo "bases=$BASES"
  echo "window_size=$WINDOW_SIZE"
  echo "max_window_depth=$MAX_WINDOW_DEPTH"
  echo "max_post_impg_candidates=$MAX_POST_IMPG_CANDIDATES"
  echo "impg=$IMPG"
  "$IMPG" --version
  echo "query_fasta=$QUERY_FASTA"
  echo "target_fasta=$TARGET_FASTA"
  echo "scratch=$SCRATCH"
} | tee "$OUT/metadata/homolog_vs_interchrom_runtime.${SLURM_JOB_ID:-manual}.txt"

MANIFEST="$OUT/summaries/homolog_vs_interchrom_manifest.${COMPARISON_ID}.${SLURM_JOB_ID:-manual}.tsv"
printf "comparison_id\tbasis\tfiltered_paf\tquery_bed\tclass_winner_tsv_gz\tskip_report\tstatus\n" > "$MANIFEST"

for BASIS in $BASES; do
  BASIS_ID="${BASIS/:/to}"
  FILTERED_PAF="$OUT/filtered_paf/${COMPARISON_ID}.sweepga_f32.${BASIS_ID}.noscaffold.ani.paf.gz"
  QUERY_BED="$OUT/beds/${COMPARISON_ID}.query_${WINDOW_SIZE}bp.${BASIS_ID}.maxdepth${MAX_WINDOW_DEPTH}.bed"
  OUTPUT_GZ="$OUT/outputs/${COMPARISON_ID}.sweepga_f32.${BASIS_ID}.query_${WINDOW_SIZE}bp.predepth_class_winners.impg_similarity.tsv.gz"
  SKIP_REPORT="$OUT/summaries/${COMPARISON_ID}.${BASIS_ID}.query_${WINDOW_SIZE}bp.impg_class_winner_skip.tsv"

  for required in "$FILTERED_PAF" "$QUERY_BED"; do
    if [[ ! -s "$required" ]]; then echo "missing required input: $required" >&2; exit 1; fi
  done

  echo "[$(date -Is)] impg homolog-vs-interchrom class winners basis=$BASIS"
  "$IMPG" similarity \
    --alignment-files "$FILTERED_PAF" \
    --target-bed "$QUERY_BED" \
    --sequence-files "$QUERY_FASTA" "$TARGET_FASTA" \
    --gfa-engine poa \
    --no-merge \
    --num-mappings many:many \
    --scaffold-jump 0 \
    --threads "$THREADS" \
    | python3 "$OUT/scripts/filter_impg_similarity_class_winners.py" \
        --max-candidates "$MAX_POST_IMPG_CANDIDATES" \
        --skip-report "$SKIP_REPORT" \
    | "$PIGZ" -p "$THREADS" > "$OUTPUT_GZ.tmp"
  gzip -t "$OUTPUT_GZ.tmp"
  mv "$OUTPUT_GZ.tmp" "$OUTPUT_GZ"

  printf "%s\t%s\t%s\t%s\t%s\t%s\tOK\n" \
    "$COMPARISON_ID" "$BASIS" "$FILTERED_PAF" "$QUERY_BED" "$OUTPUT_GZ" "$SKIP_REPORT" >> "$MANIFEST"
done

date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" | tee -a "$OUT/metadata/homolog_vs_interchrom_runtime.${SLURM_JOB_ID:-manual}.txt"
