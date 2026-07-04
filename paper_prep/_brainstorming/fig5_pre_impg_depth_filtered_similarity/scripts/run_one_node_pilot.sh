#!/usr/bin/env bash
#SBATCH --job-name=fig5-pre-impg-depth
#SBATCH --partition=tux
#SBATCH --nodes=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=700G
#SBATCH --time=24:00:00
#SBATCH --output=paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/logs/slurm-%j.out
#SBATCH --error=paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/logs/slurm-%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-/moosefs/erikg/phrs}"

OUT="paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity"
COMPARISON_ID="${COMPARISON_ID:-PAN027pat_vs_PAN011_joint}"
BASES="${BASES:-1:1 4:4 10:10}"
WINDOW_SIZE="${WINDOW_SIZE:-2000}"
MAX_WINDOW_DEPTH="${MAX_WINDOW_DEPTH:-100}"
TOP_N="${TOP_N:-20}"
MAX_POST_IMPG_CANDIDATES="${MAX_POST_IMPG_CANDIDATES:-500}"
SCORING="${SWEEPGA_SCORING:-ani}"
RUN_TOPN="${RUN_TOPN:-1}"
RUN_CLASS_WINNERS="${RUN_CLASS_WINNERS:-0}"
THREADS="${SLURM_CPUS_PER_TASK:-96}"
SCRATCH_BASE="${SCRATCH_BASE:-/dev/shm}"

SWEEPGA="${SWEEPGA:-/home/erikg/.cargo/bin/sweepga}"
IMPG="${IMPG:-/home/erikg/.cargo/bin/impg}"
BGZIP="${BGZIP:-/home/erikg/.guix-profile/bin/bgzip}"
PIGZ="${PIGZ:-$(command -v pigz)}"

RAW_PAF="/moosefs/erikg/phrs/.wg-worktrees/agent-2727/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency32/raw_paf/${COMPARISON_ID}.sweepga_frequency32_many_many_j0.paf.gz"
QUERY_FASTA="/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/${COMPARISON_ID}.query.fa"
TARGET_FASTA="/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/${COMPARISON_ID}.target.fa"
QUERY_FAI="${QUERY_FASTA}.fai"
CENTROMERE_BED="data/chm13-annotations.bed"

mkdir -p "$OUT"/{jobs,logs,filtered_paf,beds,outputs,metadata,summaries}

if [[ ! -x "$SWEEPGA" ]]; then echo "missing sweepga: $SWEEPGA" >&2; exit 1; fi
if [[ ! -x "$IMPG" ]]; then echo "missing impg: $IMPG" >&2; exit 1; fi
if [[ ! -x "$BGZIP" ]]; then echo "missing bgzip: $BGZIP" >&2; exit 1; fi
if [[ ! -x "$PIGZ" ]]; then echo "missing pigz: $PIGZ" >&2; exit 1; fi
for required in "$RAW_PAF" "$QUERY_FASTA" "$TARGET_FASTA" "$QUERY_FAI" "$CENTROMERE_BED"; do
  if [[ ! -s "$required" ]]; then echo "missing required input: $required" >&2; exit 1; fi
done

export LC_ALL=C
export RAYON_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"

SCRATCH="$(mktemp -d "$SCRATCH_BASE/fig5_pre_impg.${COMPARISON_ID}.${SLURM_JOB_ID:-manual}.XXXXXX")"
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
  echo "top_n=$TOP_N"
  echo "run_topn=$RUN_TOPN"
  echo "run_class_winners=$RUN_CLASS_WINNERS"
  echo "sweepga=$SWEEPGA"
  "$SWEEPGA" --version
  echo "impg=$IMPG"
  "$IMPG" --version
  echo "raw_paf=$RAW_PAF"
  echo "query_fasta=$QUERY_FASTA"
  echo "target_fasta=$TARGET_FASTA"
  echo "scratch=$SCRATCH"
} | tee "$OUT/metadata/runtime.${SLURM_JOB_ID:-manual}.txt"

RAW_UNCOMPRESSED="$SCRATCH/${COMPARISON_ID}.raw.paf"
"$PIGZ" -dc -p "$THREADS" "$RAW_PAF" > "$RAW_UNCOMPRESSED"

MANIFEST="$OUT/summaries/pre_impg_depth_filtered_manifest.${COMPARISON_ID}.${SLURM_JOB_ID:-manual}.tsv"
printf "comparison_id\tbasis\traw_paf\tfiltered_paf\tquery_bed\tdepth_summary\ttopn_impg_similarity_tsv_gz\ttopn_skip_report\tclass_winner_tsv_gz\tclass_winner_skip_report\tstatus\n" > "$MANIFEST"

for BASIS in $BASES; do
  BASIS_ID="${BASIS/:/to}"
  FILTERED_PAF="$OUT/filtered_paf/${COMPARISON_ID}.sweepga_f32.${BASIS_ID}.noscaffold.ani.paf.gz"
  FILTERED_TMP="$SCRATCH/${COMPARISON_ID}.${BASIS_ID}.filtered.paf"
  QUERY_BED="$OUT/beds/${COMPARISON_ID}.query_${WINDOW_SIZE}bp.${BASIS_ID}.maxdepth${MAX_WINDOW_DEPTH}.bed"
  DEPTH_REPORT="$OUT/summaries/${COMPARISON_ID}.${BASIS_ID}.query_${WINDOW_SIZE}bp.depth_report.tsv"
  DEPTH_SUMMARY="$OUT/summaries/${COMPARISON_ID}.${BASIS_ID}.query_${WINDOW_SIZE}bp.depth_summary.tsv"
  TOPN_OUTPUT_GZ="$OUT/outputs/${COMPARISON_ID}.sweepga_f32.${BASIS_ID}.query_${WINDOW_SIZE}bp.predepth_top${TOP_N}.impg_similarity.tsv.gz"
  TOPN_SKIP_REPORT="$OUT/summaries/${COMPARISON_ID}.${BASIS_ID}.query_${WINDOW_SIZE}bp.impg_topn_skip.tsv"
  CLASS_OUTPUT_GZ="$OUT/outputs/${COMPARISON_ID}.sweepga_f32.${BASIS_ID}.query_${WINDOW_SIZE}bp.predepth_class_winners.impg_similarity.tsv.gz"
  CLASS_SKIP_REPORT="$OUT/summaries/${COMPARISON_ID}.${BASIS_ID}.query_${WINDOW_SIZE}bp.impg_class_winner_skip.tsv"

  echo "[$(date -Is)] sweepga filter basis=$BASIS"
  "$SWEEPGA" \
    --num-mappings "$BASIS" \
    --scaffold-jump 0 \
    --scoring "$SCORING" \
    --temp-dir "$SCRATCH" \
    --output-file "$FILTERED_TMP" \
    "$RAW_UNCOMPRESSED"
  "$BGZIP" -@ "$THREADS" -c "$FILTERED_TMP" > "$FILTERED_PAF.tmp"
  "$BGZIP" -t "$FILTERED_PAF.tmp"
  mv "$FILTERED_PAF.tmp" "$FILTERED_PAF"
  rm -f "$FILTERED_TMP"

  echo "[$(date -Is)] pre-IMPG query-window depth filter basis=$BASIS"
  python3 "$OUT/scripts/filter_query_windows_by_paf_depth.py" \
    --paf "$FILTERED_PAF" \
    --query-fai "$QUERY_FAI" \
    --centromere-bed "$CENTROMERE_BED" \
    --out-bed "$QUERY_BED" \
    --report "$DEPTH_REPORT" \
    --summary "$DEPTH_SUMMARY" \
    --window-size "$WINDOW_SIZE" \
    --min-depth 1 \
    --max-depth "$MAX_WINDOW_DEPTH" \
    --interchrom-only \
    --pigz-threads "$THREADS"

  if [[ "$RUN_TOPN" == "1" ]]; then
    echo "[$(date -Is)] impg top-N similarity basis=$BASIS"
    "$IMPG" similarity \
      --alignment-files "$FILTERED_PAF" \
      --target-bed "$QUERY_BED" \
      --sequence-files "$QUERY_FASTA" "$TARGET_FASTA" \
      --gfa-engine poa \
      --no-merge \
      --num-mappings many:many \
      --scaffold-jump 0 \
      --threads "$THREADS" \
      | python3 "$OUT/scripts/filter_impg_similarity_topn.py" \
          --top-n "$TOP_N" \
          --max-candidates "$MAX_POST_IMPG_CANDIDATES" \
          --interchrom-only \
          --skip-report "$TOPN_SKIP_REPORT" \
      | "$PIGZ" -p "$THREADS" > "$TOPN_OUTPUT_GZ.tmp"
    gzip -t "$TOPN_OUTPUT_GZ.tmp"
    mv "$TOPN_OUTPUT_GZ.tmp" "$TOPN_OUTPUT_GZ"
  else
    TOPN_OUTPUT_GZ="NA"
    TOPN_SKIP_REPORT="NA"
  fi

  if [[ "$RUN_CLASS_WINNERS" == "1" ]]; then
    echo "[$(date -Is)] impg same-chrom/interchrom class winners basis=$BASIS"
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
          --skip-report "$CLASS_SKIP_REPORT" \
      | "$PIGZ" -p "$THREADS" > "$CLASS_OUTPUT_GZ.tmp"
    gzip -t "$CLASS_OUTPUT_GZ.tmp"
    mv "$CLASS_OUTPUT_GZ.tmp" "$CLASS_OUTPUT_GZ"
  else
    CLASS_OUTPUT_GZ="NA"
    CLASS_SKIP_REPORT="NA"
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tOK\n" \
    "$COMPARISON_ID" "$BASIS" "$RAW_PAF" "$FILTERED_PAF" "$QUERY_BED" "$DEPTH_SUMMARY" \
    "$TOPN_OUTPUT_GZ" "$TOPN_SKIP_REPORT" "$CLASS_OUTPUT_GZ" "$CLASS_SKIP_REPORT" >> "$MANIFEST"
done

date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" | tee -a "$OUT/metadata/runtime.${SLURM_JOB_ID:-manual}.txt"
