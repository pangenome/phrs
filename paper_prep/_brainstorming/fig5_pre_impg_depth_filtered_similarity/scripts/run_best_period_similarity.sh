#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

OUT="paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity"
COMPARISON_ID="${COMPARISON_ID:-PAN027pat_vs_PAN011_joint}"
BASES="${BASES:-1:1 4:4 10:10}"
WINDOW_SIZE="${WINDOW_SIZE:-2000}"
MAX_WINDOW_DEPTH="${MAX_WINDOW_DEPTH:-100}"
MAX_POST_IMPG_CANDIDATES="${MAX_POST_IMPG_CANDIDATES:-5000}"
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-48}}"

SOURCE_ROOT="${SOURCE_ROOT:-/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity}"
INPUT_ROOT="${INPUT_ROOT:-/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs}"
IMPG="${IMPG:-/home/erikg/.cargo/bin/impg}"
PIGZ="${PIGZ:-$(command -v pigz)}"

QUERY_FASTA="${QUERY_FASTA:-$INPUT_ROOT/${COMPARISON_ID}.query.fa}"
TARGET_FASTA="${TARGET_FASTA:-$INPUT_ROOT/${COMPARISON_ID}.target.fa}"

mkdir -p "$OUT"/{outputs,summaries,metadata}

for required in "$IMPG" "$PIGZ" "$QUERY_FASTA" "$TARGET_FASTA"; do
  if [[ ! -e "$required" ]]; then
    echo "missing required input/tool: $required" >&2
    exit 1
  fi
done
if [[ ! -x "$IMPG" ]]; then echo "missing executable impg: $IMPG" >&2; exit 1; fi
if [[ ! -x "$PIGZ" ]]; then echo "missing executable pigz: $PIGZ" >&2; exit 1; fi

export LC_ALL=C
export RAYON_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"

MANIFEST="$OUT/summaries/pre_impg_best_all_manifest.${COMPARISON_ID}.$(date -u +%Y%m%dT%H%M%SZ).tsv"
printf "comparison_id\tbasis\tfiltered_paf\tquery_bed\timpg_similarity_tsv_gz\tskip_report\tstatus\n" > "$MANIFEST"

{
  date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ"
  echo "comparison_id=$COMPARISON_ID"
  echo "bases=$BASES"
  echo "window_size=$WINDOW_SIZE"
  echo "max_window_depth=$MAX_WINDOW_DEPTH"
  echo "max_post_impg_candidates=$MAX_POST_IMPG_CANDIDATES"
  echo "threads=$THREADS"
  echo "source_root=$SOURCE_ROOT"
  echo "query_fasta=$QUERY_FASTA"
  echo "target_fasta=$TARGET_FASTA"
  echo "impg=$IMPG"
  "$IMPG" --version
} | tee "$OUT/metadata/best_all_runtime.$(date -u +%Y%m%dT%H%M%SZ).txt"

for BASIS in $BASES; do
  BASIS_ID="${BASIS/:/to}"
  FILTERED_PAF="$SOURCE_ROOT/filtered_paf/${COMPARISON_ID}.sweepga_f32.${BASIS_ID}.noscaffold.ani.paf.gz"
  QUERY_BED="$SOURCE_ROOT/beds/${COMPARISON_ID}.query_${WINDOW_SIZE}bp.${BASIS_ID}.maxdepth${MAX_WINDOW_DEPTH}.bed"
  OUTPUT_GZ="$OUT/outputs/${COMPARISON_ID}.sweepga_f32.${BASIS_ID}.query_${WINDOW_SIZE}bp.predepth_best_all.impg_similarity.tsv.gz"
  SKIP_REPORT="$OUT/summaries/${COMPARISON_ID}.${BASIS_ID}.query_${WINDOW_SIZE}bp.impg_best_all_skip.tsv"

  for required in "$FILTERED_PAF" "$QUERY_BED"; do
    if [[ ! -s "$required" ]]; then
      echo "missing required existing SweepGA/pre-IMPG artifact: $required" >&2
      exit 1
    fi
  done

  echo "[$(date -Is)] impg best-all similarity basis=$BASIS"
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

date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" | tee -a "$OUT/metadata/best_all_runtime.latest.txt"
