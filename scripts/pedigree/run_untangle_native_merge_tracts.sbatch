#!/usr/bin/env bash
#SBATCH --job-name=phr-native-untangle
#SBATCH --partition=workers
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/slurm-%j.out
#SBATCH --error=/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/slurm-%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/moosefs/erikg/phrs/.wg-worktrees/agent-2556}"
EXTERNAL_OUT="${EXTERNAL_OUT:-/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm}"
THREADS="${SLURM_CPUS_PER_TASK:-24}"
ODGI="${ODGI:-/home/erikg/.guix-profile/bin/odgi}"
GRAPH="${GRAPH:-/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/pggb/washu.1Mb.telo_500kb_trimmed.fa.gz.6e0e250.11fba48.13f423a.smooth.final.og}"
UNTANGLE_DIR="${UNTANGLE_DIR:-/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle}"
CUT_EVERY="${CUT_EVERY:-50000}"
MIN_JACCARD="${MIN_JACCARD:-0.8}"
EMIT_N_BEST="${EMIT_N_BEST:-4}"

mkdir -p "$EXTERNAL_OUT"
cd "$REPO_ROOT"

make_queries() {
  local pair="$1"
  local child="$2"
  local child_hap="$3"
  awk -v prefix="${child}#${child_hap}#" 'index($0, prefix) == 1 { print }' \
    "$UNTANGLE_DIR/all_paths.txt" > "$EXTERNAL_OUT/queries_${pair}.txt"
}

run_bedpe() {
  local pair="$1"
  local merge_dist="$2"
  local query_list="$EXTERNAL_OUT/queries_${pair}.txt"
  local target_list="$UNTANGLE_DIR/targets_${pair}.txt"
  local out="$EXTERNAL_OUT/${pair}.e${CUT_EVERY}.m${merge_dist}.j${MIN_JACCARD}.n${EMIT_N_BEST}.bedpe.gz"
  local log="${out}.log"

  {
    printf 'command\t'
    printf '%q ' "$ODGI" untangle -t "$THREADS" -i "$GRAPH" -Q "$query_list" -R "$target_list" \
      -e "$CUT_EVERY" -m "$merge_dist" -j "$MIN_JACCARD" -n "$EMIT_N_BEST"
    printf '| gzip -c > %q\n' "$out"
  } > "$log"

  "$ODGI" untangle \
    -t "$THREADS" \
    -i "$GRAPH" \
    -Q "$query_list" \
    -R "$target_list" \
    -e "$CUT_EVERY" \
    -m "$merge_dist" \
    -j "$MIN_JACCARD" \
    -n "$EMIT_N_BEST" \
    2>> "$log" | gzip -c > "$out"
  printf '\nexit_status\t0\n' >> "$log"
}

run_paf() {
  local pair="$1"
  local merge_dist="$2"
  local query_list="$EXTERNAL_OUT/queries_${pair}.txt"
  local target_list="$UNTANGLE_DIR/targets_${pair}.txt"
  local out="$EXTERNAL_OUT/${pair}.e${CUT_EVERY}.m${merge_dist}.j${MIN_JACCARD}.n${EMIT_N_BEST}.paf.gz"
  local log="${out}.log"

  {
    printf 'command\t'
    printf '%q ' "$ODGI" untangle -t "$THREADS" -i "$GRAPH" -Q "$query_list" -R "$target_list" \
      -e "$CUT_EVERY" -m "$merge_dist" -j "$MIN_JACCARD" -n "$EMIT_N_BEST" -p
    printf '| gzip -c > %q\n' "$out"
  } > "$log"

  "$ODGI" untangle \
    -t "$THREADS" \
    -i "$GRAPH" \
    -Q "$query_list" \
    -R "$target_list" \
    -e "$CUT_EVERY" \
    -m "$merge_dist" \
    -j "$MIN_JACCARD" \
    -n "$EMIT_N_BEST" \
    -p \
    2>> "$log" | gzip -c > "$out"
  printf '\nexit_status\t0\n' >> "$log"
}

make_queries PAN027_vs_PAN010 PAN027 1
make_queries PAN027_vs_PAN011 PAN027 2
make_queries PAN028_vs_PAN027 PAN028 1

for pair in PAN027_vs_PAN010 PAN027_vs_PAN011 PAN028_vs_PAN027; do
  for merge_dist in 0 1000; do
    run_bedpe "$pair" "$merge_dist"
    run_paf "$pair" "$merge_dist"
  done
done
