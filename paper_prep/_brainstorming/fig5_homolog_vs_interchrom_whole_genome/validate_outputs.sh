#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
OUT_DIR="$ROOT/paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome"

Rscript "$OUT_DIR/scripts/plot_homolog_vs_interchrom_whole_genome.R" "$ROOT" "$OUT_DIR"

test -s "$OUT_DIR/fig5_homolog_vs_interchrom_whole_genome.pdf"
test -s "$OUT_DIR/fig5_homolog_vs_interchrom_whole_genome.png"
test -s "$OUT_DIR/fig5_homolog_vs_interchrom_whole_genome.svg"
test -s "$OUT_DIR/context_intervals.tsv"
test -s "$OUT_DIR/top_pair_summary_for_plot.tsv"

grep -q 'PAR1_XY_positive_control' "$OUT_DIR/context_intervals.tsv"
grep -q 'PAN027_chr9q_chr3q_PHR_candidate' "$OUT_DIR/context_intervals.tsv"
grep -q 'PAN028_chr9q_chr3q_PHR_candidate' "$OUT_DIR/context_intervals.tsv"
grep -q 'PHR-acro' "$OUT_DIR/context_intervals.tsv"
grep -q 'chr9' "$OUT_DIR/context_intervals.tsv"
grep -q '^10:10' "$OUT_DIR/top_pair_summary_for_plot.tsv"
grep -q 'chrX	chrY' "$OUT_DIR/top_pair_summary_for_plot.tsv"
