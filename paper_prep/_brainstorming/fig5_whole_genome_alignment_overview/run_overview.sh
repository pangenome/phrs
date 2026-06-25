#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

python3 scripts/build_whole_genome_alignment_overview.py
Rscript scripts/plot_whole_genome_alignment_overview.R .
