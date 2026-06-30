#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

python3 paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft/scripts/plot_whole_genome_ribbon_draft.py
