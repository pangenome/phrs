#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

/home/erikg/.guix-profile/bin/Rscript \
  paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome/scripts/plot_whole_genome_homolog_vs_interchrom.R

