#!/usr/bin/env bash
#SBATCH --job-name=fig5_f16_multiway
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --output=logs/fig5_f16_multiway.%j.out
#SBATCH --error=logs/fig5_f16_multiway.%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PANEL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SOURCE_PACKAGE="${SOURCE_PACKAGE:-/moosefs/erikg/phrs/.wg-worktrees/agent-2649/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16}"

mkdir -p "$PANEL_DIR/logs" "$PANEL_DIR/preview_png"

python3 "$PANEL_DIR/scripts/extract_multiway_support.py" \
  --panel-dir "$PANEL_DIR" \
  --source-package "$SOURCE_PACKAGE"

Rscript "$PANEL_DIR/scripts/plot_multiway_panels.R" "$PANEL_DIR"
