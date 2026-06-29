#!/usr/bin/env bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

SCRIPT="paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/scripts/run_one_node_pilot.sh"
sbatch --parsable "$SCRIPT"

