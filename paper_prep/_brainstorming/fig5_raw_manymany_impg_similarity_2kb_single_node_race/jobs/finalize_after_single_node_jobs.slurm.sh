#!/usr/bin/env bash
#SBATCH --job-name=fig5_1node_finalize
#SBATCH --partition=tux
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --output=/moosefs/erikg/phrs/.wg-worktrees/agent-2859/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_single_node_race/logs/finalize_after_single_node.%j.out
#SBATCH --error=/moosefs/erikg/phrs/.wg-worktrees/agent-2859/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_single_node_race/logs/finalize_after_single_node.%j.err
set -euo pipefail
export LC_ALL=C
python3 /moosefs/erikg/phrs/.wg-worktrees/agent-2859/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_single_node_race/scripts/finalize_single_node_race.py
