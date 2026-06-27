#!/usr/bin/env bash
#SBATCH --job-name=fig5_impg_finalize_2kb
#SBATCH --partition=workers
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --output=/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/finalize_after_arrays.%j.out
#SBATCH --error=/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/finalize_after_arrays.%j.err
set -euo pipefail
export LC_ALL=C

python3 /moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/scripts/finalize_after_slurm_arrays.py \
  --live-base /moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded \
  --main-base /moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded \
  --array-ids 1706840 1706841 1706842 1706843 1706844 1706845
