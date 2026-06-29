#!/usr/bin/env bash
#SBATCH --job-name=fig5_topn_finalize
#SBATCH --partition=workers
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --output=/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_topn_depthcapped/logs/finalize_topn_depthcapped.%j.out
#SBATCH --error=/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_topn_depthcapped/logs/finalize_topn_depthcapped.%j.err
set -euo pipefail
export LC_ALL=C

python3 /moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_topn_depthcapped/scripts/finalize_topn_depthcapped_impg.py
