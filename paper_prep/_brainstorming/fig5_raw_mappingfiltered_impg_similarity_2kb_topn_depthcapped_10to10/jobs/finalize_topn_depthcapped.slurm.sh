#!/usr/bin/env bash
#SBATCH --job-name=fig5_topn_10to10_finalize
#SBATCH --partition=workers
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --output=/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_mappingfiltered_impg_similarity_2kb_topn_depthcapped_10to10/logs/finalize_topn_depthcapped.%j.out
#SBATCH --error=/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_mappingfiltered_impg_similarity_2kb_topn_depthcapped_10to10/logs/finalize_topn_depthcapped.%j.err
set -euo pipefail
export LC_ALL=C

python3 /moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_mappingfiltered_impg_similarity_2kb_topn_depthcapped_10to10/scripts/finalize_topn_depthcapped_impg.py
