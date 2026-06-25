#!/usr/bin/env bash
#SBATCH --job-name=bgzf_PAN027pat_vs_P
#SBATCH --partition=workers
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --output=/moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/logs/normalize_bgzf.PAN027pat_vs_PAN011_joint.%j.out
#SBATCH --error=/moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/logs/normalize_bgzf.PAN027pat_vs_PAN011_joint.%j.err
set -euo pipefail
export LC_ALL=C
source_paf=/moosefs/erikg/phrs/.wg-worktrees/agent-2727/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency32/raw_paf/PAN027pat_vs_PAN011_joint.sweepga_frequency32_many_many_j0.paf.gz
dest_paf=/moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/bgzf_raw_paf/PAN027pat_vs_PAN011_joint.sweepga_frequency32_many_many_j0.bgzf.paf.gz
tmp_paf="${dest_paf}.tmp.${SLURM_JOB_ID:-manual}"
if [ -s "$dest_paf" ]; then
  /home/erikg/.guix-profile/bin/bgzip -t "$dest_paf"
  echo "BGZF already exists: $dest_paf"
  exit 0
fi
gzip -dc "$source_paf" | /home/erikg/.guix-profile/bin/bgzip -@ ${SLURM_CPUS_PER_TASK} -c > "$tmp_paf"
/home/erikg/.guix-profile/bin/bgzip -t "$tmp_paf"
mv "$tmp_paf" "$dest_paf"
