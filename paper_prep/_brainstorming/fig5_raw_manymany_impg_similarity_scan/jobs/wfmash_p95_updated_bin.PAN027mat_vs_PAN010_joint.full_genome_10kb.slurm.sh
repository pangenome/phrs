#!/usr/bin/env bash
#SBATCH --job-name=impg_wfmash_PAN027mat_vs_P
#SBATCH --partition=workers
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --output=/moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/logs/wfmash_p95_updated_bin.PAN027mat_vs_PAN010_joint.full_genome_10kb.%j.out
#SBATCH --error=/moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/logs/wfmash_p95_updated_bin.PAN027mat_vs_PAN010_joint.full_genome_10kb.%j.err
set -euo pipefail
export LC_ALL=C
mkdir -p /moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/outputs
tmp_tsv=/moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/outputs/wfmash_p95_updated_bin.PAN027mat_vs_PAN010_joint.full_genome_10kb.impg_similarity.tsv
out_gz=/moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/outputs/wfmash_p95_updated_bin.PAN027mat_vs_PAN010_joint.full_genome_10kb.impg_similarity.tsv.gz
metadata=/moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/outputs/wfmash_p95_updated_bin.PAN027mat_vs_PAN010_joint.full_genome_10kb.metadata.json
command_text='/home/erikg/.cargo/bin/impg similarity --alignment-files /moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/raw_paf/updated_bin_v0.24.2-12-ge040aa10/PAN027mat_vs_PAN010_joint.literal_p95.wfmash-v0.24.2-12-ge040aa10.paf.gz --target-bed /moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/beds/PAN027mat_vs_PAN010_joint.full_genome_10kb.bed --sequence-files /moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027mat_vs_PAN010_joint.query.fa /moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027mat_vs_PAN010_joint.target.fa --gfa-engine poa --no-merge --num-mappings many:many --scaffold-jump 0 --threads ${SLURM_CPUS_PER_TASK}'
start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
impg_version=$(/home/erikg/.cargo/bin/impg --version)
node=${SLURMD_NODENAME:-$(hostname)}
partition=${SLURM_JOB_PARTITION:-unknown}
cat > "$metadata" <<JSON
{
  "method": "wfmash_p95_updated_bin",
  "comparison_id": "PAN027mat_vs_PAN010_joint",
  "source_raw_paf": "/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/raw_paf/updated_bin_v0.24.2-12-ge040aa10/PAN027mat_vs_PAN010_joint.literal_p95.wfmash-v0.24.2-12-ge040aa10.paf.gz",
  "impg_alignment_paf": "/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/raw_paf/updated_bin_v0.24.2-12-ge040aa10/PAN027mat_vs_PAN010_joint.literal_p95.wfmash-v0.24.2-12-ge040aa10.paf.gz",
  "query_fasta": "/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027mat_vs_PAN010_joint.query.fa",
  "target_fasta": "/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027mat_vs_PAN010_joint.target.fa",
  "target_bed": "/moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/beds/PAN027mat_vs_PAN010_joint.full_genome_10kb.bed",
  "output_tsv_gz": "/moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/outputs/wfmash_p95_updated_bin.PAN027mat_vs_PAN010_joint.full_genome_10kb.impg_similarity.tsv.gz",
  "impg_path": "/home/erikg/.cargo/bin/impg",
  "impg_version": "$impg_version",
  "slurm_job_id": "${SLURM_JOB_ID:-manual}",
  "slurm_cpus_per_task": "${SLURM_CPUS_PER_TASK:-unset}",
  "node": "$node",
  "partition": "$partition",
  "start_utc": "$start_utc",
  "command": "$command_text"
}
JSON
echo "$command_text"
/home/erikg/.cargo/bin/impg similarity --alignment-files /moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/raw_paf/updated_bin_v0.24.2-12-ge040aa10/PAN027mat_vs_PAN010_joint.literal_p95.wfmash-v0.24.2-12-ge040aa10.paf.gz --target-bed /moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/beds/PAN027mat_vs_PAN010_joint.full_genome_10kb.bed --sequence-files /moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027mat_vs_PAN010_joint.query.fa /moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027mat_vs_PAN010_joint.target.fa --gfa-engine poa --no-merge --num-mappings many:many --scaffold-jump 0 --threads ${SLURM_CPUS_PER_TASK} > "$tmp_tsv"
gzip -f "$tmp_tsv"
finish_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
python3 - <<'PY' "$metadata" "$finish_utc"
import json, sys
path, finish = sys.argv[1:]
with open(path) as fh:
    data = json.load(fh)
data['finish_utc'] = finish
data['status'] = 'OK'
with open(path, 'w') as fh:
    json.dump(data, fh, indent=2, sort_keys=True)
    fh.write('\n')
PY
