#!/usr/bin/env bash
#SBATCH --job-name=fig5_1node_wfmash_PAN028mat_vs_P
#SBATCH --partition=tux
#SBATCH --cpus-per-task=96
#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --output=/moosefs/erikg/phrs/.wg-worktrees/agent-2859/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_single_node_race/logs/wfmash_p95_updated_bin.PAN028mat_vs_PAN027_joint.%j.out
#SBATCH --error=/moosefs/erikg/phrs/.wg-worktrees/agent-2859/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_single_node_race/logs/wfmash_p95_updated_bin.PAN028mat_vs_PAN027_joint.%j.err
set -euo pipefail
export LC_ALL=C
method=wfmash_p95_updated_bin
comparison_id=PAN028mat_vs_PAN027_joint
source_raw_paf=/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/raw_paf/updated_bin_v0.24.2-12-ge040aa10/PAN028mat_vs_PAN027_joint.literal_p95.wfmash-v0.24.2-12-ge040aa10.paf.gz
impg_alignment_paf=/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/raw_paf/updated_bin_v0.24.2-12-ge040aa10/PAN028mat_vs_PAN027_joint.literal_p95.wfmash-v0.24.2-12-ge040aa10.paf.gz
query_fasta=/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN028mat_vs_PAN027_joint.query.fa
target_fasta=/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN028mat_vs_PAN027_joint.target.fa
full_target_bed=/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/beds/PAN028mat_vs_PAN027_joint.full_genome_2kb.bed
output_tsv_gz=/moosefs/erikg/phrs/.wg-worktrees/agent-2859/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_single_node_race/outputs/all_hits/wfmash_p95_updated_bin/wfmash_p95_updated_bin.PAN028mat_vs_PAN027_joint.full_genome_2kb.impg_similarity.tsv.gz
metadata_json=/moosefs/erikg/phrs/.wg-worktrees/agent-2859/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_single_node_race/metadata/jobs/wfmash_p95_updated_bin/wfmash_p95_updated_bin.PAN028mat_vs_PAN027_joint.full_genome_2kb.metadata.json
impg=/home/erikg/.cargo/bin/impg
mkdir -p "$(dirname "$output_tsv_gz")" "$(dirname "$metadata_json")"
tmp_tsv="${output_tsv_gz%.gz}.tmp.${SLURM_JOB_ID:-manual}"
start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
node="${SLURMD_NODENAME:-$(hostname)}"
partition="${SLURM_JOB_PARTITION:-unknown}"
impg_version=$("$impg" --version)
command_text="$impg similarity --alignment-files $impg_alignment_paf --target-bed $full_target_bed --sequence-files $query_fasta $target_fasta --gfa-engine poa --no-merge --num-mappings many:many --scaffold-jump 0 --threads ${SLURM_CPUS_PER_TASK}"
python3 - "$metadata_json" "$method" "$comparison_id" "$source_raw_paf" "$impg_alignment_paf" "$query_fasta" "$target_fasta" "$full_target_bed" "$output_tsv_gz" "$impg" "$impg_version" "$command_text" "$start_utc" "$node" "$partition" <<'PY'
import json, os, sys
keys = ['metadata_json','method','comparison_id','source_raw_paf','impg_alignment_paf','query_fasta','target_fasta','full_target_bed','output_tsv_gz','impg_path','impg_version','command','start_utc','node','partition']
data = dict(zip(keys, sys.argv[1:]))
data.update({
    'slurm_job_id': os.environ.get('SLURM_JOB_ID', 'manual'),
    'slurm_cpus_per_task': os.environ.get('SLURM_CPUS_PER_TASK', 'unset'),
    'status': 'RUNNING',
})
with open(data.pop('metadata_json'), 'w') as handle:
    json.dump(data, handle, indent=2, sort_keys=True)
    handle.write('\n')
PY
echo "$command_text"
"$impg" similarity \
  --alignment-files "$impg_alignment_paf" \
  --target-bed "$full_target_bed" \
  --sequence-files "$query_fasta" "$target_fasta" \
  --gfa-engine poa \
  --no-merge \
  --num-mappings many:many \
  --scaffold-jump 0 \
  --threads "${SLURM_CPUS_PER_TASK}" > "$tmp_tsv"
gzip -f "$tmp_tsv"
mv "${tmp_tsv}.gz" "$output_tsv_gz"
finish_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
python3 - "$metadata_json" "$finish_utc" <<'PY'
import json, sys
path, finish = sys.argv[1:]
with open(path) as handle:
    data = json.load(handle)
data['finish_utc'] = finish
data['status'] = 'OK'
with open(path, 'w') as handle:
    json.dump(data, handle, indent=2, sort_keys=True)
    handle.write('\n')
PY
