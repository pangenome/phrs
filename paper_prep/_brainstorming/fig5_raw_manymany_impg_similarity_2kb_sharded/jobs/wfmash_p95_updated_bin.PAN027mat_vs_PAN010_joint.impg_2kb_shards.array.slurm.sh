#!/usr/bin/env bash
#SBATCH --job-name=fig5_impg_wfmash_PAN027mat_vs_P
#SBATCH --partition=workers
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --array=0-151%6
#SBATCH --output=/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/wfmash_p95_updated_bin.PAN027mat_vs_PAN010_joint.shard_%a.%A.out
#SBATCH --error=/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/wfmash_p95_updated_bin.PAN027mat_vs_PAN010_joint.shard_%a.%A.err
set -euo pipefail
export LC_ALL=C
task_manifest=/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/manifests/array_tasks/wfmash_p95_updated_bin.PAN027mat_vs_PAN010_joint.tasks.tsv
impg=/home/erikg/.cargo/bin/impg
row=$(awk -F '\t' -v idx="${SLURM_ARRAY_TASK_ID}" 'NR > 1 && $4 == idx {print; exit}' "$task_manifest")
if [ -z "$row" ]; then
  echo "No task row for shard ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi
IFS=$'\t' read -r method comparison_id source_raw_paf shard_index impg_alignment_paf query_fasta target_fasta shard_bed output_tsv_gz metadata_json command_text <<< "$row"
mkdir -p "$(dirname "$output_tsv_gz")" "$(dirname "$metadata_json")"
tmp_tsv="${output_tsv_gz%.gz}.tmp.${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}"
start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
node="${SLURMD_NODENAME:-$(hostname)}"
partition="${SLURM_JOB_PARTITION:-unknown}"
impg_version=$("$impg" --version)
python3 - "$metadata_json" "$method" "$comparison_id" "$source_raw_paf" "$impg_alignment_paf" "$query_fasta" "$target_fasta" "$shard_bed" "$output_tsv_gz" "$impg" "$impg_version" "$command_text" "$start_utc" "$node" "$partition" <<'PY'
import json, sys
keys = ['metadata_json','method','comparison_id','source_raw_paf','impg_alignment_paf','query_fasta','target_fasta','bed_shard','output_tsv_gz','impg_path','impg_version','command','start_utc','node','partition']
data = dict(zip(keys, sys.argv[1:]))
data.update({
    'slurm_job_id': __import__('os').environ.get('SLURM_ARRAY_JOB_ID', __import__('os').environ.get('SLURM_JOB_ID', 'manual')),
    'slurm_array_task_id': __import__('os').environ.get('SLURM_ARRAY_TASK_ID', '0'),
    'slurm_cpus_per_task': __import__('os').environ.get('SLURM_CPUS_PER_TASK', 'unset'),
    'status': 'RUNNING',
})
with open(data.pop('metadata_json'), 'w') as handle:
    json.dump(data, handle, indent=2, sort_keys=True)
    handle.write('\n')
PY
echo "$command_text"
"$impg" similarity \
  --alignment-files "$impg_alignment_paf" \
  --target-bed "$shard_bed" \
  --sequence-files "$query_fasta" "$target_fasta" \
  --gfa-engine poa \
  --no-merge \
  --num-mappings many:many \
  --scaffold-jump 0 \
  --threads "${SLURM_CPUS_PER_TASK}" > "$tmp_tsv"
gzip -f "$tmp_tsv"
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
