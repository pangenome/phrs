#!/usr/bin/env bash
#SBATCH --job-name=fig5_topn_sweepga_fastga_frequ
#SBATCH --partition=workers
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --array=0-151%1
#SBATCH --output=/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_mappingfiltered_impg_similarity_2kb_topn_depthcapped_10to10/logs/sweepga_fastga_frequency32.PAN028mat_vs_PAN027_joint.topn_depthcapped.array.slurm.%a.%A.out
#SBATCH --error=/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_mappingfiltered_impg_similarity_2kb_topn_depthcapped_10to10/logs/sweepga_fastga_frequency32.PAN028mat_vs_PAN027_joint.topn_depthcapped.array.slurm.%a.%A.err
set -euo pipefail
export LC_ALL=C
task_manifest=/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_mappingfiltered_impg_similarity_2kb_topn_depthcapped_10to10/manifests/array_tasks/sweepga_fastga_frequency32.PAN028mat_vs_PAN027_joint.tasks.tsv
impg=/home/erikg/.cargo/bin/impg
filter=/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_mappingfiltered_impg_similarity_2kb_topn_depthcapped_10to10/scripts/filter_impg_similarity_topn.py
pigz=/usr/bin/pigz
row=$(awk -F '\t' -v idx="${SLURM_ARRAY_TASK_ID}" 'NR > 1 && $4 == idx {print; exit}' "$task_manifest")
if [ -z "$row" ]; then echo "No task row for shard ${SLURM_ARRAY_TASK_ID}" >&2; exit 2; fi
IFS=$'\t' read -r method comparison_id source_raw_paf shard_index impg_alignment_paf query_fasta target_fasta full_target_bed filtered_bed output_tsv_gz skip_report metadata_json command_text <<< "$row"
mkdir -p "$(dirname "$output_tsv_gz")" "$(dirname "$skip_report")" "$(dirname "$metadata_json")"
tmp_gz="${output_tsv_gz}.tmp.${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}"
tmp_skip="${skip_report}.tmp.${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}"
start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
node="${SLURMD_NODENAME:-$(hostname)}"
partition="${SLURM_JOB_PARTITION:-unknown}"
impg_version=$("$impg" --version)
python3 - "$metadata_json" "$method" "$comparison_id" "$source_raw_paf" "$impg_alignment_paf" "$query_fasta" "$target_fasta" "$filtered_bed" "$output_tsv_gz" "$skip_report" "$impg" "$impg_version" "$command_text" "$start_utc" "$node" "$partition" <<'PY'
import json, os, sys
keys = ['metadata_json','method','comparison_id','source_raw_paf','impg_alignment_paf','query_fasta','target_fasta','filtered_bed','output_tsv_gz','skip_report','impg_path','impg_version','command','start_utc','node','partition']
data = dict(zip(keys, sys.argv[1:]))
data.update({'slurm_job_id': os.environ.get('SLURM_ARRAY_JOB_ID', os.environ.get('SLURM_JOB_ID', 'manual')), 'slurm_array_task_id': os.environ.get('SLURM_ARRAY_TASK_ID', '0'), 'slurm_cpus_per_task': os.environ.get('SLURM_CPUS_PER_TASK', 'unset'), 'status': 'RUNNING'})
with open(data.pop('metadata_json'), 'w') as handle:
    json.dump(data, handle, indent=2, sort_keys=True); handle.write('\n')
PY
if [ ! -s "$filtered_bed" ]; then
  printf 'chrom\tstart\tend\tgroup.a\tgroup.b\tgroup.a.length\tgroup.b.length\tintersection\tjaccard.similarity\tcosine.similarity\tdice.similarity\testimated.identity\n' | "$pigz" -p "${SLURM_CPUS_PER_TASK}" > "$tmp_gz"
  printf 'chrom\tstart\tend\traw_candidate_count\tretained_count\treason\n' > "$tmp_skip"
else
  echo "$command_text"
  "$impg" similarity --alignment-files "$impg_alignment_paf" --target-bed "$filtered_bed" --sequence-files "$query_fasta" "$target_fasta" --gfa-engine poa --no-merge --num-mappings 10:10 --scaffold-jump 0 --threads "${SLURM_CPUS_PER_TASK}" | python3 "$filter" --top-n 20 --max-candidates 500 --interchrom-only --skip-report "$tmp_skip" | "$pigz" -p "${SLURM_CPUS_PER_TASK}" > "$tmp_gz"
fi
mv "$tmp_gz" "$output_tsv_gz"
mv "$tmp_skip" "$skip_report"
finish_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
python3 - "$metadata_json" "$finish_utc" <<'PY'
import json, sys
path, finish = sys.argv[1:]
with open(path) as handle: data = json.load(handle)
data['finish_utc'] = finish; data['status'] = 'OK'
with open(path, 'w') as handle: json.dump(data, handle, indent=2, sort_keys=True); handle.write('\n')
PY
