#!/usr/bin/env python3
"""Build and optionally submit top-N/depth-capped Fig5 IMPG shard jobs."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


OUT = Path(__file__).resolve().parents[1]
REPO = Path("/moosefs/erikg/phrs")
SOURCE = REPO / ".wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded"
CENTROMERES = REPO / "data/chm13-annotations.bed"
IMPG = Path("/home/erikg/.cargo/bin/impg")
PIGZ = Path("/usr/bin/pigz")
FILTER = OUT / "scripts/filter_impg_similarity_topn.py"
CHR_RE = re.compile(r"(?:^|[_#])chr([0-9XY]+)$")


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def seq_chrom(seq: str) -> str:
    if "_chr" in seq:
        return "chr" + seq.rsplit("_chr", 1)[1]
    return seq


def load_centromeres(pad: int) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with CENTROMERES.open() as handle:
        for line in handle:
            chrom, start, end, name, *_ = line.rstrip("\n").split("\t")
            if name != "Centromere":
                continue
            out[chrom].append((max(0, int(start) - pad), int(end) + pad))
    return out


def overlaps(intervals: list[tuple[int, int]], start: int, end: int) -> bool:
    return any(start < stop and begin < end for begin, stop in intervals)


def filter_bed(source: Path, dest: Path, centromeres: dict[str, list[tuple[int, int]]]) -> dict[str, object]:
    total = kept = excluded = 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    with source.open() as inp, dest.open("w") as out:
        for line in inp:
            if not line.strip():
                continue
            total += 1
            seq, start_text, end_text, *rest = line.rstrip("\n").split("\t")
            start = int(start_text)
            end = int(end_text)
            chrom = seq_chrom(seq)
            if overlaps(centromeres.get(chrom, []), start, end):
                excluded += 1
                continue
            out.write("\t".join([seq, start_text, end_text, *rest]) + "\n")
            kept += 1
    return {"source_bed": str(source), "filtered_bed": str(dest), "total_windows": total, "kept_windows": kept, "excluded_centromere_windows": excluded}


def job_script(
    path: Path,
    task_manifest: Path,
    top_n: int,
    max_candidates: int,
    num_mappings: str,
    time_limit: str,
    cpus: int,
    max_parallel: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name=fig5_topn_{path.name[:20]}",
        "#SBATCH --partition=workers",
        f"#SBATCH --cpus-per-task={cpus}",
        "#SBATCH --nodes=1",
        f"#SBATCH --time={time_limit}",
        f"#SBATCH --array=0-PLACEHOLDER%{max_parallel}",
        f"#SBATCH --output={OUT}/logs/{path.stem}.%a.%A.out",
        f"#SBATCH --error={OUT}/logs/{path.stem}.%a.%A.err",
        "set -euo pipefail",
        "export LC_ALL=C",
        f"task_manifest={shlex.quote(str(task_manifest))}",
        f"impg={shlex.quote(str(IMPG))}",
        f"filter={shlex.quote(str(FILTER))}",
        f"pigz={shlex.quote(str(PIGZ))}",
        "row=$(awk -F '\\t' -v idx=\"${SLURM_ARRAY_TASK_ID}\" 'NR > 1 && $4 == idx {print; exit}' \"$task_manifest\")",
        "if [ -z \"$row\" ]; then echo \"No task row for shard ${SLURM_ARRAY_TASK_ID}\" >&2; exit 2; fi",
        "IFS=$'\\t' read -r method comparison_id source_raw_paf shard_index impg_alignment_paf query_fasta target_fasta full_target_bed filtered_bed output_tsv_gz skip_report metadata_json command_text <<< \"$row\"",
        "mkdir -p \"$(dirname \"$output_tsv_gz\")\" \"$(dirname \"$skip_report\")\" \"$(dirname \"$metadata_json\")\"",
        "tmp_gz=\"${output_tsv_gz}.tmp.${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}\"",
        "tmp_skip=\"${skip_report}.tmp.${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}\"",
        "start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)",
        "node=\"${SLURMD_NODENAME:-$(hostname)}\"",
        "partition=\"${SLURM_JOB_PARTITION:-unknown}\"",
        "impg_version=$(\"$impg\" --version)",
        "python3 - \"$metadata_json\" \"$method\" \"$comparison_id\" \"$source_raw_paf\" \"$impg_alignment_paf\" \"$query_fasta\" \"$target_fasta\" \"$filtered_bed\" \"$output_tsv_gz\" \"$skip_report\" \"$impg\" \"$impg_version\" \"$command_text\" \"$start_utc\" \"$node\" \"$partition\" <<'PY'",
        "import json, os, sys",
        "keys = ['metadata_json','method','comparison_id','source_raw_paf','impg_alignment_paf','query_fasta','target_fasta','filtered_bed','output_tsv_gz','skip_report','impg_path','impg_version','command','start_utc','node','partition']",
        "data = dict(zip(keys, sys.argv[1:]))",
        "data.update({'slurm_job_id': os.environ.get('SLURM_ARRAY_JOB_ID', os.environ.get('SLURM_JOB_ID', 'manual')), 'slurm_array_task_id': os.environ.get('SLURM_ARRAY_TASK_ID', '0'), 'slurm_cpus_per_task': os.environ.get('SLURM_CPUS_PER_TASK', 'unset'), 'status': 'RUNNING'})",
        "with open(data.pop('metadata_json'), 'w') as handle:",
        "    json.dump(data, handle, indent=2, sort_keys=True); handle.write('\\n')",
        "PY",
        "if [ ! -s \"$filtered_bed\" ]; then",
        "  printf 'chrom\\tstart\\tend\\tgroup.a\\tgroup.b\\tgroup.a.length\\tgroup.b.length\\tintersection\\tjaccard.similarity\\tcosine.similarity\\tdice.similarity\\testimated.identity\\n' | \"$pigz\" -p \"${SLURM_CPUS_PER_TASK}\" > \"$tmp_gz\"",
        "  printf 'chrom\\tstart\\tend\\traw_candidate_count\\tretained_count\\treason\\n' > \"$tmp_skip\"",
        "else",
        "  echo \"$command_text\"",
        "  \"$impg\" similarity --alignment-files \"$impg_alignment_paf\" --target-bed \"$filtered_bed\" --sequence-files \"$query_fasta\" \"$target_fasta\" --gfa-engine poa --no-merge --num-mappings " + shlex.quote(num_mappings) + " --scaffold-jump 0 --threads \"${SLURM_CPUS_PER_TASK}\" | python3 \"$filter\" --top-n " + str(top_n) + " --max-candidates " + str(max_candidates) + " --interchrom-only --skip-report \"$tmp_skip\" | \"$pigz\" -p \"${SLURM_CPUS_PER_TASK}\" > \"$tmp_gz\"",
        "fi",
        "mv \"$tmp_gz\" \"$output_tsv_gz\"",
        "mv \"$tmp_skip\" \"$skip_report\"",
        "finish_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)",
        "python3 - \"$metadata_json\" \"$finish_utc\" <<'PY'",
        "import json, sys",
        "path, finish = sys.argv[1:]",
        "with open(path) as handle: data = json.load(handle)",
        "data['finish_utc'] = finish; data['status'] = 'OK'",
        "with open(path, 'w') as handle: json.dump(data, handle, indent=2, sort_keys=True); handle.write('\\n')",
        "PY",
    ]
    text = "\n".join(lines) + "\n"
    shard_count = sum(1 for _ in task_manifest.open()) - 1
    text = text.replace("0-PLACEHOLDER", f"0-{shard_count - 1}")
    path.write_text(text)
    path.chmod(0o755)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--max-candidates", type=int, default=500)
    parser.add_argument("--num-mappings", default="many:many")
    parser.add_argument("--centromere-pad", type=int, default=0)
    parser.add_argument("--cpus", type=int, default=48)
    parser.add_argument("--max-parallel", type=int, default=6)
    parser.add_argument("--time", default="7-00:00:00")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args()

    centromeres = load_centromeres(args.centromere_pad)
    rows = read_tsv(SOURCE / "manifests/shard_manifest.tsv")
    bed_reports = []
    new_rows = []
    for row in rows:
        filtered_bed = OUT / "bed_shards_no_centromere" / row["comparison_id"] / Path(row["bed_shard"]).name
        if not filtered_bed.exists():
            bed_reports.append(filter_bed(Path(row["bed_shard"]), filtered_bed, centromeres))
        output = OUT / "outputs/shards" / row["method"] / row["comparison_id"] / f"{row['method']}.{row['comparison_id']}.full_genome_2kb.shard_{int(row['shard_index']):04d}.topn{args.top_n}.impg_similarity.tsv.gz"
        skip = OUT / "outputs/skipped_windows" / row["method"] / row["comparison_id"] / f"{row['method']}.{row['comparison_id']}.full_genome_2kb.shard_{int(row['shard_index']):04d}.skipped.tsv"
        meta = OUT / "metadata/shards" / row["method"] / row["comparison_id"] / f"{row['method']}.{row['comparison_id']}.full_genome_2kb.shard_{int(row['shard_index']):04d}.metadata.json"
        command = f"{IMPG} similarity --alignment-files {row['impg_alignment_paf']} --target-bed {filtered_bed} --sequence-files {row['query_fasta']} {row['target_fasta']} --gfa-engine poa --no-merge --num-mappings {args.num_mappings} --scaffold-jump 0 --threads ${{SLURM_CPUS_PER_TASK}} | filter top {args.top_n}, skip >{args.max_candidates}"
        new_rows.append({**row, "filtered_bed": str(filtered_bed), "output_tsv_gz": str(output), "skip_report": str(skip), "metadata_json": str(meta), "command": command})

    fields = ["method", "comparison_id", "source_raw_paf", "shard_index", "impg_alignment_paf", "query_fasta", "target_fasta", "full_target_bed", "filtered_bed", "output_tsv_gz", "skip_report", "metadata_json", "command"]
    write_tsv(OUT / "manifests/shard_manifest.tsv", new_rows, fields)
    write_tsv(OUT / "manifests/centromere_filter_report.tsv", bed_reports, ["source_bed", "filtered_bed", "total_windows", "kept_windows", "excluded_centromere_windows"])

    submitted = []
    for key, pair_rows in sorted(defaultdict(list, {k: [r for r in new_rows if (r["method"], r["comparison_id"]) == k] for k in {(r["method"], r["comparison_id"]) for r in new_rows}}).items()):
        method, comparison_id = key
        task_manifest = OUT / "manifests/array_tasks" / f"{method}.{comparison_id}.tasks.tsv"
        write_tsv(task_manifest, sorted(pair_rows, key=lambda r: int(r["shard_index"])), fields)
        script = OUT / "jobs" / f"{method}.{comparison_id}.topn_depthcapped.array.slurm.sh"
        job_script(script, task_manifest, args.top_n, args.max_candidates, args.num_mappings, args.time, args.cpus, args.max_parallel)
        job_id = ""
        state = "CREATED"
        if args.submit:
            proc = subprocess.run(["sbatch", "--parsable", str(script)], check=True, text=True, stdout=subprocess.PIPE)
            job_id = proc.stdout.strip()
            state = "SUBMITTED"
        submitted.append({"method": method, "comparison_id": comparison_id, "job_script": str(script), "task_manifest": str(task_manifest), "shard_count": len(pair_rows), "array_range": f"0-{len(pair_rows)-1}%{args.max_parallel}", "partition": "workers", "cpus_per_task": args.cpus, "submitted_array_job_id": job_id, "state": state})
    write_tsv(OUT / "manifests/slurm_array_manifest.tsv", submitted, ["method", "comparison_id", "job_script", "task_manifest", "shard_count", "array_range", "partition", "cpus_per_task", "submitted_array_job_id", "state"])
    metadata = {"created_utc": utc_now(), "top_n": args.top_n, "max_candidates": args.max_candidates, "num_mappings": args.num_mappings, "centromere_pad": args.centromere_pad, "source": str(SOURCE), "submitted": args.submit}
    (OUT / "pipeline_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
