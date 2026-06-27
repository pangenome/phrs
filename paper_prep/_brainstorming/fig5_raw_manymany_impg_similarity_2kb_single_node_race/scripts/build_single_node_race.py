#!/usr/bin/env python3
"""Build and optionally submit Fig5 raw many:many IMPG 2 kb single-node race jobs."""

from __future__ import annotations

import argparse
import csv
import json
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_single_node_race"
SHARDED = ROOT / "paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded"
SHARDED_MANIFEST = SHARDED / "manifests/shard_manifest.tsv"
IMPG = Path("/home/erikg/.cargo/bin/impg")
COMPARISONS = (
    "PAN027mat_vs_PAN010_joint",
    "PAN027pat_vs_PAN011_joint",
    "PAN028mat_vs_PAN027_joint",
)
METHODS = (
    "sweepga_fastga_frequency32",
    "wfmash_p95_updated_bin",
)


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


def sharded_source_rows() -> list[dict[str, str]]:
    rows = read_tsv(SHARDED_MANIFEST)
    by_pair: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        method = row["method"]
        comparison_id = row["comparison_id"]
        if method not in METHODS or comparison_id not in COMPARISONS:
            continue
        key = (method, comparison_id)
        previous = by_pair.get(key)
        if previous is None or int(row["shard_index"]) < int(previous["shard_index"]):
            by_pair[key] = row
    expected = {(method, comparison) for method in METHODS for comparison in COMPARISONS}
    found = set(by_pair)
    if found != expected:
        raise SystemExit(f"Missing sharded manifest method/comparison rows: {sorted(expected - found)}")
    return [by_pair[(method, comparison)] for method in METHODS for comparison in COMPARISONS]


def write_job_script(row: dict[str, object], cpus: int, partition: str, time_limit: str) -> Path:
    method = str(row["method"])
    comparison_id = str(row["comparison_id"])
    script = OUT / "jobs" / f"{method}.{comparison_id}.full_bed.single_node.slurm.sh"
    script.parent.mkdir(parents=True, exist_ok=True)
    output_tsv_gz = Path(str(row["output_tsv_gz"]))
    metadata_json = Path(str(row["metadata_json"]))
    script.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                f"#SBATCH --job-name=fig5_1node_{method[:6]}_{comparison_id[:14]}",
                f"#SBATCH --partition={partition}",
                f"#SBATCH --cpus-per-task={cpus}",
                "#SBATCH --nodes=1",
                f"#SBATCH --time={time_limit}",
                f"#SBATCH --output={OUT / 'logs' / (method + '.' + comparison_id + '.%j.out')}",
                f"#SBATCH --error={OUT / 'logs' / (method + '.' + comparison_id + '.%j.err')}",
                "set -euo pipefail",
                "export LC_ALL=C",
                f"method={shlex.quote(method)}",
                f"comparison_id={shlex.quote(comparison_id)}",
                f"source_raw_paf={shlex.quote(str(row['source_raw_paf']))}",
                f"impg_alignment_paf={shlex.quote(str(row['impg_alignment_paf']))}",
                f"query_fasta={shlex.quote(str(row['query_fasta']))}",
                f"target_fasta={shlex.quote(str(row['target_fasta']))}",
                f"full_target_bed={shlex.quote(str(row['full_target_bed']))}",
                f"output_tsv_gz={shlex.quote(str(output_tsv_gz))}",
                f"metadata_json={shlex.quote(str(metadata_json))}",
                f"impg={shlex.quote(str(IMPG))}",
                "mkdir -p \"$(dirname \"$output_tsv_gz\")\" \"$(dirname \"$metadata_json\")\"",
                'tmp_tsv="${output_tsv_gz%.gz}.tmp.${SLURM_JOB_ID:-manual}"',
                'start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)',
                'node="${SLURMD_NODENAME:-$(hostname)}"',
                'partition="${SLURM_JOB_PARTITION:-unknown}"',
                'impg_version=$("$impg" --version)',
                'command_text="$impg similarity --alignment-files $impg_alignment_paf --target-bed $full_target_bed --sequence-files $query_fasta $target_fasta --gfa-engine poa --no-merge --num-mappings many:many --scaffold-jump 0 --threads ${SLURM_CPUS_PER_TASK}"',
                'python3 - "$metadata_json" "$method" "$comparison_id" "$source_raw_paf" "$impg_alignment_paf" "$query_fasta" "$target_fasta" "$full_target_bed" "$output_tsv_gz" "$impg" "$impg_version" "$command_text" "$start_utc" "$node" "$partition" <<\'PY\'',
                "import json, os, sys",
                "keys = ['metadata_json','method','comparison_id','source_raw_paf','impg_alignment_paf','query_fasta','target_fasta','full_target_bed','output_tsv_gz','impg_path','impg_version','command','start_utc','node','partition']",
                "data = dict(zip(keys, sys.argv[1:]))",
                "data.update({",
                "    'slurm_job_id': os.environ.get('SLURM_JOB_ID', 'manual'),",
                "    'slurm_cpus_per_task': os.environ.get('SLURM_CPUS_PER_TASK', 'unset'),",
                "    'status': 'RUNNING',",
                "})",
                "with open(data.pop('metadata_json'), 'w') as handle:",
                "    json.dump(data, handle, indent=2, sort_keys=True)",
                "    handle.write('\\n')",
                "PY",
                "echo \"$command_text\"",
                '"$impg" similarity \\',
                '  --alignment-files "$impg_alignment_paf" \\',
                '  --target-bed "$full_target_bed" \\',
                '  --sequence-files "$query_fasta" "$target_fasta" \\',
                "  --gfa-engine poa \\",
                "  --no-merge \\",
                "  --num-mappings many:many \\",
                "  --scaffold-jump 0 \\",
                '  --threads "${SLURM_CPUS_PER_TASK}" > "$tmp_tsv"',
                'gzip -f "$tmp_tsv"',
                "mv \"${tmp_tsv}.gz\" \"$output_tsv_gz\"",
                'finish_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)',
                'python3 - "$metadata_json" "$finish_utc" <<\'PY\'',
                "import json, sys",
                "path, finish = sys.argv[1:]",
                "with open(path) as handle:",
                "    data = json.load(handle)",
                "data['finish_utc'] = finish",
                "data['status'] = 'OK'",
                "with open(path, 'w') as handle:",
                "    json.dump(data, handle, indent=2, sort_keys=True)",
                "    handle.write('\\n')",
                "PY",
            ]
        )
        + "\n"
    )
    script.chmod(0o755)
    return script


def write_finalizer_script(partition: str) -> Path:
    script = OUT / "jobs/finalize_after_single_node_jobs.slurm.sh"
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "#SBATCH --job-name=fig5_1node_finalize",
                f"#SBATCH --partition={partition}",
                "#SBATCH --cpus-per-task=4",
                "#SBATCH --nodes=1",
                "#SBATCH --time=04:00:00",
                f"#SBATCH --output={OUT / 'logs/finalize_after_single_node.%j.out'}",
                f"#SBATCH --error={OUT / 'logs/finalize_after_single_node.%j.err'}",
                "set -euo pipefail",
                "export LC_ALL=C",
                f"python3 {shlex.quote(str(OUT / 'scripts/finalize_single_node_race.py'))}",
            ]
        )
        + "\n"
    )
    script.chmod(0o755)
    return script


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--partition", default="tux")
    parser.add_argument("--cpus", type=int, default=96)
    parser.add_argument("--time", default="7-00:00:00")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    for subdir in ("jobs", "logs", "manifests", "metadata/jobs", "outputs/all_hits", "summaries"):
        (OUT / subdir).mkdir(parents=True, exist_ok=True)

    impg_version = subprocess.run([str(IMPG), "--version"], check=True, text=True, stdout=subprocess.PIPE).stdout.strip()
    help_text = subprocess.run([str(IMPG), "similarity", "--help"], check=True, text=True, stdout=subprocess.PIPE).stdout
    (OUT / "impg_similarity_help.txt").write_text(help_text)

    rows: list[dict[str, object]] = []
    for source in sharded_source_rows():
        method = source["method"]
        comparison_id = source["comparison_id"]
        output_tsv_gz = OUT / "outputs/all_hits" / method / f"{method}.{comparison_id}.full_genome_2kb.impg_similarity.tsv.gz"
        metadata_json = OUT / "metadata/jobs" / method / f"{method}.{comparison_id}.full_genome_2kb.metadata.json"
        full_target_bed = Path(source["full_target_bed"])
        for required in (
            Path(source["source_raw_paf"]),
            Path(source["impg_alignment_paf"]),
            Path(source["query_fasta"]),
            Path(source["target_fasta"]),
            full_target_bed,
        ):
            if not required.exists():
                raise FileNotFoundError(required)
        command = " ".join(
            [
                shlex.quote(str(IMPG)),
                "similarity",
                "--alignment-files",
                shlex.quote(source["impg_alignment_paf"]),
                "--target-bed",
                shlex.quote(str(full_target_bed)),
                "--sequence-files",
                shlex.quote(source["query_fasta"]),
                shlex.quote(source["target_fasta"]),
                "--gfa-engine poa --no-merge --num-mappings many:many --scaffold-jump 0 --threads ${SLURM_CPUS_PER_TASK}",
            ]
        )
        row = {
            "method": method,
            "comparison_id": comparison_id,
            "source_raw_paf": source["source_raw_paf"],
            "impg_alignment_paf": source["impg_alignment_paf"],
            "query_fasta": source["query_fasta"],
            "target_fasta": source["target_fasta"],
            "full_target_bed": str(full_target_bed),
            "command": command,
            "job_script": "",
            "partition": args.partition,
            "cpus_per_task": args.cpus,
            "submitted_job_id": "",
            "output_tsv_gz": str(output_tsv_gz),
            "metadata_json": str(metadata_json),
            "raw_source_manifest": source["raw_source_manifest"],
            "source_raw_paf_size_bytes": source["source_raw_paf_size_bytes"],
            "impg_path": str(IMPG),
            "impg_version": impg_version,
            "state": "PENDING_SUBMISSION",
        }
        row["job_script"] = str(write_job_script(row, args.cpus, args.partition, args.time))
        rows.append(row)

    submitted_ids: list[str] = []
    if args.submit:
        for row in rows:
            proc = subprocess.run(["sbatch", "--parsable", str(row["job_script"])], check=True, text=True, stdout=subprocess.PIPE)
            job_id = proc.stdout.strip().split(";", 1)[0]
            row["submitted_job_id"] = job_id
            row["state"] = "SUBMITTED"
            submitted_ids.append(job_id)

    finalizer_script = write_finalizer_script(args.partition)
    finalizer_job_id = ""
    finalizer_dependency = ""
    if args.submit:
        finalizer_dependency = "afterany:" + ":".join(submitted_ids)
        proc = subprocess.run(
            ["sbatch", "--parsable", f"--dependency={finalizer_dependency}", str(finalizer_script)],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        )
        finalizer_job_id = proc.stdout.strip().split(";", 1)[0]

    job_fields = [
        "method",
        "comparison_id",
        "source_raw_paf",
        "impg_alignment_paf",
        "query_fasta",
        "target_fasta",
        "full_target_bed",
        "command",
        "job_script",
        "partition",
        "cpus_per_task",
        "submitted_job_id",
        "output_tsv_gz",
        "metadata_json",
        "raw_source_manifest",
        "source_raw_paf_size_bytes",
        "impg_path",
        "impg_version",
        "state",
    ]
    write_tsv(OUT / "manifests/single_node_job_manifest.tsv", rows, job_fields)
    write_tsv(
        OUT / "manifests/slurm_submission_manifest.tsv",
        [
            {
                "submitted_job_count": len(submitted_ids),
                "submitted_job_ids": ",".join(submitted_ids),
                "finalizer_script": str(finalizer_script),
                "finalizer_dependency": finalizer_dependency,
                "submitted_finalizer_job_id": finalizer_job_id,
                "partition": args.partition,
                "cpus_per_task": args.cpus,
                "state": "SUBMITTED" if args.submit else "PENDING_SUBMISSION",
            }
        ],
        [
            "submitted_job_count",
            "submitted_job_ids",
            "finalizer_script",
            "finalizer_dependency",
            "submitted_finalizer_job_id",
            "partition",
            "cpus_per_task",
            "state",
        ],
    )
    metadata = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "scope": "raw unfiltered many:many PAF-backed IMPG similarity over existing full-genome 2 kb target BEDs",
        "single_node_non_array_jobs": True,
        "no_new_alignment_or_graph_construction": True,
        "source_sharded_manifest": str(SHARDED_MANIFEST),
        "impg_path": str(IMPG),
        "impg_version": impg_version,
        "partition": args.partition,
        "cpus_per_task": args.cpus,
        "methods": list(METHODS),
        "comparisons": list(COMPARISONS),
        "finalizer_is_slurm_dependency_job": bool(finalizer_dependency),
        "finalizer_dependency": finalizer_dependency,
        "finalizer_job_id": finalizer_job_id,
    }
    (OUT / "pipeline_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(f"Generated {len(rows)} single-node full-BED IMPG jobs in {OUT}")
    if args.submit:
        print(f"Submitted jobs: {','.join(submitted_ids)}")
        print(f"Submitted finalizer: {finalizer_job_id} ({finalizer_dependency})")


if __name__ == "__main__":
    main()
