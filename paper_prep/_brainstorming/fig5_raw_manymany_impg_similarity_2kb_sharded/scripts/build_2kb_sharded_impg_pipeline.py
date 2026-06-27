#!/usr/bin/env python3
"""Build and optionally submit Fig5 raw many:many IMPG 2 kb shard jobs.

This pipeline consumes existing raw PAFs and query/target FASTAs only.  It does
not run WFMASH, SweepGA, FastGA, minimap2, seqwish, odgi, or graph construction.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded"
WFMASH_MANIFEST = ROOT / (
    "paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/"
    "summaries/query_grid_filter_manifest.tsv"
)
SWEEPGA_MANIFEST = ROOT / (
    "paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency32/"
    "summaries/query_grid_chop_filter_manifest.tsv"
)
INPUT_MANIFEST = ROOT / (
    "paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/"
    "summaries/input_manifest.tsv"
)
PREVIOUS_BGZF_MANIFEST = ROOT.parent / "agent-2762/paper_prep/_brainstorming/" / (
    "fig5_raw_manymany_impg_similarity_scan/bgzf_manifest.tsv"
)
IMPG = Path("/home/erikg/.cargo/bin/impg")
BGZIP = Path("/home/erikg/.guix-profile/bin/bgzip")

COMPARISONS = (
    "PAN027mat_vs_PAN010_joint",
    "PAN027pat_vs_PAN011_joint",
    "PAN028mat_vs_PAN027_joint",
)
METHODS = (
    ("wfmash_p95_updated_bin", WFMASH_MANIFEST),
    ("sweepga_fastga_frequency32", SWEEPGA_MANIFEST),
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


def bgzip_valid(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        subprocess.run([str(BGZIP), "-t", str(path)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        return False
    return True


def raw_rows_for_method(method: str, manifest: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for row in read_tsv(manifest):
        comparison_id = row["comparison_id"]
        if comparison_id not in COMPARISONS:
            continue
        if row.get("status") != "OK":
            continue
        if int(row.get("chop_length_bp", "0")) != 2000:
            continue
        key = (comparison_id, row["raw_paf"])
        if key in seen:
            continue
        seen.add(key)
        rows.append(
            {
                "method": method,
                "comparison_id": comparison_id,
                "source_manifest": str(manifest),
                "raw_paf": row["raw_paf"],
                "raw_row_count": row.get("raw_row_count", ""),
            }
        )
    return sorted(rows, key=lambda item: item["comparison_id"])


def previous_bgzf_by_source() -> dict[str, str]:
    out: dict[str, str] = {}
    if not PREVIOUS_BGZF_MANIFEST.exists():
        return out
    for row in read_tsv(PREVIOUS_BGZF_MANIFEST):
        source = row.get("source_raw_paf", "")
        bgzf = Path(row.get("bgzf_paf", ""))
        if source and bgzip_valid(bgzf):
            out[source] = str(bgzf)
    return out


def write_full_bed_and_shards(
    comparison_id: str,
    target_fasta: Path,
    window_size: int,
    shard_windows: int,
) -> tuple[Path, list[dict[str, object]]]:
    fai = Path(f"{target_fasta}.fai")
    if not fai.exists():
        raise FileNotFoundError(f"Missing target FASTA index: {fai}")
    bed_dir = OUT / "beds"
    shard_dir = OUT / "bed_shards" / comparison_id
    bed_dir.mkdir(parents=True, exist_ok=True)
    shard_dir.mkdir(parents=True, exist_ok=True)
    full_bed = bed_dir / f"{comparison_id}.full_genome_2kb.bed"

    shard_rows: list[dict[str, object]] = []
    shard_index = 0
    shard_line_count = 0
    shard_start_window = 0
    window_index = 0
    current_shard: Path | None = None
    shard_handle = None

    def open_shard() -> None:
        nonlocal current_shard, shard_handle, shard_line_count, shard_start_window
        shard_start_window = window_index
        shard_line_count = 0
        current_shard = shard_dir / f"{comparison_id}.full_genome_2kb.shard_{shard_index:04d}.bed"
        shard_handle = current_shard.open("w")

    def close_shard() -> None:
        nonlocal shard_handle
        if shard_handle is None or current_shard is None:
            return
        shard_handle.close()
        shard_rows.append(
            {
                "comparison_id": comparison_id,
                "target_fasta": str(target_fasta),
                "target_fai": str(fai),
                "full_bed": str(full_bed),
                "shard_index": shard_index,
                "shard_bed": str(current_shard),
                "window_size_bp": window_size,
                "shard_window_count": shard_line_count,
                "first_window_index": shard_start_window,
                "last_window_index": window_index - 1,
            }
        )
        shard_handle = None

    with fai.open() as src, full_bed.open("w") as full:
        open_shard()
        for line in src:
            seq, length_text, *_ = line.rstrip("\n").split("\t")
            length = int(length_text)
            for start in range(0, length, window_size):
                end = min(start + window_size, length)
                record = f"{seq}\t{start}\t{end}\n"
                full.write(record)
                assert shard_handle is not None
                shard_handle.write(record)
                shard_line_count += 1
                window_index += 1
                if shard_line_count == shard_windows:
                    close_shard()
                    shard_index += 1
                    open_shard()
        close_shard()
    return full_bed, shard_rows


def write_bgzf_job(row: dict[str, object], cpus: int, partition: str) -> Path:
    script = Path(str(row["job_script"]))
    script.parent.mkdir(parents=True, exist_ok=True)
    Path(str(row["bgzf_paf"])).parent.mkdir(parents=True, exist_ok=True)
    script.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                f"#SBATCH --job-name=fig5_bgzf_{str(row['comparison_id'])[:14]}",
                f"#SBATCH --partition={partition}",
                f"#SBATCH --cpus-per-task={cpus}",
                "#SBATCH --nodes=1",
                "#SBATCH --time=12:00:00",
                f"#SBATCH --output={row['stdout_log']}",
                f"#SBATCH --error={row['stderr_log']}",
                "set -euo pipefail",
                "export LC_ALL=C",
                f"source_paf={shlex.quote(str(row['source_raw_paf']))}",
                f"dest_paf={shlex.quote(str(row['bgzf_paf']))}",
                'tmp_paf="${dest_paf}.tmp.${SLURM_JOB_ID:-manual}"',
                f"mkdir -p {shlex.quote(str(Path(str(row['bgzf_paf'])).parent))}",
                "if [ -s \"$dest_paf\" ]; then",
                f"  {shlex.quote(str(BGZIP))} -t \"$dest_paf\"",
                "  exit 0",
                "fi",
                f"gzip -dc \"$source_paf\" | {shlex.quote(str(BGZIP))} -@ ${{SLURM_CPUS_PER_TASK}} -c > \"$tmp_paf\"",
                f"{shlex.quote(str(BGZIP))} -t \"$tmp_paf\"",
                "mv \"$tmp_paf\" \"$dest_paf\"",
            ]
        )
        + "\n"
    )
    script.chmod(0o755)
    return script


def write_array_script(method: str, comparison_id: str, task_manifest: Path, task_count: int, cpus: int, partition: str, max_parallel: int) -> Path:
    script = OUT / "jobs" / f"{method}.{comparison_id}.impg_2kb_shards.array.slurm.sh"
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                f"#SBATCH --job-name=fig5_impg_{method[:6]}_{comparison_id[:14]}",
                f"#SBATCH --partition={partition}",
                f"#SBATCH --cpus-per-task={cpus}",
                "#SBATCH --nodes=1",
                "#SBATCH --time=24:00:00",
                f"#SBATCH --array=0-{task_count - 1}%{max_parallel}",
                f"#SBATCH --output={OUT / 'logs' / (method + '.' + comparison_id + '.shard_%a.%A.out')}",
                f"#SBATCH --error={OUT / 'logs' / (method + '.' + comparison_id + '.shard_%a.%A.err')}",
                "set -euo pipefail",
                "export LC_ALL=C",
                f"task_manifest={shlex.quote(str(task_manifest))}",
                f"impg={shlex.quote(str(IMPG))}",
                "row=$(awk -F '\\t' -v idx=\"${SLURM_ARRAY_TASK_ID}\" 'NR > 1 && $4 == idx {print; exit}' \"$task_manifest\")",
                "if [ -z \"$row\" ]; then",
                "  echo \"No task row for shard ${SLURM_ARRAY_TASK_ID}\" >&2",
                "  exit 2",
                "fi",
                "IFS=$'\\t' read -r method comparison_id source_raw_paf shard_index impg_alignment_paf query_fasta target_fasta shard_bed output_tsv_gz metadata_json command_text <<< \"$row\"",
                "mkdir -p \"$(dirname \"$output_tsv_gz\")\" \"$(dirname \"$metadata_json\")\"",
                'tmp_tsv="${output_tsv_gz%.gz}.tmp.${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}"',
                'start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)',
                'node="${SLURMD_NODENAME:-$(hostname)}"',
                'partition="${SLURM_JOB_PARTITION:-unknown}"',
                'impg_version=$("$impg" --version)',
                'python3 - "$metadata_json" "$method" "$comparison_id" "$source_raw_paf" "$impg_alignment_paf" "$query_fasta" "$target_fasta" "$shard_bed" "$output_tsv_gz" "$impg" "$impg_version" "$command_text" "$start_utc" "$node" "$partition" <<\'PY\'',
                "import json, sys",
                "keys = ['metadata_json','method','comparison_id','source_raw_paf','impg_alignment_paf','query_fasta','target_fasta','bed_shard','output_tsv_gz','impg_path','impg_version','command','start_utc','node','partition']",
                "data = dict(zip(keys, sys.argv[1:]))",
                "data.update({",
                "    'slurm_job_id': __import__('os').environ.get('SLURM_ARRAY_JOB_ID', __import__('os').environ.get('SLURM_JOB_ID', 'manual')),",
                "    'slurm_array_task_id': __import__('os').environ.get('SLURM_ARRAY_TASK_ID', '0'),",
                "    'slurm_cpus_per_task': __import__('os').environ.get('SLURM_CPUS_PER_TASK', 'unset'),",
                "    'status': 'RUNNING',",
                "})",
                "with open(data.pop('metadata_json'), 'w') as handle:",
                "    json.dump(data, handle, indent=2, sort_keys=True)",
                "    handle.write('\\n')",
                "PY",
                "echo \"$command_text\"",
                '"$impg" similarity \\',
                '  --alignment-files "$impg_alignment_paf" \\',
                '  --target-bed "$shard_bed" \\',
                '  --sequence-files "$query_fasta" "$target_fasta" \\',
                "  --gfa-engine poa \\",
                "  --no-merge \\",
                "  --num-mappings many:many \\",
                "  --scaffold-jump 0 \\",
                '  --threads "${SLURM_CPUS_PER_TASK}" > "$tmp_tsv"',
                'gzip -f "$tmp_tsv"',
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--window-size", type=int, default=2000)
    parser.add_argument("--shard-windows", type=int, default=20000)
    parser.add_argument("--cpus", type=int, default=48)
    parser.add_argument("--partition", default="workers")
    parser.add_argument("--max-parallel", type=int, default=18)
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "logs").mkdir(parents=True, exist_ok=True)
    inputs = {row["comparison_id"]: row for row in read_tsv(INPUT_MANIFEST) if row["comparison_id"] in COMPARISONS}
    if set(inputs) != set(COMPARISONS):
        missing = sorted(set(COMPARISONS) - set(inputs))
        raise SystemExit(f"Missing input manifest rows: {missing}")

    previous_bgzf = previous_bgzf_by_source()
    raw_rows: list[dict[str, str]] = []
    for method, manifest in METHODS:
        raw_rows.extend(raw_rows_for_method(method, manifest))
    if len(raw_rows) != 6:
        raise SystemExit(f"Expected six raw method/comparison rows, found {len(raw_rows)}")

    bed_rows_by_comparison: dict[str, list[dict[str, object]]] = {}
    full_beds: dict[str, Path] = {}
    for comparison_id in COMPARISONS:
        target_fasta = Path(inputs[comparison_id]["target_fasta"])
        full_bed, shard_rows = write_full_bed_and_shards(
            comparison_id, target_fasta, args.window_size, args.shard_windows
        )
        full_beds[comparison_id] = full_bed
        bed_rows_by_comparison[comparison_id] = shard_rows

    bgzf_rows: list[dict[str, object]] = []
    shard_rows: list[dict[str, object]] = []
    task_rows_by_pair: dict[tuple[str, str], list[dict[str, object]]] = {}
    impg_version = subprocess.run([str(IMPG), "--version"], check=True, text=True, stdout=subprocess.PIPE).stdout.strip()
    help_text = subprocess.run([str(IMPG), "similarity", "--help"], check=True, text=True, stdout=subprocess.PIPE).stdout
    (OUT / "impg_similarity_help.txt").write_text(help_text)

    for raw in raw_rows:
        method = raw["method"]
        comparison_id = raw["comparison_id"]
        raw_paf = Path(raw["raw_paf"])
        if not raw_paf.exists():
            raise FileNotFoundError(raw_paf)
        input_row = inputs[comparison_id]
        query_fasta = Path(input_row["query_fasta"])
        target_fasta = Path(input_row["target_fasta"])
        if not query_fasta.exists() or not target_fasta.exists():
            raise FileNotFoundError(f"Missing FASTA for {comparison_id}")
        impg_paf = raw_paf
        bgzf_status = "not_required"
        bgzf_source = ""
        if method == "sweepga_fastga_frequency32":
            if str(raw_paf) in previous_bgzf:
                impg_paf = Path(previous_bgzf[str(raw_paf)])
                bgzf_status = "reused_validated_previous_task_bgzf"
                bgzf_source = str(PREVIOUS_BGZF_MANIFEST)
            else:
                impg_paf = OUT / "bgzf_raw_paf" / raw_paf.name.replace(".paf.gz", ".bgzf.paf.gz")
                bgzf_status = "needs_bgzip_normalization"
            bgzf_row = {
                "comparison_id": comparison_id,
                "source_raw_paf": str(raw_paf),
                "bgzf_paf": str(impg_paf),
                "bgzf_status": bgzf_status,
                "bgzf_source_manifest": bgzf_source,
                "job_script": str(OUT / "jobs" / f"normalize_bgzf.{comparison_id}.slurm.sh"),
                "stdout_log": str(OUT / "logs" / f"normalize_bgzf.{comparison_id}.%j.out"),
                "stderr_log": str(OUT / "logs" / f"normalize_bgzf.{comparison_id}.%j.err"),
                "partition": args.partition,
                "cpus_per_task": args.cpus,
                "submitted_job_id": "",
            }
            if bgzf_status == "needs_bgzip_normalization":
                write_bgzf_job(bgzf_row, args.cpus, args.partition)
            bgzf_rows.append(bgzf_row)

        for bed_row in bed_rows_by_comparison[comparison_id]:
            shard_index = int(bed_row["shard_index"])
            stem = f"{method}.{comparison_id}.full_genome_2kb.shard_{shard_index:04d}"
            output_tsv_gz = OUT / "outputs" / "shards" / method / comparison_id / f"{stem}.impg_similarity.tsv.gz"
            metadata_json = OUT / "metadata" / "shards" / method / comparison_id / f"{stem}.metadata.json"
            command = " ".join(
                [
                    shlex.quote(str(IMPG)),
                    "similarity",
                    "--alignment-files",
                    shlex.quote(str(impg_paf)),
                    "--target-bed",
                    shlex.quote(str(bed_row["shard_bed"])),
                    "--sequence-files",
                    shlex.quote(str(query_fasta)),
                    shlex.quote(str(target_fasta)),
                    "--gfa-engine poa --no-merge --num-mappings many:many --scaffold-jump 0 --threads ${SLURM_CPUS_PER_TASK}",
                ]
            )
            row = {
                "method": method,
                "comparison_id": comparison_id,
                "source_raw_paf": str(raw_paf),
                "shard_index": shard_index,
                "impg_alignment_paf": str(impg_paf),
                "query_fasta": str(query_fasta),
                "target_fasta": str(target_fasta),
                "bed_shard": str(bed_row["shard_bed"]),
                "output_tsv_gz": str(output_tsv_gz),
                "metadata_json": str(metadata_json),
                "command": command,
                "full_target_bed": str(full_beds[comparison_id]),
                "raw_source_manifest": raw["source_manifest"],
                "source_raw_paf_size_bytes": raw_paf.stat().st_size,
                "impg_path": str(IMPG),
                "impg_version": impg_version,
                "slurm_array_job_id": "",
                "slurm_array_task_id": shard_index,
                "node": "",
                "partition": args.partition,
                "slurm_cpus_per_task": args.cpus,
                "state": "PENDING_SUBMISSION",
            }
            shard_rows.append(row)
            task_rows_by_pair.setdefault((method, comparison_id), []).append(row)

    array_rows: list[dict[str, object]] = []
    bgzf_job_by_comparison: dict[str, str] = {}
    if args.submit:
        for row in bgzf_rows:
            if row["bgzf_status"] == "needs_bgzip_normalization":
                proc = subprocess.run(["sbatch", "--parsable", str(row["job_script"])], check=True, text=True, stdout=subprocess.PIPE)
                row["submitted_job_id"] = proc.stdout.strip()
                bgzf_job_by_comparison[str(row["comparison_id"])] = str(row["submitted_job_id"])

    for (method, comparison_id), rows in sorted(task_rows_by_pair.items()):
        task_manifest = OUT / "manifests" / "array_tasks" / f"{method}.{comparison_id}.tasks.tsv"
        task_fields = [
            "method",
            "comparison_id",
            "source_raw_paf",
            "shard_index",
            "impg_alignment_paf",
            "query_fasta",
            "target_fasta",
            "bed_shard",
            "output_tsv_gz",
            "metadata_json",
            "command",
        ]
        write_tsv(task_manifest, rows, task_fields)
        script = write_array_script(method, comparison_id, task_manifest, len(rows), args.cpus, args.partition, args.max_parallel)
        dependency = ""
        submitted = ""
        if args.submit:
            cmd = ["sbatch", "--parsable"]
            if method == "sweepga_fastga_frequency32" and comparison_id in bgzf_job_by_comparison:
                dependency = f"afterok:{bgzf_job_by_comparison[comparison_id]}"
                cmd.append(f"--dependency={dependency}")
            cmd.append(str(script))
            proc = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE)
            submitted = proc.stdout.strip()
            for shard in rows:
                shard["slurm_array_job_id"] = submitted
                shard["state"] = "SUBMITTED"
        array_rows.append(
            {
                "method": method,
                "comparison_id": comparison_id,
                "array_script": str(script),
                "task_manifest": str(task_manifest),
                "shard_count": len(rows),
                "array_range": f"0-{len(rows) - 1}%{args.max_parallel}",
                "partition": args.partition,
                "cpus_per_task": args.cpus,
                "dependency": dependency,
                "submitted_array_job_id": submitted,
                "state": "SUBMITTED" if submitted else "PENDING_SUBMISSION",
            }
        )

    write_tsv(
        OUT / "manifests/target_bed_shards.tsv",
        [row for rows in bed_rows_by_comparison.values() for row in rows],
        [
            "comparison_id",
            "target_fasta",
            "target_fai",
            "full_bed",
            "shard_index",
            "shard_bed",
            "window_size_bp",
            "shard_window_count",
            "first_window_index",
            "last_window_index",
        ],
    )
    write_tsv(
        OUT / "manifests/raw_paf_bgzf_manifest.tsv",
        bgzf_rows,
        [
            "comparison_id",
            "source_raw_paf",
            "bgzf_paf",
            "bgzf_status",
            "bgzf_source_manifest",
            "job_script",
            "stdout_log",
            "stderr_log",
            "partition",
            "cpus_per_task",
            "submitted_job_id",
        ],
    )
    write_tsv(
        OUT / "manifests/shard_manifest.tsv",
        shard_rows,
        [
            "method",
            "comparison_id",
            "source_raw_paf",
            "impg_alignment_paf",
            "query_fasta",
            "target_fasta",
            "full_target_bed",
            "bed_shard",
            "shard_index",
            "command",
            "slurm_array_job_id",
            "slurm_array_task_id",
            "node",
            "partition",
            "slurm_cpus_per_task",
            "impg_path",
            "impg_version",
            "output_tsv_gz",
            "metadata_json",
            "raw_source_manifest",
            "source_raw_paf_size_bytes",
            "state",
        ],
    )
    write_tsv(
        OUT / "manifests/slurm_array_manifest.tsv",
        array_rows,
        [
            "method",
            "comparison_id",
            "array_script",
            "task_manifest",
            "shard_count",
            "array_range",
            "partition",
            "cpus_per_task",
            "dependency",
            "submitted_array_job_id",
            "state",
        ],
    )
    metadata = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "scope": "raw unfiltered many:many PAF-backed IMPG similarity over exact full-genome 2 kb target windows",
        "window_size_bp": args.window_size,
        "shard_windows": args.shard_windows,
        "impg_path": str(IMPG),
        "impg_version": impg_version,
        "no_new_alignment_or_graph_construction": True,
        "methods": [method for method, _ in METHODS],
        "comparisons": list(COMPARISONS),
    }
    (OUT / "pipeline_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(f"Generated {len(shard_rows)} IMPG shard tasks in {OUT}")


if __name__ == "__main__":
    main()
