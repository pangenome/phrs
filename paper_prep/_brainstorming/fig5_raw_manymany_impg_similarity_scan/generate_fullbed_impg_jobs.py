#!/usr/bin/env python3
"""Generate full-BED IMPG similarity jobs for Fig5 raw many:many PAFs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import shlex
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan"
WFMASH_MANIFEST = ROOT / (
    "paper_prep/_brainstorming/"
    "pedigree_whole_genome_wfmash_p95_updated_bin/summaries/"
    "query_grid_filter_manifest.tsv"
)
SWEEPGA_MANIFEST = ROOT / (
    "paper_prep/_brainstorming/"
    "pedigree_whole_genome_sweepga_fastga_frequency32/summaries/"
    "query_grid_chop_filter_manifest.tsv"
)
INPUT_MANIFEST = ROOT / (
    "paper_prep/_brainstorming/"
    "pedigree_whole_genome_wfmash_p95_updated_bin/summaries/input_manifest.tsv"
)
IMPG = Path("/home/erikg/.cargo/bin/impg")


@dataclass(frozen=True)
class Job:
    method: str
    comparison_id: str
    source_raw_paf: Path
    impg_alignment_paf: Path
    query_fasta: Path
    target_fasta: Path
    target_bed: Path
    output_tsv_gz: Path
    metadata_json: Path
    stdout_log: Path
    stderr_log: Path
    job_script: Path
    cpus: int
    partition: str

    @property
    def command(self) -> list[str]:
        return [
            str(IMPG),
            "similarity",
            "--alignment-files",
            str(self.impg_alignment_paf),
            "--target-bed",
            str(self.target_bed),
            "--sequence-files",
            str(self.query_fasta),
            str(self.target_fasta),
            "--gfa-engine",
            "poa",
            "--no-merge",
            "--num-mappings",
            "many:many",
            "--scaffold-jump",
            "0",
            "--threads",
            "${SLURM_CPUS_PER_TASK}",
        ]

    @property
    def command_text(self) -> str:
        rendered = []
        for part in self.command:
            if part == "${SLURM_CPUS_PER_TASK}":
                rendered.append(part)
            else:
                rendered.append(shlex.quote(part))
        return " ".join(rendered)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha256_file(path: Path) -> str:
    import hashlib

    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_bed_from_fai(target_fasta: Path, out_bed: Path, tile_size: int) -> int:
    fai = Path(f"{target_fasta}.fai")
    if not fai.exists():
        raise FileNotFoundError(f"Missing FASTA index: {fai}")
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with fai.open() as src, out_bed.open("w") as dst:
        for line in src:
            seq, length_text, *_ = line.rstrip("\n").split("\t")
            length = int(length_text)
            for start in range(0, length, tile_size):
                end = min(start + tile_size, length)
                dst.write(f"{seq}\t{start}\t{end}\n")
                n += 1
    return n


def unique_raw_rows(rows: list[dict[str, str]], comparison_col: str) -> list[dict[str, str]]:
    seen: set[tuple[str, str]] = set()
    out: list[dict[str, str]] = []
    for row in rows:
        key = (row[comparison_col], row["raw_paf"])
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return sorted(out, key=lambda r: r[comparison_col])


def build_jobs(cpus: int, partition: str, tile_size: int) -> tuple[list[Job], list[dict[str, str]]]:
    inputs = {row["comparison_id"]: row for row in read_tsv(INPUT_MANIFEST)}
    raw_specs = [
        ("wfmash_p95_updated_bin", read_tsv(WFMASH_MANIFEST), "comparison_id"),
        ("sweepga_fastga_frequency32", read_tsv(SWEEPGA_MANIFEST), "comparison_id"),
    ]

    jobs: list[Job] = []
    bed_rows: list[dict[str, str]] = []
    for method, rows, comparison_col in raw_specs:
        for row in unique_raw_rows(rows, comparison_col):
            comparison_id = row[comparison_col]
            if comparison_id not in inputs:
                raise KeyError(f"{comparison_id} missing from {INPUT_MANIFEST}")
            input_row = inputs[comparison_id]
            target_fasta = Path(input_row["target_fasta"])
            query_fasta = Path(input_row["query_fasta"])
            bed = OUT / "beds" / f"{comparison_id}.full_genome_{tile_size // 1000}kb.bed"
            n_tiles = write_bed_from_fai(target_fasta, bed, tile_size)
            bed_rows.append(
                {
                    "comparison_id": comparison_id,
                    "target_fasta": str(target_fasta),
                    "target_fai": f"{target_fasta}.fai",
                    "bed": str(bed),
                    "tile_size_bp": str(tile_size),
                    "tile_count": str(n_tiles),
                }
            )

            stem = f"{method}.{comparison_id}.full_genome_{tile_size // 1000}kb"
            source_raw_paf = Path(row["raw_paf"])
            impg_alignment_paf = source_raw_paf
            if method == "sweepga_fastga_frequency32":
                impg_alignment_paf = OUT / "bgzf_raw_paf" / source_raw_paf.name.replace(
                    ".paf.gz", ".bgzf.paf.gz"
                )
            jobs.append(
                Job(
                    method=method,
                    comparison_id=comparison_id,
                    source_raw_paf=source_raw_paf,
                    impg_alignment_paf=impg_alignment_paf,
                    query_fasta=query_fasta,
                    target_fasta=target_fasta,
                    target_bed=bed,
                    output_tsv_gz=OUT / "outputs" / f"{stem}.impg_similarity.tsv.gz",
                    metadata_json=OUT / "outputs" / f"{stem}.metadata.json",
                    stdout_log=OUT / "logs" / f"{stem}.%j.out",
                    stderr_log=OUT / "logs" / f"{stem}.%j.err",
                    job_script=OUT / "jobs" / f"{stem}.slurm.sh",
                    cpus=cpus,
                    partition=partition,
                )
            )
    return jobs, bed_rows


def write_job_script(job: Job) -> None:
    job.job_script.parent.mkdir(parents=True, exist_ok=True)
    job.output_tsv_gz.parent.mkdir(parents=True, exist_ok=True)
    job.stdout_log.parent.mkdir(parents=True, exist_ok=True)
    out_text = str(job.output_tsv_gz)
    tmp_tsv = out_text[:-3] if out_text.endswith(".gz") else f"{out_text}.tmp"
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name=impg_{job.method[:6]}_{job.comparison_id[:14]}",
        f"#SBATCH --partition={job.partition}",
        f"#SBATCH --cpus-per-task={job.cpus}",
        "#SBATCH --nodes=1",
        "#SBATCH --time=24:00:00",
        f"#SBATCH --output={job.stdout_log}",
        f"#SBATCH --error={job.stderr_log}",
        "set -euo pipefail",
        "export LC_ALL=C",
        f"mkdir -p {shlex.quote(str(job.output_tsv_gz.parent))}",
        f"tmp_tsv={shlex.quote(tmp_tsv)}",
        f"out_gz={shlex.quote(str(job.output_tsv_gz))}",
        f"metadata={shlex.quote(str(job.metadata_json))}",
        f"command_text={shlex.quote(job.command_text)}",
        "start_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)",
        "impg_version=$(/home/erikg/.cargo/bin/impg --version)",
        "node=${SLURMD_NODENAME:-$(hostname)}",
        "partition=${SLURM_JOB_PARTITION:-unknown}",
        "cat > \"$metadata\" <<JSON",
        "{",
        f"  \"method\": {json.dumps(job.method)},",
        f"  \"comparison_id\": {json.dumps(job.comparison_id)},",
        f"  \"source_raw_paf\": {json.dumps(str(job.source_raw_paf))},",
        f"  \"impg_alignment_paf\": {json.dumps(str(job.impg_alignment_paf))},",
        f"  \"query_fasta\": {json.dumps(str(job.query_fasta))},",
        f"  \"target_fasta\": {json.dumps(str(job.target_fasta))},",
        f"  \"target_bed\": {json.dumps(str(job.target_bed))},",
        f"  \"output_tsv_gz\": {json.dumps(str(job.output_tsv_gz))},",
        f"  \"impg_path\": {json.dumps(str(IMPG))},",
        "  \"impg_version\": \"$impg_version\",",
        "  \"slurm_job_id\": \"${SLURM_JOB_ID:-manual}\",",
        "  \"slurm_cpus_per_task\": \"${SLURM_CPUS_PER_TASK:-unset}\",",
        "  \"node\": \"$node\",",
        "  \"partition\": \"$partition\",",
        "  \"start_utc\": \"$start_utc\",",
        "  \"command\": \"$command_text\"",
        "}",
        "JSON",
        "echo \"$command_text\"",
        f"{job.command_text} > \"$tmp_tsv\"",
        "gzip -f \"$tmp_tsv\"",
        "finish_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)",
        "python3 - <<'PY' \"$metadata\" \"$finish_utc\"",
        "import json, sys",
        "path, finish = sys.argv[1:]",
        "with open(path) as fh:",
        "    data = json.load(fh)",
        "data['finish_utc'] = finish",
        "data['status'] = 'OK'",
        "with open(path, 'w') as fh:",
        "    json.dump(data, fh, indent=2, sort_keys=True)",
        "    fh.write('\\n')",
        "PY",
    ]
    job.job_script.write_text("\n".join(lines) + "\n")
    job.job_script.chmod(0o755)


def write_bgzf_scripts(jobs: list[Job], cpus: int, partition: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for job in jobs:
        if job.source_raw_paf == job.impg_alignment_paf:
            continue
        script = OUT / "jobs" / f"normalize_bgzf.{job.comparison_id}.slurm.sh"
        stdout = OUT / "logs" / f"normalize_bgzf.{job.comparison_id}.%j.out"
        stderr = OUT / "logs" / f"normalize_bgzf.{job.comparison_id}.%j.err"
        script.parent.mkdir(parents=True, exist_ok=True)
        job.impg_alignment_paf.parent.mkdir(parents=True, exist_ok=True)
        lines = [
            "#!/usr/bin/env bash",
            f"#SBATCH --job-name=bgzf_{job.comparison_id[:14]}",
            f"#SBATCH --partition={partition}",
            f"#SBATCH --cpus-per-task={cpus}",
            "#SBATCH --nodes=1",
            "#SBATCH --time=12:00:00",
            f"#SBATCH --output={stdout}",
            f"#SBATCH --error={stderr}",
            "set -euo pipefail",
            "export LC_ALL=C",
            f"source_paf={shlex.quote(str(job.source_raw_paf))}",
            f"dest_paf={shlex.quote(str(job.impg_alignment_paf))}",
            'tmp_paf="${dest_paf}.tmp.${SLURM_JOB_ID:-manual}"',
            "if [ -s \"$dest_paf\" ]; then",
            "  /home/erikg/.guix-profile/bin/bgzip -t \"$dest_paf\"",
            "  echo \"BGZF already exists: $dest_paf\"",
            "  exit 0",
            "fi",
            "gzip -dc \"$source_paf\" | /home/erikg/.guix-profile/bin/bgzip -@ ${SLURM_CPUS_PER_TASK} -c > \"$tmp_paf\"",
            "/home/erikg/.guix-profile/bin/bgzip -t \"$tmp_paf\"",
            "mv \"$tmp_paf\" \"$dest_paf\"",
        ]
        script.write_text("\n".join(lines) + "\n")
        script.chmod(0o755)
        rows.append(
            {
                "comparison_id": job.comparison_id,
                "source_raw_paf": str(job.source_raw_paf),
                "bgzf_paf": str(job.impg_alignment_paf),
                "job_script": str(script),
                "stdout_log": str(stdout),
                "stderr_log": str(stderr),
                "partition": partition,
                "cpus_per_task": str(cpus),
                "submitted_job_id": "",
            }
        )
    return rows


def write_table(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def submit(job: Job) -> str:
    proc = subprocess.run(
        ["sbatch", "--parsable", str(job.job_script)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return proc.stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--submit", action="store_true", help="Submit generated Slurm scripts")
    parser.add_argument("--cpus", type=int, default=48)
    parser.add_argument("--partition", default="workers")
    parser.add_argument("--tile-size", type=int, default=10_000)
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    jobs, bed_rows = build_jobs(args.cpus, args.partition, args.tile_size)
    if len(jobs) != 6:
        raise RuntimeError(f"Expected 6 primary jobs, generated {len(jobs)}")
    for job in jobs:
        for path in (job.source_raw_paf, job.query_fasta, job.target_fasta, job.target_bed):
            if not path.exists():
                raise FileNotFoundError(path)
        write_job_script(job)
    bgzf_rows = write_bgzf_scripts(jobs, min(args.cpus, 48), args.partition)

    version = subprocess.run([str(IMPG), "--version"], check=True, text=True, stdout=subprocess.PIPE)
    help_text = subprocess.run(
        [str(IMPG), "similarity", "--help"], check=True, text=True, stdout=subprocess.PIPE
    )
    (OUT / "impg_similarity_help.txt").write_text(help_text.stdout)

    run_rows = []
    for job in jobs:
        run_rows.append(
            {
                "method": job.method,
                "comparison_id": job.comparison_id,
                "source_raw_paf": str(job.source_raw_paf),
                "impg_alignment_paf": str(job.impg_alignment_paf),
                "source_raw_paf_exists": str(job.source_raw_paf.exists()),
                "source_raw_paf_sha256": sha256_file(job.source_raw_paf),
                "impg_alignment_paf_exists": str(job.impg_alignment_paf.exists()),
                "query_fasta": str(job.query_fasta),
                "target_fasta": str(job.target_fasta),
                "target_bed": str(job.target_bed),
                "job_script": str(job.job_script),
                "output_tsv_gz": str(job.output_tsv_gz),
                "metadata_json": str(job.metadata_json),
                "stdout_log": str(job.stdout_log),
                "stderr_log": str(job.stderr_log),
                "partition": job.partition,
                "cpus_per_task": str(job.cpus),
                "impg_path": str(IMPG),
                "impg_version": version.stdout.strip(),
                "command": job.command_text,
                "submitted_job_id": "",
            }
        )

    if args.submit:
        bgzf_by_comparison: dict[str, str] = {}
        for row in bgzf_rows:
            proc = subprocess.run(
                ["sbatch", "--parsable", row["job_script"]],
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            row["submitted_job_id"] = proc.stdout.strip()
            bgzf_by_comparison[row["comparison_id"]] = row["submitted_job_id"]
        for row, job in zip(run_rows, jobs):
            if job.comparison_id in bgzf_by_comparison and job.method == "sweepga_fastga_frequency32":
                proc = subprocess.run(
                    [
                        "sbatch",
                        "--parsable",
                        f"--dependency=afterok:{bgzf_by_comparison[job.comparison_id]}",
                        str(job.job_script),
                    ],
                    check=True,
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                row["submitted_job_id"] = proc.stdout.strip()
                row["dependency"] = f"afterok:{bgzf_by_comparison[job.comparison_id]}"
            else:
                row["submitted_job_id"] = submit(job)

    write_table(
        OUT / "run_manifest.tsv",
        run_rows,
        [
            "method",
            "comparison_id",
            "source_raw_paf",
            "impg_alignment_paf",
            "source_raw_paf_exists",
            "source_raw_paf_sha256",
            "impg_alignment_paf_exists",
            "query_fasta",
            "target_fasta",
            "target_bed",
            "job_script",
            "output_tsv_gz",
            "metadata_json",
            "stdout_log",
            "stderr_log",
            "partition",
            "cpus_per_task",
            "impg_path",
            "impg_version",
            "command",
            "submitted_job_id",
            "dependency",
        ],
    )
    write_table(
        OUT / "bgzf_manifest.tsv",
        bgzf_rows,
        [
            "comparison_id",
            "source_raw_paf",
            "bgzf_paf",
            "job_script",
            "stdout_log",
            "stderr_log",
            "partition",
            "cpus_per_task",
            "submitted_job_id",
        ],
    )
    write_table(
        OUT / "bed_manifest.tsv",
        bed_rows,
        ["comparison_id", "target_fasta", "target_fai", "bed", "tile_size_bp", "tile_count"],
    )
    probe = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "impg_path": str(IMPG),
        "impg_version": version.stdout.strip(),
        "merge_distance_probe": "IMPG 0.4.1 rejects --merge-distance 0 with --no-merge; jobs use --no-merge.",
        "sequence_file_requirement": (
            "IMPG 0.4.1 similarity requires --sequence-files for POA similarity output, "
            "so query_fasta and target_fasta from input_manifest.tsv are passed explicitly."
        ),
    }
    (OUT / "pipeline_metadata.json").write_text(json.dumps(probe, indent=2, sort_keys=True) + "\n")
    print(f"Generated {len(jobs)} jobs under {OUT}")


if __name__ == "__main__":
    main()
