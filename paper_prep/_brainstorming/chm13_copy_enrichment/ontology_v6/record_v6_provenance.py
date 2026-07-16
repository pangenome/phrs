#!/usr/bin/env python3
"""Freeze Slurm accounting, logs, code stages, RNG, and completion evidence."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Sequence

import v6_engine as V


ROLE = {
    "initial": ("slurm_v6_initial.sbatch", "initial"),
    "finalizer": ("slurm_v6_finalize.sbatch", "finalize"),
    "targeted": ("slurm_v6_targeted_extension.sbatch", "targeted"),
}


def collect(specifications: Sequence[str], output: Path = V.RELEASE) -> dict[str, object]:
    jobs = []
    for specification in specifications:
        role, job_id = specification.split("=", 1)
        if role not in ROLE or not job_id.isdigit():
            raise ValueError(f"invalid job specification: {specification}")
        jobs.append((role, job_id))
    fields = [
        "role", "job_id", "scheduler_job_id", "job_name", "partition", "state", "exit_code",
        "elapsed", "allocated_cpus", "requested_memory", "max_rss", "node_list", "script",
        "seed", "permutations", "stdout_log", "stderr_log", "stdout_bytes", "stdout_sha256",
        "stderr_bytes", "stderr_sha256", "query_utc",
    ]
    rows = []
    roots = {}
    commands = []
    for role, job_id in jobs:
        command = [
            "sacct", "-j", job_id, "--parsable2", "--noheader",
            "--format=JobIDRaw,JobName,Partition,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList",
        ]
        completed = subprocess.run(command, text=True, capture_output=True, check=True)
        commands.append(command)
        script, log_stem = ROLE[role]
        stdout = V.RESULTS / "logs" / f"{log_stem}-{job_id}.out"
        stderr = V.RESULTS / "logs" / f"{log_stem}-{job_id}.err"
        for line in completed.stdout.splitlines():
            values = line.split("|")
            if len(values) < 10 or not values[0]:
                continue
            row = {
                "role": role, "job_id": job_id, "scheduler_job_id": values[0],
                "job_name": values[1], "partition": values[2], "state": values[3],
                "exit_code": values[4], "elapsed": values[5], "allocated_cpus": values[6],
                "requested_memory": values[7], "max_rss": values[8], "node_list": values[9],
                "script": str((V.HERE / script).relative_to(V.REPO)),
                "seed": V.MASTER_SEED,
                "permutations": V.INITIAL_PERMUTATIONS if role != "targeted" else "staged_as_required",
                "stdout_log": str(stdout.relative_to(V.REPO)),
                "stderr_log": str(stderr.relative_to(V.REPO)),
                "stdout_bytes": stdout.stat().st_size if stdout.is_file() else "",
                "stdout_sha256": V.sha256(stdout) if stdout.is_file() else "",
                "stderr_bytes": stderr.stat().st_size if stderr.is_file() else "",
                "stderr_sha256": V.sha256(stderr) if stderr.is_file() else "",
                "query_utc": V.utcnow(),
            }
            rows.append(row)
            if values[0] == job_id:
                roots[role] = row
    V.write_rows(output / "SLURM_JOBS.tsv", fields, rows)
    complete = all(
        roots.get(role, {}).get("state", "").startswith("COMPLETED")
        and roots.get(role, {}).get("exit_code") == "0:0"
        for role, _job_id in jobs
    )
    completion = {
        "schema_version": V.SCHEMA_VERSION,
        "jobs": {role: {
            "job_id": job_id,
            "state": roots.get(role, {}).get("state", "MISSING"),
            "exit_code": roots.get(role, {}).get("exit_code", "MISSING"),
            "elapsed": roots.get(role, {}).get("elapsed", "MISSING"),
            "requested_resources": {
                "partition": "workers", "cpus": 1,
                "memory": "16G" if role == "targeted" else "32G", "time": "12:00:00",
            },
        } for role, job_id in jobs},
        "completion_status": "COMPLETE" if complete else "FAILED_OR_INCOMPLETE",
        "accounting_commands": commands,
        "recorded_utc": V.utcnow(),
    }
    V.atomic_json(output / "SLURM_COMPLETION.json", completion)

    initial_manifest = V.WORK / "primary/run_manifest.json"
    selective_manifest = V.WORK / "selective_extension/run_manifest.json"
    provenance = {
        "schema_version": V.SCHEMA_VERSION,
        "array_manifest_path": str(initial_manifest.relative_to(V.REPO)),
        "array_manifest_sha256": V.sha256(initial_manifest),
        "array_manifest": json.loads(initial_manifest.read_text()),
        "array_engine_sha256": json.loads(initial_manifest.read_text())["immutable_configuration"]["input_checksums"]["v6_engine"],
        "release_engine_sha256": V.sha256(V.HERE / "v6_engine.py"),
        "targeted_engine_sha256": V.sha256(V.HERE / "v6_targeted_extension.py"),
        "current_git_commit": V.git_commit(),
        "initial_seed": V.MASTER_SEED,
        "bit_generator": "PCG64DXSM",
        "spawn_key": [0],
        "rng_state_saved_after_every_batch": True,
        "selective_extension_manifest_path": (
            str(selective_manifest.relative_to(V.REPO)) if selective_manifest.is_file() else "NOT_REQUIRED"
        ),
        "selective_extension_manifest_sha256": (
            V.sha256(selective_manifest) if selective_manifest.is_file() else "NOT_REQUIRED"
        ),
        "slurm_completion_status": completion["completion_status"],
        "scripts": {role: str((V.HERE / ROLE[role][0]).relative_to(V.REPO)) for role, _ in jobs},
        "recorded_utc": V.utcnow(),
    }
    V.atomic_json(output / "COMPUTATIONAL_PROVENANCE.json", provenance)
    return completion


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", action="append", required=True,
                        help="role=job_id; roles: initial, finalizer, targeted")
    parser.add_argument("--output", type=Path, default=V.RELEASE)
    args = parser.parse_args(argv)
    report = collect(args.job, args.output)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["completion_status"] == "COMPLETE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
