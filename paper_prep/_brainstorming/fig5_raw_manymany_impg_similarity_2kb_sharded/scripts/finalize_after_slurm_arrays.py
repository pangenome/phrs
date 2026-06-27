#!/usr/bin/env python3
"""Finalize Fig5 IMPG shards after Slurm array jobs have reached terminal state."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


SUCCESS_STATES = {"COMPLETED"}


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def check_arrays(array_ids: list[str]) -> None:
    cmd = [
        "sacct",
        "-n",
        "-P",
        "-j",
        ",".join(array_ids),
        "--format=JobID,State,ExitCode",
    ]
    proc = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE)
    bad: list[str] = []
    seen_tasks = 0
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        job_id, state, exit_code = line.split("|")[:3]
        if "." in job_id:
            continue
        if "_" not in job_id:
            continue
        if "[" in job_id:
            bad.append(f"{job_id} {state} {exit_code}")
            continue
        seen_tasks += 1
        if state not in SUCCESS_STATES or not exit_code.startswith("0:"):
            bad.append(f"{job_id} {state} {exit_code}")
    if bad:
        sample = "\n".join(bad[:50])
        raise SystemExit(f"Slurm arrays are not all successful; first bad states:\n{sample}")
    if seen_tasks == 0:
        raise SystemExit("No completed Slurm array task records found in sacct")


def normalize_tmp_outputs(live_base: Path) -> tuple[int, list[str]]:
    rows = read_tsv(live_base / "manifests" / "shard_manifest.tsv")
    renamed = 0
    missing: list[str] = []
    for row in rows:
        output = Path(row["output_tsv_gz"])
        if output.exists() and output.stat().st_size > 0:
            continue
        candidates = sorted(
            output.parent.glob(output.name.removesuffix(".gz") + ".tmp.*.gz"),
            key=lambda path: (path.stat().st_mtime, path.name),
        )
        if not candidates:
            missing.append(str(output))
            continue
        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(candidates[-1]), str(output))
        renamed += 1
    return renamed, missing


def rsync_tree(live_base: Path, main_base: Path) -> None:
    if live_base.resolve() == main_base.resolve():
        return
    main_base.mkdir(parents=True, exist_ok=True)
    for rel in ["outputs", "metadata", "summaries"]:
        src = live_base / rel
        if src.exists():
            subprocess.run(["rsync", "-a", f"{src}/", f"{main_base / rel}/"], check=True)
    for rel in [
        "manifests/shard_completion_manifest.tsv",
        "manifests/assembled_outputs.tsv",
        "REPORT.md",
    ]:
        src = live_base / rel
        dst = main_base / rel
        if src.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)


def append_report(live_base: Path, text: str) -> None:
    report = live_base / "REPORT.md"
    with report.open("a") as handle:
        handle.write("\n\n## Slurm dependency finalization\n\n")
        handle.write(text.rstrip() + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--live-base", required=True, type=Path)
    parser.add_argument("--main-base", required=True, type=Path)
    parser.add_argument("--array-ids", required=True, nargs="+")
    args = parser.parse_args()

    live_base = args.live_base.resolve()
    main_base = args.main_base.resolve()
    check_arrays(args.array_ids)
    renamed, missing = normalize_tmp_outputs(live_base)
    if missing:
        sample = "\n".join(missing[:50])
        raise SystemExit(f"{len(missing)} shard outputs still missing after tmp normalization; first missing:\n{sample}")

    finalizer = live_base / "scripts" / "finalize_2kb_sharded_impg.py"
    subprocess.run(["python3", str(finalizer)], check=True)
    append_report(
        live_base,
        "\n".join(
            [
                f"- Finalized at: `{utc_now()}`",
                f"- Slurm array dependency IDs: `{','.join(args.array_ids)}`",
                f"- Tmp shard outputs renamed into manifest paths: `{renamed}`",
                "- Plotting summaries select one best interchromosomal hit per 2 kb target window.",
                "- Tie-breaking: estimated.identity, then intersection length, dice, cosine, jaccard, then stable lexical target coordinates.",
                f"- Live output tree: `{live_base}`",
                f"- Synced output tree: `{main_base}`",
            ]
        ),
    )
    rsync_tree(live_base, main_base)


if __name__ == "__main__":
    main()
