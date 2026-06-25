#!/usr/bin/env python3
"""Summarize focused odgi untangle support for the Fig5 graph pilot."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def maybe_int(value: str) -> int:
    try:
        return int(value)
    except ValueError:
        return 0


def parse_bedpe(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rows.append(
                {
                    "query_path": parts[0],
                    "query_start": parts[1],
                    "query_end": parts[2],
                    "target_path": parts[3],
                    "target_start": parts[4],
                    "target_end": parts[5],
                    "raw_columns": "|".join(parts),
                }
            )
    return rows


def classify_segments(
    run_id: str, bedpe_rows: list[dict[str, str]], focus_rows: list[dict[str, str]]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    focus_compiled = [
        (
            row,
            re.compile(row["query_regex"]),
            re.compile(row["target_regex"]),
        )
        for row in focus_rows
    ]
    segments: list[dict[str, object]] = []
    summary: dict[str, dict[str, object]] = {}
    for row in focus_rows:
        summary[row["focus_id"]] = {
            "run_id": run_id,
            "focus_id": row["focus_id"],
            "class": row["class"],
            "query_regex": row["query_regex"],
            "target_regex": row["target_regex"],
            "segment_count": 0,
            "total_query_bp": 0,
            "max_segment_bp": 0,
            "query_paths": set(),
            "target_paths": set(),
            "interpretation": row["interpretation"],
        }

    for row in bedpe_rows:
        query_bp = max(0, maybe_int(row["query_end"]) - maybe_int(row["query_start"]))
        for focus, query_re, target_re in focus_compiled:
            if not query_re.search(row["query_path"]):
                continue
            if not target_re.search(row["target_path"]):
                continue
            focus_id = focus["focus_id"]
            segments.append(
                {
                    "run_id": run_id,
                    "focus_id": focus_id,
                    "class": focus["class"],
                    "query_path": row["query_path"],
                    "query_start": row["query_start"],
                    "query_end": row["query_end"],
                    "query_bp": query_bp,
                    "target_path": row["target_path"],
                    "target_start": row["target_start"],
                    "target_end": row["target_end"],
                    "raw_columns": row["raw_columns"],
                }
            )
            acc = summary[focus_id]
            acc["segment_count"] = int(acc["segment_count"]) + 1
            acc["total_query_bp"] = int(acc["total_query_bp"]) + query_bp
            acc["max_segment_bp"] = max(int(acc["max_segment_bp"]), query_bp)
            acc["query_paths"].add(row["query_path"])
            acc["target_paths"].add(row["target_path"])

    summary_rows: list[dict[str, object]] = []
    for focus in focus_rows:
        acc = summary[focus["focus_id"]]
        summary_rows.append(
            {
                **acc,
                "query_path_count": len(acc["query_paths"]),
                "target_path_count": len(acc["target_paths"]),
                "query_paths": ";".join(sorted(acc["query_paths"])),
                "target_paths": ";".join(sorted(acc["target_paths"])),
            }
        )
    return summary_rows, segments


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--paths", required=True, type=Path)
    parser.add_argument("--bedpe", required=True, type=Path)
    parser.add_argument("--focus", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--segments", required=True, type=Path)
    args = parser.parse_args()

    focus_rows = read_tsv(args.focus)
    bedpe_rows = parse_bedpe(args.bedpe)
    summary_rows, segment_rows = classify_segments(args.run_id, bedpe_rows, focus_rows)

    write_tsv(
        args.summary,
        summary_rows,
        [
            "run_id",
            "focus_id",
            "class",
            "query_regex",
            "target_regex",
            "segment_count",
            "total_query_bp",
            "max_segment_bp",
            "query_path_count",
            "target_path_count",
            "query_paths",
            "target_paths",
            "interpretation",
        ],
    )
    write_tsv(
        args.segments,
        segment_rows,
        [
            "run_id",
            "focus_id",
            "class",
            "query_path",
            "query_start",
            "query_end",
            "query_bp",
            "target_path",
            "target_start",
            "target_end",
            "raw_columns",
        ],
    )


if __name__ == "__main__":
    main()
