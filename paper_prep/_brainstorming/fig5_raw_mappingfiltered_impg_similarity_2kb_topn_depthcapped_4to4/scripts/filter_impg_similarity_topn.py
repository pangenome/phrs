#!/usr/bin/env python3
"""Keep bounded top-N IMPG similarity rows per target window."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


CHR_RE = re.compile(r"(?:^|[_#])chr([0-9XY]+)(?::|$)")


def chrom_name(seq: str) -> str:
    match = CHR_RE.search(seq)
    return f"chr{match.group(1)}" if match else "unknown"


def split_group(group: str) -> tuple[str, int, int]:
    seq, coords = group.rsplit(":", 1)
    start, end = coords.split("-", 1)
    return seq, int(start), int(end)


def target_other(row: dict[str, str]) -> tuple[str, str] | None:
    target_seq = row["chrom"]
    target_start = int(row["start"])
    target_end = int(row["end"])
    a_seq, a_start, a_end = split_group(row["group.a"])
    b_seq, b_start, b_end = split_group(row["group.b"])
    if (a_seq, a_start, a_end) == (target_seq, target_start, target_end):
        return target_seq, b_seq
    if (b_seq, b_start, b_end) == (target_seq, target_start, target_end):
        return target_seq, a_seq
    return None


def rank(row: dict[str, str]) -> tuple[float, int, float, float, float, str, str]:
    return (
        float(row["estimated.identity"]),
        int(float(row["intersection"])),
        float(row["dice.similarity"]),
        float(row["cosine.similarity"]),
        float(row["jaccard.similarity"]),
        row["group.b"],
        row["group.a"],
    )


def emit_group(
    writer: csv.DictWriter,
    skip_writer: csv.writer,
    rows: list[dict[str, str]],
    top_n: int,
    max_candidates: int,
    interchrom_only: bool,
) -> tuple[int, int]:
    if not rows:
        return 0, 0
    first = rows[0]
    raw_count = len(rows)
    if raw_count > max_candidates:
        skip_writer.writerow([first["chrom"], first["start"], first["end"], raw_count, 0, "too_many_candidates"])
        return 0, 1

    kept: list[dict[str, str]] = []
    for row in rows:
        sides = target_other(row)
        if sides is None:
            continue
        target_seq, other_seq = sides
        if row["group.a"] == row["group.b"]:
            continue
        if interchrom_only and chrom_name(target_seq) == chrom_name(other_seq):
            continue
        kept.append(row)

    for row in sorted(kept, key=rank, reverse=True)[:top_n]:
        writer.writerow(row)
    return min(len(kept), top_n), 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--max-candidates", type=int, default=500)
    parser.add_argument("--skip-report", required=True, type=Path)
    parser.add_argument("--interchrom-only", action="store_true")
    args = parser.parse_args()

    args.skip_report.parent.mkdir(parents=True, exist_ok=True)
    with args.skip_report.open("w", newline="") as skip_handle:
        skip_writer = csv.writer(skip_handle, delimiter="\t", lineterminator="\n")
        skip_writer.writerow(["chrom", "start", "end", "raw_candidate_count", "retained_count", "reason"])

        reader = csv.DictReader(sys.stdin, delimiter="\t")
        if reader.fieldnames is None:
            return
        writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=reader.fieldnames, lineterminator="\n")
        writer.writeheader()

        current_key: tuple[str, str, str] | None = None
        group_rows: list[dict[str, str]] = []
        for row in reader:
            key = (row["chrom"], row["start"], row["end"])
            if current_key is not None and key != current_key:
                emit_group(writer, skip_writer, group_rows, args.top_n, args.max_candidates, args.interchrom_only)
                group_rows = []
            current_key = key
            group_rows.append(row)
        emit_group(writer, skip_writer, group_rows, args.top_n, args.max_candidates, args.interchrom_only)


if __name__ == "__main__":
    main()
