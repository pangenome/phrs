#!/usr/bin/env python3
"""Keep same-chromosome and interchromosomal IMPG winners per query window."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


CHR_RE = re.compile(r"(?:^|[#_/])chr([0-9]+|X|Y|M)(?:[._#:/-]|$)")


def chrom_name(seq: str) -> str:
    match = CHR_RE.search(seq)
    if match:
        return f"chr{match.group(1)}"
    match = re.search(r"chr([0-9]+|X|Y|M)", seq)
    if match:
        return f"chr{match.group(1)}"
    return seq


def split_group(group: str) -> tuple[str, int, int]:
    seq, coords = group.rsplit(":", 1)
    start, end = coords.split("-", 1)
    return seq, int(start), int(end)


def query_and_other(row: dict[str, str]) -> tuple[str, str] | None:
    query_seq = row["chrom"]
    query_start = int(row["start"])
    query_end = int(row["end"])
    a_seq, a_start, a_end = split_group(row["group.a"])
    b_seq, b_start, b_end = split_group(row["group.b"])
    if (a_seq, a_start, a_end) == (query_seq, query_start, query_end):
        return query_seq, b_seq
    if (b_seq, b_start, b_end) == (query_seq, query_start, query_end):
        return query_seq, a_seq
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
    writer: csv.DictWriter[str],
    skip_writer: csv.writer,
    rows: list[dict[str, str]],
    max_candidates: int,
) -> None:
    if not rows:
        return
    first = rows[0]
    raw_count = len(rows)
    if raw_count > max_candidates:
        skip_writer.writerow([first["chrom"], first["start"], first["end"], raw_count, 0, 0, "too_many_candidates"])
        return

    best: dict[str, tuple[tuple[object, ...], dict[str, str], str, str]] = {}
    class_counts = {"same_chrom": 0, "interchrom": 0}
    for row in rows:
        sides = query_and_other(row)
        if sides is None:
            continue
        query_seq, other_seq = sides
        if row["group.a"] == row["group.b"]:
            continue
        query_chrom = chrom_name(query_seq)
        other_chrom = chrom_name(other_seq)
        winner_class = "same_chrom" if query_chrom == other_chrom else "interchrom"
        class_counts[winner_class] += 1
        rr = rank(row)
        if winner_class not in best or rr > best[winner_class][0]:
            best[winner_class] = (rr, row, query_seq, other_seq)

    if not best:
        skip_writer.writerow([first["chrom"], first["start"], first["end"], raw_count, 0, 0, "no_usable_candidates"])
        return

    for winner_class in ("same_chrom", "interchrom"):
        if winner_class not in best:
            continue
        _rr, row, query_seq, other_seq = best[winner_class]
        out = dict(row)
        out["winner_class"] = winner_class
        out["query_seq"] = query_seq
        out["other_seq"] = other_seq
        out["query_chrom"] = chrom_name(query_seq)
        out["other_chrom"] = chrom_name(other_seq)
        out["raw_candidate_count"] = raw_count
        out["class_candidate_count"] = class_counts[winner_class]
        writer.writerow(out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-candidates", type=int, default=5000)
    parser.add_argument("--skip-report", required=True, type=Path)
    args = parser.parse_args()

    args.skip_report.parent.mkdir(parents=True, exist_ok=True)
    with args.skip_report.open("w", newline="") as skip_handle:
        skip_writer = csv.writer(skip_handle, delimiter="\t", lineterminator="\n")
        skip_writer.writerow([
            "chrom",
            "start",
            "end",
            "raw_candidate_count",
            "same_chrom_candidate_count",
            "interchrom_candidate_count",
            "reason",
        ])

        reader = csv.DictReader(sys.stdin, delimiter="\t")
        if reader.fieldnames is None:
            return
        extra_fields = [
            "winner_class",
            "query_seq",
            "other_seq",
            "query_chrom",
            "other_chrom",
            "raw_candidate_count",
            "class_candidate_count",
        ]
        writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=reader.fieldnames + extra_fields, lineterminator="\n")
        writer.writeheader()

        current_key: tuple[str, str, str] | None = None
        group_rows: list[dict[str, str]] = []
        for row in reader:
            key = (row["chrom"], row["start"], row["end"])
            if current_key is not None and key != current_key:
                emit_group(writer, skip_writer, group_rows, args.max_candidates)
                group_rows = []
            current_key = key
            group_rows.append(row)
        emit_group(writer, skip_writer, group_rows, args.max_candidates)


if __name__ == "__main__":
    main()
