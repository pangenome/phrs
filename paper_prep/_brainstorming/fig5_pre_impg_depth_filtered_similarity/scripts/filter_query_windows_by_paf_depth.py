#!/usr/bin/env python3
"""Build a query BED after pre-IMPG PAF-depth screening."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
import subprocess
from array import array
from pathlib import Path
from typing import Iterable, TextIO


CHR_RE = re.compile(r"(?:^|[#_/])chr([0-9]+|X|Y|M)(?:[._#:/-]|$)")


def parse_chrom(name: str) -> str:
    match = CHR_RE.search(name)
    if match:
        return f"chr{match.group(1)}"
    match = re.search(r"chr([0-9]+|X|Y|M)", name)
    if match:
        return f"chr{match.group(1)}"
    return name


def open_maybe_gzip(path: Path, pigz_threads: int) -> tuple[subprocess.Popen[str] | None, TextIO]:
    if str(path).endswith(".gz") and shutil.which("pigz"):
        proc = subprocess.Popen(
            ["pigz", "-dc", "-p", str(max(1, pigz_threads)), str(path)],
            stdout=subprocess.PIPE,
            text=True,
        )
        assert proc.stdout is not None
        return proc, proc.stdout
    if str(path).endswith(".gz"):
        return None, gzip.open(path, "rt")
    return None, path.open()


def read_fai(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            seq, length, *_ = line.rstrip("\n").split("\t")
            rows.append((seq, int(length)))
    return rows


def load_centromeres(path: Path) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *rest = line.rstrip("\n").split("\t")
            label = "\t".join(rest).lower()
            if "centromere" not in label:
                continue
            out.setdefault(chrom, []).append((int(start), int(end)))
    return out


def overlaps_any(intervals: Iterable[tuple[int, int]], start: int, end: int) -> bool:
    for iv_start, iv_end in intervals:
        if start < iv_end and end > iv_start:
            return True
    return False


def increment_depth(depths: dict[str, array], seq: str, start: int, end: int, window: int, cap: int) -> int:
    arr = depths.get(seq)
    if arr is None or end <= start:
        return 0
    first = start // window
    last = (end - 1) // window
    touched = 0
    for idx in range(first, min(last + 1, len(arr))):
        if arr[idx] <= cap:
            arr[idx] += 1
        touched += 1
    return touched


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf", required=True, type=Path)
    parser.add_argument("--query-fai", required=True, type=Path)
    parser.add_argument("--centromere-bed", required=True, type=Path)
    parser.add_argument("--out-bed", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--window-size", type=int, default=2000)
    parser.add_argument("--min-depth", type=int, default=1)
    parser.add_argument("--max-depth", type=int, default=100)
    parser.add_argument("--interchrom-only", action="store_true")
    parser.add_argument("--pigz-threads", type=int, default=1)
    args = parser.parse_args()

    seqs = read_fai(args.query_fai)
    centromeres = load_centromeres(args.centromere_bed)
    cap = args.max_depth
    depths: dict[str, array] = {
        seq: array("H", [0]) * ((length + args.window_size - 1) // args.window_size)
        for seq, length in seqs
    }

    paf_rows = 0
    counted_rows = 0
    touched_windows = 0
    proc, handle = open_maybe_gzip(args.paf, args.pigz_threads)
    try:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            paf_rows += 1
            q_name = fields[0]
            t_name = fields[5]
            if args.interchrom_only and parse_chrom(q_name) == parse_chrom(t_name):
                continue
            try:
                q_start = int(fields[2])
                q_end = int(fields[3])
            except ValueError:
                continue
            hit = increment_depth(depths, q_name, q_start, q_end, args.window_size, cap)
            if hit:
                counted_rows += 1
                touched_windows += hit
    finally:
        handle.close()
        if proc is not None:
            rc = proc.wait()
            if rc != 0:
                raise RuntimeError(f"pigz failed with exit {rc} for {args.paf}")

    args.out_bed.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)

    total_windows = 0
    kept_windows = 0
    zero_depth_windows = 0
    high_depth_windows = 0
    centromere_windows = 0
    max_observed_depth = 0
    depth_sum_kept = 0
    with args.out_bed.open("w") as bed, args.report.open("w", newline="") as rep:
        writer = csv.writer(rep, delimiter="\t", lineterminator="\n")
        writer.writerow(["seq", "start", "end", "depth", "decision", "reason"])
        for seq, length in seqs:
            arr = depths[seq]
            chrom = parse_chrom(seq)
            for idx, depth in enumerate(arr):
                start = idx * args.window_size
                end = min(start + args.window_size, length)
                total_windows += 1
                max_observed_depth = max(max_observed_depth, int(depth))
                reason = "keep"
                decision = "KEEP"
                if overlaps_any(centromeres.get(chrom, []), start, end):
                    centromere_windows += 1
                    decision = "DROP"
                    reason = "centromere"
                elif depth < args.min_depth:
                    zero_depth_windows += 1
                    decision = "DROP"
                    reason = "below_min_depth"
                elif depth > args.max_depth:
                    high_depth_windows += 1
                    decision = "DROP"
                    reason = "above_max_depth"
                else:
                    kept_windows += 1
                    depth_sum_kept += int(depth)
                    bed.write(f"{seq}\t{start}\t{end}\n")
                writer.writerow([seq, start, end, int(depth), decision, reason])

    summary_fields = [
        "paf",
        "query_fai",
        "window_size",
        "min_depth",
        "max_depth",
        "interchrom_only",
        "paf_rows",
        "counted_rows",
        "touched_window_events",
        "total_windows",
        "kept_windows",
        "zero_depth_windows",
        "high_depth_windows",
        "centromere_windows",
        "max_observed_depth",
        "mean_kept_depth",
        "out_bed",
        "report",
    ]
    with args.summary.open("w", newline="") as handle_out:
        writer = csv.DictWriter(handle_out, delimiter="\t", fieldnames=summary_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "paf": args.paf,
                "query_fai": args.query_fai,
                "window_size": args.window_size,
                "min_depth": args.min_depth,
                "max_depth": args.max_depth,
                "interchrom_only": int(args.interchrom_only),
                "paf_rows": paf_rows,
                "counted_rows": counted_rows,
                "touched_window_events": touched_windows,
                "total_windows": total_windows,
                "kept_windows": kept_windows,
                "zero_depth_windows": zero_depth_windows,
                "high_depth_windows": high_depth_windows,
                "centromere_windows": centromere_windows,
                "max_observed_depth": max_observed_depth,
                "mean_kept_depth": f"{(depth_sum_kept / kept_windows) if kept_windows else 0:.6f}",
                "out_bed": args.out_bed,
                "report": args.report,
            }
        )


if __name__ == "__main__":
    main()

