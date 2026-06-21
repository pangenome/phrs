#!/usr/bin/env python3
import collections
import glob
import os
import re

from common import CANDIDATE_WINDOWS, PACKAGE_DIR, paf_glob, paf_records, read_tsv, write_tsv


def chrom_from_name(name):
    match = re.search(r"#chr([0-9XYM]+)", name)
    if match:
        return "chr" + match.group(1)
    match = re.search(r"chr([0-9XYM]+)", name)
    if match:
        return "chr" + match.group(1)
    return name


def overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def merge_coverage(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    merged = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def summarize_window(paf_path, window):
    by_chrom = collections.defaultdict(lambda: {"rows": 0, "overlap_sum": 0, "intervals": []})
    qname = window["query_name"]
    qstart = int(window["query_start"])
    qend = int(window["query_end"])
    for f in paf_records(paf_path):
        if f[0] != qname:
            continue
        row_qstart = int(f[2])
        row_qend = int(f[3])
        ov = overlap(row_qstart, row_qend, qstart, qend)
        if ov <= 0:
            continue
        chrom = chrom_from_name(f[5])
        entry = by_chrom[chrom]
        entry["rows"] += 1
        entry["overlap_sum"] += ov
        entry["intervals"].append((max(row_qstart, qstart), min(row_qend, qend)))
    base = os.path.basename(paf_path).split(".")
    run_label = os.path.basename(os.path.dirname(paf_path))
    rows = []
    for chrom, entry in sorted(by_chrom.items()):
        rows.append({
            "event_id": window["event_id"],
            "run_label": run_label,
            "comparison_id": window["comparison_id"],
            "parameter_set": base[1],
            "query_name": qname,
            "query_start": str(qstart),
            "query_end": str(qend),
            "target_chrom": chrom,
            "expected_target_chrom": window["expected_target_chrom"],
            "paf_rows_overlapping_window": str(entry["rows"]),
            "query_overlap_bp_sum": str(entry["overlap_sum"]),
            "query_covered_bp_union": str(merge_coverage(entry["intervals"])),
            "chr3_support": "yes" if chrom == "chr3" and entry["rows"] > 0 else "no",
            "paf_path": paf_path,
        })
    if not rows:
        rows.append({
            "event_id": window["event_id"],
            "run_label": run_label,
            "comparison_id": window["comparison_id"],
            "parameter_set": base[1],
            "query_name": qname,
            "query_start": str(qstart),
            "query_end": str(qend),
            "target_chrom": "NO_OVERLAP",
            "expected_target_chrom": window["expected_target_chrom"],
            "paf_rows_overlapping_window": "0",
            "query_overlap_bp_sum": "0",
            "query_covered_bp_union": "0",
            "chr3_support": "no",
            "paf_path": paf_path,
        })
    return rows


def main():
    windows = read_tsv(CANDIDATE_WINDOWS)
    pafs = sorted(glob.glob(paf_glob(), recursive=True))
    rows = []
    for window in windows:
        for paf in pafs:
            base = os.path.basename(paf).split(".")
            if base[0] != window["comparison_id"]:
                continue
            rows.extend(summarize_window(paf, window))
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "candidate_window_support.tsv"), rows, [
        "event_id",
        "run_label",
        "comparison_id",
        "parameter_set",
        "query_name",
        "query_start",
        "query_end",
        "target_chrom",
        "expected_target_chrom",
        "paf_rows_overlapping_window",
        "query_overlap_bp_sum",
        "query_covered_bp_union",
        "chr3_support",
        "paf_path",
    ])


if __name__ == "__main__":
    main()
