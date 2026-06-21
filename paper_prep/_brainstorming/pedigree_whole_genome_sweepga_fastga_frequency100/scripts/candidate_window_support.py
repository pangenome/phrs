#!/usr/bin/env python3
import glob
import os
from collections import defaultdict

from common import PACKAGE_DIR, paf_records, write_tsv


WINDOWS = [
    {
        "candidate_id": "PAN027_chr9q_chr3q_PHR_candidate",
        "comparison_id": "PAN027pat_vs_PAN011_joint",
        "query_name": "PAN027#2#chr9",
        "query_start": 135704825,
        "query_end": 136204825,
    },
    {
        "candidate_id": "PAN028_chr9q_chr3q_PHR_candidate",
        "comparison_id": "PAN028mat_vs_PAN027_joint",
        "query_name": "PAN028#1#chr9",
        "query_start": 134380985,
        "query_end": 134880985,
    },
]


def target_chrom(target_name):
    if "#chr" not in target_name:
        if "_chr" not in target_name:
            return "unknown"
        return "chr" + target_name.rsplit("_chr", 1)[1].split("_", 1)[0].split("#", 1)[0]
    return "chr" + target_name.rsplit("#chr", 1)[1].split("_", 1)[0].split("#", 1)[0]


def union_bp(intervals):
    if not intervals:
        return 0
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return sum(end - start for start, end in merged)


def file_label(path):
    rel = os.path.relpath(path, PACKAGE_DIR)
    if rel.startswith("raw_paf/"):
        return "raw_many_many_j0"
    if rel.startswith("chopped_paf/"):
        base = os.path.basename(path)
        return "chopped_many_many_" + base.split(".chopped_", 1)[1].split(".paf.gz", 1)[0]
    if rel.startswith("chopped_paf_l"):
        dirname = rel.split("/", 1)[0]
        return "chopped_many_many_" + dirname.replace("chopped_paf_", "")
    if rel.startswith("filtered_paf/"):
        base = os.path.basename(path)
        parts = base.split(".")
        return parts[1] if len(parts) > 2 else "filtered"
    return rel


def summarize_window(path, window):
    by_target = defaultdict(lambda: {"rows": 0, "sum_overlap_bp": 0, "intervals": []})
    for f in paf_records(path):
        if f[0] != window["query_name"]:
            continue
        q_start = int(f[2])
        q_end = int(f[3])
        ov_start = max(q_start, window["query_start"])
        ov_end = min(q_end, window["query_end"])
        if ov_start >= ov_end:
            continue
        chrom = target_chrom(f[5])
        by_target[chrom]["rows"] += 1
        by_target[chrom]["sum_overlap_bp"] += ov_end - ov_start
        by_target[chrom]["intervals"].append((ov_start, ov_end))
    rows = []
    for chrom, vals in sorted(by_target.items()):
        rows.append({
            "candidate_id": window["candidate_id"],
            "comparison_id": window["comparison_id"],
            "layer": file_label(path),
            "paf_path": path,
            "query_name": window["query_name"],
            "query_start": str(window["query_start"]),
            "query_end": str(window["query_end"]),
            "target_chrom": chrom,
            "overlap_rows": str(vals["rows"]),
            "sum_overlap_bp": str(vals["sum_overlap_bp"]),
            "query_union_bp": str(union_bp(vals["intervals"])),
            "emits_chr3_target": "yes" if chrom == "chr3" and vals["rows"] else "no",
        })
    if not rows:
        rows.append({
            "candidate_id": window["candidate_id"],
            "comparison_id": window["comparison_id"],
            "layer": file_label(path),
            "paf_path": path,
            "query_name": window["query_name"],
            "query_start": str(window["query_start"]),
            "query_end": str(window["query_end"]),
            "target_chrom": "none",
            "overlap_rows": "0",
            "sum_overlap_bp": "0",
            "query_union_bp": "0",
            "emits_chr3_target": "no",
        })
    return rows


def main():
    rows = []
    paths = []
    paths.extend(glob.glob(os.path.join(PACKAGE_DIR, "raw_paf", "*.paf.gz")))
    paths.extend(glob.glob(os.path.join(PACKAGE_DIR, "chopped_paf_l*_o*", "*.paf.gz")))
    paths.extend(glob.glob(os.path.join(PACKAGE_DIR, "filtered_paf", "*.paf.gz")))
    for window in WINDOWS:
        for path in sorted(paths):
            if os.path.basename(path).split(".")[0] != window["comparison_id"]:
                continue
            rows.extend(summarize_window(path, window))
    write_tsv(
        os.path.join(PACKAGE_DIR, "summaries", "candidate_window_support.tsv"),
        rows,
        [
            "candidate_id",
            "comparison_id",
            "layer",
            "paf_path",
            "query_name",
            "query_start",
            "query_end",
            "target_chrom",
            "overlap_rows",
            "sum_overlap_bp",
            "query_union_bp",
            "emits_chr3_target",
        ],
    )


if __name__ == "__main__":
    main()
