#!/usr/bin/env python3
import collections
import glob
import os
import re

from common import PACKAGE_DIR, paf_glob, paf_records, write_tsv


def chrom_from_name(name):
    match = re.search(r"#chr([0-9XYM]+)", name)
    if match:
        return "chr" + match.group(1)
    match = re.search(r"chr([0-9XYM]+)", name)
    if match:
        return "chr" + match.group(1)
    return name


def summarize_file(path):
    records = 0
    target_counts = collections.Counter()
    for f in paf_records(path):
        records += 1
        target_counts[chrom_from_name(f[5])] += 1
    rel_parent = os.path.basename(os.path.dirname(path))
    base = os.path.basename(path).split(".")
    return [
        {
            "paf_path": path,
            "run_label": rel_parent,
            "comparison_id": base[0],
            "parameter_set": base[1],
            "target_chrom": chrom,
            "rows": str(count),
            "total_rows": str(records),
        }
        for chrom, count in sorted(target_counts.items())
    ] or [{
        "paf_path": path,
        "run_label": rel_parent,
        "comparison_id": base[0],
        "parameter_set": base[1],
        "target_chrom": "NO_ROWS",
        "rows": "0",
        "total_rows": "0",
    }]


def main():
    rows = []
    for path in sorted(glob.glob(paf_glob(), recursive=True)):
        rows.extend(summarize_file(path))
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "paf_file_summary.tsv"), rows, [
        "paf_path",
        "run_label",
        "comparison_id",
        "parameter_set",
        "target_chrom",
        "rows",
        "total_rows",
    ])


if __name__ == "__main__":
    main()
