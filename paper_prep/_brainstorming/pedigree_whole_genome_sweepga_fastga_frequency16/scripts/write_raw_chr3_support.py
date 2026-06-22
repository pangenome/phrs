#!/usr/bin/env python3
import csv
import glob
import os

from candidate_window_support import WINDOWS, summarize_window
from common import PACKAGE_DIR, write_tsv


FIELDS = [
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
    "raw_paf_status",
    "note",
]


def read_existing(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    rows = []
    raw_paths = sorted(glob.glob(os.path.join(PACKAGE_DIR, "raw_paf", "*.paf.gz")))
    for window in WINDOWS:
        matched = [p for p in raw_paths if os.path.basename(p).split(".")[0] == window["comparison_id"]]
        if matched:
            for path in matched:
                for row in summarize_window(path, window):
                    row["layer"] = "raw_frequency16"
                    row["raw_paf_status"] = "OK"
                    row["note"] = "Raw PAF produced by explicit --fastga-frequency 16 many:many scaffold-jump 0 run."
                    rows.append(row)
        else:
            rows.append({
                "candidate_id": window["candidate_id"],
                "comparison_id": window["comparison_id"],
                "layer": "raw_frequency16",
                "paf_path": "none",
                "query_name": window["query_name"],
                "query_start": str(window["query_start"]),
                "query_end": str(window["query_end"]),
                "target_chrom": "chr3",
                "overlap_rows": "0",
                "sum_overlap_bp": "0",
                "query_union_bp": "0",
                "emits_chr3_target": "no",
                "raw_paf_status": "MISSING_RAW_PAF",
                "note": "No raw PAF found for this comparison.",
            })

    pathological = read_existing(os.path.join(PACKAGE_DIR, "summaries", "pathological_runtime.tsv"))
    if pathological:
        by_comparison = {row["comparison_id"]: row for row in pathological}
        for row in rows:
            p = by_comparison.get(row["comparison_id"])
            if not p:
                continue
            row["raw_paf_status"] = "PATHOLOGICAL_NO_RAW_PAF"
            row["note"] = (
                "No raw PAF was emitted before the frequency-16 Slurm job was cancelled; "
                "FastGA remained active with zero/near-zero .1aln output. "
                + p.get("note", "")
            ).strip()

    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "raw_chr3_support.tsv"), rows, FIELDS)


if __name__ == "__main__":
    main()
