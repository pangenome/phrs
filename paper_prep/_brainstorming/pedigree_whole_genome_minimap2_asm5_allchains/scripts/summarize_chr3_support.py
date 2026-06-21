#!/usr/bin/env python3
import csv
import os

from common import PACKAGE_DIR, write_tsv


def read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    support_path = os.path.join(PACKAGE_DIR, "summaries", "candidate_window_support.tsv")
    rows = read_rows(support_path)
    by_event = {}
    for row in rows:
        event = row["event_id"]
        entry = by_event.setdefault(event, {
            "event_id": event,
            "comparison_id": row["comparison_id"],
            "query_name": row["query_name"],
            "query_start": row["query_start"],
            "query_end": row["query_end"],
            "expected_target_chrom": row["expected_target_chrom"],
            "minimap2_chr3_support": "no",
            "chr3_rows": "0",
            "chr3_query_overlap_bp_sum": "0",
            "chr3_query_covered_bp_union": "0",
            "all_target_chromosomes_overlapping_window": [],
            "answer": "no",
            "comparison_to_other_aligners": "",
        })
        chrom = row["target_chrom"]
        if chrom not in entry["all_target_chromosomes_overlapping_window"]:
            entry["all_target_chromosomes_overlapping_window"].append(chrom)
        if chrom == "chr3" and int(row["paf_rows_overlapping_window"]) > 0:
            entry["minimap2_chr3_support"] = "yes"
            entry["answer"] = "yes"
            entry["chr3_rows"] = row["paf_rows_overlapping_window"]
            entry["chr3_query_overlap_bp_sum"] = row["query_overlap_bp_sum"]
            entry["chr3_query_covered_bp_union"] = row["query_covered_bp_union"]

    out_rows = []
    for event in sorted(by_event):
        entry = by_event[event]
        entry["all_target_chromosomes_overlapping_window"] = ",".join(sorted(entry["all_target_chromosomes_overlapping_window"]))
        entry["comparison_to_other_aligners"] = (
            "Updated wfmash -p 95 was chr3-positive for the PAN027/PAN028 Fig5 windows; "
            "sweepGA/FastGA default raw PAF was chr3-negative; this minimap2 asm5 all-chain row is %s."
        ) % ("chr3-positive" if entry["answer"] == "yes" else "chr3-negative")
        out_rows.append(entry)

    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "minimap2_chr3_support_summary.tsv"), out_rows, [
        "event_id",
        "comparison_id",
        "query_name",
        "query_start",
        "query_end",
        "expected_target_chrom",
        "minimap2_chr3_support",
        "chr3_rows",
        "chr3_query_overlap_bp_sum",
        "chr3_query_covered_bp_union",
        "all_target_chromosomes_overlapping_window",
        "answer",
        "comparison_to_other_aligners",
    ])


if __name__ == "__main__":
    main()
