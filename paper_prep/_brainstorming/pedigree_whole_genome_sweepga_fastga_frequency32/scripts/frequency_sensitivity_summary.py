#!/usr/bin/env python3
import csv
import os
from collections import defaultdict

from common import PACKAGE_DIR, write_tsv


PRIOR_SWEEPGA = os.path.abspath(os.path.join(
    PACKAGE_DIR, "..", "pedigree_whole_genome_sweepga_updated_bin", "summaries", "candidate_window_support.tsv"
))
WFMASH = os.path.abspath(os.path.join(
    PACKAGE_DIR, "..", "pedigree_whole_genome_wfmash_p95_updated_bin", "summaries", "candidate_window_support.tsv"
))
CURRENT = os.path.join(PACKAGE_DIR, "summaries", "raw_chr3_support.tsv")


def read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def yes_no(rows, target_key, target_value, support_key):
    return "yes" if any(row.get(target_key) == target_value and row.get(support_key) == "yes" for row in rows) else "no"


def main():
    current_rows = read_rows(CURRENT)
    prior_rows = [r for r in read_rows(PRIOR_SWEEPGA) if r.get("layer") == "raw_many_many_j0"]
    wfmash_rows = read_rows(WFMASH)

    by_candidate = defaultdict(lambda: {"candidate_id": "", "comparison_id": ""})
    for row in current_rows:
        key = row["candidate_id"]
        by_candidate[key]["candidate_id"] = row["candidate_id"]
        by_candidate[key]["comparison_id"] = row["comparison_id"]
        by_candidate[key]["query_name"] = row["query_name"]
        by_candidate[key]["query_start"] = row["query_start"]
        by_candidate[key]["query_end"] = row["query_end"]

    rows = []
    for candidate_id in sorted(by_candidate):
        cur = [r for r in current_rows if r["candidate_id"] == candidate_id]
        prior = [r for r in prior_rows if r["candidate_id"] == candidate_id]
        wf = [r for r in wfmash_rows if r["event_id"] == candidate_id]
        freq32_chr3 = yes_no(cur, "target_chrom", "chr3", "emits_chr3_target")
        prior_chr3 = yes_no(prior, "target_chrom", "chr3", "emits_chr3_target")
        wfmash_chr3 = yes_no(wf, "target_chrom", "chr3", "chr3_support")
        rows.append({
            **by_candidate[candidate_id],
            "explicit_fastga_frequency32_raw_chr3": freq32_chr3,
            "prior_updated_bin_raw_no_explicit_frequency_chr3": prior_chr3,
            "updated_wfmash_p95_chr3": wfmash_chr3,
            "frequency32_recovers_chr3_vs_prior": "yes" if freq32_chr3 == "yes" and prior_chr3 == "no" else "no",
            "wfmash_positive_sweepga_frequency32_negative": "yes" if wfmash_chr3 == "yes" and freq32_chr3 == "no" else "no",
        })

    write_tsv(
        os.path.join(PACKAGE_DIR, "summaries", "frequency_sensitivity_summary.tsv"),
        rows,
        [
            "candidate_id",
            "comparison_id",
            "query_name",
            "query_start",
            "query_end",
            "explicit_fastga_frequency32_raw_chr3",
            "prior_updated_bin_raw_no_explicit_frequency_chr3",
            "updated_wfmash_p95_chr3",
            "frequency32_recovers_chr3_vs_prior",
            "wfmash_positive_sweepga_frequency32_negative",
        ],
    )


if __name__ == "__main__":
    main()
