#!/usr/bin/env python3
"""Summarize Fig5 candidate/control support in 2 kb query-space bins.

This is intentionally PAF-only: it reuses existing SweepGA and WFMASH evidence
and avoids seqwish/ODGI graph construction.
"""

from __future__ import annotations

import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path


PILOT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = PILOT_DIR.parents[2]
WINDOWS_TSV = PILOT_DIR / "config" / "evaluation_windows.tsv"
SUMMARY_DIR = PILOT_DIR / "summaries"

SWEEPGA_DIR = REPO_ROOT / "paper_prep/_brainstorming/pedigree_direct_sweepga_joint_parent"
WFMASH_DIR = REPO_ROOT / "paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin"

BIN_SIZE = 2000
SWEEPGA_FILTERS = [
    ("sweepga_many_many_noscaffold", "filtered_paf/{comparison_id}.many_many_noscaffold.paf.gz"),
    ("sweepga_four_many_noscaffold", "filtered_paf/{comparison_id}.four_many_noscaffold.paf.gz"),
    ("sweepga_one_one_noscaffold", "filtered_paf/{comparison_id}.one_one_noscaffold.paf.gz"),
    ("sweepga_simple_i95_l1k_q80", "filtered_paf/{comparison_id}.simple_i95_l1k_q80.paf.gz"),
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def parse_chrom(name: str) -> str:
    match = re.search(r"#chr([0-9XYM]+)", name)
    if match:
        return f"chr{match.group(1)}"
    match = re.search(r"chr([0-9XYM]+)", name)
    if match:
        return f"chr{match.group(1)}"
    return name


def parse_arm(name: str) -> str:
    chrom = parse_chrom(name)
    if name.endswith("_parm"):
        return f"{chrom}p"
    if name.endswith("_qarm"):
        return f"{chrom}q"
    return chrom


def target_class(target_chrom: str, expected: str, control_class: str, query_chrom: str) -> str:
    acros = {"chr13", "chr14", "chr15", "chr21", "chr22"}
    if target_chrom == expected:
        return "expected"
    if target_chrom == query_chrom:
        return "same_chromosome"
    if control_class == "PAR_positive" and target_chrom in {"chrX", "chrY"}:
        return "PAR_partner"
    if control_class == "acrocentric_p_control" and target_chrom in acros:
        return "cross_acrocentric"
    return "other"


def identity_from_paf(fields: list[str]) -> float:
    try:
        matches = float(fields[9])
        aln_len = float(fields[10])
        if aln_len > 0:
            return matches / aln_len
    except (IndexError, ValueError):
        pass
    for field in fields[12:]:
        if field.startswith("dv:f:"):
            return max(0.0, min(1.0, 1.0 - float(field.split(":")[-1])))
    return math.nan


def summarize_paf_for_windows(method: str, paf_path: Path, windows: list[dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    by_query = defaultdict(list)
    for window in windows:
        if window["comparison_id"] in paf_path.name:
            by_query[window["query_name"]].append(window)

    bin_acc = defaultdict(lambda: {"rows": 0, "aligned_bp_sum": 0, "identity_bp_sum": 0.0})
    with open_text(paf_path) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            query_name = fields[0]
            if query_name not in by_query:
                continue
            q_start = int(fields[2])
            q_end = int(fields[3])
            target_chrom = parse_chrom(fields[5])
            target_arm = parse_arm(fields[5])
            ident = identity_from_paf(fields)
            for window in by_query[query_name]:
                w_start = int(window["window_start"])
                w_end = int(window["window_end"])
                if overlap(q_start, q_end, w_start, w_end) == 0:
                    continue
                query_chrom = parse_chrom(query_name)
                first_bin = max(0, (max(q_start, w_start) - w_start) // BIN_SIZE)
                last_bin = max(0, (min(q_end, w_end) - 1 - w_start) // BIN_SIZE)
                for bin_index in range(first_bin, last_bin + 1):
                    b_start = w_start + bin_index * BIN_SIZE
                    b_end = min(w_end, b_start + BIN_SIZE)
                    ov = overlap(q_start, q_end, b_start, b_end)
                    if ov <= 0:
                        continue
                    key = (
                        method,
                        window["event_id"],
                        window["comparison_id"],
                        query_name,
                        bin_index,
                        b_start,
                        b_end,
                        target_chrom,
                        target_arm,
                        target_class(target_chrom, window["expected_target_chrom"], window["control_class"], query_chrom),
                    )
                    acc = bin_acc[key]
                    acc["rows"] += 1
                    acc["aligned_bp_sum"] += ov
                    if not math.isnan(ident):
                        acc["identity_bp_sum"] += ident * ov

    bin_rows = []
    for key, acc in sorted(bin_acc.items()):
        (
            method_id,
            event_id,
            comparison_id,
            query_name,
            bin_index,
            bin_start,
            bin_end,
            target_chrom,
            target_arm,
            cls,
        ) = key
        aligned = acc["aligned_bp_sum"]
        mean_identity = acc["identity_bp_sum"] / aligned if aligned else math.nan
        bin_rows.append({
            "method_id": method_id,
            "event_id": event_id,
            "comparison_id": comparison_id,
            "query_name": query_name,
            "bin_size_bp": BIN_SIZE,
            "bin_index": bin_index,
            "bin_start": bin_start,
            "bin_end": bin_end,
            "target_chrom": target_chrom,
            "target_arm": target_arm,
            "target_class": cls,
            "paf_rows": acc["rows"],
            "aligned_bp_sum": aligned,
            "bin_fraction_aligned": round(aligned / max(1, bin_end - bin_start), 6),
            "mean_identity_weighted": round(mean_identity, 6) if not math.isnan(mean_identity) else "",
            "mean_match_distance": round(1.0 - mean_identity, 6) if not math.isnan(mean_identity) else "",
            "paf_path": str(paf_path),
        })

    summary_acc = defaultdict(lambda: {"bins": set(), "aligned_bp_sum": 0, "identity_bp_sum": 0.0, "rows": 0})
    for row in bin_rows:
        key = (
            row["method_id"],
            row["event_id"],
            row["comparison_id"],
            row["target_chrom"],
            row["target_class"],
        )
        acc = summary_acc[key]
        acc["bins"].add(row["bin_index"])
        aligned = int(row["aligned_bp_sum"])
        acc["aligned_bp_sum"] += aligned
        acc["rows"] += int(row["paf_rows"])
        if row["mean_identity_weighted"] != "":
            acc["identity_bp_sum"] += float(row["mean_identity_weighted"]) * aligned

    summary_rows = []
    for key, acc in sorted(summary_acc.items()):
        method_id, event_id, comparison_id, target_chrom, cls = key
        aligned = acc["aligned_bp_sum"]
        mean_identity = acc["identity_bp_sum"] / aligned if aligned else math.nan
        summary_rows.append({
            "method_id": method_id,
            "event_id": event_id,
            "comparison_id": comparison_id,
            "target_chrom": target_chrom,
            "target_class": cls,
            "bins_with_support": len(acc["bins"]),
            "aligned_bp_sum_across_bins": aligned,
            "paf_rows_across_bins": acc["rows"],
            "mean_identity_weighted": round(mean_identity, 6) if not math.isnan(mean_identity) else "",
            "mean_match_distance": round(1.0 - mean_identity, 6) if not math.isnan(mean_identity) else "",
        })
    present = {(row["event_id"], row["comparison_id"]) for row in summary_rows}
    for window in windows:
        if window["comparison_id"] not in paf_path.name:
            continue
        key = (window["event_id"], window["comparison_id"])
        if key in present:
            continue
        summary_rows.append({
            "method_id": method,
            "event_id": window["event_id"],
            "comparison_id": window["comparison_id"],
            "target_chrom": "NO_SUPPORT",
            "target_class": "none",
            "bins_with_support": 0,
            "aligned_bp_sum_across_bins": 0,
            "paf_rows_across_bins": 0,
            "mean_identity_weighted": "",
            "mean_match_distance": "",
        })
    return bin_rows, summary_rows


def load_wfmash_2kb_rows() -> list[dict[str, object]]:
    path = WFMASH_DIR / "summaries/query_grid_filter_candidate_window_support.tsv"
    rows = []
    if not path.exists():
        return rows
    for row in read_tsv(path):
        if row.get("chop_length_bp") != "2000":
            continue
        rows.append({
            "method_id": "wfmash_p95_one_to_one_2kb",
            "event_id": row["event_id"],
            "comparison_id": row["comparison_id"],
            "target_chrom": row["target_chrom"],
            "target_class": "expected" if row["target_chrom"] == row["expected_target_chrom"] else "same_chromosome" if row["target_chrom"] == "chr9" else "other",
            "bins_with_support": row["retained_rows_overlapping_window"],
            "aligned_bp_sum_across_bins": row["query_overlap_bp_sum"],
            "paf_rows_across_bins": row["retained_rows_overlapping_window"],
            "mean_identity_weighted": "",
            "mean_match_distance": "",
        })
    return rows


def make_method_comparison(window_summary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_event_method = defaultdict(list)
    for row in window_summary_rows:
        by_event_method[(row["event_id"], row["comparison_id"], row["method_id"])].append(row)

    rows = []
    for (event_id, comparison_id, method_id), group in sorted(by_event_method.items()):
        expected = [r for r in group if r["target_class"] == "expected"]
        same = [r for r in group if r["target_class"] == "same_chromosome"]
        cross_acro = [r for r in group if r["target_class"] == "cross_acrocentric"]
        par_partner = [r for r in group if r["target_class"] == "PAR_partner"]
        top = max(group, key=lambda r: int(r["aligned_bp_sum_across_bins"]))
        expected_bp = sum(int(r["aligned_bp_sum_across_bins"]) for r in expected)
        same_bp = sum(int(r["aligned_bp_sum_across_bins"]) for r in same)
        rows.append({
            "event_id": event_id,
            "comparison_id": comparison_id,
            "method_id": method_id,
            "top_target_chrom": top["target_chrom"],
            "top_target_class": top["target_class"],
            "top_aligned_bp_sum_across_bins": top["aligned_bp_sum_across_bins"],
            "expected_target_aligned_bp_sum": expected_bp,
            "expected_target_bins_with_support": sum(int(r["bins_with_support"]) for r in expected),
            "same_chromosome_aligned_bp_sum": same_bp,
            "expected_vs_same_bp_ratio": round(expected_bp / same_bp, 6) if same_bp else "",
            "par_partner_aligned_bp_sum": sum(int(r["aligned_bp_sum_across_bins"]) for r in par_partner),
            "cross_acrocentric_aligned_bp_sum": sum(int(r["aligned_bp_sum_across_bins"]) for r in cross_acro),
            "supports_expected_target": "yes" if expected_bp > 0 else "no",
        })
    return rows


def main() -> None:
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    windows = read_tsv(WINDOWS_TSV)
    all_bin_rows = []
    all_summary_rows = []

    comparisons = sorted({row["comparison_id"] for row in windows})
    for comparison_id in comparisons:
        for method, template in SWEEPGA_FILTERS:
            paf_path = SWEEPGA_DIR / template.format(comparison_id=comparison_id)
            if not paf_path.exists():
                raise FileNotFoundError(paf_path)
            bin_rows, summary_rows = summarize_paf_for_windows(method, paf_path, windows)
            all_bin_rows.extend(bin_rows)
            all_summary_rows.extend(summary_rows)

    wfmash_rows = load_wfmash_2kb_rows()
    all_summary_rows.extend(wfmash_rows)
    comparison_rows = make_method_comparison(all_summary_rows)

    write_tsv(SUMMARY_DIR / "bin_target_support.tsv", all_bin_rows, [
        "method_id",
        "event_id",
        "comparison_id",
        "query_name",
        "bin_size_bp",
        "bin_index",
        "bin_start",
        "bin_end",
        "target_chrom",
        "target_arm",
        "target_class",
        "paf_rows",
        "aligned_bp_sum",
        "bin_fraction_aligned",
        "mean_identity_weighted",
        "mean_match_distance",
        "paf_path",
    ])
    write_tsv(SUMMARY_DIR / "window_target_summary.tsv", all_summary_rows, [
        "method_id",
        "event_id",
        "comparison_id",
        "target_chrom",
        "target_class",
        "bins_with_support",
        "aligned_bp_sum_across_bins",
        "paf_rows_across_bins",
        "mean_identity_weighted",
        "mean_match_distance",
    ])
    write_tsv(SUMMARY_DIR / "method_comparison.tsv", comparison_rows, [
        "event_id",
        "comparison_id",
        "method_id",
        "top_target_chrom",
        "top_target_class",
        "top_aligned_bp_sum_across_bins",
        "expected_target_aligned_bp_sum",
        "expected_target_bins_with_support",
        "same_chromosome_aligned_bp_sum",
        "expected_vs_same_bp_ratio",
        "par_partner_aligned_bp_sum",
        "cross_acrocentric_aligned_bp_sum",
        "supports_expected_target",
    ])


if __name__ == "__main__":
    main()
