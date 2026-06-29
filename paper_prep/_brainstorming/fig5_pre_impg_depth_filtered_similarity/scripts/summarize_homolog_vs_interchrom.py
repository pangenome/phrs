#!/usr/bin/env python3
"""Summarize whole-genome homolog-vs-interchrom IMPG class winners."""

from __future__ import annotations

import csv
import gzip
import statistics
from collections import defaultdict
from pathlib import Path


BASE_DIR = Path("paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity")
BASIS_OUTPUTS = {
    "1:1": BASE_DIR / "outputs/PAN027pat_vs_PAN011_joint.sweepga_f32.1to1.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz",
    "4:4": BASE_DIR / "outputs/PAN027pat_vs_PAN011_joint.sweepga_f32.4to4.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz",
    "10:10": BASE_DIR / "outputs/PAN027pat_vs_PAN011_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz",
}
SUMMARY_DIR = BASE_DIR / "summaries"
WINDOW_BP = 2000


def parse_class_winners(path: Path) -> dict[tuple[str, int, int], dict[str, dict[str, str]]]:
    groups: dict[tuple[str, int, int], dict[str, dict[str, str]]] = {}
    with gzip.open(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["chrom"], int(row["start"]), int(row["end"]))
            groups.setdefault(key, {})[row["winner_class"]] = row
    return groups


def merge_intervals(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_pair: dict[tuple[str, str], list[tuple[int, int, dict[str, object]]]] = defaultdict(list)
    for row in rows:
        by_pair[(str(row["query_chrom"]), str(row["other_chrom"]))].append((int(row["start"]), int(row["end"]), row))

    tracts: list[dict[str, object]] = []
    for (query_chrom, target_chrom), intervals in by_pair.items():
        intervals.sort()
        cur_start: int | None = None
        cur_end: int | None = None
        cur_rows: list[dict[str, object]] = []
        for start, end, row in intervals:
            if cur_start is None:
                cur_start, cur_end, cur_rows = start, end, [row]
            elif start <= int(cur_end):
                cur_end = max(int(cur_end), end)
                cur_rows.append(row)
            else:
                tracts.append(make_tract(query_chrom, target_chrom, int(cur_start), int(cur_end), cur_rows))
                cur_start, cur_end, cur_rows = start, end, [row]
        if cur_start is not None:
            tracts.append(make_tract(query_chrom, target_chrom, int(cur_start), int(cur_end), cur_rows))
    return tracts


def make_tract(query_chrom: str, target_chrom: str, start: int, end: int, rows: list[dict[str, object]]) -> dict[str, object]:
    deltas = [float(row["delta"]) for row in rows]
    return {
        "query_chrom": query_chrom,
        "target_chrom": target_chrom,
        "query_start": start,
        "query_end": end,
        "bp": end - start,
        "windows": len(rows),
        "mean_delta": sum(deltas) / len(deltas),
        "max_delta": max(deltas),
        "mean_inter_identity": sum(float(row["inter_identity"]) for row in rows) / len(rows),
        "mean_same_identity": sum(float(row["same_identity"]) for row in rows) / len(rows),
    }


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = {}
            for field in fields:
                value = row.get(field, "")
                if isinstance(value, float):
                    value = f"{value:.6f}"
                out[field] = value
            writer.writerow(out)


def main() -> None:
    overall_rows: list[dict[str, object]] = []
    pair_rows: list[dict[str, object]] = []
    tract_rows: list[dict[str, object]] = []

    for basis, path in BASIS_OUTPUTS.items():
        if not path.exists():
            raise FileNotFoundError(path)
        groups = parse_class_winners(path)
        counts: dict[str, int] = defaultdict(int)
        inter_win_rows: list[dict[str, object]] = []

        for (_seq, start, end), rows in groups.items():
            same = rows.get("same_chrom")
            inter = rows.get("interchrom")
            if same is not None and inter is not None:
                same_identity = float(same["estimated.identity"])
                inter_identity = float(inter["estimated.identity"])
                if inter_identity > same_identity:
                    counts["inter_beats_same"] += 1
                    out = dict(inter)
                    out["delta"] = inter_identity - same_identity
                    out["same_identity"] = same_identity
                    out["inter_identity"] = inter_identity
                    out["start"] = start
                    out["end"] = end
                    inter_win_rows.append(out)
                elif same_identity > inter_identity:
                    counts["same_beats_inter"] += 1
                else:
                    counts["tie"] += 1
                counts["both"] += 1
            elif inter is not None:
                counts["inter_only"] += 1
            elif same is not None:
                counts["same_only"] += 1

        overall_rows.append(
            {
                "basis": basis,
                "windows_total": len(groups),
                "both_same_and_inter": counts["both"],
                "inter_beats_same_windows": counts["inter_beats_same"],
                "inter_beats_same_bp": counts["inter_beats_same"] * WINDOW_BP,
                "same_beats_inter_windows": counts["same_beats_inter"],
                "tie_windows": counts["tie"],
                "inter_only_windows": counts["inter_only"],
                "same_only_windows": counts["same_only"],
            }
        )

        by_pair: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
        for row in inter_win_rows:
            by_pair[(str(row["query_chrom"]), str(row["other_chrom"]))].append(row)

        for (query_chrom, target_chrom), rows in by_pair.items():
            deltas = [float(row["delta"]) for row in rows]
            pair_rows.append(
                {
                    "basis": basis,
                    "query_chrom": query_chrom,
                    "target_chrom": target_chrom,
                    "windows": len(rows),
                    "bp": len(rows) * WINDOW_BP,
                    "mean_delta": sum(deltas) / len(deltas),
                    "median_delta": statistics.median(deltas),
                    "max_delta": max(deltas),
                    "mean_inter_identity": sum(float(row["inter_identity"]) for row in rows) / len(rows),
                    "mean_same_identity": sum(float(row["same_identity"]) for row in rows) / len(rows),
                }
            )

        for tract in merge_intervals(inter_win_rows):
            tract["basis"] = basis
            tract_rows.append(tract)

    pair_rows.sort(key=lambda row: (row["basis"] != "10:10", -int(row["bp"]), row["query_chrom"], row["target_chrom"]))
    tract_rows.sort(key=lambda row: (row["basis"] != "10:10", -int(row["bp"]), row["query_chrom"], int(row["query_start"]), row["target_chrom"]))

    write_tsv(
        SUMMARY_DIR / "homolog_vs_interchrom_overall.tsv",
        overall_rows,
        [
            "basis",
            "windows_total",
            "both_same_and_inter",
            "inter_beats_same_windows",
            "inter_beats_same_bp",
            "same_beats_inter_windows",
            "tie_windows",
            "inter_only_windows",
            "same_only_windows",
        ],
    )
    write_tsv(
        SUMMARY_DIR / "homolog_vs_interchrom_pair_summary.tsv",
        pair_rows,
        [
            "basis",
            "query_chrom",
            "target_chrom",
            "windows",
            "bp",
            "mean_delta",
            "median_delta",
            "max_delta",
            "mean_inter_identity",
            "mean_same_identity",
        ],
    )
    write_tsv(
        SUMMARY_DIR / "homolog_vs_interchrom_top_tracts.tsv",
        tract_rows,
        [
            "basis",
            "query_chrom",
            "target_chrom",
            "query_start",
            "query_end",
            "bp",
            "windows",
            "mean_delta",
            "max_delta",
            "mean_inter_identity",
            "mean_same_identity",
        ],
    )


if __name__ == "__main__":
    main()
