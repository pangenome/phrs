#!/usr/bin/env python3
"""Extract telomeric interchrom-over-homolog windows for Fig5 zoom panels."""

from __future__ import annotations

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


OUT_DIR = Path("paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels")
CLASS_WINNERS = Path(
    "paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs/"
    "PAN027pat_vs_PAN011_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz"
)
QUERY_FAI = Path(
    "/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/"
    "pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027pat_vs_PAN011_joint.query.fa.fai"
)
ZOOM_BP = 500_000
WINDOW_BP = 2_000
PANEL_ORDER = [
    ("chrX", "p", "PAR1 X/Y control"),
    ("chr22", "p", "chr22p acrocentric"),
    ("chr5", "q", "chr5q subtelomeric"),
    ("chr15", "p", "chr15p acrocentric"),
    ("chr9", "q", "PAN027 chr9q -> chr3q"),
]

CHR_RE = re.compile(r"(chr(?:[0-9]+|X|Y|M))")


def chrom_name(value: str) -> str:
    match = CHR_RE.search(value)
    return match.group(1) if match else value


def read_fai(path: Path) -> dict[str, tuple[str, int]]:
    out: dict[str, tuple[str, int]] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            seq, length, *_ = line.rstrip("\n").split("\t")
            out[chrom_name(seq)] = (seq, int(length))
    return out


def relative_interval(chrom: str, arm: str, start: int, end: int, lengths: dict[str, tuple[str, int]]) -> tuple[int, int] | None:
    _seq, length = lengths[chrom]
    if arm == "p":
        if start >= ZOOM_BP:
            return None
        return max(0, start), min(ZOOM_BP, end)
    if end <= length - ZOOM_BP:
        return None
    return max(0, length - end), min(ZOOM_BP, length - start)


def target_bucket(target_chrom: str) -> str:
    if target_chrom in {"chr3", "chrY", "chr1", "chr13", "chr14", "chr15", "chr21", "chr22"}:
        return target_chrom
    return "other"


def main() -> None:
    lengths = read_fai(QUERY_FAI)
    panel_by_chrom_arm = {(chrom, arm): (idx + 1, label) for idx, (chrom, arm, label) in enumerate(PANEL_ORDER)}
    groups: dict[tuple[str, int, int], dict[str, dict[str, str]]] = defaultdict(dict)

    with gzip.open(CLASS_WINNERS, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["chrom"], int(row["start"]), int(row["end"]))
            groups[key][row["winner_class"]] = row

    segments: list[dict[str, object]] = []
    for (query_name, start, end), rows in groups.items():
        if "same_chrom" not in rows or "interchrom" not in rows:
            continue
        same = rows["same_chrom"]
        inter = rows["interchrom"]
        same_identity = float(same["estimated.identity"])
        inter_identity = float(inter["estimated.identity"])
        if inter_identity <= same_identity:
            continue
        chrom = chrom_name(query_name)
        for arm in ("p", "q"):
            panel = panel_by_chrom_arm.get((chrom, arm))
            if panel is None:
                continue
            rel = relative_interval(chrom, arm, start, end, lengths)
            if rel is None:
                continue
            rel_start, rel_end = rel
            if rel_end <= rel_start:
                continue
            panel_order, panel_label = panel
            target_chrom = inter["other_chrom"]
            segments.append(
                {
                    "panel_order": panel_order,
                    "panel_id": f"{chrom}_{arm}",
                    "panel_label": panel_label,
                    "query_name": query_name,
                    "query_chrom": chrom,
                    "arm": arm,
                    "query_length": lengths[chrom][1],
                    "zoom_bp": ZOOM_BP,
                    "query_start": start,
                    "query_end": end,
                    "relative_start": rel_start,
                    "relative_end": rel_end,
                    "window_overlap_bp": rel_end - rel_start,
                    "target_chrom": target_chrom,
                    "target_bucket": target_bucket(target_chrom),
                    "target_name": inter["other_seq"],
                    "same_identity": f"{same_identity:.6f}",
                    "inter_identity": f"{inter_identity:.6f}",
                    "delta_identity": f"{inter_identity - same_identity:.6f}",
                    "intersection": inter["intersection"],
                    "inter_group_a": inter["group.a"],
                    "inter_group_b": inter["group.b"],
                }
            )

    summaries: list[dict[str, object]] = []
    by_panel: dict[str, list[dict[str, object]]] = defaultdict(list)
    by_panel_target: dict[tuple[str, str], int] = defaultdict(int)
    for row in segments:
        by_panel[str(row["panel_id"])].append(row)
        by_panel_target[(str(row["panel_id"]), str(row["target_chrom"]))] += int(row["window_overlap_bp"])

    for panel_order, (chrom, arm, label) in enumerate(PANEL_ORDER, start=1):
        panel_id = f"{chrom}_{arm}"
        rows = by_panel.get(panel_id, [])
        top_targets = [
            f"{target}:{bp}"
            for (pid, target), bp in sorted(by_panel_target.items(), key=lambda kv: (kv[0][0] != panel_id, -kv[1], kv[0][1]))
            if pid == panel_id
        ]
        summaries.append(
            {
                "panel_order": panel_order,
                "panel_id": panel_id,
                "panel_label": label,
                "query_chrom": chrom,
                "arm": arm,
                "query_length": lengths[chrom][1],
                "zoom_bp": ZOOM_BP,
                "winning_windows": len(rows),
                "winning_bp": sum(int(row["window_overlap_bp"]) for row in rows),
                "top_targets": ";".join(top_targets),
            }
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    segment_fields = [
        "panel_order",
        "panel_id",
        "panel_label",
        "query_name",
        "query_chrom",
        "arm",
        "query_length",
        "zoom_bp",
        "query_start",
        "query_end",
        "relative_start",
        "relative_end",
        "window_overlap_bp",
        "target_chrom",
        "target_bucket",
        "target_name",
        "same_identity",
        "inter_identity",
        "delta_identity",
        "intersection",
        "inter_group_a",
        "inter_group_b",
    ]
    summary_fields = [
        "panel_order",
        "panel_id",
        "panel_label",
        "query_chrom",
        "arm",
        "query_length",
        "zoom_bp",
        "winning_windows",
        "winning_bp",
        "top_targets",
    ]
    with (OUT_DIR / "zoom_window_segments.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=segment_fields, lineterminator="\n")
        writer.writeheader()
        for row in sorted(segments, key=lambda r: (int(r["panel_order"]), int(r["relative_start"]), str(r["target_chrom"]))):
            writer.writerow(row)
    with (OUT_DIR / "zoom_panel_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=summary_fields, lineterminator="\n")
        writer.writeheader()
        for row in summaries:
            writer.writerow(row)


if __name__ == "__main__":
    main()

