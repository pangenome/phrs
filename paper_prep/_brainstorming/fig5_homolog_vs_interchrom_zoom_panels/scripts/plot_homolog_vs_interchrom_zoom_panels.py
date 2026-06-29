#!/usr/bin/env python3
"""Render Fig5 telomeric homolog-vs-interchrom zoom panels.

The input class-winner table has one best same-chromosome row and one best
interchromosomal row per 2 kb query window. This script keeps the 2 kb query
windows where the interchromosomal winner ranks above the same-chromosome
winner, ranks distal 500 kb p/q arm windows by support, and renders the
strongest telomeric arms including PAR1 and the chr9q/chr3q candidate.
"""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


OUT_DIR = Path("paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels")
PREFIX = "fig5_homolog_vs_interchrom_zoom_panels"
LOCAL_SOURCE = Path(
    "paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs/"
    "PAN027pat_vs_PAN011_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz"
)
MIRROR_SOURCE = Path(
    "/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs/"
    "PAN027pat_vs_PAN011_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz"
)
CHROM_LENGTHS = Path(
    "paper_prep/_brainstorming/fig5_whole_genome_length_scaled_tracks/length_scaled_track_chromosomes.tsv"
)
COMPARISON_ID = "PAN027pat_vs_PAN011_joint"
WINDOW_BP = 500_000
REQUIRED_ARMS = {("chrX", "p"), ("chr9", "q")}
TOP_ARM_COUNT = 5

@dataclass(frozen=True)
class WinnerWindow:
    query_chrom: str
    arm: str
    query_start: int
    query_end: int
    target_chrom: str
    target_name: str
    target_start: int
    target_end: int
    inter_identity: float
    same_identity: float
    inter_minus_same_identity: float
    intersection: int
    raw_candidate_count: str
    class_candidate_count: str


def source_path() -> Path:
    if LOCAL_SOURCE.exists():
        return LOCAL_SOURCE
    if MIRROR_SOURCE.exists():
        return MIRROR_SOURCE
    raise FileNotFoundError(f"missing input: {LOCAL_SOURCE} or {MIRROR_SOURCE}")


def read_chrom_lengths() -> dict[str, int]:
    lengths: dict[str, int] = {}
    with CHROM_LENGTHS.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["comparison_id"] == COMPARISON_ID:
                length = int(row["query_chrom_length"])
                if length > 0:
                    lengths[row["query_chrom"]] = length
    if not lengths:
        raise RuntimeError(f"no chromosome lengths found for {COMPARISON_ID}")
    return lengths


def rank(row: dict[str, str]) -> tuple[float, int, float, float, str]:
    return (
        float(row["estimated.identity"]),
        int(float(row["intersection"])),
        float(row["dice.similarity"]),
        float(row["jaccard.similarity"]),
        row["other_seq"],
    )


def similarity_rank(row: dict[str, str]) -> tuple[float]:
    return (float(row["estimated.identity"]),)


def iter_telomeric_winners(path: Path, lengths: dict[str, int]) -> list[WinnerWindow]:
    by_query_window: dict[tuple[str, int, int], dict[str, tuple[tuple[object, ...], dict[str, str]]]] = defaultdict(dict)
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            query_chrom = row["query_chrom"]
            if query_chrom not in lengths:
                continue
            winner_class = row["winner_class"]
            if winner_class not in {"same_chrom", "interchrom"}:
                continue
            query_start = int(row["start"])
            query_end = int(row["end"])
            key = (query_chrom, query_start, query_end)
            row_rank = rank(row)
            if winner_class not in by_query_window[key] or row_rank > by_query_window[key][winner_class][0]:
                by_query_window[key][winner_class] = (row_rank, row)

    winners: list[WinnerWindow] = []
    for (query_chrom, query_start, query_end), class_rows in sorted(by_query_window.items()):
        if "same_chrom" not in class_rows or "interchrom" not in class_rows:
            continue
        inter_rank, inter = class_rows["interchrom"]
        same_rank, same = class_rows["same_chrom"]
        if similarity_rank(inter) <= similarity_rank(same):
            continue
        chrom_length = lengths[query_chrom]
        if query_start < WINDOW_BP:
            arm = "p"
        elif query_end > chrom_length - WINDOW_BP:
            arm = "q"
        else:
            continue
        inter_identity = float(inter["estimated.identity"])
        same_identity = float(same["estimated.identity"])
        winners.append(
            WinnerWindow(
                query_chrom=query_chrom,
                arm=arm,
                query_start=query_start,
                query_end=query_end,
                target_chrom=inter["other_chrom"],
                target_name=inter["other_seq"],
                target_start=_target_start(inter),
                target_end=_target_end(inter),
                inter_identity=inter_identity,
                same_identity=same_identity,
                inter_minus_same_identity=inter_identity - same_identity,
                intersection=int(float(inter["intersection"])),
                raw_candidate_count=inter.get("raw_candidate_count", ""),
                class_candidate_count=inter.get("class_candidate_count", ""),
            )
        )
    return winners


def _target_start(row: dict[str, str]) -> int:
    field = row["group.a"] if row["group.a"].startswith(row["other_seq"] + ":") else row["group.b"]
    return int(field.rsplit(":", 1)[1].split("-", 1)[0])


def _target_end(row: dict[str, str]) -> int:
    field = row["group.a"] if row["group.a"].startswith(row["other_seq"] + ":") else row["group.b"]
    return int(field.rsplit(":", 1)[1].split("-", 1)[1])


def arm_summary(winners: list[WinnerWindow]) -> list[dict[str, object]]:
    by_arm: dict[tuple[str, str], list[WinnerWindow]] = defaultdict(list)
    for winner in winners:
        by_arm[(winner.query_chrom, winner.arm)].append(winner)

    rows: list[dict[str, object]] = []
    for (query_chrom, arm), arm_rows in by_arm.items():
        targets: dict[str, int] = defaultdict(int)
        for row in arm_rows:
            targets[row.target_chrom] += row.query_end - row.query_start
        rows.append(
            {
                "query_chrom": query_chrom,
                "arm": arm,
                "window_bp": WINDOW_BP,
                "winner_windows": len(arm_rows),
                "winner_bp": sum(row.query_end - row.query_start for row in arm_rows),
                "mean_inter_minus_same_identity": sum(row.inter_minus_same_identity for row in arm_rows) / len(arm_rows),
                "max_inter_minus_same_identity": max(row.inter_minus_same_identity for row in arm_rows),
                "target_bp": ";".join(f"{chrom}:{bp}" for chrom, bp in sorted(targets.items(), key=lambda item: (-item[1], item[0]))),
            }
        )
    return sorted(rows, key=lambda row: (-int(row["winner_bp"]), row["query_chrom"], row["arm"]))


def selected_arms(summary_rows: list[dict[str, object]]) -> list[tuple[str, str]]:
    selected = [(str(row["query_chrom"]), str(row["arm"])) for row in summary_rows[:TOP_ARM_COUNT]]
    for arm in REQUIRED_ARMS:
        if arm not in selected:
            selected.append(arm)
    order = {(str(row["query_chrom"]), str(row["arm"])): i for i, row in enumerate(summary_rows)}
    return sorted(selected, key=lambda arm: order.get(arm, 999))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def serialize_segments(winners: list[WinnerWindow], arms: list[tuple[str, str]], lengths: dict[str, int], source: Path) -> list[dict[str, object]]:
    keep = set(arms)
    rows: list[dict[str, object]] = []
    for winner in winners:
        if (winner.query_chrom, winner.arm) not in keep:
            continue
        chrom_length = lengths[winner.query_chrom]
        display_start = 0 if winner.arm == "p" else chrom_length - WINDOW_BP
        display_end = WINDOW_BP if winner.arm == "p" else chrom_length
        rows.append(
            {
                "comparison_id": COMPARISON_ID,
                "basis": "10:10",
                "query_chrom": winner.query_chrom,
                "arm": winner.arm,
                "display_start": display_start,
                "display_end": display_end,
                "query_start": winner.query_start,
                "query_end": winner.query_end,
                "target_chrom": winner.target_chrom,
                "target_name": winner.target_name,
                "target_start": winner.target_start,
                "target_end": winner.target_end,
                "inter_identity": f"{winner.inter_identity:.7f}",
                "same_identity": f"{winner.same_identity:.7f}",
                "inter_minus_same_identity": f"{winner.inter_minus_same_identity:.7f}",
                "intersection": winner.intersection,
                "raw_candidate_count": winner.raw_candidate_count,
                "class_candidate_count": winner.class_candidate_count,
                "source_tsv_gz": source,
            }
        )
    return rows


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    src = source_path()
    lengths = read_chrom_lengths()
    winners = iter_telomeric_winners(src, lengths)
    summary_rows = arm_summary(winners)
    arms = selected_arms(summary_rows)
    segments = serialize_segments(winners, arms, lengths, src)

    write_tsv(
        OUT_DIR / "telomeric_interchrom_winner_arm_summary.tsv",
        summary_rows,
        [
            "query_chrom",
            "arm",
            "window_bp",
            "winner_windows",
            "winner_bp",
            "mean_inter_minus_same_identity",
            "max_inter_minus_same_identity",
            "target_bp",
        ],
    )
    write_tsv(
        OUT_DIR / "telomeric_interchrom_winner_segments.tsv",
        segments,
        [
            "comparison_id",
            "basis",
            "query_chrom",
            "arm",
            "display_start",
            "display_end",
            "query_start",
            "query_end",
            "target_chrom",
            "target_name",
            "target_start",
            "target_end",
            "inter_identity",
            "same_identity",
            "inter_minus_same_identity",
            "intersection",
            "raw_candidate_count",
            "class_candidate_count",
            "source_tsv_gz",
        ],
    )
    write_tsv(
        OUT_DIR / "artifact_manifest.tsv",
        [
            {"artifact": f"{PREFIX}.pdf", "role": "figure", "status": "GENERATED", "source": src},
            {"artifact": f"{PREFIX}.png", "role": "figure", "status": "GENERATED", "source": src},
            {"artifact": f"{PREFIX}.svg", "role": "figure", "status": "GENERATED", "source": src},
            {"artifact": "telomeric_interchrom_winner_segments.tsv", "role": "data", "status": "GENERATED", "source": src},
            {"artifact": "telomeric_interchrom_winner_arm_summary.tsv", "role": "summary", "status": "GENERATED", "source": src},
            {"artifact": "selected_telomeric_arms.tsv", "role": "manifest", "status": "GENERATED", "source": src},
        ],
        ["artifact", "role", "status", "source"],
    )
    write_tsv(
        OUT_DIR / "selected_telomeric_arms.tsv",
        [
            {
                "plot_order": i,
                "query_chrom": chrom,
                "arm": arm,
                "reason": "required_control" if (chrom, arm) == ("chrX", "p") else "required_candidate" if (chrom, arm) == ("chr9", "q") else "top_telomeric_arm",
            }
            for i, (chrom, arm) in enumerate(arms, start=1)
        ],
        ["plot_order", "query_chrom", "arm", "reason"],
    )


if __name__ == "__main__":
    main()
