#!/usr/bin/env python3
"""Adapt pre-IMPG similarity rows to the established Fig5 query-grid schema."""

from __future__ import annotations

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


PANEL_DIR = Path("paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_query_grid_panels")
SOURCE_DIR = Path("paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity")
WINDOW_CONFIG = Path("paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_query_grid_chop_filter_panels/config/panel_windows.tsv")

BASIS_OUTPUTS = {
    "1:1": SOURCE_DIR / "outputs/PAN027pat_vs_PAN011_joint.sweepga_f32.1to1.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz",
    "4:4": SOURCE_DIR / "outputs/PAN027pat_vs_PAN011_joint.sweepga_f32.4to4.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz",
    "10:10": SOURCE_DIR / "outputs/PAN027pat_vs_PAN011_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz",
}

CHR_RE = re.compile(r"(?:^|[#_/])chr([0-9]+|X|Y|M)(?:[._#:/-]|$)")


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


def chrom_from_name(name: str) -> str:
    match = CHR_RE.search(name)
    if match:
        return f"chr{match.group(1)}"
    match = re.search(r"chr([0-9]+|X|Y|M)", name)
    if match:
        return f"chr{match.group(1)}"
    return name


def overlap_bp(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def union_bp(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    total = 0
    cur_start, cur_end = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start > cur_end:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
        else:
            cur_end = max(cur_end, end)
    total += cur_end - cur_start
    return total


def split_group(group: str) -> tuple[str, int, int]:
    seq, coords = group.rsplit(":", 1)
    start, end = coords.split("-", 1)
    return seq, int(start), int(end)


def other_group(row: dict[str, str]) -> tuple[str, int, int] | None:
    query_seq = row["chrom"]
    query_start = int(row["start"])
    query_end = int(row["end"])
    a_seq, a_start, a_end = split_group(row["group.a"])
    b_seq, b_start, b_end = split_group(row["group.b"])
    if (a_seq, a_start, a_end) == (query_seq, query_start, query_end):
        return b_seq, b_start, b_end
    if (b_seq, b_start, b_end) == (query_seq, query_start, query_end):
        return a_seq, a_start, a_end
    return None


def target_bucket(chrom: str) -> str:
    if chrom in {"chr3", "chr9", "chrX", "chrY"}:
        return chrom
    return "other"


def rank(row: dict[str, str]) -> tuple[float, int, float, float, str]:
    return (
        float(row["estimated.identity"]),
        int(float(row["intersection"])),
        float(row["dice.similarity"]),
        float(row["jaccard.similarity"]),
        row["group.b"],
    )


def similarity_class(row: dict[str, str], query_name: str, target_name: str) -> str:
    winner_class = row.get("winner_class", "")
    if winner_class in {"same_chrom", "interchrom"}:
        return winner_class
    return "same_chrom" if chrom_from_name(query_name) == chrom_from_name(target_name) else "interchrom"


def keep_windows() -> list[dict[str, str]]:
    rows = []
    for row in read_tsv(WINDOW_CONFIG):
        if row["comparison_id"] == "PAN027pat_vs_PAN011_joint":
            rows.append(row)
    return rows


def extract_segments(windows: list[dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    segments: list[dict[str, object]] = []
    manifest: list[dict[str, object]] = []
    windows_by_name = defaultdict(list)
    for window in windows:
        windows_by_name[window["query_name"]].append(window)

    for basis, path in BASIS_OUTPUTS.items():
        if not path.exists():
            raise FileNotFoundError(path)
        seen: set[tuple[object, ...]] = set()
        best_by_query_window: dict[tuple[object, ...], dict[str, tuple[tuple[object, ...], dict[str, str], tuple[str, int, int]]]] = defaultdict(dict)
        with gzip.open(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                query_name = row["chrom"]
                if query_name not in windows_by_name:
                    continue
                query_start = int(row["start"])
                query_end = int(row["end"])
                other = other_group(row)
                if other is None:
                    continue
                target_name, target_start, target_end = other
                target_chrom = chrom_from_name(target_name)
                query_chrom = chrom_from_name(query_name)
                winner_class = similarity_class(row, query_name, target_name)
                for window in windows_by_name[query_name]:
                    w_start = int(window["query_start"])
                    w_end = int(window["query_end"])
                    ov = overlap_bp(query_start, query_end, w_start, w_end)
                    if ov <= 0:
                        continue
                    key = (basis, window["event_id"], query_name, query_start, query_end)
                    rr = rank(row)
                    if winner_class not in best_by_query_window[key] or rr > best_by_query_window[key][winner_class][0]:
                        best_by_query_window[key][winner_class] = (rr, row, other)
        for (basis, event_id, query_name, query_start, query_end), class_hits in sorted(best_by_query_window.items()):
            same_hit = class_hits.get("same_chrom")
            inter_hit = class_hits.get("interchrom")
            if same_hit is not None and inter_hit is not None:
                if inter_hit[0] > same_hit[0]:
                    winner_class = "interchrom"
                    inter_beats_same = "yes"
                    comparison_status = "interchrom_beats_same_chrom"
                else:
                    winner_class = "same_chrom"
                    inter_beats_same = "no"
                    comparison_status = "same_chrom_beats_interchrom"
            elif inter_hit is not None:
                winner_class = "interchrom"
                inter_beats_same = "no_same_chrom"
                comparison_status = "no_same_chrom_candidate"
            elif same_hit is not None:
                winner_class = "same_chrom"
                inter_beats_same = "no_interchrom"
                comparison_status = "no_interchrom_candidate"
            else:
                continue

            _rr, row, other = class_hits[winner_class]
            target_name, target_start, target_end = other
            target_chrom = chrom_from_name(target_name)
            window = next(w for w in windows_by_name[query_name] if w["event_id"] == event_id)
            w_start = int(window["query_start"])
            w_end = int(window["query_end"])
            clip_start = max(int(query_start), w_start)
            clip_end = min(int(query_end), w_end)
            dedup = (basis, event_id, query_name, clip_start, clip_end, target_name, target_start, target_end)
            if dedup in seen:
                continue
            seen.add(dedup)
            expected = set(window["expected_target_chroms"].split(","))
            intersection = int(float(row["intersection"]))
            identity = float(row["estimated.identity"])
            same_identity = float(same_hit[1]["estimated.identity"]) if same_hit is not None else None
            inter_identity = float(inter_hit[1]["estimated.identity"]) if inter_hit is not None else None
            inter_minus_same = inter_identity - same_identity if same_identity is not None and inter_identity is not None else None
            same_target = chrom_from_name(other_group(same_hit[1])[0]) if same_hit is not None and other_group(same_hit[1]) is not None else ""
            inter_target = chrom_from_name(other_group(inter_hit[1])[0]) if inter_hit is not None and other_group(inter_hit[1]) is not None else ""
            is_expected_interchrom_win = (
                winner_class == "interchrom"
                and inter_beats_same == "yes"
                and target_chrom in expected
            )
            segments.append(
                {
                    "event_id": event_id,
                    "comparison_id": window["comparison_id"],
                    "panel_label": window["panel_label"],
                    "query_name": query_name,
                    "query_chrom": chrom_from_name(query_name),
                    "window_start": w_start,
                    "window_end": w_end,
                    "basis": basis,
                    "chunk_mode": "pre-impg-query-window",
                    "num_mappings": basis,
                    "scaffold_jump": "0",
                    "scoring": "impg-estimated-identity",
                    "filter_overlap": "na",
                    "query_start": query_start,
                    "query_end": query_end,
                    "query_clip_start": clip_start,
                    "query_clip_end": clip_end,
                    "window_overlap_bp": clip_end - clip_start,
                    "target_name": target_name,
                    "target_chrom": target_chrom,
                    "target_bucket": target_bucket(target_chrom),
                    "target_start": target_start,
                    "target_end": target_end,
                    "strand": ".",
                    "matches": round(intersection * identity),
                    "alignment_length": intersection,
                    "identity": f"{identity:.6f}",
                    "winner_class": winner_class,
                    "comparison_status": comparison_status,
                    "has_same_chrom": "yes" if same_hit is not None else "no",
                    "has_interchrom": "yes" if inter_hit is not None else "no",
                    "same_target_chrom": same_target,
                    "inter_target_chrom": inter_target,
                    "same_identity": f"{same_identity:.6f}" if same_identity is not None else "",
                    "inter_identity": f"{inter_identity:.6f}" if inter_identity is not None else "",
                    "inter_minus_same_identity": f"{inter_minus_same:.6f}" if inter_minus_same is not None else "",
                    "inter_beats_same": inter_beats_same,
                    "expected_target_chroms": window["expected_target_chroms"],
                    "is_expected_target": "yes" if target_chrom in expected else "no",
                    "is_expected_interchrom_win": "yes" if is_expected_interchrom_win else "no",
                    "source_tsv_gz": path,
                }
            )
        for window in windows:
            manifest.append(
                {
                    "event_id": window["event_id"],
                    "panel_label": window["panel_label"],
                    "comparison_id": window["comparison_id"],
                    "query_name": window["query_name"],
                    "window_start": window["query_start"],
                    "window_end": window["query_end"],
                    "basis": basis,
                    "source_tsv_gz": path,
                    "status": "OK",
                }
            )

    by_key = defaultdict(list)
    for row in segments:
        by_key[(row["event_id"], row["basis"])].append(row)

    summary: list[dict[str, object]] = []
    for window in windows:
        for basis in BASIS_OUTPUTS:
            rows = by_key[(window["event_id"], basis)]
            expected_rows = [row for row in rows if row["is_expected_target"] == "yes"]
            interchrom_winner_rows = [row for row in rows if row["winner_class"] == "interchrom"]
            inter_beats_same_rows = [row for row in rows if row["inter_beats_same"] == "yes"]
            expected_interchrom_win_rows = [row for row in rows if row["is_expected_interchrom_win"] == "yes"]
            by_target_sum = defaultdict(int)
            by_target_intervals = defaultdict(list)
            for row in rows:
                by_target_sum[row["target_chrom"]] += int(row["window_overlap_bp"])
                by_target_intervals[row["target_chrom"]].append((int(row["query_clip_start"]), int(row["query_clip_end"])))
            summary.append(
                {
                    "event_id": window["event_id"],
                    "comparison_id": window["comparison_id"],
                    "query_name": window["query_name"],
                    "query_chrom": chrom_from_name(window["query_name"]),
                    "window_start": window["query_start"],
                    "window_end": window["query_end"],
                    "basis": basis,
                    "chunk_mode": "pre-impg-query-window",
                    "num_mappings": basis,
                    "scaffold_jump": "0",
                    "scoring": "impg-estimated-identity",
                    "expected_target_chroms": window["expected_target_chroms"],
                    "segment_rows": len(rows),
                    "expected_target_rows": len(expected_rows),
                    "sum_expected_overlap_bp": sum(int(r["window_overlap_bp"]) for r in expected_rows),
                    "union_expected_overlap_bp": union_bp([(int(r["query_clip_start"]), int(r["query_clip_end"])) for r in expected_rows]),
                    "same_chrom_winner_rows": len([row for row in rows if row["winner_class"] == "same_chrom"]),
                    "interchrom_winner_rows": len(interchrom_winner_rows),
                    "interchrom_winner_union_bp": union_bp([(int(r["query_clip_start"]), int(r["query_clip_end"])) for r in interchrom_winner_rows]),
                    "inter_beats_same_rows": len(inter_beats_same_rows),
                    "inter_beats_same_union_bp": union_bp([(int(r["query_clip_start"]), int(r["query_clip_end"])) for r in inter_beats_same_rows]),
                    "expected_interchrom_win_rows": len(expected_interchrom_win_rows),
                    "expected_interchrom_win_union_bp": union_bp([(int(r["query_clip_start"]), int(r["query_clip_end"])) for r in expected_interchrom_win_rows]),
                    "target_sum_overlap_bp": ";".join(f"{k}:{v}" for k, v in sorted(by_target_sum.items())),
                    "target_union_overlap_bp": ";".join(f"{k}:{union_bp(v)}" for k, v in sorted(by_target_intervals.items())),
                    "status": "OK" if rows else "NO_ROWS",
                }
            )
    return segments, summary, manifest


def main() -> None:
    windows = keep_windows()
    segments, summary, manifest = extract_segments(windows)
    segment_fields = [
        "event_id", "comparison_id", "panel_label", "query_name", "query_chrom", "window_start", "window_end",
        "basis", "chunk_mode", "num_mappings", "scaffold_jump", "scoring", "filter_overlap",
        "query_start", "query_end", "query_clip_start", "query_clip_end", "window_overlap_bp",
        "target_name", "target_chrom", "target_bucket", "target_start", "target_end", "strand",
        "matches", "alignment_length", "identity", "winner_class", "comparison_status",
        "has_same_chrom", "has_interchrom", "same_target_chrom", "inter_target_chrom",
        "same_identity", "inter_identity", "inter_minus_same_identity", "inter_beats_same",
        "expected_target_chroms", "is_expected_target", "is_expected_interchrom_win", "source_tsv_gz",
    ]
    summary_fields = [
        "event_id", "comparison_id", "query_name", "query_chrom", "window_start", "window_end",
        "basis", "chunk_mode", "num_mappings", "scaffold_jump", "scoring", "expected_target_chroms",
        "segment_rows", "expected_target_rows", "sum_expected_overlap_bp", "union_expected_overlap_bp",
        "same_chrom_winner_rows", "interchrom_winner_rows", "interchrom_winner_union_bp",
        "inter_beats_same_rows", "inter_beats_same_union_bp",
        "expected_interchrom_win_rows", "expected_interchrom_win_union_bp",
        "target_sum_overlap_bp", "target_union_overlap_bp", "status",
    ]
    manifest_fields = ["event_id", "panel_label", "comparison_id", "query_name", "window_start", "window_end", "basis", "source_tsv_gz", "status"]
    write_tsv(PANEL_DIR / "pre_impg_query_grid_panel_segments.tsv", segments, segment_fields)
    write_tsv(PANEL_DIR / "pre_impg_query_grid_panel_summary.tsv", summary, summary_fields)
    write_tsv(PANEL_DIR / "pre_impg_query_grid_panel_manifest.tsv", manifest, manifest_fields)


if __name__ == "__main__":
    main()
