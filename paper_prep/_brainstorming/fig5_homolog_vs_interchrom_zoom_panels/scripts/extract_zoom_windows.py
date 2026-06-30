#!/usr/bin/env python3
"""Extract telomeric interchrom-over-homolog windows for Fig5 zoom panels."""

from __future__ import annotations

import csv
import gzip
import argparse
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
PHR_TABLE = Path(
    "/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/"
    "all-vs-all.1Mb.p95.id95.len.tsv"
)
ZOOM_BP = 500_000
WINDOW_BP = 2_000
PANEL_ORDER = [
    ("chrX", "p", "chrXp PAR1"),
    ("chr5", "q", "chr5q"),
    ("chr9", "q", "chr9q"),
]

CHR_RE = re.compile(r"(chr(?:[0-9]+|X|Y|M))")
TARGET_HAP_RE = re.compile(r"#(h[0-9]+)_")
PHR_SEQ_RE = re.compile(
    r"^(?P<query_prefix>[^:]+):(?P<seq_start>[0-9]+)-(?P<seq_end>[0-9]+)_"
    r"(?P<label_chrom>chr(?:[0-9]+|X|Y|M))_(?P<arm>[pq])arm$"
)


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


def target_haplotype(target_name: str) -> str:
    match = TARGET_HAP_RE.search(target_name)
    return match.group(1) if match else "NA"


def project_phr_interval(
    row: dict[str, str],
    lengths: dict[str, tuple[str, int]],
    source_path: Path,
) -> dict[str, object] | None:
    match = PHR_SEQ_RE.match(row["seq"])
    if match is None:
        return None

    query_name = match.group("query_prefix").split(".")[0]
    if not query_name.startswith("PAN027#2#"):
        return None

    chrom = chrom_name(query_name)
    if chrom not in lengths:
        return None

    arm = "p" if match.group("arm") == "p" else "q"
    panel_lookup = {(chrom, arm): (idx + 1, label) for idx, (chrom, arm, label) in enumerate(PANEL_ORDER)}
    panel = panel_lookup.get((chrom, arm))
    if panel is None:
        return None

    if row["region_start"] == "." or row["region_end"] == ".":
        return None

    seq_start = int(match.group("seq_start"))
    seq_end = int(match.group("seq_end"))
    region_start = int(row["region_start"])
    region_end = int(row["region_end"])
    if region_end <= region_start:
        return None

    full_start = seq_start + region_start
    full_end = seq_start + region_end
    query_length = lengths[chrom][1]

    if arm == "p":
        plot_start = max(0, min(ZOOM_BP, full_start))
        plot_end = max(0, min(ZOOM_BP, full_end))
    else:
        window_start = query_length - ZOOM_BP
        plot_start = max(0, min(ZOOM_BP, full_start - window_start))
        plot_end = max(0, min(ZOOM_BP, full_end - window_start))

    if plot_end <= plot_start:
        return None

    panel_order, panel_label = panel
    return {
        "panel_order": panel_order,
        "panel_id": f"{chrom}_{arm}",
        "panel_label": panel_label,
        "query_name": query_name,
        "query_chrom": chrom,
        "arm": arm,
        "query_length": query_length,
        "zoom_bp": ZOOM_BP,
        "phr_source": str(source_path),
        "phr_seq": row["seq"],
        "flank_start": seq_start,
        "flank_end": seq_end + 1,
        "region_start": region_start,
        "region_end": region_end,
        "query_full_start": full_start,
        "query_full_end": full_end,
        "plot_start": plot_start,
        "plot_end": plot_end,
        "chrs_involved": row.get("chrs_involved", "NA"),
        "arms_involved": row.get("arms_involved", "NA"),
    }


def read_phr_table(path: Path, lengths: dict[str, tuple[str, int]]) -> dict[tuple[str, str], list[dict[str, object]]]:
    out: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            projected = project_phr_interval(row, lengths, path)
            if projected is None:
                continue
            out[(str(projected["query_chrom"]), str(projected["arm"]))].append(projected)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--phr-only",
        action="store_true",
        help="Refresh only projected WashU PHR intervals; keep committed window snapshots untouched.",
    )
    args = parser.parse_args()

    lengths = read_fai(QUERY_FAI)
    phr_by_chrom_arm = read_phr_table(PHR_TABLE, lengths)
    panel_by_chrom_arm = {(chrom, arm): (idx + 1, label) for idx, (chrom, arm, label) in enumerate(PANEL_ORDER)}
    groups: dict[tuple[str, int, int], dict[str, dict[str, str]]] = defaultdict(dict)

    if not args.phr_only:
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
                    "target_haplotype": target_haplotype(inter["other_seq"]),
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
    by_panel_target: dict[tuple[str, str, str], int] = defaultdict(int)
    for row in segments:
        by_panel[str(row["panel_id"])].append(row)
        by_panel_target[
            (str(row["panel_id"]), str(row["target_chrom"]), str(row["target_haplotype"]))
        ] += int(row["window_overlap_bp"])

    for panel_order, (chrom, arm, label) in enumerate(PANEL_ORDER, start=1):
        panel_id = f"{chrom}_{arm}"
        rows = by_panel.get(panel_id, [])
        top_targets = [
            f"{target} {hap}:{bp}"
            for (pid, target, hap), bp in sorted(
                by_panel_target.items(),
                key=lambda kv: (kv[0][0] != panel_id, -kv[1], kv[0][1], kv[0][2]),
            )
            if pid == panel_id
        ]
        phr_rows = phr_by_chrom_arm.get((chrom, arm), [])
        first_phr = phr_rows[0] if phr_rows else {}
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
                "phr_interval_count": len(phr_rows),
                "phr_source": first_phr.get("phr_source", "NA"),
                "phr_seq": first_phr.get("phr_seq", "NA"),
                "phr_region_start": first_phr.get("region_start", "NA"),
                "phr_region_end": first_phr.get("region_end", "NA"),
                "phr_full_start": first_phr.get("query_full_start", "NA"),
                "phr_full_end": first_phr.get("query_full_end", "NA"),
            }
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    phr_fields = [
        "panel_order",
        "panel_id",
        "panel_label",
        "query_name",
        "query_chrom",
        "arm",
        "query_length",
        "zoom_bp",
        "phr_source",
        "phr_seq",
        "flank_start",
        "flank_end",
        "region_start",
        "region_end",
        "query_full_start",
        "query_full_end",
        "plot_start",
        "plot_end",
        "chrs_involved",
        "arms_involved",
    ]
    phr_rows_flat = [
        row
        for rows in phr_by_chrom_arm.values()
        for row in rows
    ]
    with (OUT_DIR / "zoom_phr_intervals.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=phr_fields, lineterminator="\n")
        writer.writeheader()
        for row in sorted(phr_rows_flat, key=lambda r: (int(r["panel_order"]), int(r["plot_start"]), int(r["plot_end"]))):
            writer.writerow(row)

    if args.phr_only:
        return

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
        "target_haplotype",
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
        "phr_interval_count",
        "phr_source",
        "phr_seq",
        "phr_region_start",
        "phr_region_end",
        "phr_full_start",
        "phr_full_end",
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
