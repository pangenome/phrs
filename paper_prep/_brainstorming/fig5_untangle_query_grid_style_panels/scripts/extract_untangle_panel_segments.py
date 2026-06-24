#!/usr/bin/env python3
"""Extract strict untangle primary-path rows for query-grid-style Fig5 panels."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


EVENT_ORDER = [
    "PAR1_XY_positive_control",
    "PAN027_chr9q_chr3q_PHR_candidate",
    "PAN028_chr9q_chr3q_PHR_candidate",
]

PANEL_LABELS = {
    "PAR1_XY_positive_control": "PAR1 X/Y control",
    "PAN027_chr9q_chr3q_PHR_candidate": "PAN027 chr9q -> chr3q",
    "PAN028_chr9q_chr3q_PHR_candidate": "PAN028 chr9q -> chr3q",
}

ROLE_ORDER = {
    "same-chromosome context": 1,
    "PAR positive control": 2,
    "primary donor": 3,
    "side fragment": 4,
    "low-confidence tail": 5,
}

SOURCE_DIR = Path("paper_prep/_brainstorming/fig5_synteny_recombination_schematic")
SELECTED_SEGMENTS = SOURCE_DIR / "selected_segments.tsv"
EVENT_MANIFEST = SOURCE_DIR / "event_manifest.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def split_interval(value: str, default_chrom: str | None = None) -> tuple[str, int, int]:
    if ":" in value:
        chrom, rest = value.split(":", 1)
    elif default_chrom is not None:
        chrom, rest = default_chrom, value
    else:
        raise ValueError(f"missing chromosome in interval: {value}")
    start, end = rest.split("-", 1)
    return chrom, int(start), int(end)


def strip_suffix(value: str, suffix: str) -> str:
    if value.endswith(suffix):
        return value[: -len(suffix)]
    return value


def parse_bp_map(value: str) -> dict[str, int]:
    if not value or value == "none":
        return {}
    parsed: dict[str, int] = {}
    for item in value.split(";"):
        key, size = item.split(":", 1)
        parsed[key] = int(strip_suffix(size, "bp"))
    return parsed


def fmt_bp(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f} Mb"
    if value >= 1_000:
        return f"{value / 1_000:.1f} kb"
    return f"{value} bp"


def hap_code(haplotype_label: str) -> str:
    if haplotype_label in {"haplotype1", "maternal"}:
        return "h1_or_maternal"
    if haplotype_label in {"haplotype2", "paternal"}:
        return "h2_or_paternal"
    return haplotype_label


def target_color_key(target_arm: str, event_role: str, haplotype_label: str) -> str:
    if event_role in {"side fragment", "low-confidence tail"}:
        return f"side_caveat_{target_arm}"
    if target_arm in {"chr3q", "chr9q"}:
        return f"{target_arm}_{hap_code(haplotype_label)}"
    return target_arm


def build_segments(selected_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    selected = [row for row in selected_rows if row["event_id"] in EVENT_ORDER]
    order_index = {event: i + 1 for i, event in enumerate(EVENT_ORDER)}

    out: list[dict[str, object]] = []
    for row in selected:
        query_chrom, native_start, native_end = split_interval(row["native_query_interval_0based_half_open"])
        target_chrom = strip_suffix(strip_suffix(row["target_arm"], "p"), "q")
        if row["target_native_interval_0based_half_open_if_recovered"] not in ("", "not_available"):
            target_chrom, target_start, target_end = split_interval(
                row["target_native_interval_0based_half_open_if_recovered"],
                default_chrom=target_chrom,
            )
        else:
            target_start = target_end = ""

        event_role = row["event_role"]
        side_caveat = "yes" if event_role in {"side fragment", "low-confidence tail"} else "no"
        out.append(
            {
                "event_id": row["event_id"],
                "event_order": order_index[row["event_id"]],
                "panel_label": PANEL_LABELS[row["event_id"]],
                "pair": row["pair"],
                "transmission": row["transmission"],
                "query_name": row["query_name"],
                "query_arm": row["query_arm"],
                "query_sample": row["query_sample"],
                "query_haplotype": row["query_haplotype"],
                "query_haplotype_label": row["query_haplotype_label"],
                "query_source_window_native_0based_half_open": row[
                    "query_source_window_native_0based_half_open"
                ],
                "native_query_chrom": query_chrom,
                "native_query_start_0based": native_start,
                "native_query_end_0based_exclusive": native_end,
                "native_query_interval_0based_half_open": row["native_query_interval_0based_half_open"],
                "local_query_start_0based": row["local_query_start_0based"],
                "local_query_end_0based_exclusive": row["local_query_end_0based_exclusive"],
                "segment_length_bp": row["segment_length_bp"],
                "target_name": row["target_name"],
                "target_arm": row["target_arm"],
                "target_chrom": target_chrom,
                "target_sample": row["target_sample"],
                "target_haplotype": row["target_haplotype"],
                "target_haplotype_label": row["target_haplotype_label"],
                "target_color_key": target_color_key(
                    row["target_arm"], event_role, row["target_haplotype_label"]
                ),
                "target_source_window_native_0based_half_open": row[
                    "target_source_window_native_0based_half_open"
                ],
                "target_native_interval_0based_half_open_if_recovered": row[
                    "target_native_interval_0based_half_open_if_recovered"
                ],
                "target_native_start_0based_if_recovered": target_start,
                "target_native_end_0based_exclusive_if_recovered": target_end,
                "strand": row["strand"],
                "identity": row["identity"],
                "jaccard": row["jaccard"],
                "nb": row["nb"],
                "interchromosomal": row["interchromosomal"],
                "event_role": event_role,
                "event_role_order": ROLE_ORDER.get(event_role, 99),
                "side_fragment_caveat": side_caveat,
                "query_community": row["query_community"],
                "target_community": row["target_community"],
                "community_status": row["community_status"],
                "patch_pattern_annotation": row["patch_pattern_annotation"],
                "overlaps_phr": row["overlaps_phr"],
                "has_phr": row["has_phr"],
                "panel_geometry_source": "paper_prep/_brainstorming/fig5_synteny_recombination_schematic/selected_segments.tsv",
                "geometry_source": row["geometry_source"],
                "target_coordinate_source": row["target_coordinate_source"],
                "annotation_source": row["annotation_source"],
            }
        )

    out.sort(
        key=lambda row: (
            int(row["event_order"]),
            int(row["native_query_start_0based"]),
            int(row["native_query_end_0based_exclusive"]),
            int(row["event_role_order"]),
            str(row["target_arm"]),
            str(row["target_haplotype"]),
        )
    )
    return out


def build_summary(events: list[dict[str, str]], segments: list[dict[str, object]]) -> list[dict[str, object]]:
    by_event: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in segments:
        by_event[str(row["event_id"])].append(row)

    event_by_id = {row["event_id"]: row for row in events}
    summary_rows: list[dict[str, object]] = []
    for event_id in EVENT_ORDER:
        event = event_by_id[event_id]
        rows = by_event[event_id]
        window_chrom, window_start, window_end = split_interval(event["query_native_window_0based_half_open"])
        role_bp: dict[str, int] = defaultdict(int)
        target_bp: dict[str, int] = defaultdict(int)
        target_haps: set[str] = set()
        side_labels: list[str] = []
        for row in rows:
            size = int(row["segment_length_bp"])
            role_bp[str(row["event_role"])] += size
            target_bp[str(row["target_arm"])] += size
            target_haps.add(f"{row['target_arm']}:{row['target_haplotype']}")
            if row["side_fragment_caveat"] == "yes":
                side_labels.append(f"{row['target_arm']} {fmt_bp(size)}")

        summary_rows.append(
            {
                "event_id": event_id,
                "event_order": event["event_order"],
                "panel_label": PANEL_LABELS[event_id],
                "event_class": event["event_class"],
                "pair": event["pair"],
                "transmission": event["transmission"],
                "query_arm": event["query_arm"],
                "query_window_chrom": window_chrom,
                "query_window_start_0based": window_start,
                "query_window_end_0based_exclusive": window_end,
                "query_native_window_0based_half_open": event["query_native_window_0based_half_open"],
                "coordinate_system": event["coordinate_system"],
                "strict_support_row_count": event["strict_support_row_count"],
                "strict_support_total_bp": event["strict_support_total_bp"],
                "same_chromosome_context_bp_by_arm": event["same_chromosome_context_bp_by_arm"],
                "primary_donor_arms_bp": event["primary_donor_arms_bp"],
                "side_fragment_arms_bp": event["side_fragment_arms_bp"],
                "low_confidence_tail_bp": event["low_confidence_tail_bp"],
                "segment_rows": len(rows),
                "role_bp_summary": ";".join(f"{role}:{bp}" for role, bp in sorted(role_bp.items())),
                "target_arm_bp_summary": ";".join(f"{arm}:{bp}" for arm, bp in sorted(target_bp.items())),
                "target_haplotypes": ";".join(sorted(target_haps)),
                "side_fragment_caveat_labels": ";".join(side_labels) if side_labels else "none",
                "interpretation_boundary": event["interpretation_boundary"],
                "drawing_geometry_source": event["drawing_geometry_source"],
                "annotation_source": event["annotation_source"],
                "do_not_use_for_geometry": event["do_not_use_for_geometry"],
                "status": "OK",
            }
        )
    return summary_rows


def build_manifest(events: list[dict[str, str]], segment_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    event_counts: dict[str, int] = defaultdict(int)
    for row in segment_rows:
        event_counts[str(row["event_id"])] += 1

    manifest_rows: list[dict[str, object]] = []
    for event in events:
        if event["event_id"] not in EVENT_ORDER:
            continue
        manifest_rows.append(
            {
                "event_id": event["event_id"],
                "event_order": event["event_order"],
                "panel_label": PANEL_LABELS[event["event_id"]],
                "pair": event["pair"],
                "query_arm": event["query_arm"],
                "query_native_window_0based_half_open": event["query_native_window_0based_half_open"],
                "coordinate_system": event["coordinate_system"],
                "strict_support_row_count_manifest": event["strict_support_row_count"],
                "strict_support_row_count_written": event_counts[event["event_id"]],
                "geometry_source": "paper_prep/_brainstorming/fig5_synteny_recombination_schematic/selected_segments.tsv",
                "event_source": "paper_prep/_brainstorming/fig5_synteny_recombination_schematic/event_manifest.tsv",
                "optional_annotation_source": "/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv labels only when exactly joined upstream",
                "raw_fasta_query_grid_source": "not_used",
                "excluded_geometry": event["do_not_use_for_geometry"],
                "status": "OK",
            }
        )
    return manifest_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-dir", type=Path, required=True)
    parser.add_argument("--segments", type=Path, default=SELECTED_SEGMENTS)
    parser.add_argument("--events", type=Path, default=EVENT_MANIFEST)
    args = parser.parse_args()

    events = read_tsv(args.events)
    selected_rows = read_tsv(args.segments)
    segments = build_segments(selected_rows)
    summary = build_summary(events, segments)
    manifest = build_manifest(events, segments)

    segment_fields = [
        "event_id",
        "event_order",
        "panel_label",
        "pair",
        "transmission",
        "query_name",
        "query_arm",
        "query_sample",
        "query_haplotype",
        "query_haplotype_label",
        "query_source_window_native_0based_half_open",
        "native_query_chrom",
        "native_query_start_0based",
        "native_query_end_0based_exclusive",
        "native_query_interval_0based_half_open",
        "local_query_start_0based",
        "local_query_end_0based_exclusive",
        "segment_length_bp",
        "target_name",
        "target_arm",
        "target_chrom",
        "target_sample",
        "target_haplotype",
        "target_haplotype_label",
        "target_color_key",
        "target_source_window_native_0based_half_open",
        "target_native_interval_0based_half_open_if_recovered",
        "target_native_start_0based_if_recovered",
        "target_native_end_0based_exclusive_if_recovered",
        "strand",
        "identity",
        "jaccard",
        "nb",
        "interchromosomal",
        "event_role",
        "event_role_order",
        "side_fragment_caveat",
        "query_community",
        "target_community",
        "community_status",
        "patch_pattern_annotation",
        "overlaps_phr",
        "has_phr",
        "panel_geometry_source",
        "geometry_source",
        "target_coordinate_source",
        "annotation_source",
    ]
    summary_fields = [
        "event_id",
        "event_order",
        "panel_label",
        "event_class",
        "pair",
        "transmission",
        "query_arm",
        "query_window_chrom",
        "query_window_start_0based",
        "query_window_end_0based_exclusive",
        "query_native_window_0based_half_open",
        "coordinate_system",
        "strict_support_row_count",
        "strict_support_total_bp",
        "same_chromosome_context_bp_by_arm",
        "primary_donor_arms_bp",
        "side_fragment_arms_bp",
        "low_confidence_tail_bp",
        "segment_rows",
        "role_bp_summary",
        "target_arm_bp_summary",
        "target_haplotypes",
        "side_fragment_caveat_labels",
        "interpretation_boundary",
        "drawing_geometry_source",
        "annotation_source",
        "do_not_use_for_geometry",
        "status",
    ]
    manifest_fields = [
        "event_id",
        "event_order",
        "panel_label",
        "pair",
        "query_arm",
        "query_native_window_0based_half_open",
        "coordinate_system",
        "strict_support_row_count_manifest",
        "strict_support_row_count_written",
        "geometry_source",
        "event_source",
        "optional_annotation_source",
        "raw_fasta_query_grid_source",
        "excluded_geometry",
        "status",
    ]

    write_tsv(args.panel_dir / "untangle_panel_segments.tsv", segments, segment_fields)
    write_tsv(args.panel_dir / "untangle_panel_summary.tsv", summary, summary_fields)
    write_tsv(args.panel_dir / "untangle_panel_manifest.tsv", manifest, manifest_fields)


if __name__ == "__main__":
    main()
