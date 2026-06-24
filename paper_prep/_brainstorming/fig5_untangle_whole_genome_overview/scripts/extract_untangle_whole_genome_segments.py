#!/usr/bin/env python3
"""Build whole-genome strict untangle overview TSVs for Fig5."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


STRICT_SEGMENTS = Path("paper_prep/_brainstorming/fig5_sweepga_1to1_redraw/conservative_segments.tsv")
EVENT_DIR = Path("paper_prep/_brainstorming/fig5_synteny_recombination_schematic")
SELECTED_SEGMENTS = EVENT_DIR / "selected_segments.tsv"
EVENT_MANIFEST = EVENT_DIR / "event_manifest.tsv"
PATCHES = Path("/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv")

EVENT_LABELS = {
    "PAR1_XY_positive_control": "PAR1",
    "PAN027_chr9q_chr3q_PHR_candidate": "PAN027 chr9q->chr3q",
    "PAN028_chr9q_chr3q_PHR_candidate": "PAN028 chr9q->chr3q",
}

ARM_SUFFIX = {"parm": "p", "qarm": "q"}


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


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    if chrom.startswith("chr"):
        rest = chrom[3:]
    else:
        rest = chrom
    if rest.isdigit():
        return (int(rest), "")
    if rest == "X":
        return (23, "")
    if rest == "Y":
        return (24, "")
    return (99, rest)


def arm_sort_key(arm: str) -> tuple[int, int, str]:
    chrom = arm[:-1]
    side = arm[-1]
    side_order = 0 if side == "p" else 1
    ckey = chrom_sort_key(chrom)
    return (ckey[0], side_order, ckey[1])


def parse_window_name(name: str) -> dict[str, object]:
    pattern = re.compile(
        r"^(?P<sample>[^#]+)#(?P<hapnum>[^#]+)#(?P<chrom>chr[^.:]+)\."
        r"(?P<label>[^:]+):(?P<start>\d+)-(?P<end>\d+)_(?P<armchrom>chr[^_]+)_(?P<armside>parm|qarm)$"
    )
    match = pattern.match(name)
    if not match:
        raise ValueError(f"could not parse untangle window name: {name}")
    values = match.groupdict()
    source_start = int(values["start"])
    source_end_inclusive = int(values["end"])
    source_end = source_end_inclusive + 1
    arm = f"{values['armchrom']}{ARM_SUFFIX[values['armside']]}"
    return {
        "sample": values["sample"],
        "haplotype": f"{values['sample']}#{values['hapnum']}",
        "haplotype_number": values["hapnum"],
        "haplotype_label": values["label"],
        "chrom": values["chrom"],
        "arm": arm,
        "arm_side": ARM_SUFFIX[values["armside"]],
        "source_start": source_start,
        "source_end": source_end,
        "source_interval": f"{values['chrom']}:{source_start}-{source_end}",
    }


def parse_chr_from_arm(arm: str) -> str:
    return arm[:-1]


def fmt_bp(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}Mb"
    if value >= 1_000:
        return f"{value / 1_000:.1f}kb"
    return f"{value}bp"


def make_event_join(selected_rows: list[dict[str, str]], event_rows: list[dict[str, str]]) -> dict[tuple[str, ...], dict[str, str]]:
    events = {row["event_id"]: row for row in event_rows}
    joined: dict[tuple[str, ...], dict[str, str]] = {}
    for row in selected_rows:
        event = events[row["event_id"]]
        key = (
            row["pair"],
            row["query_name"],
            row["local_query_start_0based"],
            row["local_query_end_0based_exclusive"],
            row["target_name"],
            row["target_arm"],
        )
        joined[key] = {
            "event_id": row["event_id"],
            "event_order": event["event_order"],
            "event_label": EVENT_LABELS.get(row["event_id"], row["event_id"]),
            "event_class": event["event_class"],
            "event_role": row["event_role"],
            "event_interpretation": event["interpretation_boundary"],
            "candidate_window_native_0based_half_open": event["query_native_window_0based_half_open"],
            "side_fragment_caveat": "yes"
            if row["event_role"] in {"side fragment", "low-confidence tail"}
            else "no",
        }
    return joined


def load_patch_metadata(path: Path) -> dict[tuple[str, str, str, str], str]:
    if not path.exists():
        return {}
    out: dict[tuple[str, str, str, str], str] = {}
    for row in read_tsv(path):
        key = (row["label"], row["query"], row["patch_start"], row["patch_end"])
        out[key] = ";".join(
            [
                f"pattern={row.get('pattern', '')}",
                f"query_community={row.get('query_community', '')}",
                f"ref_community={row.get('ref_community', '')}",
                f"community_status={row.get('community_status', '')}",
                f"overlaps_phr={row.get('overlaps_phr', '')}",
            ]
        )
    return out


def build_segments(strict_rows: list[dict[str, str]], event_join: dict[tuple[str, ...], dict[str, str]]) -> list[dict[str, object]]:
    arms = sorted({row["query_arm"] for row in strict_rows}, key=arm_sort_key)
    arm_index = {arm: idx + 1 for idx, arm in enumerate(arms)}
    out: list[dict[str, object]] = []

    for row_id, row in enumerate(strict_rows, start=1):
        query = parse_window_name(row["query_name"])
        target = parse_window_name(row["target_name"])
        local_start = int(row["query_start"])
        local_end = int(row["query_end"])
        native_start = int(query["source_start"]) + local_start
        native_end = int(query["source_start"]) + local_end
        segment_length = max(0, local_end - local_start)
        target_chrom = parse_chr_from_arm(row["target_arm"])
        query_chrom = parse_chr_from_arm(row["query_arm"])
        key = (
            row["pair"],
            row["query_name"],
            row["query_start"],
            row["query_end"],
            row["target_name"],
            row["target_arm"],
        )
        event = event_join.get(key, {})
        target_hap_short = target["haplotype_number"]
        color_group = row["target_arm"]
        if row["interchromosomal"] == "0":
            color_group = "same_chromosome_primary_path"
        elif event.get("side_fragment_caveat") == "yes":
            color_group = f"side_caveat_{row['target_arm']}"

        out.append(
            {
                "row_id": row_id,
                "pair": row["pair"],
                "transmission": row["transmission"],
                "query_name": row["query_name"],
                "query_sample": query["sample"],
                "query_haplotype": query["haplotype"],
                "query_haplotype_label": query["haplotype_label"],
                "query_chrom": query_chrom,
                "query_arm": row["query_arm"],
                "query_arm_order": arm_index[row["query_arm"]],
                "query_source_window_native_0based_half_open": query["source_interval"],
                "query_source_start_0based": query["source_start"],
                "query_source_end_0based_exclusive": query["source_end"],
                "local_query_start_0based": local_start,
                "local_query_end_0based_exclusive": local_end,
                "native_query_start_0based": native_start,
                "native_query_end_0based_exclusive": native_end,
                "native_query_interval_0based_half_open": f"{query_chrom}:{native_start}-{native_end}",
                "query_length": row["query_length"],
                "segment_length_bp": segment_length,
                "target_name": row["target_name"],
                "target_sample": target["sample"],
                "target_haplotype": row["target_hap"],
                "target_haplotype_number": target_hap_short,
                "target_haplotype_label": target["haplotype_label"],
                "target_chrom": target_chrom,
                "target_arm": row["target_arm"],
                "target_source_window_native_0based_half_open": target["source_interval"],
                "target_color_group": color_group,
                "strand": row["strand"],
                "identity": row["identity"],
                "jaccard": row["jaccard"],
                "nb": row["nb"],
                "interchromosomal": row["interchromosomal"],
                "strict_primary_path_only": "yes",
                "event_id": event.get("event_id", ""),
                "event_order": event.get("event_order", ""),
                "event_label": event.get("event_label", ""),
                "event_class": event.get("event_class", ""),
                "event_role": event.get("event_role", ""),
                "side_fragment_caveat": event.get("side_fragment_caveat", "no"),
                "candidate_window_native_0based_half_open": event.get(
                    "candidate_window_native_0based_half_open", ""
                ),
                "event_interpretation": event.get("event_interpretation", ""),
                "coordinate_system": "native sample assembly window coordinates; not CHM13-projected",
                "geometry_source": str(STRICT_SEGMENTS),
                "annotation_source": str(EVENT_MANIFEST) if event else "",
                "do_not_use_for_geometry": "permissive multimap/nth-best rows; patches.tsv patch extents",
            }
        )

    out.sort(
        key=lambda row: (
            str(row["pair"]),
            int(row["query_arm_order"]),
            int(row["local_query_start_0based"]),
            int(row["local_query_end_0based_exclusive"]),
            str(row["target_arm"]),
            str(row["target_haplotype"]),
        )
    )
    return out


def build_summary(segment_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in segment_rows:
        grouped[
            (
                str(row["pair"]),
                str(row["transmission"]),
                str(row["query_arm"]),
                str(row["target_arm"]),
            )
        ].append(row)

    out: list[dict[str, object]] = []
    for (pair, transmission, query_arm, target_arm), rows in grouped.items():
        total_bp = sum(int(row["segment_length_bp"]) for row in rows)
        inter_bp = sum(int(row["segment_length_bp"]) for row in rows if row["interchromosomal"] == "1")
        candidate_bp = sum(int(row["segment_length_bp"]) for row in rows if row["event_id"])
        side_bp = sum(int(row["segment_length_bp"]) for row in rows if row["side_fragment_caveat"] == "yes")
        event_ids = sorted({str(row["event_id"]) for row in rows if row["event_id"]})
        roles = sorted({str(row["event_role"]) for row in rows if row["event_role"]})
        target_haps = sorted({str(row["target_haplotype"]) for row in rows})
        query_windows = sorted({str(row["query_source_window_native_0based_half_open"]) for row in rows})
        out.append(
            {
                "pair": pair,
                "transmission": transmission,
                "query_arm": query_arm,
                "query_arm_order": rows[0]["query_arm_order"],
                "target_arm": target_arm,
                "target_chrom": rows[0]["target_chrom"],
                "strict_segment_count": len(rows),
                "strict_total_bp": total_bp,
                "interchromosomal_bp": inter_bp,
                "same_chromosome_bp": total_bp - inter_bp,
                "candidate_callout_bp": candidate_bp,
                "side_fragment_caveat_bp": side_bp,
                "event_ids": ";".join(event_ids),
                "event_roles": ";".join(roles),
                "target_haplotypes": ";".join(target_haps),
                "query_windows_native": ";".join(query_windows),
                "coordinate_system": "native sample assembly window coordinates; not CHM13-projected",
                "geometry_source": str(STRICT_SEGMENTS),
            }
        )

    out.sort(
        key=lambda row: (
            str(row["pair"]),
            int(row["query_arm_order"]),
            0 if str(row["target_arm"]) == str(row["query_arm"]) else 1,
            -int(row["strict_total_bp"]),
            str(row["target_arm"]),
        )
    )
    return out


def write_manifest(panel_dir: Path, segments: list[dict[str, object]], summary: list[dict[str, object]]) -> None:
    event_rows = [row for row in segments if row["event_id"]]
    caveat_bp = sum(int(row["segment_length_bp"]) for row in segments if row["side_fragment_caveat"] == "yes")
    lines = [
        "# Fig5 untangle whole-genome overview manifest",
        "",
        f"- strict primary-path rows: {len(segments)}",
        f"- summary query-arm/target-arm rows: {len(summary)}",
        f"- candidate/callout strict rows: {len(event_rows)}",
        f"- side-fragment caveat bp retained: {fmt_bp(caveat_bp)}",
        f"- strict geometry source: `{STRICT_SEGMENTS}`",
        f"- corrected event labels: `{SELECTED_SEGMENTS}` and `{EVENT_MANIFEST}`",
        f"- optional metadata labels only: `{PATCHES}`",
        "",
        "Coordinates are native sample assembly query windows parsed from untangle",
        "path names. They are not CHM13-projected coordinates.",
    ]
    (panel_dir / "MANIFEST.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-dir", type=Path, default=Path("."))
    args = parser.parse_args()

    strict_rows = read_tsv(STRICT_SEGMENTS)
    selected_rows = read_tsv(SELECTED_SEGMENTS)
    event_rows = read_tsv(EVENT_MANIFEST)
    event_join = make_event_join(selected_rows, event_rows)
    load_patch_metadata(PATCHES)

    segments = build_segments(strict_rows, event_join)
    summary = build_summary(segments)

    segment_fields = [
        "row_id",
        "pair",
        "transmission",
        "query_name",
        "query_sample",
        "query_haplotype",
        "query_haplotype_label",
        "query_chrom",
        "query_arm",
        "query_arm_order",
        "query_source_window_native_0based_half_open",
        "query_source_start_0based",
        "query_source_end_0based_exclusive",
        "local_query_start_0based",
        "local_query_end_0based_exclusive",
        "native_query_start_0based",
        "native_query_end_0based_exclusive",
        "native_query_interval_0based_half_open",
        "query_length",
        "segment_length_bp",
        "target_name",
        "target_sample",
        "target_haplotype",
        "target_haplotype_number",
        "target_haplotype_label",
        "target_chrom",
        "target_arm",
        "target_source_window_native_0based_half_open",
        "target_color_group",
        "strand",
        "identity",
        "jaccard",
        "nb",
        "interchromosomal",
        "strict_primary_path_only",
        "event_id",
        "event_order",
        "event_label",
        "event_class",
        "event_role",
        "side_fragment_caveat",
        "candidate_window_native_0based_half_open",
        "event_interpretation",
        "coordinate_system",
        "geometry_source",
        "annotation_source",
        "do_not_use_for_geometry",
    ]
    summary_fields = [
        "pair",
        "transmission",
        "query_arm",
        "query_arm_order",
        "target_arm",
        "target_chrom",
        "strict_segment_count",
        "strict_total_bp",
        "interchromosomal_bp",
        "same_chromosome_bp",
        "candidate_callout_bp",
        "side_fragment_caveat_bp",
        "event_ids",
        "event_roles",
        "target_haplotypes",
        "query_windows_native",
        "coordinate_system",
        "geometry_source",
    ]
    write_tsv(args.panel_dir / "untangle_whole_genome_segments.tsv", segments, segment_fields)
    write_tsv(args.panel_dir / "untangle_whole_genome_summary.tsv", summary, summary_fields)
    write_manifest(args.panel_dir, segments, summary)


if __name__ == "__main__":
    main()
