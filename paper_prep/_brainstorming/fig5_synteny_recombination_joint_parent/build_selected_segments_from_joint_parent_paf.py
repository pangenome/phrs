#!/usr/bin/env python3
"""Build Fig5 schematic geometry from joint-parent direct sweepGA PAFs.

Regenerate from the repository root with:

    python3 paper_prep/_brainstorming/fig5_synteny_recombination_joint_parent/build_selected_segments_from_joint_parent_paf.py --filter-id one_one_noscaffold

The output intentionally matches the column schema consumed by the existing
Fig5 synteny schematic renderer. Query and target coordinates are 0-500 kb
local offsets in the telomeric source windows parsed from PAF sequence names;
native assembly coordinates are reconstructed only for labels. No whole-genome
alignment or CHM13/reference projection is performed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
SOURCE_DIR = HERE.parent / "fig5_synteny_recombination_schematic"
JOINT_DIR = HERE.parent / "pedigree_direct_sweepga_joint_parent"
EVENTS_TSV = HERE / "event_manifest.tsv"
SOURCE_SEGMENTS_TSV = SOURCE_DIR / "selected_segments.tsv"
FILTERED_PAF_DIR = JOINT_DIR / "filtered_paf"
RAW_PAF_DIR = JOINT_DIR / "raw_paf"

FIELDNAMES = [
    "event_id",
    "pair",
    "transmission",
    "query_name",
    "query_arm",
    "query_sample",
    "query_haplotype",
    "query_haplotype_label",
    "query_source_window_native_0based_half_open",
    "local_query_start_0based",
    "local_query_end_0based_exclusive",
    "native_query_start_0based",
    "native_query_end_0based_exclusive",
    "native_query_interval_0based_half_open",
    "segment_length_bp",
    "target_name",
    "target_arm",
    "target_sample",
    "target_haplotype",
    "target_haplotype_label",
    "target_source_window_native_0based_half_open",
    "target_local_interval_0based_if_recovered",
    "target_native_interval_0based_half_open_if_recovered",
    "target_interval_recovery",
    "raw_strict_paf_row_count",
    "raw_strict_paf_query_coverage_check",
    "strand",
    "identity",
    "jaccard",
    "nb",
    "interchromosomal",
    "event_role",
    "query_community",
    "target_community",
    "community_status",
    "patch_pattern_annotation",
    "overlaps_phr",
    "has_phr",
    "geometry_source",
    "target_coordinate_source",
    "annotation_source",
]


@dataclass(frozen=True)
class SeqName:
    raw: str
    sample: str
    haplotype: str
    haplotype_label: str
    chrom: str
    native_start: int
    native_end_inclusive: int
    arm: str

    @property
    def window_half_open(self) -> str:
        return f"{self.chrom}:{self.native_start}-{self.native_end_inclusive + 1}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def parse_seq_name(value: str) -> SeqName:
    match = re.match(
        r"^(?P<sample>[^#]+)#(?P<hap>[^#]+)#(?P<chrom>[^.:]+)\.(?P<label>[^:]+):"
        r"(?P<start>\d+)-(?P<end>\d+)_(?P<arm_chrom>chr[^_]+)_(?P<arm_kind>[pq]arm)$",
        value,
    )
    if not match:
        raise ValueError(f"Cannot parse PAF sequence name: {value}")
    arm_suffix = "p" if match.group("arm_kind") == "parm" else "q"
    return SeqName(
        raw=value,
        sample=match.group("sample"),
        haplotype=match.group("hap"),
        haplotype_label=match.group("label"),
        chrom=match.group("chrom"),
        native_start=int(match.group("start")),
        native_end_inclusive=int(match.group("end")),
        arm=f"{match.group('arm_chrom')}{arm_suffix}",
    )


def displayed_haplotype(seq: SeqName) -> str:
    if seq.haplotype != "joint":
        return f"{seq.sample}#{seq.haplotype}"
    label_to_hap = {
        "haplotype1": "1",
        "maternal": "1",
        "haplotype2": "2",
        "paternal": "2",
    }
    return f"{seq.sample}#{label_to_hap.get(seq.haplotype_label, seq.haplotype_label)}"


def comparison_ids(event: dict[str, str]) -> list[str]:
    if event["pair"] == "PAN027_vs_PAN011":
        return ["PAN027pat_vs_PAN011_joint"]
    if event["pair"] == "PAN027_vs_PAN010":
        return ["PAN027mat_vs_PAN010_joint"]
    if event["pair"] == "PAN028_vs_PAN027":
        return ["PAN028mat_vs_PAN027_joint"]
    raise ValueError(f"No direct sweepGA comparison mapping for {event['event_id']}")


def role_for(event: dict[str, str], target: SeqName) -> str:
    if event["event_id"] == "PAR1_XY_positive_control":
        if target.arm == event["query_arm"]:
            return "same-chromosome context"
        if target.arm == "chrYp":
            return "PAR positive control"
        return "side fragment"
    if target.arm == event["query_arm"]:
        return "same-chromosome context"
    if target.arm == "chr3q":
        return "primary donor"
    if target.arm == "chr20q":
        return "low-confidence tail"
    return "side fragment"


def annotation_lookup() -> dict[tuple[str, str, str], dict[str, str]]:
    lookup: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in read_tsv(SOURCE_SEGMENTS_TSV):
        key = (row["event_id"], row["target_arm"], row["event_role"])
        lookup.setdefault(key, row)
    return lookup


def paf_rows(path: Path, query_name: str) -> list[list[str]]:
    rows: list[list[str]] = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 12 and fields[0] == query_name:
                rows.append(fields)
    return rows


def make_row(
    event: dict[str, str],
    fields: list[str],
    paf_path: Path,
    filter_id: str,
    annotations: dict[tuple[str, str, str], dict[str, str]],
) -> dict[str, str]:
    query = parse_seq_name(fields[0])
    target = parse_seq_name(fields[5])
    q_start = int(fields[2])
    q_end = int(fields[3])
    t_start = int(fields[7])
    t_end = int(fields[8])
    matches = int(fields[9])
    aln_len = int(fields[10])
    role = role_for(event, target)
    native_q_start = query.native_start + q_start
    native_q_end = query.native_start + q_end
    native_t_start = target.native_start + t_start
    native_t_end = target.native_start + t_end
    identity = matches / aln_len if aln_len else 0.0
    query_len = max(1, q_end - q_start)
    target_len = max(1, t_end - t_start)
    jaccard = matches / max(1, query_len + target_len - matches)
    ann = annotations.get((event["event_id"], target.arm, role), {})

    return {
        "event_id": event["event_id"],
        "pair": event["pair"],
        "transmission": event["transmission"],
        "query_name": query.raw,
        "query_arm": event["query_arm"],
        "query_sample": query.sample,
        "query_haplotype": query.haplotype,
        "query_haplotype_label": query.haplotype_label,
        "query_source_window_native_0based_half_open": query.window_half_open,
        "local_query_start_0based": str(q_start),
        "local_query_end_0based_exclusive": str(q_end),
        "native_query_start_0based": str(native_q_start),
        "native_query_end_0based_exclusive": str(native_q_end),
        "native_query_interval_0based_half_open": f"{query.chrom}:{native_q_start}-{native_q_end}",
        "segment_length_bp": str(q_end - q_start),
        "target_name": target.raw,
        "target_arm": target.arm,
        "target_sample": target.sample,
        "target_haplotype": displayed_haplotype(target),
        "target_haplotype_label": target.haplotype_label,
        "target_source_window_native_0based_half_open": target.window_half_open,
        "target_local_interval_0based_if_recovered": f"{t_start}-{t_end}",
        "target_native_interval_0based_half_open_if_recovered": f"{target.chrom}:{native_t_start}-{native_t_end}",
        "target_interval_recovery": f"joint_parent_direct_sweepga_{filter_id}_paf",
        "raw_strict_paf_row_count": "1",
        "raw_strict_paf_query_coverage_check": f"{q_start}-{q_end}",
        "strand": fields[4],
        "identity": f"{identity * 100:.4f}".rstrip("0").rstrip("."),
        "jaccard": f"{jaccard:.6f}".rstrip("0").rstrip("."),
        "nb": f"joint_parent_{filter_id}",
        "interchromosomal": "0" if target.chrom == query.chrom else "1",
        "event_role": role,
        "query_community": ann.get("query_community", "not_available"),
        "target_community": ann.get("target_community", "not_available"),
        "community_status": ann.get("community_status", "not_available"),
        "patch_pattern_annotation": ann.get("patch_pattern_annotation", "not_available"),
        "overlaps_phr": ann.get("overlaps_phr", "not_joined_to_graph_patch"),
        "has_phr": ann.get("has_phr", "not_joined_to_graph_patch"),
        "geometry_source": (
            f"joint-parent direct sweepGA PAF; filter_id={filter_id}; "
            "parent hap1+hap2 targets filtered jointly with --scaffold-jump 0"
        ),
        "target_coordinate_source": rel(paf_path),
        "annotation_source": (
            f"{rel(SOURCE_SEGMENTS_TSV)} for optional community/PHR labels by "
            "event_id,target_arm,event_role only; geometry is direct PAF-native"
        ),
    }


def paf_path_for(comparison_id: str, filter_id: str) -> Path:
    if filter_id == "many_many_noscaffold":
        return RAW_PAF_DIR / f"{comparison_id}.sweepga_many_many_j0.paf.gz"
    return FILTERED_PAF_DIR / f"{comparison_id}.{filter_id}.paf.gz"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filter-id",
        required=True,
        choices=[
            "one_one_noscaffold",
            "one_many_noscaffold",
            "two_many_noscaffold",
            "four_many_noscaffold",
            "many_many_noscaffold",
        ],
    )
    args = parser.parse_args()
    out_tsv = HERE / f"selected_segments.joint_parent_{args.filter_id}.tsv"
    events = read_tsv(EVENTS_TSV)
    annotations = annotation_lookup()
    output_rows: list[dict[str, str]] = []
    for event in events:
        query_name = event["child_query_source"].split(" [", 1)[0]
        for comparison_id in comparison_ids(event):
            paf_path = paf_path_for(comparison_id, args.filter_id)
            if not paf_path.exists():
                raise FileNotFoundError(paf_path)
            for fields in paf_rows(paf_path, query_name):
                output_rows.append(make_row(event, fields, paf_path, args.filter_id, annotations))

    output_rows.sort(
        key=lambda row: (
            int(next(e["event_order"] for e in events if e["event_id"] == row["event_id"])),
            int(row["local_query_start_0based"]),
            row["target_haplotype"],
            row["target_arm"],
            int(row["target_local_interval_0based_if_recovered"].split("-", 1)[0]),
        )
    )

    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Wrote {len(output_rows)} joint-parent {args.filter_id} segments to {rel(out_tsv)}")


if __name__ == "__main__":
    main()
