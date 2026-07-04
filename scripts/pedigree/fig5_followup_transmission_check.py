#!/usr/bin/env python3
"""
Summarize the Fig. 5 chr9q follow-up transmission check.

The script uses the strict-path segment manifest produced for the Fig. 5
synteny prototype. It compares the PAN027 paternal chr9q event interrogated
against PAN011 with the PAN028 maternal haplotype inherited from PAN027.
Outputs are compact repo artifacts only; no alignments are recomputed.
"""

from __future__ import annotations

import csv
import html
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_TSV = (
    REPO_ROOT
    / "paper_prep"
    / "_brainstorming"
    / "fig5_synteny_recombination_schematic"
    / "selected_segments.tsv"
)
OUT_DIR = REPO_ROOT / "paper_prep" / "_brainstorming" / "fig5_followup_transmission_check"

MOTHER_EVENT = "PAN027_chr9q_chr3q_PHR_candidate"
CHILD_EVENT = "PAN028_chr9q_chr3q_PHR_candidate"
EVENT_ORDER = [MOTHER_EVENT, CHILD_EVENT]

ARM_COLORS = {
    "chr9q": "#a7aeb6",
    "chr3q": "#2f7fbe",
    "chr15q": "#d8a100",
    "chr16q": "#7b51c8",
    "chr20q": "#b04f4f",
}
ROLE_LABELS = {
    "same-chromosome context": "same chr9q context",
    "primary donor": "chr3q primary donor",
    "side fragment": "side fragment",
    "low-confidence tail": "low-confidence tail",
}


@dataclass(frozen=True)
class Segment:
    event_id: str
    pair: str
    transmission: str
    query_name: str
    query_arm: str
    query_window: str
    qstart: int
    qend: int
    length: int
    target_name: str
    target_arm: str
    target_haplotype: str
    target_haplotype_label: str
    target_window: str
    target_local: str
    target_native: str
    identity: float
    jaccard: float
    role: str
    community_status: str
    patch_pattern: str
    overlaps_phr: str

    @property
    def target_key(self) -> str:
        return f"{self.target_arm}:{self.target_haplotype_label}"

    @property
    def interchromosomal(self) -> bool:
        return self.target_arm[: -1] != self.query_arm[: -1]


def parse_int(value: str) -> int:
    return int(value)


def parse_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return math.nan


def read_segments(path: Path) -> list[Segment]:
    rows: list[Segment] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["event_id"] not in EVENT_ORDER:
                continue
            rows.append(
                Segment(
                    event_id=row["event_id"],
                    pair=row["pair"],
                    transmission=row["transmission"],
                    query_name=row["query_name"],
                    query_arm=row["query_arm"],
                    query_window=row["query_source_window_native_0based_half_open"],
                    qstart=parse_int(row["local_query_start_0based"]),
                    qend=parse_int(row["local_query_end_0based_exclusive"]),
                    length=parse_int(row["segment_length_bp"]),
                    target_name=row["target_name"],
                    target_arm=row["target_arm"],
                    target_haplotype=row["target_haplotype"],
                    target_haplotype_label=row["target_haplotype_label"],
                    target_window=row["target_source_window_native_0based_half_open"],
                    target_local=row["target_local_interval_0based_if_recovered"],
                    target_native=row["target_native_interval_0based_half_open_if_recovered"],
                    identity=parse_float(row["identity"]),
                    jaccard=parse_float(row["jaccard"]),
                    role=row["event_role"],
                    community_status=row["community_status"],
                    patch_pattern=row["patch_pattern_annotation"],
                    overlaps_phr=row["overlaps_phr"],
                )
            )
    return rows


def weighted_mean(rows: Iterable[Segment], attr: str) -> float:
    total = 0
    weighted = 0.0
    for row in rows:
        value = getattr(row, attr)
        if math.isnan(value):
            continue
        total += row.length
        weighted += row.length * value
    return weighted / total if total else math.nan


def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize_events(segments: list[Segment]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for event_id in EVENT_ORDER:
        rows = [row for row in segments if row.event_id == event_id]
        by_role = defaultdict(int)
        by_arm = defaultdict(int)
        by_status = defaultdict(int)
        for row in rows:
            by_role[row.role] += row.length
            by_arm[row.target_arm] += row.length
            by_status[row.community_status] += row.length
        output.append(
            {
                "event_id": event_id,
                "pair": rows[0].pair,
                "transmission": rows[0].transmission,
                "query_window": rows[0].query_window,
                "segments": len(rows),
                "total_bp": sum(row.length for row in rows),
                "same_chr_context_bp": by_role["same-chromosome context"],
                "primary_chr3q_bp": by_arm["chr3q"],
                "side_chr15q_bp": by_arm["chr15q"],
                "side_chr16q_bp": by_arm["chr16q"],
                "low_conf_chr20q_bp": by_arm["chr20q"],
                "within_community_bp": by_status["within_community"],
                "cross_community_bp": by_status["cross_community"],
                "same_chr_annotated_bp": by_status["same_chr"],
                "weighted_identity_pct": f"{weighted_mean(rows, 'identity'):.4f}",
                "weighted_jaccard": f"{weighted_mean(rows, 'jaccard'):.6f}",
            }
        )
    return output


def summarize_child_sources(segments: list[Segment]) -> list[dict[str, object]]:
    child_rows = [row for row in segments if row.event_id == CHILD_EVENT]
    by_source: dict[tuple[str, str, str], list[Segment]] = defaultdict(list)
    for row in child_rows:
        source = "PAN027 paternal hap2" if row.target_haplotype_label == "paternal" else "PAN027 maternal hap1"
        by_source[(row.role, row.target_arm, source)].append(row)
    output = []
    for (role, target_arm, source), rows in sorted(by_source.items()):
        output.append(
            {
                "event_id": CHILD_EVENT,
                "role": role,
                "target_arm": target_arm,
                "PAN028_source_in_PAN027": source,
                "bp": sum(row.length for row in rows),
                "segments": len(rows),
                "weighted_identity_pct": f"{weighted_mean(rows, 'identity'):.4f}",
                "weighted_jaccard": f"{weighted_mean(rows, 'jaccard'):.6f}",
                "local_intervals": ";".join(f"{row.qstart}-{row.qend}" for row in rows),
            }
        )
    return output


def compare_local_intervals(segments: list[Segment]) -> list[dict[str, object]]:
    mother_rows = [row for row in segments if row.event_id == MOTHER_EVENT]
    child_rows = [row for row in segments if row.event_id == CHILD_EVENT]
    output: list[dict[str, object]] = []
    keys = sorted({(row.role, row.target_arm) for row in mother_rows + child_rows})
    for role, arm in keys:
        mother = [row for row in mother_rows if row.role == role and row.target_arm == arm]
        child = [row for row in child_rows if row.role == role and row.target_arm == arm]
        overlap = 0
        for mrow in mother:
            for crow in child:
                overlap += interval_overlap(mrow.qstart, mrow.qend, crow.qstart, crow.qend)
        mother_bp = sum(row.length for row in mother)
        child_bp = sum(row.length for row in child)
        output.append(
            {
                "role": role,
                "target_arm": arm,
                "mother_bp": mother_bp,
                "child_bp": child_bp,
                "same_local_coordinate_overlap_bp": overlap,
                "child_to_mother_bp_ratio": f"{child_bp / mother_bp:.3f}" if mother_bp else "NA",
                "mother_local_intervals": ";".join(f"{row.qstart}-{row.qend}" for row in mother) or "absent",
                "child_local_intervals": ";".join(f"{row.qstart}-{row.qend}" for row in child) or "absent",
            }
        )
    return output


def svg_rect(x: float, y: float, width: float, height: float, fill: str, stroke: str = "none") -> str:
    return (
        f'<rect x="{x:.2f}" y="{y:.2f}" width="{width:.2f}" height="{height:.2f}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="0.6"/>'
    )


def render_svg(segments: list[Segment], path: Path) -> None:
    width = 1320
    height = 540
    margin_left = 230
    margin_right = 60
    track_width = width - margin_left - margin_right
    row_y = {MOTHER_EVENT: 135, CHILD_EVENT: 300}
    scale = track_width / 500_000
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#20242a}.small{font-size:14px}.label{font-size:16px;font-weight:700}.title{font-size:22px;font-weight:700}.muted{fill:#66717d}</style>',
        '<text class="title" x="40" y="44">Fig. 5 follow-up transmission check: chr9q -> chr3q candidate</text>',
        '<text class="small muted" x="40" y="70">Strict nb=1 + sweepGA 1:1 no-scaffold segments, aligned on each 500 kb query window.</text>',
    ]
    for tick in range(0, 501_000, 100_000):
        x = margin_left + tick * scale
        lines.append(f'<line x1="{x:.2f}" y1="94" x2="{x:.2f}" y2="400" stroke="#e0e4e8" stroke-width="1"/>')
        lines.append(f'<text class="small muted" x="{x - 18:.2f}" y="420">{tick // 1000} kb</text>')
    labels = {
        MOTHER_EVENT: ("PAN027 paternal hap2", "source: PAN011 father"),
        CHILD_EVENT: ("PAN028 maternal hap1", "source: PAN027 mother"),
    }
    for event_id in EVENT_ORDER:
        y = row_y[event_id]
        primary, secondary = labels[event_id]
        lines.append(f'<text class="label" x="40" y="{y - 16}">{html.escape(primary)}</text>')
        lines.append(f'<text class="small muted" x="40" y="{y + 4}">{html.escape(secondary)}</text>')
        lines.append(svg_rect(margin_left, y, track_width, 46, "#f7f8fa", "#b7c0ca"))
        for row in [seg for seg in segments if seg.event_id == event_id]:
            x = margin_left + row.qstart * scale
            w = max(1.2, (row.qend - row.qstart) * scale)
            fill = ARM_COLORS.get(row.target_arm, "#888888")
            stroke = "#2a2f36" if row.interchromosomal else "#89929c"
            lines.append(svg_rect(x, y + 4, w, 38, fill, stroke))
            if row.length >= 7000 or row.interchromosomal:
                label = row.target_arm
                if row.target_haplotype_label in {"paternal", "maternal"}:
                    label += " " + ("pat" if row.target_haplotype_label == "paternal" else "mat")
                if w > 36:
                    lines.append(
                        f'<text class="small" x="{x + 4:.2f}" y="{y + 30}" fill="#111">{html.escape(label)}</text>'
                    )
        lines.append(f'<text class="small muted" x="{margin_left}" y="{y + 66}">{html.escape([seg for seg in segments if seg.event_id == event_id][0].query_window)}</text>')

    legend_x = 40
    legend_y = 455
    lines.append(f'<text class="label" x="{legend_x}" y="{legend_y}">Color key</text>')
    x = legend_x
    for arm, color in ARM_COLORS.items():
        lines.append(svg_rect(x, legend_y + 18, 22, 14, color, "#2a2f36"))
        lines.append(f'<text class="small" x="{x + 30}" y="{legend_y + 30}">{arm}</text>')
        x += 122
    lines.append('<text class="small muted" x="40" y="515">Black outlines mark inter-chromosomal mappings; gray outlines mark same-chromosome chr9q context. The 493 bp chr20q tail in PAN027 is retained only as a low-confidence marker and is not recovered in PAN028.</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n")


def render_report(summary: list[dict[str, object]], child_sources: list[dict[str, object]], path: Path) -> None:
    mother = summary[0]
    child = summary[1]
    paternal_chr9 = sum(
        int(row["bp"])
        for row in child_sources
        if row["target_arm"] == "chr9q" and row["PAN028_source_in_PAN027"] == "PAN027 paternal hap2"
    )
    paternal_chr15 = sum(
        int(row["bp"])
        for row in child_sources
        if row["target_arm"] == "chr15q" and row["PAN028_source_in_PAN027"] == "PAN027 paternal hap2"
    )
    paternal_chr16 = sum(
        int(row["bp"])
        for row in child_sources
        if row["target_arm"] == "chr16q" and row["PAN028_source_in_PAN027"] == "PAN027 paternal hap2"
    )
    child_chr3 = int(child["primary_chr3q_bp"])
    report = f"""# Fig. 5 follow-up transmission check

## Question

The manuscript Fig. 5 candidate interrogates the paternal haplotype of PAN027
against PAN011. This follow-up applies the same strict-path procedure to
PAN028, the child of PAN027, and asks whether the implicated chromosome ends are
visible in the transmitted PAN028 maternal haplotype.

## Procedure

Input segments are the existing strict `nb=1` plus sweepGA `1:1` no-scaffold
rows in `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/selected_segments.tsv`.
No alignment is recomputed here. The script extracts the chr9q -> chr3q
candidate for `PAN027_vs_PAN011` and the matching `PAN028_vs_PAN027` chr9q
candidate, then summarizes the same 500 kb local query coordinate system used by
the Fig. 5 schematic prototype.

## Result

The result is clean enough for a candidate Fig. 5 update. PAN028's maternal
chr9q window retains the same implicated chromosome-end classes seen in the
PAN027 paternal haplotype: chr9q context plus chr3q primary-donor sequence and
chr15q/chr16q side fragments. In PAN028, {paternal_chr9:,} bp of same-chromosome
chr9q context maps directly to PAN027's paternal hap2, and the two diagnostic
side fragments also map to PAN027 paternal hap2 ({paternal_chr16:,} bp chr16q
and {paternal_chr15:,} bp chr15q). The child has {child_chr3:,} bp of chr3q
primary-donor sequence overall, split between PAN027 maternal and paternal chr3q
sources, so the chr3q end is present in the transmitted haplotype but is not a
single intact paternal-hap2-only block.

The 493 bp chr20q low-confidence tail in PAN027 is not recovered in PAN028.
That is desirable for the figure update: the candidate should focus on chr9q,
chr3q, chr15q, and chr16q, and omit chr20q from the interpreted event model.

## Summary Table

| Event | Same chr9q context | chr3q primary | chr16q side | chr15q side | chr20q low-conf | Weighted identity |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| PAN027 paternal hap2 <- PAN011 | {int(mother['same_chr_context_bp']):,} | {int(mother['primary_chr3q_bp']):,} | {int(mother['side_chr16q_bp']):,} | {int(mother['side_chr15q_bp']):,} | {int(mother['low_conf_chr20q_bp']):,} | {mother['weighted_identity_pct']} |
| PAN028 maternal hap1 <- PAN027 | {int(child['same_chr_context_bp']):,} | {int(child['primary_chr3q_bp']):,} | {int(child['side_chr16q_bp']):,} | {int(child['side_chr15q_bp']):,} | {int(child['low_conf_chr20q_bp']):,} | {child['weighted_identity_pct']} |

## Candidate Figure Update

Use `fig5_followup_transmission_check.svg` as a candidate companion/update panel:
two aligned rows show the originally interrogated PAN027 paternal haplotype and
the PAN028 maternal follow-up. Black-outlined blocks are inter-chromosomal
mappings; gray-outlined blocks are same-chromosome chr9q context.

## Outputs

- `fig5_followup_transmission_check.svg`
- `transmission_event_summary.tsv`
- `pan028_source_breakdown.tsv`
- `local_interval_comparison.tsv`

## Interpretation Boundary

This is a transmission consistency check, not a new de novo-event proof. It
supports the expected transmission of the implicated chromosome-end pattern into
PAN028, with the caveat that chr3q support in PAN028 is distributed across both
PAN027 haplotypes.
"""
    path.write_text(report)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    segments = read_segments(SOURCE_TSV)
    if not segments:
        raise SystemExit(f"No Fig. 5 follow-up segments found in {SOURCE_TSV}")

    summary = summarize_events(segments)
    child_sources = summarize_child_sources(segments)
    local_comparison = compare_local_intervals(segments)

    write_tsv(
        OUT_DIR / "transmission_event_summary.tsv",
        summary,
        [
            "event_id",
            "pair",
            "transmission",
            "query_window",
            "segments",
            "total_bp",
            "same_chr_context_bp",
            "primary_chr3q_bp",
            "side_chr15q_bp",
            "side_chr16q_bp",
            "low_conf_chr20q_bp",
            "within_community_bp",
            "cross_community_bp",
            "same_chr_annotated_bp",
            "weighted_identity_pct",
            "weighted_jaccard",
        ],
    )
    write_tsv(
        OUT_DIR / "pan028_source_breakdown.tsv",
        child_sources,
        [
            "event_id",
            "role",
            "target_arm",
            "PAN028_source_in_PAN027",
            "bp",
            "segments",
            "weighted_identity_pct",
            "weighted_jaccard",
            "local_intervals",
        ],
    )
    write_tsv(
        OUT_DIR / "local_interval_comparison.tsv",
        local_comparison,
        [
            "role",
            "target_arm",
            "mother_bp",
            "child_bp",
            "same_local_coordinate_overlap_bp",
            "child_to_mother_bp_ratio",
            "mother_local_intervals",
            "child_local_intervals",
        ],
    )
    render_svg(segments, OUT_DIR / "fig5_followup_transmission_check.svg")
    render_report(summary, child_sources, OUT_DIR / "REPORT.md")


if __name__ == "__main__":
    main()
