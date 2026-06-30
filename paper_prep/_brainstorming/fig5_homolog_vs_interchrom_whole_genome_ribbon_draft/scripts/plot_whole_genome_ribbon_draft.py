#!/usr/bin/env python3
"""Whole-genome ribbon draft for PAN027 paternal interchrom-over-homolog calls."""

from __future__ import annotations

import csv
import gzip
import html
import os
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
OUT_DIR = ROOT / "paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft"
CLASS_WINNERS_NAME = (
    "PAN027pat_vs_PAN011_joint.sweepga_f32.10to10.query_2000bp."
    "predepth_class_winners.impg_similarity.tsv.gz"
)
LOCAL_CLASS_WINNERS = (
    ROOT
    / "paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs/"
    / CLASS_WINNERS_NAME
)
MOOSEFS_CLASS_WINNERS = (
    Path("/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs")
    / CLASS_WINNERS_NAME
)
CLASS_WINNERS = LOCAL_CLASS_WINNERS if LOCAL_CLASS_WINNERS.exists() else MOOSEFS_CLASS_WINNERS
QUERY_FAI = Path(
    "/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/"
    "pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027pat_vs_PAN011_joint.query.fa.fai"
)
TARGET_FAI = Path(
    "/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/"
    "pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027pat_vs_PAN011_joint.target.fa.fai"
)

SVG_OUT = OUT_DIR / "fig5_homolog_vs_interchrom_whole_genome_ribbon_draft.svg"
HOMOLOG_SVG_OUT = OUT_DIR / "fig5_homologous_recombination_context_ribbon_draft.svg"
RUNS_OUT = OUT_DIR / "whole_genome_ribbon_runs.tsv"
SUMMARY_OUT = OUT_DIR / "whole_genome_ribbon_summary.tsv"
HOMOLOG_RUNS_OUT = OUT_DIR / "whole_genome_homologous_context_runs.tsv"
HOMOLOG_SUMMARY_OUT = OUT_DIR / "whole_genome_homologous_context_summary.tsv"
CONVERSION_STATUS = OUT_DIR / "conversion_status.txt"

CHROM_ORDER = [f"chr{i}" for i in range(1, 22 + 1)] + ["chrX", "chrY"]
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}
CHR_RE = re.compile(r"(chr(?:[0-9]+|X|Y|M))")
TARGET_RE = re.compile(r"PAN011#joint#(?P<hap>h[12])_(?P<chrom>chr(?:[0-9]+|X|Y|M))")
LOC_RE = re.compile(r"^(?P<seq>.+):(?P<start>[0-9]+)-(?P<end>[0-9]+)$")

PAGE_W = 3000
PAGE_H = 900
TRACK_X0 = 360
TRACK_W = 2360
TRACK_H = 28
Y_H1 = 235
Y_QUERY = 465
Y_H2 = 675
GRID_Y0 = 145
GRID_Y1 = 745
LEGEND_Y = 810
FOOTNOTE_Y1 = 855
FOOTNOTE_Y2 = 883
TEXT = "#202124"
MUTED = "#5f6368"
GRID = "#e8eaed"
HOMOLOG_COLOR = "#b8bdc3"
HOMOLOG_RIBBON = "#cfd3d7"
HOMOLOG_MIN_BP = 50_000
HOMOLOG_MIN_IDENTITY = 0.95

COLORS = {
    "PAR_XY": "#E7298A",
    "chr5_chr1_candidate": "#1F77B4",
    "chr9_chr3_candidate": "#D95F02",
    "acro_acro": "#6f6f6f",
    "acro_other": "#9e9e9e",
    "other_nonacro": "#2C7FB8",
}


@dataclass
class GenomeLayout:
    label: str
    offsets: dict[str, int]
    lengths: dict[str, int]
    total: int


@dataclass
class Segment:
    query_seq: str
    query_chrom: str
    query_start: int
    query_end: int
    donor_seq: str
    donor_haplotype: str
    target_chrom: str
    donor_start: int
    donor_end: int
    bp: int
    same_identity: float
    inter_identity: float


@dataclass
class Run:
    run_id: str
    query_seq: str
    query_chrom: str
    query_start: int
    query_end: int
    donor_seq: str
    donor_haplotype: str
    target_chrom: str
    donor_start: int
    donor_end: int
    bp: int
    windows: int
    weighted_same_identity: float
    weighted_inter_identity: float

    @property
    def mean_same_identity(self) -> float:
        return self.weighted_same_identity / max(1, self.bp)

    @property
    def mean_inter_identity(self) -> float:
        return self.weighted_inter_identity / max(1, self.bp)

    @property
    def category(self) -> str:
        if {self.query_chrom, self.target_chrom} == {"chrX", "chrY"}:
            return "PAR_XY"
        if self.query_chrom == "chr5" and self.target_chrom == "chr1":
            return "chr5_chr1_candidate"
        if self.query_chrom == "chr9" and self.target_chrom == "chr3":
            return "chr9_chr3_candidate"
        if self.query_chrom in ACRO and self.target_chrom in ACRO:
            return "acro_acro"
        if self.query_chrom in ACRO or self.target_chrom in ACRO:
            return "acro_other"
        return "other_nonacro"


class SVG:
    def __init__(self, width: int, height: int) -> None:
        self.width = width
        self.height = height
        self.items: list[str] = []

    def add(self, item: str) -> None:
        self.items.append(item)

    def text(
        self,
        x: float,
        y: float,
        text: str,
        size: float = 14,
        weight: str = "400",
        fill: str = TEXT,
        anchor: str = "start",
        rotate: float | None = None,
    ) -> None:
        transform = ""
        if rotate is not None:
            transform = f' transform="rotate({rotate:.1f} {x:.1f} {y:.1f})"'
        self.add(
            f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size:.1f}" '
            f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}"{transform}>'
            f"{html.escape(str(text))}</text>"
        )

    def rect(
        self,
        x: float,
        y: float,
        w: float,
        h: float,
        fill: str,
        stroke: str = "none",
        sw: float = 1.0,
        opacity: float = 1.0,
        rx: float = 0.0,
    ) -> None:
        self.add(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" '
            f'rx="{rx:.2f}" fill="{fill}" stroke="{stroke}" stroke-width="{sw:.2f}" '
            f'opacity="{opacity:.3f}"/>'
        )

    def line(self, x1: float, y1: float, x2: float, y2: float, stroke: str, sw: float = 1.0, opacity: float = 1.0) -> None:
        self.add(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="{stroke}" stroke-width="{sw:.2f}" opacity="{opacity:.3f}"/>'
        )

    def path(self, d: str, fill: str, stroke: str = "none", sw: float = 0.0, opacity: float = 1.0) -> None:
        self.add(
            f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw:.2f}" opacity="{opacity:.3f}"/>'
        )

    def write(self, path: Path) -> None:
        path.write_text(
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" height="{self.height}" '
            f'viewBox="0 0 {self.width} {self.height}">\n'
            '<rect width="100%" height="100%" fill="white"/>\n'
            '<style>text{font-family:Arial,Helvetica,sans-serif;dominant-baseline:alphabetic}</style>\n'
            + "\n".join(self.items)
            + "\n</svg>\n"
        )


def chrom_name(value: str) -> str:
    match = CHR_RE.search(value)
    return match.group(1) if match else value


def parse_loc(value: str) -> tuple[str, int, int]:
    match = LOC_RE.match(value)
    if match is None:
        raise ValueError(f"could not parse interval: {value}")
    return match.group("seq"), int(match.group("start")), int(match.group("end"))


def target_meta(seq: str) -> tuple[str, str]:
    match = TARGET_RE.match(seq)
    if match is None:
        return "NA", chrom_name(seq)
    return match.group("hap"), match.group("chrom")


def read_query_fai(path: Path) -> tuple[dict[str, str], dict[str, int]]:
    seq_by_chrom: dict[str, str] = {}
    length_by_chrom: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            seq, length, *_ = line.rstrip("\n").split("\t")
            chrom = chrom_name(seq)
            seq_by_chrom[chrom] = seq
            length_by_chrom[chrom] = int(length)
    return seq_by_chrom, length_by_chrom


def read_target_fai(path: Path) -> tuple[dict[tuple[str, str], str], dict[tuple[str, str], int], dict[str, int]]:
    seq_by_key: dict[tuple[str, str], str] = {}
    length_by_key: dict[tuple[str, str], int] = {}
    length_by_seq: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            seq, length_text, *_ = line.rstrip("\n").split("\t")
            length = int(length_text)
            hap, chrom = target_meta(seq)
            seq_by_key[(hap, chrom)] = seq
            length_by_key[(hap, chrom)] = length
            length_by_seq[seq] = length
    return seq_by_key, length_by_key, length_by_seq


def layout_for_lengths(label: str, lengths: dict[str, int]) -> GenomeLayout:
    offsets: dict[str, int] = {}
    offset = 0
    for chrom in CHROM_ORDER:
        length = lengths.get(chrom)
        if length is None:
            continue
        offsets[chrom] = offset
        offset += length
    return GenomeLayout(label=label, offsets=offsets, lengths=lengths, total=offset)


def target_layout(hap: str, length_by_key: dict[tuple[str, str], int]) -> GenomeLayout:
    lengths = {
        chrom: length_by_key[(hap, chrom)]
        for chrom in CHROM_ORDER
        if (hap, chrom) in length_by_key
    }
    label = f"PAN011 father {hap}"
    return layout_for_lengths(label, lengths)


def x_for(layout: GenomeLayout, chrom: str, pos: int) -> float:
    offset = layout.offsets[chrom] + max(0, min(pos, layout.lengths[chrom]))
    return TRACK_X0 + offset / layout.total * TRACK_W


def interval_x_native(layout: GenomeLayout, chrom: str, start: int, end: int) -> tuple[float, float]:
    x0 = x_for(layout, chrom, start)
    x1 = x_for(layout, chrom, end)
    if x1 < x0:
        x0, x1 = x1, x0
    return x0, x1


def interval_x(layout: GenomeLayout, chrom: str, start: int, end: int) -> tuple[float, float]:
    x0, x1 = interval_x_native(layout, chrom, start, end)
    if x1 - x0 < 3.0:
        mid = (x0 + x1) / 2
        x0, x1 = mid - 1.5, mid + 1.5
    return x0, x1


def donor_interval(row: dict[str, str]) -> tuple[str, int, int]:
    other_seq = row["other_seq"]
    for field in ("group.a", "group.b"):
        seq, start, end = parse_loc(row[field])
        if seq == other_seq:
            return seq, start, end
    seq, start, end = parse_loc(row["group.a"])
    return seq, start, end


def same_chrom_interval(row: dict[str, str]) -> tuple[str, int, int]:
    other_seq = row["other_seq"]
    for field in ("group.a", "group.b"):
        seq, start, end = parse_loc(row[field])
        if seq == other_seq:
            return seq, start, end
    return parse_loc(row["group.a"])


def read_segments(path: Path) -> list[Segment]:
    grouped: dict[tuple[str, int, int], dict[str, dict[str, str]]] = defaultdict(dict)
    with gzip.open(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["chrom"], int(row["start"]), int(row["end"]))
            grouped[key][row["winner_class"]] = row

    segments: list[Segment] = []
    for (query_seq, start, end), rows in grouped.items():
        if "same_chrom" not in rows or "interchrom" not in rows:
            continue
        same = rows["same_chrom"]
        inter = rows["interchrom"]
        same_identity = float(same["estimated.identity"])
        inter_identity = float(inter["estimated.identity"])
        if inter_identity <= same_identity:
            continue
        donor_seq, donor_start, donor_end = donor_interval(inter)
        donor_hap, target_chrom = target_meta(donor_seq)
        segments.append(
            Segment(
                query_seq=query_seq,
                query_chrom=inter["query_chrom"],
                query_start=start,
                query_end=end,
                donor_seq=donor_seq,
                donor_haplotype=donor_hap,
                target_chrom=target_chrom,
                donor_start=donor_start,
                donor_end=donor_end,
                bp=end - start,
                same_identity=same_identity,
                inter_identity=inter_identity,
            )
        )
    return segments


def read_homolog_segments(path: Path) -> list[Segment]:
    segments: list[Segment] = []
    with gzip.open(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["winner_class"] != "same_chrom":
                continue
            identity = float(row["estimated.identity"])
            if identity < HOMOLOG_MIN_IDENTITY:
                continue
            donor_seq, donor_start, donor_end = same_chrom_interval(row)
            donor_hap, target_chrom = target_meta(donor_seq)
            if donor_hap not in {"h1", "h2"}:
                continue
            query_chrom = row["query_chrom"]
            if target_chrom != query_chrom:
                continue
            start = int(row["start"])
            end = int(row["end"])
            segments.append(
                Segment(
                    query_seq=row["chrom"],
                    query_chrom=query_chrom,
                    query_start=start,
                    query_end=end,
                    donor_seq=donor_seq,
                    donor_haplotype=donor_hap,
                    target_chrom=target_chrom,
                    donor_start=donor_start,
                    donor_end=donor_end,
                    bp=end - start,
                    same_identity=identity,
                    inter_identity=identity,
                )
            )
    return segments


def group_runs(segments: list[Segment]) -> list[Run]:
    runs: list[Run] = []
    current: Run | None = None
    current_key: tuple[str, str] | None = None
    run_counts: dict[str, int] = defaultdict(int)

    ordered = sorted(
        segments,
        key=lambda s: (CHROM_ORDER.index(s.query_chrom), s.query_start, s.donor_seq, s.donor_start),
    )
    for segment in ordered:
        key = (segment.query_seq, segment.donor_seq)
        same_run = (
            current is not None
            and current_key == key
            and segment.query_start <= current.query_end + 2_000
            and (
                abs(segment.donor_start - current.donor_end) <= 10_000
                or abs(current.donor_start - segment.donor_end) <= 10_000
            )
        )
        if same_run:
            current.query_end = max(current.query_end, segment.query_end)
            current.donor_start = min(current.donor_start, segment.donor_start)
            current.donor_end = max(current.donor_end, segment.donor_end)
            current.bp += segment.bp
            current.windows += 1
            current.weighted_same_identity += segment.same_identity * segment.bp
            current.weighted_inter_identity += segment.inter_identity * segment.bp
            continue

        panel_id = segment.query_chrom
        run_counts[panel_id] += 1
        current = Run(
            run_id=f"{panel_id}_run{run_counts[panel_id]:04d}",
            query_seq=segment.query_seq,
            query_chrom=segment.query_chrom,
            query_start=segment.query_start,
            query_end=segment.query_end,
            donor_seq=segment.donor_seq,
            donor_haplotype=segment.donor_haplotype,
            target_chrom=segment.target_chrom,
            donor_start=segment.donor_start,
            donor_end=segment.donor_end,
            bp=segment.bp,
            windows=1,
            weighted_same_identity=segment.same_identity * segment.bp,
            weighted_inter_identity=segment.inter_identity * segment.bp,
        )
        runs.append(current)
        current_key = key
    return runs


def group_homolog_runs(segments: list[Segment]) -> list[Run]:
    runs: list[Run] = []
    current: Run | None = None
    current_key: tuple[str, str] | None = None
    run_counts: dict[str, int] = defaultdict(int)

    ordered = sorted(
        segments,
        key=lambda s: (CHROM_ORDER.index(s.query_chrom), s.query_start, s.donor_seq, s.donor_start),
    )
    for segment in ordered:
        key = (segment.query_seq, segment.donor_seq)
        same_run = (
            current is not None
            and current_key == key
            and segment.query_start <= current.query_end + 10_000
            and (
                abs(segment.donor_start - current.donor_end) <= 50_000
                or abs(current.donor_start - segment.donor_end) <= 50_000
            )
        )
        if same_run:
            current.query_end = max(current.query_end, segment.query_end)
            current.donor_start = min(current.donor_start, segment.donor_start)
            current.donor_end = max(current.donor_end, segment.donor_end)
            current.bp += segment.bp
            current.windows += 1
            current.weighted_same_identity += segment.same_identity * segment.bp
            current.weighted_inter_identity += segment.inter_identity * segment.bp
            continue

        run_counts[segment.query_chrom] += 1
        current = Run(
            run_id=f"{segment.query_chrom}_homolog_run{run_counts[segment.query_chrom]:04d}",
            query_seq=segment.query_seq,
            query_chrom=segment.query_chrom,
            query_start=segment.query_start,
            query_end=segment.query_end,
            donor_seq=segment.donor_seq,
            donor_haplotype=segment.donor_haplotype,
            target_chrom=segment.target_chrom,
            donor_start=segment.donor_start,
            donor_end=segment.donor_end,
            bp=segment.bp,
            windows=1,
            weighted_same_identity=segment.same_identity * segment.bp,
            weighted_inter_identity=segment.inter_identity * segment.bp,
        )
        runs.append(current)
        current_key = key
    return runs


def high_confidence_runs(runs: list[Run]) -> list[Run]:
    out = [
        run
        for run in runs
        if run.bp >= 10_000 and run.mean_inter_identity >= 0.95 and run.donor_haplotype in {"h1", "h2"}
    ]
    return sorted(
        out,
        key=lambda r: (
            {"acro_acro": 0, "acro_other": 1, "other_nonacro": 2, "chr5_chr1_candidate": 3, "chr9_chr3_candidate": 4, "PAR_XY": 5}[r.category],
            r.bp,
        ),
    )


def high_confidence_homolog_runs(runs: list[Run]) -> list[Run]:
    return sorted(
        [run for run in runs if run.bp >= HOMOLOG_MIN_BP and run.donor_haplotype in {"h1", "h2"}],
        key=lambda r: (-r.bp, CHROM_ORDER.index(r.query_chrom), r.query_start, r.donor_haplotype),
    )


def write_runs(path: Path, runs: list[Run]) -> None:
    fields = [
        "run_id",
        "query_chrom",
        "query_start",
        "query_end",
        "target_chrom",
        "donor_haplotype",
        "donor_seq",
        "donor_start",
        "donor_end",
        "bp",
        "windows",
        "mean_same_identity",
        "mean_inter_identity",
        "delta_identity",
        "category",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for run in sorted(runs, key=lambda r: (-r.bp, r.query_chrom, r.query_start, r.target_chrom)):
            writer.writerow(
                {
                    "run_id": run.run_id,
                    "query_chrom": run.query_chrom,
                    "query_start": run.query_start,
                    "query_end": run.query_end,
                    "target_chrom": run.target_chrom,
                    "donor_haplotype": run.donor_haplotype,
                    "donor_seq": run.donor_seq,
                    "donor_start": run.donor_start,
                    "donor_end": run.donor_end,
                    "bp": run.bp,
                    "windows": run.windows,
                    "mean_same_identity": f"{run.mean_same_identity:.6f}",
                    "mean_inter_identity": f"{run.mean_inter_identity:.6f}",
                    "delta_identity": f"{run.mean_inter_identity - run.mean_same_identity:.6f}",
                    "category": run.category,
                }
            )


def write_homolog_runs(path: Path, runs: list[Run]) -> None:
    fields = [
        "run_id",
        "query_chrom",
        "query_start",
        "query_end",
        "target_chrom",
        "donor_haplotype",
        "donor_seq",
        "donor_start",
        "donor_end",
        "bp",
        "windows",
        "mean_identity",
        "category",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for run in runs:
            writer.writerow(
                {
                    "run_id": run.run_id,
                    "query_chrom": run.query_chrom,
                    "query_start": run.query_start,
                    "query_end": run.query_end,
                    "target_chrom": run.target_chrom,
                    "donor_haplotype": run.donor_haplotype,
                    "donor_seq": run.donor_seq,
                    "donor_start": run.donor_start,
                    "donor_end": run.donor_end,
                    "bp": run.bp,
                    "windows": run.windows,
                    "mean_identity": f"{run.mean_same_identity:.6f}",
                    "category": "homologous_same_chrom",
                }
            )


def write_summary(path: Path, runs: list[Run], all_runs: list[Run]) -> None:
    by_category: dict[str, tuple[int, int]] = {}
    for category in COLORS:
        selected = [run for run in runs if run.category == category]
        by_category[category] = (len(selected), sum(run.bp for run in selected))
    with path.open("w", newline="") as handle:
        fields = ["metric", "value"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerow({"metric": "source", "value": str(CLASS_WINNERS)})
        writer.writerow({"metric": "all_inter_beats_same_runs", "value": str(len(all_runs))})
        writer.writerow({"metric": "drawn_high_conf_runs", "value": str(len(runs))})
        writer.writerow({"metric": "drawn_high_conf_bp", "value": str(sum(run.bp for run in runs))})
        for category, (count, bp) in by_category.items():
            writer.writerow({"metric": f"{category}_runs", "value": str(count)})
            writer.writerow({"metric": f"{category}_bp", "value": str(bp)})


def write_homolog_summary(path: Path, homolog_runs: list[Run], all_homolog_runs: list[Run], inter_runs: list[Run]) -> None:
    with path.open("w", newline="") as handle:
        fields = ["metric", "value"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerow({"metric": "source", "value": str(CLASS_WINNERS)})
        writer.writerow({"metric": "homolog_layer_definition", "value": "grouped same_chrom father-child donor intervals drawn to child intervals from the 10:10 IMPG class-winner table"})
        writer.writerow({"metric": "homolog_min_identity", "value": str(HOMOLOG_MIN_IDENTITY)})
        writer.writerow({"metric": "homolog_min_bp", "value": str(HOMOLOG_MIN_BP)})
        writer.writerow({"metric": "all_homolog_runs", "value": str(len(all_homolog_runs))})
        writer.writerow({"metric": "drawn_homolog_runs", "value": str(len(homolog_runs))})
        writer.writerow({"metric": "drawn_homolog_bp", "value": str(sum(run.bp for run in homolog_runs))})
        writer.writerow({"metric": "drawn_interchrom_runs", "value": str(len(inter_runs))})
        writer.writerow({"metric": "drawn_interchrom_bp", "value": str(sum(run.bp for run in inter_runs))})


def ribbon_path(xa0: float, xa1: float, ya: float, xb0: float, xb1: float, yb: float) -> str:
    c = abs(yb - ya) * 0.48
    return (
        f"M {xa0:.2f} {ya:.2f} "
        f"C {xa0:.2f} {ya + c:.2f}, {xb0:.2f} {yb - c:.2f}, {xb0:.2f} {yb:.2f} "
        f"L {xb1:.2f} {yb:.2f} "
        f"C {xb1:.2f} {yb - c:.2f}, {xa1:.2f} {ya + c:.2f}, {xa1:.2f} {ya:.2f} Z"
    )


def draw_genome_track(svg: SVG, layout: GenomeLayout, y: float) -> None:
    svg.text(TRACK_X0 - 26, y + TRACK_H - 2, layout.label, 30, "700", TEXT, "end")
    for idx, chrom in enumerate(CHROM_ORDER):
        if chrom not in layout.lengths:
            continue
        x0 = x_for(layout, chrom, 0)
        x1 = x_for(layout, chrom, layout.lengths[chrom])
        fill = "#eef0f2" if idx % 2 == 0 else "#dfe3e6"
        svg.rect(x0, y, x1 - x0, TRACK_H, fill, "#ffffff", 1.1, 1.0, rx=0)
        if x1 - x0 >= 38:
            label = chrom.replace("chr", "")
            svg.text((x0 + x1) / 2, y - 9, label, 19, "600", MUTED, "middle")
        if chrom in {"chr1", "chr3", "chr5", "chr9", "chr13", "chr15", "chr21", "chr22", "chrX", "chrY"}:
            svg.line(x0, GRID_Y0, x0, GRID_Y1, GRID, 0.7, 0.55)
    svg.rect(TRACK_X0, y, TRACK_W, TRACK_H, "none", "#bdc1c6", 1.1, 1.0, rx=0)


def ribbon_opacity(category: str) -> float:
    if category in {"PAR_XY", "chr5_chr1_candidate", "chr9_chr3_candidate"}:
        return 0.42
    if category == "acro_acro":
        return 0.16
    return 0.18


def draw_interval(svg: SVG, x0: float, x1: float, y: float, color: str, opacity: float = 0.95) -> None:
    svg.rect(x0, y - 5, x1 - x0, TRACK_H + 10, color, "none", 0, opacity, rx=1.5)


def ribbon_y_for(source_y: float, target_y: float) -> tuple[float, float]:
    if source_y < target_y:
        return source_y + TRACK_H + 8, target_y - 8
    return source_y - 8, target_y + TRACK_H + 8


def label_text(run: Run) -> str:
    if run.category == "PAR_XY":
        return "PAR1 chrX ~ chrY"
    if run.category == "chr5_chr1_candidate":
        return "chr5q ~ chr1p"
    if run.category == "chr9_chr3_candidate":
        return "chr9q ~ chr3q"
    return f"{run.query_chrom} ~ {run.target_chrom}"


def draw_callouts(svg: SVG, runs: list[Run], query_layout: GenomeLayout) -> None:
    selected: list[Run] = []
    for category in ("PAR_XY", "chr5_chr1_candidate", "chr9_chr3_candidate"):
        rows = [run for run in runs if run.category == category]
        rows = sorted(rows, key=lambda r: -r.bp)
        selected.extend(rows[:1])

    lanes: list[float] = []
    for run in sorted(selected, key=lambda r: x_for(query_layout, r.query_chrom, (r.query_start + r.query_end) // 2)):
        x0, x1 = interval_x(query_layout, run.query_chrom, run.query_start, run.query_end)
        center = (x0 + x1) / 2
        width = max(170, len(label_text(run)) * 15)
        lane = 0
        for idx, last_right in enumerate(lanes):
            if center - width / 2 > last_right + 24:
                lane = idx
                lanes[idx] = center + width / 2
                break
        else:
            lane = len(lanes)
            lanes.append(center + width / 2)
        y = Y_QUERY - 100 - lane * 42
        color = COLORS[run.category]
        svg.line(center, Y_QUERY - 9, center, y + 8, color, 1.4, 0.9)
        svg.text(center, y, label_text(run), 27, "700", color, "middle")
        svg.text(center, y + 29, f"{run.bp / 1000:.0f} kb, {run.donor_haplotype}", 20, "400", MUTED, "middle")


def draw_legend(svg: SVG, runs: list[Run]) -> None:
    x = TRACK_X0
    y = LEGEND_Y
    svg.text(x, y, "Ribbon classes", 24, "700", MUTED)
    x += 210
    for category, label in [
        ("PAR_XY", "PAR1 positive control"),
        ("chr5_chr1_candidate", "chr5q/chr1p candidate"),
        ("chr9_chr3_candidate", "chr9q/chr3q candidate"),
        ("acro_acro", "acrocentric"),
        ("acro_other", "acro-other"),
        ("other_nonacro", "other"),
    ]:
        count = sum(1 for run in runs if run.category == category)
        if count == 0:
            continue
        svg.rect(x, y - 20, 27, 21, COLORS[category], "none", 0, 0.9)
        svg.text(x + 39, y, f"{label} ({count})", 21, "400", TEXT)
        x += 390 if category != "other_nonacro" else 220


def draw_homolog_legend(svg: SVG, runs: list[Run], homolog_runs: list[Run]) -> None:
    x = TRACK_X0
    y = LEGEND_Y
    svg.text(x, y, "Ribbon classes", 24, "700", MUTED)
    x += 210
    svg.rect(x, y - 20, 27, 21, HOMOLOG_COLOR, "none", 0, 0.9)
    svg.text(x + 39, y, f"same chromosome homologous ({len(homolog_runs)})", 21, "400", TEXT)
    x += 560
    for category, label in [
        ("PAR_XY", "PAR1 non-homologous"),
        ("chr5_chr1_candidate", "chr5q/chr1p"),
        ("chr9_chr3_candidate", "chr9q/chr3q"),
        ("acro_acro", "acrocentric"),
    ]:
        count = sum(1 for run in runs if run.category == category)
        if count == 0:
            continue
        svg.rect(x, y - 20, 27, 21, COLORS[category], "none", 0, 0.9)
        svg.text(x + 39, y, f"{label} ({count})", 21, "400", TEXT)
        x += 360


def render(runs: list[Run], query_layout: GenomeLayout, hap1_layout: GenomeLayout, hap2_layout: GenomeLayout) -> None:
    svg = SVG(PAGE_W, PAGE_H)
    svg.text(TRACK_X0, 58, "Whole-genome 10:10 SweepGA/F32 IMPG class winners with donor ribbons", 42, "700", TEXT)
    svg.text(
        TRACK_X0,
        98,
        "PAN027 paternal child query; ribbons mark high-confidence 2 kb runs where the best interchromosomal match beats the best same-chromosome match.",
        24,
        "400",
        MUTED,
    )

    for y, label in [
        (Y_H1 - 16, "father donor h1"),
        (Y_QUERY + TRACK_H + 8, "child query"),
        (Y_H2 - 16, "father donor h2"),
    ]:
        svg.line(TRACK_X0, y, TRACK_X0 + TRACK_W, y, "#f1f3f4", 0.8, 0.8)
        svg.text(TRACK_X0 + TRACK_W + 22, y + 7, label, 20, "400", MUTED)

    draw_genome_track(svg, query_layout, Y_QUERY)
    draw_genome_track(svg, hap1_layout, Y_H1)
    draw_genome_track(svg, hap2_layout, Y_H2)

    target_layouts = {"h1": hap1_layout, "h2": hap2_layout}
    target_y = {"h1": Y_H1, "h2": Y_H2}

    for run in runs:
        target_layout_obj = target_layouts.get(run.donor_haplotype)
        if target_layout_obj is None or run.target_chrom not in target_layout_obj.lengths:
            continue
        qx0, qx1 = interval_x_native(query_layout, run.query_chrom, run.query_start, run.query_end)
        dx0, dx1 = interval_x_native(target_layout_obj, run.target_chrom, run.donor_start, run.donor_end)
        color = COLORS[run.category]
        donor_edge_y, child_edge_y = ribbon_y_for(target_y[run.donor_haplotype], Y_QUERY)
        d = ribbon_path(dx0, dx1, donor_edge_y, qx0, qx1, child_edge_y)
        svg.path(d, color, "none", 0, ribbon_opacity(run.category))

    for run in runs:
        target_layout_obj = target_layouts.get(run.donor_haplotype)
        if target_layout_obj is None or run.target_chrom not in target_layout_obj.lengths:
            continue
        color = COLORS[run.category]
        qx0, qx1 = interval_x(query_layout, run.query_chrom, run.query_start, run.query_end)
        dx0, dx1 = interval_x(target_layout_obj, run.target_chrom, run.donor_start, run.donor_end)
        emph = 1.0 if run.category in {"PAR_XY", "chr5_chr1_candidate", "chr9_chr3_candidate"} else 0.55
        draw_interval(svg, qx0, qx1, Y_QUERY, color, 0.88 * emph)
        draw_interval(svg, dx0, dx1, target_y[run.donor_haplotype], color, 0.88 * emph)

    draw_callouts(svg, runs, query_layout)
    draw_legend(svg, runs)

    total_bp = sum(run.bp for run in runs)
    svg.text(
        TRACK_X0,
        FOOTNOTE_Y1,
        f"Drawn: {len(runs)} high-confidence runs, {total_bp / 1e6:.2f} Mb total; threshold: run >=10 kb and mean interchrom identity >=0.95.",
        19,
        "400",
        MUTED,
    )
    svg.text(
        TRACK_X0,
        FOOTNOTE_Y2,
        "All coordinates are native whole-genome coordinates collapsed into length-scaled chromosome-order tracks.",
        19,
        "400",
        MUTED,
    )
    svg.write(SVG_OUT)


def render_homolog_context(
    inter_runs: list[Run],
    homolog_runs: list[Run],
    query_layout: GenomeLayout,
    hap1_layout: GenomeLayout,
    hap2_layout: GenomeLayout,
) -> None:
    svg = SVG(PAGE_W, PAGE_H)
    svg.text(TRACK_X0, 58, "Whole-genome homologous inheritance context with non-homologous winners", 42, "700", TEXT)
    svg.text(
        TRACK_X0,
        98,
        "Light-gray ribbons are full same-chromosome father-child homologous chains; colored ribbons are interchromosomal winners from the same 10:10 IMPG scan.",
        24,
        "400",
        MUTED,
    )

    for y, label in [
        (Y_H1 - 16, "father donor h1"),
        (Y_QUERY + TRACK_H + 8, "child query"),
        (Y_H2 - 16, "father donor h2"),
    ]:
        svg.line(TRACK_X0, y, TRACK_X0 + TRACK_W, y, "#f1f3f4", 0.8, 0.8)
        svg.text(TRACK_X0 + TRACK_W + 22, y + 7, label, 20, "400", MUTED)

    draw_genome_track(svg, query_layout, Y_QUERY)
    draw_genome_track(svg, hap1_layout, Y_H1)
    draw_genome_track(svg, hap2_layout, Y_H2)

    target_layouts = {"h1": hap1_layout, "h2": hap2_layout}
    target_y = {"h1": Y_H1, "h2": Y_H2}

    for run in homolog_runs:
        target_layout_obj = target_layouts.get(run.donor_haplotype)
        if target_layout_obj is None or run.target_chrom not in target_layout_obj.lengths:
            continue
        qx0, qx1 = interval_x_native(query_layout, run.query_chrom, run.query_start, run.query_end)
        dx0, dx1 = interval_x_native(target_layout_obj, run.target_chrom, run.donor_start, run.donor_end)
        donor_edge_y, child_edge_y = ribbon_y_for(target_y[run.donor_haplotype], Y_QUERY)
        d = ribbon_path(dx0, dx1, donor_edge_y, qx0, qx1, child_edge_y)
        svg.path(d, HOMOLOG_RIBBON, "none", 0, 0.13)

    for run in homolog_runs:
        target_layout_obj = target_layouts.get(run.donor_haplotype)
        if target_layout_obj is None or run.target_chrom not in target_layout_obj.lengths:
            continue
        qx0, qx1 = interval_x(query_layout, run.query_chrom, run.query_start, run.query_end)
        dx0, dx1 = interval_x(target_layout_obj, run.target_chrom, run.donor_start, run.donor_end)
        draw_interval(svg, qx0, qx1, Y_QUERY, HOMOLOG_COLOR, 0.22)
        draw_interval(svg, dx0, dx1, target_y[run.donor_haplotype], HOMOLOG_COLOR, 0.22)

    for run in inter_runs:
        target_layout_obj = target_layouts.get(run.donor_haplotype)
        if target_layout_obj is None or run.target_chrom not in target_layout_obj.lengths:
            continue
        qx0, qx1 = interval_x_native(query_layout, run.query_chrom, run.query_start, run.query_end)
        dx0, dx1 = interval_x_native(target_layout_obj, run.target_chrom, run.donor_start, run.donor_end)
        color = COLORS[run.category]
        donor_edge_y, child_edge_y = ribbon_y_for(target_y[run.donor_haplotype], Y_QUERY)
        d = ribbon_path(dx0, dx1, donor_edge_y, qx0, qx1, child_edge_y)
        opacity = 0.55 if run.category in {"PAR_XY", "chr5_chr1_candidate", "chr9_chr3_candidate"} else 0.28
        svg.path(d, color, "none", 0, opacity)

    for run in inter_runs:
        target_layout_obj = target_layouts.get(run.donor_haplotype)
        if target_layout_obj is None or run.target_chrom not in target_layout_obj.lengths:
            continue
        color = COLORS[run.category]
        qx0, qx1 = interval_x(query_layout, run.query_chrom, run.query_start, run.query_end)
        dx0, dx1 = interval_x(target_layout_obj, run.target_chrom, run.donor_start, run.donor_end)
        emph = 1.0 if run.category in {"PAR_XY", "chr5_chr1_candidate", "chr9_chr3_candidate"} else 0.65
        draw_interval(svg, qx0, qx1, Y_QUERY, color, 0.95 * emph)
        draw_interval(svg, dx0, dx1, target_y[run.donor_haplotype], color, 0.95 * emph)

    draw_callouts(svg, inter_runs, query_layout)
    draw_homolog_legend(svg, inter_runs, homolog_runs)

    svg.text(
        TRACK_X0,
        FOOTNOTE_Y1,
        f"Light gray: {len(homolog_runs)} grouped same-chromosome chains >=50 kb drawn as donor-interval to child-interval ribbons.",
        19,
        "400",
        MUTED,
    )
    svg.text(
        TRACK_X0,
        FOOTNOTE_Y2,
        "All tracks use native whole-genome coordinates collapsed into length-scaled chromosome-order tracks.",
        19,
        "400",
        MUTED,
    )
    svg.write(HOMOLOG_SVG_OUT)


def find_rsvg() -> str | None:
    found = shutil.which("rsvg-convert")
    if found:
        return found
    if os.environ.get("GUIX_ENVIRONMENT"):
        candidate = Path(os.environ["GUIX_ENVIRONMENT"]) / "bin/rsvg-convert"
        if candidate.exists():
            return str(candidate)
    for candidate in [
        Path("/gnu/store/kawnjzdr5wi6x77psqsvmaqqni359df5-profile/bin/rsvg-convert"),
        Path("/gnu/store/bb8ijpv1y3wpppfqd7r0pkk25xckag19-librsvg-2.54.5/bin/rsvg-convert"),
        Path("/gnu/store/42q62mvxp7hhixhavpchm7pzjrl29630-librsvg-2.58.5/bin/rsvg-convert"),
    ]:
        if candidate.exists():
            return str(candidate)
    return None


def convert_one(svg_path: Path) -> str | None:
    rsvg = find_rsvg()
    if rsvg is None:
        return None
    subprocess.run([rsvg, "-f", "pdf", "-o", str(svg_path.with_suffix(".pdf")), str(svg_path)], check=True)
    subprocess.run([rsvg, "-f", "png", "-o", str(svg_path.with_suffix(".png")), str(svg_path)], check=True)
    return subprocess.run([rsvg, "--version"], check=True, text=True, capture_output=True).stdout.strip()


def convert_outputs() -> list[str]:
    messages: list[str] = []
    version = convert_one(SVG_OUT)
    if version is None:
        return ["SVG only: no rsvg-convert found."]
    homolog_version = convert_one(HOMOLOG_SVG_OUT)
    messages.append(f"converted PDF and PNG with {version}")
    if homolog_version is not None:
        messages.append(f"converted homologous-context PDF and PNG with {homolog_version}")
    return messages


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    _query_seqs, query_lengths = read_query_fai(QUERY_FAI)
    _target_seqs, target_lengths_by_key, _target_lengths_by_seq = read_target_fai(TARGET_FAI)
    query_layout = layout_for_lengths("PAN027 pat child", query_lengths)
    hap1_layout = target_layout("h1", target_lengths_by_key)
    hap2_layout = target_layout("h2", target_lengths_by_key)

    segments = read_segments(CLASS_WINNERS)
    all_runs = group_runs(segments)
    runs = high_confidence_runs(all_runs)
    homolog_segments = read_homolog_segments(CLASS_WINNERS)
    all_homolog_runs = group_homolog_runs(homolog_segments)
    homolog_runs = high_confidence_homolog_runs(all_homolog_runs)
    write_runs(RUNS_OUT, runs)
    write_summary(SUMMARY_OUT, runs, all_runs)
    write_homolog_runs(HOMOLOG_RUNS_OUT, homolog_runs)
    write_homolog_summary(HOMOLOG_SUMMARY_OUT, homolog_runs, all_homolog_runs, runs)
    render(runs, query_layout, hap1_layout, hap2_layout)
    render_homolog_context(runs, homolog_runs, query_layout, hap1_layout, hap2_layout)
    messages = convert_outputs()
    CONVERSION_STATUS.write_text("\n".join(messages) + "\n")
    for message in messages:
        print(message)
    print(
        f"segments={len(segments)} all_runs={len(all_runs)} drawn_runs={len(runs)} "
        f"homolog_segments={len(homolog_segments)} homolog_runs={len(all_homolog_runs)} "
        f"drawn_homolog_runs={len(homolog_runs)}"
    )
    print(f"wrote {SVG_OUT}")
    print(f"wrote {RUNS_OUT}")


if __name__ == "__main__":
    main()
