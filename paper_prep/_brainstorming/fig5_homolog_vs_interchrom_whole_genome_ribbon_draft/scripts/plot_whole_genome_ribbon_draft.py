#!/usr/bin/env python3
"""Whole-genome ribbon renderer for interchrom-over-homolog class-winner calls."""

from __future__ import annotations

import argparse
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
DEFAULT_OUT_DIR = ROOT / "paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft"
CLASS_WINNERS_NAME = (
    "PAN027pat_vs_PAN011_joint.sweepga_f32.10to10.query_2000bp."
    "predepth_class_winners.impg_similarity.tsv.gz"
)
def _first_existing(paths: list[Path]) -> Path:
    for p in paths:
        if p.exists():
            return p
    return paths[-1]

# Prefer the repo-vendored copies in data/ (self-contained; no moosefs needed).
# Fall back to the original moosefs sources when absent. Override on the CLI with
# --class-winners / --query-fai / --target-fai.
DATA_CLASS_WINNERS = ROOT / "data/fig5_PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz"
LOCAL_CLASS_WINNERS = (
    ROOT
    / "paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs/"
    / CLASS_WINNERS_NAME
)
MOOSEFS_CLASS_WINNERS = (
    Path("/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs")
    / CLASS_WINNERS_NAME
)
CLASS_WINNERS = _first_existing([DATA_CLASS_WINNERS, LOCAL_CLASS_WINNERS, MOOSEFS_CLASS_WINNERS])
DEFAULT_QUERY_FAI = _first_existing([
    ROOT / "data/fig5_PAN027pat_vs_PAN011_joint.query.fa.fai",
    Path("/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/"
         "pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027pat_vs_PAN011_joint.query.fa.fai"),
])
DEFAULT_TARGET_FAI = _first_existing([
    ROOT / "data/fig5_PAN027pat_vs_PAN011_joint.target.fa.fai",
    Path("/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/"
         "pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027pat_vs_PAN011_joint.target.fa.fai"),
])

CHROM_ORDER = [f"chr{i}" for i in range(1, 22 + 1)] + ["chrX", "chrY"]
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}
CHR_RE = re.compile(r"(chr(?:[0-9]+|X|Y|M))")
TARGET_RE = re.compile(r".+#joint#(?P<hap>h[12])_(?P<chrom>chr(?:[0-9]+|X|Y|M))")
LOC_RE = re.compile(r"^(?P<seq>.+):(?P<start>[0-9]+)-(?P<end>[0-9]+)$")

PAGE_W = 3600
PAGE_H = 840
TRACK_X0 = 850
TRACK_W = 2580
TRACK_H = 28
TRACK_LABEL_X = TRACK_X0 - 36
Y_H1 = 115
Y_QUERY = 345
Y_H2 = 555
GRID_Y0 = 25
GRID_Y1 = 625
LEGEND_Y = 690
FOOTNOTE_Y1 = 955
FOOTNOTE_Y2 = 990
TEXT = "#202124"
MUTED = "#5f6368"
FONT_SCALE = 1.9   # global multiplier on every text size (bump to enlarge fonts)
GRID = "#e8eaed"
CHROM_BORDER = "#111111"
HOMOLOG_COLOR = "#b8bdc3"
HOMOLOG_RIBBON = "#cfd3d7"
HOMOLOG_RIBBON_OPACITY = 0.07
HOMOLOG_MIN_BP = 10_000
HOMOLOG_MIN_IDENTITY = 0.95
CONTIGUOUS_MERGE_GAP_BP = 0
INTERCHROM_RIBBON_MIN_W = 7.0

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


@dataclass
class RenderConfig:
    comparison_id: str
    class_winners: Path
    query_fai: Path
    target_fai: Path
    out_dir: Path
    query_label: str
    target_h1_label: str
    target_h2_label: str
    layer_label: str
    svg_out: Path
    homolog_svg_out: Path
    runs_out: Path
    summary_out: Path
    homolog_runs_out: Path
    homolog_summary_out: Path
    merge_audit_out: Path
    conversion_status: Path


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
            f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size * FONT_SCALE:.1f}" '
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


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path)


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


def target_layout(hap: str, length_by_key: dict[tuple[str, str], int], label: str) -> GenomeLayout:
    lengths = {
        chrom: length_by_key[(hap, chrom)]
        for chrom in CHROM_ORDER
        if (hap, chrom) in length_by_key
    }
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


def interval_x_with_min_width(layout: GenomeLayout, chrom: str, start: int, end: int, min_w: float) -> tuple[float, float]:
    x0, x1 = interval_x_native(layout, chrom, start, end)
    if x1 - x0 < min_w:
        mid = (x0 + x1) / 2
        half = min_w / 2
        x0, x1 = mid - half, mid + half
    return x0, x1


def homolog_visual_width(bp: int) -> float:
    return max(6.0, min(34.0, 4.0 + bp / 90_000.0))


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


def touches(left: int, right: int) -> bool:
    return abs(left - right) <= CONTIGUOUS_MERGE_GAP_BP


def renumber_runs(runs: list[Run], run_name: str) -> list[Run]:
    counts: dict[str, int] = defaultdict(int)
    ordered = sorted(
        runs,
        key=lambda r: (
            CHROM_ORDER.index(r.query_chrom),
            r.query_start,
            r.query_end,
            r.donor_haplotype,
            r.target_chrom,
            r.donor_start,
            r.donor_end,
        ),
    )
    for run in ordered:
        counts[run.query_chrom] += 1
        run.run_id = f"{run.query_chrom}_{run_name}{counts[run.query_chrom]:04d}"
    return ordered


def donor_step_direction(previous: Segment, segment: Segment, current_direction: int | None) -> int | None:
    if not touches(segment.query_start, previous.query_end):
        return None
    if current_direction in {None, 1} and touches(segment.donor_start, previous.donor_end):
        return 1
    if current_direction in {None, -1} and touches(segment.donor_end, previous.donor_start):
        return -1
    return None


def start_run(segment: Segment, run_name: str, run_counts: dict[str, int]) -> Run:
    run_counts[segment.query_chrom] += 1
    donor_start, donor_end = sorted((segment.donor_start, segment.donor_end))
    return Run(
        run_id=f"{segment.query_chrom}_{run_name}{run_counts[segment.query_chrom]:04d}",
        query_seq=segment.query_seq,
        query_chrom=segment.query_chrom,
        query_start=segment.query_start,
        query_end=segment.query_end,
        donor_seq=segment.donor_seq,
        donor_haplotype=segment.donor_haplotype,
        target_chrom=segment.target_chrom,
        donor_start=donor_start,
        donor_end=donor_end,
        bp=segment.bp,
        windows=1,
        weighted_same_identity=segment.same_identity * segment.bp,
        weighted_inter_identity=segment.inter_identity * segment.bp,
    )


def extend_run(run: Run, segment: Segment) -> None:
    donor_start, donor_end = sorted((segment.donor_start, segment.donor_end))
    run.query_start = min(run.query_start, segment.query_start)
    run.query_end = max(run.query_end, segment.query_end)
    run.donor_start = min(run.donor_start, donor_start)
    run.donor_end = max(run.donor_end, donor_end)
    run.bp += segment.bp
    run.windows += 1
    run.weighted_same_identity += segment.same_identity * segment.bp
    run.weighted_inter_identity += segment.inter_identity * segment.bp


def group_end_to_end_runs(segments: list[Segment], run_name: str) -> list[Run]:
    buckets: dict[tuple[str, str], list[Segment]] = defaultdict(list)
    for segment in segments:
        buckets[(segment.query_seq, segment.donor_seq)].append(segment)

    runs: list[Run] = []
    run_counts: dict[str, int] = defaultdict(int)
    for _key, key_segments in sorted(
        buckets.items(),
        key=lambda item: (
            CHROM_ORDER.index(item[1][0].query_chrom),
            min(segment.query_start for segment in item[1]),
            item[0][1],
        ),
    ):
        current: Run | None = None
        previous: Segment | None = None
        direction: int | None = None
        for segment in sorted(
            key_segments,
            key=lambda s: (CHROM_ORDER.index(s.query_chrom), s.query_start, s.query_end, s.donor_start, s.donor_end),
        ):
            next_direction = donor_step_direction(previous, segment, direction) if previous is not None else None
            if current is not None and next_direction is not None:
                direction = next_direction if direction is None else direction
                extend_run(current, segment)
            else:
                current = start_run(segment, run_name, run_counts)
                runs.append(current)
                direction = None
            previous = segment
    return renumber_runs(runs, run_name)


def group_runs(segments: list[Segment]) -> list[Run]:
    return group_end_to_end_runs(segments, "run")


def group_homolog_runs(segments: list[Segment]) -> list[Run]:
    return group_end_to_end_runs(segments, "homolog_run")


def query_span(run: Run) -> int:
    return run.query_end - run.query_start


def donor_span(run: Run) -> int:
    return abs(run.donor_end - run.donor_start)


def validate_run_spans(runs: list[Run], label: str, allow_query_gaps: bool = False) -> None:
    first_window_only = [
        run
        for run in runs
        if run.windows > 1 and query_span(run) <= max(2000, CONTIGUOUS_MERGE_GAP_BP)
    ]
    query_span_invalid = [
        run for run in runs if query_span(run) < run.bp or (not allow_query_gaps and query_span(run) != run.bp)
    ]
    if first_window_only or query_span_invalid:
        example = (first_window_only or query_span_invalid)[0]
        raise AssertionError(
            f"{label} merged span validation failed for {example.run_id}: "
            f"windows={example.windows} bp={example.bp} "
            f"query={example.query_start}-{example.query_end} "
            f"donor={example.donor_start}-{example.donor_end}"
        )


def span_audit_metrics(runs: list[Run]) -> dict[str, int]:
    return {
        "multi_window_runs_with_first_window_only_query_span": sum(
            1 for run in runs if run.windows > 1 and query_span(run) <= max(2000, CONTIGUOUS_MERGE_GAP_BP)
        ),
        "query_span_shorter_than_bp_runs": sum(1 for run in runs if query_span(run) < run.bp),
        "query_span_longer_than_bp_runs": sum(1 for run in runs if query_span(run) > run.bp),
        "max_query_span_bp": max((query_span(run) for run in runs), default=0),
        "max_donor_span_bp": max((donor_span(run) for run in runs), default=0),
        "max_windows_per_run": max((run.windows for run in runs), default=0),
    }


def merge_audit_metrics(segments: list[Segment], runs: list[Run]) -> dict[str, int]:
    buckets: dict[tuple[str, str], list[Segment]] = defaultdict(list)
    for segment in segments:
        buckets[(segment.query_seq, segment.donor_seq)].append(segment)

    absorbed_fragments = 0
    max_query_gap = 0
    max_donor_gap = 0
    for _key, key_segments in buckets.items():
        previous: Segment | None = None
        direction: int | None = None
        for segment in sorted(
            key_segments,
            key=lambda s: (CHROM_ORDER.index(s.query_chrom), s.query_start, s.query_end, s.donor_start, s.donor_end),
        ):
            next_direction = donor_step_direction(previous, segment, direction) if previous is not None else None
            if previous is not None and next_direction is not None:
                max_query_gap = max(max_query_gap, abs(segment.query_start - previous.query_end))
                if next_direction == 1:
                    max_donor_gap = max(max_donor_gap, abs(segment.donor_start - previous.donor_end))
                else:
                    max_donor_gap = max(max_donor_gap, abs(segment.donor_end - previous.donor_start))
                direction = next_direction if direction is None else direction
                absorbed_fragments += 1
            else:
                direction = None
            previous = segment
    return {
        "raw_segments_2kb": len(segments),
        "end_to_end_runs": len(runs),
        "absorbed_fragments": absorbed_fragments,
        "multi_window_runs": sum(1 for run in runs if run.windows > 1),
        "max_absorbed_query_endpoint_gap_bp": max_query_gap,
        "max_absorbed_donor_endpoint_gap_bp": max_donor_gap,
        "merged_windows": sum(run.windows for run in runs),
    }


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
        "query_span_bp",
        "donor_span_bp",
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
                    "query_span_bp": query_span(run),
                    "donor_span_bp": donor_span(run),
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
        "query_span_bp",
        "donor_span_bp",
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
                    "query_span_bp": query_span(run),
                    "donor_span_bp": donor_span(run),
                    "windows": run.windows,
                    "mean_identity": f"{run.mean_same_identity:.6f}",
                    "category": "homologous_same_chrom",
                }
            )


def write_summary(path: Path, runs: list[Run], all_runs: list[Run], segments: list[Segment], source: Path) -> None:
    by_category: dict[str, tuple[int, int]] = {}
    for category in COLORS:
        selected = [run for run in runs if run.category == category]
        by_category[category] = (len(selected), sum(run.bp for run in selected))
    all_span_audit = span_audit_metrics(all_runs)
    drawn_span_audit = span_audit_metrics(runs)
    with path.open("w", newline="") as handle:
        fields = ["metric", "value"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerow({"metric": "source", "value": display_path(source)})
        writer.writerow({"metric": "end_to_end_merge_gap_bp", "value": str(CONTIGUOUS_MERGE_GAP_BP)})
        writer.writerow({"metric": "inter_beats_same_segments_2kb", "value": str(len(segments))})
        writer.writerow({"metric": "inter_beats_same_end_to_end_runs", "value": str(len(all_runs))})
        writer.writerow({"metric": "drawn_high_conf_runs", "value": str(len(runs))})
        writer.writerow({"metric": "drawn_high_conf_bp", "value": str(sum(run.bp for run in runs))})
        for metric, value in all_span_audit.items():
            writer.writerow({"metric": f"all_{metric}", "value": str(value)})
        for metric, value in drawn_span_audit.items():
            writer.writerow({"metric": f"drawn_{metric}", "value": str(value)})
        for category, (count, bp) in by_category.items():
            writer.writerow({"metric": f"{category}_runs", "value": str(count)})
            writer.writerow({"metric": f"{category}_bp", "value": str(bp)})


def write_homolog_summary(
    path: Path,
    homolog_runs: list[Run],
    all_homolog_runs: list[Run],
    homolog_segments: list[Segment],
    inter_runs: list[Run],
    source: Path,
) -> None:
    all_span_audit = span_audit_metrics(all_homolog_runs)
    drawn_span_audit = span_audit_metrics(homolog_runs)
    with path.open("w", newline="") as handle:
        fields = ["metric", "value"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerow({"metric": "source", "value": display_path(source)})
        writer.writerow({"metric": "homolog_layer_definition", "value": "grouped same_chrom father-child donor intervals drawn to child intervals from the 10:10 IMPG class-winner table"})
        writer.writerow({"metric": "homolog_min_identity", "value": str(HOMOLOG_MIN_IDENTITY)})
        writer.writerow({"metric": "homolog_min_bp", "value": str(HOMOLOG_MIN_BP)})
        writer.writerow({"metric": "end_to_end_merge_gap_bp", "value": str(CONTIGUOUS_MERGE_GAP_BP)})
        writer.writerow({"metric": "homolog_segments_2kb", "value": str(len(homolog_segments))})
        writer.writerow({"metric": "homolog_end_to_end_runs", "value": str(len(all_homolog_runs))})
        writer.writerow({"metric": "drawn_homolog_runs", "value": str(len(homolog_runs))})
        writer.writerow({"metric": "drawn_homolog_bp", "value": str(sum(run.bp for run in homolog_runs))})
        for metric, value in all_span_audit.items():
            writer.writerow({"metric": f"all_{metric}", "value": str(value)})
        for metric, value in drawn_span_audit.items():
            writer.writerow({"metric": f"drawn_{metric}", "value": str(value)})
        writer.writerow({"metric": "drawn_interchrom_runs", "value": str(len(inter_runs))})
        writer.writerow({"metric": "drawn_interchrom_bp", "value": str(sum(run.bp for run in inter_runs))})


def write_merge_audit(
    path: Path,
    inter_segments: list[Segment],
    inter_runs: list[Run],
    inter_drawn_runs: list[Run],
    homolog_segments: list[Segment],
    homolog_runs: list[Run],
    homolog_drawn_runs: list[Run],
    source: Path,
) -> None:
    fields = ["layer", "metric", "value"]
    rows: list[dict[str, str]] = []
    for layer, segments, runs, drawn_runs in [
        ("interchrom", inter_segments, inter_runs, inter_drawn_runs),
        ("homolog", homolog_segments, homolog_runs, homolog_drawn_runs),
    ]:
        audit = merge_audit_metrics(segments, runs)
        span_audit = span_audit_metrics(runs)
        drawn_span_audit = span_audit_metrics(drawn_runs)
        values = {
            "source": display_path(source),
            "merge_rule": "same query_seq/donor_seq with query-adjacent and donor-adjacent endpoints in a consistent donor direction",
            "end_to_end_merge_gap_bp": str(CONTIGUOUS_MERGE_GAP_BP),
            "raw_segments_2kb": str(audit["raw_segments_2kb"]),
            "end_to_end_runs": str(audit["end_to_end_runs"]),
            "drawn_runs": str(len(drawn_runs)),
            "drawn_bp": str(sum(run.bp for run in drawn_runs)),
            "absorbed_fragments": str(audit["absorbed_fragments"]),
            "multi_window_runs": str(audit["multi_window_runs"]),
            "max_absorbed_query_endpoint_gap_bp": str(audit["max_absorbed_query_endpoint_gap_bp"]),
            "max_absorbed_donor_endpoint_gap_bp": str(audit["max_absorbed_donor_endpoint_gap_bp"]),
            "merged_windows": str(audit["merged_windows"]),
        }
        for metric, value in span_audit.items():
            values[f"all_{metric}"] = str(value)
        for metric, value in drawn_span_audit.items():
            values[f"drawn_{metric}"] = str(value)
        for metric, value in values.items():
            rows.append({"layer": layer, "metric": metric, "value": value})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def ribbon_path(xa0: float, xa1: float, ya: float, xb0: float, xb1: float, yb: float) -> str:
    c = (yb - ya) * 0.48
    return (
        f"M {xa0:.2f} {ya:.2f} "
        f"C {xa0:.2f} {ya + c:.2f}, {xb0:.2f} {yb - c:.2f}, {xb0:.2f} {yb:.2f} "
        f"L {xb1:.2f} {yb:.2f} "
        f"C {xb1:.2f} {yb - c:.2f}, {xa1:.2f} {ya + c:.2f}, {xa1:.2f} {ya:.2f} Z"
    )


def draw_genome_track(svg: SVG, layout: GenomeLayout, y: float) -> None:
    svg.text(TRACK_LABEL_X, y + TRACK_H - 2, layout.label, 30, "700", TEXT, "end")
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


def draw_chromosome_boundary_bars(svg: SVG, layout: GenomeLayout, y: float) -> None:
    xs: list[float] = []
    last_chrom: str | None = None
    for chrom in CHROM_ORDER:
        if chrom not in layout.lengths:
            continue
        xs.append(x_for(layout, chrom, 0))
        last_chrom = chrom
    if last_chrom is not None:
        xs.append(x_for(layout, last_chrom, layout.lengths[last_chrom]))

    for x in xs:
        svg.line(x, y - 7, x, y + TRACK_H + 7, CHROM_BORDER, 1.8, 1.0)


def draw_chromosome_labels(svg: SVG, layout: GenomeLayout, y: float) -> None:
    for chrom in CHROM_ORDER:
        if chrom not in layout.lengths:
            continue
        x0 = x_for(layout, chrom, 0)
        x1 = x_for(layout, chrom, layout.lengths[chrom])
        if x1 - x0 >= 38:
            label = chrom.replace("chr", "")
            svg.text((x0 + x1) / 2, y - 9, label, 20, "700", TEXT, "middle")


def draw_chromosome_track_overlay(svg: SVG, layout: GenomeLayout, y: float) -> None:
    draw_chromosome_boundary_bars(svg, layout, y)
    draw_chromosome_labels(svg, layout, y)


def draw_all_chromosome_track_overlays(svg: SVG, query_layout: GenomeLayout, hap1_layout: GenomeLayout, hap2_layout: GenomeLayout) -> None:
    draw_chromosome_track_overlay(svg, hap1_layout, Y_H1)
    draw_chromosome_track_overlay(svg, query_layout, Y_QUERY)
    draw_chromosome_track_overlay(svg, hap2_layout, Y_H2)


def ribbon_opacity(category: str) -> float:
    if category in {"PAR_XY", "chr5_chr1_candidate", "chr9_chr3_candidate"}:
        return 0.75
    if category == "acro_acro":
        return 0.34
    return 0.42


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
        y = Y_QUERY - 140 - lane * 42
        color = COLORS[run.category]
        svg.line(center, Y_QUERY - 9, center, y + 8, color, 1.4, 0.9)
        svg.text(center, y, label_text(run), 27, "700", color, "middle")
        svg.text(center, y + int(29 * FONT_SCALE), f"{run.bp / 1000:.0f} kb, {run.donor_haplotype}", 20, "400", MUTED, "middle")


def draw_legend(svg: SVG, runs: list[Run]) -> None:
    x = TRACK_X0
    y = LEGEND_Y
    svg.text(x, y, "Ribbon classes", 24, "700", MUTED)
    x += int(210 * FONT_SCALE)
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
        x += int((390 if category != "other_nonacro" else 220) * FONT_SCALE)


def draw_homolog_legend(svg: SVG, runs: list[Run], homolog_runs: list[Run]) -> None:
    sw = 27 * FONT_SCALE   # swatch scales with the font
    sh = 21 * FONT_SCALE
    # row 1: title + the same-chromosome homologous class
    x = TRACK_X0
    y = LEGEND_Y
    svg.text(x, y, "Ribbon classes", 24, "700", MUTED)
    x += int(230 * FONT_SCALE)
    svg.rect(x, y - sh + 3, sw, sh, HOMOLOG_COLOR, "none", 0, 0.9)
    svg.text(x + sw + 16, y, f"same chromosome ({len(homolog_runs)})", 21, "400", TEXT)
    x += int(560 * FONT_SCALE)
    acro_n = sum(1 for run in runs if run.category == "acro_acro")
    svg.rect(x, y - sh + 3, sw, sh, COLORS["acro_acro"], "none", 0, 0.9)
    svg.text(x + sw + 16, y, f"acrocentric ({acro_n})", 21, "400", TEXT)
    # row 2: highlighted candidate classes
    x = TRACK_X0
    y = LEGEND_Y + int(54 * FONT_SCALE)
    for category, label in [
        ("PAR_XY", "PAR1"),
        ("chr5_chr1_candidate", "chr5q/chr1p"),
        ("chr9_chr3_candidate", "chr9q/chr3q"),
    ]:
        count = sum(1 for run in runs if run.category == category)
        if count == 0:
            continue
        svg.rect(x, y - sh + 3, sw, sh, COLORS[category], "none", 0, 0.9)
        svg.text(x + sw + 16, y, f"{label} ({count})", 21, "400", TEXT)
        x += int(360 * FONT_SCALE)


def render(
    runs: list[Run],
    query_layout: GenomeLayout,
    hap1_layout: GenomeLayout,
    hap2_layout: GenomeLayout,
    config: RenderConfig,
) -> None:
    svg = SVG(PAGE_W, PAGE_H)

    for y, label in [
        (Y_H1 - 16, config.target_h1_label),
        (Y_QUERY + TRACK_H + 8, config.query_label),
        (Y_H2 - 16, config.target_h2_label),
    ]:
        svg.line(TRACK_X0, y, TRACK_X0 + TRACK_W, y, "#f1f3f4", 0.8, 0.8)

    draw_genome_track(svg, query_layout, Y_QUERY)
    draw_genome_track(svg, hap1_layout, Y_H1)
    draw_genome_track(svg, hap2_layout, Y_H2)

    target_layouts = {"h1": hap1_layout, "h2": hap2_layout}
    target_y = {"h1": Y_H1, "h2": Y_H2}

    for run in runs:
        target_layout_obj = target_layouts.get(run.donor_haplotype)
        if target_layout_obj is None or run.target_chrom not in target_layout_obj.lengths:
            continue
        qx0, qx1 = interval_x_with_min_width(query_layout, run.query_chrom, run.query_start, run.query_end, INTERCHROM_RIBBON_MIN_W)
        dx0, dx1 = interval_x_with_min_width(target_layout_obj, run.target_chrom, run.donor_start, run.donor_end, INTERCHROM_RIBBON_MIN_W)
        color = COLORS[run.category]
        donor_edge_y, child_edge_y = ribbon_y_for(target_y[run.donor_haplotype], Y_QUERY)
        d = ribbon_path(dx0, dx1, donor_edge_y, qx0, qx1, child_edge_y)
        svg.path(d, color, "none", 0, ribbon_opacity(run.category))

    for run in runs:
        target_layout_obj = target_layouts.get(run.donor_haplotype)
        if target_layout_obj is None or run.target_chrom not in target_layout_obj.lengths:
            continue
        color = COLORS[run.category]
        qx0, qx1 = interval_x_with_min_width(query_layout, run.query_chrom, run.query_start, run.query_end, INTERCHROM_RIBBON_MIN_W)
        dx0, dx1 = interval_x_with_min_width(target_layout_obj, run.target_chrom, run.donor_start, run.donor_end, INTERCHROM_RIBBON_MIN_W)
        emph = 1.0 if run.category in {"PAR_XY", "chr5_chr1_candidate", "chr9_chr3_candidate"} else 0.55
        draw_interval(svg, qx0, qx1, Y_QUERY, color, 0.88 * emph)
        draw_interval(svg, dx0, dx1, target_y[run.donor_haplotype], color, 0.88 * emph)

    draw_callouts(svg, runs, query_layout)
    draw_all_chromosome_track_overlays(svg, query_layout, hap1_layout, hap2_layout)
    draw_legend(svg, runs)

    total_bp = sum(run.bp for run in runs)
    svg.write(config.svg_out)


def render_homolog_context(
    inter_runs: list[Run],
    homolog_runs: list[Run],
    query_layout: GenomeLayout,
    hap1_layout: GenomeLayout,
    hap2_layout: GenomeLayout,
    config: RenderConfig,
) -> None:
    svg = SVG(PAGE_W, PAGE_H)

    for y, label in [
        (Y_H1 - 16, config.target_h1_label),
        (Y_QUERY + TRACK_H + 8, config.query_label),
        (Y_H2 - 16, config.target_h2_label),
    ]:
        svg.line(TRACK_X0, y, TRACK_X0 + TRACK_W, y, "#f1f3f4", 0.8, 0.8)

    draw_genome_track(svg, query_layout, Y_QUERY)
    draw_genome_track(svg, hap1_layout, Y_H1)
    draw_genome_track(svg, hap2_layout, Y_H2)

    target_layouts = {"h1": hap1_layout, "h2": hap2_layout}
    target_y = {"h1": Y_H1, "h2": Y_H2}

    for run in homolog_runs:
        target_layout_obj = target_layouts.get(run.donor_haplotype)
        if target_layout_obj is None or run.target_chrom not in target_layout_obj.lengths:
            continue
        min_w = homolog_visual_width(run.bp)
        qx0, qx1 = interval_x_with_min_width(query_layout, run.query_chrom, run.query_start, run.query_end, min_w)
        dx0, dx1 = interval_x_with_min_width(target_layout_obj, run.target_chrom, run.donor_start, run.donor_end, min_w)
        donor_edge_y, child_edge_y = ribbon_y_for(target_y[run.donor_haplotype], Y_QUERY)
        d = ribbon_path(dx0, dx1, donor_edge_y, qx0, qx1, child_edge_y)
        svg.path(d, HOMOLOG_RIBBON, "none", 0, HOMOLOG_RIBBON_OPACITY)

    for run in homolog_runs:
        target_layout_obj = target_layouts.get(run.donor_haplotype)
        if target_layout_obj is None or run.target_chrom not in target_layout_obj.lengths:
            continue
        min_w = homolog_visual_width(run.bp)
        qx0, qx1 = interval_x_with_min_width(query_layout, run.query_chrom, run.query_start, run.query_end, min_w)
        dx0, dx1 = interval_x_with_min_width(target_layout_obj, run.target_chrom, run.donor_start, run.donor_end, min_w)
        draw_interval(svg, qx0, qx1, Y_QUERY, HOMOLOG_COLOR, 0.22)
        draw_interval(svg, dx0, dx1, target_y[run.donor_haplotype], HOMOLOG_COLOR, 0.22)

    for run in inter_runs:
        target_layout_obj = target_layouts.get(run.donor_haplotype)
        if target_layout_obj is None or run.target_chrom not in target_layout_obj.lengths:
            continue
        qx0, qx1 = interval_x_with_min_width(query_layout, run.query_chrom, run.query_start, run.query_end, INTERCHROM_RIBBON_MIN_W)
        dx0, dx1 = interval_x_with_min_width(target_layout_obj, run.target_chrom, run.donor_start, run.donor_end, INTERCHROM_RIBBON_MIN_W)
        color = COLORS[run.category]
        donor_edge_y, child_edge_y = ribbon_y_for(target_y[run.donor_haplotype], Y_QUERY)
        d = ribbon_path(dx0, dx1, donor_edge_y, qx0, qx1, child_edge_y)
        opacity = 0.82 if run.category in {"PAR_XY", "chr5_chr1_candidate", "chr9_chr3_candidate"} else 0.48
        svg.path(d, color, "none", 0, opacity)

    for run in inter_runs:
        target_layout_obj = target_layouts.get(run.donor_haplotype)
        if target_layout_obj is None or run.target_chrom not in target_layout_obj.lengths:
            continue
        color = COLORS[run.category]
        qx0, qx1 = interval_x_with_min_width(query_layout, run.query_chrom, run.query_start, run.query_end, INTERCHROM_RIBBON_MIN_W)
        dx0, dx1 = interval_x_with_min_width(target_layout_obj, run.target_chrom, run.donor_start, run.donor_end, INTERCHROM_RIBBON_MIN_W)
        emph = 1.0 if run.category in {"PAR_XY", "chr5_chr1_candidate", "chr9_chr3_candidate"} else 0.65
        draw_interval(svg, qx0, qx1, Y_QUERY, color, 0.95 * emph)
        draw_interval(svg, dx0, dx1, target_y[run.donor_haplotype], color, 0.95 * emph)

    draw_callouts(svg, inter_runs, query_layout)
    draw_all_chromosome_track_overlays(svg, query_layout, hap1_layout, hap2_layout)
    draw_homolog_legend(svg, inter_runs, homolog_runs)

    svg.write(config.homolog_svg_out)


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
    subprocess.run([rsvg, "-f", "png", "-z", "2", "-o", str(svg_path.with_suffix(".png")), str(svg_path)], check=True)
    return subprocess.run([rsvg, "--version"], check=True, text=True, capture_output=True).stdout.strip()


def convert_outputs(config: RenderConfig) -> list[str]:
    messages: list[str] = []
    version = convert_one(config.svg_out)
    if version is None:
        return ["SVG only: no rsvg-convert found."]
    homolog_version = convert_one(config.homolog_svg_out)
    messages.append(f"converted PDF and PNG with {version}")
    if homolog_version is not None:
        messages.append(f"converted homologous-context PDF and PNG with {homolog_version}")
    return messages


def build_config(args: argparse.Namespace) -> RenderConfig:
    out_dir = args.output_dir
    prefix = args.output_prefix or args.comparison_id
    legacy_default_names = (
        args.output_prefix is None
        and args.comparison_id == "PAN027pat_vs_PAN011_joint"
        and out_dir == DEFAULT_OUT_DIR
    )
    if legacy_default_names:
        svg_out = out_dir / "fig5_homolog_vs_interchrom_whole_genome_ribbon_draft.svg"
        homolog_svg_out = out_dir / "fig5_homologous_recombination_context_ribbon_draft.svg"
        runs_out = out_dir / "whole_genome_ribbon_runs.tsv"
        summary_out = out_dir / "whole_genome_ribbon_summary.tsv"
        homolog_runs_out = out_dir / "whole_genome_homologous_context_runs.tsv"
        homolog_summary_out = out_dir / "whole_genome_homologous_context_summary.tsv"
        merge_audit_out = out_dir / "whole_genome_ribbon_merge_audit.tsv"
        conversion_status = out_dir / "conversion_status.txt"
    else:
        svg_out = out_dir / f"{prefix}.whole_genome_ribbon.svg"
        homolog_svg_out = out_dir / f"{prefix}.whole_genome_homologous_context_ribbon.svg"
        runs_out = out_dir / f"{prefix}.whole_genome_ribbon_runs.tsv"
        summary_out = out_dir / f"{prefix}.whole_genome_ribbon_summary.tsv"
        homolog_runs_out = out_dir / f"{prefix}.whole_genome_homologous_context_runs.tsv"
        homolog_summary_out = out_dir / f"{prefix}.whole_genome_homologous_context_summary.tsv"
        merge_audit_out = out_dir / f"{prefix}.whole_genome_ribbon_merge_audit.tsv"
        conversion_status = out_dir / f"{prefix}.conversion_status.txt"
    return RenderConfig(
        comparison_id=args.comparison_id,
        class_winners=args.class_winners,
        query_fai=args.query_fai,
        target_fai=args.target_fai,
        out_dir=out_dir,
        query_label=args.query_label,
        target_h1_label=args.target_h1_label,
        target_h2_label=args.target_h2_label,
        layer_label=args.layer_label,
        svg_out=svg_out,
        homolog_svg_out=homolog_svg_out,
        runs_out=runs_out,
        summary_out=summary_out,
        homolog_runs_out=homolog_runs_out,
        homolog_summary_out=homolog_summary_out,
        merge_audit_out=merge_audit_out,
        conversion_status=conversion_status,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render whole-genome homologous-competition ribbons from a 10:10 "
            "IMPG class-winner TSV.GZ. Filtering thresholds are fixed to the "
            "accepted Fig. 5 draft values."
        )
    )
    parser.add_argument("--comparison-id", default="PAN027pat_vs_PAN011_joint")
    parser.add_argument("--class-winners", type=Path, default=CLASS_WINNERS)
    parser.add_argument("--query-fai", type=Path, default=DEFAULT_QUERY_FAI)
    parser.add_argument("--target-fai", type=Path, default=DEFAULT_TARGET_FAI)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--output-prefix", default=None)
    parser.add_argument("--query-label", default="PAN027 paternal child query")
    parser.add_argument("--target-h1-label", default="PAN011 father donor h1")
    parser.add_argument("--target-h2-label", default="PAN011 father donor h2")
    parser.add_argument("--layer-label", default="father-child")
    return parser.parse_args()


def main() -> None:
    config = build_config(parse_args())
    config.out_dir.mkdir(parents=True, exist_ok=True)
    _query_seqs, query_lengths = read_query_fai(config.query_fai)
    _target_seqs, target_lengths_by_key, _target_lengths_by_seq = read_target_fai(config.target_fai)
    query_layout = layout_for_lengths(config.query_label, query_lengths)
    hap1_layout = target_layout("h1", target_lengths_by_key, config.target_h1_label)
    hap2_layout = target_layout("h2", target_lengths_by_key, config.target_h2_label)

    segments = read_segments(config.class_winners)
    all_runs = group_runs(segments)
    validate_run_spans(all_runs, "interchrom all")
    runs = high_confidence_runs(all_runs)
    validate_run_spans(runs, "interchrom drawn")
    homolog_segments = read_homolog_segments(config.class_winners)
    all_homolog_runs = group_homolog_runs(homolog_segments)
    validate_run_spans(all_homolog_runs, "homolog all")
    homolog_runs = high_confidence_homolog_runs(all_homolog_runs)
    validate_run_spans(homolog_runs, "homolog drawn")
    write_runs(config.runs_out, runs)
    write_summary(config.summary_out, runs, all_runs, segments, config.class_winners)
    write_homolog_runs(config.homolog_runs_out, homolog_runs)
    write_homolog_summary(config.homolog_summary_out, homolog_runs, all_homolog_runs, homolog_segments, runs, config.class_winners)
    write_merge_audit(config.merge_audit_out, segments, all_runs, runs, homolog_segments, all_homolog_runs, homolog_runs, config.class_winners)
    render(runs, query_layout, hap1_layout, hap2_layout, config)
    render_homolog_context(runs, homolog_runs, query_layout, hap1_layout, hap2_layout, config)
    messages = convert_outputs(config)
    config.conversion_status.write_text("\n".join(messages) + "\n")
    for message in messages:
        print(message)
    print(
        f"segments={len(segments)} all_runs={len(all_runs)} drawn_runs={len(runs)} "
        f"homolog_segments={len(homolog_segments)} homolog_runs={len(all_homolog_runs)} "
        f"drawn_homolog_runs={len(homolog_runs)}"
    )
    print(f"wrote {config.svg_out}")
    print(f"wrote {config.runs_out}")


if __name__ == "__main__":
    main()
