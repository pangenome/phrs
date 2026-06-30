#!/usr/bin/env python3
"""Draft donor/recipient ribbon view for PAN027 paternal Fig5 candidates."""

from __future__ import annotations

import csv
import html
import math
import os
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
HERE = ROOT / "paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft"
ZOOM_DIR = ROOT / "paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels"
SEGMENTS_TSV = ZOOM_DIR / "zoom_window_segments.tsv"
PHR_TSV = ZOOM_DIR / "zoom_phr_intervals.tsv"
PHR_TABLE = Path(
    "/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/"
    "all-vs-all.1Mb.p95.id95.len.tsv"
)
TARGET_FAI = (
    ROOT
    / ".wg-worktrees/agent-2636/paper_prep/_brainstorming/"
    / "pedigree_whole_genome_wfmash_p95_updated_bin/inputs/"
    / "PAN027pat_vs_PAN011_joint.target.fa.fai"
)

SVG_OUT = HERE / "fig5_donor_recipient_ribbon_draft.svg"
RUNS_OUT = HERE / "donor_recipient_runs.tsv"
CONVERSION_STATUS = HERE / "conversion_status.txt"

PANEL_ORDER = ["chrX_p", "chr5_q", "chr9_q"]
PANEL_TITLES = {
    "chrX_p": "chrXp child recipient <- chrYp father donor",
    "chr5_q": "chr5q child recipient <- chr1p father donor",
    "chr9_q": "chr9q child recipient <- chr3q father donor",
}
DOMINANT_BY_PANEL = {
    "chrX_p": {"chrY"},
    "chr5_q": {"chr1"},
    "chr9_q": {"chr3"},
}
COLORS = {
    "chrY": "#E7298A",
    "chr1": "#4E79A7",
    "chr3": "#D95F02",
    "other": "#9E9E9E",
}
PHR_FILL = "#111111"
CONTINUATION_FILL = "#9aa0a6"

PAGE_W = 2600
TRACK_X0 = 540
TRACK_W = 1580
TRACK_H = 20
PANEL_GAP = 78
DONOR_ROW_GAP = 54
TOP = 152
TEXT = "#222222"
MUTED = "#5f6368"
LIGHT = "#d7dadd"
FAINT = "#f4f5f6"

LOC_RE = re.compile(r"^(?P<seq>.+):(?P<start>[0-9]+)-(?P<end>[0-9]+)$")
TARGET_RE = re.compile(r"PAN011#joint#(?P<hap>h[12])_(?P<chrom>chr(?:[0-9]+|X|Y|M))")
PHR_SEQ_RE = re.compile(
    r"^PAN011#(?P<hap_no>[12])#(?P<chrom>chr(?:[0-9]+|X|Y|M))\.[^:]+:"
    r"(?P<seq_start>[0-9]+)-(?P<seq_end>[0-9]+)_(?P<label_chrom>chr(?:[0-9]+|X|Y|M))_(?P<arm>[pq])arm$"
)


@dataclass
class Segment:
    panel_order: int
    panel_id: str
    panel_label: str
    query_name: str
    query_chrom: str
    arm: str
    query_length: int
    zoom_bp: int
    query_start: int
    query_end: int
    recipient_start: int
    recipient_end: int
    target_chrom: str
    target_bucket: str
    target_haplotype: str
    donor_seq: str
    donor_start: int
    donor_end: int
    bp: int
    same_identity: float
    inter_identity: float


@dataclass
class Run:
    run_id: str
    panel_order: int
    panel_id: str
    panel_label: str
    query_chrom: str
    arm: str
    query_length: int
    zoom_bp: int
    target_chrom: str
    target_haplotype: str
    donor_seq: str
    query_start: int
    query_end: int
    recipient_start: int
    recipient_end: int
    donor_start: int
    donor_end: int
    bp: int
    windows: int
    drawn: bool = False
    donor_window_start: int = 0
    donor_window_end: int = 0
    donor_window_label: str = ""
    suppressed: bool = False

    @property
    def target_bucket(self) -> str:
        return self.target_chrom if self.target_chrom in COLORS else "other"

    @property
    def label(self) -> str:
        prefix = self.target_chrom
        if self.target_bucket == "other":
            prefix = f"other: {self.target_chrom}"
        return f"{prefix} {self.target_haplotype}"


@dataclass(frozen=True)
class PhrInterval:
    seq: str
    haplotype: str
    chrom: str
    arm: str
    full_start: int
    full_end: int


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


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
    ) -> None:
        self.add(
            f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size:.1f}" '
            f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">{esc(text)}</text>'
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

    def polygon(
        self,
        points: list[tuple[float, float]],
        fill: str,
        stroke: str = "none",
        sw: float = 0.0,
        opacity: float = 1.0,
    ) -> None:
        pts = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
        self.add(
            f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" '
            f'stroke-width="{sw:.2f}" opacity="{opacity:.3f}"/>'
        )

    def line(
        self,
        x1: float,
        y1: float,
        x2: float,
        y2: float,
        stroke: str = TEXT,
        sw: float = 1.0,
        opacity: float = 1.0,
    ) -> None:
        self.add(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="{stroke}" stroke-width="{sw:.2f}" opacity="{opacity:.3f}"/>'
        )

    def path(self, d: str, fill: str, stroke: str = "none", sw: float = 0.0, opacity: float = 1.0) -> None:
        self.add(
            f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw:.2f}" opacity="{opacity:.3f}"/>'
        )

    def write(self, path: Path) -> None:
        body = "\n".join(self.items)
        path.write_text(
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" height="{self.height}" '
            f'viewBox="0 0 {self.width} {self.height}">\n'
            '<rect width="100%" height="100%" fill="white"/>\n'
            '<style>text{font-family:Arial,Helvetica,sans-serif;dominant-baseline:alphabetic}</style>\n'
            f"{body}\n</svg>\n"
        )


def read_fai(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            seq, length, *_ = line.rstrip("\n").split("\t")
            lengths[seq] = int(length)
    return lengths


def parse_loc(value: str) -> tuple[str, int, int]:
    match = LOC_RE.match(value)
    if match is None:
        raise ValueError(f"could not parse interval: {value}")
    return match.group("seq"), int(match.group("start")), int(match.group("end"))


def read_segments(path: Path) -> list[Segment]:
    segments: list[Segment] = []
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            donor_seq, donor_start, donor_end = parse_loc(row["inter_group_a"])
            segments.append(
                Segment(
                    panel_order=int(row["panel_order"]),
                    panel_id=row["panel_id"],
                    panel_label=row["panel_label"],
                    query_name=row["query_name"],
                    query_chrom=row["query_chrom"],
                    arm=row["arm"],
                    query_length=int(row["query_length"]),
                    zoom_bp=int(row["zoom_bp"]),
                    query_start=int(row["query_start"]),
                    query_end=int(row["query_end"]),
                    recipient_start=int(row["relative_start"]),
                    recipient_end=int(row["relative_end"]),
                    target_chrom=row["target_chrom"],
                    target_bucket=row["target_bucket"],
                    target_haplotype=row["target_haplotype"],
                    donor_seq=donor_seq,
                    donor_start=donor_start,
                    donor_end=donor_end,
                    bp=int(row["window_overlap_bp"]),
                    same_identity=float(row["same_identity"]),
                    inter_identity=float(row["inter_identity"]),
                )
            )
    return segments


def read_phr_intervals(path: Path) -> dict[str, list[dict[str, str]]]:
    out: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            out[row["panel_id"]].append(row)
    return out


def read_population_phrs(path: Path) -> dict[tuple[str, str, str], list[PhrInterval]]:
    out: dict[tuple[str, str, str], list[PhrInterval]] = defaultdict(list)
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["region_start"] == "." or row["region_end"] == ".":
                continue
            match = PHR_SEQ_RE.match(row["seq"])
            if match is None:
                continue
            haplotype = f"h{match.group('hap_no')}"
            chrom = match.group("chrom")
            arm = "p" if row["arm"] == "parm" else "q"
            flank_start = int(match.group("seq_start"))
            out[(haplotype, chrom, arm)].append(
                PhrInterval(
                    seq=row["seq"],
                    haplotype=haplotype,
                    chrom=chrom,
                    arm=arm,
                    full_start=flank_start + int(row["region_start"]),
                    full_end=flank_start + int(row["region_end"]),
                )
            )
    return out


def group_runs(segments: list[Segment]) -> list[Run]:
    runs: list[Run] = []
    by_panel: dict[str, list[Segment]] = defaultdict(list)
    for segment in segments:
        by_panel[segment.panel_id].append(segment)

    for panel_id in PANEL_ORDER:
        current: Run | None = None
        ordered = sorted(
            by_panel.get(panel_id, []),
            key=lambda s: (s.query_start, s.target_chrom, s.target_haplotype, s.donor_seq, s.donor_start),
        )
        panel_runs: list[Run] = []
        for segment in ordered:
            same = (
                current is not None
                and current.target_chrom == segment.target_chrom
                and current.target_haplotype == segment.target_haplotype
                and current.donor_seq == segment.donor_seq
            )
            query_adjacent = current is not None and segment.query_start <= current.query_end + 1
            donor_adjacent = current is not None and (
                abs(segment.donor_start - current.donor_end) <= 10_000
                or abs(current.donor_start - segment.donor_end) <= 10_000
            )
            if current is not None and same and query_adjacent and donor_adjacent:
                current.query_end = max(current.query_end, segment.query_end)
                current.recipient_start = min(current.recipient_start, segment.recipient_start)
                current.recipient_end = max(current.recipient_end, segment.recipient_end)
                current.donor_start = min(current.donor_start, segment.donor_start)
                current.donor_end = max(current.donor_end, segment.donor_end)
                current.bp += segment.bp
                current.windows += 1
                continue

            run_id = f"{panel_id}_run{len(panel_runs) + 1:02d}"
            current = Run(
                run_id=run_id,
                panel_order=segment.panel_order,
                panel_id=panel_id,
                panel_label=segment.panel_label,
                query_chrom=segment.query_chrom,
                arm=segment.arm,
                query_length=segment.query_length,
                zoom_bp=segment.zoom_bp,
                target_chrom=segment.target_chrom,
                target_haplotype=segment.target_haplotype,
                donor_seq=segment.donor_seq,
                query_start=segment.query_start,
                query_end=segment.query_end,
                recipient_start=segment.recipient_start,
                recipient_end=segment.recipient_end,
                donor_start=segment.donor_start,
                donor_end=segment.donor_end,
                bp=segment.bp,
                windows=1,
            )
            panel_runs.append(current)
        runs.extend(panel_runs)

    return runs


def should_suppress(run: Run) -> bool:
    if run.panel_id == "chr9_q" and run.target_chrom == "chr7":
        return True
    if run.panel_id == "chr5_q" and run.target_chrom == "chr1" and run.target_haplotype == "h1":
        return True
    return False


def donor_window(run: Run, target_lengths: dict[str, int]) -> tuple[int, int, str]:
    length = target_lengths.get(run.donor_seq)
    if length is None:
        pad = 250_000
        center = (run.donor_start + run.donor_end) // 2
        return max(0, center - pad), center + pad, "local"

    width = run.zoom_bp
    if run.donor_end <= width:
        return 0, width, "p tip"
    if run.donor_start >= max(0, length - width):
        return max(0, length - width), length, "q tip"
    center = (run.donor_start + run.donor_end) // 2
    start = max(0, min(length - width, center - width // 2))
    return start, start + width, "internal"


def mark_drawn_runs(runs: list[Run], target_lengths: dict[str, int]) -> None:
    for run in runs:
        run.suppressed = should_suppress(run)
        dominant = run.target_chrom in DOMINANT_BY_PANEL.get(run.panel_id, set())
        run.drawn = not run.suppressed and (run.bp >= 4_000 or dominant)
        start, end, label = donor_window(run, target_lengths)
        run.donor_window_start = start
        run.donor_window_end = end
        run.donor_window_label = label


def write_runs(path: Path, runs: list[Run]) -> None:
    fields = [
        "run_id",
        "panel_order",
        "panel_id",
        "panel_label",
        "query_chrom",
        "arm",
        "query_length",
        "zoom_bp",
        "target_chrom",
        "target_haplotype",
        "donor_seq",
        "query_start",
        "query_end",
        "recipient_start",
        "recipient_end",
        "donor_start",
        "donor_end",
        "bp",
        "windows",
        "drawn",
        "donor_window_start",
        "donor_window_end",
        "donor_window_label",
        "suppressed",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for run in runs:
            writer.writerow({field: getattr(run, field) for field in fields})


def x_for_window(start: int, end: int, window_start: int, window_end: int) -> tuple[float, float]:
    denom = max(1, window_end - window_start)
    x0 = TRACK_X0 + (start - window_start) / denom * TRACK_W
    x1 = TRACK_X0 + (end - window_start) / denom * TRACK_W
    if x1 < x0:
        x0, x1 = x1, x0
    min_w = 3.0
    if x1 - x0 < min_w:
        mid = (x0 + x1) / 2
        x0, x1 = mid - min_w / 2, mid + min_w / 2
    return x0, x1


def draw_continuation(svg: SVG, y: float, side: str) -> None:
    mid = y + TRACK_H / 2
    if side == "right":
        x = TRACK_X0 + TRACK_W
        points = [(x + 8, y + 2), (x + 8, y + TRACK_H - 2), (x + 25, mid)]
    else:
        x = TRACK_X0
        points = [(x - 8, y + 2), (x - 8, y + TRACK_H - 2), (x - 25, mid)]
    svg.polygon(points, CONTINUATION_FILL, "none", 0, 0.95)


def draw_cut_glyphs(svg: SVG, y: float, anchor: str) -> None:
    if anchor == "left":
        draw_continuation(svg, y, "right")
    elif anchor == "right":
        draw_continuation(svg, y, "left")
    else:
        draw_continuation(svg, y, "left")
        draw_continuation(svg, y, "right")


def draw_phr_bar(svg: SVG, y: float, start: int, end: int, window_start: int, window_end: int) -> None:
    px0, px1 = x_for_window(start, end, window_start, window_end)
    svg.rect(px0, y - 17, px1 - px0, 6, PHR_FILL, "none", 0, 1.0, rx=0)


def ribbon_path(xa0: float, xa1: float, ya: float, xb0: float, xb1: float, yb: float) -> str:
    c = abs(yb - ya) * 0.45
    return (
        f"M {xa0:.2f} {ya:.2f} "
        f"C {xa0:.2f} {ya + c:.2f}, {xb0:.2f} {yb - c:.2f}, {xb0:.2f} {yb:.2f} "
        f"L {xb1:.2f} {yb:.2f} "
        f"C {xb1:.2f} {yb - c:.2f}, {xa1:.2f} {ya + c:.2f}, {xa1:.2f} {ya:.2f} Z"
    )


def draw_ticks(svg: SVG, y: float, start: int, end: int, anchor: str) -> None:
    width = end - start
    for offset in range(0, width + 1, 100_000):
        x = TRACK_X0 + offset / width * TRACK_W
        svg.line(x, y + TRACK_H, x, y + TRACK_H + 10, "#777777", 0.8)
        if offset in (0, width) or offset % 200_000 == 0:
            svg.text(x, y + TRACK_H + 26, f"{(start + offset) / 1e6:.3f}", 12, "400", MUTED, "middle")
    svg.text(TRACK_X0 + TRACK_W + 15, y + TRACK_H + 26, "Mb", 12, "400", MUTED)
    svg.text(TRACK_X0 if anchor == "left" else TRACK_X0 + TRACK_W, y + TRACK_H + 44, anchor, 11, "400", MUTED, "middle")


def draw_track(svg: SVG, y: float, label: str, start: int, end: int, anchor: str) -> None:
    svg.text(TRACK_X0 - 28, y + TRACK_H - 1, label, 18, "700", TEXT, "end")
    svg.rect(TRACK_X0, y, TRACK_W, TRACK_H, FAINT, "#c8ccd0", 1.0, rx=2)
    if anchor in {"left", "right"}:
        cap_x = TRACK_X0 if anchor == "left" else TRACK_X0 + TRACK_W
        svg.line(cap_x, y - 3, cap_x, y + TRACK_H + 3, "#111111", 2.2)
    draw_cut_glyphs(svg, y, anchor)
    draw_ticks(
        svg,
        y,
        start,
        end,
        "chr start" if anchor == "left" else "chr end" if anchor == "right" else "local window",
    )


def label_lanes(items: list[tuple[float, float, str]]) -> list[tuple[float, float, str, int]]:
    lanes: list[float] = []
    out: list[tuple[float, float, str, int]] = []
    for x0, x1, label in sorted(items, key=lambda item: (item[0] + item[1]) / 2):
        center = (x0 + x1) / 2
        width = max(42, len(label) * 6.2)
        placed = False
        for lane, last_right in enumerate(lanes):
            if center - width / 2 > last_right + 8:
                lanes[lane] = center + width / 2
                out.append((center, width, label, lane))
                placed = True
                break
        if not placed:
            lane = len(lanes)
            lanes.append(center + width / 2)
            out.append((center, width, label, lane))
    return out


def target_label(run: Run) -> str:
    if run.target_bucket == "other":
        return f"other: {run.target_chrom} {run.target_haplotype}"
    return f"{run.target_chrom} {run.target_haplotype}"


def arm_from_window_label(label: str) -> str:
    if label.startswith("p"):
        return "p"
    if label.startswith("q"):
        return "q"
    return ""


def donor_label(run: Run) -> str:
    match = TARGET_RE.match(run.donor_seq)
    if match:
        hap_number = match.group("hap")[1:]
        return (
            f"PAN011 haplotype {hap_number} (father) "
            f"{match.group('chrom')}{arm_from_window_label(run.donor_window_label)}"
        )
    return run.donor_seq.replace("PAN011#joint#", "PAN011 ")


def draw_panel(
    svg: SVG,
    panel_id: str,
    runs: list[Run],
    phr_by_panel: dict[str, list[dict[str, str]]],
    donor_phrs: dict[tuple[str, str, str], list[PhrInterval]],
    y: float,
) -> float:
    panel_runs = [run for run in runs if run.panel_id == panel_id]
    drawn = [run for run in panel_runs if run.drawn and not run.suppressed]
    if not panel_runs:
        return y

    donor_rows: list[tuple[tuple[str, int, int, str], list[Run]]] = []
    donor_row_lookup: dict[tuple[str, int, int, str], list[Run]] = {}
    for run in sorted(drawn, key=lambda item: (item.recipient_start, item.donor_seq, item.donor_start)):
        key = (run.donor_seq, run.donor_window_start, run.donor_window_end, run.donor_window_label)
        if key not in donor_row_lookup:
            donor_row_lookup[key] = []
            donor_rows.append((key, donor_row_lookup[key]))
        donor_row_lookup[key].append(run)

    first = panel_runs[0]
    rec_start = 0 if first.arm == "p" else first.query_length - first.zoom_bp
    rec_end = first.zoom_bp if first.arm == "p" else first.query_length
    rec_anchor = "left" if first.arm == "p" else "right"

    panel_h = 122 + max(1, len(donor_rows)) * DONOR_ROW_GAP
    svg.rect(44, y - 48, PAGE_W - 88, panel_h, "#ffffff", "#e0e2e4", 0.8, rx=3)
    svg.text(68, y - 18, PANEL_TITLES.get(panel_id, panel_id), 25, "700")

    rec_y = y + 44
    draw_track(
        svg,
        rec_y,
        f"PAN027 paternal haplotype (child) {first.query_chrom}{first.arm}",
        rec_start,
        rec_end,
        rec_anchor,
    )

    phr_rows = phr_by_panel.get(panel_id, [])
    for phr in phr_rows:
        draw_phr_bar(svg, rec_y, int(phr["query_full_start"]), int(phr["query_full_end"]), rec_start, rec_end)

    for run in drawn:
        rx0, rx1 = x_for_window(run.query_start, run.query_end, rec_start, rec_end)
        color = COLORS[run.target_bucket]
        svg.rect(rx0, rec_y - 5, rx1 - rx0, TRACK_H + 10, color, "none", 0, 0.95, rx=1)

    donor_y0 = rec_y + 86
    for idx, (_key, row_runs) in enumerate(donor_rows):
        run = row_runs[0]
        donor_y = donor_y0 + idx * DONOR_ROW_GAP
        draw_track(
            svg,
            donor_y,
            donor_label(run),
            run.donor_window_start,
            run.donor_window_end,
            "left" if run.donor_window_label == "p tip" else "right" if run.donor_window_label == "q tip" else "local",
        )

        donor_arm = arm_from_window_label(run.donor_window_label)
        if donor_arm:
            for phr in donor_phrs.get((run.target_haplotype, run.target_chrom, donor_arm), []):
                draw_phr_bar(
                    svg,
                    donor_y,
                    phr.full_start,
                    phr.full_end,
                    run.donor_window_start,
                    run.donor_window_end,
                )

        for item in row_runs:
            color = COLORS[item.target_bucket]
            dx0, dx1 = x_for_window(item.donor_start, item.donor_end, item.donor_window_start, item.donor_window_end)
            svg.rect(dx0, donor_y - 4, dx1 - dx0, TRACK_H + 8, color, "none", 0, 0.95, rx=1)
            rx0, rx1 = x_for_window(item.query_start, item.query_end, rec_start, rec_end)
            d = ribbon_path(rx0, rx1, rec_y + TRACK_H + 7, dx0, dx1, donor_y - 7)
            svg.path(d, color, "none", 0, 0.24 if item.target_bucket == "other" else 0.30)

    return y + panel_h + PANEL_GAP


def draw_legend(svg: SVG, y: float) -> None:
    svg.text(68, y, "Recipient color / ribbon target:", 13, "700", MUTED)
    x = 310
    labels = [("chrY", "target chrY"), ("chr1", "target chr1"), ("chr3", "target chr3"), ("other", "other target, labeled")]
    for bucket, label in labels:
        svg.rect(x, y - 13, 14, 14, COLORS[bucket], "none")
        svg.text(x + 20, y, label, 13, "400", TEXT)
        x += 190 if bucket != "other" else 260


def render(
    runs: list[Run],
    phr_by_panel: dict[str, list[dict[str, str]]],
    donor_phrs: dict[tuple[str, str, str], list[PhrInterval]],
) -> None:
    panel_heights = []
    for panel_id in PANEL_ORDER:
        row_keys = {
            (run.donor_seq, run.donor_window_start, run.donor_window_end, run.donor_window_label)
            for run in runs
            if run.panel_id == panel_id and run.drawn and not run.suppressed
        }
        n = len(row_keys)
        panel_heights.append(122 + max(1, n) * DONOR_ROW_GAP)
    height = TOP + sum(panel_heights) + PANEL_GAP * len(panel_heights) + 40
    svg = SVG(PAGE_W, height)
    svg.text(
        68,
        54,
        "Fig5 donor-recipient ribbon draft: PAN027 paternal haplotype (child) vs PAN011 (father)",
        31,
        "700",
    )
    y = TOP
    for panel_id in PANEL_ORDER:
        y = draw_panel(svg, panel_id, runs, phr_by_panel, donor_phrs, y)
    svg.write(SVG_OUT)


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


def convert_outputs() -> list[str]:
    messages: list[str] = []
    rsvg = find_rsvg()
    if rsvg is None:
        return ["SVG only: no rsvg-convert found."]
    pdf = SVG_OUT.with_suffix(".pdf")
    png = SVG_OUT.with_suffix(".png")
    subprocess.run([rsvg, "-f", "pdf", "-o", str(pdf), str(SVG_OUT)], check=True)
    subprocess.run([rsvg, "-f", "png", "-o", str(png), str(SVG_OUT)], check=True)
    version = subprocess.run([rsvg, "--version"], check=True, text=True, capture_output=True).stdout.strip()
    messages.append(f"converted PDF and PNG with {version} ({rsvg})")
    return messages


def main() -> None:
    HERE.mkdir(parents=True, exist_ok=True)
    segments = read_segments(SEGMENTS_TSV)
    target_lengths = read_fai(TARGET_FAI)
    runs = group_runs(segments)
    mark_drawn_runs(runs, target_lengths)
    write_runs(RUNS_OUT, runs)
    phr_by_panel = read_phr_intervals(PHR_TSV)
    donor_phrs = read_population_phrs(PHR_TABLE)
    render(runs, phr_by_panel, donor_phrs)
    messages = convert_outputs()
    CONVERSION_STATUS.write_text("\n".join(messages) + "\n")
    for message in messages:
        print(message)
    print(f"wrote {SVG_OUT}")
    print(f"wrote {RUNS_OUT}")


if __name__ == "__main__":
    main()
