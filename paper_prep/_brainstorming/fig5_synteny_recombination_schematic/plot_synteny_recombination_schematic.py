#!/usr/bin/env python3
"""Render Fig5 synteny/recombination schematic prototypes.

Regenerate from the repository root with:

    python3 paper_prep/_brainstorming/fig5_synteny_recombination_schematic/plot_synteny_recombination_schematic.py

The renderer intentionally uses only the Python standard library. Geometry is
read from event_manifest.tsv and selected_segments.tsv in this directory.
"""

from __future__ import annotations

import csv
import html
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable


ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
EVENTS_TSV = HERE / "event_manifest.tsv"
SEGMENTS_TSV = HERE / "selected_segments.tsv"
CHROM_SIZES = ROOT / "data" / "chm13.chrom.sizes"

FULL_SVG = HERE / "fig5_synteny_recombination_full.svg"
FOCUS_SVG = HERE / "fig5_synteny_recombination_focus.svg"

PAGE_W = 1680
LEFT_LABEL_W = 360
TRACK_X0 = 430
TRACK_X1 = 1520
TRACK_W = TRACK_X1 - TRACK_X0

TRACK_H = 18
LANE_GAP = 11
TRACK_GAP = 56
EVENT_GAP = 54
TOP = 86

QUERY_KEY = "child"


ROLE_STYLE = {
    "same-chromosome context": {
        "fill": "#bfc5ca",
        "stroke": "#68717a",
        "opacity": 0.34,
        "label": "same-chromosome context",
    },
    "PAR positive control": {
        "fill": "#2076b4",
        "stroke": "#0d4d7c",
        "opacity": 0.58,
        "label": "PAR1 positive-control flow",
    },
    "primary donor": {
        "fill": "#bd3d2a",
        "stroke": "#812516",
        "opacity": 0.58,
        "label": "primary PHR donor flow",
    },
    "side fragment": {
        "fill": "#a87935",
        "stroke": "#76511e",
        "opacity": 0.34,
        "label": "side fragment",
    },
    "low-confidence tail": {
        "fill": "#777777",
        "stroke": "#555555",
        "opacity": 0.24,
        "label": "low-confidence tail",
    },
}

TRACK_FILL = "#f4f5f5"
TRACK_STROKE = "#31363a"
TEXT = "#16191c"
MUTED = "#5f686f"
LIGHT = "#d8dcdf"


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start

    def label(self) -> str:
        return f"{self.chrom}:{self.start:,}-{self.end:,}"


@dataclass
class Segment:
    row: dict[str, str]
    event_id: str
    role: str
    query: Interval
    target: Interval
    query_window: Interval
    target_window: Interval
    target_name: str
    target_arm: str
    target_sample: str
    target_haplotype_label: str
    length: int


@dataclass
class Track:
    key: str
    label: str
    sublabel: str
    kind: str
    lanes: list[str]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_interval(value: str, default_chrom: str | None = None) -> Interval:
    value = value.strip()
    if ":" in value:
        chrom, rest = value.split(":", 1)
    else:
        if default_chrom is None:
            raise ValueError(f"Missing chromosome in interval: {value}")
        chrom, rest = default_chrom, value
    start_s, end_s = rest.split("-", 1)
    return Interval(chrom=chrom, start=int(start_s), end=int(end_s))


def arm_to_chrom(arm: str) -> str:
    if arm.endswith("p") or arm.endswith("q"):
        return arm[:-1]
    return arm


def read_chrom_sizes() -> dict[str, int]:
    sizes: dict[str, int] = {}
    with CHROM_SIZES.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, size = line.split()[:2]
            sizes[chrom] = int(size)
    return sizes


def fmt_bp(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f} Mb"
    if value >= 1_000:
        return f"{value / 1_000:.0f} kb"
    return f"{value} bp"


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


class SVG:
    def __init__(self, width: int, height: int) -> None:
        self.width = width
        self.height = height
        self.parts: list[str] = []

    def add(self, line: str) -> None:
        self.parts.append(line)

    def text(
        self,
        x: float,
        y: float,
        text: str,
        size: int = 13,
        weight: str = "400",
        fill: str = TEXT,
        anchor: str = "start",
        extra: str = "",
    ) -> None:
        self.add(
            f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
            f'font-family="Arial, Helvetica, sans-serif" font-weight="{weight}" '
            f'fill="{fill}" text-anchor="{anchor}" {extra}>{esc(text)}</text>'
        )

    def rect(
        self,
        x: float,
        y: float,
        w: float,
        h: float,
        fill: str,
        stroke: str = "none",
        sw: float = 1,
        rx: float = 0,
        opacity: float | None = None,
        extra: str = "",
    ) -> None:
        op = "" if opacity is None else f' opacity="{opacity:.3f}"'
        self.add(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" '
            f'rx="{rx:.2f}" fill="{fill}" stroke="{stroke}" stroke-width="{sw:.2f}"{op} {extra}/>'
        )

    def line(
        self,
        x1: float,
        y1: float,
        x2: float,
        y2: float,
        stroke: str = TEXT,
        sw: float = 1,
        opacity: float | None = None,
        extra: str = "",
    ) -> None:
        op = "" if opacity is None else f' opacity="{opacity:.3f}"'
        self.add(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="{stroke}" stroke-width="{sw:.2f}"{op} {extra}/>'
        )

    def path(
        self,
        d: str,
        fill: str,
        stroke: str = "none",
        sw: float = 1,
        opacity: float | None = None,
        extra: str = "",
    ) -> None:
        op = "" if opacity is None else f' opacity="{opacity:.3f}"'
        self.add(
            f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw:.2f}"{op} {extra}/>'
        )

    def write(self, path: Path) -> None:
        style = """
        <style>
          .small { font-size: 11px; fill: #5f686f; font-family: Arial, Helvetica, sans-serif; }
          .tiny { font-size: 9px; fill: #5f686f; font-family: Arial, Helvetica, sans-serif; }
        </style>
        """.strip()
        body = "\n".join(self.parts)
        path.write_text(
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" '
            f'height="{self.height}" viewBox="0 0 {self.width} {self.height}">\n'
            f"<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n{style}\n{body}\n</svg>\n"
        )


def load_data() -> tuple[list[dict[str, str]], list[Segment], dict[str, int]]:
    events = sorted(read_tsv(EVENTS_TSV), key=lambda r: int(r["event_order"]))
    chrom_sizes = read_chrom_sizes()
    segments: list[Segment] = []
    for row in read_tsv(SEGMENTS_TSV):
        query_window = parse_interval(row["query_source_window_native_0based_half_open"])
        target_window = parse_interval(row["target_source_window_native_0based_half_open"])
        query = parse_interval(row["native_query_interval_0based_half_open"])
        target = parse_interval(
            row["target_native_interval_0based_half_open_if_recovered"],
            default_chrom=target_window.chrom,
        )
        segments.append(
            Segment(
                row=row,
                event_id=row["event_id"],
                role=row["event_role"],
                query=query,
                target=target,
                query_window=query_window,
                target_window=target_window,
                target_name=row["target_name"],
                target_arm=row["target_arm"],
                target_sample=row["target_sample"],
                target_haplotype_label=row["target_haplotype_label"],
                length=int(row["segment_length_bp"]),
            )
        )
    return events, segments, chrom_sizes


def event_segments(segments: Iterable[Segment], event_id: str) -> list[Segment]:
    return [s for s in segments if s.event_id == event_id]


def track_for_role(role: str) -> str:
    if role == "same-chromosome context":
        return "same"
    if role == "PAR positive control":
        return "par"
    if role == "primary donor":
        return "primary"
    if role in {"side fragment", "low-confidence tail"}:
        return "side"
    return "other"


def lane_label(seg: Segment) -> str:
    return (
        f"{seg.target_sample} {seg.target_haplotype_label} "
        f"{seg.target_arm} {seg.target_window.label()}"
    )


def build_tracks(event: dict[str, str], segs: list[Segment]) -> list[Track]:
    query_seg = segs[0]
    child_label = event["child_query_source"].split(" [", 1)[0]
    tracks = [
        Track(
            key=QUERY_KEY,
            label=f"recombinant/child: {child_label}",
            sublabel=event["query_native_window_0based_half_open"],
            kind="query",
            lanes=["child"],
        )
    ]

    order = ["same", "par", "primary", "side"]
    labels = {
        "same": "parental same-chromosome context",
        "par": "non-homologous PAR1 source",
        "primary": "primary non-homologous PHR donor",
        "side": "secondary side fragments / caveats",
    }
    for key in order:
        source = [s for s in segs if track_for_role(s.role) == key]
        if not source:
            continue
        seen: dict[str, str] = {}
        for seg in source:
            lane = lane_label(seg)
            seen[lane] = lane
        tracks.append(
            Track(
                key=key,
                label=labels[key],
                sublabel=", ".join(sorted({s.target_arm for s in source})),
                kind="source",
                lanes=list(seen.values()),
            )
        )
    return tracks


def lane_y(base_y: float, track: Track, lane: str) -> float:
    index = track.lanes.index(lane)
    total_h = len(track.lanes) * TRACK_H + (len(track.lanes) - 1) * LANE_GAP
    return base_y - total_h / 2 + index * (TRACK_H + LANE_GAP)


def track_height(track: Track) -> float:
    return len(track.lanes) * TRACK_H + (len(track.lanes) - 1) * LANE_GAP


def draw_ideogram(
    svg: SVG,
    x0: float,
    y: float,
    width: float,
    height: float,
    domain_start: int,
    domain_end: int,
    window: Interval | None,
    exact_bands: bool,
) -> None:
    svg.rect(x0, y, width, height, TRACK_FILL, TRACK_STROKE, 1.1, rx=height / 2)
    if exact_bands:
        band_n = 12
    else:
        band_n = 8
    for i in range(band_n):
        if i % 2:
            bx = x0 + width * i / band_n
            bw = width / band_n
            svg.rect(bx, y + 1.2, bw, height - 2.4, "#e3e6e8", "none", rx=0)
    if not exact_bands:
        cx = x0 + width * 0.5
        d = (
            f"M {cx - 9:.2f} {y + 1:.2f} L {cx:.2f} {y + height / 2:.2f} "
            f"L {cx - 9:.2f} {y + height - 1:.2f} Z"
            f"M {cx + 9:.2f} {y + 1:.2f} L {cx:.2f} {y + height / 2:.2f} "
            f"L {cx + 9:.2f} {y + height - 1:.2f} Z"
        )
        svg.path(d, "#cfd4d7", "#9ba3a8", 0.5, 0.7)
    if window is not None:
        denom = max(1, domain_end - domain_start)
        wx0 = x0 + (window.start - domain_start) / denom * width
        wx1 = x0 + (window.end - domain_start) / denom * width
        if wx1 < x0 or wx0 > x0 + width:
            return
        wx0 = max(x0, wx0)
        wx1 = min(x0 + width, wx1)
        min_w = 5.0
        if wx1 - wx0 < min_w:
            mid = (wx0 + wx1) / 2
            wx0 = max(x0, mid - min_w / 2)
            wx1 = min(x0 + width, mid + min_w / 2)
        svg.rect(wx0, y - 3, wx1 - wx0, height + 6, "none", "#111111", 1.4, rx=2)


def ribbon_path(xa0: float, xa1: float, ya: float, xb0: float, xb1: float, yb: float) -> str:
    min_w = 1.4
    if abs(xa1 - xa0) < min_w:
        mid = (xa0 + xa1) / 2
        xa0, xa1 = mid - min_w / 2, mid + min_w / 2
    if abs(xb1 - xb0) < min_w:
        mid = (xb0 + xb1) / 2
        xb0, xb1 = mid - min_w / 2, mid + min_w / 2
    c = abs(yb - ya) * 0.42
    if yb > ya:
        c1y = ya + c
        c2y = yb - c
    else:
        c1y = ya - c
        c2y = yb + c
    return (
        f"M {xa0:.2f} {ya:.2f} "
        f"C {xa0:.2f} {c1y:.2f}, {xb0:.2f} {c2y:.2f}, {xb0:.2f} {yb:.2f} "
        f"L {xb1:.2f} {yb:.2f} "
        f"C {xb1:.2f} {c2y:.2f}, {xa1:.2f} {c1y:.2f}, {xa1:.2f} {ya:.2f} Z"
    )


def full_mapper(domain_start: int, domain_end: int) -> Callable[[Interval], tuple[float, float]]:
    denom = max(1, domain_end - domain_start)

    def mapper(interval: Interval) -> tuple[float, float]:
        return (
            TRACK_X0 + (interval.start - domain_start) / denom * TRACK_W,
            TRACK_X0 + (interval.end - domain_start) / denom * TRACK_W,
        )

    return mapper


def focus_mapper(window: Interval) -> Callable[[Interval], tuple[float, float]]:
    denom = max(1, window.end - window.start)

    def mapper(interval: Interval) -> tuple[float, float]:
        return (
            TRACK_X0 + (interval.start - window.start) / denom * TRACK_W,
            TRACK_X0 + (interval.end - window.start) / denom * TRACK_W,
        )

    return mapper


def role_totals(segs: list[Segment]) -> str:
    totals: dict[str, int] = {}
    for seg in segs:
        totals[seg.role] = totals.get(seg.role, 0) + seg.length
    ordered = [
        ("PAR positive control", "PAR"),
        ("primary donor", "primary donor"),
        ("side fragment", "side fragments"),
        ("low-confidence tail", "low-conf. tail"),
        ("same-chromosome context", "same-chr"),
    ]
    pieces = [f"{label} {fmt_bp(totals[role])}" for role, label in ordered if role in totals]
    return "; ".join(pieces)


def draw_event(
    svg: SVG,
    event: dict[str, str],
    segs: list[Segment],
    chrom_sizes: dict[str, int],
    y0: float,
    mode: str,
) -> float:
    tracks = build_tracks(event, segs)
    row_h = sum(max(track_height(t), TRACK_H) for t in tracks) + TRACK_GAP * (len(tracks) - 1)
    y_center0 = y0 + 94

    title = {
        "PAR1_XY_positive_control": "PAR1 positive control",
        "PAN027_chr9q_chr3q_PHR_candidate": "PHR candidate 1: PAN027 chr9q with chr3q donor",
        "PAN028_chr9q_chr3q_PHR_candidate": "PHR candidate 2: PAN028 strict chr9q path",
    }.get(event["event_id"], event["event_id"])

    svg.text(42, y0 + 26, f"{event['event_order']}. {title}", 20, "700")
    svg.text(42, y0 + 49, event["transmission"], 13, "400", MUTED)
    svg.text(42, y0 + 68, role_totals(segs), 12, "400", MUTED)

    if mode == "full":
        svg.text(
            TRACK_X0,
            y0 + 25,
            "Full chromosome/arm context; bands are neutral schematic bands from chromosome-size context, not exact cytobands.",
            12,
            "400",
            MUTED,
        )
    else:
        svg.text(
            TRACK_X0,
            y0 + 25,
            "Focus view uses native assembly window coordinates and one physical scale across tracks.",
            12,
            "400",
            MUTED,
        )

    track_y: dict[str, float] = {}
    cursor = y_center0
    for track in tracks:
        track_y[track.key] = cursor
        cursor += max(track_height(track), TRACK_H) + TRACK_GAP

    lane_positions: dict[tuple[str, str], float] = {}
    lane_maps: dict[tuple[str, str], Callable[[Interval], tuple[float, float]]] = {}

    # Draw tracks and build mappers.
    for track in tracks:
        cy = track_y[track.key]
        svg.text(42, cy - max(8, track_height(track) / 2 + 7), track.label, 13, "700")
        svg.text(42, cy - max(8, track_height(track) / 2 - 10), track.sublabel, 11, "400", MUTED)
        for lane in track.lanes:
            ly = lane_y(cy, track, lane)
            lane_positions[(track.key, lane)] = ly

            if track.key == QUERY_KEY:
                window = segs[0].query_window
                chrom = window.chrom
            else:
                lane_seg = next(s for s in segs if track_for_role(s.role) == track.key and lane_label(s) == lane)
                window = lane_seg.target_window
                chrom = window.chrom

            if mode == "full":
                max_seen = max(
                    [window.end]
                    + [s.query.end for s in segs if track.key == QUERY_KEY and s.query.chrom == chrom]
                    + [s.target.end for s in segs if track.key != QUERY_KEY and lane_label(s) == lane]
                )
                domain_end = max(chrom_sizes.get(chrom, max_seen), max_seen)
                mapper = full_mapper(0, domain_end)
                draw_ideogram(svg, TRACK_X0, ly, TRACK_W, TRACK_H, 0, domain_end, window, exact_bands=False)
                svg.text(TRACK_X1 + 13, ly + 13, f"{chrom} length context {fmt_bp(domain_end)}", 10, "400", MUTED)
            else:
                mapper = focus_mapper(window)
                draw_ideogram(svg, TRACK_X0, ly, TRACK_W, TRACK_H, window.start, window.end, None, exact_bands=True)
                svg.text(TRACK_X1 + 13, ly + 13, window.label(), 10, "400", MUTED)
                svg.text(TRACK_X0, ly - 6, f"{window.start:,}", 9, "400", MUTED)
                svg.text(TRACK_X1, ly - 6, f"{window.end:,}", 9, "400", MUTED, anchor="end")

            if track.key != QUERY_KEY:
                svg.text(TRACK_X0 - 10, ly + 13, lane, 10, "400", MUTED, anchor="end")
            lane_maps[(track.key, lane)] = mapper

    # Ribbons first, then block overlays.
    for seg in segs:
        style = ROLE_STYLE[seg.role]
        source_key = track_for_role(seg.role)
        source_lane = lane_label(seg)
        source_y = lane_positions[(source_key, source_lane)] + TRACK_H + 1
        query_y = lane_positions[(QUERY_KEY, "child")] - 1
        sx0, sx1 = lane_maps[(source_key, source_lane)](seg.target)
        qx0, qx1 = lane_maps[(QUERY_KEY, "child")](seg.query)
        d = ribbon_path(sx0, sx1, source_y, qx0, qx1, query_y)
        svg.path(d, style["fill"], style["stroke"], 0.45, style["opacity"])

    for seg in segs:
        style = ROLE_STYLE[seg.role]
        qx0, qx1 = lane_maps[(QUERY_KEY, "child")](seg.query)
        qy = lane_positions[(QUERY_KEY, "child")]
        sx0, sx1 = lane_maps[(track_for_role(seg.role), lane_label(seg))](seg.target)
        sy = lane_positions[(track_for_role(seg.role), lane_label(seg))]
        min_block = 1.8 if mode == "focus" else 3.2
        for x0, x1, y, stroke_extra in [
            (qx0, qx1, qy, "#202428"),
            (sx0, sx1, sy, style["stroke"]),
        ]:
            w = max(min_block, x1 - x0)
            if w != x1 - x0:
                mid = (x0 + x1) / 2
                x0 = mid - w / 2
            svg.rect(x0, y - 2.5, w, TRACK_H + 5, style["fill"], stroke_extra, 0.6, rx=1.4, opacity=0.86)

        if mode == "focus" and seg.role != "same-chromosome context":
            mid = (qx0 + qx1) / 2
            if seg.length >= 900:
                svg.line(mid, qy - 5, mid, qy - 19, style["stroke"], 0.8)
                svg.text(
                    mid,
                    qy - 24,
                    f"{seg.role}: {fmt_bp(seg.length)}",
                    9,
                    "400",
                    style["stroke"],
                    anchor="middle",
                )

    # Native-coordinate event-window labels.
    qwin = segs[0].query_window
    svg.text(TRACK_X0, y0 + row_h + 116, f"Child/query native window: {qwin.label()}", 11, "700", TEXT)
    if mode == "focus":
        bar_y = y0 + row_h + 138
        bar_w = TRACK_W * 100_000 / max(1, qwin.length)
        svg.line(TRACK_X0, bar_y, TRACK_X0 + bar_w, bar_y, TEXT, 2.2)
        svg.line(TRACK_X0, bar_y - 5, TRACK_X0, bar_y + 5, TEXT, 1.2)
        svg.line(TRACK_X0 + bar_w, bar_y - 5, TRACK_X0 + bar_w, bar_y + 5, TEXT, 1.2)
        svg.text(TRACK_X0 + bar_w / 2, bar_y + 17, "100 kb", 10, "700", TEXT, anchor="middle")

    return y0 + row_h + EVENT_GAP + 120


def legend(svg: SVG, y: float) -> None:
    x = 42
    svg.text(x, y, "Flow legend", 12, "700")
    x += 90
    for role in [
        "PAR positive control",
        "primary donor",
        "side fragment",
        "low-confidence tail",
        "same-chromosome context",
    ]:
        style = ROLE_STYLE[role]
        svg.rect(x, y - 11, 30, 9, style["fill"], style["stroke"], 0.4, rx=1.2, opacity=style["opacity"] + 0.2)
        svg.text(x + 37, y - 3, style["label"], 10, "400", MUTED)
        x += 210 if role != "same-chromosome context" else 260


def render(mode: str, output: Path) -> None:
    events, segments, chrom_sizes = load_data()
    # Estimate page height from the actual number of lanes.
    y = TOP
    for event in events:
        tracks = build_tracks(event, event_segments(segments, event["event_id"]))
        y += sum(max(track_height(t), TRACK_H) for t in tracks) + TRACK_GAP * (len(tracks) - 1) + EVENT_GAP + 120
    height = math.ceil(y + 48)

    svg = SVG(PAGE_W, height)
    svg.text(42, 38, "Fig5 prototype: strict-path synteny/recombination schematic", 24, "700")
    svg.text(
        42,
        61,
        "Geometry comes from selected_segments.tsv; autosomal PHR rows are candidate exchange-compatible patches, not event-level validation.",
        13,
        "400",
        MUTED,
    )
    legend(svg, 82)

    y = TOP
    for event in events:
        segs = event_segments(segments, event["event_id"])
        y = draw_event(svg, event, segs, chrom_sizes, y, mode)
        svg.line(42, y - 30, PAGE_W - 42, y - 30, LIGHT, 0.8)
    svg.write(output)


def convert_pdf(svg_path: Path) -> tuple[bool, str]:
    pdf_path = svg_path.with_suffix(".pdf")
    rsvg = shutil.which("rsvg-convert")
    if rsvg:
        subprocess.run([rsvg, "-f", "pdf", "-o", str(pdf_path), str(svg_path)], check=True)
        return True, f"Converted with {rsvg}"
    inkscape = shutil.which("inkscape")
    if inkscape:
        subprocess.run([inkscape, str(svg_path), "--export-type=pdf", f"--export-filename={pdf_path}"], check=True)
        return True, f"Converted with {inkscape}"
    try:
        import cairosvg  # type: ignore
    except Exception:
        return False, "No rsvg-convert, inkscape, or Python cairosvg was available."
    cairosvg.svg2pdf(url=str(svg_path), write_to=str(pdf_path))
    return True, "Converted with Python cairosvg"


def main() -> None:
    render("full", FULL_SVG)
    render("focus", FOCUS_SVG)
    messages = []
    for svg_path in [FULL_SVG, FOCUS_SVG]:
        ok, message = convert_pdf(svg_path)
        messages.append(f"{svg_path.name}: {'PDF written' if ok else 'SVG only'} ({message})")
    status_path = HERE / "pdf_conversion_status.txt"
    status_path.write_text("\n".join(messages) + "\n")
    for message in messages:
        print(message)


if __name__ == "__main__":
    main()
