#!/usr/bin/env python3
"""
Render compact candidate Figure 5 panels for PAR1 and autosomal PHR events.

Usage:
    python3 paper_prep/_brainstorming/fig5_par1_phr_candidate_panels/plot_fig5_par1_phr_candidate_panels.py

The script reads existing WashU odgi-untangle recombination patch output and
writes the figure PDF/SVG plus panel_event_summary.tsv next to this file. It
does not run odgi, sweepGA, or any heavy discovery pipeline. Rendering uses only
the Python standard library so this can run on lightweight worktrees.
"""

from __future__ import annotations

import argparse
import csv
import html
import re
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


HERE = Path(__file__).resolve().parent
DEFAULT_PATCHES = Path(
    "/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv"
)
DEFAULT_SEGMENTS = (
    HERE.parent / "fig5_sweepga_1to1_redraw" / "conservative_segments.tsv"
)

FIG_BASE = "fig5_par1_phr_candidate_panels"
FULL_SPAN = 500_000
WIDTH = 1100
HEIGHT = 760
LEFT = 86
RIGHT = 32
PANEL_W = WIDTH - LEFT - RIGHT

PANELS = [
    {
        "panel": "A",
        "title": "A  PAR1 positive control: paternal chrX p maps to chrY p",
        "short_title": "PAR1 positive control",
        "label": "PAN027 paternal (hap2) vs PAN011 (father)",
        "pair": "PAN027_vs_PAN011",
        "query": "PAN027#2#chrX.paternal:12265-512264_chrX_parm",
        "primary_arms": {"chrYp"},
        "secondary_arms": set(),
        "low_conf_arms": set(),
        "interpretation": "Known male PAR1 X/Y recombination sanity check; not autosomal PHR evidence.",
        "phr_status": "PAR1 positive-control",
        "xlim": (0, FULL_SPAN),
        "callout": "Main chrYp block at chrX:12,265-155,863; strict nb=1 1:1 path; C15/C15.",
    },
    {
        "panel": "B",
        "title": "B  PAN027 autosomal candidate: chr9 q terminal patch mostly chr3 q",
        "short_title": "PAN027 chr9q candidate",
        "label": "PAN027 paternal (hap2) vs PAN011 (father)",
        "pair": "PAN027_vs_PAN011",
        "query": "PAN027#2#chr9.paternal:135704825-136204824_chr9_qarm",
        "primary_arms": {"chr3q"},
        "secondary_arms": {"chr15q", "chr16q"},
        "low_conf_arms": {"chr20q"},
        "interpretation": "Candidate terminal PHR exchange patch; not a clean full crossover.",
        "phr_status": "autosomal PHR candidate",
        "xlim": (280_000, FULL_SPAN),
        "callout": "Strict primary path switches near chr9:136,151,769; terminal tract is mostly chr3q with small side fragments.",
    },
    {
        "panel": "C",
        "title": "C  Independent PAN028 support: chr3 q terminal patch maps to chr9 q",
        "short_title": "PAN028 chr3q candidate",
        "label": "PAN028 maternal (hap1) vs PAN027 (mother)",
        "pair": "PAN028_vs_PAN027",
        "query": "PAN028#1#chr3.haplotype1:199233840-199733839_chr3_qarm",
        "primary_arms": {"chr9q"},
        "secondary_arms": {"chr7p", "chr16q", "chr20q"},
        "low_conf_arms": set(),
        "interpretation": "Independent candidate compatible with chr3q/chr9q C3 exchange.",
        "phr_status": "autosomal PHR candidate",
        "xlim": (220_000, FULL_SPAN),
        "callout": "Strict primary path has candidate side fragments across chr3:199,496,793-199,733,150; permissive chr9q patch calls are annotation only.",
    },
]

COLORS = {
    "background": "#d9d9d9",
    "same_chr_alt_hap": "#c4c4c4",
    "query_outline": "#8a8a8a",
    "primary_chrYp": "#2f6f9f",
    "primary_chr3q_chr9q": "#008a78",
    "secondary": "#a98a3a",
    "cross_community": "#9a7a8e",
    "low_conf": "#bdb7ad",
    "phr_band": "#f3f0e8",
    "axis": "#4a4a4a",
    "text": "#1f1f1f",
}


@dataclass(frozen=True)
class SourceInterval:
    sample: str
    hap: str
    chrom: str
    hap_label: str
    start: int
    end: int
    arm_name: str

    @property
    def chrom_short(self) -> str:
        return self.chrom

    @property
    def window(self) -> str:
        return f"{self.chrom_short}:{self.start:,}-{self.end:,}"

    def project(self, local_start: int, local_end: int) -> tuple[int, int]:
        return self.start + local_start, self.start + local_end

    def interval_label(self, local_start: int, local_end: int) -> str:
        start, end = self.project(local_start, local_end)
        return f"{self.chrom_short}:{start:,}-{end:,}"


@dataclass(frozen=True)
class Patch:
    label: str
    query: str
    query_chr: str
    query_arm: str
    patch_start: int
    patch_end: int
    patch_size: int
    ref_chr: str
    ref_arm: str
    ref_chrarm: str
    ref_hap: str
    mean_score: str
    min_score: str
    max_score: str
    n_segments: str
    is_interchr: bool
    pattern: str
    query_community: str
    ref_community: str
    community_status: str
    overlaps_phr: str
    has_phr: str
    phr_start: str
    phr_end: str

    @property
    def interval(self) -> str:
        return f"{self.patch_start}-{self.patch_end}"


@dataclass(frozen=True)
class Segment:
    pair: str
    transmission: str
    query_name: str
    query_arm: str
    query_start: int
    query_end: int
    query_length: int
    target_name: str
    target_hap: str
    target_arm: str
    strand: str
    identity: str
    jaccard: str
    nb: str
    interchromosomal: bool
    annotation: Patch | None = None

    @property
    def patch_start(self) -> int:
        return self.query_start

    @property
    def patch_end(self) -> int:
        return self.query_end

    @property
    def patch_size(self) -> int:
        return self.query_end - self.query_start

    @property
    def is_interchr(self) -> bool:
        return self.interchromosomal

    @property
    def ref_chrarm(self) -> str:
        return self.target_arm

    @property
    def ref_hap(self) -> str:
        return self.target_hap.split("#")[-1]

    @property
    def interval(self) -> str:
        return f"{self.query_start}-{self.query_end}"

    @property
    def mean_score(self) -> str:
        return self.annotation.mean_score if self.annotation else self.jaccard

    @property
    def min_score(self) -> str:
        return self.annotation.min_score if self.annotation else self.jaccard

    @property
    def max_score(self) -> str:
        return self.annotation.max_score if self.annotation else self.jaccard

    @property
    def n_segments(self) -> str:
        return self.annotation.n_segments if self.annotation else "1"

    @property
    def pattern(self) -> str:
        return self.annotation.pattern if self.annotation else "strict_1to1_primary_path"

    @property
    def query_community(self) -> str:
        return self.annotation.query_community if self.annotation else "not_available"

    @property
    def ref_community(self) -> str:
        return self.annotation.ref_community if self.annotation else "not_available"

    @property
    def community_status(self) -> str:
        return self.annotation.community_status if self.annotation else "not_available"

    @property
    def overlaps_phr(self) -> str:
        return self.annotation.overlaps_phr if self.annotation else "not_available"

    @property
    def has_phr(self) -> str:
        return self.annotation.has_phr if self.annotation else "not_available"

    @property
    def phr_start(self) -> str:
        return self.annotation.phr_start if self.annotation else "None"

    @property
    def phr_end(self) -> str:
        return self.annotation.phr_end if self.annotation else "None"


class Canvas:
    def rect(self, x: float, y: float, w: float, h: float, fill: str, stroke: str = "none", sw: float = 0.0, opacity: float = 1.0) -> None:
        raise NotImplementedError

    def line(self, x1: float, y1: float, x2: float, y2: float, stroke: str, sw: float = 1.0, opacity: float = 1.0) -> None:
        raise NotImplementedError

    def text(self, x: float, y: float, value: str, size: float, fill: str = "#222222", weight: str = "normal", anchor: str = "start") -> None:
        raise NotImplementedError


class SvgCanvas(Canvas):
    def __init__(self, width: int, height: int) -> None:
        self.width = width
        self.height = height
        self.ops: list[str] = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
            '<style>text{font-family:DejaVu Sans,Arial,sans-serif;letter-spacing:0}</style>',
        ]

    def rect(self, x: float, y: float, w: float, h: float, fill: str, stroke: str = "none", sw: float = 0.0, opacity: float = 1.0) -> None:
        self.ops.append(
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="{fill}" stroke="{stroke}" stroke-width="{sw:.2f}" opacity="{opacity:.3f}"/>'
        )

    def line(self, x1: float, y1: float, x2: float, y2: float, stroke: str, sw: float = 1.0, opacity: float = 1.0) -> None:
        self.ops.append(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{stroke}" stroke-width="{sw:.2f}" opacity="{opacity:.3f}"/>'
        )

    def text(self, x: float, y: float, value: str, size: float, fill: str = "#222222", weight: str = "normal", anchor: str = "start") -> None:
        self.ops.append(
            f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size:.2f}" fill="{fill}" font-weight="{weight}" text-anchor="{anchor}">{html.escape(value)}</text>'
        )

    def write(self, path: Path) -> None:
        path.write_text("\n".join(self.ops + ["</svg>\n"]))


class PdfCanvas(Canvas):
    def __init__(self, width: int, height: int) -> None:
        self.width = float(width)
        self.height = float(height)
        self.ops: list[str] = []

    def _y(self, y: float) -> float:
        return self.height - y

    def rect(self, x: float, y: float, w: float, h: float, fill: str, stroke: str = "none", sw: float = 0.0, opacity: float = 1.0) -> None:
        # PDF opacity needs an ExtGState; use solid colors for portability.
        if fill != "none":
            self.ops.append(f"{pdf_color(fill)} rg")
        if stroke != "none" and sw > 0:
            self.ops.append(f"{pdf_color(stroke)} RG {sw:.3f} w")
        self.ops.append(f"{x:.3f} {self._y(y + h):.3f} {w:.3f} {h:.3f} re")
        if fill != "none" and stroke != "none" and sw > 0:
            self.ops.append("B")
        elif fill != "none":
            self.ops.append("f")
        elif stroke != "none" and sw > 0:
            self.ops.append("S")

    def line(self, x1: float, y1: float, x2: float, y2: float, stroke: str, sw: float = 1.0, opacity: float = 1.0) -> None:
        self.ops.append(f"{pdf_color(stroke)} RG {sw:.3f} w {x1:.3f} {self._y(y1):.3f} m {x2:.3f} {self._y(y2):.3f} l S")

    def text(self, x: float, y: float, value: str, size: float, fill: str = "#222222", weight: str = "normal", anchor: str = "start") -> None:
        safe = value.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")
        approx_w = len(value) * size * 0.52
        if anchor == "middle":
            x -= approx_w / 2
        elif anchor == "end":
            x -= approx_w
        font = "/F2" if weight == "bold" else "/F1"
        self.ops.append(f"BT {font} {size:.2f} Tf {pdf_color(fill)} rg {x:.3f} {self._y(y):.3f} Td ({safe}) Tj ET")

    def write(self, path: Path) -> None:
        stream = "\n".join(self.ops).encode("latin-1", "replace")
        objects = [
            b"<< /Type /Catalog /Pages 2 0 R >>",
            b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
            (
                f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {self.width:.3f} {self.height:.3f}] "
                f"/Resources << /Font << /F1 4 0 R /F2 5 0 R >> >> /Contents 6 0 R >>"
            ).encode(),
            b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
            b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>",
            b"<< /Length " + str(len(stream)).encode() + b" >>\nstream\n" + stream + b"\nendstream",
        ]
        data = bytearray(b"%PDF-1.4\n")
        offsets = [0]
        for idx, obj in enumerate(objects, start=1):
            offsets.append(len(data))
            data.extend(f"{idx} 0 obj\n".encode())
            data.extend(obj)
            data.extend(b"\nendobj\n")
        xref = len(data)
        data.extend(f"xref\n0 {len(objects) + 1}\n0000000000 65535 f \n".encode())
        for offset in offsets[1:]:
            data.extend(f"{offset:010d} 00000 n \n".encode())
        data.extend(
            f"trailer << /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n".encode()
        )
        path.write_bytes(bytes(data))


def pdf_color(hex_color: str) -> str:
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16) / 255
    g = int(hex_color[2:4], 16) / 255
    b = int(hex_color[4:6], 16) / 255
    return f"{r:.4f} {g:.4f} {b:.4f}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patches", type=Path, default=DEFAULT_PATCHES)
    parser.add_argument("--segments", type=Path, default=DEFAULT_SEGMENTS)
    parser.add_argument("--out-dir", type=Path, default=HERE)
    return parser.parse_args()


QUERY_RE = re.compile(
    r"^(?P<sample>[^#]+)#(?P<hap>[^#]+)#(?P<chrom>[^.:_]+)\.(?P<hap_label>[^:]+):(?P<start>\d+)-(?P<end>\d+)_(?P<arm>.+)$"
)


def parse_source_interval(name: str) -> SourceInterval:
    match = QUERY_RE.match(name)
    if not match:
        raise ValueError(f"cannot parse source interval from name: {name}")
    return SourceInterval(
        sample=match.group("sample"),
        hap=match.group("hap"),
        chrom=match.group("chrom"),
        hap_label=match.group("hap_label"),
        start=int(match.group("start")),
        end=int(match.group("end")),
        arm_name=match.group("arm"),
    )


def read_patches(path: Path) -> list[Patch]:
    if not path.exists():
        raise FileNotFoundError(f"patch table not found: {path}")
    rows: list[Patch] = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(
                Patch(
                    label=row["label"],
                    query=row["query"],
                    query_chr=row["query_chr"],
                    query_arm=row["query_arm"],
                    patch_start=int(row["patch_start"]),
                    patch_end=int(row["patch_end"]),
                    patch_size=int(row["patch_size"]),
                    ref_chr=row["ref_chr"],
                    ref_arm=row["ref_arm"],
                    ref_chrarm=row["ref_chrarm"],
                    ref_hap=row["ref_hap"],
                    mean_score=row["mean_score"],
                    min_score=row["min_score"],
                    max_score=row["max_score"],
                    n_segments=row["n_segments"],
                    is_interchr=row["is_interchr"] == "True",
                    pattern=row["pattern"],
                    query_community=row["query_community"],
                    ref_community=row["ref_community"],
                    community_status=row["community_status"],
                    overlaps_phr=row["overlaps_phr"],
                    has_phr=row["has_phr"],
                    phr_start=row["phr_start"],
                    phr_end=row["phr_end"],
                )
            )
    return rows


def patch_index(patches: list[Patch]) -> dict[tuple[str, str, int, int, str, str], Patch]:
    return {
        (p.query, p.patch_start, p.patch_end, p.ref_chrarm, p.ref_hap): p
        for p in patches
    }


def read_segments(path: Path, patches: list[Patch]) -> list[Segment]:
    if not path.exists():
        raise FileNotFoundError(f"conservative segment table not found: {path}")
    idx = patch_index(patches)
    rows: list[Segment] = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["nb"] != "1":
                continue
            target_hap_num = row["target_hap"].split("#")[-1]
            # patches.tsv uses the human-readable transmission label rather than
            # pair id, so match on query/interval/donor when available.
            ann = idx.get(
                (
                    row["query_name"],
                    int(row["query_start"]),
                    int(row["query_end"]),
                    row["target_arm"],
                    target_hap_num,
                )
            )
            rows.append(
                Segment(
                    pair=row["pair"],
                    transmission=row["transmission"],
                    query_name=row["query_name"],
                    query_arm=row["query_arm"],
                    query_start=int(row["query_start"]),
                    query_end=int(row["query_end"]),
                    query_length=int(row["query_length"]),
                    target_name=row["target_name"],
                    target_hap=row["target_hap"],
                    target_arm=row["target_arm"],
                    strand=row["strand"],
                    identity=row["identity"],
                    jaccard=row["jaccard"],
                    nb=row["nb"],
                    interchromosomal=row["interchromosomal"] == "1",
                    annotation=ann,
                )
            )
    return rows


def panel_rows(segments: list[Segment], panel: dict) -> list[Segment]:
    rows = [p for p in segments if p.pair == panel["pair"] and p.query_name == panel["query"]]
    if not rows:
        raise RuntimeError(f"no patch rows for panel {panel['panel']}: {panel['query']}")
    return sorted(rows, key=lambda p: p.patch_start)


def event_class(patch: Segment, panel: dict) -> str:
    if not patch.is_interchr:
        return "same-chromosome background"
    if patch.ref_chrarm in panel["low_conf_arms"] or patch.patch_size < 1_000 or float(patch.min_score) < 0.90:
        return "low-confidence"
    if panel["phr_status"].startswith("PAR1"):
        return "positive-control"
    if patch.ref_chrarm in panel["primary_arms"] and patch.community_status in {"within_community", "not_available"}:
        return "autosomal-candidate"
    if patch.community_status != "within_community" or patch.ref_chrarm in panel["secondary_arms"]:
        return "secondary fragment"
    return "autosomal-candidate"


def patch_style(patch: Segment, panel: dict) -> tuple[str, str]:
    cls = event_class(patch, panel)
    if cls == "same-chromosome background":
        fill = COLORS["same_chr_alt_hap"] if patch.ref_hap != "2" else COLORS["background"]
        return fill, "#ffffff"
    if cls == "positive-control":
        return COLORS["primary_chrYp"], "#1d4668"
    if cls == "autosomal-candidate":
        return COLORS["primary_chr3q_chr9q"], "#00594e"
    if cls == "secondary fragment":
        fill = COLORS["cross_community"] if patch.community_status != "within_community" else COLORS["secondary"]
        return fill, "#6e5a2a"
    return COLORS["low_conf"], "#8e877d"


def selected_for_summary(patch: Segment, panel: dict) -> bool:
    if patch.is_interchr:
        return True
    if panel["panel"] == "A":
        return patch.patch_start in {143431, 194850, 331663}
    return False


def write_summary(segments: list[Segment], out_path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for panel in PANELS:
        query_src = parse_source_interval(panel["query"])
        for patch in panel_rows(segments, panel):
            if not selected_for_summary(patch, panel):
                continue
            target_src = parse_source_interval(patch.target_name)
            rows.append(
                {
                    "panel": panel["panel"],
                    "panel_title": panel["short_title"],
                    "transmission": patch.transmission,
                    "query": patch.query_name,
                    "coordinate_convention": "0-based half-open native assembly coordinates parsed from source sequence names",
                    "query_arm": patch.query_arm,
                    "local_interval_bp": patch.interval,
                    "display_query_interval_native_0based": query_src.interval_label(patch.query_start, patch.query_end),
                    "display_query_start_native_0based": str(query_src.project(patch.query_start, patch.query_end)[0]),
                    "display_query_end_native_0based": str(query_src.project(patch.query_start, patch.query_end)[1]),
                    "donor_arm": patch.ref_chrarm,
                    "donor_hap": patch.target_hap,
                    "donor_source_window_native_0based": target_src.window,
                    "donor_segment_interval_native_0based": "not_recovered_from_conservative_segment_table",
                    "length_bp": str(patch.patch_size),
                    "identity_or_mean_score": patch.identity,
                    "jaccard_or_min_score": patch.jaccard,
                    "max_score_from_patch_annotation": patch.max_score,
                    "n_segments": patch.n_segments,
                    "drawing_source": "conservative_segments.tsv nb=1 sweepGA 1:1 no-scaffold primary path",
                    "pattern": patch.pattern,
                    "query_community": patch.query_community,
                    "ref_community": patch.ref_community,
                    "community_status": patch.community_status,
                    "overlaps_phr": patch.overlaps_phr,
                    "has_phr": patch.has_phr,
                    "phr_or_par_status": panel["phr_status"],
                    "event_class": event_class(patch, panel),
                    "interpretation_boundary": panel["interpretation"],
                }
            )
    fieldnames = list(rows[0].keys())
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def xmap(value: int, xmin: int, xmax: int) -> float:
    return LEFT + (value - xmin) / (xmax - xmin) * PANEL_W


def wrap_text(c: Canvas, x: float, y: float, text: str, width: int, size: float, fill: str = "#555555", line_height: float = 10) -> None:
    for i, line in enumerate(textwrap.wrap(text, width=width)):
        c.text(x, y + i * line_height, line, size=size, fill=fill)


def draw_axis(c: Canvas, y: float, xmin: int, xmax: int, source: SourceInterval) -> None:
    c.line(LEFT, y + 98, LEFT + PANEL_W, y + 98, COLORS["axis"], 0.7)
    tick_start = (xmin // 50_000) * 50_000
    for tick in range(tick_start, xmax + 1, 50_000):
        if tick < xmin:
            continue
        x = xmap(tick, xmin, xmax)
        c.line(x, y + 94, x, y + 101, COLORS["axis"], 0.6)
        c.text(x, y + 113, f"{source.start + tick:,}", 7.2, "#555555", anchor="middle")
    c.text(LEFT + PANEL_W, y + 126, f"{source.chrom_short} native bp (0-based half-open)", 7.2, "#555555", anchor="end")


def draw_panel(c: Canvas, rows: list[Segment], panel: dict, y: float) -> None:
    xmin, xmax = panel["xlim"]
    source = parse_source_interval(panel["query"])
    c.text(LEFT, y, panel["title"], 11, COLORS["text"], weight="bold")
    c.text(LEFT, y + 18, panel["label"], 8.2, "#555555")
    c.text(LEFT, y + 30, f"Displayed query window: {source.interval_label(xmin, xmax)}", 7.4, "#555555")
    wrap_text(c, LEFT + 545, y + 18, panel["callout"], 72, 7.7, "#555555", 9)

    phr_starts = [int(r.phr_start) for r in rows if r.phr_start != "None"]
    phr_ends = [int(r.phr_end) for r in rows if r.phr_end != "None"]
    if phr_starts and phr_ends:
        phr_start = max(xmin, min(phr_starts))
        phr_end = min(xmax, max(phr_ends))
        if phr_end > phr_start:
            c.rect(xmap(phr_start, xmin, xmax), y + 58, xmap(phr_end, xmin, xmax) - xmap(phr_start, xmin, xmax), 30, COLORS["phr_band"])
            c.text(xmap(phr_start, xmin, xmax), y + 54, "PHR/PAR interval", 7, "#766f63")

    c.rect(LEFT, y + 67, PANEL_W, 14, "#eeeeee", COLORS["query_outline"], 0.5)
    label_toggle = False
    for patch in rows:
        start = max(patch.patch_start, xmin)
        end = min(patch.patch_end, xmax)
        if end <= start:
            continue
        fill, edge = patch_style(patch, panel)
        px = xmap(start, xmin, xmax)
        pw = max(0.9, xmap(end, xmin, xmax) - px)
        if patch.is_interchr:
            py, ph = y + 63, 22
        else:
            py, ph = y + 68, 12
        c.rect(px, py, pw, ph, fill, edge, 0.5)
        cls = event_class(patch, panel)
        if patch.is_interchr and (patch.patch_size >= 1_000 or cls != "low-confidence"):
            cx = px + pw / 2
            ty = y + (48 if not label_toggle else 38)
            label_toggle = not label_toggle
            label = f"{patch.ref_chrarm} h{patch.ref_hap}"
            if cls == "low-confidence":
                label += " low"
            elif patch.community_status != "within_community":
                label += " side"
            c.line(cx, py, cx, ty + 4, edge, 0.5)
            c.text(cx, ty, label, 7, edge, anchor="middle")

    wrap_text(c, LEFT, y + 126, panel["interpretation"], 105, 7.8, "#555555", 9)
    draw_axis(c, y, xmin, xmax, source)


def grouped_summary_line(rows: Iterable[dict[str, str]], event_class_name: str) -> str:
    entries = [r for r in rows if r["event_class"] == event_class_name]
    if not entries:
        return ""
    parts = [
        f"{r['donor_arm']} {r['display_query_interval_native_0based']} ({r['community_status']})"
        for r in entries
    ]
    return "; ".join(parts)


def draw_overview(c: Canvas, y: float, summary_rows: list[dict[str, str]]) -> None:
    c.text(LEFT, y, "D  Selected finite candidate set", 11, COLORS["text"], weight="bold")
    c.text(LEFT, y + 18, "Deliberately selected review set; not a genome-wide rediscovery run.", 8.2, "#555555")
    col_x = [LEFT, LEFT + 54, LEFT + 205, LEFT + 590, LEFT + 840]
    col_w = [48, 145, 375, 240, 190]
    headers = ["panel", "role", "primary evidence", "side fragments", "low-confidence"]
    row_h = 40
    table_y = y + 36
    for i, header in enumerate(headers):
        c.rect(col_x[i], table_y, col_w[i], 22, "#efefef", "#cfcfcf", 0.4)
        c.text(col_x[i] + 4, table_y + 15, header, 7.3, COLORS["text"], weight="bold")

    table_rows = []
    for panel in PANELS:
        rows = [r for r in summary_rows if r["panel"] == panel["panel"]]
        primary = grouped_summary_line(rows, "positive-control") or grouped_summary_line(rows, "autosomal-candidate")
        secondary = grouped_summary_line(rows, "secondary fragment") or "none emphasized"
        low = grouped_summary_line(rows, "low-confidence") or "none emphasized"
        table_rows.append([panel["panel"], panel["short_title"], primary, secondary, low])

    for ridx, row in enumerate(table_rows):
        ry = table_y + 22 + ridx * row_h
        for cidx, value in enumerate(row):
            c.rect(col_x[cidx], ry, col_w[cidx], row_h, "#ffffff", "#d6d6d6", 0.35)
            lines = textwrap.wrap(value, width=max(6, int(col_w[cidx] / 5.0)))[:4]
            for lidx, line in enumerate(lines):
                c.text(col_x[cidx] + 4, ry + 11 + lidx * 8.5, line, 6.7, "#333333")

    legend = [
        ("same-chromosome background", COLORS["background"]),
        ("PAR1 chrY donor", COLORS["primary_chrYp"]),
        ("chr3/chr9 candidate donor", COLORS["primary_chr3q_chr9q"]),
        ("side/cross-community", COLORS["cross_community"]),
        ("low-confidence tail", COLORS["low_conf"]),
    ]
    lx = LEFT
    ly = y + 188
    for label, color in legend:
        c.rect(lx, ly, 15, 12, color, "#666666", 0.3)
        c.text(lx + 20, ly + 10, label, 7.2, "#444444")
        lx += 188


def render(canvas: Canvas, selected: dict[str, list[Segment]], summary_rows: list[dict[str, str]]) -> None:
    canvas.text(
        16,
        26,
        "Candidate Fig. 5 asset pack: PAR1 control and autosomal PHR-compatible pedigree patches",
        13,
        COLORS["text"],
        weight="bold",
    )
    canvas.text(
        16,
        44,
        "Strict nb=1 sweepGA 1:1 no-scaffold primary path; axes are native assembly coordinates parsed from source names.",
        8,
        "#666666",
    )
    panel_ys = [76, 222, 368]
    for panel, y in zip(PANELS, panel_ys):
        draw_panel(canvas, selected[panel["panel"]], panel, y)
    draw_overview(canvas, 532, summary_rows)


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(Path.cwd().resolve()))
    except ValueError:
        return str(path)


def write_readme(out_dir: Path, patch_path: Path, segment_path: Path) -> None:
    patch_display = display_path(patch_path)
    segment_display = display_path(segment_path)
    readme = f"""# Fig5 PAR1/PHR Candidate Panels

This directory contains an experimental candidate asset pack for reviewing a
compact Figure 5-style pedigree recombination panel. It does not modify
`submission/paper.tex`, `submission/fig/MainFigures/Fig5_pedigree_untangle.pdf`,
or any bibliography file.

## Outputs

- `{FIG_BASE}.pdf`
- `{FIG_BASE}.svg`
- `panel_event_summary.tsv`
- `plot_fig5_par1_phr_candidate_panels.py`

## Provenance

Drawing input:

`{segment_display}`

The plotted rectangles come from `conservative_segments.tsv`, the strict
`nb=1` sweepGA 1:1 no-scaffold primary-path table produced by the companion
redraw task. This correction does not plot nth-best, multimap, permissive
secondary, or non-primary PAF rows.

Annotation/provenance input:

`{patch_display}`

This table is the existing WashU `odgi untangle` recombination patch output.
The script uses it only to recover community labels, score summaries, PHR/PAR
status, and interpretive labels when a conservative segment has the same query,
interval, donor arm, and donor haplotype.

The script filters three pre-selected, review-facing examples:

- PAN027 paternal hap2 vs PAN011 father, chrX p PAR1 positive control.
- PAN027 paternal hap2 vs PAN011 father, chr9 q terminal autosomal PHR candidate.
- PAN028 maternal hap1 vs PAN027 mother, chr3 q independent autosomal PHR candidate.

## Coordinate Convention

Axes, callouts, and `panel_event_summary.tsv` use 0-based half-open native
assembly coordinates parsed from source names such as
`PAN027#2#chr9.paternal:135704825-136204824_chr9_qarm`. A local segment
`[a,b)` is displayed as `chr:(start+a)-(start+b)`. For example, local
`[446944,472441)` in the PAN027 paternal chr9q source window is displayed as
`chr9:136,151,769-136,177,266`.

No CHM13 projection or liftover table is used here. The source names are native
sample assembly windows, so the figure deliberately labels the coordinates as
native assembly coordinates, not CHM13 coordinates.

Donor/source target names are also native assembly windows. The conservative
summary table does not carry exact target-side segment start/end fields, so
`panel_event_summary.tsv` records each donor source window but marks exact donor
segment intervals as not recovered in this lightweight presentation correction.

## Regeneration

From the repository root:

```bash
python3 paper_prep/_brainstorming/fig5_par1_phr_candidate_panels/plot_fig5_par1_phr_candidate_panels.py
```

Optional explicit input/output:

```bash
python3 paper_prep/_brainstorming/fig5_par1_phr_candidate_panels/plot_fig5_par1_phr_candidate_panels.py \\
  --segments {segment_display} \\
  --patches {patch_display} \\
  --out-dir paper_prep/_brainstorming/fig5_par1_phr_candidate_panels
```

## Interpretation Boundaries

Panel A is a known male PAR1 X/Y positive control and should be kept separate
from the autosomal PHR interpretation.

Panels B and C are candidate event-level examples compatible with chr3q/chr9q
C3 PHR exchange. They are not presented as clean full crossovers. In Panel B,
the terminal tract is mostly chr3q, while the chr15q segment is marked as a
smaller side fragment within the single selected 1:1 path, and the tiny chr20q
tail is treated as low-confidence. In Panel C, strict primary-path drawing
shows chr7p, chr16q, and chr20q side fragments in the selected chr3q window;
previous permissive chr9q patch calls are not drawn as alternate alignments.

The optional acrocentric/known-system panel was intentionally omitted here:
the finite four-panel asset keeps PAR1 plus the two autosomal chr3q/chr9q
examples legible and avoids a repetitive p-arm-dominated panel.

The asset is for review and presentation triage before any manuscript
integration.
"""
    (out_dir / "README.md").write_text(readme)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    patches = read_patches(args.patches)
    segments = read_segments(args.segments, patches)
    selected = {panel["panel"]: panel_rows(segments, panel) for panel in PANELS}
    summary_rows = write_summary(segments, args.out_dir / "panel_event_summary.tsv")

    svg = SvgCanvas(WIDTH, HEIGHT)
    render(svg, selected, summary_rows)
    svg.write(args.out_dir / f"{FIG_BASE}.svg")

    pdf = PdfCanvas(WIDTH, HEIGHT)
    render(pdf, selected, summary_rows)
    pdf.write(args.out_dir / f"{FIG_BASE}.pdf")

    write_readme(args.out_dir, args.patches, args.segments)
    print(f"wrote {args.out_dir / (FIG_BASE + '.svg')}")
    print(f"wrote {args.out_dir / (FIG_BASE + '.pdf')}")
    print(f"wrote {args.out_dir / 'panel_event_summary.tsv'}")
    print(f"wrote {args.out_dir / 'README.md'}")


if __name__ == "__main__":
    main()
