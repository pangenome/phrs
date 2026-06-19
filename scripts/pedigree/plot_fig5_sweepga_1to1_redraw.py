#!/usr/bin/env python3
"""
Regenerate a Fig. 5-style inspection redraw from conservative native PAF mappings.

Pipeline per WashU pedigree pair:

  native odgi untangle m1000 n4 PAF
    -> keep only rows tagged nb:i:1
    -> sweepga --num-mappings 1:1 --scaffold-jump 0
    -> plot conservative mappings

Large PAF intermediates are written outside the repository by default. Compact
tables and author-facing SVG/PDF outputs are written under paper_prep/_brainstorming.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import math
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]
NATIVE_DIR = Path("/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm")
EXTERNAL_DIR = Path("/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_sweepga_1to1_noscaffold")
OUT_DIR = REPO_ROOT / "paper_prep" / "_brainstorming" / "fig5_sweepga_1to1_redraw"
SWEEPGA = Path("/moosefs/erikg/sweepga/target/release/sweepga")
PHR_INTERVALS = Path("/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/all-vs-all.1Mb.p95.id95.len.tsv")

PAIR_INFO = {
    "PAN027_vs_PAN010": {
        "label": "PAN027 maternal hap1 <- PAN010 mother",
        "query_sample": "PAN027",
        "query_hap": "1",
        "target_sample": "PAN010",
    },
    "PAN027_vs_PAN011": {
        "label": "PAN027 paternal hap2 <- PAN011 father",
        "query_sample": "PAN027",
        "query_hap": "2",
        "target_sample": "PAN011",
    },
    "PAN028_vs_PAN027": {
        "label": "PAN028 maternal hap1 <- PAN027 mother",
        "query_sample": "PAN028",
        "query_hap": "1",
        "target_sample": "PAN027",
    },
}

ARM_ORDER = [f"chr{i}{arm}" for i in range(1, 23) for arm in ("p", "q")] + [
    "chrXp",
    "chrXq",
    "chrYp",
    "chrYq",
]
CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
PLOT_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
PATH_RE = re.compile(r"^(PAN\d+)#([12])#.*:(\d+)-(\d+)_(chr(?:[0-9]+|X|Y))_([pq])arm$")
TAG_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]):[A-Za-z]:(.*)$")
INTERCHR_PALETTE = {
    "chr10p": "#e73535",
    "chr13p": "#f05a28",
    "chr15p": "#3f4bad",
    "chr15q": "#c8d431",
    "chr19p": "#2e91d6",
    "chr21p": "#ff8a66",
    "chr22p": "#a88f86",
    "chr9q": "#008a78",
}


@dataclass(frozen=True)
class Segment:
    pair: str
    qname: str
    qarm: str
    qlen: int
    paf_qlen: int
    qstart: int
    qend: int
    strand: str
    tname: str
    target_hap: str
    target_arm: str
    target_sample: str
    matches: int
    block_len: int
    identity: float
    jc: float
    nb: int

    @property
    def length(self) -> int:
        return max(0, self.qend - self.qstart)

    @property
    def interchromosomal(self) -> bool:
        return chrom_of(self.qarm) != chrom_of(self.target_arm)


@dataclass(frozen=True)
class PhrInterval:
    qname: str
    start: int
    end: int
    chrs_involved: str
    arms_involved: str


class PdfWriter:
    """Tiny vector PDF writer for rectangles, lines, and simple text."""

    def __init__(self, width: float, height: float) -> None:
        self.width = width
        self.height = height
        self.ops: List[str] = []

    def _y(self, y: float) -> float:
        return self.height - y

    def rect(
        self,
        x: float,
        y: float,
        w: float,
        h: float,
        fill: str | None = None,
        stroke: str | None = None,
        sw: float = 0.4,
    ) -> None:
        if fill:
            self.ops.append(f"{pdf_color(fill)} rg")
        if stroke:
            self.ops.append(f"{pdf_color(stroke)} RG {sw:.3f} w")
        self.ops.append(f"{x:.3f} {self._y(y + h):.3f} {w:.3f} {h:.3f} re")
        if fill and stroke:
            self.ops.append("B")
        elif fill:
            self.ops.append("f")
        elif stroke:
            self.ops.append("S")

    def line(self, x1: float, y1: float, x2: float, y2: float, stroke: str = "#999999", sw: float = 0.3) -> None:
        self.ops.append(f"{pdf_color(stroke)} RG {sw:.3f} w {x1:.3f} {self._y(y1):.3f} m {x2:.3f} {self._y(y2):.3f} l S")

    def text(self, x: float, y: float, value: str, size: float = 8, fill: str = "#222222") -> None:
        safe = value.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")
        self.ops.append(f"BT /F1 {size:.2f} Tf {pdf_color(fill)} rg {x:.3f} {self._y(y):.3f} Td ({safe}) Tj ET")

    def write(self, path: Path) -> None:
        stream = "\n".join(self.ops).encode("latin-1", "replace")
        objects = [
            b"<< /Type /Catalog /Pages 2 0 R >>",
            b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
            (
                f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {self.width:.3f} {self.height:.3f}] "
                f"/Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>"
            ).encode(),
            b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-dir", type=Path, default=NATIVE_DIR)
    parser.add_argument("--external-dir", type=Path, default=EXTERNAL_DIR)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--sweepga", type=Path, default=SWEEPGA)
    parser.add_argument("--phr-intervals", type=Path, default=PHR_INTERVALS)
    parser.add_argument("--force", action="store_true", help="Regenerate external sweepGA PAF outputs even if present")
    parser.add_argument("--skip-sweepga", action="store_true", help="Use existing external PAF outputs only")
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def path_info(name: str) -> Tuple[str, str, str, int] | None:
    match = PATH_RE.match(name)
    if not match:
        return None
    sample, hap, start, end, chrom, arm = match.groups()
    span = int(end) - int(start) + 1
    if span <= 0:
        return None
    return sample, hap, f"{chrom}{arm}", span


def chrom_of(chrarm: str) -> str:
    return chrarm[:-1]


def arm_sort_key(chrarm: str) -> Tuple[int, int]:
    chrom = chrom_of(chrarm)
    if chrom == "chrX":
        c = 23
    elif chrom == "chrY":
        c = 24
    else:
        c = int(chrom[3:])
    return c, 0 if chrarm.endswith("p") else 1


def parse_tags(fields: Sequence[str]) -> Dict[str, str]:
    tags: Dict[str, str] = {}
    for field in fields:
        match = TAG_RE.match(field)
        if match:
            tags[match.group(1)] = match.group(2)
    return tags


def parse_float(value: str, default: float = math.nan) -> float:
    try:
        return float(value)
    except ValueError:
        return default


def read_paf(path: Path, pair: str, require_nb1: bool = False) -> List[Segment]:
    info = PAIR_INFO[pair]
    rows: List[Segment] = []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qinfo = path_info(fields[0])
            tinfo = path_info(fields[5])
            if not qinfo or not tinfo:
                continue
            qsample, qhap, qarm, query_span = qinfo
            tsample, thap, tarm, _target_span = tinfo
            if qsample != info["query_sample"] or qhap != info["query_hap"] or tsample != info["target_sample"]:
                continue
            tags = parse_tags(fields[12:])
            nb = int(tags.get("nb", "1"))
            if require_nb1 and nb != 1:
                continue
            qstart = int(fields[2])
            qend = int(fields[3])
            if qstart < 0 or qend <= qstart or qend > query_span:
                continue
            rows.append(
                Segment(
                    pair=pair,
                    qname=fields[0],
                    qarm=qarm,
                    qlen=query_span,
                    paf_qlen=int(fields[1]),
                    qstart=qstart,
                    qend=qend,
                    strand=fields[4],
                    tname=fields[5],
                    target_hap=f"{tsample}#{thap}",
                    target_arm=tarm,
                    target_sample=tsample,
                    matches=int(fields[9]),
                    block_len=max(1, int(fields[10])),
                    identity=parse_float(tags.get("id", str(100 * int(fields[9]) / max(1, int(fields[10]))))),
                    jc=parse_float(tags.get("jc", "nan")),
                    nb=nb,
                )
            )
    return rows


def load_phr_intervals(path: Path, query_names: set[str]) -> Dict[str, List[PhrInterval]]:
    if not path.exists():
        raise FileNotFoundError(path)
    intervals: Dict[str, List[PhrInterval]] = {name: [] for name in query_names}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"seq", "region_start", "region_end", "chrs_involved", "arms_involved"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise RuntimeError(f"PHR interval table missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            qname = row["seq"]
            if qname not in query_names:
                continue
            if row["region_start"] == "." or row["region_end"] == ".":
                continue
            start = int(row["region_start"])
            end = int(row["region_end"])
            if end <= start:
                continue
            intervals[qname].append(
                PhrInterval(
                    qname=qname,
                    start=max(0, start),
                    end=min(500000, end),
                    chrs_involved=row["chrs_involved"],
                    arms_involved=row["arms_involved"],
                )
            )
    return intervals


def native_paf_path(native_dir: Path, pair: str) -> Path:
    return native_dir / f"{pair}.e50000.m1000.j0.8.n4.paf.gz"


def plain_sweep_path(external_dir: Path, pair: str) -> Path:
    return external_dir / f"{pair}.e50000.m1000.j0.8.n4.sweepga_1to1_noscaffold.paf"


def native_uncompressed_path(external_dir: Path, pair: str) -> Path:
    return external_dir / f"{pair}.e50000.m1000.j0.8.n4.native.paf"


def nb1_input_path(external_dir: Path, pair: str) -> Path:
    return external_dir / f"{pair}.e50000.m1000.j0.8.n4.native_nb1.paf"


def conservative_sweep_path(external_dir: Path, pair: str) -> Path:
    return external_dir / f"{pair}.e50000.m1000.j0.8.n4.native_nb1.sweepga_1to1_noscaffold.paf"


def write_nb1_input(native_paf: Path, out_paf: Path) -> int:
    out_paf.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open_text(native_paf) as src, out_paf.open("w") as dst:
        for line in src:
            if "nb:i:1" in line.rstrip("\n").split("\t")[12:]:
                dst.write(line)
                count += 1
    return count


def write_uncompressed_input(native_paf: Path, out_paf: Path) -> int:
    out_paf.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open_text(native_paf) as src, out_paf.open("w") as dst:
        for line in src:
            dst.write(line)
            if line.strip() and not line.startswith("#"):
                count += 1
    return count


def run_sweepga(sweepga: Path, input_paf: Path, output_paf: Path, force: bool) -> None:
    if output_paf.exists() and not force:
        return
    if not sweepga.exists():
        raise FileNotFoundError(f"sweepga binary not found: {sweepga}")
    output_paf.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(sweepga),
        "--num-mappings",
        "1:1",
        "--scaffold-jump",
        "0",
        "--output-file",
        str(output_paf),
        str(input_paf),
    ]
    subprocess.run(cmd, check=True)


def ensure_outputs(args: argparse.Namespace) -> Tuple[List[dict], Dict[str, Path]]:
    summary_rows: List[dict] = []
    conservative_paths: Dict[str, Path] = {}
    for pair in PAIR_INFO:
        native = native_paf_path(args.native_dir, pair)
        if not native.exists():
            raise FileNotFoundError(native)
        native_uncompressed = native_uncompressed_path(args.external_dir, pair)
        nb1_in = nb1_input_path(args.external_dir, pair)
        plain = plain_sweep_path(args.external_dir, pair)
        conservative = conservative_sweep_path(args.external_dir, pair)
        if args.force or not native_uncompressed.exists():
            write_uncompressed_input(native, native_uncompressed)
        if args.force or not nb1_in.exists():
            write_nb1_input(native, nb1_in)
        if not args.skip_sweepga:
            run_sweepga(args.sweepga, native_uncompressed, plain, args.force)
            run_sweepga(args.sweepga, nb1_in, conservative, args.force)
        conservative_paths[pair] = conservative

        native_rows = read_paf(native, pair)
        nb1_rows = read_paf(nb1_in, pair)
        plain_rows = read_paf(plain, pair) if plain.exists() else []
        conservative_rows = read_paf(conservative, pair) if conservative.exists() else []
        for stage, rows, path in [
            ("native_n4", native_rows, native),
            ("native_nb1_prefilter", nb1_rows, nb1_in),
            ("plain_sweepga_n4_1to1_noscaffold", plain_rows, plain),
            ("conservative_nb1_sweepga_1to1_noscaffold", conservative_rows, conservative),
        ]:
            nb_values = sorted({r.nb for r in rows})
            summary_rows.append(
                {
                    "pair": pair,
                    "transmission": PAIR_INFO[pair]["label"],
                    "stage": stage,
                    "rows": len(rows),
                    "nb_values": ",".join(str(v) for v in nb_values),
                    "interchromosomal_rows": sum(1 for r in rows if r.interchromosomal),
                    "path": str(path),
                    "command_or_filter": command_for_stage(
                        stage,
                        args.sweepga,
                        native,
                        native_uncompressed,
                        nb1_in,
                        plain,
                        conservative,
                    ),
                }
            )
    return summary_rows, conservative_paths


def command_for_stage(
    stage: str,
    sweepga: Path,
    native: Path,
    native_uncompressed: Path,
    nb1_in: Path,
    plain: Path,
    conservative: Path,
) -> str:
    if stage == "native_n4":
        return "existing native odgi untangle -e 50000 -m 1000 -j 0.8 -n 4 PAF"
    if stage == "native_nb1_prefilter":
        return f"filter {native} to rows whose optional tag list contains nb:i:1"
    if stage == "plain_sweepga_n4_1to1_noscaffold":
        return f"materialize uncompressed native PAF from {native} to {native_uncompressed}; {sweepga} --num-mappings 1:1 --scaffold-jump 0 --output-file {plain} {native_uncompressed}"
    return f"{sweepga} --num-mappings 1:1 --scaffold-jump 0 --output-file {conservative} {nb1_in}"


def coalesce(rows: Iterable[Segment]) -> List[Segment]:
    merged: List[Segment] = []
    for row in sorted(rows, key=lambda r: (r.pair, arm_sort_key(r.qarm), r.qstart, r.target_arm, r.target_hap, r.strand)):
        if (
            merged
            and row.pair == merged[-1].pair
            and row.qname == merged[-1].qname
            and row.target_arm == merged[-1].target_arm
            and row.target_hap == merged[-1].target_hap
            and row.strand == merged[-1].strand
            and row.qstart <= merged[-1].qend + 1
        ):
            prev = merged[-1]
            merged[-1] = Segment(
                prev.pair,
                prev.qname,
                prev.qarm,
                prev.qlen,
                prev.paf_qlen,
                prev.qstart,
                max(prev.qend, row.qend),
                prev.strand,
                prev.tname,
                prev.target_hap,
                prev.target_arm,
                prev.target_sample,
                prev.matches + row.matches,
                prev.block_len + row.block_len,
                max(prev.identity, row.identity),
                max(prev.jc, row.jc),
                prev.nb,
            )
        else:
            merged.append(row)
    return merged


def mapping_style(row: Segment) -> Tuple[str, str, float, str]:
    if row.target_arm == row.qarm:
        return "#cfcfcf", "#9a9a9a", 0.62, "same-arm inherited/background mapping"
    if chrom_of(row.target_arm) == chrom_of(row.qarm):
        return "#b9b1a6", "#8d8378", 0.72, "same-chromosome off-arm mapping"
    return "#b3345b", "#6f1d46", 0.93, "inter-chromosomal candidate mapping"


def submitted_style(row: Segment) -> Tuple[str, str, float, str]:
    if chrom_of(row.target_arm) == chrom_of(row.qarm):
        return "#bfbfbf", "#333333", 0.95, "same chromosome"
    fill = INTERCHR_PALETTE.get(row.target_arm, "#d95f02")
    return fill, "#222222", 0.96, "inter-chromosomal"


def target_lane(row: Segment) -> int:
    if row.target_hap.endswith("#2"):
        return 1
    return 0


def pdf_color(hex_color: str) -> str:
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16) / 255
    g = int(hex_color[2:4], 16) / 255
    b = int(hex_color[4:6], 16) / 255
    return f"{r:.4f} {g:.4f} {b:.4f}"


def panel_geometry() -> Tuple[int, int, int, int, int, int]:
    left = 64
    top = 78
    track_w = 190
    row_h = 9
    pair_gap = 38
    arm_gap = 16
    return left, top, track_w, row_h, pair_gap, arm_gap


@dataclass(frozen=True)
class DrawnRect:
    pair: str
    qarm: str
    x: float
    width: float
    panel_start: float
    panel_end: float


def clip_to_track(row: Segment, tx: float, track_w: float) -> Tuple[float, float]:
    start = max(0, min(row.qstart, row.qlen))
    end = max(0, min(row.qend, row.qlen))
    if end <= start:
        return tx, 0.0
    x = tx + (start / row.qlen) * track_w
    w = ((end - start) / row.qlen) * track_w
    panel_end = tx + track_w
    if x < tx:
        w -= tx - x
        x = tx
    if x + w > panel_end:
        w = panel_end - x
    return x, max(0.0, w)


def validate_segments(outputs: Dict[str, List[Segment]]) -> None:
    bad = [
        row
        for rows in outputs.values()
        for row in rows
        if row.qstart < 0 or row.qend <= row.qstart or row.qend > row.qlen
    ]
    if bad:
        examples = ", ".join(f"{r.pair}:{r.qarm}:{r.qstart}-{r.qend}/{r.qlen}" for r in bad[:5])
        raise RuntimeError(f"invalid conservative segment coordinates after coalescing: {len(bad)} rows; {examples}")


def validate_drawn_rectangles(rectangles: Sequence[DrawnRect], tolerance: float = 1e-6) -> None:
    off_panel = [
        r
        for r in rectangles
        if r.width < -tolerance or r.x < r.panel_start - tolerance or r.x + r.width > r.panel_end + tolerance
    ]
    if off_panel:
        examples = ", ".join(
            f"{r.pair}:{r.qarm}:x={r.x:.3f},w={r.width:.3f},panel={r.panel_start:.3f}-{r.panel_end:.3f}"
            for r in off_panel[:5]
        )
        raise RuntimeError(f"off-panel rectangle geometry: {len(off_panel)} rectangles; {examples}")


def render(
    outputs: Dict[str, List[Segment]],
    phr_intervals: Dict[str, List[PhrInterval]],
    svg_path: Path,
    pdf_path: Path,
) -> List[DrawnRect]:
    left, top, track_w, row_h, pair_gap, arm_gap = panel_geometry()
    panel_w = track_w * 2 + arm_gap
    height = top + len(CHROM_ORDER) * row_h + 100
    width = left + len(PAIR_INFO) * panel_w + (len(PAIR_INFO) - 1) * pair_gap + 44
    title = "Fig. 5-style inspection redraw: conservative native PAF nb:i:1 -> sweepGA 1:1, no scaffold"
    subtitle = "Rows are child/query chromosome ends; p-arm tracks are left and q-arm tracks right. Teal bars mark PHR spans from the WashU all-vs-all p95/id95 table."
    rectangles: List[DrawnRect] = []

    svg: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#222}.title{font-size:16px;font-weight:700}.small{font-size:9px}.label{font-size:8px}.panel{fill:#fbfbfb;stroke:#d7d7d7;stroke-width:.8}.axis{stroke:#a8a8a8;stroke-width:.45}.hit{stroke:#333;stroke-width:.1}.inter{stroke:#5a2ca0;stroke-width:.7}.trackbg{fill:#f2f2f2;stroke:#e1e1e1;stroke-width:.25}",
        "</style>",
        f'<text class="title" x="{left}" y="24">{html.escape(title)}</text>',
        f'<text class="small" x="{left}" y="43">{html.escape(subtitle)}</text>',
    ]
    pdf = PdfWriter(width, height)
    pdf.text(left, 24, title, 13)
    pdf.text(left, 43, subtitle[:160], 7)

    for pidx, (pair, info) in enumerate(PAIR_INFO.items()):
        x0 = left + pidx * (panel_w + pair_gap)
        svg.append(f'<rect class="panel" x="{x0 - 7}" y="{top - 31}" width="{panel_w + 14}" height="{len(CHROM_ORDER) * row_h + 52}" rx="2"/>')
        pdf.rect(x0 - 7, top - 31, panel_w + 14, len(CHROM_ORDER) * row_h + 52, fill="#fbfbfb", stroke="#d7d7d7")
        svg.append(f'<text class="small" x="{x0}" y="{top - 17}">{html.escape(info["label"])}</text>')
        pdf.text(x0, top - 17, info["label"], 7)
        for side, arm_suffix in enumerate(("p", "q")):
            tx = x0 + side * (track_w + arm_gap)
            svg.append(f'<text class="label" x="{tx}" y="{top - 5}">{arm_suffix}-arm 0</text>')
            svg.append(f'<text class="label" x="{tx + track_w - 28}" y="{top - 5}">500 kb</text>')
            pdf.text(tx, top - 5, f"{arm_suffix}-arm 0", 6)
            pdf.text(tx + track_w - 28, top - 5, "500 kb", 6)
        for ridx, chrom in enumerate(CHROM_ORDER):
            y = top + ridx * row_h
            svg.append(f'<text class="label" x="{x0 - 11}" y="{y + 6}" text-anchor="end">{chrom[3:]}</text>')
            pdf.text(x0 - 23, y + 6, chrom[3:], 6)
            for side, arm_suffix in enumerate(("p", "q")):
                tx = x0 + side * (track_w + arm_gap)
                svg.append(f'<rect class="trackbg" x="{tx}" y="{y + 1}" width="{track_w}" height="{row_h - 2}"/>')
                pdf.rect(tx, y + 1, track_w, row_h - 2, fill="#f2f2f2", stroke="#e1e1e1", sw=0.2)
                qname = qname_for_track(outputs[pair], chrom, arm_suffix)
                if qname:
                    for phr in phr_intervals.get(qname, []):
                        start = max(0, min(phr.start, 500000))
                        end = max(0, min(phr.end, 500000))
                        if end <= start:
                            continue
                        px = tx + (start / 500000) * track_w
                        pw = ((end - start) / 500000) * track_w
                        svg.append(
                            f'<rect x="{px:.2f}" y="{y + 1.25:.2f}" width="{pw:.2f}" height="1.55" fill="#287f8f" opacity=".82">'
                            f'<title>PHR span {html.escape(qname)} {start}-{end}; chromosomes {html.escape(phr.chrs_involved)}; arms {html.escape(phr.arms_involved)}</title></rect>'
                        )
                        pdf.rect(px, y + 1.25, pw, 1.55, fill="#287f8f")
        for row in outputs[pair]:
            chrom = chrom_of(row.qarm)
            if chrom not in CHROM_ORDER:
                continue
            side = 0 if row.qarm.endswith("p") else 1
            tx = x0 + side * (track_w + arm_gap)
            y = top + CHROM_ORDER.index(chrom) * row_h + 1.8
            x, w = clip_to_track(row, tx, track_w)
            if w <= 0:
                continue
            fill, stroke, opacity, category = mapping_style(row)
            klass = "hit inter" if row.interchromosomal else "hit"
            rectangles.append(DrawnRect(pair, row.qarm, x, w, tx, tx + track_w))
            svg.append(
                f'<rect class="{klass}" x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{row_h - 3.6:.2f}" fill="{fill}" opacity="{opacity}">'
                f'<title>{html.escape(pair)} {html.escape(row.qarm)} {row.qstart}-{row.qend} / {row.qlen} -> {html.escape(row.target_hap)}:{html.escape(row.target_arm)}; {category}; id={row.identity:.4g} nb={row.nb}</title></rect>'
            )
            pdf.rect(x, y, w, row_h - 3.6, fill=fill, stroke=stroke, sw=0.45 if row.interchromosomal else 0.12)
            if row.interchromosomal and row.length >= 4500 and w >= 2.0:
                label = row.target_arm
                label_x = min(x + w + 1.5, tx + track_w - 18)
                svg.append(f'<text class="label" x="{label_x:.2f}" y="{y + row_h - 4.1:.2f}" fill="#5b1538">{html.escape(label)}</text>')
                pdf.text(label_x, y + row_h - 4.1, label, 5.5, fill="#5b1538")

    legend_y = height - 62
    legend_items = [
        ("#cfcfcf", "#9a9a9a", "same-arm/background"),
        ("#b9b1a6", "#8d8378", "same chromosome, off-arm"),
        ("#b3345b", "#6f1d46", "inter-chromosomal candidate"),
        ("#287f8f", "#287f8f", "PHR span"),
    ]
    lx = left
    for fill, stroke, label in legend_items:
        svg.append(f'<rect x="{lx}" y="{legend_y - 8}" width="16" height="8" fill="{fill}" stroke="{stroke}" stroke-width=".5"/>')
        svg.append(f'<text class="small" x="{lx + 21}" y="{legend_y - 1}">{html.escape(label)}</text>')
        pdf.rect(lx, legend_y - 8, 16, 8, fill=fill, stroke=stroke, sw=0.35)
        pdf.text(lx + 21, legend_y - 1, label, 7)
        lx += 168
    svg.append(f'<text class="small" x="{left}" y="{height - 27}">Inspection output only. This redraw does not replace submission/fig/MainFigures/Fig5_pedigree_untangle.pdf or edit the manuscript.</text>')
    pdf.text(left, height - 27, "Inspection output only. This redraw does not replace the manuscript figure.", 7)
    svg.append("</svg>")
    validate_drawn_rectangles(rectangles)
    svg_path.write_text("\n".join(svg) + "\n")
    pdf.write(pdf_path)
    return rectangles


def qname_for_track(rows: Sequence[Segment], chrom: str, arm_suffix: str) -> str | None:
    qarm = f"{chrom}{arm_suffix}"
    for row in rows:
        if row.qarm == qarm:
            return row.qname
    return None


def submitted_style_render_pair(
    pair: str,
    rows: Sequence[Segment],
    phr_intervals: Dict[str, List[PhrInterval]],
    svg_path: Path,
    pdf_path: Path,
) -> List[DrawnRect]:
    left = 55
    top = 78
    track_w = 560
    track_h = 34
    row_gap = 6
    col_gap = 88
    label_w = 42
    lane_h = 10.5
    lane_gap = 2.5
    phr_h = 4.0
    col_w = label_w + track_w
    width = left + 2 * col_w + col_gap + 40
    height = top + len(PLOT_CHROMS) * (track_h + row_gap) + 58
    rectangles: List[DrawnRect] = []
    by_arm = {(row.qarm, row.qname) for row in rows}
    rows_by_arm: Dict[str, List[Segment]] = {}
    for row in rows:
        rows_by_arm.setdefault(row.qarm, []).append(row)
    title = f"{PAIR_INFO[pair]['label']} - conservative sweepGA 1:1, submitted-style view"
    subtitle = "h1/h2 = parent haplotypes | gray = same chromosome | color = non-homologous chromosome | purple = PHR"

    svg: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#222}.title{font-size:16px;font-weight:700}.small{font-size:9px;fill:#777}.arm{font-size:14px;font-weight:700}.chr{font-size:9px;font-weight:700}.lane{font-size:9px;font-weight:700}.tick{font-size:8px;fill:#666}.label{font-size:8px;font-weight:700}.track{fill:#ffffff;stroke:#222;stroke-width:.8}.nophr{fill:#f7dddd;opacity:.55}.grid{stroke:#e4e4e4;stroke-width:.55}.seg{stroke:#222;stroke-width:.22}.phr{fill:#7e1fa2;opacity:.92}",
        "</style>",
        f'<text class="title" x="{width / 2:.1f}" y="23" text-anchor="middle">{html.escape(title)}</text>',
        f'<text class="small" x="{width / 2:.1f}" y="42" text-anchor="middle">{html.escape(subtitle)}</text>',
    ]
    pdf = PdfWriter(width, height)
    pdf.text(width / 2 - len(title) * 2.2, 23, title, 12)
    pdf.text(width / 2 - len(subtitle) * 1.7, 42, subtitle, 7)

    for col, arm_suffix in enumerate(("p", "q")):
        col_x = left + col * (col_w + col_gap)
        tx0 = col_x + label_w
        svg.append(f'<text class="arm" x="{tx0 + track_w / 2:.1f}" y="{top - 15}" text-anchor="middle">{arm_suffix}-arm</text>')
        pdf.text(tx0 + track_w / 2 - 16, top - 15, f"{arm_suffix}-arm", 10)
        for chrom_idx, chrom in enumerate(PLOT_CHROMS):
            y = top + chrom_idx * (track_h + row_gap)
            qarm = f"{chrom}{arm_suffix}"
            qname = qname_for_track(rows, chrom, arm_suffix)
            has_phr = bool(qname and phr_intervals.get(qname))
            svg.append(f'<rect x="{col_x}" y="{y + 4}" width="{label_w - 7}" height="{track_h - 8}" fill="#e6e6e6" stroke="#777" stroke-width=".7"/>')
            svg.append(f'<text class="chr" x="{col_x + label_w - 12}" y="{y + track_h / 2 + 3}" text-anchor="end">{html.escape(chrom)}</text>')
            pdf.rect(col_x, y + 4, label_w - 7, track_h - 8, fill="#e6e6e6", stroke="#777777", sw=0.5)
            pdf.text(col_x + 5, y + track_h / 2 + 3, chrom, 6)
            svg.append(f'<rect class="track" x="{tx0}" y="{y}" width="{track_w}" height="{track_h}"/>')
            pdf.rect(tx0, y, track_w, track_h, fill="#ffffff", stroke="#222222", sw=0.55)
            if not has_phr:
                svg.append(f'<rect class="nophr" x="{tx0}" y="{y}" width="{track_w}" height="{track_h}"/>')
                pdf.rect(tx0, y, track_w, track_h, fill="#f7dddd")
            for tick in range(0, 501, 100):
                x = tx0 + (tick / 500) * track_w
                svg.append(f'<line class="grid" x1="{x:.2f}" y1="{y}" x2="{x:.2f}" y2="{y + track_h}"/>')
                pdf.line(x, y, x, y + track_h, stroke="#e4e4e4", sw=0.2)
            for lane, label in enumerate(("h1", "h2")):
                ly = y + 4 + lane * (lane_h + lane_gap)
                svg.append(f'<text class="lane" x="{tx0 - 8}" y="{ly + 8}" text-anchor="end">{label}</text>')
                svg.append(f'<line class="grid" x1="{tx0}" y1="{ly + lane_h + 1}" x2="{tx0 + track_w}" y2="{ly + lane_h + 1}"/>')
                pdf.text(tx0 - 16, ly + 8, label, 6)
                pdf.line(tx0, ly + lane_h + 1, tx0 + track_w, ly + lane_h + 1, stroke="#e9e9e9", sw=0.2)
            for phr in phr_intervals.get(qname or "", []):
                start = max(0, min(phr.start, 500000))
                end = max(0, min(phr.end, 500000))
                if end <= start:
                    continue
                px = tx0 + (start / 500000) * track_w
                pw = ((end - start) / 500000) * track_w
                py = y + track_h - phr_h - 2
                svg.append(
                    f'<rect class="phr" x="{px:.2f}" y="{py:.2f}" width="{pw:.2f}" height="{phr_h:.2f}">'
                    f'<title>PHR {html.escape(qname or qarm)} {start}-{end}; arms {html.escape(phr.arms_involved)}</title></rect>'
                )
                pdf.rect(px, py, pw, phr_h, fill="#7e1fa2")

            for row in rows_by_arm.get(qarm, []):
                lane = target_lane(row)
                ly = y + 4 + lane * (lane_h + lane_gap)
                x, w = clip_to_track(row, tx0, track_w)
                if w <= 0:
                    continue
                fill, stroke, opacity, category = submitted_style(row)
                rectangles.append(DrawnRect(pair, row.qarm, x, w, tx0, tx0 + track_w))
                svg.append(
                    f'<rect class="seg" x="{x:.2f}" y="{ly:.2f}" width="{w:.2f}" height="{lane_h:.2f}" fill="{fill}" opacity="{opacity}">'
                    f'<title>{html.escape(row.qarm)} {row.qstart}-{row.qend} -> {html.escape(row.target_hap)}:{html.escape(row.target_arm)}; {category}; id={row.identity:.4g}</title></rect>'
                )
                pdf.rect(x, ly, w, lane_h, fill=fill, stroke=stroke, sw=0.15)
                if row.interchromosomal and (row.length >= 4500 or row.target_arm == "chr15q"):
                    label = row.target_arm[3:]
                    lx = min(max(x + w / 2 - 8, tx0), tx0 + track_w - 20)
                    ly_label = ly - 1
                    svg.append(f'<rect x="{lx - 1:.2f}" y="{ly_label - 8:.2f}" width="{len(label) * 5 + 5}" height="10" fill="#ffffcc" stroke="#444" stroke-width=".4"/>')
                    svg.append(f'<text class="label" x="{lx + 1:.2f}" y="{ly_label:.2f}">{html.escape(label)}</text>')
                    pdf.rect(lx - 1, ly_label - 8, len(label) * 5 + 5, 10, fill="#ffffcc", stroke="#444444", sw=0.2)
                    pdf.text(lx + 1, ly_label, label, 6)

        axis_y = top + len(PLOT_CHROMS) * (track_h + row_gap) - row_gap + 17
        for tick in range(0, 501, 100):
            x = tx0 + (tick / 500) * track_w
            svg.append(f'<text class="tick" x="{x:.2f}" y="{axis_y}" text-anchor="middle">{tick}</text>')
            pdf.text(x - 5, axis_y, str(tick), 6)
        svg.append(f'<text class="tick" x="{tx0 + track_w / 2:.2f}" y="{axis_y + 14}" text-anchor="middle">Position (kb)</text>')
        pdf.text(tx0 + track_w / 2 - 25, axis_y + 14, "Position (kb)", 7)

    legend_y = height - 20
    lx = left + 55
    legend = [
        ("chr10p", INTERCHR_PALETTE["chr10p"]),
        ("chr13p", INTERCHR_PALETTE["chr13p"]),
        ("chr15p", INTERCHR_PALETTE["chr15p"]),
        ("chr15q", INTERCHR_PALETTE["chr15q"]),
        ("chr19p", INTERCHR_PALETTE["chr19p"]),
        ("chr21p", INTERCHR_PALETTE["chr21p"]),
        ("chr22p", INTERCHR_PALETTE["chr22p"]),
        ("chr9q", INTERCHR_PALETTE["chr9q"]),
        ("same chr", "#bfbfbf"),
        ("PHR", "#7e1fa2"),
    ]
    for label, color in legend:
        svg.append(f'<rect x="{lx}" y="{legend_y - 9}" width="52" height="10" fill="{color}" stroke="#222" stroke-width=".45"/>')
        svg.append(f'<text class="label" x="{lx + 26}" y="{legend_y - 12}" text-anchor="middle">{html.escape(label[3:] if label.startswith("chr") else label)}</text>')
        pdf.rect(lx, legend_y - 9, 52, 10, fill=color, stroke="#222222", sw=0.2)
        pdf.text(lx + 16, legend_y - 12, label[3:] if label.startswith("chr") else label, 6)
        lx += 78
    svg.append("</svg>")
    validate_drawn_rectangles(rectangles)
    svg_path.write_text("\n".join(svg) + "\n")
    pdf.write(pdf_path)
    return rectangles


def render_submitted_style_outputs(
    outputs: Dict[str, List[Segment]],
    phr_intervals: Dict[str, List[PhrInterval]],
    out_dir: Path,
) -> List[DrawnRect]:
    rectangles: List[DrawnRect] = []
    for pair, rows in outputs.items():
        rectangles.extend(
            submitted_style_render_pair(
                pair,
                rows,
                phr_intervals,
                out_dir / f"fig5_sweepga_1to1_submitted_style_{pair}.svg",
                out_dir / f"fig5_sweepga_1to1_submitted_style_{pair}.pdf",
            )
        )
    return rectangles


def write_tsv(path: Path, rows: Sequence[dict], fields: Sequence[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_readme(out_dir: Path, summary_rows: Sequence[dict], phr_interval_path: Path) -> None:
    counts = {(r["pair"], r["stage"]): r for r in summary_rows}
    lines = [
        "# Fig5 sweepGA 1:1 inspection redraw",
        "",
        "This directory contains an author-facing inspection redraw of the pedigree/Fig. 5-style untangle panel.",
        "It is exploratory output only and does not update `submission/paper.tex` or `submission/fig/MainFigures/Fig5_pedigree_untangle.pdf`.",
        "",
        "## Inputs and filtering",
        "",
        "For each WashU transmission, the script starts from the existing native odgi untangle m1000 n4 PAF in `/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/`.",
        "The conservative plotted view first filters the native PAF to rows whose optional tags contain `nb:i:1`, then runs:",
        "",
        "```bash",
        "/moosefs/erikg/sweepga/target/release/sweepga --num-mappings 1:1 --scaffold-jump 0 --output-file <conservative.paf> <native_nb1.paf>",
        "```",
        "",
        "This differs from plain sweepGA 1:1 because plain sweepGA is run directly on the native n4 PAF and can retain `nb:i:2`, `nb:i:3`, or `nb:i:4` alternates when those are equivalent under its filter. The conservative redraw removes those alternates before sweepGA, so the plotted PAF contains only `nb:i:1` rows.",
        "",
        "It also differs from native first-best alone because sweepGA applies an additional reciprocal 1:1 plane-sweep filter with scaffolding disabled (`--scaffold-jump 0`).",
        "",
        "## Row counts",
        "",
        "| Pair | Native n4 | Native nb1 | Plain sweepGA n4 1:1 no scaffold | Conservative nb1 -> sweepGA 1:1 no scaffold | Plain nb values | Conservative nb values |",
        "| --- | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for pair in PAIR_INFO:
        native = counts[(pair, "native_n4")]
        nb1 = counts[(pair, "native_nb1_prefilter")]
        plain = counts[(pair, "plain_sweepga_n4_1to1_noscaffold")]
        cons = counts[(pair, "conservative_nb1_sweepga_1to1_noscaffold")]
        lines.append(
            f"| `{pair}` | {native['rows']} | {nb1['rows']} | {plain['rows']} | {cons['rows']} | {plain['nb_values']} | {cons['nb_values']} |"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `fig5_sweepga_1to1_redraw.svg` and `fig5_sweepga_1to1_redraw.pdf`: compact author-facing redraw with teal PHR-span overlays.",
            "- `fig5_sweepga_1to1_submitted_style_<PAIR>.svg` and `.pdf`: sweepGA-filtered data rendered in the same p-arm/q-arm, h1/h2-track format as the raw `odgi untangle` figure.",
            "- `summary_counts.tsv`: row counts, `nb` values, inter-chromosomal row counts, source paths, and exact command/filter per stage.",
            "- `conservative_segments.tsv`: compact plotted segment table after coalescing adjacent same-target rows.",
            "- `phr_intervals.tsv`: PHR spans plotted on the child/query 500 kb tracks.",
        "- `validation_report.tsv`: coordinate, `nb`, query-length, and SVG/PDF rectangle-bound checks from the last regeneration.",
        "",
        "## PHR overlay",
        "",
        f"PHR spans are read from `{phr_interval_path}`. Coordinates are already in the 0-500 kb flank coordinate system used by the child/query tracks. Rows with `.` region coordinates are treated as no plotted PHR span.",
        "",
        "## Coordinate correction",
        "",
            "The superseded first redraw incorrectly used PAF column 2 (`fields[1]`) as the full plotting denominator. In these odgi PAFs that field can be 41,248 even when `qstart`/`qend` span hundreds of kilobases, which created invalid coalesced rows and off-panel SVG rectangles.",
            "",
            "The corrected redraw derives `query_length` from the child/query path name interval, for example `PAN027#1#chr3.maternal:9503-509502_chr3_parm` gives 500,000 bp. The raw PAF column 2 value is retained only as `paf_query_length` for audit. Rows with invalid child/query coordinates are dropped before coalescing, and rendering validates that all rectangles stay within their p/q-arm track bounds.",
            "",
            "Regenerate from the repository root with:",
            "",
            "```bash",
            "python3 scripts/pedigree/plot_fig5_sweepga_1to1_redraw.py",
            "```",
        ]
    )
    (out_dir / "README.md").write_text("\n".join(lines) + "\n")


def validation_report(
    outputs: Dict[str, List[Segment]],
    rectangles: Sequence[DrawnRect],
    phr_intervals: Dict[str, List[PhrInterval]],
    submitted_style_rectangles: Sequence[DrawnRect],
) -> List[dict]:
    rows = [row for pair_rows in outputs.values() for row in pair_rows]
    query_names = {row.qname for row in rows}
    invalid_coordinates = sum(1 for row in rows if row.qstart < 0 or row.qend <= row.qstart or row.qend > row.qlen)
    invalid_phrs = sum(
        1
        for intervals in phr_intervals.values()
        for interval in intervals
        if interval.start < 0 or interval.end <= interval.start or interval.end > 500000
    )
    off_panel = sum(
        1
        for rect in rectangles
        if rect.width < -1e-6 or rect.x < rect.panel_start - 1e-6 or rect.x + rect.width > rect.panel_end + 1e-6
    )
    query_lengths = ",".join(str(v) for v in sorted({row.qlen for row in rows}))
    paf_query_lengths = ",".join(str(v) for v in sorted({row.paf_qlen for row in rows}))
    nb_values = ",".join(str(v) for v in sorted({row.nb for row in rows}))
    return [
        {"check": "plotted_segments", "value": str(len(rows)), "status": "PASS"},
        {"check": "invalid_coordinate_rows", "value": str(invalid_coordinates), "status": "PASS" if invalid_coordinates == 0 else "FAIL"},
        {"check": "query_lengths", "value": query_lengths, "status": "PASS" if query_lengths == "500000" else "FAIL"},
        {"check": "raw_paf_query_lengths_audit", "value": paf_query_lengths, "status": "PASS"},
        {"check": "nb_values", "value": nb_values, "status": "PASS" if nb_values == "1" else "FAIL"},
        {"check": "query_names_with_phr_rows", "value": str(sum(1 for name in query_names if phr_intervals.get(name))), "status": "PASS"},
        {"check": "plotted_phr_intervals", "value": str(sum(len(v) for v in phr_intervals.values())), "status": "PASS"},
        {"check": "invalid_phr_intervals", "value": str(invalid_phrs), "status": "PASS" if invalid_phrs == 0 else "FAIL"},
        {"check": "drawn_rectangles", "value": str(len(rectangles)), "status": "PASS"},
        {"check": "off_panel_rectangles", "value": str(off_panel), "status": "PASS" if off_panel == 0 else "FAIL"},
        {"check": "submitted_style_rectangles", "value": str(len(submitted_style_rectangles)), "status": "PASS"},
    ]


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary_rows, conservative_paths = ensure_outputs(args)

    outputs: Dict[str, List[Segment]] = {}
    compact_rows: List[dict] = []
    query_names: set[str] = set()
    for pair, path in conservative_paths.items():
        rows = read_paf(path, pair, require_nb1=True)
        if any(r.nb != 1 for r in rows):
            raise RuntimeError(f"non-nb1 row survived in conservative plotted PAF for {pair}")
        compact = coalesce(rows)
        outputs[pair] = compact
        for row in compact:
            query_names.add(row.qname)
            compact_rows.append(
                {
                    "pair": pair,
                    "transmission": PAIR_INFO[pair]["label"],
                    "query_name": row.qname,
                    "query_arm": row.qarm,
                    "query_start": row.qstart,
                    "query_end": row.qend,
                    "query_length": row.qlen,
                    "paf_query_length": row.paf_qlen,
                    "target_name": row.tname,
                    "target_hap": row.target_hap,
                    "target_arm": row.target_arm,
                    "strand": row.strand,
                    "identity": f"{row.identity:.6g}",
                    "jaccard": f"{row.jc:.6g}",
                    "nb": row.nb,
                    "interchromosomal": int(row.interchromosomal),
                }
            )
    phr_intervals = load_phr_intervals(args.phr_intervals, query_names)
    phr_rows = []
    for qname, intervals in sorted(phr_intervals.items()):
        info = path_info(qname)
        qarm = info[2] if info else ""
        for interval in intervals:
            phr_rows.append(
                {
                    "query_name": qname,
                    "query_arm": qarm,
                    "phr_start": interval.start,
                    "phr_end": interval.end,
                    "chrs_involved": interval.chrs_involved,
                    "arms_involved": interval.arms_involved,
                }
            )
    write_tsv(
        args.out_dir / "summary_counts.tsv",
        summary_rows,
        ["pair", "transmission", "stage", "rows", "nb_values", "interchromosomal_rows", "path", "command_or_filter"],
    )
    write_tsv(
        args.out_dir / "conservative_segments.tsv",
        compact_rows,
        [
            "pair",
            "transmission",
            "query_name",
            "query_arm",
            "query_start",
            "query_end",
            "query_length",
            "paf_query_length",
            "target_name",
            "target_hap",
            "target_arm",
            "strand",
            "identity",
            "jaccard",
            "nb",
            "interchromosomal",
        ],
    )
    write_tsv(
        args.out_dir / "phr_intervals.tsv",
        phr_rows,
        ["query_name", "query_arm", "phr_start", "phr_end", "chrs_involved", "arms_involved"],
    )
    validate_segments(outputs)
    rectangles = render(
        outputs,
        phr_intervals,
        args.out_dir / "fig5_sweepga_1to1_redraw.svg",
        args.out_dir / "fig5_sweepga_1to1_redraw.pdf",
    )
    validate_drawn_rectangles(rectangles)
    submitted_style_rectangles = render_submitted_style_outputs(outputs, phr_intervals, args.out_dir)
    validate_drawn_rectangles(submitted_style_rectangles)
    write_tsv(
        args.out_dir / "validation_report.tsv",
        validation_report(outputs, rectangles, phr_intervals, submitted_style_rectangles),
        ["check", "value", "status"],
    )
    write_readme(args.out_dir, summary_rows, args.phr_intervals)


if __name__ == "__main__":
    main()
