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
PATH_RE = re.compile(r"^(PAN\d+)#([12])#.*_(chr(?:[0-9]+|X|Y))_([pq])arm$")
TAG_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]):[A-Za-z]:(.*)$")


@dataclass(frozen=True)
class Segment:
    pair: str
    qname: str
    qarm: str
    qlen: int
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
    parser.add_argument("--force", action="store_true", help="Regenerate external sweepGA PAF outputs even if present")
    parser.add_argument("--skip-sweepga", action="store_true", help="Use existing external PAF outputs only")
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def path_info(name: str) -> Tuple[str, str, str] | None:
    match = PATH_RE.match(name)
    if not match:
        return None
    sample, hap, chrom, arm = match.groups()
    return sample, hap, f"{chrom}{arm}"


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
            qsample, qhap, qarm = qinfo
            tsample, thap, tarm = tinfo
            if qsample != info["query_sample"] or qhap != info["query_hap"] or tsample != info["target_sample"]:
                continue
            tags = parse_tags(fields[12:])
            nb = int(tags.get("nb", "1"))
            if require_nb1 and nb != 1:
                continue
            rows.append(
                Segment(
                    pair=pair,
                    qname=fields[0],
                    qarm=qarm,
                    qlen=int(fields[1]),
                    qstart=max(0, int(fields[2])),
                    qend=min(int(fields[1]), int(fields[3])),
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


def target_color(target_arm: str) -> str:
    c, side = arm_sort_key(target_arm)
    hue = (c * 47 + side * 21) % 360
    return hsl_to_hex(hue, 0.58 if side == 0 else 0.66, 0.42 if side == 0 else 0.52)


def hsl_to_hex(h: float, s: float, l: float) -> str:
    c = (1 - abs(2 * l - 1)) * s
    hp = h / 60
    x = c * (1 - abs(hp % 2 - 1))
    if hp < 1:
        r, g, b = c, x, 0
    elif hp < 2:
        r, g, b = x, c, 0
    elif hp < 3:
        r, g, b = 0, c, x
    elif hp < 4:
        r, g, b = 0, x, c
    elif hp < 5:
        r, g, b = x, 0, c
    else:
        r, g, b = c, 0, x
    m = l - c / 2
    return "#{:02x}{:02x}{:02x}".format(round((r + m) * 255), round((g + m) * 255), round((b + m) * 255))


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


def render(outputs: Dict[str, List[Segment]], svg_path: Path, pdf_path: Path) -> None:
    left, top, track_w, row_h, pair_gap, arm_gap = panel_geometry()
    panel_w = track_w * 2 + arm_gap
    height = top + len(CHROM_ORDER) * row_h + 100
    width = left + len(PAIR_INFO) * panel_w + (len(PAIR_INFO) - 1) * pair_gap + 44
    title = "Fig. 5-style inspection redraw: conservative native PAF nb:i:1 -> sweepGA 1:1, no scaffold"
    subtitle = "Rows are child/query chromosome ends; p-arm tracks are left and q-arm tracks right within each pedigree pair. Fill color encodes target/parent arm; purple outline marks inter-chromosomal mappings."

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
        for row in outputs[pair]:
            chrom = chrom_of(row.qarm)
            if chrom not in CHROM_ORDER:
                continue
            side = 0 if row.qarm.endswith("p") else 1
            tx = x0 + side * (track_w + arm_gap)
            y = top + CHROM_ORDER.index(chrom) * row_h + 1.8
            qdenom = max(row.qlen, row.qend, 1)
            x = tx + (row.qstart / qdenom) * track_w
            w = max(0.45, (row.length / qdenom) * track_w)
            fill = target_color(row.target_arm)
            stroke = "#5a2ca0" if row.interchromosomal else "#333333"
            klass = "hit inter" if row.interchromosomal else "hit"
            opacity = "0.92" if row.interchromosomal else "0.42"
            svg.append(
                f'<rect class="{klass}" x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{row_h - 3.6:.2f}" fill="{fill}" opacity="{opacity}">'
                f'<title>{html.escape(pair)} {html.escape(row.qarm)} {row.qstart}-{row.qend} -> {html.escape(row.target_hap)}:{html.escape(row.target_arm)} id={row.identity:.4g} nb={row.nb}</title></rect>'
            )
            pdf.rect(x, y, w, row_h - 3.6, fill=fill, stroke=stroke, sw=0.45 if row.interchromosomal else 0.12)

    legend_y = height - 62
    svg.append(f'<text class="small" x="{left}" y="{legend_y}">Target-arm color examples</text>')
    pdf.text(left, legend_y, "Target-arm color examples", 7)
    lx = left + 126
    for arm in ["chr3p", "chr4q", "chr10p", "chr13p", "chr14p", "chr15p", "chr21p", "chr22p", "chrXq"]:
        svg.append(f'<rect x="{lx}" y="{legend_y - 8}" width="12" height="8" fill="{target_color(arm)}"/>')
        svg.append(f'<text class="label" x="{lx + 15}" y="{legend_y - 1}">{arm}</text>')
        pdf.rect(lx, legend_y - 8, 12, 8, fill=target_color(arm))
        pdf.text(lx + 15, legend_y - 1, arm, 6)
        lx += 52
    svg.append(f'<text class="small" x="{left}" y="{height - 27}">Inspection output only. This redraw does not replace submission/fig/MainFigures/Fig5_pedigree_untangle.pdf or edit the manuscript.</text>')
    pdf.text(left, height - 27, "Inspection output only. This redraw does not replace the manuscript figure.", 7)
    svg.append("</svg>")
    svg_path.write_text("\n".join(svg) + "\n")
    pdf.write(pdf_path)


def write_tsv(path: Path, rows: Sequence[dict], fields: Sequence[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_readme(out_dir: Path, summary_rows: Sequence[dict]) -> None:
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
            "- `fig5_sweepga_1to1_redraw.svg` and `fig5_sweepga_1to1_redraw.pdf`: compact author-facing redraw.",
            "- `summary_counts.tsv`: row counts, `nb` values, inter-chromosomal row counts, source paths, and exact command/filter per stage.",
            "- `conservative_segments.tsv`: compact plotted segment table after coalescing adjacent same-target rows.",
            "",
            "Regenerate from the repository root with:",
            "",
            "```bash",
            "python3 scripts/pedigree/plot_fig5_sweepga_1to1_redraw.py",
            "```",
        ]
    )
    (out_dir / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary_rows, conservative_paths = ensure_outputs(args)

    outputs: Dict[str, List[Segment]] = {}
    compact_rows: List[dict] = []
    for pair, path in conservative_paths.items():
        rows = read_paf(path, pair, require_nb1=True)
        if any(r.nb != 1 for r in rows):
            raise RuntimeError(f"non-nb1 row survived in conservative plotted PAF for {pair}")
        compact = coalesce(rows)
        outputs[pair] = compact
        for row in compact:
            compact_rows.append(
                {
                    "pair": pair,
                    "transmission": PAIR_INFO[pair]["label"],
                    "query_name": row.qname,
                    "query_arm": row.qarm,
                    "query_start": row.qstart,
                    "query_end": row.qend,
                    "query_length": row.qlen,
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
    render(
        outputs,
        args.out_dir / "fig5_sweepga_1to1_redraw.svg",
        args.out_dir / "fig5_sweepga_1to1_redraw.pdf",
    )
    write_readme(args.out_dir, summary_rows)


if __name__ == "__main__":
    main()
