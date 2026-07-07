#!/usr/bin/env python3
"""Compose manuscript Fig. 5 without embedding panel PNGs through librsvg."""

from __future__ import annotations

import os
import shutil
import struct
import subprocess
import zlib
from binascii import crc32
from html import escape
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
MAIN_DIR = ROOT / "paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft"
MAT_DIR = ROOT / "paper_prep/_brainstorming/fig5_extended_maternal_whole_genome_ribbons"

PANEL_W = 7200
PANEL_CROP_H = 1300
PEDIGREE_H = 610
PANEL_GAP = 40
LEGEND_H = 320

TEXT = "#202124"
MUTED = "#5f6368"
GRID = "#d9dde1"
HOMOLOG = "#b8bdc3"
ACRO = "#7f7f7f"
PAR1 = "#E7298A"
CHR5_CHR1 = "#2C7FB8"
CHR9_CHR3 = "#E26D0A"
OTHER = "#969696"


def text(x: float, y: float, value: str, size: int, weight: int = 400, fill: str = TEXT,
         anchor: str = "start") -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">{escape(value)}</text>'
    )


def line(x1: float, y1: float, x2: float, y2: float, stroke: str = TEXT, width: float = 5.0,
         opacity: float = 1.0) -> str:
    return (
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
        f'stroke="{stroke}" stroke-width="{width:.1f}" opacity="{opacity:.2f}"/>'
    )


def svg_doc(width: int, height: int, body: str) -> str:
    return "\n".join([
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,Helvetica,sans-serif;dominant-baseline:alphabetic}</style>',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
        body,
        "</svg>",
        "",
    ])


def pedigree_svg() -> str:
    parts: list[str] = []
    parts.append(text(80, 132, "A", 118, 700))
    parts.append(text(270, 96, "Pedigree comparisons", 78, 700))

    father_x, mother_x = 1180, 2350
    child_x, pan028_x = 1765, 3600
    parent_y, child_y = 250, 455
    r = 72

    parts.extend([
        f'<rect x="{father_x - r}" y="{parent_y - r}" width="{2 * r}" height="{2 * r}" '
        f'fill="#ffffff" stroke="{TEXT}" stroke-width="9"/>',
        f'<circle cx="{mother_x}" cy="{parent_y}" r="{r}" fill="#ffffff" stroke="{TEXT}" stroke-width="9"/>',
        line(father_x + r, parent_y, mother_x - r, parent_y, TEXT, 6),
        line((father_x + mother_x) / 2, parent_y, (father_x + mother_x) / 2, child_y - r, TEXT, 6),
        f'<circle cx="{child_x}" cy="{child_y}" r="{r}" fill="#ffffff" stroke="{TEXT}" stroke-width="9"/>',
        line((father_x + mother_x) / 2, child_y - r, child_x, child_y - r, TEXT, 6),
        f'<circle cx="{pan028_x}" cy="{child_y}" r="{r}" fill="#ffffff" stroke="{TEXT}" stroke-width="9"/>',
        line(child_x + r, child_y, pan028_x - r, child_y, TEXT, 6),
        text(father_x, parent_y + 158, "PAN011 father", 48, 700, TEXT, "middle"),
        text(mother_x, parent_y + 158, "PAN010 mother", 48, 700, TEXT, "middle"),
        text(child_x, child_y + 158, "PAN027", 48, 700, TEXT, "middle"),
        text(pan028_x, child_y + 158, "PAN028", 48, 700, TEXT, "middle"),
    ])

    x0 = 4550
    parts.extend([
        text(x0, 188, "B  PAN027 paternal haplotype vs PAN011 father", 52, 700, TEXT),
        text(x0, 305, "C  PAN027 maternal haplotype vs PAN010 mother", 52, 700, TEXT),
        text(x0, 422, "D  PAN028 maternal haplotype vs PAN027 mother", 52, 700, TEXT),
        f'<rect x="{x0 - 72}" y="145" width="34" height="34" fill="{CHR9_CHR3}" opacity="0.95"/>',
        f'<rect x="{x0 - 72}" y="262" width="34" height="34" fill="{ACRO}" opacity="0.95"/>',
        f'<rect x="{x0 - 72}" y="379" width="34" height="34" fill="{CHR5_CHR1}" opacity="0.95"/>',
        line(0, PEDIGREE_H - 22, PANEL_W, PEDIGREE_H - 22, GRID, 3, 1),
    ])
    return svg_doc(PANEL_W, PEDIGREE_H, "\n".join(parts))


def legend_svg() -> str:
    entries = [
        (HOMOLOG, "same chromosome"),
        (ACRO, "acrocentric"),
        (PAR1, "PAR1"),
        (CHR5_CHR1, "chr5q/chr1p"),
        (CHR9_CHR3, "chr9q/chr3q"),
        (OTHER, "other non-homologous"),
    ]
    parts = [text(620, 100, "Ribbon classes", 62, 700, MUTED)]
    x = 1500
    for idx, (color, label) in enumerate(entries):
        row = idx // 3
        col = idx % 3
        xx = x + col * 1500
        yy = 66 + row * 126
        parts.append(f'<rect x="{xx}" y="{yy}" width="70" height="54" fill="{color}" opacity="0.9"/>')
        parts.append(text(xx + 96, yy + 47, label, 54, 400, TEXT))
    return svg_doc(PANEL_W, LEGEND_H, "\n".join(parts))


def find_rsvg() -> str:
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
    raise RuntimeError("rsvg-convert not found")


def render_strip(svg: str, stem: str) -> Path:
    svg_path = MAIN_DIR / f"{stem}.svg"
    png_path = MAIN_DIR / f"{stem}.png"
    svg_path.write_text(svg)
    subprocess.run([find_rsvg(), "-f", "png", "-o", str(png_path), str(svg_path)], check=True)
    return png_path


def paeth(a: int, b: int, c: int) -> int:
    p = a + b - c
    pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
    if pa <= pb and pa <= pc:
        return a
    if pb <= pc:
        return b
    return c


def read_png_rgb(path: Path) -> tuple[int, int, list[bytes]]:
    data = path.read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"not a PNG: {path}")
    pos = 8
    width = height = color_type = None
    idat: list[bytes] = []
    while pos < len(data):
        n = struct.unpack(">I", data[pos:pos + 4])[0]
        pos += 4
        chunk_type = data[pos:pos + 4]
        pos += 4
        chunk = data[pos:pos + n]
        pos += n + 4
        if chunk_type == b"IHDR":
            width, height, bit_depth, color_type, comp, filt, interlace = struct.unpack(">IIBBBBB", chunk)
            if bit_depth != 8 or comp != 0 or filt != 0 or interlace != 0:
                raise ValueError(f"unsupported PNG format: {path}")
            if color_type not in {2, 6}:
                raise ValueError(f"unsupported PNG color type {color_type}: {path}")
        elif chunk_type == b"IDAT":
            idat.append(chunk)
        elif chunk_type == b"IEND":
            break
    if width is None or height is None or color_type is None:
        raise ValueError(f"missing PNG header: {path}")

    bpp = 3 if color_type == 2 else 4
    stride = width * bpp
    raw = zlib.decompress(b"".join(idat))
    rows: list[bytes] = []
    prev = bytearray(stride)
    i = 0
    for _ in range(height):
        filter_type = raw[i]
        i += 1
        row = bytearray(raw[i:i + stride])
        i += stride
        if filter_type == 1:
            for x in range(bpp, stride):
                row[x] = (row[x] + row[x - bpp]) & 255
        elif filter_type == 2:
            for x in range(stride):
                row[x] = (row[x] + prev[x]) & 255
        elif filter_type == 3:
            for x in range(stride):
                left = row[x - bpp] if x >= bpp else 0
                row[x] = (row[x] + ((left + prev[x]) // 2)) & 255
        elif filter_type == 4:
            for x in range(stride):
                left = row[x - bpp] if x >= bpp else 0
                up = prev[x]
                up_left = prev[x - bpp] if x >= bpp else 0
                row[x] = (row[x] + paeth(left, up, up_left)) & 255
        elif filter_type != 0:
            raise ValueError(f"unsupported PNG filter {filter_type}: {path}")

        if color_type == 2:
            rows.append(bytes(row))
        else:
            rgb = bytearray(width * 3)
            for px in range(width):
                r, g, b, a = row[4 * px:4 * px + 4]
                rgb[3 * px] = (r * a + 255 * (255 - a)) // 255
                rgb[3 * px + 1] = (g * a + 255 * (255 - a)) // 255
                rgb[3 * px + 2] = (b * a + 255 * (255 - a)) // 255
            rows.append(bytes(rgb))
        prev = row
    return width, height, rows


def png_chunk(kind: bytes, payload: bytes) -> bytes:
    return (
        struct.pack(">I", len(payload))
        + kind
        + payload
        + struct.pack(">I", crc32(kind + payload) & 0xFFFFFFFF)
    )


def write_png_rgb(path: Path, width: int, rows: list[bytes]) -> None:
    height = len(rows)
    raw = b"".join(b"\x00" + row for row in rows)
    payload = [
        b"\x89PNG\r\n\x1a\n",
        png_chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)),
        png_chunk(b"IDAT", zlib.compress(raw, level=6)),
        png_chunk(b"IEND", b""),
    ]
    path.write_bytes(b"".join(payload))


def checked_rows(path: Path, crop: int | None = None) -> list[bytes]:
    width, height, rows = read_png_rgb(path)
    if width != PANEL_W:
        raise ValueError(f"unexpected width in {path}: {width}")
    if crop is not None:
        if height < crop:
            raise ValueError(f"cannot crop {crop} rows from {path} with height {height}")
        return rows[:crop]
    return rows


def main() -> None:
    panel_paths = [
        MAIN_DIR / "fig5_homologous_recombination_context_ribbon_draft.png",
        MAT_DIR / "PAN027mat_vs_PAN010_joint/PAN027mat_vs_PAN010_joint.whole_genome_homologous_context_ribbon.png",
        MAT_DIR / "PAN028mat_vs_PAN027_joint/PAN028mat_vs_PAN027_joint.whole_genome_homologous_context_ribbon.png",
    ]
    for panel_path in panel_paths:
        if not panel_path.exists():
            raise FileNotFoundError(panel_path)

    pedigree_png = render_strip(pedigree_svg(), "fig5_pedigree_key")
    legend_png = render_strip(legend_svg(), "fig5_ribbon_legend")
    gap = bytes([255]) * PANEL_W * 3

    rows: list[bytes] = []
    rows.extend(checked_rows(pedigree_png))
    rows.extend([gap] * PANEL_GAP)
    for idx, panel_path in enumerate(panel_paths):
        if idx:
            rows.extend([gap] * PANEL_GAP)
        rows.extend(checked_rows(panel_path, PANEL_CROP_H))
    rows.extend(checked_rows(legend_png))

    out_png = MAIN_DIR / "fig5_pedigree_recombination_combined.png"
    write_png_rgb(out_png, PANEL_W, rows)

    for stale in [
        MAIN_DIR / "fig5_pedigree_recombination_combined.svg",
        MAIN_DIR / "fig5_pedigree_recombination_combined.pdf",
    ]:
        if stale.exists():
            stale.unlink()


if __name__ == "__main__":
    main()
