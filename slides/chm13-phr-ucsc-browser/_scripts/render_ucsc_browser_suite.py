#!/usr/bin/env python3
"""Render terminal-anchored CHM13 PHR UCSC browser panels and a Typst deck.

The script intentionally uses UCSC's real hgTracks image output.  For each BED
row it fetches hgTracks HTML, verifies the chm13.phrs.bed custom track, downloads
the emitted trash sprites, crops the visible side-label and data panels by the
CSS offsets, appends them, and writes a manifest plus Typst source/PDF.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from urllib.parse import urlencode, urljoin


DB = "hub_3671779_hs1"
HGSID = "3966979908_lGks6rs34CqmdawD8iwY2YCYKVd3"
PIX = 1800
UCSC_BASE = "https://genome.ucsc.edu/cgi-bin/hgTracks"
CHROM_SIZES_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes"
)
USER_AGENT = "Mozilla/5.0"

REPO_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = REPO_ROOT / "slides" / "chm13-phr-ucsc-browser"
ASSET_DIR = OUT_DIR / "_assets" / "ucsc"
HTML_DIR = ASSET_DIR / "html"
SPRITE_DIR = ASSET_DIR / "sprites"
PANEL_DIR = ASSET_DIR / "panels"
CHROM_SIZES_PATH = ASSET_DIR / "hs1.chrom.sizes"
MANIFEST_PATH = OUT_DIR / "manifest.tsv"
AUDIT_PATH = OUT_DIR / "audit_37_vs_41.tsv"
TYP_PATH = OUT_DIR / "chm13_phr_ucsc_browser_suite.typ"
PDF_PATH = OUT_DIR / "CHM13_PHR_UCSC_browser_suite.pdf"

BED_PATH = REPO_ROOT / "chm13.phrs.bed"
ARCH_PATH = REPO_ROOT / "paper_prep" / "figures" / "fig1" / "architecture_per_arm.tsv"


@dataclass(frozen=True)
class Slice:
    src: str
    x: int
    y: int
    width: int
    height: int


@dataclass(frozen=True)
class BrowserRow:
    row_id: str
    side: Slice
    data: Slice


@dataclass
class PhrPanel:
    bed_index: int
    label: str
    chrom_arm: str
    chrom: str
    arm: str
    bed_start0: int
    bed_end0: int
    bed_field4: str
    phr_inclusive_bp: int
    chrom_size_bp: int
    browser_start1: int
    browser_end1: int
    browser_window_bp: int
    terminal_gap_bp: int
    p_terminal_gap_bp: int
    q_terminal_gap_bp: int
    ucsc_url: str
    html_path: Path
    image_path: Path
    html_track_confirmed: bool = False
    ucsc_position_display: str = ""
    image_width_px: int = 0
    image_height_px: int = 0
    image_stdev_gray: float = 0.0


def run(cmd: list[str], *, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    printable = " ".join(cmd)
    print(f"[run] {printable}", flush=True)
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def tool_path(name: str, fallback: str | None = None) -> str:
    if fallback and Path(fallback).exists():
        return fallback
    found = shutil.which(name)
    if found:
        return found
    raise SystemExit(f"Required tool not found on PATH: {name}")


def fetch_url(url: str, dest: Path, force: bool) -> None:
    if dest.exists() and not force:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    cmd = [
        "curl",
        "-fsSL",
        "-A",
        USER_AGENT,
        "--compressed",
        "--retry",
        "3",
        "--retry-delay",
        "2",
        "--connect-timeout",
        "20",
        "--max-time",
        "180",
        url,
        "-o",
        str(tmp),
    ]
    run(cmd)
    tmp.replace(dest)
    if dest.suffix.lower() in {".html", ".htm"}:
        strip_trailing_whitespace(dest)


def strip_trailing_whitespace(path: Path) -> None:
    text = path.read_text(errors="replace")
    cleaned = "\n".join(line.rstrip(" \t") for line in text.splitlines())
    if text.endswith("\n"):
        cleaned += "\n"
    path.write_text(cleaned)


def read_chrom_sizes(force: bool) -> dict[str, int]:
    fetch_url(CHROM_SIZES_URL, CHROM_SIZES_PATH, force=force)
    sizes: dict[str, int] = {}
    with CHROM_SIZES_PATH.open() as fh:
        for line in fh:
            chrom, size = line.rstrip("\n").split("\t")[:2]
            sizes[chrom] = int(size)
    return sizes


def ucsc_url(chrom: str, start1: int, end1: int) -> str:
    query = urlencode(
        {
            "db": DB,
            "hgsid": HGSID,
            "position": f"{chrom}:{start1}-{end1}",
            "pix": str(PIX),
        }
    )
    return f"{UCSC_BASE}?{query}"


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def read_bed_panels(sizes: dict[str, int]) -> list[PhrPanel]:
    panels: list[PhrPanel] = []
    seen: dict[str, int] = {}
    with BED_PATH.open() as fh:
        for idx, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"BED row {idx} has fewer than 3 columns")
            chrom = fields[0]
            if chrom not in sizes:
                raise KeyError(f"Missing hs1 chromosome size for {chrom}")
            bed_start0 = int(fields[1])
            bed_end0 = int(fields[2])
            bed_field4 = fields[3] if len(fields) > 3 else ""
            chrom_size = sizes[chrom]
            phr_len = bed_end0 - bed_start0 + 1
            if phr_len <= 0:
                raise ValueError(f"BED row {idx} has non-positive inclusive length: {line}")
            view_len = (3 * phr_len + 1) // 2
            p_gap = bed_start0
            q_gap = chrom_size - bed_end0
            if p_gap <= q_gap:
                arm = "p"
                start1 = 1
                end1 = min(chrom_size, view_len)
                terminal_gap = p_gap
            else:
                arm = "q"
                start1 = max(1, chrom_size - view_len + 1)
                end1 = chrom_size
                terminal_gap = q_gap
            chrom_arm = f"{chrom}_{arm}"
            base_label = f"{chrom}{arm}"
            seen[base_label] = seen.get(base_label, 0) + 1
            label = base_label if seen[base_label] == 1 else f"{base_label}_{seen[base_label]}"
            position_name = safe_name(f"{label}_{chrom}_{start1}_{end1}")
            panel_path = PANEL_DIR / f"{position_name}.png"
            html_path = HTML_DIR / f"{position_name}.html"
            panels.append(
                PhrPanel(
                    bed_index=idx,
                    label=label,
                    chrom_arm=chrom_arm,
                    chrom=chrom,
                    arm=arm,
                    bed_start0=bed_start0,
                    bed_end0=bed_end0,
                    bed_field4=bed_field4,
                    phr_inclusive_bp=phr_len,
                    chrom_size_bp=chrom_size,
                    browser_start1=start1,
                    browser_end1=end1,
                    browser_window_bp=end1 - start1 + 1,
                    terminal_gap_bp=terminal_gap,
                    p_terminal_gap_bp=p_gap,
                    q_terminal_gap_bp=q_gap,
                    ucsc_url=ucsc_url(chrom, start1, end1),
                    html_path=html_path,
                    image_path=panel_path,
                )
            )
    return panels


def extract_img_table(html_text: str) -> str:
    start = html_text.find("<TABLE id='imgTbl'")
    if start == -1:
        raise ValueError("UCSC HTML does not contain TABLE id='imgTbl'")
    end = html_text.find("</TABLE>", start)
    if end == -1:
        raise ValueError("UCSC HTML contains imgTbl start but no closing table")
    return html_text[start : end + len("</TABLE>")]


def parse_slice(block: str, kind: str) -> Slice | None:
    div_match = re.search(
        rf"<div style='width:(\d+)px; height:(\d+)px;' class='sliceDiv {kind}[^']*'>(.*?)</div>",
        block,
        re.S,
    )
    if not div_match:
        return None
    width = int(div_match.group(1))
    height = int(div_match.group(2))
    div_body = div_match.group(3)
    img_match = re.search(r"<IMG\b[^>]*\bsrc='([^']+)'[^>]*\bstyle='([^']+)'", div_body, re.S)
    if not img_match:
        return None
    src = html.unescape(img_match.group(1))
    style = html.unescape(img_match.group(2))
    left_match = re.search(r"left:\s*(-?\d+)px", style)
    top_match = re.search(r"top:\s*(-?\d+)px", style)
    if not left_match or not top_match:
        return None
    left = int(left_match.group(1))
    top = int(top_match.group(1))
    return Slice(src=src, x=-left, y=-top, width=width, height=height)


def parse_ucsc_html(html_text: str, page_url: str) -> tuple[list[BrowserRow], bool, str]:
    decoded = html.unescape(html_text)
    track_confirmed = "chm13.phrs.bed" in decoded and re.search(
        r"<TR id='tr_ct_[^']*chm13phrsbed", html_text
    )
    pos_match = re.search(r"<span class='positionDisplay [^']*'[^>]*>(.*?)</span>", html_text)
    position_display = ""
    if pos_match:
        position_display = html.unescape(re.sub(r"<[^>]+>", "", pos_match.group(1))).strip()
    table = extract_img_table(html_text)
    rows: list[BrowserRow] = []
    for row_match in re.finditer(r"<TR id='tr_([^']+)'[^>]*>(.*?)</TR>", table, re.S):
        row_id = row_match.group(1)
        block = row_match.group(2)
        side = parse_slice(block, "sideLab")
        data = parse_slice(block, "dataImg")
        if side and data:
            rows.append(
                BrowserRow(
                    row_id=row_id,
                    side=Slice(
                        src=urljoin(page_url, side.src),
                        x=side.x,
                        y=side.y,
                        width=side.width,
                        height=side.height,
                    ),
                    data=Slice(
                        src=urljoin(page_url, data.src),
                        x=data.x,
                        y=data.y,
                        width=data.width,
                        height=data.height,
                    ),
                )
            )
    if not rows:
        raise ValueError("UCSC HTML had imgTbl but no crop-able browser rows")
    if not any("chm13phrsbed" in row.row_id for row in rows):
        track_confirmed = False
    return rows, bool(track_confirmed), position_display


def one_value(values: Iterable[int], label: str) -> int:
    unique = sorted(set(values))
    if len(unique) != 1:
        raise ValueError(f"Expected one {label}; observed {unique}")
    return unique[0]


def one_src(values: Iterable[str], label: str) -> str:
    unique = sorted(set(values))
    if len(unique) != 1:
        raise ValueError(f"Expected one {label}; observed {unique}")
    return unique[0]


def identify_dimensions(magick: str, image: Path) -> tuple[int, int]:
    result = run([magick, "identify", "-format", "%w\t%h", str(image)])
    width, height = result.stdout.strip().split("\t")
    return int(width), int(height)


def image_stdev(magick: str, image: Path) -> float:
    result = run(
        [
            magick,
            str(image),
            "-colorspace",
            "Gray",
            "-format",
            "%[fx:standard_deviation]",
            "info:",
        ]
    )
    return float(result.stdout.strip())


def crop_ucsc_panel(panel: PhrPanel, rows: list[BrowserRow], force: bool, magick: str) -> None:
    if panel.image_path.exists() and not force:
        panel.image_width_px, panel.image_height_px = identify_dimensions(magick, panel.image_path)
        panel.image_stdev_gray = image_stdev(magick, panel.image_path)
        return

    data_src = one_src((row.data.src for row in rows), "data sprite source")
    side_src = one_src((row.side.src for row in rows), "side sprite source")

    SPRITE_DIR.mkdir(parents=True, exist_ok=True)
    PANEL_DIR.mkdir(parents=True, exist_ok=True)

    data_sprite = SPRITE_DIR / f"{safe_name(panel.label)}_data.png"
    side_sprite = SPRITE_DIR / f"{safe_name(panel.label)}_side.png"
    fetch_url(data_src, data_sprite, force=True)
    fetch_url(side_src, side_sprite, force=True)

    side_x = one_value((row.side.x for row in rows), "side sprite x offset")
    side_width = one_value((row.side.width for row in rows), "side crop width")
    data_x = one_value((row.data.x for row in rows), "data sprite x offset")
    data_width = one_value((row.data.width for row in rows), "data crop width")
    crop_height = max(
        max(row.side.y + row.side.height for row in rows),
        max(row.data.y + row.data.height for row in rows),
    )

    side_sprite_width, side_sprite_height = identify_dimensions(magick, side_sprite)
    data_sprite_width, data_sprite_height = identify_dimensions(magick, data_sprite)
    if side_x + side_width > side_sprite_width or crop_height > side_sprite_height:
        raise ValueError(
            f"Side crop exceeds sprite for {panel.label}: "
            f"{side_width}x{crop_height}+{side_x}+0 vs {side_sprite_width}x{side_sprite_height}"
        )
    if data_x + data_width > data_sprite_width or crop_height > data_sprite_height:
        raise ValueError(
            f"Data crop exceeds sprite for {panel.label}: "
            f"{data_width}x{crop_height}+{data_x}+0 vs {data_sprite_width}x{data_sprite_height}"
        )

    with tempfile.TemporaryDirectory(prefix=f"{panel.label}_ucsc_") as tmp_dir:
        tmp = Path(tmp_dir)
        side_crop = tmp / "side.png"
        data_crop = tmp / "data.png"
        run(
            [
                magick,
                str(side_sprite),
                "-crop",
                f"{side_width}x{crop_height}+{side_x}+0",
                "+repage",
                "-background",
                "white",
                "-alpha",
                "remove",
                str(side_crop),
            ]
        )
        run(
            [
                magick,
                str(data_sprite),
                "-crop",
                f"{data_width}x{crop_height}+{data_x}+0",
                "+repage",
                "-background",
                "white",
                "-alpha",
                "remove",
                str(data_crop),
            ]
        )
        run(
            [
                magick,
                str(side_crop),
                str(data_crop),
                "+append",
                "-bordercolor",
                "#cfcfcf",
                "-border",
                "1",
                str(panel.image_path),
            ]
        )
    panel.image_width_px, panel.image_height_px = identify_dimensions(magick, panel.image_path)
    panel.image_stdev_gray = image_stdev(magick, panel.image_path)
    if panel.image_width_px < 1700 or panel.image_stdev_gray <= 0.001:
        raise ValueError(
            f"Panel image failed resolution/nonblank check for {panel.label}: "
            f"{panel.image_width_px}x{panel.image_height_px}, stdev={panel.image_stdev_gray}"
        )


def repo_rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def out_rel(path: Path) -> str:
    return path.relative_to(OUT_DIR).as_posix()


def write_manifest(panels: list[PhrPanel]) -> None:
    columns = [
        "label",
        "chrom",
        "arm",
        "chrom_arm",
        "bed_start0",
        "bed_end0",
        "phr_inclusive_bp",
        "chrom_size_bp",
        "browser_start1",
        "browser_end1",
        "browser_window_bp",
        "terminal_gap_bp",
        "p_terminal_gap_bp",
        "q_terminal_gap_bp",
        "bed_field4",
        "image_path",
        "image_width_px",
        "image_height_px",
        "image_stdev_gray",
        "html_path",
        "html_track_confirmed",
        "ucsc_position_display",
        "ucsc_url",
    ]
    with MANIFEST_PATH.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for panel in panels:
            writer.writerow(
                {
                    "label": panel.label,
                    "chrom": panel.chrom,
                    "arm": panel.arm,
                    "chrom_arm": panel.chrom_arm,
                    "bed_start0": panel.bed_start0,
                    "bed_end0": panel.bed_end0,
                    "phr_inclusive_bp": panel.phr_inclusive_bp,
                    "chrom_size_bp": panel.chrom_size_bp,
                    "browser_start1": panel.browser_start1,
                    "browser_end1": panel.browser_end1,
                    "browser_window_bp": panel.browser_window_bp,
                    "terminal_gap_bp": panel.terminal_gap_bp,
                    "p_terminal_gap_bp": panel.p_terminal_gap_bp,
                    "q_terminal_gap_bp": panel.q_terminal_gap_bp,
                    "bed_field4": panel.bed_field4,
                    "image_path": repo_rel(panel.image_path),
                    "image_width_px": panel.image_width_px,
                    "image_height_px": panel.image_height_px,
                    "image_stdev_gray": f"{panel.image_stdev_gray:.6f}",
                    "html_path": repo_rel(panel.html_path),
                    "html_track_confirmed": str(panel.html_track_confirmed).lower(),
                    "ucsc_position_display": panel.ucsc_position_display,
                    "ucsc_url": panel.ucsc_url,
                }
            )


def read_architecture_arms() -> list[str]:
    arms: list[str] = []
    with ARCH_PATH.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            arms.append(row["ChromArm"])
    return arms


def write_audit(panels: list[PhrPanel]) -> tuple[list[str], list[str]]:
    bed_arms = {panel.chrom_arm for panel in panels}
    arch_arms = read_architecture_arms()
    arch_set = set(arch_arms)
    missing = [arm for arm in arch_arms if arm not in bed_arms]
    extra = sorted(bed_arms - arch_set)
    with AUDIT_PATH.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["arm", "in_chm13_phrs_bed", "in_architecture_per_arm", "status", "note"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for arm in arch_arms:
            in_bed = arm in bed_arms
            writer.writerow(
                {
                    "arm": arm,
                    "in_chm13_phrs_bed": str(in_bed).lower(),
                    "in_architecture_per_arm": "true",
                    "status": "rendered" if in_bed else "missing_from_chm13_phrs_bed",
                    "note": (
                        "Rendered from chm13.phrs.bed"
                        if in_bed
                        else "Present in architecture_per_arm.tsv but absent from chm13.phrs.bed; no UCSC panel rendered"
                    ),
                }
            )
        for arm in extra:
            writer.writerow(
                {
                    "arm": arm,
                    "in_chm13_phrs_bed": "true",
                    "in_architecture_per_arm": "false",
                    "status": "extra_in_chm13_phrs_bed",
                    "note": "Rendered from chm13.phrs.bed but not listed in architecture_per_arm.tsv",
                }
            )
    return missing, extra


def typ_str(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def write_typst(panels: list[PhrPanel], missing: list[str], extra: list[str]) -> None:
    lines: list[str] = []
    lines.extend(
        [
            "// Standalone CHM13/hs1 UCSC browser slide suite for PHR review.",
            "// Generated by _scripts/render_ucsc_browser_suite.py.",
            "",
            "#set page(",
            "  width: 13.33in,",
            "  height: 7.5in,",
            "  margin: (x: 0.18in, y: 0.14in),",
            ")",
            '#set text(font: ("DejaVu Sans", "Liberation Sans"), size: 10pt, lang: "en")',
            "#set par(justify: false, leading: 0.64em)",
            "",
            '#let col-title = rgb("#1a3a6b")',
            '#let col-muted = rgb("#555555")',
            '#let col-rule = rgb("#ccd6e0")',
            "",
            "#let browser-slide(num, label, region, phr, gap, img) = {",
            "  grid(",
            "    rows: (0.52in, 1fr, 0.16in),",
            "    row-gutter: 0.02in,",
            "    [",
            "      #grid(",
            "        columns: (1fr, 1fr),",
            "        align: (left, right),",
            "        [#text(size: 12pt, weight: \"semibold\", fill: col-title)[#num #h(0.10in) #label]],",
            "        [#text(size: 9.2pt, fill: col-muted)[#region]],",
            "      )",
            "      #line(length: 100%, stroke: 0.6pt + col-rule)",
            "      #text(size: 7.6pt, fill: col-muted)[#phr #h(0.12in) #gap]",
            "    ],",
            "    box(width: 100%, height: 100%)[",
            "      #image(img, width: 100%, height: 100%, fit: \"contain\")",
            "    ],",
            "    align(center)[#text(size: 5.8pt, fill: col-muted)[UCSC Genome Browser CHM13/hs1; chm13.phrs.bed custom track visible; terminal end included in window]],",
            "  )",
            "}",
            "",
        ]
    )

    missing_text = ", ".join(missing) if missing else "none"
    extra_text = ", ".join(extra) if extra else "none"
    lines.extend(
        [
            "#align(center + horizon)[",
            "  #block(width: 11.2in)[",
            "    #align(center)[#text(size: 28pt, weight: \"bold\", fill: col-title)[CHM13/hs1 PHR UCSC Browser Suite]]",
            "    #v(0.28in)",
            f"    #text(size: 14pt)[Rendered {len(panels)} main UCSC browser-image slides from `chm13.phrs.bed`. Every panel uses a terminal-anchored 1.5x zoom-out around the PHR interval and keeps the chromosome end in view.]",
            "    #v(0.18in)",
            f"    #text(size: 12pt, fill: col-muted)[Audit: `architecture_per_arm.tsv` lists {len(read_architecture_arms())} arms, while `chm13.phrs.bed` contains {len(panels)} intervals. Missing from the CHM13 BED: {missing_text}. Extra in the CHM13 BED relative to the architecture audit: {extra_text}.]",
            "    #v(0.18in)",
            "    #text(size: 10.5pt, fill: col-muted)[Manifest: `manifest.tsv`; audit: `audit_37_vs_41.tsv`; cached UCSC HTML and cropped panels: `_assets/ucsc/`.]",
            "  ]",
            "]",
            "",
        ]
    )
    for i, panel in enumerate(panels, start=1):
        lines.append("#pagebreak()")
        region = f"{panel.chrom}:{panel.browser_start1}-{panel.browser_end1}"
        phr = (
            f"PHR BED {panel.chrom}:{panel.bed_start0}-{panel.bed_end0}; "
            f"inclusive PHR {panel.phr_inclusive_bp:,} bp; browser window {panel.browser_window_bp:,} bp"
        )
        gap = f"{panel.arm}-arm terminal gap {panel.terminal_gap_bp:,} bp"
        lines.append(
            "#browser-slide("
            + ", ".join(
                [
                    typ_str(f"{i:02d}"),
                    typ_str(panel.label),
                    typ_str(region),
                    typ_str(phr),
                    typ_str(gap),
                    typ_str(out_rel(panel.image_path)),
                ]
            )
            + ")"
        )
        lines.append("")
    TYP_PATH.write_text("\n".join(lines))


def compile_typst(typst: str) -> None:
    run([typst, "compile", "--root", str(OUT_DIR), str(TYP_PATH), str(PDF_PATH)])


def render(force: bool) -> list[PhrPanel]:
    os.environ["PATH"] = "/home/erikg/micromamba/bin:" + os.environ.get("PATH", "")
    magick = tool_path("magick", "/home/erikg/micromamba/bin/magick")
    typst = tool_path("typst", "/home/erikg/.local/bin/typst")

    for directory in [OUT_DIR, ASSET_DIR, HTML_DIR, SPRITE_DIR, PANEL_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    sizes = read_chrom_sizes(force=force)
    panels = read_bed_panels(sizes)

    for panel in panels:
        print(
            f"[panel] {panel.label}: {panel.chrom}:{panel.browser_start1}-{panel.browser_end1}",
            flush=True,
        )
        fetch_url(panel.ucsc_url, panel.html_path, force=force or not panel.image_path.exists())
        html_text = panel.html_path.read_text(errors="replace")
        rows, track_confirmed, position_display = parse_ucsc_html(html_text, panel.ucsc_url)
        if not track_confirmed:
            raise ValueError(f"UCSC HTML did not confirm visible chm13.phrs.bed for {panel.label}")
        panel.html_track_confirmed = track_confirmed
        panel.ucsc_position_display = position_display
        crop_ucsc_panel(panel, rows, force=force, magick=magick)

    write_manifest(panels)
    missing, extra = write_audit(panels)
    write_typst(panels, missing, extra)
    compile_typst(typst)
    return panels


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force",
        action="store_true",
        help="refetch UCSC HTML/sprites and rebuild panel PNGs even if cached files exist",
    )
    args = parser.parse_args(argv)
    panels = render(force=args.force)
    chr9q = [p for p in panels if p.label == "chr9q"]
    if not chr9q:
        raise SystemExit("No chr9q panel was generated")
    got = f"{chr9q[0].chrom}:{chr9q[0].browser_start1}-{chr9q[0].browser_end1}"
    if got != "chr9:150279748-150617247":
        raise SystemExit(f"chr9q window mismatch: got {got}")
    print(f"[done] rendered {len(panels)} browser panels plus 1 audit/title slide")
    print(f"[done] {PDF_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
