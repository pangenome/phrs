#!/usr/bin/env python3
"""
Render a PAN027 maternal-vs-PAN010 untangle source/native visual diff.

The script only writes compact repo artifacts. Native odgi outputs remain in
their external Slurm directory, and sweepGA PAFs are generated in the external
sweepGA test directory when absent.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import math
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = REPO_ROOT / "paper_prep" / "_brainstorming" / "fig5_maternal_native_diff"

OLD_BED = Path(
    "/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/"
    "PAN027_vs_PAN010.e50000.m1000.bed.gz"
)
NATIVE_DIR = Path("/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm")
NATIVE_BEDPE = NATIVE_DIR / "PAN027_vs_PAN010.e50000.m1000.j0.8.n4.bedpe.gz"
NATIVE_PAF = NATIVE_DIR / "PAN027_vs_PAN010.e50000.m1000.j0.8.n4.paf.gz"
SWEEPGA_DIR = Path("/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_sweepga_test")
SWEEPGA = Path("/home/erikg/.cargo/bin/sweepga")

PAIR = "PAN027_vs_PAN010"
QUERY_SAMPLE = "PAN027"
QUERY_HAP = "1"
TARGET_SAMPLE = "PAN010"
BIN_SIZE = 1000
QUERY_BP = 500_000

VIEW_ORDER = [
    "old_bed_n1",
    "native_bedpe_n1",
    "native_bedpe_n4",
    "native_paf_n1",
    "native_paf_n4",
    "sweepga_paf_1_many",
    "sweepga_paf_2_many",
    "sweepga_paf_4_many",
]

VIEW_LABELS = {
    "old_bed_n1": "old BED first-best",
    "native_bedpe_n1": "native BEDPE n1",
    "native_bedpe_n4": "native BEDPE n4",
    "native_paf_n1": "native PAF n1",
    "native_paf_n4": "native PAF n4",
    "sweepga_paf_1_many": "sweepGA PAF 1:many",
    "sweepga_paf_2_many": "sweepGA PAF 2:many",
    "sweepga_paf_4_many": "sweepGA PAF 4:many",
}

SWEEPGA_PREFIX = "sweepga_paf_"

ARM_ORDER = [
    f"chr{i}{arm}" for i in range(1, 23) for arm in ("p", "q")
] + ["chrXp", "chrXq", "chrYp", "chrYq"]

ARM_RE = re.compile(r"^(PAN\d+)#([12])#.*_(chr(?:[0-9]+|X|Y))_([pq])arm$")
TAG_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]):[A-Za-z]:(.+)$")


@dataclass(frozen=True)
class Segment:
    source: str
    qname: str
    qarm: str
    qstart: int
    qend: int
    target_name: str
    target_hap: str
    target_arm: str
    score: float
    nth: int
    strand: str

    @property
    def length(self) -> int:
        return max(0, self.qend - self.qstart)

    @property
    def interchromosomal(self) -> bool:
        return self.qarm[: self.qarm.rfind(self.qarm[-1])] != self.target_arm[: self.target_arm.rfind(self.target_arm[-1])]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--old-bed", type=Path, default=OLD_BED)
    parser.add_argument("--native-bedpe", type=Path, default=NATIVE_BEDPE)
    parser.add_argument("--native-paf", type=Path, default=NATIVE_PAF)
    parser.add_argument("--sweepga-dir", type=Path, default=SWEEPGA_DIR)
    parser.add_argument("--skip-sweepga", action="store_true", help="Do not generate missing external sweepGA PAFs")
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def path_info(name: str) -> tuple[str, str, str] | None:
    match = ARM_RE.search(name)
    if not match:
        return None
    sample, hap, chrom, arm = match.groups()
    return sample, hap, f"{chrom}{arm}"


def chrom_of(chrarm: str) -> str:
    return chrarm[:-1]


def parse_float(value: str, default: float = math.nan) -> float:
    try:
        return float(value)
    except ValueError:
        return default


def parse_int(value: str, default: int = 0) -> int:
    try:
        return int(value)
    except ValueError:
        return default


def read_bed_like(path: Path, source: str, max_nth: int = 4) -> list[Segment]:
    rows: list[Segment] = []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            qinfo = path_info(fields[0])
            tinfo = path_info(fields[3])
            if not qinfo or not tinfo:
                continue
            qsample, qhap, qarm = qinfo
            tsample, thap, tarm = tinfo
            if qsample != QUERY_SAMPLE or qhap != QUERY_HAP or tsample != TARGET_SAMPLE:
                continue
            nth = parse_int(fields[9], 1)
            if nth > max_nth:
                continue
            rows.append(
                Segment(
                    source=source,
                    qname=fields[0],
                    qarm=qarm,
                    qstart=max(0, parse_int(fields[1])),
                    qend=min(QUERY_BP, parse_int(fields[2])),
                    target_name=fields[3],
                    target_hap=f"{tsample}#{thap}",
                    target_arm=tarm,
                    score=parse_float(fields[6]),
                    nth=nth,
                    strand=fields[7],
                )
            )
    return rows


def read_paf(path: Path, source: str, max_nth: int = 4) -> list[Segment]:
    rows: list[Segment] = []
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
            if qsample != QUERY_SAMPLE or qhap != QUERY_HAP or tsample != TARGET_SAMPLE:
                continue
            tags = {}
            for tag in fields[12:]:
                match = TAG_RE.match(tag)
                if match:
                    tags[match.group(1)] = match.group(2)
            nth = parse_int(tags.get("nb", "1"), 1)
            if nth > max_nth:
                continue
            score = parse_float(tags.get("jc", tags.get("id", "nan")))
            rows.append(
                Segment(
                    source=source,
                    qname=fields[0],
                    qarm=qarm,
                    qstart=max(0, parse_int(fields[2])),
                    qend=min(QUERY_BP, parse_int(fields[3])),
                    target_name=fields[5],
                    target_hap=f"{tsample}#{thap}",
                    target_arm=tarm,
                    score=score,
                    nth=nth,
                    strand=fields[4],
                )
            )
    return rows


def ensure_sweepga_outputs(native_paf: Path, sweepga_dir: Path, skip: bool) -> dict[str, Path]:
    sweepga_dir.mkdir(parents=True, exist_ok=True)
    native_uncompressed = sweepga_dir / "PAN027_vs_PAN010.e50000.m1000.j0.8.n4.native.paf"
    outputs = {
        "sweepga_paf_1_many": sweepga_dir / "PAN027_vs_PAN010.e50000.m1000.j0.8.n4.uncompressed.1_many.paf",
        "sweepga_paf_2_many": sweepga_dir / "PAN027_vs_PAN010.e50000.m1000.j0.8.n4.uncompressed.2_many.paf",
        "sweepga_paf_4_many": sweepga_dir / "PAN027_vs_PAN010.e50000.m1000.j0.8.n4.uncompressed.4_many.paf",
    }
    if skip:
        return {k: v for k, v in outputs.items() if v.exists()}
    if not SWEEPGA.exists():
        return {k: v for k, v in outputs.items() if v.exists()}
    if not native_uncompressed.exists():
        with gzip.open(native_paf, "rb") as src, native_uncompressed.open("wb") as dst:
            dst.write(src.read())
    for view, out_path in outputs.items():
        if out_path.exists():
            continue
        num = view[len(SWEEPGA_PREFIX) :].replace("_", ":")
        subprocess.run(
            [
                str(SWEEPGA),
                "--num-mappings",
                num,
                "--scaffold-jump",
                "0",
                "--output-file",
                str(out_path),
                str(native_uncompressed),
            ],
            check=True,
        )
    return {k: v for k, v in outputs.items() if v.exists()}


def coalesce_segments(rows: list[Segment]) -> list[Segment]:
    grouped = defaultdict(list)
    for row in sorted(rows, key=lambda r: (r.source, r.qname, r.nth, r.target_arm, r.target_hap, r.qstart, r.qend)):
        key = (row.source, row.qname, row.qarm, row.nth, row.target_arm, row.target_hap, row.strand)
        bucket = grouped[key]
        if bucket and row.qstart <= bucket[-1].qend + 1:
            prev = bucket[-1]
            bucket[-1] = Segment(
                prev.source,
                prev.qname,
                prev.qarm,
                prev.qstart,
                max(prev.qend, row.qend),
                prev.target_name,
                prev.target_hap,
                prev.target_arm,
                max(prev.score, row.score),
                prev.nth,
                prev.strand,
            )
        else:
            bucket.append(row)
    merged: list[Segment] = []
    for bucket in grouped.values():
        merged.extend(bucket)
    return sorted(merged, key=lambda r: (VIEW_ORDER.index(r.source), ARM_ORDER.index(r.qarm), r.qstart, r.nth, r.target_arm))


def bin_sets(rows: list[Segment]) -> dict[tuple[str, int], set[str]]:
    bins: dict[tuple[str, int], set[str]] = defaultdict(set)
    for row in rows:
        for idx in range(row.qstart // BIN_SIZE, max(row.qstart // BIN_SIZE, (row.qend - 1) // BIN_SIZE) + 1):
            lo = idx * BIN_SIZE
            hi = lo + BIN_SIZE
            if row.qstart < hi and row.qend > lo:
                bins[(row.qarm, idx)].add(row.target_arm)
    return bins


def summarize_views(rows_by_view: dict[str, list[Segment]]) -> list[dict[str, str]]:
    old_bins = bin_sets(rows_by_view.get("old_bed_n1", []))
    summaries: list[dict[str, str]] = []
    for view in VIEW_ORDER:
        rows = rows_by_view.get(view, [])
        by_arm = defaultdict(list)
        for row in rows:
            by_arm[row.qarm].append(row)
        bins = bin_sets(rows)
        for qarm in ARM_ORDER:
            arm_rows = by_arm.get(qarm, [])
            if not arm_rows and qarm not in {arm for arm, _ in bins}:
                continue
            target_arms = sorted({r.target_arm for r in arm_rows}, key=arm_sort_key)
            target_haps = sorted({f"{r.target_hap}:{r.target_arm}" for r in arm_rows}, key=lambda x: (arm_sort_key(x.split(":")[-1]), x))
            rank2plus_bins = bin_sets([r for r in arm_rows if r.nth > 1])
            rank2plus_bp = binned_bp(rank2plus_bins, lambda values: bool(values))
            multi_target_bp = binned_bp({k: v for k, v in bins.items() if k[0] == qarm}, lambda values: len(values) > 1)
            missing_old_bp = 0
            added_targets: set[str] = set()
            if view != "old_bed_n1":
                for idx in range(QUERY_BP // BIN_SIZE):
                    old_targets = old_bins.get((qarm, idx), set())
                    new_targets = bins.get((qarm, idx), set())
                    if old_targets and new_targets and old_targets.isdisjoint(new_targets):
                        missing_old_bp += BIN_SIZE
                    added_targets.update(new_targets - old_targets)
            summaries.append(
                {
                    "view": view,
                    "view_label": VIEW_LABELS[view],
                    "query_arm": qarm,
                    "n_segments": str(len(arm_rows)),
                    "covered_bp_binned": str(binned_bp({k: v for k, v in bins.items() if k[0] == qarm}, lambda values: bool(values))),
                    "interchromosomal_bp_binned": str(
                        binned_bp(
                            {k: v for k, v in bins.items() if k[0] == qarm},
                            lambda values, q=qarm: any(chrom_of(v) != chrom_of(q) for v in values),
                        )
                    ),
                    "rank2plus_present_bp_binned": str(rank2plus_bp),
                    "multi_target_bp_binned": str(multi_target_bp),
                    "old_n1_target_absent_bp_binned": str(missing_old_bp),
                    "target_arms": ",".join(target_arms),
                    "target_hap_arms": ",".join(target_haps[:20]) + (";..." if len(target_haps) > 20 else ""),
                    "added_target_arms_vs_old": ",".join(sorted(added_targets, key=arm_sort_key)),
                }
            )
    return summaries


def binned_bp(bins: dict[tuple[str, int], set[str]], pred) -> int:
    return sum(BIN_SIZE for values in bins.values() if pred(values))


def arm_sort_key(arm: str) -> tuple[int, int]:
    chrom = chrom_of(arm)
    suffix = arm[-1]
    if chrom == "chrX":
        c = 23
    elif chrom == "chrY":
        c = 24
    else:
        c = int(chrom[3:])
    return c, 0 if suffix == "p" else 1


def color_for_arm(arm: str) -> str:
    c, side = arm_sort_key(arm)
    hue = (c * 47 + side * 17) % 360
    sat = 54 if side == 0 else 70
    light = 42 if side == 0 else 55
    return hsl_to_hex(hue, sat / 100, light / 100)


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


def render_svg(rows_by_view: dict[str, list[Segment]], path: Path) -> None:
    panel_w = 205
    row_h = 8
    top = 86
    left = 72
    gap = 16
    width = left + len(VIEW_ORDER) * panel_w + (len(VIEW_ORDER) - 1) * gap + 26
    height = top + len(ARM_ORDER) * row_h + 72
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#222}.title{font-size:16px;font-weight:700}.small{font-size:9px}.label{font-size:8px}.axis{stroke:#999;stroke-width:.5}.panel{fill:#fafafa;stroke:#ddd}.hit{stroke:#111;stroke-width:.08}.alt{stroke:#111;stroke-width:.25;stroke-dasharray:1 1}.miss{fill:#fff;stroke:#d73027;stroke-width:.7}",
        "</style>",
        f'<text class="title" x="{left}" y="24">PAN027 maternal haplotype vs PAN010: old/source vs native odgi and sweepGA views</text>',
        f'<text class="small" x="{left}" y="43">Rows are PAN027#1 query arms; x-axis is each 500 kb telomere-anchored flank. Color encodes PAN010 target arm. Thin sublanes show nth-best ranks; red outlines flag bins where a regenerated view lacks the old first-best target arm.</text>',
    ]
    old_bins = bin_sets(rows_by_view.get("old_bed_n1", []))
    for i, view in enumerate(VIEW_ORDER):
        x0 = left + i * (panel_w + gap)
        lines.append(f'<rect class="panel" x="{x0-2}" y="{top-24}" width="{panel_w+4}" height="{len(ARM_ORDER)*row_h+34}" rx="2"/>')
        lines.append(f'<text class="small" x="{x0}" y="{top-11}">{html.escape(VIEW_LABELS[view])}</text>')
        lines.append(f'<line class="axis" x1="{x0}" y1="{top-5}" x2="{x0+panel_w}" y2="{top-5}"/>')
        lines.append(f'<text class="label" x="{x0}" y="{top-29}">0</text><text class="label" x="{x0+panel_w-22}" y="{top-29}">500 kb</text>')
        for idx, arm in enumerate(ARM_ORDER):
            y = top + idx * row_h
            if i == 0:
                lines.append(f'<text class="label" x="{left-8}" y="{y+5}" text-anchor="end">{arm}</text>')
            lines.append(f'<line class="axis" x1="{x0}" y1="{y+row_h-1}" x2="{x0+panel_w}" y2="{y+row_h-1}" opacity=".20"/>')
        for row in rows_by_view.get(view, []):
            if row.qarm not in ARM_ORDER:
                continue
            ybase = top + ARM_ORDER.index(row.qarm) * row_h
            lane_h = 1.45
            lane_y = ybase + min(row.nth - 1, 3) * lane_h + 0.7
            x = x0 + (row.qstart / QUERY_BP) * panel_w
            w = max(0.45, (row.length / QUERY_BP) * panel_w)
            klass = "hit" if row.nth == 1 else "hit alt"
            opacity = "0.92" if row.nth == 1 else "0.55"
            lines.append(
                f'<rect class="{klass}" x="{x:.2f}" y="{lane_y:.2f}" width="{w:.2f}" height="{lane_h:.2f}" '
                f'fill="{color_for_arm(row.target_arm)}" opacity="{opacity}"><title>{html.escape(view)} {html.escape(row.qarm)} {row.qstart}-{row.qend} -> {html.escape(row.target_hap)}:{html.escape(row.target_arm)} rank {row.nth} score {row.score:.4g}</title></rect>'
            )
        if view != "old_bed_n1":
            bins = bin_sets(rows_by_view.get(view, []))
            for arm in ARM_ORDER:
                for idx in range(QUERY_BP // BIN_SIZE):
                    old_targets = old_bins.get((arm, idx), set())
                    new_targets = bins.get((arm, idx), set())
                    if old_targets and new_targets and old_targets.isdisjoint(new_targets):
                        y = top + ARM_ORDER.index(arm) * row_h + 0.4
                        x = x0 + (idx * BIN_SIZE / QUERY_BP) * panel_w
                        w = max(0.5, (BIN_SIZE / QUERY_BP) * panel_w)
                        lines.append(f'<rect class="miss" x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{row_h-1:.2f}" opacity=".8"/>')
    legend_y = height - 43
    legend = ["chr3p", "chr4p", "chr13p", "chr14p", "chr15p", "chr21p", "chr22p", "chrXp"]
    lines.append(f'<text class="small" x="{left}" y="{legend_y}">Target-arm color examples:</text>')
    lx = left + 126
    for arm in legend:
        lines.append(f'<rect x="{lx}" y="{legend_y-8}" width="12" height="8" fill="{color_for_arm(arm)}"/>')
        lines.append(f'<text class="label" x="{lx+15}" y="{legend_y-1}">{arm}</text>')
        lx += 54
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n")


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    sweepga_paths = ensure_sweepga_outputs(args.native_paf, args.sweepga_dir, args.skip_sweepga)

    raw_by_view = {
        "old_bed_n1": [r for r in read_bed_like(args.old_bed, "old_bed_n1", max_nth=1) if r.nth == 1],
        "native_bedpe_n1": [r for r in read_bed_like(args.native_bedpe, "native_bedpe_n1", max_nth=1) if r.nth == 1],
        "native_bedpe_n4": read_bed_like(args.native_bedpe, "native_bedpe_n4", max_nth=4),
        "native_paf_n1": [r for r in read_paf(args.native_paf, "native_paf_n1", max_nth=1) if r.nth == 1],
        "native_paf_n4": read_paf(args.native_paf, "native_paf_n4", max_nth=4),
    }
    for view, path in sweepga_paths.items():
        raw_by_view[view] = read_paf(path, view, max_nth=4)
    rows_by_view = {view: coalesce_segments(raw_by_view.get(view, [])) for view in VIEW_ORDER}

    source_rows = [
        {"view": "old_bed_n1", "role": "current/source first-best filtered to PAN027#1", "path": str(args.old_bed), "rows_loaded": str(len(raw_by_view["old_bed_n1"]))},
        {"view": "native_bedpe_n1", "role": "native odgi BEDPE nth.best=1", "path": str(args.native_bedpe), "rows_loaded": str(len(raw_by_view["native_bedpe_n1"]))},
        {"view": "native_bedpe_n4", "role": "native odgi BEDPE nth.best<=4", "path": str(args.native_bedpe), "rows_loaded": str(len(raw_by_view["native_bedpe_n4"]))},
        {"view": "native_paf_n1", "role": "native odgi PAF nb=1", "path": str(args.native_paf), "rows_loaded": str(len(raw_by_view["native_paf_n1"]))},
        {"view": "native_paf_n4", "role": "native odgi PAF nb<=4", "path": str(args.native_paf), "rows_loaded": str(len(raw_by_view["native_paf_n4"]))},
    ]
    for view in ("sweepga_paf_1_many", "sweepga_paf_2_many", "sweepga_paf_4_many"):
        source_rows.append(
            {
                "view": view,
                "role": f"sweepGA-filtered native PAF {view[len(SWEEPGA_PREFIX):].replace('_', ':')}",
                "path": str(sweepga_paths.get(view, "")),
                "rows_loaded": str(len(raw_by_view.get(view, []))),
            }
        )
    write_tsv(args.out_dir / "source_manifest.tsv", source_rows, ["view", "role", "path", "rows_loaded"])

    seg_rows = []
    for view in VIEW_ORDER:
        for row in rows_by_view.get(view, []):
            seg_rows.append(
                {
                    "view": view,
                    "query_name": row.qname,
                    "query_arm": row.qarm,
                    "query_start": str(row.qstart),
                    "query_end": str(row.qend),
                    "target_hap": row.target_hap,
                    "target_arm": row.target_arm,
                    "nth_best": str(row.nth),
                    "score": f"{row.score:.6g}",
                    "strand": row.strand,
                }
            )
    write_tsv(
        args.out_dir / "segments_compact.tsv",
        seg_rows,
        ["view", "query_name", "query_arm", "query_start", "query_end", "target_hap", "target_arm", "nth_best", "score", "strand"],
    )

    summary = summarize_views(rows_by_view)
    write_tsv(
        args.out_dir / "arm_comparison.tsv",
        summary,
        [
            "view",
            "view_label",
            "query_arm",
            "n_segments",
            "covered_bp_binned",
            "interchromosomal_bp_binned",
            "rank2plus_present_bp_binned",
            "multi_target_bp_binned",
            "old_n1_target_absent_bp_binned",
            "target_arms",
            "target_hap_arms",
            "added_target_arms_vs_old",
        ],
    )
    render_svg(rows_by_view, args.out_dir / "fig5_maternal_native_diff.svg")


if __name__ == "__main__":
    main()
