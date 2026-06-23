#!/usr/bin/env python3
"""Build Fig5 raw-FASTA f16 chop/filter sensitivity summaries.

The heavy PAF layers are intentionally ignored by git. This script reads the
whole-genome chopped PAFs, applies sweepGA filtering genome-wide, and only then
subsets the Fig5 candidate query windows for the committed summary/panel data.
"""

import argparse
import csv
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[2]
DEFAULT_SOURCE_PACKAGE = (
    Path("/moosefs/erikg/phrs")
    / "paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16"
)
CHOP_LENGTHS = (2000, 5000, 10000)


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def chrom_from_name(name: str) -> str:
    matches = re.findall(r"chr(?:[0-9]+|X|Y|M)", name)
    return matches[-1] if matches else name


def overlap_bp(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def union_bp(intervals) -> int:
    return sum(end - start for start, end in merge_intervals(intervals))


def existing_ani_path(source_package: Path, comparison_id: str, chop_length: int):
    path = (
        source_package
        / "filtered_paf_chop_sensitivity"
        / f"{comparison_id}.l{chop_length}.one_one_chopped.paf.gz"
    )
    return path if path.exists() else None


def filtered_path(panel_dir: Path, comparison_id: str, chop_length: int, filter_mode: str) -> Path:
    return panel_dir / "filtered_paf" / f"{comparison_id}.l{chop_length}.{filter_mode}.paf.gz"


def run_sweepga_filter(
    *,
    sweepga: str,
    src_gz: Path,
    out_gz: Path,
    mode: dict,
    scratch_base: Path,
) -> None:
    out_gz.parent.mkdir(parents=True, exist_ok=True)
    if out_gz.exists():
        return
    if not src_gz.exists():
        raise FileNotFoundError(f"missing chopped PAF: {src_gz}")
    if not scratch_base.is_dir() or not os.access(scratch_base, os.W_OK):
        raise SystemExit(f"scratch base is not writable: {scratch_base}")
    work = Path(tempfile.mkdtemp(prefix=f"sweepga.{out_gz.stem}.", dir=str(scratch_base)))
    try:
        input_paf = work / "input.paf"
        output_paf = work / "filtered.paf"
        with gzip.open(src_gz, "rb") as inp, input_paf.open("wb") as out:
            shutil.copyfileobj(inp, out)
        cmd = [
            sweepga,
            "--num-mappings",
            mode["num_mappings"],
            "--scaffold-jump",
            mode["scaffold_jump"],
            "--scoring",
            mode["scoring"],
            "--output-file",
            str(output_paf),
            str(input_paf),
        ]
        subprocess.run(cmd, check=True)
        with output_paf.open("rb") as inp, gzip.open(out_gz, "wb") as out:
            shutil.copyfileobj(inp, out)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def ensure_filtered_layers(args, comparisons, modes):
    manifest = []
    source_package = args.source_package
    scratch = Path(args.scratch)
    for comparison_id in comparisons:
        for chop_length in CHOP_LENGTHS:
            src = (
                source_package
                / f"chopped_paf_l{chop_length}_o0"
                / f"{comparison_id}.chopped_l{chop_length}_o0.paf.gz"
            )
            for mode_id, mode in modes.items():
                out = filtered_path(args.panel_dir, comparison_id, chop_length, mode_id)
                reused_from = ""
                existing = existing_ani_path(source_package, comparison_id, chop_length)
                if mode_id == "no_merge_ani" and not out.exists() and existing:
                    out.parent.mkdir(parents=True, exist_ok=True)
                    os.symlink(str(existing), str(out))
                    reused_from = str(existing)
                if not out.exists():
                    run_sweepga_filter(
                        sweepga=args.sweepga,
                        src_gz=src,
                        out_gz=out,
                        mode=mode,
                        scratch_base=scratch,
                    )
                manifest.append(
                    {
                        "comparison_id": comparison_id,
                        "chop_length_bp": str(chop_length),
                        "filter_mode": mode_id,
                        "filter_label": mode["label"],
                        "num_mappings": mode["num_mappings"],
                        "scaffold_jump": mode["scaffold_jump"],
                        "scoring": mode["scoring"],
                        "source_chopped_paf": str(src),
                        "filtered_paf": str(out),
                        "filtered_paf_sha256": sha256(out.resolve()),
                        "reused_from": reused_from,
                        "status": "OK",
                    }
                )
                sys.stderr.write(
                    f"[ok] {comparison_id} l{chop_length} {mode_id}: {out}\n"
                )
    return manifest


def summarize(panel_dir: Path, windows, modes):
    segments = []
    summary = []
    target_breakdown = []
    for window in windows:
        event_id = window["event_id"]
        comparison_id = window["comparison_id"]
        query_name = window["query_name"]
        q0 = int(window["query_start"])
        q1 = int(window["query_end"])
        expected = window.get("expected_target_chroms") or window.get("expected_target_chrom", "chr3")
        expected_set = {x for x in expected.split(",") if x}
        for chop_length in CHOP_LENGTHS:
            for mode_id, mode in modes.items():
                paf = filtered_path(panel_dir, comparison_id, chop_length, mode_id)
                rows = []
                by_target = defaultdict(list)
                by_target_sum = defaultdict(int)
                with gzip.open(paf, "rt") as handle:
                    for line_no, line in enumerate(handle, start=1):
                        if not line.strip() or line.startswith("#"):
                            continue
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) < 12 or fields[0] != query_name:
                            continue
                        qs = int(fields[2])
                        qe = int(fields[3])
                        ov = overlap_bp(qs, qe, q0, q1)
                        if ov <= 0:
                            continue
                        clipped_start = max(qs, q0)
                        clipped_end = min(qe, q1)
                        target_chrom = chrom_from_name(fields[5])
                        matches = int(fields[9])
                        aln_len = int(fields[10])
                        identity = matches / aln_len if aln_len else 0.0
                        row = {
                            "event_id": event_id,
                            "comparison_id": comparison_id,
                            "filter_mode": mode_id,
                            "filter_label": mode["label"],
                            "chop_length_bp": str(chop_length),
                            "query_name": query_name,
                            "query_chrom": chrom_from_name(query_name),
                            "window_start_abs": str(q0),
                            "window_end_abs": str(q1),
                            "query_start_abs": str(clipped_start),
                            "query_end_abs": str(clipped_end),
                            "window_overlap_bp": str(ov),
                            "target_name": fields[5],
                            "target_chrom": target_chrom,
                            "target_start": fields[7],
                            "target_end": fields[8],
                            "strand": fields[4],
                            "matches": str(matches),
                            "alignment_length": str(aln_len),
                            "identity": f"{identity:.6f}",
                            "expected_target_chroms": expected,
                            "is_expected_target": "yes" if target_chrom in expected_set else "no",
                            "filtered_paf": str(paf),
                            "filtered_paf_line": str(line_no),
                        }
                        rows.append(row)
                        by_target[target_chrom].append((clipped_start, clipped_end))
                        by_target_sum[target_chrom] += ov
                segments.extend(rows)
                chr3_rows = [r for r in rows if r["target_chrom"] == "chr3"]
                chr9_rows = [r for r in rows if r["target_chrom"] == "chr9"]
                other_intervals = []
                other_sum = 0
                for chrom, intervals in by_target.items():
                    if chrom not in {"chr3", "chr9"}:
                        other_intervals.extend(intervals)
                        other_sum += by_target_sum[chrom]
                    target_breakdown.append(
                        {
                            "event_id": event_id,
                            "comparison_id": comparison_id,
                            "filter_mode": mode_id,
                            "filter_label": mode["label"],
                            "chop_length_bp": str(chop_length),
                            "query_name": query_name,
                            "window_start_abs": str(q0),
                            "window_end_abs": str(q1),
                            "target_chrom": chrom,
                            "overlap_rows": str(sum(1 for r in rows if r["target_chrom"] == chrom)),
                            "sum_overlap_bp": str(by_target_sum[chrom]),
                            "query_union_bp": str(union_bp(intervals)),
                            "paf_path": str(paf),
                        }
                    )
                summary.append(
                    {
                        "event_id": event_id,
                        "comparison_id": comparison_id,
                        "filter_mode": mode_id,
                        "filter_label": mode["label"],
                        "chop_length_bp": str(chop_length),
                        "query_name": query_name,
                        "query_chrom": chrom_from_name(query_name),
                        "window_start_abs": str(q0),
                        "window_end_abs": str(q1),
                        "query_window_bp": str(q1 - q0),
                        "chr3_rows": str(len(chr3_rows)),
                        "chr3_sum_overlap_bp": str(sum(int(r["window_overlap_bp"]) for r in chr3_rows)),
                        "chr3_query_union_bp": str(union_bp((int(r["query_start_abs"]), int(r["query_end_abs"])) for r in chr3_rows)),
                        "chr9_rows": str(len(chr9_rows)),
                        "chr9_sum_overlap_bp": str(sum(int(r["window_overlap_bp"]) for r in chr9_rows)),
                        "chr9_query_union_bp": str(union_bp((int(r["query_start_abs"]), int(r["query_end_abs"])) for r in chr9_rows)),
                        "other_sum_overlap_bp": str(other_sum),
                        "other_query_union_bp": str(union_bp(other_intervals)),
                        "chr3_survives": "yes" if chr3_rows else "no",
                        "status": "OK" if chr3_rows else "NO_CHR3_ROWS",
                    }
                )
    return segments, summary, target_breakdown


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-dir", type=Path, default=HERE)
    parser.add_argument("--source-package", type=Path, default=DEFAULT_SOURCE_PACKAGE)
    parser.add_argument("--sweepga", default=os.environ.get("SWEEPGA", "/home/erikg/.cargo/bin/sweepga"))
    parser.add_argument("--scratch", default=os.environ.get("SWEEPGA_DEVSHM_BASE", "/dev/shm"))
    parser.add_argument("--skip-filter", action="store_true")
    args = parser.parse_args()

    windows = [
        row
        for row in read_tsv(args.panel_dir / "config/panel_windows.tsv")
        if "chr3" in (row.get("expected_target_chroms") or row.get("expected_target_chrom", "")).split(",")
    ]
    modes = {row["filter_mode"]: row for row in read_tsv(args.panel_dir / "config/filter_modes.tsv")}
    comparisons = sorted({row["comparison_id"] for row in windows})
    if not args.skip_filter:
        manifest = ensure_filtered_layers(args, comparisons, modes)
    else:
        manifest = []
    segments, summary, target_breakdown = summarize(args.panel_dir, windows, modes)

    write_tsv(
        args.panel_dir / "chop_filter_panel_segments.tsv",
        segments,
        [
            "event_id",
            "comparison_id",
            "filter_mode",
            "filter_label",
            "chop_length_bp",
            "query_name",
            "query_chrom",
            "window_start_abs",
            "window_end_abs",
            "query_start_abs",
            "query_end_abs",
            "window_overlap_bp",
            "target_name",
            "target_chrom",
            "target_start",
            "target_end",
            "strand",
            "matches",
            "alignment_length",
            "identity",
            "expected_target_chroms",
            "is_expected_target",
            "filtered_paf",
            "filtered_paf_line",
        ],
    )
    write_tsv(
        args.panel_dir / "chop_filter_panel_summary.tsv",
        summary,
        [
            "event_id",
            "comparison_id",
            "filter_mode",
            "filter_label",
            "chop_length_bp",
            "query_name",
            "query_chrom",
            "window_start_abs",
            "window_end_abs",
            "query_window_bp",
            "chr3_rows",
            "chr3_sum_overlap_bp",
            "chr3_query_union_bp",
            "chr9_rows",
            "chr9_sum_overlap_bp",
            "chr9_query_union_bp",
            "other_sum_overlap_bp",
            "other_query_union_bp",
            "chr3_survives",
            "status",
        ],
    )
    write_tsv(
        args.panel_dir / "chop_filter_target_breakdown.tsv",
        target_breakdown,
        [
            "event_id",
            "comparison_id",
            "filter_mode",
            "filter_label",
            "chop_length_bp",
            "query_name",
            "window_start_abs",
            "window_end_abs",
            "target_chrom",
            "overlap_rows",
            "sum_overlap_bp",
            "query_union_bp",
            "paf_path",
        ],
    )
    if manifest:
        write_tsv(
            args.panel_dir / "filter_manifest.tsv",
            manifest,
            [
                "comparison_id",
                "chop_length_bp",
                "filter_mode",
                "filter_label",
                "num_mappings",
                "scaffold_jump",
                "scoring",
                "source_chopped_paf",
                "filtered_paf",
                "filtered_paf_sha256",
                "reused_from",
                "status",
            ],
        )


if __name__ == "__main__":
    main()
