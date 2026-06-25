#!/usr/bin/env python3
"""Summarize existing WFMASH/SweepGA PAFs as IMPG-like query-space support.

The full scan is intended for Slurm. Work is sharded by PAF/evidence layer and
helper counts are derived from SLURM_CPUS_PER_TASK. Each worker streams a PAF,
assigns aligned query bases to configurable bins, and emits target support.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[4]
SCAN_DIR = Path(__file__).resolve().parents[1]
ACROS = {"chr13", "chr14", "chr15", "chr21", "chr22"}


@dataclass(frozen=True)
class InputRow:
    method_id: str
    evidence_layer: str
    comparison_id: str
    bin_size_bp: int
    paf_path: str
    row_count_manifest: str


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str], gzip_output: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if gzip_output else open
    mode = "wt" if gzip_output else "w"
    with opener(path, mode, newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def parse_chrom(name: str) -> str:
    match = re.search(r"(?:^|[#_/])chr([0-9]+|X|Y|M)(?:[._#:/-]|$)", name)
    if match:
        return f"chr{match.group(1)}"
    match = re.search(r"chr([0-9]+|X|Y|M)", name)
    if match:
        return f"chr{match.group(1)}"
    return name


def tag_value(fields: list[str], tag: str) -> str:
    prefix = f"{tag}:"
    for field in fields[12:]:
        if field.startswith(prefix):
            return field.split(":", 2)[-1]
    return ""


def chrom_arm(chrom: str, start: int, end: int, chrom_sizes: dict[str, int], subtelomere_bp: int) -> str:
    size = chrom_sizes.get(chrom)
    if size is None:
        return f"{chrom}_unknown"
    if end <= subtelomere_bp:
        return f"{chrom}p"
    if start >= max(0, size - subtelomere_bp):
        return f"{chrom}q"
    return f"{chrom}_internal"


def arm_side(arm: str) -> str:
    if arm.endswith("p"):
        return "p"
    if arm.endswith("q"):
        return "q"
    return "internal"


def identity_from_paf(fields: list[str]) -> float:
    try:
        matches = float(fields[9])
        block = float(fields[10])
        if block > 0:
            return matches / block
    except (IndexError, ValueError):
        pass
    dv = tag_value(fields, "dv")
    if dv:
        try:
            return max(0.0, min(1.0, 1.0 - float(dv)))
        except ValueError:
            return math.nan
    return math.nan


def open_paf(path: str, pigz_threads: int):
    if path.endswith(".gz") and shutil.which("pigz"):
        proc = subprocess.Popen(["pigz", "-dc", "-p", str(max(1, pigz_threads)), path], stdout=subprocess.PIPE, text=True)
        assert proc.stdout is not None
        return proc, proc.stdout
    if path.endswith(".gz"):
        return None, gzip.open(path, "rt")
    return None, open(path)


def iter_bin_overlaps(q_start: int, q_end: int, bin_size: int):
    first = q_start // bin_size
    last = (q_end - 1) // bin_size
    for bin_index in range(first, last + 1):
        b_start = bin_index * bin_size
        b_end = b_start + bin_size
        overlap = max(0, min(q_end, b_end) - max(q_start, b_start))
        if overlap:
            yield bin_index, b_start, b_end, overlap


def classify_support(query_chrom: str, target_chrom: str) -> str:
    if query_chrom == target_chrom:
        return "same_chromosome"
    return "interchromosomal"


def process_one(args: tuple[InputRow, dict[str, int], int, int, Path]) -> dict[str, object]:
    row, chrom_sizes, subtelomere_bp, pigz_threads, tmp_dir = args
    started = time.time()
    bin_size = row.bin_size_bp
    streaming_reduce = row.method_id == "sweepga_fastga_f32" and row.evidence_layer == "raw_many_many"
    support: dict[tuple[str, ...], dict[str, float]] = defaultdict(lambda: {"paf_rows": 0, "aligned_bp": 0, "identity_bp": 0.0})
    total_lines = 0
    stem = f"{row.method_id}.{row.evidence_layer}.{row.comparison_id}.bin_target_support.tsv"
    out = tmp_dir / stem
    obs = tmp_dir / stem.replace(".tsv", ".observations.tsv")
    obs_sorted = tmp_dir / stem.replace(".tsv", ".observations.sorted.tsv")
    obs_handle = obs.open("w") if streaming_reduce else None
    proc, handle = open_paf(row.paf_path, pigz_threads)
    try:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            total_lines += 1
            query_name = fields[0]
            target_name = fields[5]
            q_start = int(fields[2])
            q_end = int(fields[3])
            t_start = int(fields[7])
            t_end = int(fields[8])
            query_chrom = parse_chrom(query_name)
            target_chrom = parse_chrom(target_name)
            query_arm = chrom_arm(query_chrom, q_start, q_end, chrom_sizes, subtelomere_bp)
            target_arm = chrom_arm(target_chrom, t_start, t_end, chrom_sizes, subtelomere_bp)
            ident = identity_from_paf(fields)
            for bin_index, b_start, b_end, ov in iter_bin_overlaps(q_start, q_end, bin_size):
                key = (
                    row.method_id,
                    row.evidence_layer,
                    row.comparison_id,
                    query_name,
                    query_chrom,
                    query_arm,
                    str(bin_index),
                    str(b_start),
                    str(b_end),
                    target_chrom,
                    target_arm,
                    classify_support(query_chrom, target_chrom),
                )
                if obs_handle is not None:
                    ident_bp = ident * ov if not math.isnan(ident) else 0.0
                    obs_handle.write("\t".join((*key, str(ov), f"{ident_bp:.6f}", "1")) + "\n")
                    continue
                acc = support[key]
                acc["paf_rows"] += 1
                acc["aligned_bp"] += ov
                if not math.isnan(ident):
                    acc["identity_bp"] += ident * ov
    finally:
        handle.close()
        if proc is not None:
            rc = proc.wait()
            if rc != 0:
                raise RuntimeError(f"pigz failed with exit {rc} for {row.paf_path}")
        if obs_handle is not None:
            obs_handle.close()

    fields_out = [
        "method_id", "evidence_layer", "comparison_id", "query_name", "query_chrom", "query_arm",
        "bin_size_bp", "bin_index", "bin_start", "bin_end", "target_chrom", "target_arm",
        "support_class", "paf_rows", "aligned_bp_sum", "bin_fraction_aligned",
        "mean_identity_weighted", "mean_match_distance", "paf_path",
    ]
    output_rows = 0
    if streaming_reduce:
        sort_threads = max(1, pigz_threads)
        subprocess.run(
            [
                "sort",
                "--parallel",
                str(sort_threads),
                "-T",
                str(tmp_dir),
                "-o",
                str(obs_sorted),
                str(obs),
            ],
            check=True,
        )
        with obs_sorted.open() as sorted_in, out.open("w", newline="") as handle_out:
            writer = csv.DictWriter(handle_out, delimiter="\t", fieldnames=fields_out, lineterminator="\n")
            writer.writeheader()
            last_key: tuple[str, ...] | None = None
            rows_count = 0
            aligned_sum = 0
            identity_sum = 0.0

            def flush() -> None:
                nonlocal output_rows, rows_count, aligned_sum, identity_sum, last_key
                if last_key is None:
                    return
                (
                    method_id, evidence_layer, comparison_id, query_name, query_chrom, query_arm, bin_index,
                    b_start, b_end, target_chrom, target_arm, support_class,
                ) = last_key
                mean_id = identity_sum / aligned_sum if aligned_sum else math.nan
                writer.writerow({
                    "method_id": method_id,
                    "evidence_layer": evidence_layer,
                    "comparison_id": comparison_id,
                    "query_name": query_name,
                    "query_chrom": query_chrom,
                    "query_arm": query_arm,
                    "bin_size_bp": bin_size,
                    "bin_index": bin_index,
                    "bin_start": b_start,
                    "bin_end": b_end,
                    "target_chrom": target_chrom,
                    "target_arm": target_arm,
                    "support_class": support_class,
                    "paf_rows": rows_count,
                    "aligned_bp_sum": aligned_sum,
                    "bin_fraction_aligned": round(aligned_sum / bin_size, 6),
                    "mean_identity_weighted": round(mean_id, 6) if not math.isnan(mean_id) else "",
                    "mean_match_distance": round(1.0 - mean_id, 6) if not math.isnan(mean_id) else "",
                    "paf_path": row.paf_path,
                })
                output_rows += 1

            for obs_line in sorted_in:
                parts = obs_line.rstrip("\n").split("\t")
                key = tuple(parts[:12])
                aligned = int(parts[12])
                ident_bp = float(parts[13])
                nrows = int(parts[14])
                if last_key is not None and key != last_key:
                    flush()
                    rows_count = 0
                    aligned_sum = 0
                    identity_sum = 0.0
                last_key = key
                rows_count += nrows
                aligned_sum += aligned
                identity_sum += ident_bp
            flush()
        if obs.exists():
            obs.unlink()
        if obs_sorted.exists():
            obs_sorted.unlink()
        support_rows = output_rows
    else:
        with out.open("w", newline="") as handle_out:
            writer = csv.DictWriter(handle_out, delimiter="\t", fieldnames=fields_out, lineterminator="\n")
            writer.writeheader()
            for key, acc in sorted(support.items()):
                (
                    method_id, evidence_layer, comparison_id, query_name, query_chrom, query_arm, bin_index,
                    b_start, b_end, target_chrom, target_arm, support_class,
                ) = key
                aligned = int(acc["aligned_bp"])
                mean_id = acc["identity_bp"] / aligned if aligned else math.nan
                writer.writerow({
                    "method_id": method_id,
                    "evidence_layer": evidence_layer,
                    "comparison_id": comparison_id,
                    "query_name": query_name,
                    "query_chrom": query_chrom,
                    "query_arm": query_arm,
                    "bin_size_bp": bin_size,
                    "bin_index": bin_index,
                    "bin_start": b_start,
                    "bin_end": b_end,
                    "target_chrom": target_chrom,
                    "target_arm": target_arm,
                    "support_class": support_class,
                    "paf_rows": int(acc["paf_rows"]),
                    "aligned_bp_sum": aligned,
                    "bin_fraction_aligned": round(aligned / bin_size, 6),
                    "mean_identity_weighted": round(mean_id, 6) if not math.isnan(mean_id) else "",
                    "mean_match_distance": round(1.0 - mean_id, 6) if not math.isnan(mean_id) else "",
                    "paf_path": row.paf_path,
                })
        support_rows = len(support)
    return {
        "method_id": row.method_id,
        "evidence_layer": row.evidence_layer,
        "comparison_id": row.comparison_id,
        "paf_path": row.paf_path,
        "input_lines": total_lines,
        "bin_target_rows": support_rows,
        "worker_seconds": round(time.time() - started, 3),
        "pigz_threads": pigz_threads,
        "reducer": "stream_sort" if streaming_reduce else "in_memory",
        "tmp_output": str(out),
    }


def load_chrom_sizes(path: Path) -> dict[str, int]:
    sizes = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, size = line.rstrip("\n").split("\t")[:2]
            sizes[chrom] = int(size)
    return sizes


def combine_worker_outputs(worker_rows: list[dict[str, object]], out_gz: Path) -> None:
    fields = [
        "method_id", "evidence_layer", "comparison_id", "query_name", "query_chrom", "query_arm",
        "bin_size_bp", "bin_index", "bin_start", "bin_end", "target_chrom", "target_arm",
        "support_class", "paf_rows", "aligned_bp_sum", "bin_fraction_aligned",
        "mean_identity_weighted", "mean_match_distance", "paf_path",
    ]
    with gzip.open(out_gz, "wt", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in worker_rows:
            with open(row["tmp_output"], newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for rec in reader:
                    writer.writerow(rec)


def summarize_targets(bin_support_gz: Path, bin_best_gz: Path, totals_tsv: Path, focal_tsv: Path, aggregate_tsv: Path, aggregates: list[int]) -> None:
    totals: dict[tuple[str, ...], dict[str, float]] = defaultdict(lambda: {"bins": 0, "paf_rows": 0, "aligned_bp": 0, "identity_bp": 0.0})
    aggregate_acc: dict[tuple[str, ...], dict[str, float]] = defaultdict(lambda: {"bins": 0, "aligned_bp": 0, "identity_bp": 0.0})
    focal: dict[tuple[str, ...], dict[str, float]] = defaultdict(lambda: {"bins": 0, "aligned_bp": 0, "identity_bp": 0.0})
    best_obs = bin_best_gz.with_suffix(".observations.tsv")
    best_sorted = bin_best_gz.with_suffix(".observations.sorted.tsv")

    with gzip.open(bin_support_gz, "rt", newline="") as handle, best_obs.open("w") as best_handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            aligned = int(row["aligned_bp_sum"])
            paf_rows = int(row["paf_rows"])
            mean_id = float(row["mean_identity_weighted"]) if row["mean_identity_weighted"] else math.nan
            identity_bp = mean_id * aligned if not math.isnan(mean_id) else 0.0
            total_key = (
                row["method_id"], row["evidence_layer"], row["comparison_id"],
                row["query_chrom"], arm_side(row["query_arm"]), row["target_chrom"], row["target_arm"], row["support_class"],
            )
            totals[total_key]["bins"] += 1
            totals[total_key]["paf_rows"] += paf_rows
            totals[total_key]["aligned_bp"] += aligned
            totals[total_key]["identity_bp"] += identity_bp

            best_handle.write("\t".join([
                row["method_id"], row["evidence_layer"], row["comparison_id"], row["query_name"],
                row["query_chrom"], row["query_arm"], row["bin_size_bp"], row["bin_index"],
                row["bin_start"], row["bin_end"], row["target_chrom"], row["target_arm"],
                row["support_class"], row["aligned_bp_sum"], row["mean_identity_weighted"] or "nan",
                row["paf_rows"],
            ]) + "\n")

            for agg in aggregates:
                b_start = int(row["bin_start"])
                agg_start = (b_start // agg) * agg
                agg_key = (
                    row["method_id"], row["evidence_layer"], row["comparison_id"], row["query_chrom"], row["query_arm"],
                    str(agg), str(agg_start), str(agg_start + agg), row["target_chrom"], row["target_arm"], row["support_class"],
                )
                aggregate_acc[agg_key]["bins"] += 1
                aggregate_acc[agg_key]["aligned_bp"] += aligned
                aggregate_acc[agg_key]["identity_bp"] += identity_bp

            qside = arm_side(row["query_arm"])
            tside = arm_side(row["target_arm"])
            qchrom = row["query_chrom"]
            tchrom = row["target_chrom"]
            if qchrom == "chr9" and qside == "q" and tchrom == "chr3" and tside == "q":
                region = "chr9q_to_chr3q"
            elif {qchrom, tchrom} == {"chrX", "chrY"}:
                region = "PAR_XY"
            elif qchrom in ACROS and tchrom in ACROS and qside == "p" and tside == "p" and qchrom != tchrom:
                region = "acrocentric_p_cross"
            else:
                region = ""
            if region:
                fkey = (region, row["method_id"], row["evidence_layer"], row["comparison_id"], qchrom, tchrom)
                focal[fkey]["bins"] += 1
                focal[fkey]["aligned_bp"] += aligned
                focal[fkey]["identity_bp"] += identity_bp

    total_fields = [
        "method_id", "evidence_layer", "comparison_id", "query_chrom", "query_arm_side",
        "target_chrom", "target_arm", "support_class", "support_bins", "paf_rows",
        "aligned_bp_sum", "mean_identity_weighted", "mean_match_distance",
    ]
    total_rows = []
    for key, acc in sorted(totals.items(), key=lambda kv: (-kv[1]["aligned_bp"], kv[0])):
        mean_id = acc["identity_bp"] / acc["aligned_bp"] if acc["aligned_bp"] else math.nan
        out_row = dict(zip(total_fields[:8], key))
        out_row.update({
            "support_bins": int(acc["bins"]),
            "paf_rows": int(acc["paf_rows"]),
            "aligned_bp_sum": int(acc["aligned_bp"]),
            "mean_identity_weighted": round(mean_id, 6) if not math.isnan(mean_id) else "",
            "mean_match_distance": round(1.0 - mean_id, 6) if not math.isnan(mean_id) else "",
        })
        total_rows.append(out_row)
    write_tsv(totals_tsv, total_rows, total_fields)

    subprocess.run(
        [
            "sort",
            "--parallel",
            str(max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))),
            "-T",
            str(bin_best_gz.parent),
            "-o",
            str(best_sorted),
            str(best_obs),
        ],
        check=True,
    )

    best_fields = [
        "method_id", "evidence_layer", "comparison_id", "query_name", "query_chrom", "query_arm",
        "bin_size_bp", "bin_index", "bin_start", "bin_end", "best_target_chrom", "best_target_arm",
        "best_support_class", "best_aligned_bp", "best_mean_identity_weighted",
        "all_target_aligned_bp", "all_target_paf_rows",
    ]
    with best_sorted.open() as sorted_in, gzip.open(bin_best_gz, "wt", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=best_fields, lineterminator="\n")
        writer.writeheader()
        last_key: tuple[str, ...] | None = None
        all_aligned = 0
        all_rows = 0
        best_target: list[str] | None = None
        best_aligned = -1

        def flush_best() -> None:
            nonlocal last_key, all_aligned, all_rows, best_target, best_aligned
            if last_key is None or best_target is None:
                return
            out_row = dict(zip(best_fields[:10], last_key))
            out_row.update({
                "best_target_chrom": best_target[0],
                "best_target_arm": best_target[1],
                "best_support_class": best_target[2],
                "best_aligned_bp": best_target[3],
                "best_mean_identity_weighted": best_target[4] if best_target[4] != "nan" else "",
                "all_target_aligned_bp": all_aligned,
                "all_target_paf_rows": all_rows,
            })
            writer.writerow(out_row)

        for line in sorted_in:
            parts = line.rstrip("\n").split("\t")
            key = tuple(parts[:10])
            aligned = int(parts[13])
            paf_rows = int(parts[15])
            if last_key is not None and key != last_key:
                flush_best()
                all_aligned = 0
                all_rows = 0
                best_target = None
                best_aligned = -1
            last_key = key
            all_aligned += aligned
            all_rows += paf_rows
            if aligned > best_aligned:
                best_aligned = aligned
                best_target = [parts[10], parts[11], parts[12], parts[13], parts[14]]
        flush_best()
    if best_obs.exists():
        best_obs.unlink()
    if best_sorted.exists():
        best_sorted.unlink()

    focal_fields = [
        "region", "method_id", "evidence_layer", "comparison_id", "query_chrom", "target_chrom",
        "support_bins", "aligned_bp_sum", "mean_identity_weighted", "mean_match_distance",
    ]
    focal_rows = []
    for key, acc in sorted(focal.items()):
        mean_id = acc["identity_bp"] / acc["aligned_bp"] if acc["aligned_bp"] else math.nan
        out_row = dict(zip(focal_fields[:6], key))
        out_row.update({
            "support_bins": int(acc["bins"]),
            "aligned_bp_sum": int(acc["aligned_bp"]),
            "mean_identity_weighted": round(mean_id, 6) if not math.isnan(mean_id) else "",
            "mean_match_distance": round(1.0 - mean_id, 6) if not math.isnan(mean_id) else "",
        })
        focal_rows.append(out_row)
    write_tsv(focal_tsv, focal_rows, focal_fields)

    aggregate_fields = [
        "method_id", "evidence_layer", "comparison_id", "query_chrom", "query_arm",
        "aggregate_bp", "aggregate_start", "aggregate_end", "target_chrom", "target_arm",
        "support_class", "source_2kb_target_bins", "aligned_bp_sum",
        "mean_identity_weighted", "mean_match_distance",
    ]
    aggregate_rows = []
    for key, acc in sorted(aggregate_acc.items()):
        mean_id = acc["identity_bp"] / acc["aligned_bp"] if acc["aligned_bp"] else math.nan
        out_row = dict(zip(aggregate_fields[:11], key))
        out_row.update({
            "source_2kb_target_bins": int(acc["bins"]),
            "aligned_bp_sum": int(acc["aligned_bp"]),
            "mean_identity_weighted": round(mean_id, 6) if not math.isnan(mean_id) else "",
            "mean_match_distance": round(1.0 - mean_id, 6) if not math.isnan(mean_id) else "",
        })
        aggregate_rows.append(out_row)
    write_tsv(aggregate_tsv, aggregate_rows, aggregate_fields)


def maybe_write_parquet(tsv_path: Path, parquet_path: Path) -> str:
    try:
        import pyarrow.csv as pacsv  # type: ignore
        import pyarrow.parquet as pq  # type: ignore
    except Exception as exc:
        return f"SKIP: pyarrow unavailable ({type(exc).__name__})"
    table = pacsv.read_csv(tsv_path, parse_options=pacsv.ParseOptions(delimiter="\t"))
    pq.write_table(table, parquet_path)
    return "OK"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=SCAN_DIR / "manifests/paf_inputs.tsv")
    parser.add_argument("--chrom-sizes", type=Path, default=ROOT / "data/chm13.chrom.sizes")
    parser.add_argument("--outdir", type=Path, default=SCAN_DIR / "summaries")
    parser.add_argument("--tmpdir", type=Path, default=SCAN_DIR / "summaries/tmp_worker_bin_support")
    parser.add_argument("--subtelomere-bp", type=int, default=500_000)
    parser.add_argument("--aggregate-bp", default="10000,50000")
    parser.add_argument("--include-evidence-layer", action="append", default=[])
    parser.add_argument("--reuse-worker-outputs", action="store_true")
    args = parser.parse_args()

    cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
    rows = []
    for row in read_tsv(args.manifest):
        if row["path_exists"] != "yes":
            continue
        if args.include_evidence_layer and row["evidence_layer"] not in args.include_evidence_layer:
            continue
        rows.append(InputRow(row["method_id"], row["evidence_layer"], row["comparison_id"], int(row["bin_size_bp"]), row["paf_path"], row["row_count_manifest"]))
    if not rows:
        raise SystemExit("No readable PAF inputs selected")

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.tmpdir.mkdir(parents=True, exist_ok=True)
    chrom_sizes = load_chrom_sizes(args.chrom_sizes)
    workers = min(len(rows), max(1, cpus))
    pigz_threads = max(1, cpus // workers)

    started = time.time()
    worker_summaries = []
    if args.reuse_worker_outputs:
        for row in rows:
            stem = f"{row.method_id}.{row.evidence_layer}.{row.comparison_id}.bin_target_support.tsv"
            tmp_output = args.tmpdir / stem
            if not tmp_output.exists():
                raise SystemExit(f"Missing worker output for reuse: {tmp_output}")
            line_count = int(subprocess.check_output(["wc", "-l", str(tmp_output)], text=True).split()[0])
            worker_summaries.append({
                "method_id": row.method_id,
                "evidence_layer": row.evidence_layer,
                "comparison_id": row.comparison_id,
                "paf_path": row.paf_path,
                "input_lines": "",
                "bin_target_rows": max(0, line_count - 1),
                "worker_seconds": "",
                "pigz_threads": "",
                "reducer": "reused_worker_output",
                "tmp_output": str(tmp_output),
            })
    else:
        worker_args = [(row, chrom_sizes, args.subtelomere_bp, pigz_threads, args.tmpdir) for row in rows]
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(process_one, wa) for wa in worker_args]
            for future in as_completed(futures):
                worker_summaries.append(future.result())

    bin_support = args.outdir / "bin_target_support.tsv.gz"
    combine_worker_outputs(worker_summaries, bin_support)
    summarize_targets(
        bin_support,
        args.outdir / "bin_best_support.tsv.gz",
        args.outdir / "target_support_totals.tsv",
        args.outdir / "focal_region_summary.tsv",
        args.outdir / "aggregate_target_support.tsv",
        [int(x) for x in args.aggregate_bp.split(",") if x],
    )

    resource_rows = [{
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", "not_slurm"),
        "hostname": os.uname().nodename,
        "slurm_cpus_per_task": cpus,
        "input_paf_count": len(rows),
        "process_workers": workers,
        "pigz_threads_per_worker": pigz_threads,
        "accounted_helper_threads": workers * pigz_threads,
        "wall_seconds": round(time.time() - started, 3),
        "python_executable": sys.executable,
    }]
    resource_rows.extend(worker_summaries)
    write_tsv(
        args.outdir / "resource_usage.tsv",
        resource_rows,
        [
            "slurm_job_id", "hostname", "slurm_cpus_per_task", "input_paf_count", "process_workers",
            "pigz_threads_per_worker", "accounted_helper_threads", "wall_seconds", "python_executable",
            "method_id", "evidence_layer", "comparison_id", "paf_path", "input_lines",
            "bin_target_rows", "worker_seconds", "pigz_threads", "reducer", "tmp_output",
        ],
    )

    parquet_rows = []
    for name in ["target_support_totals.tsv", "focal_region_summary.tsv", "aggregate_target_support.tsv", "resource_usage.tsv"]:
        status = maybe_write_parquet(args.outdir / name, args.outdir / name.replace(".tsv", ".parquet"))
        parquet_rows.append({"tsv": name, "parquet": name.replace(".tsv", ".parquet"), "status": status})
    write_tsv(args.outdir / "parquet_status.tsv", parquet_rows, ["tsv", "parquet", "status"])


if __name__ == "__main__":
    main()
