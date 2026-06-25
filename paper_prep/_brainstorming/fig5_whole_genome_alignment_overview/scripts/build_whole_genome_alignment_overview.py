#!/usr/bin/env python3
"""Build lightweight whole-genome Fig5 alignment overview tables.

This script consumes only already-produced strict/filtered upstream artifacts:

* strict untangle whole-genome overview segments
* SweepGA/FastGA f16 and f32 query-grid chopped + SweepGA 1:1 ANI-filtered PAFs
* wfmash -p95 updated-bin query-grid chopped + SweepGA 1:1 ANI-filtered PAFs

It does not run raw whole-genome alignment or many-to-many filtering.
"""

from __future__ import annotations

import csv
import gzip
import math
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
OUT_DIR = Path(__file__).resolve().parents[1]
BIN_BP = 1_000_000

UNTANGLE_SEGMENTS = ROOT / "paper_prep/_brainstorming/fig5_untangle_whole_genome_overview/untangle_whole_genome_segments.tsv"
F16_MANIFEST = ROOT / "paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/summaries/query_grid_chop_filter_manifest.tsv"
F32_MANIFEST = ROOT / "paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency32/summaries/query_grid_chop_filter_manifest.tsv"
WFMASH_MANIFEST = ROOT / "paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/query_grid_filter_manifest.tsv"
PANEL_WINDOWS = ROOT / "paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_query_grid_chop_filter_panels/config/panel_windows.tsv"

CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CHR_RANK = {chrom: i for i, chrom in enumerate(CHR_ORDER)}
TARGET_RE = re.compile(r"chr(?:[0-9]{1,2}|X|Y)")


def open_text(path: str | Path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def parse_chrom(name: str) -> str:
    hits = TARGET_RE.findall(name)
    if not hits:
        return "unknown"
    return hits[-1]


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    return (CHR_RANK.get(chrom, 999), chrom)


def arm_proxy(chrom: str, start: int, end: int, length: int) -> str:
    """Assign p/q by midpoint when cytoband centromeres are not available.

    Upstream PAF coordinates are native sample assembly coordinates, not CHM13.
    The overview needs arm-level aggregation but exact centromere coordinates
    are unavailable for these haplotypes, so this proxy is documented in README.
    """

    if chrom == "unknown" or length <= 0:
        return "unknown"
    mid = (start + end) / 2
    suffix = "p" if mid < (length / 2) else "q"
    return f"{chrom}{suffix}"


def identity_from_paf(parts: list[str]) -> float:
    matches = int(parts[9])
    aln_len = int(parts[10])
    if aln_len <= 0:
        return 0.0
    return matches / aln_len


def iter_bin_overlaps(start: int, end: int, bin_bp: int = BIN_BP):
    if end <= start:
        return
    first = start // bin_bp
    last = (end - 1) // bin_bp
    for idx in range(first, last + 1):
        b0 = idx * bin_bp
        b1 = b0 + bin_bp
        ov = max(0, min(end, b1) - max(start, b0))
        if ov > 0:
            yield idx, b0, b1, ov


def add_support(support, matrix, key, qchrom, qlen, qstart, qend, tchrom, tlen, tstart, tend, ident):
    qarm = arm_proxy(qchrom, qstart, qend, qlen)
    tarm = arm_proxy(tchrom, tstart, tend, tlen)
    target_family = tchrom if qchrom != tchrom else "same_chromosome"
    for bin_idx, b0, b1, overlap in iter_bin_overlaps(qstart, qend):
        skey = key + (qchrom, bin_idx)
        support[skey]["query_bin_start"] = b0
        support[skey]["query_bin_end"] = b1
        support[skey]["query_chrom_length"] = max(support[skey].get("query_chrom_length", 0), qlen)
        support[skey]["total_support_bp"] += overlap
        target_key = (target_family, tchrom, tarm)
        support[skey]["targets"][target_key] += overlap
        support[skey]["target_identity_weight"][target_key] += overlap * ident
    mkey = key + (qarm, tarm, qchrom, tchrom)
    matrix[mkey]["support_bp"] += max(0, qend - qstart)
    matrix[mkey]["weighted_identity_sum"] += max(0, qend - qstart) * ident
    matrix[mkey]["row_count"] += 1


def build_method_manifest() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    rows.append(
        {
            "method_id": "untangle_strict",
            "method_family": "untangle",
            "method_label": "Untangle strict primary path",
            "comparison_id": "all_three_transmissions",
            "comparison_label": "strict untangle overview",
            "fastga_frequency": "NA",
            "chop_length_bp": "NA",
            "filter_id": "strict_primary_path_only",
            "chunk_mode": "terminal-window",
            "scoring": "strict corrected untangle geometry",
            "input_path": str(UNTANGLE_SEGMENTS),
            "input_exists": UNTANGLE_SEGMENTS.exists(),
            "input_rows": "",
            "status": "OK" if UNTANGLE_SEGMENTS.exists() else "MISSING",
        }
    )
    for manifest_path, family, label, freq in [
        (F16_MANIFEST, "sweepga_fastga_f16", "SweepGA/FastGA f16", "16"),
        (F32_MANIFEST, "sweepga_fastga_f32", "SweepGA/FastGA f32", "32"),
    ]:
        if not manifest_path.exists():
            continue
        for row in read_tsv(manifest_path):
            if row.get("status") != "OK":
                continue
            rows.append(
                {
                    "method_id": f"{family}_l{row['chop_length_bp']}_{row['comparison_id']}",
                    "method_family": family,
                    "method_label": label,
                    "comparison_id": row["comparison_id"],
                    "comparison_label": row["comparison_id"].replace("_joint", ""),
                    "fastga_frequency": row.get("fastga_frequency", freq),
                    "chop_length_bp": row["chop_length_bp"],
                    "filter_id": row["filter_id"],
                    "chunk_mode": row["chunk_mode"],
                    "scoring": row.get("scoring", "ani"),
                    "input_path": row["filtered_paf"],
                    "input_exists": os.path.exists(row["filtered_paf"]),
                    "input_rows": "",
                    "status": "OK" if os.path.exists(row["filtered_paf"]) else "MISSING",
                }
            )
    if WFMASH_MANIFEST.exists():
        for row in read_tsv(WFMASH_MANIFEST):
            if row.get("status") != "OK":
                continue
            rows.append(
                {
                    "method_id": f"wfmash_p95_l{row['chop_length_bp']}_{row['comparison_id']}",
                    "method_family": "wfmash_p95_updated_bin",
                    "method_label": "wfmash -p95 updated bin",
                    "comparison_id": row["comparison_id"],
                    "comparison_label": row["comparison_id"].replace("_joint", ""),
                    "fastga_frequency": "NA",
                    "chop_length_bp": row["chop_length_bp"],
                    "filter_id": row["filter_id"],
                    "chunk_mode": row["chunk_mode"],
                    "scoring": row.get("scoring", "ani"),
                    "input_path": row["filtered_paf"],
                    "input_exists": os.path.exists(row["filtered_paf"]),
                    "input_rows": row.get("filtered_row_count", ""),
                    "status": "OK" if os.path.exists(row["filtered_paf"]) else "MISSING",
                }
            )
    return rows


def aggregate_untangle(manifest_row, support, matrix, q_lengths) -> int:
    rows = read_tsv(UNTANGLE_SEGMENTS)
    key = (
        manifest_row["method_id"],
        manifest_row["method_family"],
        manifest_row["method_label"],
        manifest_row["comparison_id"],
        manifest_row["comparison_label"],
        str(manifest_row["fastga_frequency"]),
        str(manifest_row["chop_length_bp"]),
        manifest_row["filter_id"],
        manifest_row["chunk_mode"],
        manifest_row["scoring"],
    )
    count = 0
    for row in rows:
        qchrom = row["query_chrom"]
        tchrom = row["target_chrom"]
        qstart = int(float(row["native_query_start_0based"]))
        qend = int(float(row["native_query_end_0based_exclusive"]))
        qlen = max(int(float(row["query_source_end_0based_exclusive"])), qend)
        tstart = 0
        tend = int(float(row["segment_length_bp"]))
        tlen = max(tend, int(float(row["query_length"])))
        ident = float(row["identity"]) / 100 if float(row["identity"]) > 1.5 else float(row["identity"])
        q_lengths[(key, qchrom)] = max(q_lengths.get((key, qchrom), 0), qlen)
        add_support(support, matrix, key, qchrom, qlen, qstart, qend, tchrom, tlen, tstart, tend, ident)
        count += 1
    return count


def aggregate_paf(manifest_row, support, matrix, q_lengths) -> int:
    key = (
        manifest_row["method_id"],
        manifest_row["method_family"],
        manifest_row["method_label"],
        manifest_row["comparison_id"],
        manifest_row["comparison_label"],
        str(manifest_row["fastga_frequency"]),
        str(manifest_row["chop_length_bp"]),
        manifest_row["filter_id"],
        manifest_row["chunk_mode"],
        manifest_row["scoring"],
    )
    count = 0
    with open_text(manifest_row["input_path"]) as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qname = parts[0]
            qlen = int(parts[1])
            qstart = int(parts[2])
            qend = int(parts[3])
            tname = parts[5]
            tlen = int(parts[6])
            tstart = int(parts[7])
            tend = int(parts[8])
            ident = identity_from_paf(parts)
            qchrom = parse_chrom(qname)
            tchrom = parse_chrom(tname)
            q_lengths[(key, qchrom)] = max(q_lengths.get((key, qchrom), 0), qlen)
            add_support(support, matrix, key, qchrom, qlen, qstart, qend, tchrom, tlen, tstart, tend, ident)
            count += 1
    return count


def summarize_support(manifest_rows, support, q_lengths) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    method_keys = []
    for row in manifest_rows:
        if row["status"] != "OK":
            continue
        method_keys.append(
            (
                row["method_id"],
                row["method_family"],
                row["method_label"],
                row["comparison_id"],
                row["comparison_label"],
                str(row["fastga_frequency"]),
                str(row["chop_length_bp"]),
                row["filter_id"],
                row["chunk_mode"],
                row["scoring"],
            )
        )
    for key in method_keys:
        chroms = sorted({chrom for (mkey, chrom), length in q_lengths.items() if mkey == key and length > 0}, key=chrom_sort_key)
        for qchrom in chroms:
            qlen = q_lengths[(key, qchrom)]
            n_bins = math.ceil(qlen / BIN_BP)
            for bin_idx in range(n_bins):
                skey = key + (qchrom, bin_idx)
                state = support.get(skey)
                b0 = bin_idx * BIN_BP
                b1 = min((bin_idx + 1) * BIN_BP, qlen)
                if state:
                    total = state["total_support_bp"]
                    winner, winner_bp = max(state["targets"].items(), key=lambda item: (item[1], item[0][0]))
                    target_family, target_chrom, target_arm = winner
                    identity_weight = state["target_identity_weight"][winner]
                    mean_identity = identity_weight / winner_bp if winner_bp else 0
                    inter_targets = [(target, bp) for target, bp in state["targets"].items() if target[0] not in ("same_chromosome", "no_support")]
                    if inter_targets:
                        top_inter, top_inter_bp = max(inter_targets, key=lambda item: (item[1], item[0][1]))
                        top_inter_chrom = top_inter[1]
                        top_inter_arm = top_inter[2]
                    else:
                        top_inter_bp = 0
                        top_inter_chrom = "none"
                        top_inter_arm = "none"
                    support_fraction = min(1, total / max(1, b1 - b0))
                    no_support = "no"
                else:
                    total = 0
                    winner_bp = 0
                    target_family = "no_support"
                    target_chrom = "no_support"
                    target_arm = "no_support"
                    mean_identity = 0
                    top_inter_bp = 0
                    top_inter_chrom = "none"
                    top_inter_arm = "none"
                    support_fraction = 0
                    no_support = "yes"
                qarm = arm_proxy(qchrom, b0, b1, qlen)
                out.append(
                    {
                        "bin_bp": BIN_BP,
                        "method_id": key[0],
                        "method_family": key[1],
                        "method_label": key[2],
                        "comparison_id": key[3],
                        "comparison_label": key[4],
                        "fastga_frequency": key[5],
                        "chop_length_bp": key[6],
                        "filter_id": key[7],
                        "chunk_mode": key[8],
                        "scoring": key[9],
                        "query_chrom": qchrom,
                        "query_arm_proxy": qarm,
                        "query_chrom_length": qlen,
                        "query_bin_index": bin_idx,
                        "query_bin_start": b0,
                        "query_bin_end": b1,
                        "support_bp": int(total),
                        "support_fraction": f"{support_fraction:.6f}",
                        "winner_target_family": target_family,
                        "winner_target_chrom": target_chrom,
                        "winner_target_arm_proxy": target_arm,
                        "winner_support_bp": int(winner_bp),
                        "winner_mean_identity": f"{mean_identity:.6f}",
                        "top_interchrom_target_chrom": top_inter_chrom,
                        "top_interchrom_target_arm_proxy": top_inter_arm,
                        "top_interchrom_support_bp": int(top_inter_bp),
                        "no_support": no_support,
                    }
                )
    return out


def summarize_matrix(matrix) -> list[dict[str, object]]:
    rows = []
    for key, vals in matrix.items():
        (
            method_id,
            method_family,
            method_label,
            comparison_id,
            comparison_label,
            fastga_frequency,
            chop_length_bp,
            filter_id,
            chunk_mode,
            scoring,
            qarm,
            tarm,
            qchrom,
            tchrom,
        ) = key
        support_bp = vals["support_bp"]
        rows.append(
            {
                "method_id": method_id,
                "method_family": method_family,
                "method_label": method_label,
                "comparison_id": comparison_id,
                "comparison_label": comparison_label,
                "fastga_frequency": fastga_frequency,
                "chop_length_bp": chop_length_bp,
                "filter_id": filter_id,
                "chunk_mode": chunk_mode,
                "scoring": scoring,
                "query_arm_proxy": qarm,
                "target_arm_proxy": tarm,
                "query_chrom": qchrom,
                "target_chrom": tchrom,
                "support_bp": int(support_bp),
                "row_count": int(vals["row_count"]),
                "mean_identity": f"{(vals['weighted_identity_sum'] / support_bp) if support_bp else 0:.6f}",
                "interchromosomal": "yes" if qchrom != tchrom else "no",
            }
        )
    rows.sort(key=lambda r: (r["method_family"], int(r["chop_length_bp"]) if str(r["chop_length_bp"]).isdigit() else -1, r["comparison_id"], r["query_arm_proxy"], r["target_arm_proxy"]))
    return rows


def add_callouts(binned_rows: list[dict[str, object]]) -> None:
    if not PANEL_WINDOWS.exists():
        return
    windows = read_tsv(PANEL_WINDOWS)
    for row in binned_rows:
        row["callout_event_id"] = ""
        row["callout_label"] = ""
        row["callout_expected_target_chroms"] = ""
        row["callout_window_start"] = ""
        row["callout_window_end"] = ""
        if row["no_support"] == "yes":
            continue
        for win in windows:
            if row["comparison_id"] != win["comparison_id"]:
                continue
            if row["query_chrom"] != parse_chrom(win["query_name"]):
                continue
            w0 = int(win["query_start"])
            w1 = int(win["query_end"])
            if int(row["query_bin_start"]) < w1 and int(row["query_bin_end"]) > w0:
                row["callout_event_id"] = win["event_id"]
                row["callout_label"] = win["panel_label"]
                row["callout_expected_target_chroms"] = win["expected_target_chroms"]
                row["callout_window_start"] = w0
                row["callout_window_end"] = w1
                break


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest_rows = build_method_manifest()
    support = defaultdict(lambda: {"total_support_bp": 0, "targets": Counter(), "target_identity_weight": Counter()})
    matrix = defaultdict(lambda: {"support_bp": 0, "weighted_identity_sum": 0.0, "row_count": 0})
    q_lengths = {}

    for row in manifest_rows:
        if row["status"] != "OK":
            continue
        if row["method_family"] == "untangle":
            input_rows = aggregate_untangle(row, support, matrix, q_lengths)
        else:
            input_rows = aggregate_paf(row, support, matrix, q_lengths)
        row["input_rows"] = input_rows
        print(f"aggregated {row['method_id']}: {input_rows} rows", file=sys.stderr, flush=True)

    binned_rows = summarize_support(manifest_rows, support, q_lengths)
    add_callouts(binned_rows)
    matrix_rows = summarize_matrix(matrix)

    method_fields = [
        "method_id",
        "method_family",
        "method_label",
        "comparison_id",
        "comparison_label",
        "fastga_frequency",
        "chop_length_bp",
        "filter_id",
        "chunk_mode",
        "scoring",
        "input_path",
        "input_exists",
        "input_rows",
        "status",
    ]
    binned_fields = [
        "bin_bp",
        "method_id",
        "method_family",
        "method_label",
        "comparison_id",
        "comparison_label",
        "fastga_frequency",
        "chop_length_bp",
        "filter_id",
        "chunk_mode",
        "scoring",
        "query_chrom",
        "query_arm_proxy",
        "query_chrom_length",
        "query_bin_index",
        "query_bin_start",
        "query_bin_end",
        "support_bp",
        "support_fraction",
        "winner_target_family",
        "winner_target_chrom",
        "winner_target_arm_proxy",
        "winner_support_bp",
        "winner_mean_identity",
        "top_interchrom_target_chrom",
        "top_interchrom_target_arm_proxy",
        "top_interchrom_support_bp",
        "no_support",
        "callout_event_id",
        "callout_label",
        "callout_expected_target_chroms",
        "callout_window_start",
        "callout_window_end",
    ]
    matrix_fields = [
        "method_id",
        "method_family",
        "method_label",
        "comparison_id",
        "comparison_label",
        "fastga_frequency",
        "chop_length_bp",
        "filter_id",
        "chunk_mode",
        "scoring",
        "query_arm_proxy",
        "target_arm_proxy",
        "query_chrom",
        "target_chrom",
        "support_bp",
        "row_count",
        "mean_identity",
        "interchromosomal",
    ]
    write_tsv(OUT_DIR / "whole_genome_method_manifest.tsv", method_fields, manifest_rows)
    write_tsv(OUT_DIR / "whole_genome_binned_support.tsv", binned_fields, binned_rows)
    write_tsv(OUT_DIR / "whole_genome_support_matrix.tsv", matrix_fields, matrix_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
