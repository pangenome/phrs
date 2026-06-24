#!/usr/bin/env python3
import collections
import csv
import gzip
import hashlib
import os
import re
import subprocess
import sys

from common import CANDIDATE_WINDOWS, PACKAGE_DIR, read_tsv, write_tsv


METHOD = "wfmash_p95_updated_bin"
STATUS_DIR = os.path.join(PACKAGE_DIR, "summaries", "query_grid_filter_status")
MANIFEST = os.path.join(PACKAGE_DIR, "summaries", "query_grid_filter_manifest.tsv")
CANDIDATE_SUMMARY = os.path.join(PACKAGE_DIR, "summaries", "query_grid_filter_candidate_window_support.tsv")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_gzip(path):
    subprocess.check_call(["pigz", "-t", path])


def row_count(path):
    count = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                count += 1
    return count


def paf_records(path):
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 12:
                yield fields


def chrom_from_name(name):
    match = re.search(r"#chr([0-9XYM]+)", name)
    if match:
        return "chr" + match.group(1)
    match = re.search(r"chr([0-9XYM]+)", name)
    if match:
        return "chr" + match.group(1)
    return name


def overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def merge_coverage(intervals):
    if not intervals:
        return 0
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def read_status_rows():
    rows = []
    if not os.path.isdir(STATUS_DIR):
        return rows
    for name in sorted(os.listdir(STATUS_DIR)):
        if not name.endswith(".tsv"):
            continue
        rows.extend(read_tsv(os.path.join(STATUS_DIR, name)))
    return rows


def write_sidecar(path, digest):
    with open(path + ".sha256", "w") as fh:
        fh.write("%s  %s\n" % (digest, path))


def build_manifest(status_rows):
    manifest_rows = []
    failures = []
    for row in status_rows:
        raw = row["raw_paf"]
        chopped = row["chopped_paf"]
        filtered = row["filtered_paf"]
        status = row.get("status", "")
        try:
            validate_gzip(chopped)
            validate_gzip(filtered)
            chopped_sha = sha256(chopped)
            filtered_sha = sha256(filtered)
            write_sidecar(chopped, chopped_sha)
            write_sidecar(filtered, filtered_sha)
            raw_rows = row_count(raw)
            chopped_rows = row_count(chopped)
            filtered_rows = row_count(filtered)
        except (OSError, subprocess.CalledProcessError) as exc:
            failures.append("%s length %s: %s" % (row["comparison_id"], row["chop_length_bp"], exc))
            chopped_sha = ""
            filtered_sha = ""
            raw_rows = chopped_rows = filtered_rows = ""
            status = "MISSING_OR_INVALID"
        manifest_rows.append({
            "method": METHOD,
            "comparison_id": row["comparison_id"],
            "chop_length_bp": row["chop_length_bp"],
            "filter_id": row["filter_id"],
            "chunk_mode": "query-grid",
            "raw_paf": raw,
            "chopped_paf": chopped,
            "filtered_paf": filtered,
            "raw_row_count": str(raw_rows),
            "chopped_row_count": str(chopped_rows),
            "filtered_row_count": str(filtered_rows),
            "chopped_sha256": chopped_sha,
            "filtered_sha256": filtered_sha,
            "chop_summary": row["chop_summary"],
            "num_mappings": "1:1",
            "scaffold_jump": "0",
            "scoring": "ani",
            "filter_overlap": "0",
            "cg_verification": row["cg_verification"],
            "gzip_validation": "pigz -t",
            "chop_command": row["chop_command"],
            "filter_command": row["filter_command"],
            "pafchop_bin": row["pafchop_bin"],
            "pafchop_sha256": row["pafchop_sha256"],
            "sweepga_bin": row["sweepga_bin"],
            "sweepga_sha256": row["sweepga_sha256"],
            "pigz_bin": row["pigz_bin"],
            "pigz_sha256": row["pigz_sha256"],
            "started_utc": row["started_utc"],
            "finished_utc": row["finished_utc"],
            "status": status,
        })
    fields = [
        "method",
        "comparison_id",
        "chop_length_bp",
        "filter_id",
        "chunk_mode",
        "raw_paf",
        "chopped_paf",
        "filtered_paf",
        "raw_row_count",
        "chopped_row_count",
        "filtered_row_count",
        "chopped_sha256",
        "filtered_sha256",
        "chop_summary",
        "num_mappings",
        "scaffold_jump",
        "scoring",
        "filter_overlap",
        "cg_verification",
        "gzip_validation",
        "chop_command",
        "filter_command",
        "pafchop_bin",
        "pafchop_sha256",
        "sweepga_bin",
        "sweepga_sha256",
        "pigz_bin",
        "pigz_sha256",
        "started_utc",
        "finished_utc",
        "status",
    ]
    write_tsv(MANIFEST, manifest_rows, fields)
    return manifest_rows, failures


def summarize_candidate_windows(manifest_rows):
    windows = read_tsv(CANDIDATE_WINDOWS)
    rows = []
    for manifest in manifest_rows:
        filtered = manifest["filtered_paf"]
        if manifest["status"] != "OK":
            continue
        for window in windows:
            if manifest["comparison_id"] != window["comparison_id"]:
                continue
            qname = window["query_name"]
            qstart = int(window["query_start"])
            qend = int(window["query_end"])
            by_chrom = collections.defaultdict(lambda: {"rows": 0, "bp_sum": 0, "intervals": []})
            for fields in paf_records(filtered):
                if fields[0] != qname:
                    continue
                row_qstart = int(fields[2])
                row_qend = int(fields[3])
                ov = overlap(row_qstart, row_qend, qstart, qend)
                if ov <= 0:
                    continue
                chrom = chrom_from_name(fields[5])
                by_chrom[chrom]["rows"] += 1
                by_chrom[chrom]["bp_sum"] += ov
                by_chrom[chrom]["intervals"].append((max(row_qstart, qstart), min(row_qend, qend)))
            if not by_chrom:
                by_chrom["NO_OVERLAP"]
            if window["expected_target_chrom"] not in by_chrom:
                by_chrom[window["expected_target_chrom"]]
            for chrom, entry in sorted(by_chrom.items()):
                rows.append({
                    "method": METHOD,
                    "event_id": window["event_id"],
                    "comparison_id": window["comparison_id"],
                    "chop_length_bp": manifest["chop_length_bp"],
                    "filter_id": manifest["filter_id"],
                    "query_name": qname,
                    "query_start": str(qstart),
                    "query_end": str(qend),
                    "target_chrom": chrom,
                    "expected_target_chrom": window["expected_target_chrom"],
                    "retained_rows_overlapping_window": str(entry["rows"]),
                    "query_overlap_bp_sum": str(entry["bp_sum"]),
                    "query_covered_bp_union": str(merge_coverage(entry["intervals"])),
                    "chr3_support": "yes" if chrom == "chr3" and entry["rows"] > 0 else "no",
                    "filtered_paf": filtered,
                })
    fields = [
        "method",
        "event_id",
        "comparison_id",
        "chop_length_bp",
        "filter_id",
        "query_name",
        "query_start",
        "query_end",
        "target_chrom",
        "expected_target_chrom",
        "retained_rows_overlapping_window",
        "query_overlap_bp_sum",
        "query_covered_bp_union",
        "chr3_support",
        "filtered_paf",
    ]
    write_tsv(CANDIDATE_SUMMARY, rows, fields)
    return rows


def main():
    status_rows = read_status_rows()
    if not status_rows:
        print("no query-grid filter status rows found in %s" % STATUS_DIR, file=sys.stderr)
        return 1
    manifest_rows, failures = build_manifest(status_rows)
    candidate_rows = summarize_candidate_windows(manifest_rows)
    chr3_rows = [
        row for row in candidate_rows
        if row["target_chrom"] == "chr3" and row["comparison_id"] in {"PAN027pat_vs_PAN011_joint", "PAN028mat_vs_PAN027_joint"}
    ]
    if len(chr3_rows) != 6:
        failures.append("expected 6 chr3 candidate-window rows after filtering, saw %d" % len(chr3_rows))
    if failures:
        for failure in failures:
            print(failure, file=sys.stderr)
        return 1
    print("Validated %d wfmash query-grid filtered outputs" % len(manifest_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
