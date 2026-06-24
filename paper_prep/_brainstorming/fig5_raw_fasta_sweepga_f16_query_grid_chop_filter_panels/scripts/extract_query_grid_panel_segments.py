#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import os
import re
from collections import defaultdict


LENGTHS = (10000, 5000, 2000)
FILTER_ID = "one_to_one_ani_o0"


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def chrom_from_name(name):
    matches = re.findall(r"chr(?:[0-9]+|X|Y|M)", name)
    return matches[-1] if matches else name


def overlap_bp(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def union_bp(intervals):
    if not intervals:
        return 0
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return sum(end - start for start, end in merged)


def first_existing_path(path, repo_root):
    if os.path.exists(path):
        return path
    if path.startswith("/moosefs/"):
        suffix = "paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/"
        idx = path.find(suffix)
        if idx >= 0:
            local = os.path.join(repo_root, path[idx:])
            if os.path.exists(local):
                return local
    return path


def tag_value(fields, tag):
    prefix = f"{tag}:"
    for field in fields[12:]:
        if field.startswith(prefix):
            return field.split(":", 2)[-1]
    return ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-dir", required=True)
    parser.add_argument("--source-dir", default="paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16")
    parser.add_argument("--lengths", default="10000,5000,2000")
    args = parser.parse_args()

    repo_root = os.path.abspath(os.path.join(args.panel_dir, "..", "..", ".."))
    lengths = tuple(int(x) for x in args.lengths.replace(",", " ").split())
    source_dir = args.source_dir
    windows = read_tsv(os.path.join(args.panel_dir, "config", "panel_windows.tsv"))
    manifest = read_tsv(os.path.join(source_dir, "summaries", "query_grid_chop_filter_manifest.tsv"))
    audit = read_tsv(os.path.join(source_dir, "summaries", "query_grid_overlap_audit.tsv"))

    manifest_by_key = {}
    for row in manifest:
        if row["filter_id"] != FILTER_ID or row["chunk_mode"] != "query-grid":
            continue
        key = (row["comparison_id"], int(row["chop_length_bp"]))
        row = dict(row)
        row["filtered_paf_resolved"] = first_existing_path(row["filtered_paf"], repo_root)
        manifest_by_key[key] = row

    audit_by_key = {}
    for row in audit:
        key = (row["event_id"], row["comparison_id"], int(row["chop_length_bp"]))
        audit_by_key[key] = row

    by_comparison = defaultdict(list)
    for window in windows:
        by_comparison[window["comparison_id"]].append(window)

    segment_rows = []
    manifest_rows = []
    for comparison_id, comparison_windows in sorted(by_comparison.items()):
        for length in lengths:
            key = (comparison_id, length)
            if key not in manifest_by_key:
                raise SystemExit(f"missing manifest row for {comparison_id} length {length}")
            row = manifest_by_key[key]
            filtered_paf = row["filtered_paf_resolved"]
            if not os.path.exists(filtered_paf):
                raise SystemExit(f"missing filtered PAF: {filtered_paf}")
            expected_sha = row.get("filtered_sha256", "")
            actual_sha = sha256(filtered_paf)
            if expected_sha and actual_sha != expected_sha:
                raise SystemExit(f"sha256 mismatch for {filtered_paf}: {actual_sha} != {expected_sha}")
            for window in comparison_windows:
                manifest_rows.append({
                    "event_id": window["event_id"],
                    "panel_label": window["panel_label"],
                    "comparison_id": comparison_id,
                    "query_name": window["query_name"],
                    "window_start": window["query_start"],
                    "window_end": window["query_end"],
                    "chop_length_bp": str(length),
                    "filter_id": FILTER_ID,
                    "chunk_mode": row["chunk_mode"],
                    "output_family": row["output_family"],
                    "num_mappings": row["num_mappings"],
                    "scaffold_jump": row["scaffold_jump"],
                    "scoring": row["scoring"],
                    "filter_overlap": row["filter_overlap"],
                    "status": row["status"],
                    "raw_paf": row["raw_paf"],
                    "chopped_paf": row["chopped_paf"],
                    "filtered_paf": row["filtered_paf"],
                    "filtered_paf_resolved": filtered_paf,
                    "filtered_sha256": actual_sha,
                })
            with gzip.open(filtered_paf, "rt") as handle:
                for line_no, line in enumerate(handle, start=1):
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 12:
                        continue
                    q_name = fields[0]
                    q_start = int(fields[2])
                    q_end = int(fields[3])
                    if not any(q_name == w["query_name"] and q_end > int(w["query_start"]) and q_start < int(w["query_end"]) for w in comparison_windows):
                        continue
                    t_name = fields[5]
                    t_start = int(fields[7])
                    t_end = int(fields[8])
                    matches = int(fields[9])
                    aln_len = int(fields[10])
                    target_chrom = chrom_from_name(t_name)
                    identity = matches / aln_len if aln_len else 0.0
                    for window in comparison_windows:
                        if q_name != window["query_name"]:
                            continue
                        w_start = int(window["query_start"])
                        w_end = int(window["query_end"])
                        ov = overlap_bp(q_start, q_end, w_start, w_end)
                        if ov <= 0:
                            continue
                        expected = set(window["expected_target_chroms"].split(","))
                        clip_start = max(q_start, w_start)
                        clip_end = min(q_end, w_end)
                        segment_rows.append({
                            "event_id": window["event_id"],
                            "comparison_id": comparison_id,
                            "panel_label": window["panel_label"],
                            "query_name": q_name,
                            "query_chrom": chrom_from_name(q_name),
                            "window_start": str(w_start),
                            "window_end": str(w_end),
                            "chop_length_bp": str(length),
                            "filter_id": FILTER_ID,
                            "chunk_mode": "query-grid",
                            "num_mappings": row["num_mappings"],
                            "scaffold_jump": row["scaffold_jump"],
                            "scoring": row["scoring"],
                            "filter_overlap": row["filter_overlap"],
                            "query_start": str(q_start),
                            "query_end": str(q_end),
                            "query_clip_start": str(clip_start),
                            "query_clip_end": str(clip_end),
                            "window_overlap_bp": str(ov),
                            "target_name": t_name,
                            "target_chrom": target_chrom,
                            "target_start": str(t_start),
                            "target_end": str(t_end),
                            "strand": fields[4],
                            "matches": str(matches),
                            "alignment_length": str(aln_len),
                            "identity": f"{identity:.6f}",
                            "query_grid_chunk_index": tag_value(fields, "zc:i"),
                            "query_grid_chunk_length": tag_value(fields, "zl:i"),
                            "query_grid_chunk_overlap": tag_value(fields, "zo:i"),
                            "query_grid_mode": tag_value(fields, "zm:Z"),
                            "expected_target_chroms": window["expected_target_chroms"],
                            "is_expected_target": "yes" if target_chrom in expected else "no",
                            "filtered_paf": row["filtered_paf"],
                            "filtered_paf_resolved": filtered_paf,
                            "filtered_paf_line": str(line_no),
                        })

    segment_groups = defaultdict(list)
    for row in segment_rows:
        key = (row["event_id"], int(row["chop_length_bp"]))
        segment_groups[key].append(row)

    summary_rows = []
    for window in windows:
        for length in lengths:
            rows = segment_groups.get((window["event_id"], length), [])
            expected_rows = [r for r in rows if r["is_expected_target"] == "yes"]
            by_target_sum = defaultdict(int)
            by_target_intervals = defaultdict(list)
            for row in rows:
                by_target_sum[row["target_chrom"]] += int(row["window_overlap_bp"])
                by_target_intervals[row["target_chrom"]].append((int(row["query_clip_start"]), int(row["query_clip_end"])))
            audit_row = audit_by_key.get((window["event_id"], window["comparison_id"], length), {})
            summary_rows.append({
                "event_id": window["event_id"],
                "comparison_id": window["comparison_id"],
                "query_name": window["query_name"],
                "query_chrom": chrom_from_name(window["query_name"]),
                "window_start": window["query_start"],
                "window_end": window["query_end"],
                "chop_length_bp": str(length),
                "filter_id": FILTER_ID,
                "chunk_mode": "query-grid",
                "num_mappings": "1:1",
                "scaffold_jump": "0",
                "scoring": "ani",
                "filter_overlap": "0",
                "expected_target_chroms": window["expected_target_chroms"],
                "segment_rows": str(len(rows)),
                "expected_target_rows": str(len(expected_rows)),
                "sum_expected_overlap_bp": str(sum(int(r["window_overlap_bp"]) for r in expected_rows)),
                "union_expected_overlap_bp": str(union_bp([(int(r["query_clip_start"]), int(r["query_clip_end"])) for r in expected_rows])),
                "target_sum_overlap_bp": ";".join(f"{k}:{v}" for k, v in sorted(by_target_sum.items())),
                "target_union_overlap_bp": ";".join(f"{k}:{union_bp(v)}" for k, v in sorted(by_target_intervals.items())),
                "audit_boundary_violation_rows": audit_row.get("boundary_violation_rows", ""),
                "audit_query_redundant_bp": audit_row.get("query_redundant_bp", ""),
                "audit_chr3_query_redundant_bp": audit_row.get("chr3_query_redundant_bp", ""),
                "audit_conclusion": audit_row.get("audit_conclusion", ""),
                "status": "OK" if expected_rows else "NO_EXPECTED_TARGET_ROWS",
            })

    segment_fields = [
        "event_id", "comparison_id", "panel_label", "query_name", "query_chrom", "window_start", "window_end",
        "chop_length_bp", "filter_id", "chunk_mode", "num_mappings", "scaffold_jump", "scoring", "filter_overlap",
        "query_start", "query_end", "query_clip_start", "query_clip_end", "window_overlap_bp",
        "target_name", "target_chrom", "target_start", "target_end", "strand", "matches", "alignment_length",
        "identity", "query_grid_chunk_index", "query_grid_chunk_length", "query_grid_chunk_overlap",
        "query_grid_mode", "expected_target_chroms", "is_expected_target", "filtered_paf",
        "filtered_paf_resolved", "filtered_paf_line",
    ]
    summary_fields = [
        "event_id", "comparison_id", "query_name", "query_chrom", "window_start", "window_end",
        "chop_length_bp", "filter_id", "chunk_mode", "num_mappings", "scaffold_jump", "scoring",
        "filter_overlap", "expected_target_chroms", "segment_rows", "expected_target_rows",
        "sum_expected_overlap_bp", "union_expected_overlap_bp", "target_sum_overlap_bp",
        "target_union_overlap_bp", "audit_boundary_violation_rows", "audit_query_redundant_bp",
        "audit_chr3_query_redundant_bp", "audit_conclusion", "status",
    ]
    manifest_fields = [
        "event_id", "panel_label", "comparison_id", "query_name", "window_start", "window_end",
        "chop_length_bp", "filter_id", "chunk_mode", "output_family", "num_mappings",
        "scaffold_jump", "scoring", "filter_overlap", "status", "raw_paf",
        "chopped_paf", "filtered_paf", "filtered_paf_resolved", "filtered_sha256",
    ]
    write_tsv(os.path.join(args.panel_dir, "query_grid_panel_segments.tsv"), segment_rows, segment_fields)
    write_tsv(os.path.join(args.panel_dir, "query_grid_panel_summary.tsv"), summary_rows, summary_fields)
    write_tsv(os.path.join(args.panel_dir, "query_grid_panel_manifest.tsv"), manifest_rows, manifest_fields)


if __name__ == "__main__":
    main()
