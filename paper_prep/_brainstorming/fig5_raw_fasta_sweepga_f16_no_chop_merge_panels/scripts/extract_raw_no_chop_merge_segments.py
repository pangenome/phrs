#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import re
from collections import defaultdict


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


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


def parse_target_union(overlaps):
    return ";".join("%s:%d" % (chrom, union_bp(vals)) for chrom, vals in sorted(overlaps.items()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-dir", required=True)
    parser.add_argument("--source-package", required=True)
    args = parser.parse_args()

    windows = read_tsv(os.path.join(args.panel_dir, "config", "panel_windows.tsv"))
    modes = sorted(read_tsv(os.path.join(args.panel_dir, "config", "filter_modes.tsv")), key=lambda row: int(row["plot_order"]))
    by_comparison = defaultdict(list)
    for window in windows:
        by_comparison[window["comparison_id"]].append(window)

    segment_rows = []
    for mode in modes:
        for comparison_id, comparison_windows in sorted(by_comparison.items()):
            filtered_paf = os.path.join(args.panel_dir, "filtered_paf", "%s.%s.raw_no_chop.paf.gz" % (comparison_id, mode["filter_mode"]))
            if not os.path.exists(filtered_paf):
                raise SystemExit("missing filtered PAF: %s" % filtered_paf)
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
                    candidate_windows = [
                        w for w in comparison_windows
                        if q_name == w["query_name"] and q_end > int(w["query_start"]) and q_start < int(w["query_end"])
                    ]
                    if not candidate_windows:
                        continue
                    t_name = fields[5]
                    t_start = int(fields[7])
                    t_end = int(fields[8])
                    matches = int(fields[9])
                    aln_len = int(fields[10])
                    target_chrom = chrom_from_name(t_name)
                    identity = matches / aln_len if aln_len else 0.0
                    for window in candidate_windows:
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
                            "window_start": str(w_start),
                            "window_end": str(w_end),
                            "filter_mode": mode["filter_mode"],
                            "filter_label": mode["filter_label"],
                            "num_mappings": mode["num_mappings"],
                            "scaffold_jump": mode["scaffold_jump"],
                            "scoring": mode["scoring"],
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
                            "identity": "%.6f" % identity,
                            "expected_target_chroms": window["expected_target_chroms"],
                            "is_expected_target": "yes" if target_chrom in expected else "no",
                            "source_layer": "raw_fasta_sweepga_fastga_f16_no_chop_%s" % mode["filter_mode"],
                            "filtered_paf": filtered_paf,
                            "filtered_paf_line": str(line_no),
                        })

    grouped = defaultdict(list)
    for row in segment_rows:
        grouped[(row["event_id"], row["filter_mode"])].append(row)

    summary_rows = []
    for window in windows:
        for mode in modes:
            rows = grouped.get((window["event_id"], mode["filter_mode"]), [])
            expected_rows = [r for r in rows if r["is_expected_target"] == "yes"]
            by_target_sum = defaultdict(int)
            by_target_intervals = defaultdict(list)
            expected_intervals = []
            for row in rows:
                by_target_sum[row["target_chrom"]] += int(row["window_overlap_bp"])
                by_target_intervals[row["target_chrom"]].append((int(row["query_clip_start"]), int(row["query_clip_end"])))
            for row in expected_rows:
                expected_intervals.append((int(row["query_clip_start"]), int(row["query_clip_end"])))
            union_expected = union_bp(expected_intervals)
            chr3_survival = "NA"
            if window["expected_target_chroms"] == "chr3":
                chr3_survival = "CHR3_SURVIVES" if union_expected > 0 else "CHR3_DROPPED"
            summary_rows.append({
                "event_id": window["event_id"],
                "comparison_id": window["comparison_id"],
                "query_name": window["query_name"],
                "window_start": window["query_start"],
                "window_end": window["query_end"],
                "filter_mode": mode["filter_mode"],
                "filter_label": mode["filter_label"],
                "num_mappings": mode["num_mappings"],
                "scaffold_jump": mode["scaffold_jump"],
                "scoring": mode["scoring"],
                "expected_target_chroms": window["expected_target_chroms"],
                "segment_rows": str(len(rows)),
                "expected_target_rows": str(len(expected_rows)),
                "sum_expected_overlap_bp": str(sum(int(r["window_overlap_bp"]) for r in expected_rows)),
                "union_expected_overlap_bp": str(union_expected),
                "target_sum_overlap_bp": ";".join("%s:%d" % (k, v) for k, v in sorted(by_target_sum.items())),
                "target_union_overlap_bp": parse_target_union(by_target_intervals),
                "chr3_survival_status": chr3_survival,
                "status": "OK" if expected_rows else "NO_EXPECTED_TARGET_ROWS",
            })

    segment_fields = [
        "event_id", "comparison_id", "panel_label", "query_name", "window_start", "window_end",
        "filter_mode", "filter_label", "num_mappings", "scaffold_jump", "scoring",
        "query_start", "query_end", "query_clip_start", "query_clip_end", "window_overlap_bp",
        "target_name", "target_chrom", "target_start", "target_end", "strand", "matches",
        "alignment_length", "identity", "expected_target_chroms", "is_expected_target",
        "source_layer", "filtered_paf", "filtered_paf_line",
    ]
    summary_fields = [
        "event_id", "comparison_id", "query_name", "window_start", "window_end",
        "filter_mode", "filter_label", "num_mappings", "scaffold_jump", "scoring",
        "expected_target_chroms", "segment_rows", "expected_target_rows",
        "sum_expected_overlap_bp", "union_expected_overlap_bp",
        "target_sum_overlap_bp", "target_union_overlap_bp", "chr3_survival_status", "status",
    ]
    write_tsv(os.path.join(args.panel_dir, "raw_merge_panel_segments.tsv"), segment_rows, segment_fields)
    write_tsv(os.path.join(args.panel_dir, "raw_merge_panel_summary.tsv"), summary_rows, summary_fields)


if __name__ == "__main__":
    main()
