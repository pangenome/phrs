#!/usr/bin/env python3
"""Extract Fig5 window-level support from whole-genome raw f16 many:many PAFs.

The primary input is the raw SweepGA/FastGA f16 whole-genome many:many PAF.  The
filtered/chopped layers are read only as comparison layers so the resulting
summary explicitly separates pre-filter multiway support from final 1:1 support.
"""

import argparse
import csv
import gzip
import os
import re
from collections import defaultdict
from pathlib import Path


RAW_SUFFIX = ".sweepga_frequency16_many_many_j0.paf.gz"
FILTER_LAYERS = (
    ("chopped_many_many", ".many_many_chopped.paf.gz"),
    ("chopped_four_many", ".four_many_chopped.paf.gz"),
    ("chopped_one_one", ".one_one_chopped.paf.gz"),
)


def read_tsv(path):
    with open(path, newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_checksum(path):
    checksum_path = Path(str(path) + ".sha256")
    if not checksum_path.exists():
        return ""
    text = checksum_path.read_text().strip()
    return text.split()[0] if text else ""


def chrom_from_name(name):
    matches = re.findall(r"chr(?:[0-9]+|X|Y|M)", name)
    return matches[-1] if matches else name


def parse_query(name):
    parts = name.split("#")
    return {
        "query_sample": parts[0] if len(parts) > 0 else "",
        "query_hap": parts[1] if len(parts) > 1 else "",
        "query_chrom": parts[2] if len(parts) > 2 else chrom_from_name(name),
    }


def parse_target(name):
    # Examples: PAN011#joint#h2_chr3, PAN027#joint#h1_chr9.
    parts = name.split("#")
    sample = parts[0] if parts else ""
    hap = ""
    if len(parts) >= 3:
        m = re.match(r"(h[12])(?:_|$)", parts[2])
        if m:
            hap = m.group(1)
    return {
        "target_sample": sample,
        "target_hap": hap,
        "target_chrom": chrom_from_name(name),
    }


def overlap_bp(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def union_bp(intervals):
    if not intervals:
        return 0
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def target_bucket(target_chrom):
    if target_chrom == "chr3":
        return "chr3"
    if target_chrom == "chr9":
        return "chr9"
    return "other"


def load_raw_jobs(source_package):
    jobs = {}
    path = source_package / "summaries" / "slurm_jobs.tsv"
    if not path.exists():
        return jobs
    for row in read_tsv(path):
        if row.get("stage") == "raw_frequency16_primary":
            jobs[row["comparison_id"]] = row
    return jobs


def open_paf(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def row_from_paf(fields, line_no, layer, paf_path, window):
    q_name = fields[0]
    q_len = int(fields[1])
    q_start = int(fields[2])
    q_end = int(fields[3])
    strand = fields[4]
    t_name = fields[5]
    t_len = int(fields[6])
    t_start = int(fields[7])
    t_end = int(fields[8])
    matches = int(fields[9])
    aln_len = int(fields[10])
    mapq = fields[11]
    q_clip_start = max(q_start, window["query_start_i"])
    q_clip_end = min(q_end, window["query_end_i"])
    identity = matches / aln_len if aln_len else 0.0
    parsed_q = parse_query(q_name)
    parsed_t = parse_target(t_name)
    expected = set(window["expected_target_chroms"].split(","))
    return {
        "event_id": window["event_id"],
        "comparison_id": window["comparison_id"],
        "panel_label": window["panel_label"],
        "source_layer": layer,
        "query_name": q_name,
        "query_sample": parsed_q["query_sample"],
        "query_hap": parsed_q["query_hap"],
        "query_chrom": parsed_q["query_chrom"],
        "window_start": str(window["query_start_i"]),
        "window_end": str(window["query_end_i"]),
        "query_start": str(q_start),
        "query_end": str(q_end),
        "query_clip_start": str(q_clip_start),
        "query_clip_end": str(q_clip_end),
        "query_rel_start": str(q_clip_start - window["query_start_i"]),
        "query_rel_end": str(q_clip_end - window["query_start_i"]),
        "window_overlap_bp": str(q_clip_end - q_clip_start),
        "query_length": str(q_len),
        "target_name": t_name,
        "target_sample": parsed_t["target_sample"],
        "target_hap": parsed_t["target_hap"],
        "target_chrom": parsed_t["target_chrom"],
        "target_bucket": target_bucket(parsed_t["target_chrom"]),
        "target_start": str(t_start),
        "target_end": str(t_end),
        "target_length": str(t_len),
        "strand": strand,
        "matches": str(matches),
        "alignment_length": str(aln_len),
        "identity": f"{identity:.6f}",
        "block_length": str(aln_len),
        "mapq": mapq,
        "expected_target_chroms": window["expected_target_chroms"],
        "is_expected_target": "yes" if parsed_t["target_chrom"] in expected else "no",
        "source_paf": str(paf_path),
        "source_paf_line": str(line_no),
    }


def extract_layer(paf_path, layer, windows):
    rows = []
    if not paf_path.exists():
        return rows
    with open_paf(paf_path) as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            q_name = fields[0]
            q_start = int(fields[2])
            q_end = int(fields[3])
            for window in windows:
                if q_name != window["query_name"]:
                    continue
                if overlap_bp(q_start, q_end, window["query_start_i"], window["query_end_i"]) <= 0:
                    continue
                rows.append(row_from_paf(fields, line_no, layer, paf_path, window))
    return rows


def summarize(rows, windows, layer):
    summaries = []
    rows_by_event = defaultdict(list)
    for row in rows:
        if row["source_layer"] == layer:
            rows_by_event[row["event_id"]].append(row)
    for window in windows:
        event_rows = rows_by_event.get(window["event_id"], [])
        expected_rows = [row for row in event_rows if row["is_expected_target"] == "yes"]
        intervals_all = [(int(row["query_clip_start"]), int(row["query_clip_end"])) for row in event_rows]
        intervals_expected = [(int(row["query_clip_start"]), int(row["query_clip_end"])) for row in expected_rows]
        by_bucket = defaultdict(list)
        by_chrom = defaultdict(list)
        by_sample_hap = defaultdict(int)
        for row in event_rows:
            interval = (int(row["query_clip_start"]), int(row["query_clip_end"]))
            by_bucket[row["target_bucket"]].append(interval)
            by_chrom[row["target_chrom"]].append(interval)
            target_key = f"{row['target_sample']}:{row['target_hap']}:{row['target_chrom']}"
            by_sample_hap[target_key] += 1
        summaries.append(
            {
                "event_id": window["event_id"],
                "comparison_id": window["comparison_id"],
                "panel_label": window["panel_label"],
                "source_layer": layer,
                "query_name": window["query_name"],
                "window_start": str(window["query_start_i"]),
                "window_end": str(window["query_end_i"]),
                "window_bp": str(window["query_end_i"] - window["query_start_i"]),
                "expected_target_chroms": window["expected_target_chroms"],
                "row_count": str(len(event_rows)),
                "expected_row_count": str(len(expected_rows)),
                "row_multiplicity_per_query_bp": f"{(sum(int(r['window_overlap_bp']) for r in event_rows) / max(1, window['query_end_i'] - window['query_start_i'])):.4f}",
                "query_union_coverage_bp": str(union_bp(intervals_all)),
                "expected_query_union_coverage_bp": str(union_bp(intervals_expected)),
                "chr3_query_union_coverage_bp": str(union_bp(by_bucket.get("chr3", []))),
                "chr9_query_union_coverage_bp": str(union_bp(by_bucket.get("chr9", []))),
                "other_query_union_coverage_bp": str(union_bp(by_bucket.get("other", []))),
                "target_chrom_query_union_bp": ";".join(f"{k}:{union_bp(v)}" for k, v in sorted(by_chrom.items())),
                "target_sample_hap_row_count": ";".join(f"{k}:{v}" for k, v in sorted(by_sample_hap.items())),
                "status": "OK" if event_rows else "NO_ROWS",
            }
        )
    return summaries


def copy_raw_support_files(panel_dir, rows, windows):
    for window in windows:
        event_rows = [
            row for row in rows
            if row["source_layer"] == "raw_many_many_whole_genome" and row["event_id"] == window["event_id"]
        ]
        safe = window["event_id"].replace("/", "_")
        write_tsv(panel_dir / f"{safe}.raw_support.tsv", event_rows, SUPPORT_FIELDS)


def build_manifest(source_package, comparisons, raw_jobs):
    rows = []
    for comparison_id in sorted(comparisons):
        raw_paf = source_package / "raw_paf" / f"{comparison_id}{RAW_SUFFIX}"
        raw_job = raw_jobs.get(comparison_id, {})
        rows.append(
            {
                "comparison_id": comparison_id,
                "source_layer": "raw_many_many_whole_genome",
                "role": "source_of_truth_multiway_input",
                "path": str(raw_paf),
                "sha256": read_checksum(raw_paf),
                "slurm_job_id": raw_job.get("job_id", ""),
                "sacct_state": raw_job.get("sacct_state", ""),
                "elapsed": raw_job.get("elapsed", ""),
                "node": raw_job.get("node", ""),
                "command_log": raw_job.get("command_log", ""),
                "note": "Whole-genome raw f16 SweepGA/FastGA many:many PAF; streamed directly for this package.",
            }
        )
        for layer, suffix in FILTER_LAYERS:
            paf = source_package / "filtered_paf" / f"{comparison_id}{suffix}"
            rows.append(
                {
                    "comparison_id": comparison_id,
                    "source_layer": layer,
                    "role": "comparison_after_chop_filter",
                    "path": str(paf),
                    "sha256": read_checksum(paf),
                    "slurm_job_id": "",
                    "sacct_state": "",
                    "elapsed": "",
                    "node": "",
                    "command_log": "",
                    "note": "Existing chopped/filter layer from the same f16 package; not the source-of-truth multiway input.",
                }
            )
    return rows


SUPPORT_FIELDS = [
    "event_id", "comparison_id", "panel_label", "source_layer",
    "query_name", "query_sample", "query_hap", "query_chrom",
    "window_start", "window_end", "query_start", "query_end",
    "query_clip_start", "query_clip_end", "query_rel_start", "query_rel_end",
    "window_overlap_bp", "query_length",
    "target_name", "target_sample", "target_hap", "target_chrom", "target_bucket",
    "target_start", "target_end", "target_length", "strand",
    "matches", "alignment_length", "identity", "block_length", "mapq",
    "expected_target_chroms", "is_expected_target", "source_paf", "source_paf_line",
]

SUMMARY_FIELDS = [
    "event_id", "comparison_id", "panel_label", "source_layer", "query_name",
    "window_start", "window_end", "window_bp", "expected_target_chroms",
    "row_count", "expected_row_count", "row_multiplicity_per_query_bp",
    "query_union_coverage_bp", "expected_query_union_coverage_bp",
    "chr3_query_union_coverage_bp", "chr9_query_union_coverage_bp",
    "other_query_union_coverage_bp", "target_chrom_query_union_bp",
    "target_sample_hap_row_count", "status",
]

MANIFEST_FIELDS = [
    "comparison_id", "source_layer", "role", "path", "sha256", "slurm_job_id",
    "sacct_state", "elapsed", "node", "command_log", "note",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-dir", default=Path(__file__).resolve().parents[1], type=Path)
    parser.add_argument(
        "--source-package",
        default=Path("/moosefs/erikg/phrs/.wg-worktrees/agent-2649/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16"),
        type=Path,
    )
    parser.add_argument("--raw-only", action="store_true", help="Skip comparison filtered/chopped layers.")
    args = parser.parse_args()

    panel_dir = args.panel_dir.resolve()
    source_package = args.source_package.resolve()
    windows = list(read_tsv(panel_dir / "config" / "panel_windows.tsv"))
    for window in windows:
        window["query_start_i"] = int(window["query_start"])
        window["query_end_i"] = int(window["query_end"])
    by_comparison = defaultdict(list)
    for window in windows:
        by_comparison[window["comparison_id"]].append(window)

    raw_jobs = load_raw_jobs(source_package)
    support_rows = []
    for comparison_id, comparison_windows in sorted(by_comparison.items()):
        raw_paf = source_package / "raw_paf" / f"{comparison_id}{RAW_SUFFIX}"
        support_rows.extend(extract_layer(raw_paf, "raw_many_many_whole_genome", comparison_windows))
        if not args.raw_only:
            for layer, suffix in FILTER_LAYERS:
                paf = source_package / "filtered_paf" / f"{comparison_id}{suffix}"
                support_rows.extend(extract_layer(paf, layer, comparison_windows))

    summary_rows = []
    layers = ["raw_many_many_whole_genome"]
    if not args.raw_only:
        layers.extend(layer for layer, _ in FILTER_LAYERS)
    for layer in layers:
        summary_rows.extend(summarize(support_rows, windows, layer))

    raw_summary = [row for row in summary_rows if row["source_layer"] == "raw_many_many_whole_genome"]
    compact_rows = [
        {
            "event_id": row["event_id"],
            "comparison_id": row["comparison_id"],
            "panel_label": row["panel_label"],
            "query_name": row["query_name"],
            "window_start": row["window_start"],
            "window_end": row["window_end"],
            "raw_many_many_rows": row["row_count"],
            "row_multiplicity_per_query_bp": row["row_multiplicity_per_query_bp"],
            "raw_query_union_coverage_bp": row["query_union_coverage_bp"],
            "raw_chr3_union_bp": row["chr3_query_union_coverage_bp"],
            "raw_chr9_union_bp": row["chr9_query_union_coverage_bp"],
            "raw_other_union_bp": row["other_query_union_coverage_bp"],
            "raw_target_chrom_union_bp": row["target_chrom_query_union_bp"],
            "interpretation": "pre_1to1_filter_multiway_support",
        }
        for row in raw_summary
    ]

    write_tsv(panel_dir / "multiway_candidate_support.tsv", support_rows, SUPPORT_FIELDS)
    write_tsv(panel_dir / "multiway_candidate_summary.tsv", summary_rows, SUMMARY_FIELDS)
    write_tsv(
        panel_dir / "raw_chr3_chr9_other_support_summary.tsv",
        compact_rows,
        [
            "event_id", "comparison_id", "panel_label", "query_name", "window_start", "window_end",
            "raw_many_many_rows", "row_multiplicity_per_query_bp", "raw_query_union_coverage_bp",
            "raw_chr3_union_bp", "raw_chr9_union_bp", "raw_other_union_bp",
            "raw_target_chrom_union_bp", "interpretation",
        ],
    )
    write_tsv(
        panel_dir / "input_manifest.tsv",
        build_manifest(source_package, by_comparison.keys(), raw_jobs),
        MANIFEST_FIELDS,
    )
    copy_raw_support_files(panel_dir, support_rows, windows)


if __name__ == "__main__":
    main()
