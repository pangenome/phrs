#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import os
import re
from collections import defaultdict


def read_tsv(path):
    with open(path, newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


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


def load_raw_jobs(source_package):
    jobs = {}
    path = os.path.join(source_package, "summaries", "slurm_jobs.tsv")
    for row in read_tsv(path):
        if row.get("stage") == "raw_frequency16_primary":
            jobs[row["comparison_id"]] = row
    return jobs


def overlap_bp(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-dir", required=True)
    parser.add_argument("--source-package", required=True)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--chop-length", default="2000")
    parser.add_argument("--overlap", default="0")
    args = parser.parse_args()

    windows_path = os.path.join(args.panel_dir, "config", "panel_windows.tsv")
    windows = list(read_tsv(windows_path))
    raw_jobs = load_raw_jobs(args.source_package)
    by_comparison = defaultdict(list)
    for window in windows:
        by_comparison[window["comparison_id"]].append(window)

    segment_rows = []
    summary_rows = []
    provenance_rows = []

    for comparison_id, comparison_windows in sorted(by_comparison.items()):
        raw_paf = os.path.join(
            args.source_package,
            "raw_paf",
            f"{comparison_id}.sweepga_frequency16_many_many_j0.paf.gz",
        )
        chopped_paf = os.path.join(
            args.panel_dir,
            "work",
            f"chopped_paf_l{args.chop_length}_o{args.overlap}",
            f"{comparison_id}.chopped_l{args.chop_length}_o{args.overlap}.paf.gz",
        )
        raw_window_paf = os.path.join(
            args.panel_dir,
            "work",
            "raw_window_paf",
            f"{comparison_id}.raw_window_extract.paf.gz",
        )
        filtered_paf = os.path.join(
            args.panel_dir,
            "evidence_paf",
            f"{comparison_id}.window_extract.one_one_scoring_ani_chopped_l{args.chop_length}_o{args.overlap}.paf.gz",
        )
        if not os.path.exists(filtered_paf):
            raise SystemExit(f"missing filtered PAF: {filtered_paf}")

        raw_job = raw_jobs.get(comparison_id, {})
        provenance_rows.append(
            {
                "comparison_id": comparison_id,
                "panel_slurm_job_id": args.job_id,
                "raw_frequency16_slurm_job_id": raw_job.get("job_id", ""),
                "raw_frequency16_sacct_state": raw_job.get("sacct_state", ""),
                "raw_frequency16_elapsed": raw_job.get("elapsed", ""),
                "raw_frequency16_node": raw_job.get("node", ""),
                "source_package": args.source_package,
                "raw_paf": raw_paf,
                "raw_paf_sha256": sha256(raw_paf),
                "raw_window_extract_paf": raw_window_paf,
                "raw_window_extract_paf_sha256": sha256(raw_window_paf),
                "chopped_paf_l2000_o0": chopped_paf,
                "chopped_paf_sha256": sha256(chopped_paf),
                "filtered_paf_1to1_scoring_ani": filtered_paf,
                "filtered_paf_sha256": sha256(filtered_paf),
                "chop_length_bp": args.chop_length,
                "overlap_bp": args.overlap,
                "filter_command": "sweepga --num-mappings 1:1 --scaffold-jump 0 --scoring ani",
                "raw_command_log": raw_job.get("command_log", ""),
            }
        )

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
                    ov = overlap_bp(q_start, q_end, int(window["query_start"]), int(window["query_end"]))
                    if ov <= 0:
                        continue
                    expected = set(window["expected_target_chroms"].split(","))
                    segment_rows.append(
                        {
                            "event_id": window["event_id"],
                            "comparison_id": comparison_id,
                            "panel_label": window["panel_label"],
                            "query_name": q_name,
                            "window_start": window["query_start"],
                            "window_end": window["query_end"],
                            "query_start": q_start,
                            "query_end": q_end,
                            "query_rel_start": max(q_start, int(window["query_start"])) - int(window["query_start"]),
                            "query_rel_end": min(q_end, int(window["query_end"])) - int(window["query_start"]),
                            "window_overlap_bp": ov,
                            "target_name": t_name,
                            "target_chrom": target_chrom,
                            "target_start": t_start,
                            "target_end": t_end,
                            "strand": fields[4],
                            "matches": matches,
                            "alignment_length": aln_len,
                            "identity": f"{identity:.6f}",
                            "expected_target_chroms": window["expected_target_chroms"],
                            "is_expected_target": "yes" if target_chrom in expected else "no",
                            "source_layer": f"raw_fasta_sweepga_fastga_f16_pafchop_l{args.chop_length}_sweepga_1to1_scoring_ani",
                            "filtered_paf": filtered_paf,
                            "filtered_paf_line": line_no,
                        }
                    )

    for window in windows:
        rows = [r for r in segment_rows if r["event_id"] == window["event_id"]]
        expected_rows = [r for r in rows if r["is_expected_target"] == "yes"]
        by_target = defaultdict(int)
        for row in rows:
            by_target[row["target_chrom"]] += int(row["window_overlap_bp"])
        summary_rows.append(
            {
                "event_id": window["event_id"],
                "comparison_id": window["comparison_id"],
                "query_name": window["query_name"],
                "window_start": window["query_start"],
                "window_end": window["query_end"],
                "expected_target_chroms": window["expected_target_chroms"],
                "segment_rows": len(rows),
                "expected_target_rows": len(expected_rows),
                "sum_expected_overlap_bp": sum(int(r["window_overlap_bp"]) for r in expected_rows),
                "target_overlap_bp": ";".join(f"{k}:{v}" for k, v in sorted(by_target.items())),
                "status": "OK" if expected_rows else "NO_EXPECTED_TARGET_ROWS",
            }
        )

    segment_fields = [
        "event_id",
        "comparison_id",
        "panel_label",
        "query_name",
        "window_start",
        "window_end",
        "query_start",
        "query_end",
        "query_rel_start",
        "query_rel_end",
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
        "source_layer",
        "filtered_paf",
        "filtered_paf_line",
    ]
    summary_fields = [
        "event_id",
        "comparison_id",
        "query_name",
        "window_start",
        "window_end",
        "expected_target_chroms",
        "segment_rows",
        "expected_target_rows",
        "sum_expected_overlap_bp",
        "target_overlap_bp",
        "status",
    ]
    provenance_fields = [
        "comparison_id",
        "panel_slurm_job_id",
        "raw_frequency16_slurm_job_id",
        "raw_frequency16_sacct_state",
        "raw_frequency16_elapsed",
        "raw_frequency16_node",
        "source_package",
        "raw_paf",
        "raw_paf_sha256",
        "raw_window_extract_paf",
        "raw_window_extract_paf_sha256",
        "chopped_paf_l2000_o0",
        "chopped_paf_sha256",
        "filtered_paf_1to1_scoring_ani",
        "filtered_paf_sha256",
        "chop_length_bp",
        "overlap_bp",
        "filter_command",
        "raw_command_log",
    ]

    write_tsv(os.path.join(args.panel_dir, "raw_fasta_chopped_panel_segments.tsv"), segment_rows, segment_fields)
    write_tsv(os.path.join(args.panel_dir, "raw_fasta_chopped_panel_summary.tsv"), summary_rows, summary_fields)
    write_tsv(os.path.join(args.panel_dir, "slurm_jobs.tsv"), provenance_rows, provenance_fields)


if __name__ == "__main__":
    main()
