#!/usr/bin/env python3
import csv
import gzip
import hashlib
import os
import re
from collections import defaultdict
from pathlib import Path


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
    merged = []
    for start, end in sorted(intervals):
        if start >= end:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return sum(end - start for start, end in merged)


def open_paf(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def scan_paf(path, windows, cell):
    segment_rows = []
    by_query = defaultdict(list)
    for window in windows:
        by_query[window["query_name"]].append(window)
    with open_paf(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            q_name = fields[0]
            if q_name not in by_query:
                continue
            q_start = int(fields[2])
            q_end = int(fields[3])
            candidate_windows = [
                window for window in by_query[q_name]
                if q_end > int(window["query_start"]) and q_start < int(window["query_end"])
            ]
            if not candidate_windows:
                continue
            t_name = fields[5]
            target_chrom = chrom_from_name(t_name)
            matches = int(fields[9])
            aln_len = int(fields[10])
            identity = matches / aln_len if aln_len else 0.0
            for window in candidate_windows:
                w_start = int(window["query_start"])
                w_end = int(window["query_end"])
                clip_start = max(q_start, w_start)
                clip_end = min(q_end, w_end)
                ov = overlap_bp(q_start, q_end, w_start, w_end)
                if ov <= 0:
                    continue
                expected = set(window["expected_target_chroms"].split(","))
                segment_rows.append(
                    {
                        **cell,
                        "event_id": window["event_id"],
                        "event_label": window["event_label"],
                        "comparison_id": window["comparison_id"],
                        "query_name": q_name,
                        "window_start": str(w_start),
                        "window_end": str(w_end),
                        "query_start": str(q_start),
                        "query_end": str(q_end),
                        "query_clip_start": str(clip_start),
                        "query_clip_end": str(clip_end),
                        "window_overlap_bp": str(ov),
                        "target_name": t_name,
                        "target_chrom": target_chrom,
                        "target_start": fields[7],
                        "target_end": fields[8],
                        "strand": fields[4],
                        "matches": str(matches),
                        "alignment_length": str(aln_len),
                        "identity": f"{identity:.6f}",
                        "expected_target_chroms": window["expected_target_chroms"],
                        "is_expected_target": "yes" if target_chrom in expected else "no",
                        "paf_path": path,
                        "paf_line": str(line_no),
                    }
                )
    return segment_rows


def summarize(windows, cells, segment_rows):
    groups = defaultdict(list)
    for row in segment_rows:
        groups[(row["event_id"], row["cell_id"])].append(row)
    summary_rows = []
    target_rows = []
    for window in windows:
        expected_chroms = set(window["expected_target_chroms"].split(","))
        for cell in [cell for cell in cells if cell["comparison_id"] == window["comparison_id"]]:
            rows = groups.get((window["event_id"], cell["cell_id"]), [])
            expected_rows = [row for row in rows if row["target_chrom"] in expected_chroms]
            target_sum = defaultdict(int)
            target_intervals = defaultdict(list)
            for row in rows:
                target_sum[row["target_chrom"]] += int(row["window_overlap_bp"])
                target_intervals[row["target_chrom"]].append((int(row["query_clip_start"]), int(row["query_clip_end"])))
            all_union = union_bp([iv for intervals in target_intervals.values() for iv in intervals])
            expected_union = union_bp([(int(row["query_clip_start"]), int(row["query_clip_end"])) for row in expected_rows])
            chr3_union = union_bp(target_intervals.get("chr3", []))
            if chr3_union > 0 and "chr3" in expected_chroms:
                status = "chr3_expected"
            elif chr3_union > 0:
                status = "chr3_off_target"
            elif expected_rows:
                status = "expected_non_chr3"
            elif rows:
                status = "other_target_only"
            else:
                status = "no_rows"
            summary = {
                **cell,
                "event_id": window["event_id"],
                "event_label": window["event_label"],
                "comparison_id": window["comparison_id"],
                "query_name": window["query_name"],
                "window_start": window["query_start"],
                "window_end": window["query_end"],
                "expected_target_chroms": window["expected_target_chroms"],
                "row_count": str(len(rows)),
                "expected_target_rows": str(len(expected_rows)),
                "expected_target_sum_bp": str(sum(int(row["window_overlap_bp"]) for row in expected_rows)),
                "expected_target_union_bp": str(expected_union),
                "all_target_union_bp": str(all_union),
                "chr3_union_bp": str(chr3_union),
                "chr3_survives": "yes" if chr3_union > 0 else "no",
                "target_sum_overlap_bp": ";".join(f"{chrom}:{target_sum[chrom]}" for chrom in sorted(target_sum)),
                "target_union_overlap_bp": ";".join(f"{chrom}:{union_bp(target_intervals[chrom])}" for chrom in sorted(target_intervals)),
                "status": status,
            }
            summary_rows.append(summary)
            for chrom in sorted(target_intervals):
                target_rows.append(
                    {
                        **cell,
                        "event_id": window["event_id"],
                        "comparison_id": window["comparison_id"],
                        "target_chrom": chrom,
                        "target_sum_bp": str(target_sum[chrom]),
                        "target_union_bp": str(union_bp(target_intervals[chrom])),
                        "target_rows": str(sum(1 for row in rows if row["target_chrom"] == chrom)),
                        "is_expected_target": "yes" if chrom in expected_chroms else "no",
                    }
                )
    return summary_rows, target_rows


def main():
    repo = Path(os.environ.get("REPO_ROOT", "/moosefs/erikg/phrs"))
    package = repo / "paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_scaffold_jump_sensitivity"
    windows = read_tsv(package / "config/candidate_windows.tsv")
    tasks = read_tsv(package / "filter_tasks.tsv")
    source_by_comparison = {row["comparison_id"]: row["source_paf"] for row in read_tsv(package / "config/source_pafs.tsv")}

    cells = []
    segment_rows = []
    manifest_rows = []
    source_hash_cache = {}

    for comparison_id, source in sorted(source_by_comparison.items()):
        source_hash_cache[source] = sha256(source)
        cell = {
            "cell_id": "raw_many_many_unfiltered",
            "comparison_id": comparison_id,
            "source_kind": "raw_baseline",
            "scaffold_jump": "0",
            "min_aln_length": "raw",
            "scoring": "raw",
            "num_mappings": "many:many",
            "scaffold_mass": "10k",
            "overlap": "default",
            "paf_path": source,
        }
        cells.append(cell)
        manifest_rows.append(
            {
                **cell,
                "source_paf": source,
                "source_paf_sha256": source_hash_cache[source],
                "filtered_paf": source,
                "filtered_paf_sha256": source_hash_cache[source],
                "command": "raw many:many f16 source PAF; no final paf-filter applied",
            }
        )
        comparison_windows = [window for window in windows if window["comparison_id"] == comparison_id]
        segment_rows.extend(scan_paf(source, comparison_windows, cell))

    for task in tasks:
        out = task["output_paf"]
        if not os.path.exists(out):
            raise SystemExit(f"missing filtered PAF: {out}")
        source = task["source_paf"]
        if source not in source_hash_cache:
            source_hash_cache[source] = sha256(source)
        filtered_hash = sha256(out)
        cell = {
            "cell_id": task["cell_id"],
            "comparison_id": task["comparison_id"],
            "source_kind": "filtered",
            "scaffold_jump": task["scaffold_jump"],
            "min_aln_length": task["min_aln_length"],
            "scoring": task["scoring"],
            "num_mappings": task["num_mappings"],
            "scaffold_mass": task["scaffold_mass"],
            "overlap": task["overlap"],
            "paf_path": out,
        }
        cells.append(cell)
        manifest_rows.append(
            {
                **cell,
                "source_paf": source,
                "source_paf_sha256": source_hash_cache[source],
                "filtered_paf": out,
                "filtered_paf_sha256": filtered_hash,
                "command": task["command"],
            }
        )
        comparison_windows = [window for window in windows if window["comparison_id"] == task["comparison_id"]]
        segment_rows.extend(scan_paf(out, comparison_windows, cell))

    summary_rows, target_rows = summarize(windows, cells, segment_rows)

    write_tsv(package / "candidate_window_segments.tsv", segment_rows, [
        "event_id", "event_label", "comparison_id", "query_name", "window_start", "window_end",
        "cell_id", "source_kind", "scaffold_jump", "min_aln_length", "scoring", "num_mappings",
        "scaffold_mass", "overlap", "query_start", "query_end", "query_clip_start", "query_clip_end",
        "window_overlap_bp", "target_name", "target_chrom", "target_start", "target_end", "strand",
        "matches", "alignment_length", "identity", "expected_target_chroms", "is_expected_target",
        "paf_path", "paf_line",
    ])
    write_tsv(package / "candidate_window_summary.tsv", summary_rows, [
        "event_id", "event_label", "comparison_id", "query_name", "window_start", "window_end",
        "cell_id", "source_kind", "scaffold_jump", "min_aln_length", "scoring", "num_mappings",
        "scaffold_mass", "overlap", "expected_target_chroms", "row_count", "expected_target_rows",
        "expected_target_sum_bp", "expected_target_union_bp", "all_target_union_bp", "chr3_union_bp",
        "chr3_survives", "target_sum_overlap_bp", "target_union_overlap_bp", "status",
    ])
    write_tsv(package / "target_chrom_breakdown.tsv", target_rows, [
        "event_id", "comparison_id", "cell_id", "source_kind", "scaffold_jump", "min_aln_length",
        "scoring", "num_mappings", "scaffold_mass", "overlap", "target_chrom", "target_sum_bp",
        "target_union_bp", "target_rows", "is_expected_target",
    ])
    write_tsv(package / "filtered_paf_manifest.tsv", manifest_rows, [
        "comparison_id", "cell_id", "source_kind", "scaffold_jump", "min_aln_length", "scoring",
        "num_mappings", "scaffold_mass", "overlap", "source_paf", "source_paf_sha256",
        "filtered_paf", "filtered_paf_sha256", "paf_path", "command",
    ])


if __name__ == "__main__":
    main()
