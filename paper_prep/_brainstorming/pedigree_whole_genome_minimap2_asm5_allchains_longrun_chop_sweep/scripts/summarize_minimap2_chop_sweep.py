#!/usr/bin/env python3
import collections
import glob
import gzip
import hashlib
import os
import re

from common import CANDIDATE_WINDOWS, DEFAULT_RUN_LABEL, PACKAGE_DIR, paf_records, read_tsv, write_tsv


CHOP_LENGTH = os.environ.get("PAF_CHOP_LENGTH", "10000")
OVERLAP = os.environ.get("PAF_CHOP_OVERLAP", "0")
MINIMAP2_SUFFIX = os.environ.get("MINIMAP2_SUFFIX", "minimap2-v2.31-r1302")
PARAMETER_SET = os.environ.get("MINIMAP2_PARAMETER_SET", "asm5_allchains")


def chrom_from_name(name):
    if "#chr" in name:
        return "chr" + name.rsplit("#chr", 1)[1].split("_", 1)[0].split("#", 1)[0]
    if "_chr" in name:
        return "chr" + name.rsplit("_chr", 1)[1].split("_", 1)[0].split("#", 1)[0]
    match = re.search(r"chr([0-9XYM]+)", name)
    return "chr" + match.group(1) if match else name


def overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def union_bp(intervals):
    if not intervals:
        return 0
    total = 0
    cur_start, cur_end = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start > cur_end:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
        else:
            cur_end = max(cur_end, end)
    return total + cur_end - cur_start


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_status(path):
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1024 * 1024):
                pass
        return "yes"
    except Exception as exc:
        return "no:%s" % exc


def raw_path(comparison_id):
    return os.path.join(
        PACKAGE_DIR,
        "raw_paf",
        DEFAULT_RUN_LABEL,
        "%s.%s.%s.paf.gz" % (comparison_id, PARAMETER_SET, MINIMAP2_SUFFIX),
    )


def layer_paths_for_comparison(comparison_id):
    paths = [
        ("raw_minimap2", "raw", raw_path(comparison_id)),
        (
            "chopped_minimap2_l%s_o%s" % (CHOP_LENGTH, OVERLAP),
            "chopped",
            os.path.join(
                PACKAGE_DIR,
                "chopped_paf_l%s_o%s" % (CHOP_LENGTH, OVERLAP),
                "%s.chopped_l%s_o%s.paf.gz" % (comparison_id, CHOP_LENGTH, OVERLAP),
            ),
        ),
    ]
    for path in sorted(glob.glob(os.path.join(PACKAGE_DIR, "filtered_paf", "%s.*.paf.gz" % comparison_id))):
        filter_id = os.path.basename(path).split(".", 2)[1]
        paths.append(("sweepga_filtered_%s" % filter_id, "filtered", path))
    return paths


def summarize_window(path, window, layer, stage):
    if not os.path.exists(path):
        return [{
            "event_id": window["event_id"],
            "comparison_id": window["comparison_id"],
            "layer": layer,
            "stage": stage,
            "query_name": window["query_name"],
            "query_start": window["query_start"],
            "query_end": window["query_end"],
            "target_chrom": "NO_PAF",
            "expected_target_chrom": window["expected_target_chrom"],
            "paf_rows_overlapping_window": "0",
            "query_overlap_bp_sum": "0",
            "query_covered_bp_union": "0",
            "chr3_support": "not_evaluable",
            "paf_path": path,
        }]
    by_chrom = collections.defaultdict(lambda: {"rows": 0, "sum": 0, "intervals": []})
    qstart = int(window["query_start"])
    qend = int(window["query_end"])
    for fields in paf_records(path):
        if fields[0] != window["query_name"]:
            continue
        row_qstart = int(fields[2])
        row_qend = int(fields[3])
        ov = overlap(row_qstart, row_qend, qstart, qend)
        if ov <= 0:
            continue
        chrom = chrom_from_name(fields[5])
        by_chrom[chrom]["rows"] += 1
        by_chrom[chrom]["sum"] += ov
        by_chrom[chrom]["intervals"].append((max(row_qstart, qstart), min(row_qend, qend)))
    rows = []
    for chrom, vals in sorted(by_chrom.items()):
        rows.append({
            "event_id": window["event_id"],
            "comparison_id": window["comparison_id"],
            "layer": layer,
            "stage": stage,
            "query_name": window["query_name"],
            "query_start": window["query_start"],
            "query_end": window["query_end"],
            "target_chrom": chrom,
            "expected_target_chrom": window["expected_target_chrom"],
            "paf_rows_overlapping_window": str(vals["rows"]),
            "query_overlap_bp_sum": str(vals["sum"]),
            "query_covered_bp_union": str(union_bp(vals["intervals"])),
            "chr3_support": "yes" if chrom == "chr3" and vals["rows"] else "no",
            "paf_path": path,
        })
    if not rows:
        rows.append({
            "event_id": window["event_id"],
            "comparison_id": window["comparison_id"],
            "layer": layer,
            "stage": stage,
            "query_name": window["query_name"],
            "query_start": window["query_start"],
            "query_end": window["query_end"],
            "target_chrom": "NO_OVERLAP",
            "expected_target_chrom": window["expected_target_chrom"],
            "paf_rows_overlapping_window": "0",
            "query_overlap_bp_sum": "0",
            "query_covered_bp_union": "0",
            "chr3_support": "no",
            "paf_path": path,
        })
    return rows


def file_summary(path, comparison_id, layer, stage):
    row = {
        "comparison_id": comparison_id,
        "layer": layer,
        "stage": stage,
        "paf_path": path,
        "exists": "yes" if os.path.exists(path) else "no",
        "bytes": "",
        "sha256": "",
        "gzip_ok": "no",
        "records": "0",
        "query_sequences": "0",
        "target_sequences": "0",
        "query_bp": "0",
        "matches": "0",
        "aln_bp": "0",
        "mean_identity": "0",
    }
    if not os.path.exists(path):
        return row
    row["bytes"] = str(os.path.getsize(path))
    row["sha256"] = sha256(path)
    row["gzip_ok"] = gzip_status(path)
    records = query_bp = matches = aln_bp = 0
    queries = set()
    targets = set()
    for fields in paf_records(path):
        records += 1
        queries.add(fields[0])
        targets.add(fields[5])
        query_bp += max(0, int(fields[3]) - int(fields[2]))
        matches += int(fields[9])
        aln_bp += int(fields[10])
    row.update({
        "records": str(records),
        "query_sequences": str(len(queries)),
        "target_sequences": str(len(targets)),
        "query_bp": str(query_bp),
        "matches": str(matches),
        "aln_bp": str(aln_bp),
        "mean_identity": "%.6f" % (float(matches) / aln_bp) if aln_bp else "0",
    })
    return row


def summarize_answers(support_rows):
    grouped = {}
    for row in support_rows:
        key = (row["event_id"], row["comparison_id"], row["layer"], row["stage"])
        entry = grouped.setdefault(key, {
            "event_id": row["event_id"],
            "comparison_id": row["comparison_id"],
            "layer": row["layer"],
            "stage": row["stage"],
            "query_name": row["query_name"],
            "query_start": row["query_start"],
            "query_end": row["query_end"],
            "expected_target_chrom": row["expected_target_chrom"],
            "chr3_support": "no",
            "chr3_rows": "0",
            "chr3_query_overlap_bp_sum": "0",
            "chr3_query_covered_bp_union": "0",
            "all_target_chromosomes_overlapping_window": [],
            "answer": "no",
            "comparison_to_other_aligners": (
                "Updated wfmash p95 evidence is chr3-positive for the PAN027/PAN028 Fig5 windows; "
                "updated sweepGA/FastGA default evidence is chr3-negative."
            ),
            "note": "",
        })
        chrom = row["target_chrom"]
        if chrom not in entry["all_target_chromosomes_overlapping_window"]:
            entry["all_target_chromosomes_overlapping_window"].append(chrom)
        if chrom == "NO_PAF":
            entry["chr3_support"] = "not_evaluable"
            entry["answer"] = "not_evaluable"
            entry["note"] = "PAF file missing or incomplete for this layer."
        elif chrom == "chr3" and int(row["paf_rows_overlapping_window"]) > 0:
            entry["chr3_support"] = "yes"
            entry["answer"] = "yes"
            entry["chr3_rows"] = row["paf_rows_overlapping_window"]
            entry["chr3_query_overlap_bp_sum"] = row["query_overlap_bp_sum"]
            entry["chr3_query_covered_bp_union"] = row["query_covered_bp_union"]
            entry["note"] = "chr3 target rows overlap the Fig5 query candidate window."
    out = []
    for key in sorted(grouped):
        entry = grouped[key]
        entry["all_target_chromosomes_overlapping_window"] = ",".join(sorted(entry["all_target_chromosomes_overlapping_window"]))
        if not entry["note"]:
            entry["note"] = "No chr3 target rows overlap the Fig5 query candidate window at this layer."
        out.append(entry)
    return out


def write_runtime_diagnosis(file_rows):
    raw_rows = [r for r in file_rows if r["stage"] == "raw"]
    if raw_rows and all(r["exists"] == "yes" and r["gzip_ok"] == "yes" and int(r["records"]) > 0 for r in raw_rows):
        path = os.path.join(PACKAGE_DIR, "summaries", "longrun_runtime_diagnosis.tsv")
        if os.path.exists(path):
            os.unlink(path)
        return
    slurm_path = os.path.join(PACKAGE_DIR, "summaries", "slurm_jobs.tsv")
    slurm_rows = read_tsv(slurm_path) if os.path.exists(slurm_path) else []
    by_cid = {r["comparison_id"]: r for r in raw_rows}
    rows = []
    for job in slurm_rows:
        f = by_cid.get(job["comparison_id"], {})
        rows.append({
            "comparison_id": job["comparison_id"],
            "job_id": job["job_id"],
            "status": job.get("status", ""),
            "requested_time": job.get("requested_time", ""),
            "output_paf": job.get("output_paf", ""),
            "output_exists": f.get("exists", "no"),
            "output_bytes": f.get("bytes", ""),
            "records": f.get("records", "0"),
            "diagnosis": "not complete/evaluable after allowed long-run attempt; do not interpret header-only or unflushed output as chr3-negative",
        })
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "longrun_runtime_diagnosis.tsv"), rows, [
        "comparison_id",
        "job_id",
        "status",
        "requested_time",
        "output_paf",
        "output_exists",
        "output_bytes",
        "records",
        "diagnosis",
    ])


def main():
    windows = read_tsv(CANDIDATE_WINDOWS)
    comparisons = read_tsv(os.path.join(PACKAGE_DIR, "config", "comparisons.tsv"))
    file_rows = []
    support_rows = []
    for comp in comparisons:
        for layer, stage, path in layer_paths_for_comparison(comp["comparison_id"]):
            file_rows.append(file_summary(path, comp["comparison_id"], layer, stage))
    for window in windows:
        for layer, stage, path in layer_paths_for_comparison(window["comparison_id"]):
            support_rows.extend(summarize_window(path, window, layer, stage))

    file_fields = [
        "comparison_id",
        "layer",
        "stage",
        "paf_path",
        "exists",
        "bytes",
        "sha256",
        "gzip_ok",
        "records",
        "query_sequences",
        "target_sequences",
        "query_bp",
        "matches",
        "aln_bp",
        "mean_identity",
    ]
    support_fields = [
        "event_id",
        "comparison_id",
        "layer",
        "stage",
        "query_name",
        "query_start",
        "query_end",
        "target_chrom",
        "expected_target_chrom",
        "paf_rows_overlapping_window",
        "query_overlap_bp_sum",
        "query_covered_bp_union",
        "chr3_support",
        "paf_path",
    ]
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "paf_file_summary.tsv"), file_rows, file_fields)
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "raw_candidate_window_support.tsv"), [r for r in support_rows if r["stage"] == "raw"], support_fields)
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "candidate_window_support.tsv"), support_rows, support_fields)
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "minimap2_chop_sweep_chr3_support_summary.tsv"), summarize_answers(support_rows), [
        "event_id",
        "comparison_id",
        "layer",
        "stage",
        "query_name",
        "query_start",
        "query_end",
        "expected_target_chrom",
        "chr3_support",
        "chr3_rows",
        "chr3_query_overlap_bp_sum",
        "chr3_query_covered_bp_union",
        "all_target_chromosomes_overlapping_window",
        "answer",
        "comparison_to_other_aligners",
        "note",
    ])
    write_runtime_diagnosis(file_rows)


if __name__ == "__main__":
    main()
