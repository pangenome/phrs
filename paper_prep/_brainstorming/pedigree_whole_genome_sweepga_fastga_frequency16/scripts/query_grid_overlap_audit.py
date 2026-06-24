#!/usr/bin/env python3
import csv
import gzip
import os
import re
import shlex
import subprocess
from collections import Counter, defaultdict
from contextlib import contextmanager

from common import PACKAGE_DIR, write_tsv


MANIFEST = os.path.join(PACKAGE_DIR, "summaries", "query_grid_chop_filter_manifest.tsv")
STATUS_DIR = os.path.join(PACKAGE_DIR, "summaries", "query_grid_chop_filter_status")
CANDIDATES = os.path.join(PACKAGE_DIR, "config", "candidate_windows.tsv")
ROW_START_SEGMENTS = os.path.join(
    PACKAGE_DIR,
    "..",
    "fig5_raw_fasta_sweepga_f16_chop_filter_sensitivity_panels",
    "chop_filter_panel_segments.tsv",
)
OUT = os.path.join(PACKAGE_DIR, "summaries", "query_grid_overlap_audit.tsv")
EXAMPLES_OUT = os.path.join(PACKAGE_DIR, "summaries", "query_grid_overlap_audit_examples.tsv")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@contextmanager
def open_text(path):
    if path.endswith(".gz"):
        try:
            proc = subprocess.Popen(["pigz", "-dc", path], stdout=subprocess.PIPE, text=True)
        except OSError:
            with gzip.open(path, "rt") as fh:
                yield fh
        else:
            try:
                yield proc.stdout
            finally:
                if proc.stdout:
                    proc.stdout.close()
                proc.wait()
        return
    with open(path) as fh:
        yield fh


def parse_tags(fields):
    tags = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_target(target_name):
    hap = "unknown"
    chrom = target_name
    if "#" in target_name:
        tail = target_name.rsplit("#", 1)[-1]
    else:
        tail = target_name
    match = re.match(r"(h[12])_(chr[^_]+)$", tail)
    if match:
        hap, chrom = match.group(1), match.group(2)
    elif "_" in tail:
        chrom = tail.rsplit("_", 1)[-1]
    return chrom, hap


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


def max_depth_and_overlap(intervals):
    events = []
    for start, end in intervals:
        if end > start:
            events.append((start, 1))
            events.append((end, -1))
    if not events:
        return 0, 0
    events.sort(key=lambda item: (item[0], item[1]))
    depth = 0
    max_depth = 0
    overlap_bp = 0
    prev = None
    for pos, delta in events:
        if prev is not None and pos > prev and depth > 1:
            overlap_bp += pos - prev
        depth += delta
        max_depth = max(max_depth, depth)
        prev = pos
    return max_depth, overlap_bp


def boundary_audit(path, length):
    if not os.path.exists(path):
        return {
            "rows": 0,
            "violations": 0,
            "non_grid_starts": 0,
            "non_grid_ends": 0,
            "example": "",
            "status": "missing",
        }
    awk = r'''
    BEGIN{FS="\t"; OFS="\t"}
    {
      rows++;
      zs=""; ze=""; zm="";
      start_mod=(($3 % L)!=0);
      end_mod=(($4 % L)!=0);
      if (start_mod) non_grid_starts++;
      if (end_mod) non_grid_ends++;
      if (start_mod || end_mod) {
        for (i=13; i<=NF; i++) {
          if ($i ~ /^zs:i:/) { zs=substr($i,6) }
          else if ($i ~ /^ze:i:/) { ze=substr($i,6) }
          else if ($i ~ /^zm:Z:/) { zm=substr($i,6) }
        }
      }
      start_ok=(!start_mod || (zs != "" && $3 == zs));
      end_ok=(!end_mod || (ze != "" && $4 == ze));
      if (!start_ok || !end_ok) {
        violations++;
        if (example == "") {
          example=$1 ":" $3 "-" $4 ";target=" $6 ":" $8 "-" $9 ";zm=" zm ";zs=" zs ";ze=" ze;
        }
      }
    }
    END{print rows+0, violations+0, non_grid_starts+0, non_grid_ends+0, example}
    '''
    if path.endswith(".gz"):
        cmd = f"pigz -dc {shlex.quote(path)} | awk -v L={length} {shlex.quote(awk)}"
    else:
        cmd = f"awk -v L={length} {shlex.quote(awk)} {shlex.quote(path)}"
    out = subprocess.check_output(cmd, shell=True, text=True)
    parts = out.rstrip("\n").split("\t", 4)
    rows = int(parts[0])
    violations = int(parts[1])
    non_grid_starts = int(parts[2])
    non_grid_ends = int(parts[3])
    example = parts[4] if len(parts) > 4 else ""
    return {
        "rows": rows,
        "violations": violations,
        "non_grid_starts": non_grid_starts,
        "non_grid_ends": non_grid_ends,
        "example": example,
        "status": "ok",
    }


def summarize_filtered_candidate(path, candidate):
    empty = {
        "post_filter_rows": 0,
        "sum_overlap_bp": 0,
        "query_union_bp": 0,
        "query_redundant_bp": 0,
        "query_overlap_max_depth": 0,
        "target_groups": "",
        "chr3_rows": 0,
        "chr3_sum_overlap_bp": 0,
        "chr3_query_union_bp": 0,
        "chr3_query_redundant_bp": 0,
        "chr3_hap_counts": "",
        "chr3_adjacent_hap_switches": 0,
        "chr3_overlap_hap_switch_pairs": 0,
        "alternating_h1_h2_offset_persists": "no",
        "chr3_retained_chunks": "",
        "example_rows": [],
    }
    if not candidate or not os.path.exists(path):
        return empty
    qname = candidate["query_name"]
    win_start = int(candidate["query_start"])
    win_end = int(candidate["query_end"])
    expected_chrom = candidate["expected_target_chrom"]
    intervals = []
    rows = []
    group_intervals = defaultdict(list)
    group_sum = Counter()
    group_rows = Counter()
    chr3_rows = []
    awk = 'BEGIN{FS="\t"} $1 == Q && $3 < E && $4 > S {print NR "\t" $0}'
    if path.endswith(".gz"):
        cmd = f"pigz -dc {shlex.quote(path)} | awk -v Q={shlex.quote(qname)} -v S={win_start} -v E={win_end} {shlex.quote(awk)}"
    else:
        cmd = f"awk -v Q={shlex.quote(qname)} -v S={win_start} -v E={win_end} {shlex.quote(awk)} {shlex.quote(path)}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            line_no_text, paf_line = line.rstrip("\n").split("\t", 1)
            line_no = int(line_no_text)
            fields = paf_line.split("\t")
            if len(fields) < 12 or fields[0] != qname:
                continue
            q_start = int(fields[2])
            q_end = int(fields[3])
            ov_start = max(q_start, win_start)
            ov_end = min(q_end, win_end)
            if ov_end <= ov_start:
                continue
            target_chrom, target_hap = parse_target(fields[5])
            target_start = int(fields[7])
            target_end = int(fields[8])
            interval = (ov_start, ov_end)
            key = f"{target_chrom}/{target_hap}"
            ov = ov_end - ov_start
            rows.append((ov_start, ov_end, key, target_hap, target_chrom, fields, line_no))
            intervals.append(interval)
            group_intervals[key].append(interval)
            group_sum[key] += ov
            group_rows[key] += 1
            if target_chrom == expected_chrom:
                chr3_rows.append((ov_start, ov_end, target_hap, fields[5], target_start, target_end, fields[4], line_no))
    finally:
        if proc.stdout:
            proc.stdout.close()
        proc.wait()
    if not rows:
        return empty
    all_sum = sum(end - start for start, end in intervals)
    all_union = union_bp(intervals)
    max_depth, redundant_bp = max_depth_and_overlap(intervals)
    group_bits = []
    for key in sorted(group_rows):
        group_bits.append(
            f"{key}:rows={group_rows[key]},sum={group_sum[key]},union={union_bp(group_intervals[key])}"
        )
    chr3_intervals = [(row[0], row[1]) for row in chr3_rows]
    chr3_sum = sum(end - start for start, end in chr3_intervals)
    chr3_union = union_bp(chr3_intervals)
    chr3_max_depth, chr3_redundant = max_depth_and_overlap(chr3_intervals)
    hap_counts = Counter(row[2] for row in chr3_rows)
    chr3_sorted = sorted(chr3_rows)
    adjacent_switches = 0
    overlap_switches = 0
    for prev, cur in zip(chr3_sorted, chr3_sorted[1:]):
        if prev[2] != cur[2]:
            adjacent_switches += 1
            if cur[0] < prev[1]:
                overlap_switches += 1
    retained = []
    for ov_start, ov_end, hap, target_name, target_start, target_end, strand, _line_no in chr3_sorted:
        retained.append(f"{ov_start}-{ov_end}->{target_name}:{target_start}-{target_end}({strand})")
    examples = []
    for ov_start, ov_end, key, _hap, _chrom, fields, line_no in sorted(rows)[:8]:
        examples.append({
            "source": "query_grid_filtered",
            "line_no": str(line_no),
            "query_interval": f"{ov_start}-{ov_end}",
            "target_group": key,
            "target_interval": f"{fields[5]}:{fields[7]}-{fields[8]}",
            "strand": fields[4],
            "matches": fields[9],
            "alignment_length": fields[10],
            "mapq": fields[11],
        })
    return {
        "post_filter_rows": len(rows),
        "sum_overlap_bp": all_sum,
        "query_union_bp": all_union,
        "query_redundant_bp": all_sum - all_union,
        "query_overlap_max_depth": max_depth,
        "target_groups": ";".join(group_bits),
        "chr3_rows": len(chr3_rows),
        "chr3_sum_overlap_bp": chr3_sum,
        "chr3_query_union_bp": chr3_union,
        "chr3_query_redundant_bp": chr3_sum - chr3_union,
        "chr3_hap_counts": ";".join(f"{hap}:{hap_counts[hap]}" for hap in sorted(hap_counts)),
        "chr3_adjacent_hap_switches": adjacent_switches,
        "chr3_overlap_hap_switch_pairs": overlap_switches,
        "alternating_h1_h2_offset_persists": "yes" if overlap_switches > 0 and chr3_redundant > 0 else "no",
        "chr3_retained_chunks": ";".join(retained),
        "example_rows": examples,
    }


def row_start_summary():
    by_key = {}
    examples = []
    if not os.path.exists(ROW_START_SEGMENTS):
        return by_key, examples
    with open(ROW_START_SEGMENTS, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows_by_key = defaultdict(list)
        for row in reader:
            if row["filter_mode"] != "no_merge_ani":
                continue
            key = (row["event_id"], row["comparison_id"], int(row["chop_length_bp"]))
            rows_by_key[key].append(row)
    for key, rows in rows_by_key.items():
        intervals = []
        chr3_intervals = []
        chr3_rows = []
        for row in rows:
            interval = (int(row["query_clip_start"]), int(row["query_clip_end"]))
            intervals.append(interval)
            if row["target_chrom"] == "chr3":
                chr3_intervals.append(interval)
                chrom, hap = parse_target(row["target_name"])
                chr3_rows.append((interval[0], interval[1], hap, row))
        all_sum = sum(end - start for start, end in intervals)
        chr3_sum = sum(end - start for start, end in chr3_intervals)
        chr3_sorted = sorted(chr3_rows)
        switches = 0
        overlap_switches = 0
        for prev, cur in zip(chr3_sorted, chr3_sorted[1:]):
            if prev[2] != cur[2]:
                switches += 1
                if cur[0] < prev[1]:
                    overlap_switches += 1
        by_key[key] = {
            "row_start_rows": len(rows),
            "row_start_sum_overlap_bp": all_sum,
            "row_start_query_union_bp": union_bp(intervals),
            "row_start_chr3_rows": len(chr3_rows),
            "row_start_chr3_sum_overlap_bp": chr3_sum,
            "row_start_chr3_query_union_bp": union_bp(chr3_intervals),
            "row_start_chr3_adjacent_hap_switches": switches,
            "row_start_chr3_overlap_hap_switch_pairs": overlap_switches,
        }
        for ov_start, ov_end, hap, row in chr3_sorted[:8]:
            examples.append({
                "event_id": key[0],
                "comparison_id": key[1],
                "chop_length_bp": str(key[2]),
                "source": "row_start_filtered_no_merge_ani",
                "line_no": row["filtered_paf_line"],
                "query_interval": f"{ov_start}-{ov_end}",
                "target_group": f"chr3/{hap}",
                "target_interval": f"{row['target_name']}:{row['target_start']}-{row['target_end']}",
                "strand": row["strand"],
                "matches": row["matches"],
                "alignment_length": row["alignment_length"],
                "mapq": "",
            })
    return by_key, examples


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    manifest_rows = read_tsv(MANIFEST)
    candidates = {row["comparison_id"]: row for row in read_tsv(CANDIDATES)}
    row_start_by_key, row_start_examples = row_start_summary()
    rows = []
    examples = []
    completed_keys = set()
    for manifest in manifest_rows:
        comparison_id = manifest["comparison_id"]
        length = int(manifest["chop_length_bp"])
        completed_keys.add((comparison_id, length))
        candidate = candidates.get(comparison_id)
        event_id = candidate["event_id"] if candidate else "NO_FIG5_CHR9_CHR3_CANDIDATE"
        boundary = boundary_audit(manifest["chopped_paf"], length)
        filtered = summarize_filtered_candidate(manifest["filtered_paf"], candidate)
        row_start = row_start_by_key.get((event_id, comparison_id, length), {})
        conclusion = "boundary_only"
        if candidate:
            if filtered["alternating_h1_h2_offset_persists"] == "yes":
                conclusion = "filtering_still_ambiguous"
            elif filtered["chr3_rows"] and filtered["query_redundant_bp"] == 0 and filtered["chr3_query_redundant_bp"] == 0:
                conclusion = "chopper_fixed_offset_overlap_artifact"
            elif filtered["query_redundant_bp"] > 0:
                conclusion = "remaining_query_overlap_issue"
            else:
                conclusion = "no_chr3_after_filter"
        rows.append({
            "event_id": event_id,
            "comparison_id": comparison_id,
            "chop_length_bp": str(length),
            "artifact_status": manifest["status"],
            "chopped_paf": manifest["chopped_paf"],
            "filtered_paf": manifest["filtered_paf"],
            "boundary_rows": str(boundary["rows"]),
            "boundary_status": boundary["status"],
            "boundary_violation_rows": str(boundary["violations"]),
            "non_grid_start_rows": str(boundary["non_grid_starts"]),
            "non_grid_end_rows": str(boundary["non_grid_ends"]),
            "boundary_violation_example": boundary["example"],
            "post_filter_rows_in_window": str(filtered["post_filter_rows"]),
            "sum_overlap_bp": str(filtered["sum_overlap_bp"]),
            "query_union_bp": str(filtered["query_union_bp"]),
            "query_redundant_bp": str(filtered["query_redundant_bp"]),
            "query_overlap_max_depth": str(filtered["query_overlap_max_depth"]),
            "target_chrom_haplotype_metrics": filtered["target_groups"],
            "chr3_retained_rows": str(filtered["chr3_rows"]),
            "chr3_sum_overlap_bp": str(filtered["chr3_sum_overlap_bp"]),
            "chr3_query_union_bp": str(filtered["chr3_query_union_bp"]),
            "chr3_query_redundant_bp": str(filtered["chr3_query_redundant_bp"]),
            "chr3_haplotype_counts": filtered["chr3_hap_counts"],
            "chr3_adjacent_hap_switches": str(filtered["chr3_adjacent_hap_switches"]),
            "chr3_overlap_hap_switch_pairs": str(filtered["chr3_overlap_hap_switch_pairs"]),
            "alternating_h1_h2_offset_persists": filtered["alternating_h1_h2_offset_persists"],
            "chr3_retained_chunks": filtered["chr3_retained_chunks"],
            "row_start_rows": str(row_start.get("row_start_rows", "")),
            "row_start_sum_overlap_bp": str(row_start.get("row_start_sum_overlap_bp", "")),
            "row_start_query_union_bp": str(row_start.get("row_start_query_union_bp", "")),
            "row_start_chr3_rows": str(row_start.get("row_start_chr3_rows", "")),
            "row_start_chr3_sum_overlap_bp": str(row_start.get("row_start_chr3_sum_overlap_bp", "")),
            "row_start_chr3_query_union_bp": str(row_start.get("row_start_chr3_query_union_bp", "")),
            "row_start_chr3_adjacent_hap_switches": str(row_start.get("row_start_chr3_adjacent_hap_switches", "")),
            "row_start_chr3_overlap_hap_switch_pairs": str(row_start.get("row_start_chr3_overlap_hap_switch_pairs", "")),
            "audit_conclusion": conclusion,
        })
        for example in filtered.pop("example_rows"):
            example.update({
                "event_id": event_id,
                "comparison_id": comparison_id,
                "chop_length_bp": str(length),
            })
            examples.append(example)
    for status_name in sorted(os.listdir(STATUS_DIR)):
        if not status_name.endswith(".tsv"):
            continue
        for status in read_tsv(os.path.join(STATUS_DIR, status_name)):
            key = (status["comparison_id"], int(status["chop_length_bp"]))
            if key in completed_keys:
                continue
            candidate = candidates.get(status["comparison_id"])
            rows.append({
                "event_id": candidate["event_id"] if candidate else "NO_FIG5_CHR9_CHR3_CANDIDATE",
                "comparison_id": status["comparison_id"],
                "chop_length_bp": status["chop_length_bp"],
                "artifact_status": f"status_only_{status.get('status', 'missing')}",
                "chopped_paf": status.get("chopped_paf", ""),
                "filtered_paf": status.get("filtered_paf", ""),
                "boundary_rows": "0",
                "boundary_status": "not_audited_missing_manifest_output",
                "boundary_violation_rows": "",
                "non_grid_start_rows": "",
                "non_grid_end_rows": "",
                "boundary_violation_example": "",
                "post_filter_rows_in_window": "",
                "sum_overlap_bp": "",
                "query_union_bp": "",
                "query_redundant_bp": "",
                "query_overlap_max_depth": "",
                "target_chrom_haplotype_metrics": "",
                "chr3_retained_rows": "",
                "chr3_sum_overlap_bp": "",
                "chr3_query_union_bp": "",
                "chr3_query_redundant_bp": "",
                "chr3_haplotype_counts": "",
                "chr3_adjacent_hap_switches": "",
                "chr3_overlap_hap_switch_pairs": "",
                "alternating_h1_h2_offset_persists": "",
                "chr3_retained_chunks": "",
                "row_start_rows": "",
                "row_start_sum_overlap_bp": "",
                "row_start_query_union_bp": "",
                "row_start_chr3_rows": "",
                "row_start_chr3_sum_overlap_bp": "",
                "row_start_chr3_query_union_bp": "",
                "row_start_chr3_adjacent_hap_switches": "",
                "row_start_chr3_overlap_hap_switch_pairs": "",
                "audit_conclusion": "missing_query_grid_output_reported_not_rerun",
            })
    fields = [
        "event_id", "comparison_id", "chop_length_bp", "artifact_status",
        "chopped_paf", "filtered_paf", "boundary_rows", "boundary_status",
        "boundary_violation_rows", "non_grid_start_rows", "non_grid_end_rows",
        "boundary_violation_example", "post_filter_rows_in_window", "sum_overlap_bp",
        "query_union_bp", "query_redundant_bp", "query_overlap_max_depth",
        "target_chrom_haplotype_metrics", "chr3_retained_rows", "chr3_sum_overlap_bp",
        "chr3_query_union_bp", "chr3_query_redundant_bp", "chr3_haplotype_counts",
        "chr3_adjacent_hap_switches", "chr3_overlap_hap_switch_pairs",
        "alternating_h1_h2_offset_persists", "chr3_retained_chunks",
        "row_start_rows", "row_start_sum_overlap_bp", "row_start_query_union_bp",
        "row_start_chr3_rows", "row_start_chr3_sum_overlap_bp",
        "row_start_chr3_query_union_bp", "row_start_chr3_adjacent_hap_switches",
        "row_start_chr3_overlap_hap_switch_pairs", "audit_conclusion",
    ]
    rows.sort(key=lambda row: (row["event_id"], row["comparison_id"], int(row["chop_length_bp"])))
    write_tsv(OUT, rows, fields)
    example_fields = [
        "event_id", "comparison_id", "chop_length_bp", "source", "line_no",
        "query_interval", "target_group", "target_interval", "strand",
        "matches", "alignment_length", "mapq",
    ]
    examples.extend(row_start_examples)
    examples.sort(key=lambda row: (row["event_id"], row["comparison_id"], int(row["chop_length_bp"]), row["source"], row["query_interval"]))
    write_tsv(EXAMPLES_OUT, examples, example_fields)
    print(OUT)
    print(EXAMPLES_OUT)


if __name__ == "__main__":
    main()
