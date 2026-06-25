#!/usr/bin/env bash
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
prefix="fig5_whole_genome_length_scaled_tracks"

for path in \
    "$here/$prefix.pdf" \
    "$here/length_scaled_track_segments.tsv" \
    "$here/length_scaled_track_chromosomes.tsv" \
    "$here/length_scaled_track_summary.tsv" \
    "$here/length_scaled_track_manifest.tsv"; do
    [[ -s "$path" ]] || { echo "missing or empty: $path" >&2; exit 1; }
done

VALIDATE_DIR="$here" python3 - <<'PY'
import csv
import os
from pathlib import Path

here = Path(os.environ["VALIDATE_DIR"])
segments = list(csv.DictReader((here / "length_scaled_track_segments.tsv").open(), delimiter="\t"))
chroms = list(csv.DictReader((here / "length_scaled_track_chromosomes.tsv").open(), delimiter="\t"))

pages = {
    "PAN027_hap1_maternal_vs_PAN010",
    "PAN027_hap2_paternal_vs_PAN011",
    "PAN028_hap1_maternal_vs_PAN027",
}
methods = {
    "untangle_strict_primary_path",
    "wfmash_p95_qgrid2kb_1to1_ani",
    "sweepga_fastga_f32_qgrid2kb_1to1_ani",
}
seen = {(row["page_id"], row["method_id"]) for row in segments}
missing = sorted((page, method) for page in pages for method in methods if (page, method) not in seen)
if missing:
    raise SystemExit(f"missing page/method combinations: {missing}")

chrom_counts = {}
for row in chroms:
    if int(float(row["query_chrom_length"])) > 0:
        chrom_counts[row["page_id"]] = chrom_counts.get(row["page_id"], 0) + 1
bad = {page: count for page, count in chrom_counts.items() if count < 23}
if bad:
    raise SystemExit(f"expected at least 23 nonzero query chromosomes per page, got {bad}")

for page in pages:
    for method in methods:
        inter = [
            row for row in segments
            if row["page_id"] == page
            and row["method_id"] == method
            and row["display_state"] == "interchromosomal"
        ]
        if not inter:
            raise SystemExit(f"no interchromosomal rows for {page} {method}")

bad_binned = [
    row for row in segments
    if row["method_id"] in {"wfmash_p95_qgrid2kb_1to1_ani", "sweepga_fastga_f32_qgrid2kb_1to1_ani"}
    and row["source_row_type"] != "filtered_paf_record_exact_query_interval"
]
if bad_binned:
    raise SystemExit(f"PAF-derived rows must be exact PAF records, saw {bad_binned[0]['source_row_type']}")

for row in segments:
    if row["source_row_type"] == "filtered_paf_record_exact_query_interval":
        if not row["source_file"] or not row["paf_query_name"] or not row["paf_target_name"]:
            raise SystemExit("exact PAF row missing source_file/paf_query_name/paf_target_name")
        length = int(float(row["segment_end"])) - int(float(row["segment_start"]))
        if length != int(float(row["segment_length_bp"])):
            raise SystemExit("segment length does not match exact PAF query interval")

for event in ("PAN027_chr9q_chr3q_PHR_candidate", "PAN028_chr9q_chr3q_PHR_candidate"):
    hits = [
        row for row in segments
        if row["callout_event_id"] == event
        and row["display_target_chrom"] == "chr3"
        and row["display_state"] == "interchromosomal"
    ]
    if not hits:
        raise SystemExit(f"missing exact chr3 support in callout {event}")

print(f"validated {len(segments)} segments across {len(pages)} pages and {len(methods)} methods")
PY

for page in PAN027_hap1_maternal_vs_PAN010 PAN027_hap2_paternal_vs_PAN011 PAN028_hap1_maternal_vs_PAN027; do
    [[ -s "$here/$prefix.$page.pdf" ]] || { echo "missing page pdf: $page" >&2; exit 1; }
    [[ -s "$here/$prefix.$page.png" ]] || { echo "missing page png: $page" >&2; exit 1; }
    [[ -s "$here/$prefix.$page.svg" ]] || { echo "missing page svg: $page" >&2; exit 1; }
done

echo "OK: length-scaled whole-genome track outputs validated"
