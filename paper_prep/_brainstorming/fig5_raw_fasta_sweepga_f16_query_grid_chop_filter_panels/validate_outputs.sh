#!/usr/bin/env bash
set -euo pipefail

panel_dir="${1:-$(cd "$(dirname "$0")" && pwd)}"
prefix="fig5_raw_fasta_sweepga_f16_query_grid_chop_filter_panels"

required=(
  "$prefix.pdf"
  "$prefix.png"
  "$prefix.svg"
  "query_grid_panel_segments.tsv"
  "query_grid_panel_summary.tsv"
  "query_grid_panel_manifest.tsv"
  "README.md"
)

for rel in "${required[@]}"; do
  path="$panel_dir/$rel"
  if [[ ! -s "$path" ]]; then
    echo "missing or empty: $path" >&2
    exit 1
  fi
done

python3 - "$panel_dir" <<'PY'
import csv
import os
import sys

panel_dir = sys.argv[1]
summary_path = os.path.join(panel_dir, "query_grid_panel_summary.tsv")
segments_path = os.path.join(panel_dir, "query_grid_panel_segments.tsv")
manifest_path = os.path.join(panel_dir, "query_grid_panel_manifest.tsv")

with open(summary_path, newline="") as handle:
    summary = list(csv.DictReader(handle, delimiter="\t"))
with open(segments_path, newline="") as handle:
    segments = list(csv.DictReader(handle, delimiter="\t"))
with open(manifest_path, newline="") as handle:
    manifest = list(csv.DictReader(handle, delimiter="\t"))

events = {
    "PAR1_XY_positive_control",
    "PAN027_chr9q_chr3q_PHR_candidate",
    "PAN028_chr9q_chr3q_PHR_candidate",
}
lengths = {"10000", "5000", "2000"}

if len(summary) != len(events) * len(lengths):
    raise SystemExit(f"expected 9 summary rows, observed {len(summary)}")
if len(manifest) != 3 * len(lengths):
    raise SystemExit(f"expected 9 manifest rows, observed {len(manifest)}")
if not segments:
    raise SystemExit("no segment rows")

summary_keys = {(r["event_id"], r["chop_length_bp"]) for r in summary}
missing = {(event, length) for event in events for length in lengths} - summary_keys
if missing:
    raise SystemExit(f"missing summary rows: {sorted(missing)}")

bad_lengths = {r["chop_length_bp"] for r in summary} - lengths
if bad_lengths:
    raise SystemExit(f"unexpected chop lengths: {sorted(bad_lengths)}")
if "1000" in {r["chop_length_bp"] for r in summary + segments + manifest}:
    raise SystemExit("1 kb rows must not be present")

for row in summary:
    if row["status"] != "OK":
        raise SystemExit(f"non-OK row: {row['event_id']} {row['chop_length_bp']} {row['status']}")
    if row["filter_id"] != "one_to_one_ani_o0" or row["chunk_mode"] != "query-grid":
        raise SystemExit(f"unexpected filter/chunk metadata: {row}")
    if row["scaffold_jump"] != "0" or row["scoring"] != "ani" or row["filter_overlap"] != "0":
        raise SystemExit(f"unexpected SweepGA filter metadata: {row}")

for row in segments:
    q0 = int(row["query_start"])
    q1 = int(row["query_end"])
    c0 = int(row["query_clip_start"])
    c1 = int(row["query_clip_end"])
    w0 = int(row["window_start"])
    w1 = int(row["window_end"])
    if not (w0 <= c0 < c1 <= w1):
        raise SystemExit(f"clip outside window: {row}")
    if row["query_grid_mode"] != "query-grid":
        raise SystemExit(f"segment is not tagged query-grid: {row}")
    if q0 == c0 and q1 == c1 and w0 == 0 and w1 == 500000:
        continue
    if row["event_id"] != "PAR1_XY_positive_control" and q0 < 1_000_000:
        raise SystemExit(f"candidate query coordinate does not look genomic: {row}")

print("validated query-grid panel outputs")
PY
