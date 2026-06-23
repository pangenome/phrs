#!/usr/bin/env bash
set -euo pipefail

panel_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
summary="$panel_dir/chop_filter_panel_summary.tsv"
segments="$panel_dir/chop_filter_panel_segments.tsv"
pdf="$panel_dir/fig5_raw_fasta_sweepga_f16_chop_filter_sensitivity_panels.pdf"

test -s "$summary"
test -s "$segments"
test -s "$pdf"

python3 - "$summary" "$segments" <<'PY'
import csv
import sys

summary_path, segments_path = sys.argv[1:3]
with open(summary_path, newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

expected_rows = 2 * 3 * 4
if len(rows) != expected_rows:
    raise SystemExit(f"expected {expected_rows} summary rows, found {len(rows)}")

bad = [r for r in rows if r["chr3_survives"] != "yes" or int(r["chr3_query_union_bp"]) <= 0]
if bad:
    labels = ", ".join(f'{r["event_id"]}/{r["chop_length_bp"]}/{r["filter_mode"]}' for r in bad)
    raise SystemExit(f"chr3 failed to survive in: {labels}")

with open(segments_path, newline="") as handle:
    segs = list(csv.DictReader(handle, delimiter="\t"))
if not segs:
    raise SystemExit("no segment rows")
for row in segs:
    qs = int(row["query_start_abs"])
    qe = int(row["query_end_abs"])
    ws = int(row["window_start_abs"])
    we = int(row["window_end_abs"])
    if qs < ws or qe > we or qe <= qs:
        raise SystemExit(f"segment outside absolute query window: {row}")

print(f"validated {len(rows)} summary rows and {len(segs)} segment rows")
PY
