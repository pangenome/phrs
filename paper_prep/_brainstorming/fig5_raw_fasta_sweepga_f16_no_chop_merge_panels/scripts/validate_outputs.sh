#!/usr/bin/env bash
set -euo pipefail

panel_dir="${1:-paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_no_chop_merge_panels}"

required=(
  "$panel_dir/raw_merge_panel_segments.tsv"
  "$panel_dir/raw_merge_panel_summary.tsv"
  "$panel_dir/raw_merge_panel_manifest.tsv"
  "$panel_dir/fig5_raw_fasta_sweepga_f16_no_chop_merge_panels.pdf"
  "$panel_dir/fig5_raw_fasta_sweepga_f16_no_chop_merge_panels.svg"
)

for path in "${required[@]}"; do
  test -s "$path" || {
    echo "missing or empty output: $path" >&2
    exit 1
  }
done

python3 - "$panel_dir/raw_merge_panel_summary.tsv" <<'PY'
import csv
import sys

path = sys.argv[1]
rows = list(csv.DictReader(open(path, newline=""), delimiter="\t"))
expected_modes = {"no_merge_ani", "merge50k_ani", "merge50k_log_length_ani"}
expected_events = {
    "PAR1_XY_positive_control",
    "PAN027_chr9q_chr3q_PHR_candidate",
    "PAN028_chr9q_chr3q_PHR_candidate",
}
modes = {row["filter_mode"] for row in rows}
events = {row["event_id"] for row in rows}
if modes != expected_modes:
    raise SystemExit(f"unexpected filter modes: {sorted(modes)}")
if events != expected_events:
    raise SystemExit(f"unexpected events: {sorted(events)}")
if len(rows) != len(expected_modes) * len(expected_events):
    raise SystemExit(f"unexpected summary row count: {len(rows)}")
chr3_rows = [
    row for row in rows
    if row["event_id"].endswith("PHR_candidate") and row["expected_target_chroms"] == "chr3"
]
if not chr3_rows:
    raise SystemExit("no chr3 candidate rows in summary")
for row in chr3_rows:
    int(row["segment_rows"])
    int(row["expected_target_rows"])
    int(row["union_expected_overlap_bp"])
    if row["chr3_survival_status"] not in {"CHR3_SURVIVES", "CHR3_DROPPED"}:
        raise SystemExit(f"bad chr3 survival status: {row['chr3_survival_status']}")
print(f"validated {len(rows)} summary rows across {len(expected_modes)} no-chop modes")
PY
