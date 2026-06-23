#!/usr/bin/env bash
set -euo pipefail

PANEL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

for f in \
    "$PANEL_DIR/fig5_raw_fasta_sweepga_f16_chopped_panels.pdf" \
    "$PANEL_DIR/fig5_raw_fasta_sweepga_f16_chopped_panels.svg" \
    "$PANEL_DIR/raw_fasta_chopped_panel_segments.tsv" \
    "$PANEL_DIR/raw_fasta_chopped_panel_summary.tsv" \
    "$PANEL_DIR/slurm_jobs.tsv"; do
    [[ -s "$f" ]] || { echo "missing or empty: $f" >&2; exit 1; }
done

python3 - "$PANEL_DIR" <<'PY'
import csv
import os
import sys

panel_dir = sys.argv[1]
summary_path = os.path.join(panel_dir, "raw_fasta_chopped_panel_summary.tsv")
jobs_path = os.path.join(panel_dir, "slurm_jobs.tsv")
segments_path = os.path.join(panel_dir, "raw_fasta_chopped_panel_segments.tsv")

required = {
    "PAR1_XY_positive_control",
    "PAN027_chr9q_chr3q_PHR_candidate",
    "PAN028_chr9q_chr3q_PHR_candidate",
}
with open(summary_path, newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
seen = {row["event_id"] for row in rows}
missing = required - seen
if missing:
    raise SystemExit(f"missing summary events: {sorted(missing)}")
bad = [row for row in rows if row["status"] != "OK" or int(row["expected_target_rows"]) <= 0]
if bad:
    raise SystemExit(f"events without expected target support: {bad}")

with open(jobs_path, newline="") as handle:
    jobs = list(csv.DictReader(handle, delimiter="\t"))
if not jobs:
    raise SystemExit("slurm_jobs.tsv has no rows")
for row in jobs:
    if not row["panel_slurm_job_id"] or row["panel_slurm_job_id"] == "manual":
        raise SystemExit("panel provenance does not contain a Slurm job id")
    if not row["raw_paf"].endswith(".sweepga_frequency16_many_many_j0.paf.gz"):
        raise SystemExit(f"raw provenance does not point at f16 raw PAF: {row['raw_paf']}")
    if "window_extract.one_one_scoring_ani_chopped_l2000_o0" not in row["filtered_paf_1to1_scoring_ani"]:
        raise SystemExit("filtered provenance does not point at 2 kb 1:1 ANI PAF")
    if not row["raw_window_extract_paf"].endswith(".raw_window_extract.paf.gz"):
        raise SystemExit("missing raw window-extract PAF provenance")
    if not row.get("render_slurm_job_id"):
        raise SystemExit("missing render Slurm job id")

with open(segments_path, newline="") as handle:
    segments = list(csv.DictReader(handle, delimiter="\t"))
if not segments:
    raise SystemExit("segments file has no rows")
if any("untangle" in row["source_layer"].lower() or "untangle" in row["filtered_paf"].lower() for row in segments):
    raise SystemExit("segments contain untangle-derived provenance")
PY

grep -q "raw FASTA SweepGA/FastGA f16" "$PANEL_DIR/README.md"
! git -C "$PANEL_DIR/../../.." diff --name-only -- submission/ | grep -q .
