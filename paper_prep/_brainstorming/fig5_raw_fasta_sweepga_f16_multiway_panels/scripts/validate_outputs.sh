#!/usr/bin/env bash
set -euo pipefail

PANEL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

python3 - "$PANEL_DIR" <<'PY'
import csv
import sys
from pathlib import Path

panel = Path(sys.argv[1])
required = [
    "multiway_candidate_support.tsv",
    "multiway_candidate_summary.tsv",
    "raw_chr3_chr9_other_support_summary.tsv",
    "input_manifest.tsv",
    "PAR1_XY_positive_control.raw_support.tsv",
    "PAN027_chr9q_chr3q_PHR_candidate.raw_support.tsv",
    "PAN028_chr9q_chr3q_PHR_candidate.raw_support.tsv",
    "fig5_raw_fasta_sweepga_f16_multiway_panels.pdf",
    "fig5_raw_fasta_sweepga_f16_multiway_panels.svg",
    "fig5_raw_fasta_sweepga_f16_multiway_panels.png",
    "preview_png/fig5_raw_fasta_sweepga_f16_multiway_panels.png",
]
missing = [name for name in required if not (panel / name).exists() or (panel / name).stat().st_size == 0]
if missing:
    raise SystemExit("missing/empty outputs: " + ", ".join(missing))

events = {
    "PAR1_XY_positive_control",
    "PAN027_chr9q_chr3q_PHR_candidate",
    "PAN028_chr9q_chr3q_PHR_candidate",
}
with open(panel / "multiway_candidate_summary.tsv", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
raw = [r for r in rows if r["source_layer"] == "raw_many_many_whole_genome"]
if {r["event_id"] for r in raw} != events:
    raise SystemExit("raw summary does not cover exactly the three panel events")
for row in raw:
    if int(row["row_count"]) <= 0:
        raise SystemExit(f"raw row_count is zero for {row['event_id']}")
    if int(row["query_union_coverage_bp"]) <= 0:
        raise SystemExit(f"raw query union coverage is zero for {row['event_id']}")
    if row["event_id"].startswith("PAN0") and int(row["chr3_query_union_coverage_bp"]) <= 0:
        raise SystemExit(f"missing raw chr3 support for {row['event_id']}")

with open(panel / "input_manifest.tsv", newline="") as handle:
    manifest = list(csv.DictReader(handle, delimiter="\t"))
raw_manifest = [r for r in manifest if r["source_layer"] == "raw_many_many_whole_genome"]
if not raw_manifest:
    raise SystemExit("manifest lacks raw many:many source rows")
for row in raw_manifest:
    if not row["path"].endswith(".sweepga_frequency16_many_many_j0.paf.gz"):
        raise SystemExit("raw manifest path is not the whole-genome f16 many:many PAF: " + row["path"])
    if row["role"] != "source_of_truth_multiway_input":
        raise SystemExit("raw manifest role does not mark source-of-truth input")
print("validated fig5 f16 multiway panel outputs")
PY
