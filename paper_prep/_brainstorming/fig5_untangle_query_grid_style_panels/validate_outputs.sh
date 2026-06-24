#!/usr/bin/env bash
set -euo pipefail

panel_dir="${1:-$(cd "$(dirname "$0")" && pwd)}"
prefix="fig5_untangle_query_grid_style_panels"

required=(
  "$prefix.pdf"
  "$prefix.png"
  "$prefix.svg"
  "untangle_panel_segments.tsv"
  "untangle_panel_summary.tsv"
  "untangle_panel_manifest.tsv"
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
import xml.etree.ElementTree as ET

panel_dir = sys.argv[1]
prefix = "fig5_untangle_query_grid_style_panels"

with open(os.path.join(panel_dir, "untangle_panel_segments.tsv"), newline="") as handle:
    segments = list(csv.DictReader(handle, delimiter="\t"))
with open(os.path.join(panel_dir, "untangle_panel_summary.tsv"), newline="") as handle:
    summary = list(csv.DictReader(handle, delimiter="\t"))
with open(os.path.join(panel_dir, "untangle_panel_manifest.tsv"), newline="") as handle:
    manifest = list(csv.DictReader(handle, delimiter="\t"))

events = [
    "PAR1_XY_positive_control",
    "PAN027_chr9q_chr3q_PHR_candidate",
    "PAN028_chr9q_chr3q_PHR_candidate",
]

if [row["event_id"] for row in summary] != events:
    raise SystemExit("summary rows are not in required event order")
if [row["event_id"] for row in manifest] != events:
    raise SystemExit("manifest rows are not in required event order")
if len(segments) != 38:
    raise SystemExit(f"expected 38 strict selected segment rows, observed {len(segments)}")

required_roles = {
    "PAR1_XY_positive_control": {"same-chromosome context", "PAR positive control"},
    "PAN027_chr9q_chr3q_PHR_candidate": {
        "same-chromosome context",
        "primary donor",
        "side fragment",
        "low-confidence tail",
    },
    "PAN028_chr9q_chr3q_PHR_candidate": {
        "same-chromosome context",
        "primary donor",
        "side fragment",
    },
}

for event in events:
    roles = {row["event_role"] for row in segments if row["event_id"] == event}
    missing = required_roles[event] - roles
    if missing:
        raise SystemExit(f"{event} missing roles: {sorted(missing)}")

for row in segments:
    if row["nb"] != "1":
        raise SystemExit(f"non-primary nb row found: {row}")
    q0 = int(row["native_query_start_0based"])
    q1 = int(row["native_query_end_0based_exclusive"])
    if not q0 < q1:
        raise SystemExit(f"bad native query interval: {row}")
    if row["panel_geometry_source"] != "paper_prep/_brainstorming/fig5_synteny_recombination_schematic/selected_segments.tsv":
        raise SystemExit(f"unexpected panel geometry source: {row['panel_geometry_source']}")
    if "nb=1" not in row["geometry_source"]:
        raise SystemExit(f"unexpected upstream geometry source: {row['geometry_source']}")
    if row["side_fragment_caveat"] == "yes" and row["event_role"] not in {
        "side fragment",
        "low-confidence tail",
    }:
        raise SystemExit(f"bad caveat flag: {row}")

for row in manifest:
    if row["raw_fasta_query_grid_source"] != "not_used":
        raise SystemExit("raw FASTA query-grid output must not be used as geometry")
    if "permissive multimap/nth-best rows" not in row["excluded_geometry"]:
        raise SystemExit("manifest does not record excluded permissive geometry")

ET.parse(os.path.join(panel_dir, f"{prefix}.svg"))

pdf = os.path.join(panel_dir, f"{prefix}.pdf")
png = os.path.join(panel_dir, f"{prefix}.png")
if open(pdf, "rb").read(5) != b"%PDF-":
    raise SystemExit("PDF does not have a PDF header")
if open(png, "rb").read(8) != b"\x89PNG\r\n\x1a\n":
    raise SystemExit("PNG does not have a PNG header")

print("validated strict untangle query-grid-style panel outputs")
PY

if command -v file >/dev/null 2>&1; then
  file "$panel_dir/$prefix.pdf" "$panel_dir/$prefix.png" "$panel_dir/$prefix.svg" >/dev/null
fi
