#!/usr/bin/env bash
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

required=(
  "fig5_untangle_whole_genome_overview.pdf"
  "fig5_untangle_whole_genome_overview.png"
  "fig5_untangle_whole_genome_overview.svg"
  "untangle_whole_genome_segments.tsv"
  "untangle_whole_genome_summary.tsv"
  "README.md"
)

for file in "${required[@]}"; do
  path="$here/$file"
  if [[ ! -s "$path" ]]; then
    echo "missing or empty output: $file" >&2
    exit 1
  fi
done

python3 - "$here" <<'PY'
import csv
import sys
from pathlib import Path

here = Path(sys.argv[1])
segments = list(csv.DictReader((here / "untangle_whole_genome_segments.tsv").open(), delimiter="\t"))
summary = list(csv.DictReader((here / "untangle_whole_genome_summary.tsv").open(), delimiter="\t"))

if len(segments) != 1449:
    raise SystemExit(f"expected 1449 strict rows, found {len(segments)}")
if not summary:
    raise SystemExit("summary table is empty")

events = {row["event_id"] for row in segments if row["event_id"]}
required_events = {
    "PAR1_XY_positive_control",
    "PAN027_chr9q_chr3q_PHR_candidate",
    "PAN028_chr9q_chr3q_PHR_candidate",
}
missing = required_events - events
if missing:
    raise SystemExit(f"missing candidate callouts: {sorted(missing)}")

if any(row["strict_primary_path_only"] != "yes" for row in segments):
    raise SystemExit("non-strict segment marker found")

if not any(row["side_fragment_caveat"] == "yes" for row in segments):
    raise SystemExit("side-fragment caveat markers were not preserved")

if not all("not CHM13" in row["coordinate_system"] for row in segments):
    raise SystemExit("coordinate-system note is missing from segment rows")

inter_bp = sum(int(row["interchromosomal_bp"]) for row in summary)
if inter_bp <= 0:
    raise SystemExit("summary has no interchromosomal support")

print(f"validated {len(segments)} strict rows, {len(summary)} summary rows, {len(events)} event callouts")
PY

grep -q "not CHM13" "$here/README.md"
grep -q "strict primary-path" "$here/README.md"
grep -qi "side fragments" "$here/README.md"

echo "OK: whole-genome untangle overview outputs validated"
