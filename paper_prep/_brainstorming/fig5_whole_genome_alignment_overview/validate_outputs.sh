#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

required=(
  fig5_whole_genome_alignment_overview.pdf
  fig5_whole_genome_alignment_overview.png
  fig5_whole_genome_alignment_overview.svg
  whole_genome_binned_support.tsv
  whole_genome_support_matrix.tsv
  whole_genome_method_manifest.tsv
  README.md
  run_overview.sh
  scripts/build_whole_genome_alignment_overview.py
  scripts/plot_whole_genome_alignment_overview.R
)

for path in "${required[@]}"; do
  if [[ ! -s "$path" ]]; then
    echo "missing or empty: $path" >&2
    exit 1
  fi
done

python3 - <<'PY'
import csv
import os
import sys

def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

manifest = read("whole_genome_method_manifest.tsv")
binned = read("whole_genome_binned_support.tsv")
matrix = read("whole_genome_support_matrix.tsv")

ok_methods = [r for r in manifest if r["status"] == "OK"]
if len(ok_methods) < 28:
    sys.exit(f"expected at least 28 complete method/comparison/chop rows, found {len(ok_methods)}")

families = {r["method_family"] for r in ok_methods}
required_families = {"untangle", "sweepga_fastga_f16", "sweepga_fastga_f32", "wfmash_p95_updated_bin"}
missing = required_families - families
if missing:
    sys.exit(f"missing method families: {sorted(missing)}")

if len(binned) < 50_000:
    sys.exit(f"binned support row count too small: {len(binned)}")
if len(matrix) < 500:
    sys.exit(f"support matrix row count too small: {len(matrix)}")

chroms = {r["query_chrom"] for r in binned}
expected_chroms = {f"chr{i}" for i in range(1, 23)} | {"chrX"}
if not expected_chroms.issubset(chroms):
    sys.exit(f"not all query chromosomes are represented: missing {sorted(expected_chroms - chroms)}")

if not any(r["no_support"] == "yes" for r in binned):
    sys.exit("no explicit no-support bins found")
if not any(r["callout_event_id"] == "PAR1_XY_positive_control" for r in binned):
    sys.exit("PAR1 callout bins not found")
if not any(r["callout_event_id"] == "PAR1_XY_positive_control" and r["top_interchrom_target_chrom"] == "chrY" for r in binned):
    sys.exit("PAR1 callout does not expose chrY as top interchromosomal retained target")
if not any("chr9q_chr3q" in r["callout_event_id"] for r in binned):
    sys.exit("chr9q/chr3q callout bins not found")
if not any(r["interchromosomal"] == "yes" for r in matrix):
    sys.exit("support matrix has no interchromosomal rows")

for artifact in [
    "fig5_whole_genome_alignment_overview.pdf",
    "fig5_whole_genome_alignment_overview.png",
    "fig5_whole_genome_alignment_overview.svg",
]:
    if os.path.getsize(artifact) < 10_000:
        sys.exit(f"rendered artifact is suspiciously small: {artifact}")

print(f"OK: {len(ok_methods)} method rows, {len(binned)} bins, {len(matrix)} matrix rows")
PY
