#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOURCE_TSV="${ROOT}/config/pilot_sources.tsv"

command -v seqwish >/dev/null
command -v odgi >/dev/null
command -v sbatch >/dev/null
command -v pigz >/dev/null

python3 - <<'PY' "${SOURCE_TSV}"
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])
rows = list(csv.DictReader(path.open(), delimiter="\t"))
assert rows, "pilot_sources.tsv has no rows"
required = [
    "run_id",
    "query_fasta",
    "target_fasta",
    "filtered_paf",
    "filtered_paf_sha256",
    "source_manifest",
]
for row in rows:
    for key in required:
        assert row.get(key), f"{row.get('run_id', '<missing>')} missing {key}"
    assert Path(row["query_fasta"]).exists(), row["query_fasta"]
    assert Path(row["target_fasta"]).exists(), row["target_fasta"]
    assert Path(row["filtered_paf"]).exists(), row["filtered_paf"]
    assert Path(row["source_manifest"]).exists(), row["source_manifest"]
print(f"validated {len(rows)} pilot source rows")
PY

python3 "${ROOT}/scripts/summarize_untangle_support.py" --help >/dev/null
bash -n "${ROOT}/scripts/run_graph_untangle_one.sbatch"
