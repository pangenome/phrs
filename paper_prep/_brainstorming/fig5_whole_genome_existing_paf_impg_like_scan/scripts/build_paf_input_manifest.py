#!/usr/bin/env python3
"""Build the input manifest for the whole-genome existing-PAF scan."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
OUTDIR = Path(__file__).resolve().parents[1]

WFMASH_MANIFEST = ROOT / "paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/query_grid_filter_manifest.tsv"
SWEEPGA_MANIFEST = ROOT / "paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency32/summaries/query_grid_chop_filter_manifest.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def add_manifest_rows(rows: list[dict[str, object]], source_manifest: Path, method_id: str, bin_size: int) -> None:
    seen_raw: set[tuple[str, str]] = set()
    for row in read_tsv(source_manifest):
        if row.get("status") != "OK":
            continue
        if int(row["chop_length_bp"]) != bin_size:
            continue
        comparison_id = row["comparison_id"]
        common = {
            "method_id": method_id,
            "comparison_id": comparison_id,
            "bin_size_bp": bin_size,
            "source_manifest": str(source_manifest),
            "filter_id": row.get("filter_id", ""),
            "chunk_mode": row.get("chunk_mode", ""),
            "status": row.get("status", ""),
        }
        raw_key = (comparison_id, row["raw_paf"])
        if raw_key not in seen_raw:
            seen_raw.add(raw_key)
            rows.append({
                **common,
                "evidence_layer": "raw_many_many",
                "paf_path": row["raw_paf"],
                "row_count_manifest": row.get("raw_row_count", ""),
                "sha256": "",
            })
        rows.append({
            **common,
            "evidence_layer": "filtered_one_to_one",
            "paf_path": row["filtered_paf"],
            "row_count_manifest": row.get("filtered_row_count", ""),
            "sha256": row.get("filtered_sha256", ""),
        })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-size", type=int, default=2000)
    parser.add_argument("--out", type=Path, default=OUTDIR / "manifests/paf_inputs.tsv")
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    add_manifest_rows(rows, WFMASH_MANIFEST, "wfmash_p95_updated_bin", args.bin_size)
    add_manifest_rows(rows, SWEEPGA_MANIFEST, "sweepga_fastga_f32", args.bin_size)

    for row in rows:
        path = Path(str(row["paf_path"]))
        row["path_exists"] = "yes" if path.exists() else "no"
        row["path_size_bytes"] = path.stat().st_size if path.exists() else ""

    fields = [
        "method_id",
        "evidence_layer",
        "comparison_id",
        "bin_size_bp",
        "paf_path",
        "path_exists",
        "path_size_bytes",
        "row_count_manifest",
        "sha256",
        "filter_id",
        "chunk_mode",
        "source_manifest",
        "status",
    ]
    write_tsv(args.out, rows, fields)


if __name__ == "__main__":
    main()
