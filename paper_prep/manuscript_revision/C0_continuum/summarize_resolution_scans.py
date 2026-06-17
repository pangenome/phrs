#!/usr/bin/env python3
"""Summarize cached Leiden/UPGMA resolution scans for C0c/D1.

This script is intentionally lightweight: it reads already materialized TSVs
from /moosefs and writes compact revision-support tables into this directory.
It does not recompute communities.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path("/moosefs/guarracino/HPRCv2/PHR_III/similarity")
OUT = Path("paper_prep/manuscript_revision/C0_continuum")

INPUTS = {
    "arm_leiden_scan": ROOT / "hprcv2.1Mb.subtelo.arm-leiden.leiden_scan.tsv",
    "seq_leiden_scan": ROOT / "hprcv2.1Mb.subtelo.seq-leiden.leiden_scan.tsv",
    "seq_upgma_scan": ROOT / "hprcv2.1Mb.subtelo.seq-upgma.k-scan.tsv",
    "arm_leiden_assignments": ROOT / "hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv",
    "arm_upgma_assignments": ROOT / "hprcv2.1Mb.subtelo.arm-upgma-k14.assignments.tsv",
    "seq_leiden_assignments": ROOT / "hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv",
    "seq_leiden_summary": ROOT / "hprcv2.1Mb.subtelo.seq_leiden_community_summary.tsv",
    "seq_upgma_summary": ROOT / "hprcv2.1Mb.subtelo.seq-upgma-k150.summary.tsv",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value: str) -> float:
    if value in {"", "NA", "NaN", "nan"}:
        return math.nan
    return float(value)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize_arm_scan() -> None:
    rows = read_tsv(INPUTS["arm_leiden_scan"])
    by_k: dict[int, list[dict[str, str]]] = {}
    for row in rows:
        by_k.setdefault(int(row["n_communities"]), []).append(row)

    out_rows = []
    for k, group in sorted(by_k.items()):
        resolutions = [as_float(row["resolution"]) for row in group]
        silhouettes = [as_float(row["silhouette"]) for row in group]
        finite = [x for x in silhouettes if not math.isnan(x)]
        best = max(finite) if finite else "NA"
        out_rows.append(
            {
                "n_communities": k,
                "n_resolutions": len(group),
                "min_resolution": min(resolutions),
                "max_resolution": max(resolutions),
                "best_silhouette": best,
            }
        )
    write_tsv(
        OUT / "arm_leiden_resolution_by_k.tsv",
        out_rows,
        ["n_communities", "n_resolutions", "min_resolution", "max_resolution", "best_silhouette"],
    )

    finite_rows = [row for row in rows if not math.isnan(as_float(row["silhouette"]))]
    best = max(finite_rows, key=lambda row: as_float(row["silhouette"]))
    write_tsv(
        OUT / "arm_leiden_best_resolution.tsv",
        [
            {
                "selection_rule": "max_mean_silhouette",
                "resolution": best["resolution"],
                "n_communities": best["n_communities"],
                "modularity": best["modularity"],
                "silhouette": best["silhouette"],
                "source_path": INPUTS["arm_leiden_scan"],
            }
        ],
        ["selection_rule", "resolution", "n_communities", "modularity", "silhouette", "source_path"],
    )


def summarize_sequence_scans() -> None:
    leiden = read_tsv(INPUTS["seq_leiden_scan"])
    best_any = max(leiden, key=lambda row: as_float(row["modularity"]))
    target_rows = [row for row in leiden if int(row["n_communities"]) <= 50]
    best_target = max(target_rows, key=lambda row: as_float(row["modularity"]))
    write_tsv(
        OUT / "sequence_leiden_scan_summary.tsv",
        [
            {
                "selection_rule": "max_modularity_any_community_count",
                "k": best_any["k"],
                "resolution": best_any["resolution"],
                "n_communities": best_any["n_communities"],
                "modularity": best_any["modularity"],
                "source_path": INPUTS["seq_leiden_scan"],
            },
            {
                "selection_rule": "max_modularity_with_n_communities_le_50",
                "k": best_target["k"],
                "resolution": best_target["resolution"],
                "n_communities": best_target["n_communities"],
                "modularity": best_target["modularity"],
                "source_path": INPUTS["seq_leiden_scan"],
            },
        ],
        ["selection_rule", "k", "resolution", "n_communities", "modularity", "source_path"],
    )

    upgma = read_tsv(INPUTS["seq_upgma_scan"])
    best_upgma = max(upgma, key=lambda row: as_float(row["silhouette"]))
    write_tsv(
        OUT / "sequence_upgma_scan_summary.tsv",
        [
            {
                "selection_rule": "max_silhouette",
                "k": best_upgma["k"],
                "silhouette": best_upgma["silhouette"],
                "source_path": INPUTS["seq_upgma_scan"],
            }
        ],
        ["selection_rule", "k", "silhouette", "source_path"],
    )


def write_inventory() -> None:
    rows = []
    for label, path in INPUTS.items():
        rows.append(
            {
                "label": label,
                "path": path,
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else "NA",
                "n_lines": sum(1 for _ in path.open()) if path.exists() else "NA",
            }
        )
    write_tsv(OUT / "scan_file_inventory.tsv", rows, ["label", "path", "exists", "size_bytes", "n_lines"])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    write_inventory()
    summarize_arm_scan()
    summarize_sequence_scans()


if __name__ == "__main__":
    main()
