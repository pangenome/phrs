#!/usr/bin/env python3
"""Write a concise report from the IMPG-like scan summaries."""

from __future__ import annotations

import csv
from pathlib import Path


SCAN_DIR = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fmt_int(value: str) -> str:
    return f"{int(float(value)):,}" if value not in {"", None} else "0"


def top_targets(rows: list[dict[str, str]], n: int = 12) -> list[dict[str, str]]:
    return sorted(rows, key=lambda r: float(r["aligned_bp_sum"]), reverse=True)[:n]


def main() -> None:
    totals = read_tsv(SCAN_DIR / "summaries/target_support_totals.tsv")
    focal = read_tsv(SCAN_DIR / "summaries/focal_region_summary.tsv")
    resources = read_tsv(SCAN_DIR / "summaries/resource_usage.tsv")
    parquet = read_tsv(SCAN_DIR / "summaries/parquet_status.tsv")

    resource = resources[0] if resources else {}
    lines = [
        "# Fig5 whole-genome existing-PAF IMPG-like scan",
        "",
        "## Scope",
        "",
        "This scan reuses existing whole-genome WFMASH p95 and SweepGA/FastGA f32 PAFs.",
        "It does not build a seqwish/ODGI graph. Query alignments are tiled into 2 kb",
        "query-space bins, summarized by target chromosome/arm, and then aggregated into",
        "10 kb and 50 kb plotting bins.",
        "",
        "## Execution",
        "",
        f"- Slurm job id: `{resource.get('slurm_job_id', '')}`",
        f"- Host: `{resource.get('hostname', '')}`",
        f"- `SLURM_CPUS_PER_TASK`: `{resource.get('slurm_cpus_per_task', '')}`",
        f"- Process workers: `{resource.get('process_workers', '')}`",
        f"- `pigz` threads per worker: `{resource.get('pigz_threads_per_worker', '')}`",
        f"- Accounted helper threads: `{resource.get('accounted_helper_threads', '')}`",
        f"- Wall seconds: `{resource.get('wall_seconds', '')}`",
        "",
        "The worker count and decompression threads are derived from `SLURM_CPUS_PER_TASK`.",
        "On the octopus run this accounts for the full 48-CPU allocation across worker",
        "processes and `pigz` helper threads.",
        "",
        "## Key focal comparisons",
        "",
        "| region | method | evidence | comparison | query | target | bins | aligned bp | weighted identity | match distance |",
        "|---|---|---|---|---|---|---:|---:|---:|---:|",
    ]
    for row in sorted(focal, key=lambda r: (r["region"], r["method_id"], r["evidence_layer"], r["comparison_id"])):
        lines.append(
            f"| {row['region']} | {row['method_id']} | {row['evidence_layer']} | {row['comparison_id']} | "
            f"{row['query_chrom']} | {row['target_chrom']} | {fmt_int(row['support_bins'])} | "
            f"{fmt_int(row['aligned_bp_sum'])} | {row['mean_identity_weighted']} | {row['mean_match_distance']} |"
        )
    lines.extend([
        "",
        "## Whole-genome target-support patterns",
        "",
        "Top target-support totals by aligned query bp:",
        "",
        "| method | evidence | comparison | query arm | target arm | class | bins | aligned bp | weighted identity |",
        "|---|---|---|---|---|---|---:|---:|---:|",
    ])
    for row in top_targets(totals):
        lines.append(
            f"| {row['method_id']} | {row['evidence_layer']} | {row['comparison_id']} | "
            f"{row['query_chrom']}{row['query_arm_side']} | {row['target_arm']} | {row['support_class']} | "
            f"{fmt_int(row['support_bins'])} | {fmt_int(row['aligned_bp_sum'])} | {row['mean_identity_weighted']} |"
        )
    lines.extend([
        "",
        "## Outputs",
        "",
        "- `summaries/bin_target_support_manifest.tsv`: manifest of the 12 full per-bin target-support shard TSVs in `summaries/tmp_worker_bin_support/`.",
        "- `summaries/target_support_totals.tsv`: compact whole-genome support totals.",
        "- `summaries/focal_region_summary.tsv`: Fig5 chr9q->chr3q, PAR, and acrocentric controls.",
        "- `summaries/resource_usage.tsv`: Slurm allocation and helper-thread accounting.",
        "- `scripts/summarize_existing_paf_impg_like_scan.py`: full monolithic combine/best-target implementation; the cluster cancelled the multi-GB materialization step, so the worker shards are retained as the full bin-level output.",
        "",
        "## Parquet status",
        "",
        "| TSV | Parquet sidecar | status |",
        "|---|---|---|",
    ])
    for row in parquet:
        lines.append(f"| {row['tsv']} | {row['parquet']} | {row['status']} |")
    lines.extend([
        "",
        "## Arm annotation rule",
        "",
        "No centromere table is required for this PAF-only scan. Query and target arms are",
        "called from CHM13 chromosome sizes using a 500 kb subtelomeric rule: alignments in",
        "the first 500 kb are p-arm, alignments in the final 500 kb are q-arm, and all other",
        "alignments are marked internal. This is appropriate for the Fig5 subtelomeric",
        "candidate/control readout and keeps the full-genome summary auditable.",
        "",
        "## Interpretation",
        "",
        "The focal summary separates the PAN027/PAN028 chr9q-to-chr3q signal from PAR and",
        "acrocentric-p positive controls while retaining the genome-wide target landscape.",
        "The raw many:many layer captures broad direct-similarity support; the filtered",
        "one-to-one layer records the corresponding best-chain support after SweepGA filtering.",
    ])

    (SCAN_DIR / "REPORT.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
