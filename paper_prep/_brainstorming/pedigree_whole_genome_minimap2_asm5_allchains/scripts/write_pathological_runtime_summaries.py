#!/usr/bin/env python3
import os

from common import CANDIDATE_WINDOWS, PACKAGE_DIR, read_tsv, write_tsv


RUNTIME_NOTE = (
    "Primary minimap2 asm5 all-chain Slurm jobs were cancelled after ~2h35m "
    "because node-local paf.gz files remained gzip-header-sized and no complete "
    "PAF was copied back; chr3 support is therefore not evaluable from this run."
)
NEXT_COMMAND = (
    "MINIMAP2_TIME=72:00:00 MINIMAP2_MEM=192G MINIMAP2_CPUS=32 "
    "paper_prep/_brainstorming/pedigree_whole_genome_minimap2_asm5_allchains/"
    "scripts/submit_minimap2_matrix.sh"
)


def main():
    jobs = read_tsv(os.path.join(PACKAGE_DIR, "summaries", "slurm_jobs.tsv"))
    windows = read_tsv(CANDIDATE_WINDOWS)

    paf_rows = []
    for job in jobs:
        paf_rows.append({
            "paf_path": job["output_paf"],
            "run_label": job["run_label"],
            "comparison_id": job["comparison_id"],
            "parameter_set": job["parameter_set"],
            "target_chrom": "NO_COMPLETE_PAF",
            "rows": "0",
            "total_rows": "0",
            "status": job["status"],
            "note": RUNTIME_NOTE,
        })
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "paf_file_summary.tsv"), paf_rows, [
        "paf_path",
        "run_label",
        "comparison_id",
        "parameter_set",
        "target_chrom",
        "rows",
        "total_rows",
        "status",
        "note",
    ])

    support_rows = []
    for window in windows:
        matching = [job for job in jobs if job["comparison_id"] == window["comparison_id"]]
        job = matching[0] if matching else {}
        support_rows.append({
            "event_id": window["event_id"],
            "run_label": job.get("run_label", "v2.31-r1302"),
            "comparison_id": window["comparison_id"],
            "parameter_set": job.get("parameter_set", "asm5_allchains"),
            "query_name": window["query_name"],
            "query_start": window["query_start"],
            "query_end": window["query_end"],
            "target_chrom": "NO_COMPLETE_PAF",
            "expected_target_chrom": window["expected_target_chrom"],
            "paf_rows_overlapping_window": "0",
            "query_overlap_bp_sum": "0",
            "query_covered_bp_union": "0",
            "chr3_support": "not_evaluable",
            "paf_path": job.get("output_paf", ""),
            "status": job.get("status", "NO_JOB"),
            "note": RUNTIME_NOTE,
        })
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "candidate_window_support.tsv"), support_rows, [
        "event_id",
        "run_label",
        "comparison_id",
        "parameter_set",
        "query_name",
        "query_start",
        "query_end",
        "target_chrom",
        "expected_target_chrom",
        "paf_rows_overlapping_window",
        "query_overlap_bp_sum",
        "query_covered_bp_union",
        "chr3_support",
        "paf_path",
        "status",
        "note",
    ])

    chr3_rows = []
    for row in support_rows:
        chr3_rows.append({
            "event_id": row["event_id"],
            "comparison_id": row["comparison_id"],
            "query_name": row["query_name"],
            "query_start": row["query_start"],
            "query_end": row["query_end"],
            "expected_target_chrom": row["expected_target_chrom"],
            "minimap2_chr3_support": "not_evaluable",
            "chr3_rows": "0",
            "chr3_query_overlap_bp_sum": "0",
            "chr3_query_covered_bp_union": "0",
            "all_target_chromosomes_overlapping_window": "NO_COMPLETE_PAF",
            "answer": "not_evaluable",
            "comparison_to_other_aligners": (
                "Updated wfmash -p 95 was chr3-positive for this Fig5 window; "
                "sweepGA/FastGA default raw PAF was chr3-negative; minimap2 asm5 "
                "all-chain produced no evaluable complete PAF before cancellation."
            ),
            "next_command": NEXT_COMMAND,
            "note": RUNTIME_NOTE,
        })
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "minimap2_chr3_support_summary.tsv"), chr3_rows, [
        "event_id",
        "comparison_id",
        "query_name",
        "query_start",
        "query_end",
        "expected_target_chrom",
        "minimap2_chr3_support",
        "chr3_rows",
        "chr3_query_overlap_bp_sum",
        "chr3_query_covered_bp_union",
        "all_target_chromosomes_overlapping_window",
        "answer",
        "comparison_to_other_aligners",
        "next_command",
        "note",
    ])

    runtime_rows = []
    for job in jobs:
        runtime_rows.append({
            "comparison_id": job["comparison_id"],
            "job_id": job["job_id"],
            "status": job["status"],
            "node_local_paf_observation": "gzip-header-sized 4K filesystem block at ~2h35m",
            "output_copied_back": "no",
            "next_command": NEXT_COMMAND,
            "note": RUNTIME_NOTE,
        })
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "pathological_runtime.tsv"), runtime_rows, [
        "comparison_id",
        "job_id",
        "status",
        "node_local_paf_observation",
        "output_copied_back",
        "next_command",
        "note",
    ])


if __name__ == "__main__":
    main()
