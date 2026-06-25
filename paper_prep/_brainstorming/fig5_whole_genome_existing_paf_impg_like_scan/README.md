# Fig5 whole-genome existing-PAF IMPG-like scan

This directory contains a PAF-only, query-space genome-wide scan for Fig5 follow-up.
It reuses existing whole-genome WFMASH p95 and SweepGA/FastGA f32 PAFs and does not
build seqwish or ODGI graphs.

Primary entry point:

```bash
sbatch scripts/submit_whole_genome_scan.sh
```

The Slurm wrapper requests 48 CPUs on the `workers,octopus` partition by default.
The summarizer derives worker counts and decompression threads from
`SLURM_CPUS_PER_TASK`; do not run the full scan directly on the login/head node.

Important outputs:

- `manifests/paf_inputs.tsv` - existing raw and filtered PAF inputs used by the scan.
- `summaries/bin_target_support_manifest.tsv` - manifest of the full per query-bin
  and target-chromosome/arm worker shard TSVs in `summaries/tmp_worker_bin_support/`.
- `summaries/target_support_totals.tsv` - compact genome-wide target-support totals.
- `summaries/focal_region_summary.tsv` - PAN027/PAN028 chr9q->chr3q, PAR, and acrocentric controls.
- `summaries/resource_usage.tsv` - CPU allocation and helper-thread accounting.
- `REPORT.md` - concise interpretation and reproducibility notes.

If `pyarrow` is available in the Slurm environment, the TSV summaries are also written
as Parquet sidecars. This worktree currently records TSV outputs as the portable
source of truth because the default Python environment does not provide `pyarrow`.

The full monolithic `bin_target_support.tsv.gz` / `bin_best_support.tsv.gz`
materialization path is implemented in `scripts/summarize_existing_paf_impg_like_scan.py`.
On this cluster the multi-GB combine step was cancelled by Slurm after the 12 worker
shards had completed, so the shard TSVs are retained as the full bin-level output and
the committed summaries are compact derivatives from those shards.
