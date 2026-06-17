# C0b Sequence-Level Continuum Report

Date: 2026-06-17

## Scope

This report characterizes the full sequence-level Jaccard distance object for the 15,668 HPR sequences without scanning the large compressed odgi similarity TSV on the head node. The analysis was designed to run under Slurm and used the cached sequence-level distance-matrix RDS.

## Inputs

- Distance matrix RDS: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.dist_matrix.rds` (0.969 GB).
- Sequence Leiden assignments: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv` (1.411 MB).
- The inspected upstream odgi similarity TSV is large and was not decompressed or streamed on the head node; this task uses the cached RDS as the analysis substrate.

## Source Inventory Inspected

- Upstream odgi similarity TSV: `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.similarity.tsv.gz` (12,449,059,644 bytes). This file was size-inspected only; it was not decompressed or streamed on the head node.
- Cached sequence distance RDS: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.dist_matrix.rds` (969,036,139 bytes).
- Sequence Leiden k50 assignments: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv` (1,410,827 bytes).
- Upstream RDS builder inspected: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/build_dist_matrix_rds.R`; it documents the odgi TSV columns and the cached RDS construction path.

## Slurm Run

- Submitted job `1703691` with `sbatch paper_prep/manuscript_revision/C0_continuum/scripts/run_sequence_continuum.sbatch`; it completed the table-writing scan but was cancelled after an optional R PNG device hang.
- Replacement job `1703692` used the same inputs and resource request and completed successfully.
- Resource request: partition `allnodes`; 1 node; 1 task; 4 CPUs; 80G memory; 8 hours.
- Execution node and runtime: `octopus02`, 00:01:35 for job `1703692`.
- Working directory: `/moosefs/erikg/phrs/.wg-worktrees/agent-2465`.
- Stdout/stderr: `paper_prep/manuscript_revision/C0_continuum/logs/c0b_seq_continuum.1703692.out` and `paper_prep/manuscript_revision/C0_continuum/logs/c0b_seq_continuum.1703692.err`.
- Completion status: `COMPLETED`, exit code `0:0` from `sacct`.

## Outputs

- `results/sequence_similarity_distribution.tsv`: binned upper-triangle similarity distribution for all pairs and selected categories.
- `results/high_similarity_summary.tsv`: quantile and high-similarity threshold summaries.
- `results/c6_neighborhood_density.tsv`: threshold-density comparison for the C6/q-arm neighborhood and background.
- `results/sequence_similarity_peaks.tsv`: local maxima in the all-pair similarity histogram.
- The compact distributions are TSV-first so they remain stable in batch environments without relying on an R graphics device.

## C6/q-arm Definition

For this task, the arm-level C6/q-arm neighborhood is defined from the end-to-end report as: chr1_q, chr13_q, chr17_q, chr19_q, chr21_q, chr22_q.
Definable in the assignment table: TRUE.

## Distribution Summary

- Matrix size: 15668 sequences; 122,735,278 upper-triangle non-self pairs.
- All-pair approximate median similarity: 0.0125; q90: 0.3075; q99: 0.9525.
- Pairs at similarity >=0.50: 7,182,560 of 122,735,278 (0.05852).
- Within C6/q-arm pairs at similarity >=0.50: 0.5299 of within-C6 pairs; outside-C6 different-arm background: 0.03733.
- Enrichment of within-C6/q-arm density over outside-C6 different-arm background at >=0.50: 14.19x.

## Interpretation

The data support a broad continuous background of low-to-intermediate sequence similarity, with localized high-similarity peaks rather than a single cleanly separated discrete regime. The all-pair histogram and quantiles should therefore be read as evidence for a continuum that is structured by arm and sequence-community neighborhoods.

The C6/q-arm neighborhood is denser for high-similarity sequence pairs than the outside-C6 different-arm background when evaluated by fixed similarity thresholds. This supports treating the q-arm/C6 pattern as an enriched neighborhood within the continuum, not as an isolated bounded class.

## Slurm Safety

The heavy step is loading and scanning the 15,668 x 15,668 sequence-level distance object. It is intended for the sbatch wrapper in this directory. Head-node work was limited to file-size inspection, small TSV previews, script writing, submission, and result synthesis.
