# Fig5 Raw Many:Many IMPG Similarity 2 kb Single-Node Race Report

Generated: 2026-06-27

## Scope

This is a no-array, single-node Slurm race against the current sharded Fig5 raw
many:many IMPG 2 kb scan. It uses the existing raw unfiltered many:many PAFs and
the existing full-genome 2 kb target BEDs recorded by:

- `paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/manifests/shard_manifest.tsv`

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, or new alignment/graph
construction step was run. Each race job invokes `/home/erikg/.cargo/bin/impg
similarity` once with the full target BED and `--threads
"${SLURM_CPUS_PER_TASK}"`.

## Slurm Submission

The available Slurm node classes showed `tux` nodes with 96 CPUs, so the race
was submitted to partition `tux` with `--cpus-per-task=96`.

Exactly six non-array Slurm jobs were submitted:

| job_id | method | comparison | partition | cpus |
| --- | --- | --- | --- | --- |
| 1706862 | sweepga_fastga_frequency32 | PAN027mat_vs_PAN010_joint | tux | 96 |
| 1706863 | sweepga_fastga_frequency32 | PAN027pat_vs_PAN011_joint | tux | 96 |
| 1706864 | sweepga_fastga_frequency32 | PAN028mat_vs_PAN027_joint | tux | 96 |
| 1706865 | wfmash_p95_updated_bin | PAN027mat_vs_PAN010_joint | tux | 96 |
| 1706866 | wfmash_p95_updated_bin | PAN027pat_vs_PAN011_joint | tux | 96 |
| 1706867 | wfmash_p95_updated_bin | PAN028mat_vs_PAN027_joint | tux | 96 |

The dependency finalizer is Slurm job `1706868` with dependency:

```text
afterany:1706862:1706863:1706864:1706865:1706866:1706867
```

The finalizer uses `afterany` so it can write a completion manifest even if one
of the six race jobs fails; it exits nonzero unless all six outputs and metadata
are complete.

Queue state at submission check:

- `1706862` running on `tux05`
- `1706863` pending for resources
- `1706864`-`1706867` pending for priority
- `1706868` pending on dependency

The existing sharded arrays `1706840`-`1706845` were not cancelled.

## Outputs

All-hit IMPG outputs are preserved under:

- `paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_single_node_race/outputs/all_hits/`

Per-job metadata are written under:

- `paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_single_node_race/metadata/jobs/`

Finalizer summaries are written under:

- `paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_single_node_race/summaries/`

Expected summary files:

- `per_window_target_similarity_support.tsv`
- `top_interchromosomal_targets.tsv`
- `all_interchromosomal_targets.tsv`
- `chr9q_chr3q_windows.tsv`
- `par_controls.tsv`
- `acrocentric_controls.tsv`
- `full_genome_target_pattern_tracks.tsv`

The best-per-window summaries use the same rule as the sharded finalizer: one
best interchromosomal hit per 2 kb target window, with ties broken by
`estimated.identity`, `intersection`, `dice`, `cosine`, `jaccard`, then stable
lexical target/other coordinates.

## Manifests

Submission and job manifests:

- `manifests/single_node_job_manifest.tsv`
- `manifests/slurm_submission_manifest.tsv`

Completion manifest after finalizer:

- `manifests/single_node_completion_manifest.tsv`

## Supersession Status

This race is not expected to supersede the current array run until all six
single-node jobs complete, the dependency finalizer succeeds, and the generated
best-per-window summaries validate against the same plotting expectations. If
that succeeds, this simpler single-node execution shape is the intended
candidate to supersede the sharded arrays for this Fig5 raw many:many 2 kb IMPG
scan.
