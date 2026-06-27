# Fig5 Raw Many:Many IMPG Similarity 2 kb Sharded Report

Generated: 2026-06-27

## Scope

This is the corrected replacement execution for Fig5 IMPG similarity at 2 kb
resolution. It uses raw unfiltered many:many PAF-backed IMPG similarity over
full-genome target windows. No new WFMASH, SweepGA, FastGA, minimap2, seqwish,
odgi, alignment, or graph-construction step is run.

## Inputs

The submitted shards consume only:

- `raw_paf` rows from
  `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/query_grid_filter_manifest.tsv`
- `raw_paf` rows from
  `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency32/summaries/query_grid_chop_filter_manifest.tsv`
- query/target FASTAs from
  `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/input_manifest.tsv`

Required comparisons:

- `PAN027mat_vs_PAN010_joint`
- `PAN027pat_vs_PAN011_joint`
- `PAN028mat_vs_PAN027_joint`

Required methods:

- `wfmash_p95_updated_bin`
- `sweepga_fastga_frequency32`

## BED Sharding

The generator wrote exact full-genome target BEDs from each target FASTA `.fai`
using 2,000 bp windows. Windows are not expanded to fixed display widths; the
final window on each contig is allowed to be shorter, for example
`PAN010#joint#h2_chrX 153704000 153705098`.

Shard size is 20,000 target windows. This produced:

- 152 shards for `PAN027mat_vs_PAN010_joint`
- 149 shards for `PAN027pat_vs_PAN011_joint`
- 152 shards for `PAN028mat_vs_PAN027_joint`

Across two methods this is 906 IMPG shard tasks.

## SweepGA BGZF Handling

IMPG 0.4.1 requires BGZF for the SweepGA raw PAFs. The pipeline reused the
previous full-BED attempt's BGZF copies only after validating them with
`bgzip -t` in this task. The original raw source PAF remains recorded
separately from the IMPG alignment PAF in:

- `manifests/raw_paf_bgzf_manifest.tsv`
- `manifests/shard_manifest.tsv`
- per-shard `metadata/*.json`

No filtered, chopped-filtered, one-to-one, or partial timed-out SweepGA TSV is
used as evidence.

## Slurm Submission

Six Slurm arrays were submitted on `workers`, each with 48 CPUs per task and
array concurrency capped at 6 tasks:

- `1706840` sweepga_fastga_frequency32 / PAN027mat_vs_PAN010_joint
- `1706841` sweepga_fastga_frequency32 / PAN027pat_vs_PAN011_joint
- `1706842` sweepga_fastga_frequency32 / PAN028mat_vs_PAN027_joint
- `1706843` wfmash_p95_updated_bin / PAN027mat_vs_PAN010_joint
- `1706844` wfmash_p95_updated_bin / PAN027pat_vs_PAN011_joint
- `1706845` wfmash_p95_updated_bin / PAN028mat_vs_PAN027_joint

The first running SweepGA shard wrote metadata showing:

- Slurm array job ID: `1706840`
- node: `octopus09`
- partition: `workers`
- `SLURM_CPUS_PER_TASK`: `48`
- IMPG path/version: `/home/erikg/.cargo/bin/impg`, `impg 0.4.1`

## Output Finalization

The finalizer script is ready but complete finalized TSVs are not yet available
at this handoff because the Slurm arrays are still running or pending. Once all
shards for a method/comparison complete, run:

```bash
python3 paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/scripts/finalize_2kb_sharded_impg.py
```

That writes one compressed assembled IMPG TSV per method x comparison under
`outputs/assembled/` and then produces:

- `summaries/per_window_target_similarity_support.tsv`
- `summaries/top_interchromosomal_targets.tsv`
- `summaries/all_interchromosomal_targets.tsv`
- `summaries/chr9q_chr3q_windows.tsv`
- `summaries/par_controls.tsv`
- `summaries/acrocentric_controls.tsv`
- `summaries/full_genome_target_pattern_tracks.tsv`

## Finalization Attempt: 2026-06-27T13:47Z

Agent `agent-2842` attempted to finalize this run but did not run the
finalizer because the submitted Slurm arrays were still incomplete. The final
observed scheduler state was:

- `1706840_3` running on `octopus09` for `2:05:57`
- `1706840_14` running on `octopus11` for `1:53:54`
- `1706840_[15-151%6]` pending for resources
- `1706841_[0-148%6]`, `1706842_[0-151%6]`,
  `1706843_[0-151%6]`, `1706844_[0-148%6]`, and
  `1706845_[0-151%6]` pending for priority

At that point the parent Slurm worktree contained 15 shard metadata JSON files,
13 completed temp gzip outputs, and zero final-named shard outputs. This is not
sufficient for `manifests/shard_completion_manifest.tsv` to reach 906 `OK`
rows.

One recovery detail was observed: the Slurm job scripts gzip temp files named
`*.tmp.<array>_<task>.gz`, while the finalizer checks the manifest's final
`*.impg_similarity.tsv.gz` paths. After all jobs exit, normalize each completed
temp gzip file to the `output_tsv_gz` path recorded in its metadata before
running the finalizer.

For plotting summaries, preserve the raw assembled all-hit outputs and
`summaries/all_interchromosomal_targets.tsv` for audit/debug. Reduce
`summaries/per_window_target_similarity_support.tsv` and
`summaries/full_genome_target_pattern_tracks.tsv` to the single best record per
2 kb target window. Use deterministic tie-breaking by sorting within each
`method`, `comparison_id`, `target_seq`, `target_start`, `target_end` group by
descending `estimated_identity`, descending `intersection`, descending
`jaccard_similarity`, then lexicographic `other_seq`, `other_start`,
`other_end`, and `other_arm`.
