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

## Deferred Finalization Check: 2026-06-27T21:49:26Z

Finalization was checked from WG task `finalize-fig5-raw-after-arrays` in
worktree `/moosefs/erikg/phrs/.wg-worktrees/agent-2872`. The live shard tree
remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The Slurm guardrail blocked assembly because at least one of arrays
`1706840`-`1706845` was still active. Exact `squeue` state at
`2026-06-27T21:49:26Z`:

```text
1706841_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN028mat_vs_P
1706840_[50-151%6] PENDING  0:00     (Resources)   fig5_impg_sweepg_PAN027mat_vs_P
1706861            PENDING  0:00     (Dependency)  fig5_impg_finalize_2kb
1706840_49         RUNNING  1:51:00  octopus11     fig5_impg_sweepg_PAN027mat_vs_P
1706840_18         RUNNING  5:16:42  octopus09     fig5_impg_sweepg_PAN027mat_vs_P
```

The previously submitted dependency finalizer job `1706861` is still pending
with `Dependency` reason and should be harvested first when the arrays leave
RUNNING/PENDING. Do not run `scripts/finalize_2kb_sharded_impg.py` manually
against incomplete shards. A delayed WG follow-up was created to re-check Slurm
state: `finalize-fig5-raw-2`. Only after all six arrays are terminal and
successful should that follow-up harvest or run finalization.

When finalization is permitted, preserve the all-hit assembled compressed IMPG
TSVs under `outputs/assembled/` for audit. The plotting summaries must reduce to
one hit per 2 kb query window using this deterministic rule: choose the hit with
the highest similarity/ANI/support score first, then the greatest
aligned/support length, then the stable lexical target coordinate key
(`target_name`, `target_start`, `target_end`) to break any remaining tie.

## Deferred Finalization Check: 2026-06-27T23:53:10Z

Finalization was re-checked from WG task `finalize-fig5-raw-2` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2873`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobIDRaw|JobID|JobName|State|ExitCode|Elapsed|Start|End|NodeList
1706861|1706861|fig5_impg_finalize_2kb|PENDING|0:0|00:00:00|Unknown|Unknown|None assigned
```

Because job `1706861` is still pending for dependency release, it has not yet
owned tmp shard normalization, finalizer execution, or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-27T23:53:10Z`:

```text
1706841_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN028mat_vs_P
1706840_[58-151%6] PENDING  0:00     (Resources)   fig5_impg_sweepg_PAN027mat_vs_P
1706861            PENDING  0:00     (Dependency)  fig5_impg_finalize_2kb
1706840_57         RUNNING  41:28    octopus09     fig5_impg_sweepg_PAN027mat_vs_P
1706840_54         RUNNING  1:00:20  octopus11     fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_54` as `RUNNING` on `octopus11` and
`1706840_57` as `RUNNING` on `octopus09`; arrays `1706841`, `1706842`,
`1706843`, `1706844`, and `1706845` remained pending with no node assigned.
Therefore `scripts/finalize_2kb_sharded_impg.py` was not run, no partial
assembled outputs were harvested, and incomplete shards were not marked as data
failures.

A delayed WG follow-up was created to re-check Slurm state: `finalize-fig5-raw-3`.
That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates.

## Deferred Finalization Check: 2026-06-28T01:57:14Z

Finalization was re-checked from WG task `finalize-fig5-raw-3` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2874`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobIDRaw|JobID|JobName|State|ExitCode|Elapsed|Start|End|NodeList
1706861|1706861|fig5_impg_finalize_2kb|PENDING|0:0|00:00:00|Unknown|Unknown|None assigned
```

Because job `1706861` is still pending for dependency release, it has not yet
owned tmp shard normalization, finalizer execution, or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T01:57:14Z`:

```text
1706841_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN028mat_vs_P
1706840_[63-151%6] PENDING  0:00     (Resources)   fig5_impg_sweepg_PAN027mat_vs_P
1706861            PENDING  0:00     (Dependency)  fig5_impg_finalize_2kb
1706840_62         RUNNING  51:25    octopus09     fig5_impg_sweepg_PAN027mat_vs_P
1706840_60         RUNNING  1:54:07  octopus11     fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_60` as `RUNNING` on `octopus11` and
`1706840_62` as `RUNNING` on `octopus09`; arrays `1706841`, `1706842`,
`1706843`, `1706844`, and `1706845` remained pending with no node assigned.
The checked-in `manifests/shard_completion_manifest.tsv` still has 906 data
rows plus header and includes stale `MISSING_OR_INCOMPLETE` rows for shards
whose Slurm jobs have not yet run or whose completion has not yet been
harvested. Therefore the manifest is not expected to show 906 OK rows at this
guardrail stage.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-4`.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates.

## Deferred Finalization Check: 2026-06-28T04:01:07Z

Finalization was re-checked from WG task `finalize-fig5-raw-4` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2875`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobID|JobName|State|ExitCode|Elapsed|Start|End|NodeList|Reason
1706861|fig5_impg_finalize_2kb|PENDING|0:0|00:00:00|Unknown|Unknown|None assigned|Dependency
```

No live log file matching `*1706861*` was present under the live shard tree, so
job `1706861` has not yet owned tmp shard normalization, finalizer execution,
or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T04:01:07Z`:

```text
1706841_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN028mat_vs_P
1706840_[70-151%6] PENDING  0:00     (Resources)   fig5_impg_sweepg_PAN027mat_vs_P
1706861            PENDING  0:00     (Dependency)  fig5_impg_finalize_2kb
1706840_69         RUNNING  8:26     octopus09     fig5_impg_sweepg_PAN027mat_vs_P
1706840_60         RUNNING  3:58:15  octopus11     fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_60` as `RUNNING` on `octopus11` since
`2026-06-28T00:02:52` and `1706840_69` as `RUNNING` on `octopus09` since
`2026-06-28T03:52:41`; array element `1706840_[70-151%6]` and arrays
`1706841`, `1706842`, `1706843`, `1706844`, and `1706845` remain pending with
no node assigned. The live logs currently extend through
`logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_69.1706840.*`;
shard 69 is still running, and shard 60 continues to emit IMPG progress lines.

The checked-in `manifests/shard_completion_manifest.tsv` still has 906 data
rows plus header. Its state column remains `MISSING_OR_INCOMPLETE` for all 906
rows because the Slurm run has not completed and no finalizer or harvest has
updated the manifest yet. No `outputs/assembled/` or `summaries/` products are
present in this worktree at this blocked stage.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-5`.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task does not
supersede the failed `finalize-fig5-raw`; it preserves the guardrail and defers
finalization until the Slurm dependency chain reaches a terminal successful
state.

## Deferred Finalization Check: 2026-06-28T22:39:32Z

Finalization was re-checked from WG task `finalize-fig5-raw-13` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2884`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobID|JobName|State|ExitCode|Elapsed|Start|End
1706861|fig5_impg_finalize_2kb|PENDING|0:0|00:00:00|Unknown|Unknown
```

No live log file matching `*1706861*`, `*finaliz*log*`, or
`slurm-1706861*` was present under the live shard tree, so job `1706861` has
still not owned tmp shard normalization, finalizer execution, or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T22:39:32Z`:

```text
JOBID|STATE|TIME|NODES|NODELIST(REASON)|NAME
1706841_[0-148%6]|PENDING|0:00|1|(Priority)|fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]|PENDING|0:00|1|(Priority)|fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]|PENDING|0:00|1|(Priority)|fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]|PENDING|0:00|1|(Priority)|fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]|PENDING|0:00|1|(Priority)|fig5_impg_wfmash_PAN028mat_vs_P
1706840_[80-151%6]|PENDING|0:00|1|(Resources)|fig5_impg_sweepg_PAN027mat_vs_P
1706861|PENDING|0:00|1|(Dependency)|fig5_impg_finalize_2kb
1706840_79|RUNNING|11:27:01|1|octopus11|fig5_impg_sweepg_PAN027mat_vs_P
1706840_78|RUNNING|16:42:37|1|octopus09|fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_78` as `RUNNING` on `octopus09` since
`2026-06-28T05:55:27`, `1706840_79` as `RUNNING` on `octopus11` since
`2026-06-28T11:11:03`, `1706840_[80-151%6]` pending, and arrays `1706841`,
`1706842`, `1706843`, `1706844`, and `1706845` pending with no node assigned.
The concrete live log paths for the two running shards are:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.err`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.err`

Both the live and checked-in `manifests/shard_completion_manifest.tsv` files
have 906 data rows plus header. Their state columns still contain 906
`MISSING_OR_INCOMPLETE` rows because the Slurm run has not completed and no
finalizer or harvest has updated the manifest yet. The checked-in
`manifests/assembled_outputs.tsv` still has only its header. Neither the live
tree nor the checked-in mirror currently has `outputs/assembled/` or
`summaries/` products, so the six assembled compressed outputs and downstream
summary tables remain correctly blocked by RUNNING/PENDING Slurm state.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-14`. Creating the follow-up
as a direct child of `finalize-fig5-raw-13` hit the configured graph depth
limit, so it was created at the same usable monitoring level with
`finalize-fig5-raw-11` as its dependency.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task does not
supersede the failed `finalize-fig5-raw`; it preserves the guardrail and defers
finalization until the Slurm dependency chain reaches a terminal successful
state.

## Deferred Finalization Check: 2026-06-28T16:25:58Z

Finalization was re-checked from WG task `finalize-fig5-raw-10` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2881`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobID|JobName|State|ExitCode|Elapsed|Start|End|NodeList
1706861|fig5_impg_finalize_2kb|PENDING|0:0|00:00:00|Unknown|Unknown|None assigned
```

No live log file matching `*1706861*`, `finalize*1706861*`,
`*finalize*.out`, or `*finalize*.err` was present under the live shard tree, so
job `1706861` has still not owned tmp shard normalization, finalizer execution,
or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T16:25:58Z`:

```text
JOBID|STATE|REASON|TIME|TIME_LIMIT|NODES|NODELIST(REASON)|NAME
1706841_[0-148%6]|PENDING|Priority|0:00|1-00:00:00|1|(Priority)|fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]|PENDING|Priority|0:00|1-00:00:00|1|(Priority)|fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]|PENDING|Priority|0:00|1-00:00:00|1|(Priority)|fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]|PENDING|Priority|0:00|1-00:00:00|1|(Priority)|fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]|PENDING|Priority|0:00|1-00:00:00|1|(Priority)|fig5_impg_wfmash_PAN028mat_vs_P
1706840_[80-151%6]|PENDING|Resources|0:00|1-00:00:00|1|(Resources)|fig5_impg_sweepg_PAN027mat_vs_P
1706861|PENDING|Dependency|0:00|4:00:00|1|(Dependency)|fig5_impg_finalize_2kb
1706840_79|RUNNING|None|5:14:42|1-00:00:00|1|octopus11|fig5_impg_sweepg_PAN027mat_vs_P
1706840_78|RUNNING|None|10:30:18|1-00:00:00|1|octopus09|fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_78` as `RUNNING` on `octopus09` since
`2026-06-28T05:55:27`, `1706840_79` as `RUNNING` on `octopus11` since
`2026-06-28T11:11:03`, `1706840_[80-151%6]` pending, and arrays `1706841`,
`1706842`, `1706843`, `1706844`, and `1706845` pending with no node assigned.
The concrete live log paths for the two running shards are:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.err`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.err`

Both the live and checked-in `manifests/shard_completion_manifest.tsv` files
have 906 data rows plus header. The checked-in manifest state column still
contains 906 `MISSING_OR_INCOMPLETE` rows because the Slurm run has not
completed and no finalizer or harvest has updated the manifest yet. The
checked-in `manifests/assembled_outputs.tsv` still has only its header. The
checked-in tree currently has no `outputs/assembled/` or `summaries/` products,
so the six assembled compressed outputs and downstream summary tables
(`per_window_target_similarity_support.tsv`,
`full_genome_target_pattern_tracks.tsv`, `top_interchromosomal_targets.tsv`,
`all_interchromosomal_targets.tsv`, `chr9q_chr3q_windows.tsv`,
`par_controls.tsv`, and `acrocentric_controls.tsv`) remain correctly blocked by
RUNNING/PENDING Slurm state.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-11`.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task does not
supersede the failed `finalize-fig5-raw`; it preserves the guardrail and defers
finalization until the Slurm dependency chain reaches a terminal successful
state.

## Deferred Finalization Check: 2026-06-28T06:04:49Z

Finalization was re-checked from WG task `finalize-fig5-raw-5` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2876`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobID|JobName|State|ExitCode|Elapsed|Start|End|NodeList|Reason
1706861|fig5_impg_finalize_2kb|PENDING|0:0|00:00:00|Unknown|Unknown|None assigned|Dependency
```

No live log file matching `*1706861*` was present under the live shard tree, so
job `1706861` has still not owned tmp shard normalization, finalizer execution,
or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T06:04:49Z`:

```text
1706841_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN028mat_vs_P
1706840_[79-151%6] PENDING  0:00     (Resources)   fig5_impg_sweepg_PAN027mat_vs_P
1706861            PENDING  0:00     (Dependency)  fig5_impg_finalize_2kb
1706840_78         RUNNING  8:34     octopus09     fig5_impg_sweepg_PAN027mat_vs_P
1706840_60         RUNNING  6:01:09  octopus11     fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_60` as `RUNNING` on `octopus11` since
`2026-06-28T00:02:52`, `1706840_78` as `RUNNING` on `octopus09` since
`2026-06-28T05:55:27`, `1706840_[79-151%6]` pending, and arrays `1706841`,
`1706842`, `1706843`, `1706844`, and `1706845` pending with no node assigned.
The live logs currently extend through
`logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.*`.

The checked-in `manifests/shard_completion_manifest.tsv` still has 906 data
rows plus header. Its state column remains `MISSING_OR_INCOMPLETE` for all 906
rows because the Slurm run has not completed and no finalizer or harvest has
updated the manifest yet. The live output tree has shard-level outputs only
under `outputs/shards/sweepga_fastga_frequency32`; no live `outputs/assembled/`
or `summaries/` directories were present at this check.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-6`.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task does not
supersede the failed `finalize-fig5-raw`; it preserves the guardrail and defers
finalization until the Slurm dependency chain reaches a terminal successful
state.

## Deferred Finalization Check: 2026-06-28T08:09:15Z

Finalization was re-checked from WG task `finalize-fig5-raw-6` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2877`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobID|JobName|State|ExitCode|Elapsed|Start|End|NodeList|Reason
1706861|fig5_impg_finalize_2kb|PENDING|0:0|00:00:00|Unknown|Unknown|None assigned|Dependency
```

No live log file matching `*1706861*` was present under the live shard tree, so
job `1706861` has still not owned tmp shard normalization, finalizer execution,
or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T08:09:15Z`:

```text
1706841_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN028mat_vs_P
1706840_[79-151%6] PENDING  0:00     (Resources)   fig5_impg_sweepg_PAN027mat_vs_P
1706861            PENDING  0:00     (Dependency)  fig5_impg_finalize_2kb
1706840_78         RUNNING  2:13:48  octopus09     fig5_impg_sweepg_PAN027mat_vs_P
1706840_60         RUNNING  8:06:23  octopus11     fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_60` as `RUNNING` on `octopus11` since
`2026-06-28T00:02:52`, `1706840_78` as `RUNNING` on `octopus09` since
`2026-06-28T05:55:27`, `1706840_[79-151%6]` pending, and arrays `1706841`,
`1706842`, `1706843`, `1706844`, and `1706845` pending with no node assigned.
The live logs currently extend through
`logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.*`.

Both the live and checked-in `manifests/shard_completion_manifest.tsv` files
have 906 data rows plus header. Their state columns still contain 906
`MISSING_OR_INCOMPLETE` rows because the Slurm run has not completed and no
finalizer or harvest has updated the manifest yet. The live output tree still
has no `outputs/assembled/` or `summaries/` products.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-7`.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task does not
supersede the failed `finalize-fig5-raw`; it preserves the guardrail and defers
finalization until the Slurm dependency chain reaches a terminal successful
state.

## Deferred Finalization Check: 2026-06-28T10:13:29Z

Finalization was re-checked from WG task `finalize-fig5-raw-7` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2878`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobID|JobName|State|ExitCode|Start|End|Elapsed|NodeList|Reason
1706861|fig5_impg_finalize_2kb|PENDING|0:0|Unknown|Unknown|00:00:00|None assigned|Dependency
```

No live log file matching `*1706861*` was present under the live shard tree, so
job `1706861` has still not owned tmp shard normalization, finalizer execution,
or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T10:13:29Z`:

```text
1706841_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN028mat_vs_P
1706840_[79-151%6] PENDING  0:00     (Resources)   fig5_impg_sweepg_PAN027mat_vs_P
1706861            PENDING  0:00     (Dependency)  fig5_impg_finalize_2kb
1706840_78         RUNNING  4:17:51  octopus09     fig5_impg_sweepg_PAN027mat_vs_P
1706840_60         RUNNING  10:10:26 octopus11     fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_60` as `RUNNING` on `octopus11` since
`2026-06-28T00:02:52`, `1706840_78` as `RUNNING` on `octopus09` since
`2026-06-28T05:55:27`, `1706840_[79-151%6]` pending, and arrays `1706841`,
`1706842`, `1706843`, `1706844`, and `1706845` pending with no node assigned.
The concrete live log paths for the two running shards are:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_60.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_60.1706840.err`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.err`

Both the live and checked-in `manifests/shard_completion_manifest.tsv` files
have 906 data rows plus header. Their state columns still contain 906
`MISSING_OR_INCOMPLETE` rows because the Slurm run has not completed and no
finalizer or harvest has updated the manifest yet. The live output tree still
has no `outputs/assembled/` or `summaries/` products.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-8`.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task does not
supersede the failed `finalize-fig5-raw`; it preserves the guardrail and defers
finalization until the Slurm dependency chain reaches a terminal successful
state.

## Deferred Finalization Check: 2026-06-28T12:17:23Z

Finalization was re-checked from WG task `finalize-fig5-raw-8` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2879`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobID|JobName|State|ExitCode|Elapsed|Start|End|NodeList
1706861|fig5_impg_finalize_2kb|PENDING|0:0|00:00:00|Unknown|Unknown|None assigned
```

No live log file matching `*1706861*` was present under the live shard tree, so
job `1706861` has still not owned tmp shard normalization, finalizer execution,
or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T12:17:23Z`:

```text
1706841_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]  PENDING  0:00     (Priority)    fig5_impg_wfmash_PAN028mat_vs_P
1706840_[80-151%6] PENDING  0:00     (Resources)   fig5_impg_sweepg_PAN027mat_vs_P
1706861            PENDING  0:00     (Dependency)  fig5_impg_finalize_2kb
1706840_79         RUNNING  1:06:02  octopus11     fig5_impg_sweepg_PAN027mat_vs_P
1706840_78         RUNNING  6:21:38  octopus09     fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_60` as `COMPLETED` with exit code `0:0`
on `octopus11` from `2026-06-28T00:02:52` to `2026-06-28T11:11:03`,
`1706840_78` as `RUNNING` on `octopus09` since `2026-06-28T05:55:27`,
`1706840_79` as `RUNNING` on `octopus11` since `2026-06-28T11:11:03`,
`1706840_[80-151%6]` pending, and arrays `1706841`, `1706842`, `1706843`,
`1706844`, and `1706845` pending with no node assigned. The concrete live log
paths for the two running shards are:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.err`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.err`

Both the live and checked-in `manifests/shard_completion_manifest.tsv` files
have 906 data rows plus header. Their state columns still contain 906
`MISSING_OR_INCOMPLETE` rows because the Slurm run has not completed and no
finalizer or harvest has updated the manifest yet. The checked-in
`manifests/assembled_outputs.tsv` still has only its header. Neither the live
tree nor the checked-in mirror currently has `outputs/assembled/` or
`summaries/` products, so the six assembled compressed outputs and downstream
summary tables remain correctly blocked by RUNNING/PENDING Slurm state.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-9`.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task does not
supersede the failed `finalize-fig5-raw`; it preserves the guardrail and defers
finalization until the Slurm dependency chain reaches a terminal successful
state.

## Deferred Finalization Check: 2026-06-28T14:22:00Z

Finalization was re-checked from WG task `finalize-fig5-raw-9` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2880`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobID|JobName|State|ExitCode|Elapsed|Start|End
1706861|fig5_impg_finalize_2kb|PENDING|0:0|00:00:00|Unknown|Unknown
```

No live log file matching `*1706861*`, `*finaliz*log*`, or `slurm-1706861*`
was present under the live shard tree, so job `1706861` has still not owned tmp
shard normalization, finalizer execution, or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T14:22:00Z`:

```text
JOBID|STATE|TIME|NODES|NODELIST(REASON)|NAME
1706841_[0-148%6]|PENDING|0:00|1|(Priority)|fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]|PENDING|0:00|1|(Priority)|fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]|PENDING|0:00|1|(Priority)|fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]|PENDING|0:00|1|(Priority)|fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]|PENDING|0:00|1|(Priority)|fig5_impg_wfmash_PAN028mat_vs_P
1706840_[80-151%6]|PENDING|0:00|1|(Resources)|fig5_impg_sweepg_PAN027mat_vs_P
1706861|PENDING|0:00|1|(Dependency)|fig5_impg_finalize_2kb
1706840_79|RUNNING|3:10:36|1|octopus11|fig5_impg_sweepg_PAN027mat_vs_P
1706840_78|RUNNING|8:26:12|1|octopus09|fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_78` as `RUNNING` on `octopus09` since
`2026-06-28T05:55:27`, `1706840_79` as `RUNNING` on `octopus11` since
`2026-06-28T11:11:03`, `1706840_[80-151%6]` pending, and arrays `1706841`,
`1706842`, `1706843`, `1706844`, and `1706845` pending with no node assigned.
The concrete live log paths for the two running shards are:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.err`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.err`

Both the live and checked-in `manifests/shard_completion_manifest.tsv` files
have 906 data rows plus header. Their state columns still contain 906
`MISSING_OR_INCOMPLETE` rows because the Slurm run has not completed and no
finalizer or harvest has updated the manifest yet. The checked-in
`manifests/assembled_outputs.tsv` still has only its header. Neither the live
tree nor the checked-in mirror currently has `outputs/assembled/` or
`summaries/` products, so the six assembled compressed outputs and downstream
summary tables remain correctly blocked by RUNNING/PENDING Slurm state.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-10`.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task does not
supersede the failed `finalize-fig5-raw`; it preserves the guardrail and defers
finalization until the Slurm dependency chain reaches a terminal successful
state.

## Latest Guardrail Summary: 2026-06-28T16:25:58Z

The latest check is the `finalize-fig5-raw-10` section above. It supersedes the
earlier `finalize-fig5-raw-9` guardrail check only for current Slurm state, not
for data finalization. Finalization remains blocked: `1706861` is `PENDING`
with reason `Dependency`; `1706840_78` and `1706840_79` are `RUNNING`;
`1706840_[80-151%6]` is `PENDING` with reason `Resources`; and arrays
`1706841_[0-148%6]`, `1706842_[0-151%6]`, `1706843_[0-151%6]`,
`1706844_[0-148%6]`, and `1706845_[0-151%6]` are `PENDING` with reason
`Priority`.

The live source tree is
`/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`;
the checked-in target tree is
`/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`.
The manifest still has 906 data rows plus header and 906
`MISSING_OR_INCOMPLETE` states. No assembled outputs or summary tables were
generated because live Slurm jobs are still running or pending. No WFMASH,
SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard normalization,
finalizer, rsync harvest, or partial assembly command was run. Delayed follow-up
task `finalize-fig5-raw-11` was created. When finalization is eventually
permitted, plotting summaries must reduce to one best hit per 2 kb query window
using the deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task still does
not supersede failed `finalize-fig5-raw`.

## Deferred Finalization Check: 2026-06-28T18:30:45Z

Finalization was re-checked from WG task `finalize-fig5-raw-11` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2882`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobIDRaw|JobID|State|ExitCode|Elapsed|Start|End|NodeList|Reason
1706861|1706861|PENDING|0:0|00:00:00|Unknown|Unknown|None assigned|Dependency
```

No live log file matching `*1706861*`, `*finaliz*`, or `slurm-1706861*` was
present under the live shard tree, so job `1706861` has still not owned tmp
shard normalization, finalizer execution, or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T18:30:45Z`:

```text
JOBID|STATE|TIME|NODELIST(REASON)|NAME
1706841_[0-148%6]|PENDING|0:00|(Priority)|fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]|PENDING|0:00|(Priority)|fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]|PENDING|0:00|(Priority)|fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]|PENDING|0:00|(Priority)|fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]|PENDING|0:00|(Priority)|fig5_impg_wfmash_PAN028mat_vs_P
1706840_[80-151%6]|PENDING|0:00|(Resources)|fig5_impg_sweepg_PAN027mat_vs_P
1706861|PENDING|0:00|(Dependency)|fig5_impg_finalize_2kb
1706840_79|RUNNING|7:19:30|octopus11|fig5_impg_sweepg_PAN027mat_vs_P
1706840_78|RUNNING|12:35:06|octopus09|fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_78` as `RUNNING` on `octopus09` since
`2026-06-28T05:55:27`, `1706840_79` as `RUNNING` on `octopus11` since
`2026-06-28T11:11:03`, `1706840_[80-151%6]` pending, and arrays `1706841`,
`1706842`, `1706843`, `1706844`, and `1706845` pending with no node assigned.
The concrete live log paths for the two running shards are:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.err`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.err`

Both the live and checked-in `manifests/shard_completion_manifest.tsv` files
have 906 data rows plus header. Their state columns still contain 906
`MISSING_OR_INCOMPLETE` rows because the Slurm run has not completed and no
finalizer or harvest has updated the manifest yet. The checked-in
`manifests/assembled_outputs.tsv` still has only its header. Neither the live
tree nor the checked-in mirror currently has `outputs/assembled/` or
`summaries/` products, so the six assembled compressed outputs and downstream
summary tables remain correctly blocked by RUNNING/PENDING Slurm state.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-12`.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task does not
supersede the failed `finalize-fig5-raw`; it preserves the guardrail and defers
finalization until the Slurm dependency chain reaches a terminal successful
state.

## Deferred Finalization Check: 2026-06-28T20:34:26Z

Finalization was re-checked from WG task `finalize-fig5-raw-12` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2883`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobIDRaw|JobID|JobName|State|ExitCode|Elapsed|Start|End|NodeList
1706861|1706861|fig5_impg_finalize_2kb|PENDING|0:0|00:00:00|Unknown|Unknown|None assigned
```

No live log file matching `*1706861*`, `*finaliz*log*`, or
`slurm-1706861*` was present under the live shard tree, so job `1706861` has
still not owned tmp shard normalization, finalizer execution, or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T20:34:26Z`:

```text
JOBID|STATE|TIME|NODELIST(REASON)|NAME
1706841_[0-148%6]|PENDING|0:00|(Priority)|fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]|PENDING|0:00|(Priority)|fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]|PENDING|0:00|(Priority)|fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]|PENDING|0:00|(Priority)|fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]|PENDING|0:00|(Priority)|fig5_impg_wfmash_PAN028mat_vs_P
1706840_[80-151%6]|PENDING|0:00|(Resources)|fig5_impg_sweepg_PAN027mat_vs_P
1706861|PENDING|0:00|(Dependency)|fig5_impg_finalize_2kb
1706840_79|RUNNING|9:23:23|octopus11|fig5_impg_sweepg_PAN027mat_vs_P
1706840_78|RUNNING|14:38:59|octopus09|fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_78` as `RUNNING` on `octopus09` since
`2026-06-28T05:55:27`, `1706840_79` as `RUNNING` on `octopus11` since
`2026-06-28T11:11:03`, `1706840_[80-151%6]` pending, and arrays `1706841`,
`1706842`, `1706843`, `1706844`, and `1706845` pending with no node assigned.
The concrete live log paths for the two running shards are:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.err`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.err`

Both the live and checked-in `manifests/shard_completion_manifest.tsv` files
have 906 data rows plus header. Their state columns still contain 906
`MISSING_OR_INCOMPLETE` rows because the Slurm run has not completed and no
finalizer or harvest has updated the manifest yet. The checked-in
`manifests/assembled_outputs.tsv` still has only its header. Neither the live
tree nor the checked-in mirror currently has `outputs/assembled/` or
`summaries/` products, so the six assembled compressed outputs and downstream
summary tables remain correctly blocked by RUNNING/PENDING Slurm state.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-13`. Because the existing
follow-up chain reached the configured graph depth limit, that task was created
at the same graph level as this one instead of as a direct child.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task does not
supersede the failed `finalize-fig5-raw`; it preserves the guardrail and defers
finalization until the Slurm dependency chain reaches a terminal successful
state.

## Deferred Finalization Check: 2026-06-28T22:41:28Z

Finalization was re-checked from WG task `finalize-fig5-raw-14` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2885`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobIDRaw|JobID|State|ExitCode|Elapsed|Start|End|NodeList|Reason
1706861|1706861|PENDING|0:0|00:00:00|Unknown|Unknown|None assigned|Dependency
```

No live log file matching `*1706861*`, `*finaliz*log*`, or
`slurm-1706861*` was present under the live shard tree, so job `1706861` has
still not owned tmp shard normalization, finalizer execution, or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T22:41:28Z`:

```text
JOBID|STATE|TIME|NODES|NODELIST(REASON)|NAME
1706841_[0-148%6]|PENDING|0:00|1|(Priority)|fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]|PENDING|0:00|1|(Priority)|fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]|PENDING|0:00|1|(Priority)|fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]|PENDING|0:00|1|(Priority)|fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]|PENDING|0:00|1|(Priority)|fig5_impg_wfmash_PAN028mat_vs_P
1706840_[80-151%6]|PENDING|0:00|1|(Resources)|fig5_impg_sweepg_PAN027mat_vs_P
1706861|PENDING|0:00|1|(Dependency)|fig5_impg_finalize_2kb
1706840_79|RUNNING|11:30:25|1|octopus11|fig5_impg_sweepg_PAN027mat_vs_P
1706840_78|RUNNING|16:46:01|1|octopus09|fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_78` as `RUNNING` on `octopus09` since
`2026-06-28T05:55:27`, `1706840_79` as `RUNNING` on `octopus11` since
`2026-06-28T11:11:03`, `1706840_[80-151%6]` pending, and arrays `1706841`,
`1706842`, `1706843`, `1706844`, and `1706845` pending with no node assigned.
The concrete live log paths for the two running shards are:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.err`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.err`

Both the live and checked-in `manifests/shard_completion_manifest.tsv` files
have 906 data rows plus header. Their state columns still contain 906
`MISSING_OR_INCOMPLETE` rows because the Slurm run has not completed and no
finalizer or harvest has updated the manifest yet. The checked-in
`manifests/assembled_outputs.tsv` still has only its header. Neither the live
tree nor the checked-in mirror currently has `outputs/assembled/` or
`summaries/` products, so the six assembled compressed outputs and downstream
summary tables remain correctly blocked by RUNNING/PENDING Slurm state.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-15`. Because the existing
follow-up chain reached the configured graph depth limit, that task was created
at the same graph level as this one instead of as a direct child.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task does not
supersede the failed `finalize-fig5-raw`; it preserves the guardrail and defers
finalization until the Slurm dependency chain reaches a terminal successful
state.

## Deferred Finalization Check: 2026-06-28T22:45:49Z

Finalization was re-checked from WG task `finalize-fig5-raw-15` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2886`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobIDRaw|JobID|JobName|State|ExitCode|Elapsed|Start|End
1706861|1706861|fig5_impg_finalize_2kb|PENDING|0:0|00:00:00|Unknown|Unknown
```

No live log file matching `*1706861*`, `*finaliz*log*`, or
`slurm-1706861*` was present under the live shard tree, so job `1706861` has
still not owned tmp shard normalization, finalizer execution, or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T22:45:49Z`:

```text
JOBID|STATE|TIME|NODES|NODELIST(REASON)|NAME
1706841_[0-148%6]|PENDING|0:00|1|(Priority)|fig5_impg_sweepg_PAN027pat_vs_P
1706842_[0-151%6]|PENDING|0:00|1|(Priority)|fig5_impg_sweepg_PAN028mat_vs_P
1706843_[0-151%6]|PENDING|0:00|1|(Priority)|fig5_impg_wfmash_PAN027mat_vs_P
1706844_[0-148%6]|PENDING|0:00|1|(Priority)|fig5_impg_wfmash_PAN027pat_vs_P
1706845_[0-151%6]|PENDING|0:00|1|(Priority)|fig5_impg_wfmash_PAN028mat_vs_P
1706840_[80-151%6]|PENDING|0:00|1|(Resources)|fig5_impg_sweepg_PAN027mat_vs_P
1706861|PENDING|0:00|1|(Dependency)|fig5_impg_finalize_2kb
1706840_79|RUNNING|11:34:23|1|octopus11|fig5_impg_sweepg_PAN027mat_vs_P
1706840_78|RUNNING|16:49:59|1|octopus09|fig5_impg_sweepg_PAN027mat_vs_P
```

`sacct` additionally reported `1706840_78` as `RUNNING` on `octopus09` since
`2026-06-28T05:55:27`, `1706840_79` as `RUNNING` on `octopus11` since
`2026-06-28T11:11:03`, `1706840_[80-151%6]` pending, and arrays `1706841`,
`1706842`, `1706843`, `1706844`, and `1706845` pending with no start or end
time. The concrete live log paths for the two running shards are:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.err`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.err`

Both the live and checked-in `manifests/shard_completion_manifest.tsv` files
have 906 data rows plus header. Their state columns still contain 906
`MISSING_OR_INCOMPLETE` rows because the Slurm run has not completed and no
finalizer or harvest has updated the manifest yet. The checked-in and live
`manifests/assembled_outputs.tsv` files still have only their headers. Neither
the live tree nor the checked-in mirror currently has `outputs/assembled/` or
`summaries/` products, so the six assembled compressed outputs and downstream
summary tables remain correctly blocked by RUNNING/PENDING Slurm state.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-16`. Because a direct child
task would exceed the configured graph depth limit, that task was created at
the same graph level as this one instead of as a direct child.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task does not
supersede the failed `finalize-fig5-raw`; it preserves the guardrail and defers
finalization until the Slurm dependency chain reaches a terminal successful
state.

## Deferred Finalization Check: 2026-06-28T23:51:15Z

Finalization was re-checked from WG task `finalize-fig5-raw-16` in worktree
`/moosefs/erikg/phrs/.wg-worktrees/agent-2891`. The live shard tree remains:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The main target tree remains:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`

The dependency finalizer job was inspected first. `sacct -j 1706861` reported:

```text
JobID|JobName|State|ExitCode|Elapsed|Start|End|NodeList
1706861|fig5_impg_finalize_2kb|PENDING|0:0|00:00:00|Unknown|Unknown|None assigned
```

No live log file matching `*1706861*`, `*final*`, or `slurm-1706861*` was
present under the live shard tree. The harvested Slurm wrapper in
`jobs/finalize_after_slurm_arrays.slurm.sh` would write
`logs/finalize_after_arrays.%j.out` and `logs/finalize_after_arrays.%j.err` in
the main target tree, but those logs are absent because job `1706861` is still
pending on its dependency and has not owned tmp shard normalization, finalizer
execution, or rsync harvest.

The Slurm guardrail still blocks assembly because arrays `1706840`-`1706845`
remain active. Exact `squeue` state at `2026-06-28T23:51:15Z`:

```text
JOBID              PARTITION NAME                                      STATE        TIME     TIME_LIMIT NODES NODELIST(REASON)
1706841_[0-148%6] workers   fig5_impg_sweepg_PAN027pat_vs_P           PENDING      0:00     1-00:00:00 1     (Priority)
1706842_[0-151%6] workers   fig5_impg_sweepg_PAN028mat_vs_P           PENDING      0:00     1-00:00:00 1     (Priority)
1706843_[0-151%6] workers   fig5_impg_wfmash_PAN027mat_vs_P           PENDING      0:00     1-00:00:00 1     (Priority)
1706844_[0-148%6] workers   fig5_impg_wfmash_PAN027pat_vs_P           PENDING      0:00     1-00:00:00 1     (Priority)
1706845_[0-151%6] workers   fig5_impg_wfmash_PAN028mat_vs_P           PENDING      0:00     1-00:00:00 1     (Priority)
1706840_[80-151%6] workers  fig5_impg_sweepg_PAN027mat_vs_P           PENDING      0:00     1-00:00:00 1     (Resources)
1706861           workers   fig5_impg_finalize_2kb                    PENDING      0:00     4:00:00    1     (Dependency)
1706840_79        workers   fig5_impg_sweepg_PAN027mat_vs_P           RUNNING      12:39:39 1-00:00:00 1     octopus11
1706840_78        workers   fig5_impg_sweepg_PAN027mat_vs_P           RUNNING      17:55:15 1-00:00:00 1     octopus09
```

`sacct` additionally reported `1706840_78` as `RUNNING` on `octopus09` since
`2026-06-28T05:55:27`, `1706840_79` as `RUNNING` on `octopus11` since
`2026-06-28T11:11:03`, `1706840_[80-151%6]` pending, and arrays `1706841`,
`1706842`, `1706843`, `1706844`, and `1706845` pending with no node assigned.
The concrete live log paths for the two running shards are:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_78.1706840.err`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.out`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2837/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/logs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.shard_79.1706840.err`

Both the live and checked-in `manifests/shard_completion_manifest.tsv` files
have 906 data rows plus header. Their `state` columns still contain 906
`MISSING_OR_INCOMPLETE` rows because the Slurm run has not completed and no
finalizer or harvest has updated the manifest yet. The checked-in and live
`manifests/assembled_outputs.tsv` files still have only their headers. Neither
the live tree nor the checked-in mirror currently has `outputs/assembled/` or
`summaries/` products, so the six assembled compressed outputs and downstream
summary tables remain correctly blocked by RUNNING/PENDING Slurm state.

No WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, alignment, shard
normalization, finalizer, rsync harvest, or partial assembly command was run.
Incomplete shards were not marked as data failures. A delayed WG follow-up was
created to re-check Slurm state: `finalize-fig5-raw-17`. Because a direct child
task would exceed the configured graph depth limit, that task was created at
the same graph level as this one instead of as a direct child.

That follow-up should again inspect dependency finalizer job `1706861` before
manual finalization. Only after all six arrays are terminal and successful
should it normalize tmp shard filenames if needed, harvest or run the finalizer,
preserve all-hit assembled outputs under `outputs/assembled/`, and verify that
the plotting summaries reduce to one best hit per 2 kb query window using the
documented deterministic tie-break: highest similarity/ANI/support score, then
aligned/support length, then lexical target coordinates. This task does not
supersede the failed `finalize-fig5-raw`; it preserves the guardrail and defers
finalization until the Slurm dependency chain reaches a terminal successful
state.
