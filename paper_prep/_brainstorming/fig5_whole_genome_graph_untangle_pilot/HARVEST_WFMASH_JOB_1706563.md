# WFMASH graph untangle harvest: Slurm 1706563

Date: 2026-06-25

## Scope

This note harvests the WFMASH graph/untangle pilot requested for Fig5:

- run ID: `wfmash_p95_updated_bin_l2000`
- primary Slurm job: `1706563`
- configured follow-on row: `sweepga_fastga_f32_l2000`
- source pilot directory:
  `paper_prep/_brainstorming/fig5_whole_genome_graph_untangle_pilot/`

The intended success artifact pair,
`summaries/wfmash_p95_updated_bin_l2000.focus_summary.tsv` and
`summaries/wfmash_p95_updated_bin_l2000.focus_segments.tsv`, was not produced
by job `1706563` or the subsequent WFMASH retries available for inspection.

## Slurm State And Logs

`sacct` records job `1706563` as canceled:

| job | state | exit | elapsed | start | end | node |
| --- | --- | --- | --- | --- | --- | --- |
| `1706563` | `CANCELLED+` | `0:0` | `00:02:13` | `2026-06-25T16:24:30` | `2026-06-25T16:26:43` | `octopus07` |
| `1706563.batch` | `CANCELLED` | `0:15` | `00:02:15` | `2026-06-25T16:24:30` | `2026-06-25T16:26:45` | `octopus07` |

The stdout log for `1706563` is empty. The stderr log contains the decisive
failure:

```text
[seqwish::seqidx] 0.000 indexing sequences
unknown file format given to seqindex_t
[seqwish] WARNING: input FASTA file contains empty sequences, which will be ignored.
slurmstepd: error: *** JOB 1706563 ON octopus07 CANCELLED AT 2026-06-25T16:26:46 ***
```

The submission manifest already records why this happened: the initial script
passed a sequence filename list to `seqwish --seqs` instead of one concatenated
FASTA. That means `1706563` never built a graph, never ran `odgi build`, and
never reached `odgi untangle`.

Two follow-up WFMASH submissions are relevant for interpreting whether the
graph path became valid:

| job | state | relevant observation |
| --- | --- | --- |
| `1706564` | `FAILED` | Failed during staging when `/dev/shm` ran out of space while writing the concatenated whole-genome FASTA. |
| `1706565` | `CANCELLED+` | Reached `seqwish` on MooseFS-backed scratch, but was canceled during transitive closure at about 9.04 percent progress, before graph output, ODGI path listing, or untangle output. |

For `1706565`, stdout is also empty and stderr is only `seqwish` progress plus
the Slurm cancellation line. The last progress lines are still in
`seqwish::transclosure`, not in ODGI:

```text
[seqwish::transclosure] 985.130 9.04% 408000000-409000000 rank_build
slurmstepd: error: *** JOB 1706565 ON octopus07 CANCELLED AT 2026-06-25T16:46:31 ***
```

## Runtime Manifest Harvest

The available runtime command manifest is from corrected retry `1706565`, not
the failed initial `1706563` invocation. It records the effective run
configuration:

| field | value |
| --- | --- |
| `run_id` | `wfmash_p95_updated_bin_l2000` |
| `started_utc` | `2026-06-25T16:29:09Z` |
| `hostname` | `octopus07` |
| `slurm_job_id` | `1706565` |
| `threads` | `48` |
| `seqwish_bin` | `/home/erikg/.guix-profile/bin/seqwish` |
| `odgi_bin` | `/home/erikg/.guix-profile/bin/odgi` |
| `pigz_bin` | `/usr/bin/pigz` |
| `scratch_root` | `paper_prep/_brainstorming/fig5_whole_genome_graph_untangle_pilot/work/tmp` on MooseFS |
| `scratch_df` | 189T size, 94T used, 96T available at startup |
| `query_fasta` | `.../pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027pat_vs_PAN011_joint.query.fa` |
| `target_fasta` | `.../pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027pat_vs_PAN011_joint.target.fa` |
| `filtered_paf` | `.../query_grid_filter/filtered_paf_qgrid_l2000_o0/PAN027pat_vs_PAN011_joint.one_to_one_ani_o0.chopped_l2000_o0_query_grid.paf.gz` |

Because the job did not complete, the manifest has no appended `finished_utc`,
`gfa`, `og`, `paths`, `bedpe`, `paf`, `focus_summary`, or `focus_segments`
fields.

## Focus Summary Availability

No successful WFMASH graph focus outputs were available to promote:

- missing:
  `summaries/wfmash_p95_updated_bin_l2000.focus_summary.tsv`
- missing:
  `summaries/wfmash_p95_updated_bin_l2000.focus_segments.tsv`

This is expected from the job state. The summarizer runs only after:

1. `seqwish` emits the GFA,
2. `odgi build` emits the `.og`,
3. `odgi paths` writes path lists,
4. two `odgi untangle` calls write BEDPE and PAF outputs.

The inspected WFMASH graph attempts stopped before step 1 completed for the
only corrected run that got past staging.

## Focus Pattern Comparison

The graph/untangle focus bins were configured as:

| focus | role | query | target | graph result |
| --- | --- | --- | --- | --- |
| `chr9q_chr3q_candidate` | Fig5 candidate | PAN027 paternal chr9 | PAN011 chr3 | Not evaluated; no graph/untangle output. |
| `chr9q_native_context` | native control | PAN027 paternal chr9 | PAN011 chr9 | Not evaluated; no graph/untangle output. |
| `par_xy_positive_control` | positive control | PAN027 paternal chrX | PAN011 chrY | Not evaluated; no graph/untangle output. |
| `acrocentric_controls` | acrocentric controls | PAN027 paternal chr13/14/15/21/22 | PAN011 chr13/14/15/21/22 | Not evaluated; no graph/untangle output. |

The successful 2 kb direct-similarity pilot provides the best available
non-graph comparison for the same biological checks:

- chr9q/chr3q candidate: PAN027 has only `2,000 bp` of WFMASH 2 kb chr3 support
  in `1` bin, while same-chromosome chr9 support is `435,398 bp` across `223`
  bins. SweepGA direct filters recover `0 bp` chr3 for PAN027 and
  `499,208-619,724 bp` same-chromosome chr9 support.
- PAR X/Y positive control: PAN027 SweepGA direct filters recover chrY support
  of `144,103 bp` across `73` bins, alongside stronger chrX same-chromosome
  support of `358,190-362,564 bp`.
- Acrocentric controls: PAN027 chr13p, chr14p, chr15p, chr21p, and chr22p
  direct filters each recover approximately the full expected acrocentric
  chromosome target, about `499,999-500,000 bp` across `250` bins, without a
  cross-acrocentric-dominated pattern for this PAN027 comparison.

This comparison does not validate the WFMASH graph itself. It shows that the
direct PAF route has enough control behavior to use as a sanity check, while
the graph route has no completed focus evidence.

## Decision

Do not run `sweepga_fastga_f32_l2000` through the graph/untangle Slurm row yet.

Rationale:

1. The primary job `1706563` failed the graph validity gate before graph
   construction because `seqwish` was given the wrong sequence input format.
2. The corrected WFMASH retry that used the intended concatenated FASTA and
   MooseFS scratch (`1706565`) was canceled during `seqwish` transitive closure
   at about 9 percent progress.
3. No WFMASH GFA, ODGI graph, path list, BEDPE/PAF untangle output, focus
   summary, or focus segments exist for this pilot.
4. Without WFMASH focus summaries, there is no basis to compare chr9q/chr3q,
   PAR X/Y, and acrocentric-control behavior under the graph extraction logic.
5. Running the SweepGA/FastGA follow-on row now would test a different input
   before the WFMASH gate has shown that the whole-genome graph/untangle path is
   operational and interpretable.

Recommended next action is not to submit `sweepga_fastga_f32_l2000` as-is.
Either abandon this graph pilot in favor of the direct 2 kb route already
documented in `fig5_2kb_direct_similarity_pilot/REPORT.md`, or first create a
smaller bounded graph/untangle validation that can complete through ODGI and
produce focus summaries for the four configured focus classes.
