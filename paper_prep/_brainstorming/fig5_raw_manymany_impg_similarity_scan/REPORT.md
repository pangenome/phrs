# Fig5 Raw Many:Many IMPG Similarity Full-BED Run

Generated: 2026-06-25

## Scope

This directory is the only valid output area for the replacement raw many:many
IMPG similarity scan:

`paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/`

The invalid existing-PAF reducer output area was not used.

## Inputs

Primary evidence is raw PAF only:

- WFMASH updated-bin raw many:many PAFs from
  `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/query_grid_filter_manifest.tsv`
- SweepGA/FastGA f32 raw many:many PAFs from
  `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency32/summaries/query_grid_chop_filter_manifest.tsv`
- query and target FASTA names from
  `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/input_manifest.tsv`

No `filtered_paf`, `filtered_one_to_one`, chopped filtered PAFs, or PAF-overlap
reducer outputs are used as primary evidence.

## IMPG Command Validation

Installed IMPG:

- path: `/home/erikg/.cargo/bin/impg`
- version: `impg 0.4.1`
- captured help: `impg_similarity_help.txt`

Command probes found two required adjustments for this installed IMPG:

- `--merge-distance 0` and `--no-merge` are documented together in help, but
  IMPG 0.4.1 rejects using both at once. The submitted jobs use `--no-merge`,
  preserving the requested no-merging behavior.
- `impg similarity` requires `--sequence-files QUERY.fa TARGET.fa` and
  `--gfa-engine poa` for tabular similarity output from precomputed PAF
  alignments. Those required arguments are included in every job script.

The submitted execution shape is therefore:

```bash
/home/erikg/.cargo/bin/impg similarity \
  --alignment-files IMPG_ALIGNMENT.paf.gz \
  --target-bed COMPARISON.full_genome_10kb.bed \
  --sequence-files QUERY.fa TARGET.fa \
  --gfa-engine poa \
  --no-merge \
  --num-mappings many:many \
  --scaffold-jump 0 \
  --threads ${SLURM_CPUS_PER_TASK}
```

Each IMPG job requests 48 CPUs on `workers` and passes exactly
`${SLURM_CPUS_PER_TASK}` to IMPG.

## BED Generation

Generated one full-genome 10 kb target BED per comparison from the target FASTA
`.fai` files:

- `beds/PAN027mat_vs_PAN010_joint.full_genome_10kb.bed`
- `beds/PAN027pat_vs_PAN011_joint.full_genome_10kb.bed`
- `beds/PAN028mat_vs_PAN027_joint.full_genome_10kb.bed`

The 2 kb pass was not generated or submitted because the 10 kb run had not
completed and validated.

## SweepGA BGZF Normalization

The first submitted SweepGA IMPG jobs (`1706575`-`1706577`) failed before any
BED processing because IMPG reported the source PAFs are regular gzip, not
BGZF. This is a compression/indexing format issue, not evidence filtering and
not a full-BED capacity failure.

The pipeline now records both:

- `source_raw_paf`: original raw PAF from the required manifest
- `impg_alignment_paf`: BGZF-normalized copy under `bgzf_raw_paf/`

BGZF normalization jobs:

| comparison | job | state at handoff | node |
|---|---:|---|---|
| PAN027mat_vs_PAN010_joint | 1706578 | COMPLETED | octopus10 |
| PAN027pat_vs_PAN011_joint | 1706579 | COMPLETED | octopus11 |
| PAN028mat_vs_PAN027_joint | 1706580 | COMPLETED | octopus11 |

## Primary IMPG Jobs

Six primary full-BED IMPG jobs were submitted. Final retry status after the
2026-06-26 monitor:

| method | comparison | job | dependency | final state | elapsed | node | output state |
|---|---|---:|---|---|---:|---|---|
| wfmash_p95_updated_bin | PAN027mat_vs_PAN010_joint | 1706572 | none | COMPLETED | 14:35:29 | octopus07 | `outputs/wfmash_p95_updated_bin.PAN027mat_vs_PAN010_joint.full_genome_10kb.impg_similarity.tsv.gz` |
| wfmash_p95_updated_bin | PAN027pat_vs_PAN011_joint | 1706573 | none | COMPLETED | 10:48:48 | octopus08 | `outputs/wfmash_p95_updated_bin.PAN027pat_vs_PAN011_joint.full_genome_10kb.impg_similarity.tsv.gz` |
| wfmash_p95_updated_bin | PAN028mat_vs_PAN027_joint | 1706574 | none | COMPLETED | 09:45:20 | octopus09 | `outputs/wfmash_p95_updated_bin.PAN028mat_vs_PAN027_joint.full_genome_10kb.impg_similarity.tsv.gz` |
| sweepga_fastga_frequency32 | PAN027mat_vs_PAN010_joint | 1706581 | afterok:1706578 | TIMEOUT | 1-00:00:08 | octopus10 | partial uncompressed `outputs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.full_genome_10kb.impg_similarity.tsv` (6.4 GB) |
| sweepga_fastga_frequency32 | PAN027pat_vs_PAN011_joint | 1706582 | afterok:1706579 | TIMEOUT | 1-00:00:08 | octopus11 | partial uncompressed `outputs/sweepga_fastga_frequency32.PAN027pat_vs_PAN011_joint.full_genome_10kb.impg_similarity.tsv` (4.9 GB) |
| sweepga_fastga_frequency32 | PAN028mat_vs_PAN027_joint | 1706583 | afterok:1706580 | CANCELLED after blocker confirmed | 14:30:30 | octopus09 | no TSV rows; IMPG was still collecting/indexing results |

The WFMASH full-BED jobs validated the full-BED execution unit for WFMASH raw
PAFs on 48 CPU `workers` nodes. The SweepGA/FastGA f32 full-BED jobs did not
complete within the requested 24-hour walltime. Jobs `1706581` and `1706582`
were cancelled by Slurm due to time limit, with the relevant stderr endings:

```text
slurmstepd: error: *** JOB 1706581 ON octopus10 CANCELLED AT 2026-06-26T18:18:17 DUE TO TIME LIMIT ***
slurmstepd: error: *** JOB 1706582 ON octopus11 CANCELLED AT 2026-06-26T18:24:17 DUE TO TIME LIMIT ***
```

After both sibling SweepGA full-BED jobs had timed out, job `1706583` was
cancelled manually to stop the now-blocked full-BED run. Its stderr showed IMPG
was still collecting results immediately before cancellation.

## Output State

Finalized IMPG similarity outputs exist for all three WFMASH comparisons only:

- `outputs/wfmash_p95_updated_bin.PAN027mat_vs_PAN010_joint.full_genome_10kb.impg_similarity.tsv.gz`
- `outputs/wfmash_p95_updated_bin.PAN027pat_vs_PAN011_joint.full_genome_10kb.impg_similarity.tsv.gz`
- `outputs/wfmash_p95_updated_bin.PAN028mat_vs_PAN027_joint.full_genome_10kb.impg_similarity.tsv.gz`

SweepGA/FastGA f32 did not produce valid finalized IMPG outputs:

- `outputs/sweepga_fastga_frequency32.PAN027mat_vs_PAN010_joint.full_genome_10kb.impg_similarity.tsv` is a partial, uncompressed 6.4 GB file from timed-out job `1706581`.
- `outputs/sweepga_fastga_frequency32.PAN027pat_vs_PAN011_joint.full_genome_10kb.impg_similarity.tsv` is a partial, uncompressed 4.9 GB file from timed-out job `1706582`.
- `outputs/sweepga_fastga_frequency32.PAN028mat_vs_PAN027_joint.full_genome_10kb.impg_similarity.tsv` is an empty placeholder from cancelled job `1706583`.

The required six complete outputs therefore do not exist. Summary tables were
not generated, because running `summarize_impg_similarity.py` over this mixed
state would either fail on missing `.tsv.gz` files or summarize incomplete
partial evidence.

## Blocker and Proposed Next Execution

The first confirmed full-BED blocker is SweepGA/FastGA f32 runtime on 48 CPU
`workers` nodes. The full 10 kb target BED can be processed for WFMASH, but the
SweepGA raw many:many evidence layer did not finish in a 24-hour full-BED job
for two comparisons.

Recommended next attempt, in order:

1. Retry the three SweepGA/FastGA f32 comparisons on a larger allocation if
   available, preferably `tux` with 96 CPUs and a longer walltime, preserving
   the same full-BED command shape and `--threads ${SLURM_CPUS_PER_TASK}`.
2. If a larger allocation is unavailable or still times out, shard the target
   BED by chromosome/region for SweepGA only. This is now justified by the
   confirmed full-BED walltime failures in jobs `1706581` and `1706582`.
3. After all three SweepGA comparisons have complete `.tsv.gz` outputs, run:

```bash
python3 paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/summarize_impg_similarity.py
```

4. Generate the requested summary deliverables:
   `summary/per_window_target_similarity_support.tsv`,
   `summary/top_interchromosomal_targets.tsv`,
   `summary/all_interchromosomal_targets.tsv`,
   `summary/chr9q_chr3q_windows.tsv`, `summary/par_controls.tsv`,
   `summary/acrocentric_controls.tsv`, and
   `summary/full_genome_target_pattern_tracks.tsv`.
