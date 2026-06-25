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

Six primary full-BED IMPG jobs were submitted. Current Slurm snapshot at
handoff:

| method | comparison | job | dependency | state | elapsed | node |
|---|---|---:|---|---|---:|---|
| wfmash_p95_updated_bin | PAN027mat_vs_PAN010_joint | 1706572 | none | RUNNING | 01:12:30 | octopus07 |
| wfmash_p95_updated_bin | PAN027pat_vs_PAN011_joint | 1706573 | none | RUNNING | 01:12:30 | octopus08 |
| wfmash_p95_updated_bin | PAN028mat_vs_PAN027_joint | 1706574 | none | RUNNING | 01:12:30 | octopus09 |
| sweepga_fastga_frequency32 | PAN027mat_vs_PAN010_joint | 1706581 | afterok:1706578 | RUNNING | 01:05:33 | octopus10 |
| sweepga_fastga_frequency32 | PAN027pat_vs_PAN011_joint | 1706582 | afterok:1706579 | RUNNING | 00:59:33 | octopus11 |
| sweepga_fastga_frequency32 | PAN028mat_vs_PAN027_joint | 1706583 | afterok:1706580 | PENDING | 00:00:00 | Resources |

The full-BED WFMASH jobs were actively writing output at handoff. No full-BED
IMPG job had failed or exceeded limits, so no region sharding was proposed or
started.

## Current Output State

At handoff, IMPG result TSVs were still uncompressed because running job scripts
compress only after successful IMPG completion. Summaries were therefore not
generated yet.

Current WFMASH partial TSV sizes at handoff:

- PAN027mat_vs_PAN010_joint: 23 MB
- PAN027pat_vs_PAN011_joint: 39 MB
- PAN028mat_vs_PAN027_joint: 48 MB

SweepGA result TSVs were still empty at handoff while their IMPG jobs were
building/using indexes and beginning region processing.

## Next Steps

1. Wait for Slurm jobs `1706572`, `1706573`, `1706574`, `1706581`, `1706582`,
   and `1706583` to complete or fail.
2. If any full-BED IMPG job fails, stop and inspect its `logs/*.err` file. Do
   not shard until the failure is confirmed to be a full-BED processing limit.
3. If all six jobs complete, run:

```bash
python3 paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/summarize_impg_similarity.py
```

4. Commit the generated summary tables under `summary/` and update this report
   with the final methods/results.
