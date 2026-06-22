# Fig5 minimap2 asm5 all-chain long-run chop/sweep test

This package reruns the Fig5 pedigree whole-genome minimap2 sensitivity test
with `/home/erikg/bin/minimap2` v2.31-r1302 and records the downstream
10 kb chop plus sweepGA-PAF-filter plan.  It corrects the earlier cancelled
attempt by treating early header-only or unflushed node-local PAF output as
non-evaluable rather than chr3-negative.

## Direct answers

- Raw minimap2 PAF: **not evaluable yet**.  Jobs 1704357, 1704358, and 1704359
  were launched on full whole-genome query/target FASTAs and each exceeded the
  required 8 hour minimum long-run window, but no complete raw PAF had been
  copied back to `raw_paf/v2.31-r1302/` when the summaries were written.
  Therefore this run does **not** show absence of chr3 target rows; it only
  shows no complete minimap2 PAF was available after the allowed minimum
  runtime.
- 10 kb chopping plus sweepGA PAF filtering: **not run/evaluable** because the
  required complete raw minimap2 PAFs were not yet present.  The package
  contains scripts and configuration for `pafchop-rs` `l10000_o0` followed by
  sweepGA filtering with `--scaffold-jump 0`, including `many:many` and
  `4:many`, but those stages are gated on complete raw PAFs.
- Comparison context: updated wfmash p95 is chr3-positive for the PAN027/PAN028
  Fig5 candidate windows, while updated sweepGA/FastGA default evidence is
  chr3-negative.  This minimap2 long-run remains non-evaluable at the raw,
  chopped, and filtered layers until raw PAFs complete.

## Inputs and scope

`summaries/input_manifest.tsv` records full whole-genome inputs extracted from
the recovered WashU assemblies under
`/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/`.
Each comparison uses a 23-sequence child haplotype query FASTA and a
46-sequence joint-parent target FASTA.  No chromosome-only or candidate-window
FASTA was used.

Comparisons:

- `PAN027pat_vs_PAN011_joint`
- `PAN027mat_vs_PAN010_joint`
- `PAN028mat_vs_PAN027_joint`

Candidate windows:

- `PAN027#2#chr9:135704825-136204825`, expected target chromosome `chr3`
- `PAN028#1#chr9:134380985-134880985`, expected target chromosome `chr3`

## Minimap2 run

Primary command shape, as recorded in per-job command logs:

```bash
/home/erikg/bin/minimap2 -x asm5 -c --eqx -P --q-occ-frac=0 -t 32 TARGET.fa QUERY.fa | pigz -p 32 > OUT.paf.gz
```

Slurm jobs:

- 1704357: `PAN027pat_vs_PAN011_joint`
- 1704358: `PAN027mat_vs_PAN010_joint`
- 1704359: `PAN028mat_vs_PAN027_joint`

Each job requested 32 CPUs, 192G memory, and 24 hours.  Node-local scratch was
`/dev/shm`; compute nodes reported 189G tmpfs available at job start.  The
runtime policy in every command log states that jobs must not be cancelled
solely because the node-local gzip is header-sized or unflushed before
`08:00:00`.

## Summaries

Required summaries written for the no-complete-PAF state:

- `summaries/minimap2_binary.tsv`
- `summaries/sweepga_binary.tsv`
- `summaries/slurm_jobs.tsv`
- `summaries/paf_file_summary.tsv`
- `summaries/raw_candidate_window_support.tsv`
- `summaries/minimap2_chop_sweep_chr3_support_summary.tsv`
- `summaries/longrun_runtime_diagnosis.tsv`

`summaries/chop_manifest.tsv` and `summaries/filter_manifest.tsv` are
intentionally absent in this snapshot because no complete raw minimap2 PAF was
available to chop or filter.

## Logs

Binary and command provenance:

- `logs/minimap2_v2.31-r1302.help.txt`
- `logs/v2.31-r1302.PAN027pat_vs_PAN011_joint.asm5_allchains.1704357.command.log`
- `logs/v2.31-r1302.PAN027mat_vs_PAN010_joint.asm5_allchains.1704358.command.log`
- `logs/v2.31-r1302.PAN028mat_vs_PAN027_joint.asm5_allchains.1704359.command.log`

The command logs prove `/home/erikg/bin/minimap2`, realpath
`/home/erikg/bin/minimap2-v2.31-r1302`, sha256
`5a0e9d6b351f1aa5d11a5067bd29a33bc50abe70c51fc9be9e1899ec1643c949`,
version `2.31-r1302`, full FASTA inputs, `/dev/shm` scratch, `-x asm5`, `-P`,
and `--q-occ-frac=0`.

## Harvest path if raw PAFs later complete

The Slurm jobs were not cancelled by this package after the 8 hour minimum.  If
they finish under their 24 hour allocations and copy complete raw PAFs to
`raw_paf/v2.31-r1302/`, run:

```bash
paper_prep/_brainstorming/pedigree_whole_genome_minimap2_asm5_allchains_longrun_chop_sweep/scripts/run_pafchop_rs_10kb.sh
paper_prep/_brainstorming/pedigree_whole_genome_minimap2_asm5_allchains_longrun_chop_sweep/scripts/run_filter_matrix.sh
paper_prep/_brainstorming/pedigree_whole_genome_minimap2_asm5_allchains_longrun_chop_sweep/scripts/run_summaries.sh
```

`run_pafchop_rs_10kb.sh` uses `pafchop-rs` with `PAF_CHOP_LENGTH=10000` and
`PAF_CHOP_OVERLAP=0`.  `run_filter_matrix.sh` uses
`/home/erikg/.cargo/bin/sweepga` as a PAF filtering/sweeping stage, not for
realignment, with `--scaffold-jump 0`.
