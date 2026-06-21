# Fig5 whole-genome minimap2 asm5 all-chain chr3 homology

This package attempted the Fig5 chr9/chr3 whole-genome local-homology
sensitivity check with source-built minimap2:

```bash
/home/erikg/bin/minimap2 -x asm5 -c --eqx -P --q-occ-frac=0 -t 32 TARGET.fa QUERY.fa | pigz -p 32 > OUT.paf.gz
```

The run used the same full child-haplotype query and full joint-parent target
FASTA construction as the updated wfmash and sweepGA packages. It did not use
chromosome-only or window-only FASTAs.

## Outcome

Minimap2 did not produce an evaluable chr3 yes/no answer in this run.

All three Slurm jobs launched successfully and reached the minimap2 mapping
step, but they were cancelled after about 2 h 35 min because node-local
`paf.gz` files remained gzip-header-sized 4 KB filesystem blocks and no complete
PAF was copied back. The required command logs prove the exact binary, version,
options, and full-genome input paths were used; the run is therefore packaged as
a diagnosed pathological-runtime attempt rather than as chr3-positive or
chr3-negative biological evidence.

Summary answer for the PAN027 and PAN028 candidate windows:

| event | minimap2 chr3 support |
| --- | --- |
| `PAN027_chr9q_chr3q_PHR_candidate` | `not_evaluable` |
| `PAN028_chr9q_chr3q_PHR_candidate` | `not_evaluable` |

This differs from the updated wfmash `-p 95` run, which was chr3-positive for
both Fig5 windows, and from the sweepGA/FastGA default raw PAF evidence, which
was chr3-negative. The minimap2 asm5 all-chain attempt cannot currently arbitrate
between those results because it produced no complete PAF rows.

## Inputs

Input definitions are in `config/comparisons.tsv` and
`config/candidate_windows.tsv`.

Prepared full-genome FASTAs are ignored under `inputs/`. The manifest
`summaries/input_manifest.tsv` records three query FASTAs and three target
FASTAs:

- `PAN027pat_vs_PAN011_joint`: `PAN027#2` full haplotype query against a
  `PAN011#1+PAN011#2` joint-parent target.
- `PAN027mat_vs_PAN010_joint`: `PAN027#1` full haplotype query against a
  `PAN010#1+PAN010#2` joint-parent target.
- `PAN028mat_vs_PAN027_joint`: `PAN028#1` full haplotype query against a
  `PAN027#1+PAN027#2` joint-parent target.

Each manifest row states:

```text
full whole-genome haplotype FASTA records from source assembly .fai; no telomeric-window FASTA
```

## Binary provenance

`summaries/minimap2_binary.tsv` records:

- explicit binary: `/home/erikg/bin/minimap2`
- realpath: `/export/local/home/erikg/bin/minimap2-v2.31-r1302`
- version: `2.31-r1302`
- sha256: `5a0e9d6b351f1aa5d11a5067bd29a33bc50abe70c51fc9be9e1899ec1643c949`
- source checkout: `/home/erikg/minimap2`
- source tag and commit: `v2.31`, `3c28777e7e2dcc90f825de1b9f17a89cca7d4452`

The status is `PASS`.

## Slurm jobs

`summaries/slurm_jobs.tsv` records the three submitted jobs:

| comparison | job | status |
| --- | --- | --- |
| `PAN027pat_vs_PAN011_joint` | `1704346` | `CANCELLED by 1001:0:0:02:34:48` |
| `PAN027mat_vs_PAN010_joint` | `1704347` | `CANCELLED by 1001:0:0:02:34:48` |
| `PAN028mat_vs_PAN027_joint` | `1704348` | `CANCELLED by 1001:0:0:02:34:48` |

The command logs are:

- `logs/v2.31-r1302.PAN027pat_vs_PAN011_joint.asm5_allchains.1704346.command.log`
- `logs/v2.31-r1302.PAN027mat_vs_PAN010_joint.asm5_allchains.1704347.command.log`
- `logs/v2.31-r1302.PAN028mat_vs_PAN027_joint.asm5_allchains.1704348.command.log`

Each command log includes `/home/erikg/bin/minimap2`, version `2.31-r1302`,
`-x asm5`, `-P`, `--q-occ-frac=0`, full query/target FASTA paths, and the
exact `pigz` pipeline command.

## Summaries

Required outputs:

- `summaries/minimap2_binary.tsv`
- `summaries/slurm_jobs.tsv`
- `summaries/paf_file_summary.tsv`
- `summaries/candidate_window_support.tsv`
- `summaries/minimap2_chr3_support_summary.tsv`

Additional diagnostic output:

- `summaries/pathological_runtime.tsv`

Raw PAF output paths are recorded under ignored `raw_paf/v2.31-r1302/`, but no
complete raw PAF or checksum exists because the jobs were cancelled before
copy-back.

## Reproduction

Prepare inputs:

```bash
paper_prep/_brainstorming/pedigree_whole_genome_minimap2_asm5_allchains/scripts/run_prepare_inputs.sh
```

Capture binary provenance:

```bash
paper_prep/_brainstorming/pedigree_whole_genome_minimap2_asm5_allchains/scripts/capture_minimap2_binary.py
```

Submit the same primary Slurm run:

```bash
MINIMAP2_TIME=72:00:00 MINIMAP2_MEM=192G MINIMAP2_CPUS=32 \
  paper_prep/_brainstorming/pedigree_whole_genome_minimap2_asm5_allchains/scripts/submit_minimap2_matrix.sh
```

Refresh summaries after completion or cancellation:

```bash
paper_prep/_brainstorming/pedigree_whole_genome_minimap2_asm5_allchains/scripts/run_summaries.sh
```

The exact per-job command shape is emitted by `scripts/run_minimap2_one.sh` and
recorded in `logs/*.command.log`.

## Notes

No `submission/` files were modified, and no Fig5 schematic was created.
