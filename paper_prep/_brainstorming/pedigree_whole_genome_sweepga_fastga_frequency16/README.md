# Fig5 whole-genome sweepGA/FastGA frequency-16 sensitivity

Date: 2026-06-21

This package tests whether an explicit FastGA k-mer occurrence threshold of 16
rescues chr3 target homology for the Fig5 PAN027/PAN028 chr9 candidate windows
in whole-genome sweepGA/FastGA output.

Primary command shape:

```bash
/home/erikg/.cargo/bin/sweepga \
  --fastga \
  --fastga-frequency 16 \
  --num-mappings many:many \
  --scaffold-jump 0 \
  --temp-dir /dev/shm/... \
  --output-file ... \
  QUERY.fa TARGET.fa
```

The three comparisons are the same joint-parent whole-genome inputs used by
`paper_prep/_brainstorming/pedigree_whole_genome_sweepga_updated_bin/`:

- `PAN027pat_vs_PAN011_joint`
- `PAN027mat_vs_PAN010_joint`
- `PAN028mat_vs_PAN027_joint`

## Binary Provenance

`summaries/sweepga_binary.tsv` records the explicit sweepGA path, `which`,
realpath, version, sha256, and `--help` text for `/home/erikg/.cargo/bin/sweepga`.

`summaries/fastga_binary.tsv` records `/home/erikg/.cargo/bin/sweepga --check-fastga`.
Each Slurm command log also records the explicit binary path, `which`, compute-node
realpath, sha256, `--help`, `--check-fastga`, `/dev/shm` scratch path, and exact command.

## Workflow

Raw frequency-16 Slurm jobs:

- `1704349`: `PAN027pat_vs_PAN011_joint` ABORTED_TOO_EARLY 00:19:07
- `1704350`: `PAN027mat_vs_PAN010_joint` ABORTED_TOO_EARLY 00:19:07
- `1704355`: `PAN027pat_vs_PAN011_joint` RAW_PAF_OK 05:13:27
- `1704356`: `PAN027mat_vs_PAN010_joint` RAW_PAF_OK 02:37:41
- `1704351`: `PAN028mat_vs_PAN027_joint` RAW_PAF_OK 05:36:43

All three jobs are full whole-genome runs. Scratch is explicitly under `/dev/shm`;
`$SLURM_TMPDIR` is not used as sweepGA/FastGA scratch.

## Result

Direct answer: **PAN027 yes; PAN028 yes**. Explicit
`--fastga-frequency 16` made
sweepGA/FastGA emit chr3 target rows for both PAN027/PAN028 candidate windows in raw PAF.

Raw support:

- `PAN027_chr9q_chr3q_PHR_candidate`: 39 raw chr3 rows, 536932 bp summed query-window overlap, 261767 bp query-union coverage.
- `PAN028_chr9q_chr3q_PHR_candidate`: 36 raw chr3 rows, 536176 bp summed query-window overlap, 261731 bp query-union coverage.

Because raw chr3 rows appeared at frequency 16, 10 kb `pafchop-rs` and chopped
sweepGA `many:many`/`4:many` filters were run. The chr3 signal persists through
both required chopped evidence layers:

- PAN027 `many:many` chopped: 79 chr3 rows, 536932 bp summed overlap, 261767 bp query-union coverage.
- PAN027 `4:many` chopped: 27 chr3 rows, 244914 bp summed overlap, 184498 bp query-union coverage.
- PAN028 `many:many` chopped: 77 chr3 rows, 536176 bp summed overlap, 261731 bp query-union coverage.
- PAN028 `4:many` chopped: 33 chr3 rows, 294219 bp summed overlap, 202804 bp query-union coverage.

Comparator summary:

- `PAN027_chr9q_chr3q_PHR_candidate`: frequency16 raw chr3 `yes`, prior sweepGA raw `no`, wfmash p95 `yes`
- `PAN028_chr9q_chr3q_PHR_candidate`: frequency16 raw chr3 `yes`, prior sweepGA raw `no`, wfmash p95 `yes`

The updated wfmash `-p 95` comparator is treated as expected-positive evidence,
not as a filter input. The frequency-16 run supports the seed-frequency sparsification hypothesis behind the wfmash-positive / sweepGA-negative discrepancy: the prior updated-bin no-explicit-frequency sweepGA run used the much stricter effective FastGA `-f2` setting and missed raw chr3 support, while the explicit `-f16` run finished all three whole-genome comparisons and recovered chr3 rows for both candidate windows. Unlike `-f100`, `-f16` was slow but not pathological at this scale.

Validated with:

```bash
paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/scripts/validate_outputs.sh
```

## Required Summaries

- `summaries/sweepga_binary.tsv`
- `summaries/fastga_binary.tsv`
- `summaries/slurm_jobs.tsv`
- `summaries/raw_chr3_support.tsv`
- `summaries/frequency_sensitivity_summary.tsv`
- `summaries/pathological_runtime.tsv` if pathological
- `summaries/chop_manifest.tsv`, `summaries/filter_manifest.tsv`, and
  `summaries/candidate_window_support.tsv` if chopping/filtering ran

Raw, chopped, and filtered PAFs and checksums are ignored.
