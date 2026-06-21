# Fig5 whole-genome sweepGA/FastGA frequency-100 sensitivity

Date: 2026-06-21

This package tests whether making FastGA's k-mer occurrence threshold explicit
rescues chr3 target homology for the Fig5 PAN027/PAN028 chr9 candidate windows
in whole-genome sweepGA/FastGA output.

Primary command shape:

```bash
/home/erikg/.cargo/bin/sweepga \
  --fastga \
  --fastga-frequency 100 \
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

The merged prior package contains the scripts and summaries, but its ignored
full-genome input FASTAs are not present in the repository. This package
therefore uses the allowed live worker-package inputs at
`/moosefs/erikg/phrs/.wg-worktrees/agent-2639/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_updated_bin/inputs/`.
The package-local `inputs/` entries are symlinks to those live files.

## Binary Provenance

`summaries/sweepga_binary.tsv` records the explicit sweepGA path, `which`,
realpath, version, sha256, and `--help` text for:

`/home/erikg/.cargo/bin/sweepga`

`summaries/fastga_binary.tsv` records:

```bash
/home/erikg/.cargo/bin/sweepga --check-fastga
```

Each Slurm command log also records the explicit binary path, `which`,
compute-node realpath, sha256, `--help`, `--check-fastga`, scratch path, and the
exact primary command. The FastGA invocation inside sweepGA is logged with
`-f100`.

## Workflow

Raw frequency-100 Slurm jobs:

- `1704343`: `PAN027pat_vs_PAN011_joint`
- `1704344`: `PAN027mat_vs_PAN010_joint`
- `1704345`: `PAN028mat_vs_PAN027_joint`

All three jobs are full whole-genome runs. Scratch is explicitly under
`/dev/shm`; `$SLURM_TMPDIR` is not used as sweepGA/FastGA scratch.

The prior updated-bin no-explicit-frequency run completed with FastGA `-f2` and
reported no raw chr3 target rows for the PAN027 or PAN028 candidate windows. The
updated wfmash `-p 95` comparator is expected-positive and reports chr3 support
for both candidate windows.

## Result

Direct answer: **No**. Explicit `--fastga-frequency 100` did not make
sweepGA/FastGA emit chr3 target rows for the PAN027 or PAN028 candidate windows
in a usable raw PAF.

The primary frequency-100 jobs were pathological at whole-genome scale. All
three entered FastGA with `-f100`, remained CPU-active for 2:38:18 wall time,
used about 98-99 GB of `/dev/shm` scratch per job, and still had zero-byte
`.1aln` outputs with no raw PAF emitted. They were cancelled rather than
allowed to consume the full 24 hour allocation because the prior no-explicit
frequency run completed the same comparisons in 00:06:48-00:14:38 with FastGA
`-f2`.

Because no raw PAF was produced, there were no raw rows to chop or filter.
Downstream 10 kb chopping and `many:many`/`4:many` chopped filtering were
therefore intentionally skipped. This follows the task ordering: inspect raw PAF
support first, and only run chopped filters if raw chr3 rows appear.

The required raw support summary is:

- `PAN027_chr9q_chr3q_PHR_candidate`: no raw PAF; no chr3 target rows emitted
- `PAN028_chr9q_chr3q_PHR_candidate`: no raw PAF; no chr3 target rows emitted

The updated wfmash `-p 95` comparator remains expected-positive for both
candidate windows, while the prior updated sweepGA/FastGA no-explicit-frequency
run remains raw chr3-negative. The frequency-100 result therefore does **not**
support the specific hypothesis that the prior discrepancy is resolved by
raising FastGA's seed frequency threshold to 100. It is still consistent with
the broader sparsification/sensitivity issue in repetitive subtelomeric
sequence: relaxing the occurrence filter greatly expands the FastGA search
space and becomes unusable here before producing raw PAF evidence.

Validated with:

```bash
paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency100/scripts/validate_outputs.sh
```

## Required Summaries

- `summaries/sweepga_binary.tsv`
- `summaries/fastga_binary.tsv`
- `summaries/slurm_jobs.tsv`
- `summaries/raw_chr3_support.tsv`
- `summaries/frequency_sensitivity_summary.tsv`
- `summaries/pathological_runtime.tsv`

Raw PAFs and checksums are ignored under `raw_paf/`.
