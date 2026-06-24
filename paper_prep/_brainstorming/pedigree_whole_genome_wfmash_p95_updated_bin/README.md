# Fig5 whole-genome wfmash -p95 updated local binary rerun

Date: 2026-06-21

Task: `fig5-whole-genome-wfmash-p95-updated-bin`

This package is a fresh provenance-controlled rerun of the Fig5 direct
whole-genome wfmash test using the updated local binary at:

`/home/erikg/bin/wfmash`

It does not reuse the prior task-built `tools/wfmash-v0.24.2-built/wfmash`
binary from `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95/`.

## Direct answer

Yes. The updated `/home/erikg/bin/wfmash` literal whole-genome `-p 95` rerun
recovers chr3-target raw PAF rows for both Fig5 autosomal chr9 candidate
windows:

| event | comparison | chr3 rows | query overlap bp sum | query union coverage bp |
| --- | --- | ---: | ---: | ---: |
| `PAN027_chr9q_chr3q_PHR_candidate` | `PAN027pat_vs_PAN011_joint` | 2 | 73,000 | 73,000 |
| `PAN028_chr9q_chr3q_PHR_candidate` | `PAN028mat_vs_PAN027_joint` | 3 | 152,860 | 150,862 |

These values come from the raw gzipped PAFs under ignored `raw_paf/`, before
any posthoc candidate-window slicing beyond summary selection. The middle
`PAN027mat_vs_PAN010_joint` comparison also completed for provenance
completeness but is not one of the two chr9/chr3 candidate windows.

The optional permissive configuration (`-p 95 -w 1k -l 1k -n 50 -f -M`) was not
rerun, because all three required primary literal `-p 95` whole-genome runs
completed and the primary decision is based on those raw PAFs.

## Binary provenance

See `summaries/wfmash_binary.tsv` and `logs/wfmash_updated_bin.help.txt`.

The recorded binary provenance passed the task-created expectations:

- explicit path: `/home/erikg/bin/wfmash`
- realpath from the login/worktree host:
  `/export/local/home/erikg/bin/wfmash-v0.24.2-12-ge040aa10`
- version: `v0.24.2-12-ge040aa10`
- SHA-256:
  `14a6d5c7ac7be8890e904d11121341df118fad0c11193d0e91a9899e18a53d60`
- `which wfmash`: `/home/erikg/.guix-profile/bin/wfmash`, recorded to show why
  plain PATH lookup is unsafe.

The Slurm command logs record the binary actually used on the compute nodes as
`/home/erikg/bin/wfmash`, with the same version and SHA-256. On the compute
nodes `readlink -f /home/erikg/bin/wfmash` resolved under `/home/erikg/bin/`
rather than `/export/local/home/erikg/bin/`; the checksum and version match the
required updated local binary.

## Whole-genome inputs

The comparison definitions are copied from the corrected whole-genome setup and
use the recovered full WashU assemblies:

`/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/`

See:

- `config/comparisons.tsv`
- `summaries/input_manifest.tsv`

Each query FASTA contains 23 full-chromosome haplotype records. Each target
FASTA contains both parent haplotypes, 46 full-chromosome records total, with
target headers collapsed to a joint-parent prefix such as `PAN011#joint#h1_chr3`.
No window-only, chromosome-only, arm-only, or telomeric-only FASTA was used.

Large generated `inputs/` FASTAs are ignored and not committed.

## Slurm runs

All required primary jobs completed:

| comparison | job id | status | raw PAF lines | gzipped bytes | SHA-256 |
| --- | ---: | --- | ---: | ---: | --- |
| `PAN027pat_vs_PAN011_joint` | 1704325 | `COMPLETED_OUTPUT_PRESENT` | 66,591 | 2,827,411 | `53eff77a94f5510f95aa010d00f6eab985d3d7a8c2f433a73d150e907311f804` |
| `PAN027mat_vs_PAN010_joint` | 1704326 | `COMPLETED_OUTPUT_PRESENT` | 66,634 | 2,689,281 | `a6742682327ed75f362e67e3c8ca659c8289bc072087d52d0393d3a2434c8403` |
| `PAN028mat_vs_PAN027_joint` | 1704327 | `COMPLETED_OUTPUT_PRESENT` | 79,730 | 13,303,455 | `314cc5c1ca4ee2a3e25bb6934084a99d6ee4690ec192e9717fd6f6b84421c785` |

See `summaries/wfmash_jobs.tsv` for paths, checksums, stdout/stderr logs, and
Slurm status.

The exact command logs prove the required literal run shape:

```text
/home/erikg/bin/wfmash -p 95 -t 32 -B /dev/shm/.../tmp /dev/shm/.../target.fa /dev/shm/.../query.fa > /dev/shm/.../<comparison>.literal_p95.paf
```

The logs are:

- `logs/updated_bin_v0.24.2-12-ge040aa10.PAN027pat_vs_PAN011_joint.literal_p95.1704325.command.log`
- `logs/updated_bin_v0.24.2-12-ge040aa10.PAN027mat_vs_PAN010_joint.literal_p95.1704326.command.log`
- `logs/updated_bin_v0.24.2-12-ge040aa10.PAN028mat_vs_PAN027_joint.literal_p95.1704327.command.log`

Each Slurm job staged the prepared target and query FASTAs to `/dev/shm`, used
`-B /dev/shm/.../tmp`, wrote the raw PAF on `/dev/shm`, bgzipped it there, and
copied back only the gzipped raw PAF plus checksum.

Raw PAFs are intentionally ignored under:

`raw_paf/updated_bin_v0.24.2-12-ge040aa10/`

## Summaries

Required summaries:

- `summaries/wfmash_binary.tsv`
- `summaries/wfmash_jobs.tsv`
- `summaries/paf_file_summary.tsv`
- `summaries/candidate_window_support.tsv`
- `summaries/query_grid_filter_manifest.tsv`
- `summaries/query_grid_filter_candidate_window_support.tsv`

Additional useful summary:

- `summaries/input_manifest.tsv`

The candidate-window support summary is posthoc selection from raw whole-genome
PAFs by query contig/window overlap and parsed target chromosome. It reports
chr3 support for both candidate windows and also records other target
chromosomes overlapping those same windows.

## Query-grid SweepGA-compatible post-filter

The `query_grid_filter/` outputs start from the raw updated-bin wfmash `-p 95`
whole-genome PAFs and then apply the same post-alignment normalization used for
the SweepGA/FastGA query-grid comparison runs:

1. verify every raw PAF row has a `cg:Z` CIGAR tag, so exact coordinate-aware
   chopping is possible;
2. chop raw PAF rows with `pafchop-rs --chunk-mode query-grid --overlap 0` at
   10 kb, 5 kb, and 2 kb query-grid lengths;
3. filter each chopped PAF with SweepGA PAF filtering:
   `--num-mappings 1:1 --scaffold-jump 0 --scoring ani --overlap 0`.

The generated directories are:

- `query_grid_filter/chopped_paf_qgrid_l10000_o0/`
- `query_grid_filter/chopped_paf_qgrid_l5000_o0/`
- `query_grid_filter/chopped_paf_qgrid_l2000_o0/`
- `query_grid_filter/filtered_paf_qgrid_l10000_o0/`
- `query_grid_filter/filtered_paf_qgrid_l5000_o0/`
- `query_grid_filter/filtered_paf_qgrid_l2000_o0/`

All chopped and filtered `.paf.gz` outputs were validated with `pigz -t` and
have `.sha256` sidecars. The manifest records raw/chopped/filtered paths, row
counts, commands, and binary paths/checksums:

`summaries/query_grid_filter_manifest.tsv`

Candidate-window chr3 support after the query-grid + SweepGA 1:1 filter is:

| event | comparison | chop length bp | chr3 retained rows | chr3 overlap bp sum | chr3 union bp |
| --- | --- | ---: | ---: | ---: | ---: |
| `PAN027_chr9q_chr3q_PHR_candidate` | `PAN027pat_vs_PAN011_joint` | 10,000 | 0 | 0 | 0 |
| `PAN027_chr9q_chr3q_PHR_candidate` | `PAN027pat_vs_PAN011_joint` | 5,000 | 0 | 0 | 0 |
| `PAN027_chr9q_chr3q_PHR_candidate` | `PAN027pat_vs_PAN011_joint` | 2,000 | 1 | 2,000 | 2,000 |
| `PAN028_chr9q_chr3q_PHR_candidate` | `PAN028mat_vs_PAN027_joint` | 10,000 | 1 | 4,998 | 4,998 |
| `PAN028_chr9q_chr3q_PHR_candidate` | `PAN028mat_vs_PAN027_joint` | 5,000 | 3 | 15,000 | 15,000 |
| `PAN028_chr9q_chr3q_PHR_candidate` | `PAN028mat_vs_PAN027_joint` | 2,000 | 5 | 10,000 | 10,000 |

Reproduce from the repository root with:

```bash
paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/scripts/run_query_grid_filter.sh
```

## Reproduction

From the repository root:

```bash
cd paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin

scripts/capture_wfmash_binary.py
scripts/run_prepare_inputs.sh
WFMASH_PARAMETER_FILTER=literal_p95 \
  WFMASH_CPUS=32 \
  WFMASH_MEM=192G \
  WFMASH_TIME=72:00:00 \
  scripts/submit_wfmash_matrix.sh

# After Slurm completion:
scripts/run_summaries.sh
```

Do not run plain `wfmash`; the script defaults to `/home/erikg/bin/wfmash`.
The task-created PATH still resolves `wfmash` to stale Guix `0.12.5`.

## Validation notes

- All three whole-genome literal `-p 95` Slurm runs completed.
- `sha256sum -c` passed for all three ignored raw PAFs.
- `candidate_window_support.tsv` answers yes for chr3 support in both chr9
  candidate windows.
- `git status --short -- submission` was clean; no `submission/` files were
  modified.
- No Fig5 schematic was created.
