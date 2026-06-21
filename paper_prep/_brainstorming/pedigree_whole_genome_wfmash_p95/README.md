# Fig5 whole-genome wfmash `-p 95` recovery test

Date: 2026-06-21

Task: `fig5-whole-genome-wfmash-p95`

This package records a raw whole-genome wfmash homology-recovery test for the
same Fig5 pedigree comparisons used by the corrected whole-genome sweepGA
package:

- `PAN027pat_vs_PAN011_joint`
- `PAN027mat_vs_PAN010_joint`
- `PAN028mat_vs_PAN027_joint`

This is not a manuscript update. No `submission/` files were edited and no
Fig5 schematic was created.

## Inputs

Comparison definitions are copied from:

- `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/config/comparisons.tsv`
- `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/summaries/input_manifest.tsv`

The source assemblies are the recovered readable full WashU v1.1 assemblies:

- `/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/PAN010.fa.gz`
- `/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/PAN011.fa.gz`
- `/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/PAN027.fa.gz`
- `/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/PAN028.fa.gz`

`scripts/prepare_inputs.py` extracts whole-haplotype FASTAs from the source
assembly `.fai` records and collapses both target haplotypes to a joint target
prefix. Each query input has 23 full-chromosome haplotype records and each
target input has 46 full-chromosome records. No chromosome-only, arm-only,
telomeric-window, or 500 kb-window FASTA was used as a wfmash input.

Large generated input FASTAs live in ignored `inputs/`. For this run they were
hard-linked from the prior corrected sweepGA package and
`scripts/run_prepare_inputs.sh` was rerun so `summaries/input_manifest.tsv`
records paths under this package.

## Binary

The primary run uses upstream wfmash `v0.24.2`, built from source:

- Release: <https://github.com/waveygang/wfmash/releases/tag/v0.24.2>
- Source: <https://github.com/waveygang/wfmash/archive/refs/tags/v0.24.2.tar.gz>
- Tag commit: `774c01ffb9df010d6b529520033ed7dce0cb95d5`
- Reported version: `v0.24.2-0-g774c01ff`
- Built binary checksum: `00fbdb08b17ffa05142de87c0a0389c1debf317422027c4be22388ee90ef05a4`

The local Guix profile binary (`/home/erikg/.guix-profile/bin/wfmash`) was an
older `0.12.5-1+0222f7c` build and is treated only as legacy diagnostic output.
The GitHub `wfmash-v0.24.2-linux-x86_64` binary was downloaded and checksum
checked but rejected because it requires newer system glibc versions than this
host provides. See:

- `summaries/wfmash_binary.tsv`
- `logs/wfmash_v0.24.2_build_notes.md`
- `logs/wfmash_v0.24.2.version.txt`
- `logs/wfmash_v0.24.2.help.txt`
- `logs/wfmash_v0.24.2_release_binary_failure.txt`

The current help confirms the CLI distinction needed for this task:
`-s` is sketch size and `-w` is window size. Therefore the sensitive parameter
uses `-w 1k`.

## Parameters

`config/wfmash_parameters.tsv` defines:

| Parameter set | wfmash options | Status |
| --- | --- | --- |
| `literal_p95` | `-p 95` | Completed for all three comparisons and used for conclusions. |
| `permissive_p95` | `-p 95 -w 1k -l 1k -n 50 -f -M` | Submitted, then cancelled after literal `-p 95` recovered clear chr3 support at both candidate windows. |

No scaffolding filters, one-to-one filters, or downstream evidence filters were
applied before summarizing the raw PAF rows.

## Slurm Execution

Primary jobs were submitted with:

```bash
WFMASH_CPUS=32 WFMASH_MEM=192G WFMASH_TIME=72:00:00 WFMASH_SCRATCH_BASE=/dev/shm \
WFMASH_RUN_LABEL=current_v0.24.2 \
WFMASH_BIN=/moosefs/erikg/phrs/.wg-worktrees/agent-2630/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95/tools/wfmash-v0.24.2-built/wfmash \
  paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95/scripts/submit_wfmash_matrix.sh
```

Each Slurm task staged both whole-genome FASTAs under `/dev/shm`, used
`wfmash -B` under the same node-local scratch directory, wrote raw PAF on the
compute node, bgzipped it, and copied the bgzipped PAF back to ignored
`raw_paf/current_v0.24.2/`.

Command shape for completed literal jobs:

```bash
$WFMASH_BIN -p 95 -t $SLURM_CPUS_PER_TASK \
  -B /dev/shm/wfmash.$SLURM_JOB_ID.current_v0.24.2.$COMPARISON.literal_p95/tmp \
  /dev/shm/wfmash.$SLURM_JOB_ID.current_v0.24.2.$COMPARISON.literal_p95/target.fa \
  /dev/shm/wfmash.$SLURM_JOB_ID.current_v0.24.2.$COMPARISON.literal_p95/query.fa \
  > /dev/shm/wfmash.$SLURM_JOB_ID.current_v0.24.2.$COMPARISON.literal_p95/$COMPARISON.literal_p95.paf
```

Command shape for the submitted permissive jobs:

```bash
$WFMASH_BIN -p 95 -w 1k -l 1k -n 50 -f -M -t $SLURM_CPUS_PER_TASK \
  -B /dev/shm/wfmash.$SLURM_JOB_ID.current_v0.24.2.$COMPARISON.permissive_p95/tmp \
  /dev/shm/wfmash.$SLURM_JOB_ID.current_v0.24.2.$COMPARISON.permissive_p95/target.fa \
  /dev/shm/wfmash.$SLURM_JOB_ID.current_v0.24.2.$COMPARISON.permissive_p95/query.fa \
  > /dev/shm/wfmash.$SLURM_JOB_ID.current_v0.24.2.$COMPARISON.permissive_p95/$COMPARISON.permissive_p95.paf
```

Final primary job status is in `summaries/wfmash_jobs.tsv`:

| Comparison | Parameter set | Job ID | Final status |
| --- | --- | ---: | --- |
| `PAN027pat_vs_PAN011_joint` | `literal_p95` | 1704318 | completed, raw PAF present |
| `PAN027mat_vs_PAN010_joint` | `literal_p95` | 1704320 | completed, raw PAF present |
| `PAN028mat_vs_PAN027_joint` | `literal_p95` | 1704322 | completed, raw PAF present |
| `PAN027pat_vs_PAN011_joint` | `permissive_p95` | 1704319 | cancelled after literal support was confirmed |
| `PAN027mat_vs_PAN010_joint` | `permissive_p95` | 1704321 | cancelled after literal support was confirmed |
| `PAN028mat_vs_PAN027_joint` | `permissive_p95` | 1704323 | cancelled after literal support was confirmed |

The older Guix `0.12.5` diagnostic jobs and PAFs are separated under
`raw_paf/legacy_guix_0.12.5/` and `summaries/legacy_guix_0.12.5_wfmash_jobs.tsv`.
They are not used for the conclusion below.

## Raw PAF Outputs

Raw PAFs are intentionally ignored because they can be large. The manifest
records absolute paths, sizes, and checksums.

Completed primary raw PAFs:

| Comparison | Raw PAF size | SHA256 |
| --- | ---: | --- |
| `PAN027pat_vs_PAN011_joint` | 2,835,280 bytes | `96312a525e02e7e3ed0f0a73ee3fda22df2977e0c213789084952bcf8ed2c0af` |
| `PAN027mat_vs_PAN010_joint` | 2,694,821 bytes | `6e63ec7bbc172f1ec334e6de6c8eb018876a9010ee3ad3c2522522bad3da6b8e` |
| `PAN028mat_vs_PAN027_joint` | 13,284,081 bytes | `c20e9d49a958faf89ab71198ddd3d961e2d23a86803e9865f01d9e5d237db2c4` |

`summaries/paf_file_summary.tsv` gives target-chromosome distributions for
each completed raw PAF. Total raw PAF row counts are:

| Comparison | Parameter set | Total rows | chr3 target rows |
| --- | --- | ---: | ---: |
| `PAN027pat_vs_PAN011_joint` | `literal_p95` | 66,592 | 4,485 |
| `PAN027mat_vs_PAN010_joint` | `literal_p95` | 66,634 | 4,479 |
| `PAN028mat_vs_PAN027_joint` | `literal_p95` | 79,733 | 5,765 |

## Candidate Windows

`config/candidate_windows.tsv` contains the autosomal Fig5 candidate windows
used by the corrected sweepGA review:

| Event | Query interval | Expected target chromosome |
| --- | --- | --- |
| `PAN027_chr9q_chr3q_PHR_candidate` | `PAN027#2#chr9:135704825-136204825` | `chr3` |
| `PAN028_chr9q_chr3q_PHR_candidate` | `PAN028#1#chr9:134380985-134880985` | `chr3` |

Candidate-window slicing was performed only after raw whole-genome PAFs were
produced. `summaries/candidate_window_support.tsv` counts raw PAF rows whose
query interval overlaps each candidate window, grouped by target chromosome.

Result:

| Event | Comparison | Parameter set | chr3 rows overlapping candidate | Union query coverage on chr3 |
| --- | --- | --- | ---: | ---: |
| `PAN027_chr9q_chr3q_PHR_candidate` | `PAN027pat_vs_PAN011_joint` | `literal_p95` | 2 | 73,000 bp |
| `PAN028_chr9q_chr3q_PHR_candidate` | `PAN028mat_vs_PAN027_joint` | `literal_p95` | 3 | 150,862 bp |

Conclusion: current upstream wfmash `v0.24.2` with literal whole-genome
`-p 95` emits chr3-target raw PAF rows overlapping both Fig5 chr9 candidate
windows before any downstream filtering. The permissive configuration was
therefore not needed for the homology-recovery decision and was cancelled after
the literal evidence was verified.

## Reproduction

From the repository root:

```bash
cd paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95
scripts/build_wfmash_v0242.sh
scripts/run_prepare_inputs.sh
WFMASH_CPUS=32 WFMASH_MEM=192G WFMASH_TIME=72:00:00 WFMASH_SCRATCH_BASE=/dev/shm \
WFMASH_RUN_LABEL=current_v0.24.2 \
WFMASH_BIN=$PWD/tools/wfmash-v0.24.2-built/wfmash \
  scripts/submit_wfmash_matrix.sh
scripts/run_summaries.sh
```
