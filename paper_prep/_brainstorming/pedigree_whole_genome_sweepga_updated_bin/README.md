# Fig5 whole-genome sweepGA updated-binary rerun

Date: 2026-06-21

This package reruns the corrected Fig5 whole-genome joint-parent sweepGA/FastGA
workflow with byte-level provenance for the updated Cargo binary:

`/home/erikg/.cargo/bin/sweepga`

It is a full whole-genome rerun, not a telomeric-window, chromosome-only, or
arm-only test. The three comparisons match
`paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/config/comparisons.tsv`:

- `PAN027pat_vs_PAN011_joint`
- `PAN027mat_vs_PAN010_joint`
- `PAN028mat_vs_PAN027_joint`

## Binary Provenance

`summaries/sweepga_binary.tsv` records:

- explicit path: `/home/erikg/.cargo/bin/sweepga`
- `which`: `/home/erikg/.cargo/bin/sweepga`
- local realpath at package creation: `/export/local/home/erikg/.cargo/bin/sweepga`
- version: `sweepga 0.1.1`
- sha256: `a0d7ac0c3312080d67de96d85cdcad9ce0c5a7e523897109b7f598c186ab85a6`
- `--help` head

Each raw Slurm log also records the explicit path, `which`, compute-node
realpath, version, sha256, `--help`, and exact command. Compute-node `realpath`
resolved to `/home/erikg/.cargo/bin/sweepga`, while the login-node provenance
file records the expected `/export/local/home/erikg/.cargo/bin/sweepga`.

## Workflow Run

Input preparation:

- Slurm job `1704324`
- full recovered WashU assembly FASTAs under
  `/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/`
- output: `summaries/input_manifest.tsv`

Raw alignment:

- Slurm jobs `1704328`, `1704329`, `1704330`
- command shape:
  `sweepga --fastga --num-mappings many:many --scaffold-jump 0 --temp-dir /dev/shm/... --output-file ... QUERY.fa TARGET.fa`
- scratch policy: sweepGA/FastGA scratch under `/dev/shm`; `$SLURM_TMPDIR` is
  not used for sweepGA/FastGA scratch
- outputs: ignored `raw_paf/*.sweepga_many_many_j0.paf.gz`
- summary/checksums: `summaries/paf_file_summary.tsv`

Chopping:

- Final accepted chopping uses `paper_prep/_brainstorming/pafchop-rs/`
- `pafchop-rs` tests passed locally: `4/4`
- binary: `target/release/pafchop`
- binary sha256: `2cf2acbd183e61e07e5fdee1ff9093bf354f790c2f68c6df6b9e5e6b33ae09a9`
- valid Slurm jobs: `1704338`, `1704339`, `1704340`
- chunk length: `10000`
- overlap: `0`
- outputs: ignored `chopped_paf_l10000_o0/*.chopped_l10000_o0.paf.gz`
- compatibility symlinks: ignored `chopped_paf/*.chopped_l10000_o0.paf.gz`
- manifest: `summaries/chop_manifest_l10000_o0.tsv`
- required alias manifest: `summaries/chop_manifest.tsv`

The accepted 10 kb chop counts are:

| comparison | raw rows | 10 kb fragments | raw query bp | chopped query bp |
| --- | ---: | ---: | ---: | ---: |
| `PAN027pat_vs_PAN011_joint` | 260 | 313,481 | 3,133,400,532 | 3,133,400,532 |
| `PAN027mat_vs_PAN010_joint` | 409 | 316,431 | 3,161,574,384 | 3,161,574,384 |
| `PAN028mat_vs_PAN027_joint` | 16,477 | 670,629 | 6,622,012,399 | 6,622,012,399 |

Cancelled/not-used chop attempts:

- Python `scripts/chop_paf.py` Slurm job `1704331` was cancelled after
  01:56:34. Its outputs were moved to
  `quarantine/python_chop_cancelled_1704331/` and are not used as final
  evidence.
- Package-local compiled chopper jobs `1704332`-`1704334` used the rejected
  `500000` length and were cancelled before final outputs. They are not used.
- Duplicate package-level `pafchop-rs` job `1704341` was cancelled after valid
  jobs `1704338`-`1704340` completed. It is not used.

RustyBAM source was checked at `/moosefs/erikg/src/rustybam`; its PAF utilities
do not provide the exact fixed query-axis splitting operation needed here, so
the dedicated `pafchop-rs` tool was used.

Filtering:

- Slurm job `1704342`
- explicit environment: `PAF_CHOP_LENGTH=10000`, `PAF_CHOP_OVERLAP=0`,
  `SWEEPGA=/home/erikg/.cargo/bin/sweepga`, `SWEEPGA_DEVSHM_BASE=/dev/shm`
- input directory: `chopped_paf_l10000_o0/`
- matrix: `many:many`, `1:1`, `1:many`, `2:many`, `4:many`
- `many:many` chopped is a recorded copy of the chopped layer
- filtered outputs: ignored `filtered_paf/*.paf.gz`
- manifest: `summaries/filter_manifest.tsv`

## Candidate-Window Result

`summaries/candidate_window_support.tsv` covers the raw, 10 kb chopped, and
filtered support for the Fig5 PAN027/PAN028 chr9 candidate windows.

PAN027 chr9q/chr3q candidate:

- query: `PAN027#2#chr9:135704825-136204825`
- raw `many:many -j0`: chr9 support only; no chr3 target rows
- 10 kb chopped `many:many`: chr9 support only; no chr3 target rows
- filtered layers (`many:many`, `1:1`, `1:many`, `2:many`, `4:many`): chr9
  support only; no chr3 target rows

PAN028 chr9q/chr3q candidate:

- query: `PAN028#1#chr9:134380985-134880985`
- raw `many:many -j0`: chr9 support plus one chr16 side-fragment row; no chr3
  target rows
- 10 kb chopped `many:many`: chr9 support plus chr16 side-fragment rows; no
  chr3 target rows
- filtered primary layers (`many:many`, `4:many`): chr9 support plus chr16
  side-fragment rows; no chr3 target rows
- stricter diagnostic layers (`1:1`, `1:many`, `2:many`): chr9 support only;
  no chr3 target rows

Conclusion: the updated sweepGA binary does not emit chr3-target rows for the
PAN027 or PAN028 chr9 candidate windows in raw PAFs or after 10 kb chopped
filtering. This updated-binary rerun therefore does not provide direct
whole-genome sweepGA evidence for the autosomal chr9q/chr3q Fig5 candidate
structures.

## Required Summaries

- `summaries/sweepga_binary.tsv`
- `summaries/slurm_jobs.tsv`
- `summaries/chop_manifest.tsv`
- `summaries/chop_manifest_l10000_o0.tsv`
- `summaries/chop_slurm_jobs.tsv`
- `summaries/filter_manifest.tsv`
- `summaries/filter_slurm_jobs.tsv`
- `summaries/paf_file_summary.tsv`
- `summaries/candidate_window_support.tsv`

## Validation

Validated with:

```bash
SWEEPGA=/home/erikg/.cargo/bin/sweepga \
PAF_CHOP_LENGTH=10000 \
PAF_CHOP_OVERLAP=0 \
PACKAGE_DIR=$PWD/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_updated_bin \
paper_prep/_brainstorming/pedigree_whole_genome_sweepga_updated_bin/scripts/validate_outputs.sh
```

The validation regenerated PAF summaries and candidate-window support, verified
all required manifests, and ran `gzip -t` on raw, 10 kb chopped, and filtered
PAFs. No `submission/` files were modified and no Fig5 schematic was created.
