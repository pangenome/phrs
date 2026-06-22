# Fig5 updated-binary direct-alignment evidence review

Date: 2026-06-22

Task: `fig5-updated-binary-direct-alignment-review`

Reviewed updated rerun packages:

- `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/`
- `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_updated_bin/`

Prior evidence reviews used for comparison:

- `paper_prep/_brainstorming/fig5_whole_genome_wfmash_p95_evidence_review/REPORT.md`
- `paper_prep/_brainstorming/fig5_whole_genome_sweepga_evidence_review/REPORT.md`

## Direct answers

Updated wfmash chr3 support: yes. The updated `/home/erikg/bin/wfmash`
binary, run as whole-genome literal `-p 95`, emits chr3-target raw PAF rows
overlapping both PAN027/PAN028 chr9 candidate windows.

Updated sweepGA chr3 support: no. The updated
`/home/erikg/.cargo/bin/sweepga` binary, run as whole-genome FastGA
`many:many -j0`, still emits no chr3-target raw PAF rows for either candidate
window. PAN027 remains chr9-only. PAN028 remains chr9 plus a chr16 side
fragment, not chr3.

The prior wfmash-positive / sweepGA-negative discrepancy is reproduced with
the updated binaries. The updated wfmash result is consistent with the prior
source-built `v0.24.2-0-g774c01ff` review, while the updated sweepGA result is
consistent with the prior `/home/erikg/.cargo/bin/sweepga` review.

This is a technical direct-alignment review only. It does not make a final
Fig5 schematic or biological-mechanism decision.

## Summary table

Full parseable table:
`paper_prep/_brainstorming/fig5_updated_binary_direct_alignment_review/updated_binary_chr3_support_summary.tsv`

| Event | Tool | Binary/version/sha short | Raw chr3 rows | Raw chr3 union bp | Interpretation |
| --- | --- | --- | ---: | ---: | --- |
| `PAN027_chr9q_chr3q_PHR_candidate` | wfmash | `v0.24.2-12-ge040aa10;14a6d5c7` | 2 | 73,000 | Reproduces prior wfmash-positive chr3 support |
| `PAN028_chr9q_chr3q_PHR_candidate` | wfmash | `v0.24.2-12-ge040aa10;14a6d5c7` | 3 | 150,862 | Reproduces prior wfmash-positive chr3 support |
| `PAN027_chr9q_chr3q_PHR_candidate` | sweepGA/FastGA | `sweepga 0.1.1;a0d7ac0c` | 0 | 0 | Remains chr3-negative; chr9-only raw support |
| `PAN028_chr9q_chr3q_PHR_candidate` | sweepGA/FastGA | `sweepga 0.1.1;a0d7ac0c` | 0 | 0 | Remains chr3-negative; chr9 plus chr16 side-fragment raw support |

Minimal updated wfmash chr3 rows are recorded in:
`paper_prep/_brainstorming/fig5_updated_binary_direct_alignment_review/updated_wfmash_minimal_chr3_raw_paf_rows.tsv`

No minimal sweepGA chr3-row table is present because the updated sweepGA raw
PAFs have zero chr3 rows in the candidate windows.

## Binary provenance audit

### wfmash

PASS. The updated rerun uses the explicit binary requested by the task, not a
PATH-resolved stale Guix binary.

`summaries/wfmash_binary.tsv` records:

- explicit path: `/home/erikg/bin/wfmash`
- PATH `which wfmash`: `/home/erikg/.guix-profile/bin/wfmash`
- login-node realpath:
  `/export/local/home/erikg/bin/wfmash-v0.24.2-12-ge040aa10`
- version: `v0.24.2-12-ge040aa10`
- SHA-256:
  `14a6d5c7ac7be8890e904d11121341df118fad0c11193d0e91a9899e18a53d60`
- status: `PASS`

The command logs record the binary actually used on the Slurm nodes. For
example, the PAN027 paternal and PAN028 maternal logs record:

- `wfmash_bin=/home/erikg/bin/wfmash`
- `wfmash_realpath=/home/erikg/bin/wfmash-v0.24.2-12-ge040aa10`
- the same SHA-256 above
- `wfmash_version=v0.24.2-12-ge040aa10`
- command line with `/home/erikg/bin/wfmash -p 95 -t 32 -B /dev/shm/.../tmp`

Reviewed logs:

- `logs/updated_bin_v0.24.2-12-ge040aa10.PAN027pat_vs_PAN011_joint.literal_p95.1704325.command.log`
- `logs/updated_bin_v0.24.2-12-ge040aa10.PAN027mat_vs_PAN010_joint.literal_p95.1704326.command.log`
- `logs/updated_bin_v0.24.2-12-ge040aa10.PAN028mat_vs_PAN027_joint.literal_p95.1704327.command.log`

The updated run is therefore not based on stale Guix `wfmash 0.12.5`.

### sweepGA

PASS. The updated rerun uses the explicit Cargo binary requested by the task.

`summaries/sweepga_binary.tsv` records:

- explicit path: `/home/erikg/.cargo/bin/sweepga`
- PATH `which sweepga`: `/home/erikg/.cargo/bin/sweepga`
- login-node realpath: `/export/local/home/erikg/.cargo/bin/sweepga`
- version: `sweepga 0.1.1`
- SHA-256:
  `a0d7ac0c3312080d67de96d85cdcad9ce0c5a7e523897109b7f598c186ab85a6`
- help text includes the FastGA and `--temp-dir` options

The raw Slurm logs independently record:

- `sweepga_which=/home/erikg/.cargo/bin/sweepga`
- `sweepga_realpath=/home/erikg/.cargo/bin/sweepga`
- the same SHA-256 above
- command line with
  `/home/erikg/.cargo/bin/sweepga --fastga --num-mappings many:many --scaffold-jump 0 --temp-dir /dev/shm/...`

Reviewed raw logs:

- `logs/PAN027pat_vs_PAN011_joint.many_many_j0.1704328.out`
- `logs/PAN027mat_vs_PAN010_joint.many_many_j0.1704329.out`
- `logs/PAN028mat_vs_PAN027_joint.many_many_j0.1704330.out`

The sweepGA conclusion is therefore not based on an unverified path.

## Whole-genome input provenance

PASS for both tools. Both updated packages use the same three joint-parent
comparison IDs:

- `PAN027pat_vs_PAN011_joint`
- `PAN027mat_vs_PAN010_joint`
- `PAN028mat_vs_PAN027_joint`

Both `summaries/input_manifest.tsv` files point to recovered full WashU
assembly FASTAs under:

`/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/`

For each comparison, the query has 23 full-chromosome haplotype records and the
joint-parent target has 46 full-chromosome records from both target haplotypes.
The manifests explicitly label the scope as:

`full whole-genome haplotype FASTA records from source assembly .fai; no telomeric-window FASTA`

No chromosome-only, arm-only, telomeric-only, or window-only substitute FASTA is
used in either updated package.

## Scratch provenance

PASS for both tools.

wfmash `summaries/wfmash_jobs.tsv` records `scratch_base=/dev/shm` for all
three literal `-p 95` jobs. The command logs show `/dev/shm/wfmash...` job
scratch, staged query/target FASTAs on `/dev/shm`, and `-B
/dev/shm/.../tmp`.

sweepGA `summaries/slurm_jobs.tsv` records `sweepga_devshm_base=/dev/shm` for
all three raw `many:many` jobs. The raw logs show `scratch=/dev/shm/sg...` and
the command line passes `--temp-dir /dev/shm/sg...`. The filter-stage logs and
`summaries/filter_manifest.tsv` also record `/dev/shm` for filtered
sweepGA/FastGA scratch.

## Raw PAF evidence first

### Updated wfmash raw PAFs

The updated wfmash raw PAFs are the gzipped files listed in
`summaries/wfmash_jobs.tsv`; checksum validation passed for all three.

Raw candidate-window evidence from
`summaries/candidate_window_support.tsv`:

- `PAN027_chr9q_chr3q_PHR_candidate`, query
  `PAN027#2#chr9:135704825-136204825`: 2 chr3 rows, 73,000 bp summed overlap,
  73,000 bp union overlap.
- `PAN028_chr9q_chr3q_PHR_candidate`, query
  `PAN028#1#chr9:134380985-134880985`: 3 chr3 rows, 152,860 bp summed overlap,
  150,862 bp union overlap.

The minimal chr3 rows have the same coordinates and tags as the prior
source-built `v0.24.2-0-g774c01ff` wfmash review:

- PAN027 rows:
  `PAN027#2#chr9:136166000-136191000` to
  `PAN011#joint#h1_chr3:202510328-202535329`, and
  `PAN027#2#chr9:136117000-136165000` to
  `PAN011#joint#h1_chr3:202461289-202509327`.
- PAN028 rows:
  `PAN028#1#chr9:134764000-134814998` to
  `PAN027#joint#h2_chr3:202431105-202482361`,
  `PAN028#1#chr9:134813000-134863981` to
  `PAN027#joint#h2_chr3:202480363-202531372`, and
  `PAN028#1#chr9:134709119-134760000` to
  `PAN027#joint#h1_chr3:201234027-201284685`.

Interpretation: the chr3 rows do not disappear under the updated local wfmash
binary. The same minimal chr3 rows persist, with the same row counts and union
coverage as the prior source-built wfmash result.

### Updated sweepGA raw PAFs

The updated sweepGA raw PAFs are the gzipped files listed in
`summaries/paf_file_summary.tsv`; checksum validation passed for all three.

Raw candidate-window evidence from
`summaries/candidate_window_support.tsv`:

- `PAN027_chr9q_chr3q_PHR_candidate`, raw `many:many -j0`: chr9 only, 2 rows,
  500,000 bp union coverage; 0 chr3 rows.
- `PAN028_chr9q_chr3q_PHR_candidate`, raw `many:many -j0`: chr9 support, 9
  rows and 498,087 bp union coverage, plus one chr16 side-fragment row covering
  14,898 bp; 0 chr3 rows.

Interpretation: updated sweepGA does not now emit chr3 rows. It remains
consistent with the prior sweepGA review's negative chr3 result for these two
autosomal candidate windows.

## Chopped and filtered sweepGA results

This section is secondary to the raw PAF evidence above.

The updated sweepGA package chops raw PAF rows into 10 kb query-axis fragments
with `pafchop-rs` before running filtered layers:

- `summaries/chop_manifest.tsv`: `chop_length_bp=10000`, `overlap_bp=0`
- `PAN027pat_vs_PAN011_joint`: 260 raw rows to 313,481 chopped rows
- `PAN027mat_vs_PAN010_joint`: 409 raw rows to 316,431 chopped rows
- `PAN028mat_vs_PAN027_joint`: 16,477 raw rows to 670,629 chopped rows

The filter matrix uses the chopped 10 kb PAFs as inputs and records `/dev/shm`
scratch:

- `many:many`
- `1:1`
- `1:many`
- `2:many`
- `4:many`

For PAN027, the chopped `many:many` and all filtered layers remain chr9-only
for the candidate window; there are no chr3 rows.

For PAN028, chopped `many:many`, filtered `many:many`, and filtered `4:many`
retain chr9 plus chr16 side-fragment support; stricter `1:1`, `1:many`, and
`2:many` layers are chr9-only. No chopped or filtered layer emits chr3 rows.

## Biological interpretation boundary

The technical direct-alignment statement is narrow:

- updated wfmash raw whole-genome `-p 95` is chr3-positive for both candidate
  windows;
- updated sweepGA/FastGA raw whole-genome `many:many -j0`, and its 10 kb
  chopped/filtered derivatives, are chr3-negative for both candidate windows.

This discrepancy is a tool-output/provenance finding, not a final mechanism or
Fig5 drawing decision. Any downstream schematic decision should remain separate
and should account for raw many-to-many alignment behavior, same-window chr9
support, PAN028 chr16 side-fragment behavior, haplotype multiplicity, and graph
or untangle evidence outside this review.

No `submission/` files were modified, and no Fig5 schematic was created.

## Validation performed

- Read both updated package READMEs and the two prior evidence review reports.
- Audited binary provenance tables, versions, realpaths, SHA-256 values, and
  command logs for both tools.
- Confirmed both tools used the same three full whole-genome joint-parent
  comparisons and full recovered WashU FASTAs.
- Confirmed wfmash used `-B /dev/shm/.../tmp` and sweepGA/FastGA used
  `--temp-dir /dev/shm/...`.
- Recomputed SHA-256 checksums for all updated raw PAFs listed in the wfmash
  and sweepGA manifests.
- Reviewed raw PAF candidate-window summaries before considering chopped or
  filtered sweepGA layers.
- Compared updated wfmash rows against the prior source-built wfmash review:
  same chr3 support status, row counts, union bp, and minimal row coordinates.
- Compared updated sweepGA status against the prior sweepGA review: still no
  chr3 raw rows; same chr9-only PAN027 and chr9/chr16 PAN028 pattern.
- Checked that the deliverables are confined to
  `paper_prep/_brainstorming/fig5_updated_binary_direct_alignment_review/`.
