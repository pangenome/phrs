# Fig5 whole-genome sweepGA QA closeout

Date: 2026-06-21

Task: `fig5-whole-genome-sweepga-qa-closeout`

Reviewed inputs:

- Handoff: `paper_prep/_brainstorming/PEDIGREE_SWEEPGA_HANDOFF_2026-06-20.md`
- Corrected package: `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/`
- Evidence review: `paper_prep/_brainstorming/fig5_whole_genome_sweepga_evidence_review/`

## Closeout verdict

The corrected package passes the workflow/provenance QA requirements for the
whole-genome sweepGA/FastGA correction. It is a full whole-genome, joint-parent
run that preserves raw `many:many` PAFs, chops those PAFs before filtering,
uses `/dev/shm` for sweepGA/FastGA scratch, and labels strict `1:1` output as a
diagnostic control.

The corrected direct sweepGA evidence is not acceptable as a new direct Fig5
source for the autosomal chr9q/chr3q candidate structures. The current Fig5
evidence source for those autosomal structures remains the original
graph/untangle schematic artifacts:

- `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/fig5_synteny_recombination_full.pdf`
- `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/selected_segments.tsv`
- `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/event_manifest.tsv`

The corrected direct sweepGA package is acceptable as bounded direct evidence
for the paternal chrX/chrY PAR1 positive-control class only. It should not be
promoted into a direct whole-genome sweepGA schematic for the autosomal
chr9q/chr3q candidates unless a later validated analysis recovers chr3q support
in chopped primary layers.

## Required QA calls

| Issue | QA call | Evidence checked | Closeout assessment |
| --- | --- | --- | --- |
| Full whole-genome input/alignment | PASS | `summaries/input_manifest.tsv`; `config/comparisons.tsv`; raw Slurm jobs `1704307`-`1704309`; `summaries/output_file_manifest.tsv` | All three comparisons use recovered full WashU assembly FASTAs under `/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/`. Each query has 23 full-chromosome haplotype records and each target contains both parental haplotypes. No 500 kb telomeric-window, chromosome-only, or arm-only input is used for the corrected primary outputs. |
| Joint parent target | PASS | `summaries/input_manifest.tsv`; `config/comparisons.tsv`; target header fields | The required target haplotypes are combined before sweepGA filtering: `PAN011#1+PAN011#2` collapsed to `PAN011#joint`, `PAN010#1+PAN010#2` collapsed to `PAN010#joint`, and `PAN027#1+PAN027#2` collapsed to `PAN027#joint`. Filtering therefore chooses jointly across both parental haplotypes. |
| `/dev/shm` sweepGA/FastGA scratch | PASS | `summaries/slurm_jobs.tsv`; `summaries/filter_slurm_jobs.tsv`; `summaries/filter_manifest.tsv`; raw logs `PAN027pat_vs_PAN011_joint.many_many_j0.1704307.out`, `PAN027mat_vs_PAN010_joint.many_many_j0.1704308.out`, `PAN028mat_vs_PAN027_joint.many_many_j0.1704309.out`; filter log `filter_matrix.1704311.out` | Raw jobs record `sweepga_devshm_base=/dev/shm`. Their logs show `scratch=/dev/shm/sg.<job>.*`, staged `q.fa`/`t.fa` under `/dev/shm`, `sweepga --fastga --temp-dir /dev/shm/...`, and FastGA/FAtoGDB/GIXmake working with `/dev/shm` paths. The filter stage also records `/dev/shm` and logs `/dev/shm/sweepga_filter.*` input/output paths. |
| Raw whole-genome `many:many` preservation | PASS | `summaries/slurm_jobs.tsv`; `summaries/output_file_manifest.tsv`; `summaries/chop_manifest.tsv`; upstream raw PAF paths under `/moosefs/erikg/phrs/.wg-worktrees/agent-2616/.../raw_paf/` | Raw `many:many`, scaffold-jump 0 PAFs are first-class outputs for all three comparisons. `output_file_manifest.tsv` records sizes and SHA-256 hashes for the raw PAFs, and those paths were present in the upstream execution worktree during QA. |
| Chopped PAF before filtering | PASS | `summaries/chop_manifest.tsv`; `summaries/filter_manifest.tsv`; `summaries/chop_slurm_jobs.tsv`; `logs/chop_matrix.1704310.out`; `logs/filter_matrix.1704311.out` | `chop_manifest.tsv` maps each raw whole-genome PAF to a chopped PAF before any evidence filters are run. Chopping used `scripts/chop_paf.py`, `chop_length_bp=500000`, `overlap_bp=0`, and preserved total query bp per comparison. `filter_manifest.tsv` shows every filter row reads the corresponding `chopped_paf/*.chopped_l500000_o0.paf.gz` input. |
| `1:1` diagnostic labeling | PASS | `config/filter_matrix.tsv`; `summaries/filter_manifest.tsv`; evidence review `REPORT.md`; `segment_support.tsv`; `chop_filter_assessment.tsv` | Strict `1:1` is explicitly labeled `one_one_chopped` and documented as a diagnostic control expected to suppress multimap PHR signal. The evidence review uses chopped `many:many` and chopped `4:many` as primary layers; `1:1` appears only as `one_one_chopped_diagnostic` in support tables and is not used for Fig5 evidence decisions. |

## Manifest verification notes

`summaries/input_manifest.tsv` contains three corrected comparisons:

- `PAN027pat_vs_PAN011_joint`: query `PAN027#2`, target `PAN011#1+PAN011#2`, 23 query records, 46 target records.
- `PAN027mat_vs_PAN010_joint`: query `PAN027#1`, target `PAN010#1+PAN010#2`, 23 query records, 46 target records.
- `PAN028mat_vs_PAN027_joint`: query `PAN028#1`, target `PAN027#1+PAN027#2`, 23 query records, 46 target records.

All rows state the scope as full whole-genome haplotype FASTA records from
source assembly `.fai` files, with no telomeric-window FASTA. The source FASTAs
are the recovered full WashU v1.1 bgzip+faidx assemblies:

- `/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/PAN010.fa.gz`
- `/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/PAN011.fa.gz`
- `/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/PAN027.fa.gz`
- `/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/PAN028.fa.gz`

`summaries/chop_manifest.tsv` contains one raw-to-chopped mapping per
comparison:

| Comparison | Raw records | Chopped records | Raw query bp | Chopped query bp |
| --- | ---: | ---: | ---: | ---: |
| `PAN027pat_vs_PAN011_joint` | 260 | 6438 | 3133400532 | 3133400532 |
| `PAN027mat_vs_PAN010_joint` | 409 | 6570 | 3161574384 | 3161574384 |
| `PAN028mat_vs_PAN027_joint` | 16477 | 24918 | 6622012399 | 6622012399 |

The current QA worktree does not vendor the multi-GB `raw_paf/`,
`chopped_paf/`, or `filtered_paf/` directories because they are ignored from
git. `summaries/output_file_manifest.tsv` records the relative artifact paths,
byte counts, and SHA-256 hashes. During this closeout, all manifest-listed
output paths were present under the upstream execution worktree:

`/moosefs/erikg/phrs/.wg-worktrees/agent-2616/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/`

## Slurm and scratch verification

Raw alignment jobs:

| Job | Comparison | Mapping | Scratch base | Log evidence |
| --- | --- | --- | --- | --- |
| `1704307` | `PAN027pat_vs_PAN011_joint` | `many:many`, scaffold-jump 0 | `/dev/shm` | Log records `scratch=/dev/shm/sg.1704307.WsiUJL`, staged `q.fa`/`t.fa` under that directory, `sweepga --fastga --temp-dir /dev/shm/sg.1704307.WsiUJL`, FAtoGDB with `TMPDIR=/dev/shm/sg.1704307.WsiUJL`, GIXmake `-P/dev/shm/sg.1704307.WsiUJL`, and FastGA running in that directory. |
| `1704308` | `PAN027mat_vs_PAN010_joint` | `many:many`, scaffold-jump 0 | `/dev/shm` | Same `/dev/shm/sg.1704308.42ZR4D` pattern for sweepGA, staged FASTAs, FAtoGDB, GIXmake, and FastGA. |
| `1704309` | `PAN028mat_vs_PAN027_joint` | `many:many`, scaffold-jump 0 | `/dev/shm` | Same `/dev/shm/sg.1704309.hHKV5j` pattern for sweepGA, staged FASTAs, FAtoGDB, GIXmake, and FastGA. |

Filter job:

- `summaries/filter_slurm_jobs.tsv` records job `1704311` with
  `sweepga_devshm_base=/dev/shm`.
- `summaries/filter_manifest.tsv` records `/dev/shm` for all 15 filter outputs.
- `logs/filter_matrix.1704311.out` records `sweepga_filter_scratch_base=/dev/shm`
  and commands reading `/dev/shm/sweepga_filter.*/input.paf` and writing
  `/dev/shm/sweepga_filter.*/filtered.paf` for `1:1`, `1:many`, `2:many`, and
  `4:many` filters.

## Evidence-source decision for Fig5

The downstream evidence review concluded:

- PAR1 positive control is recovered in chopped whole-genome `many:many` and
  chopped joint `4:many`.
- `PAN027_chr9q_chr3q_PHR_candidate` is not recovered as a chr3q donor
  structure in chopped whole-genome `many:many` or chopped joint `4:many`; the
  primary layers show chr9 same-chromosome support only.
- `PAN028_chr9q_chr3q_PHR_candidate` is not recovered as a chr3q donor
  structure in chopped whole-genome `many:many` or chopped joint `4:many`; the
  primary layers show chr9 same-chromosome support plus one chr16 side-fragment
  row, but no chr3 support.

Therefore, no corrected direct sweepGA source is acceptable yet for the
autosomal Fig5 chr9q/chr3q structures. If Fig5 retains those autosomal events,
their provenance should stay tied to the graph/untangle evidence and they
should be labeled as candidate PHR exchange structures, not confirmed direct
whole-genome sweepGA alignments.

The corrected direct sweepGA outputs may be used only as a bounded positive
control for male paternal PAR1 exchange, with wording consistent with
`paper_prep/_brainstorming/par1_positive_control_literature/REPORT.md`.

## Validation performed for this closeout

- Read the handoff, corrected package README, corrected package configs, and
  evidence review report.
- Verified `summaries/input_manifest.tsv`, `summaries/chop_manifest.tsv`,
  `summaries/slurm_jobs.tsv`, `summaries/filter_slurm_jobs.tsv`, and
  `summaries/filter_manifest.tsv`.
- Verified raw and filter Slurm logs for `/dev/shm` scratch and FastGA/FAtoGDB
  source-adjacent temporary handling.
- Verified `summaries/output_file_manifest.tsv` paths exist in the upstream
  execution worktree.
- Verified the evidence review's `chop_filter_assessment.tsv` and
  `segment_support.tsv` label `1:1` as diagnostic and reserve primary evidence
  calls for chopped `many:many` and chopped `4:many`.
- Confirmed no `submission/` files were modified.
