# Fig5 whole-genome sweepGA evidence review

## Verdict

The corrected package is a real whole-genome sweepGA/FastGA run, not a 500 kb
window, chromosome-only, arm-only, or other substitute run. It passes the
workflow/provenance checks for full-genome inputs, `/dev/shm` scratch, raw
whole-genome PAF chopping before filtering, and chopped-filter evidence policy.

The evidence result is mixed and does not validate a new direct-alignment
schematic for Fig5:

- `PAR1_XY_positive_control`: recovered in chopped whole-genome `many:many` and
  chopped joint `4:many`.
- `PAN027_chr9q_chr3q_PHR_candidate`: not recovered as a chr3q donor structure
  in chopped whole-genome `many:many` or chopped joint `4:many`; the primary
  layers show chr9 same-chromosome support only.
- `PAN028_chr9q_chr3q_PHR_candidate`: not recovered as a chr3q donor structure
  in chopped whole-genome `many:many` or chopped joint `4:many`; the primary
  layers show chr9 same-chromosome support plus one chr16 side-fragment row, but
  no chr3 support.

No optional sibling schematic directory was created. A direct whole-genome
sweepGA schematic would currently be defensible for the paternal PAR1 positive
control only, not for the autosomal chr9q/chr3q candidate structures.

## Reviewed inputs

Corrected whole-genome package:

`paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/`

Context reviewed:

- Original schematic event definitions:
  `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/event_manifest.tsv`
- Original selected schematic segments:
  `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/selected_segments.tsv`
- Prior direct/window evidence decision:
  `paper_prep/_brainstorming/pedigree_direct_sweepga_concordance/EVIDENCE_SOURCE_DECISION.md`
- PAR1 wording/control report:
  `paper_prep/_brainstorming/par1_positive_control_literature/REPORT.md`

The current review worktree does not contain `raw_paf/`, `chopped_paf/`, or
`filtered_paf/` directories because the multi-GB PAF outputs are ignored from
git. I therefore used the exact artifact paths recorded in the package
manifests and verified that those files were still present under the upstream
artifact worktree. The package's `summaries/output_file_manifest.tsv` records
the relative artifact paths, sizes, and SHA-256 checksums.

## Whole-genome scope

`summaries/input_manifest.tsv` shows three full whole-genome comparisons:

- `PAN027pat_vs_PAN011_joint`: query `PAN027#2`, target `PAN011#1+PAN011#2`
  collapsed to `PAN011#joint`.
- `PAN027mat_vs_PAN010_joint`: query `PAN027#1`, target `PAN010#1+PAN010#2`
  collapsed to `PAN010#joint`.
- `PAN028mat_vs_PAN027_joint`: query `PAN028#1`, target `PAN027#1+PAN027#2`
  collapsed to `PAN027#joint`.

Each row records 23 full-chromosome query records and joint-parent target
records selected from recovered full WashU assemblies under:

`/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/`

The manifest explicitly labels the scope as "full whole-genome haplotype FASTA
records from source assembly .fai; no telomeric-window FASTA". The package does
not use the historical trimmed input:

`/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/washu.1Mb.telo_500kb_trimmed.fa.gz`

## `/dev/shm` scratch verification

The raw whole-genome sweepGA/FastGA Slurm jobs are recorded in
`summaries/slurm_jobs.tsv`:

- job `1704307`, `PAN027pat_vs_PAN011_joint`, `many:many`, `/dev/shm`
- job `1704308`, `PAN027mat_vs_PAN010_joint`, `many:many`, `/dev/shm`
- job `1704309`, `PAN028mat_vs_PAN027_joint`, `many:many`, `/dev/shm`

The raw-job logs verify that both sweepGA and FastGA scratch were under
`/dev/shm`. For example, the `PAN027pat_vs_PAN011_joint` log records:

- `scratch=/dev/shm/sg.1704307.WsiUJL`
- `scratch_fastga_sources=/dev/shm/sg.1704307.WsiUJL/q.fa,/dev/shm/sg.1704307.WsiUJL/t.fa`
- `sweepga --fastga ... --temp-dir /dev/shm/sg.1704307.WsiUJL ...`
- FastGA/FAtoGDB/GIXmake operating with `/dev/shm/sg.1704307.WsiUJL` paths

The same pattern appears in the `1704308` and `1704309` logs. The filter stage
also used `/dev/shm`: `summaries/filter_manifest.tsv` records
`sweepga_devshm_base=/dev/shm` for every filter row, and
`logs/filter_matrix.1704311.out` shows scratch directories such as
`/dev/shm/sweepga_filter.PAN027pat_vs_PAN011_joint.ztk68pv1/`.

## Chopping and filter order

The package uses raw whole-genome `many:many` PAFs as first-class outputs, then
chops those raw PAFs before any evidence filters are run. The exact raw-to-chop
mapping is in `summaries/chop_manifest.tsv`:

- `PAN027pat_vs_PAN011_joint`: 260 raw rows to 6438 chopped rows.
- `PAN027mat_vs_PAN010_joint`: 409 raw rows to 6570 chopped rows.
- `PAN028mat_vs_PAN027_joint`: 16477 raw rows to 24918 chopped rows.

The exact whole-genome raw PAFs are:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2616/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/raw_paf/PAN027pat_vs_PAN011_joint.sweepga_many_many_j0.paf.gz`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2616/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/raw_paf/PAN027mat_vs_PAN010_joint.sweepga_many_many_j0.paf.gz`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2616/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/raw_paf/PAN028mat_vs_PAN027_joint.sweepga_many_many_j0.paf.gz`

The exact chopped PAFs used as filter inputs are:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2616/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/chopped_paf/PAN027pat_vs_PAN011_joint.chopped_l500000_o0.paf.gz`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2616/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/chopped_paf/PAN027mat_vs_PAN010_joint.chopped_l500000_o0.paf.gz`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2616/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/chopped_paf/PAN028mat_vs_PAN027_joint.chopped_l500000_o0.paf.gz`

Chopping before filtering was required because raw whole-genome PAF rows can
span multi-Mb intervals. Running `1:1`, `4:many`, or similar filters directly
on those long rows allows merged intervals to dominate the path/similarity
metric and can erase local subtelomeric alternatives. The implemented splitter
breaks each raw PAF row into deterministic query-axis fragments before
filtering, linearly interpolating target coordinates and scaling match/aligned
counts.

The parameters are:

- `PAF_CHOP_LENGTH=500000`
- `PAF_CHOP_OVERLAP=0`

These parameters are appropriate for bounding raw whole-genome alignments and
for matching the historical 500 kb terminal-window display scale. They are
coarse for sub-50 kb donor tracts, so absence of chr3q support in this review
should be interpreted as a negative result for this corrected whole-genome run,
not as a general proof that local chr3q homology cannot exist.

`summaries/filter_manifest.tsv` confirms that every filter used the chopped PAF
input. The evidence layers reviewed here are:

- Primary: chopped `many:many`
- Primary: chopped `4:many`
- Diagnostic only: chopped `1:1`
- Diagnostic/control only: unchopped raw `many:many`

## Evidence results

The per-event support table is:

`paper_prep/_brainstorming/fig5_whole_genome_sweepga_evidence_review/segment_support.tsv`

### PAR1 positive control

The paternal chrX/chrY PAR1 positive-control class is recovered. In both
primary layers, the PAN027 paternal chrX query window has:

- chrY donor support: 1 overlapping row, 144096 bp covered.
- chrX same-chromosome support: 2 overlapping rows, 362564 bp covered.

Primary PAFs:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2616/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/filtered_paf/PAN027pat_vs_PAN011_joint.many_many_chopped.paf.gz`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2616/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/filtered_paf/PAN027pat_vs_PAN011_joint.four_many_chopped.paf.gz`

This supports the PAR1 report's bounded wording: the event is a paternal
chrX/chrY PAR1 positive-control class for detecting a known inter-chromosomal
subtelomeric recombination context. It should not be used as mechanistic proof
for autosomal PHR recombination.

### PAN027 chr9q/chr3q candidate

The corrected whole-genome primary layers do not recover the PAN027 autosomal
candidate as a chr9q/chr3q donor structure. The reviewed query interval was
`PAN027#2#chr9:135704825-136204825`.

Both chopped `many:many` and chopped `4:many` contain chr9 support across the
window:

- chr9 same-chromosome support: 3 rows, 530305 bp summed overlap, 500000 bp
  query-union coverage.

No chr3 rows overlapped the query window in either primary layer. The chopped
`1:1` diagnostic layer also showed chr9 only.

Primary PAFs:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2616/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/filtered_paf/PAN027pat_vs_PAN011_joint.many_many_chopped.paf.gz`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2616/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/filtered_paf/PAN027pat_vs_PAN011_joint.four_many_chopped.paf.gz`

Decision: do not use this corrected whole-genome direct sweepGA run as direct
alignment evidence for the PAN027 chr3q donor portion of the schematic.

### PAN028 chr9q/chr3q candidate

The corrected whole-genome primary layers do not recover the PAN028 autosomal
candidate as a chr9q/chr3q donor structure. The reviewed query interval was
`PAN028#1#chr9:134380985-134880985`.

Both chopped `many:many` and chopped `4:many` contain:

- chr9 same-chromosome support: 9 rows, 788847 bp summed overlap, 498087 bp
  query-union coverage.
- chr16 side-fragment support: 1 row, 14898 bp covered.

No chr3 rows overlapped the query window in either primary layer. The chopped
`1:1` diagnostic layer showed chr9 only.

Primary PAFs:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-2616/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/filtered_paf/PAN028mat_vs_PAN027_joint.many_many_chopped.paf.gz`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-2616/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_joint_parent/filtered_paf/PAN028mat_vs_PAN027_joint.four_many_chopped.paf.gz`

Decision: do not use this corrected whole-genome direct sweepGA run as direct
alignment evidence for the PAN028 chr3q donor portion of the schematic.

## Consequences for Fig5

The prior schematic and direct/window diagnostics remain useful context, but
the corrected whole-genome evidence changes the direct-alignment interpretation.
The PAR1 row validates the expected positive-control class. The two autosomal
rows do not validate as direct whole-genome sweepGA chr9q/chr3q structures under
the required chopped primary layers.

Recommended downstream policy:

1. Keep PAR1 as an internal positive-control benchmark if Fig5 includes a
   direct whole-genome sweepGA evidence row.
2. Do not redraw the autosomal chr9q/chr3q candidates from the corrected
   whole-genome direct PAFs unless a later validated analysis recovers chr3q
   support in chopped primary layers.
3. Treat `1:1` strictly as diagnostic; do not use it as the basis for evidence
   decisions or direct schematic geometry.
4. If autosomal candidates remain in Fig5, keep their provenance tied to the
   graph/untangle evidence and label them as candidate PHR exchange structures,
   not as confirmed direct whole-genome sweepGA alignments.
