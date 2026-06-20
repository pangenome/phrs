# Pedigree sweepGA / Fig5 handoff

Date: 2026-06-20

## Current conclusion

The direct sweepGA/Fig5 work so far is not the final analysis. The important
failure mode is now clear:

- The direct sweepGA packages that were just built used 500 kb telomeric-window
  FASTAs, not full assemblies.
- Those window-level alignments are useful as debugging controls, but they
  cannot answer the full recombination-path question.
- Strict `1:1` filtering is also the wrong display/evidence filter for PHRs. It
  suppresses the paralogous/multimap signal by design and can replace the real
  chr9q/chr3q structure with a misleading best local hit.

The next pass must use full assembly FASTAs, or at minimum whole chromosome
sequences if a full whole-genome run is computationally impossible. Do not make
another 500 kb-window schematic and treat it as biological evidence.

## Full assembly inputs found

The full WashU pedigree assemblies appear to be here:

```text
/moosefs/pangenomes/washu_pedigree/PAN010.fa.gz
/moosefs/pangenomes/washu_pedigree/PAN011.fa.gz
/moosefs/pangenomes/washu_pedigree/PAN027.fa.gz
/moosefs/pangenomes/washu_pedigree/PAN028.fa.gz
```

They are about 1.4 GB compressed each and already have `.fai` indexes:

```text
/moosefs/pangenomes/washu_pedigree/PAN010.fa.gz.fai
/moosefs/pangenomes/washu_pedigree/PAN011.fa.gz.fai
/moosefs/pangenomes/washu_pedigree/PAN027.fa.gz.fai
/moosefs/pangenomes/washu_pedigree/PAN028.fa.gz.fai
```

Avoid exploratory `zcat` on these from the head node. Use the `.fai` indexes to
inspect sequence names, and run heavy extraction/alignment through Slurm.

## What was done

### 1. Original graph/untangle schematic

Primary existing comparison figure:

```text
paper_prep/_brainstorming/fig5_synteny_recombination_schematic/fig5_synteny_recombination_full.pdf
paper_prep/_brainstorming/fig5_synteny_recombination_schematic/selected_segments.tsv
```

This remains the clearest current view. It is based on the graph/untangle path
interpretation and should not be overwritten.

### 2. Separate-haplotype direct sweepGA package

Task: `direct-sweepga-parental`

Main artifacts:

```text
paper_prep/_brainstorming/pedigree_direct_sweepga_concordance/
paper_prep/_brainstorming/pedigree_direct_sweepga_concordance/REPORT.md
paper_prep/_brainstorming/pedigree_direct_sweepga_concordance/EVIDENCE_SOURCE_DECISION.md
paper_prep/_brainstorming/pedigree_direct_sweepga_concordance/raw_paf/
paper_prep/_brainstorming/pedigree_direct_sweepga_concordance/filtered_paf/
```

This aligned 46 query telomeric windows against 46 target telomeric windows per
separate parent haplotype. That setup was conceptually wrong for strict
parent-choice display because hap1 and hap2 were filtered separately, then
combined afterward.

Important observation from this package:

- Raw `many:many` and `four_many_noscaffold` retain much more of the PHR
  multimap signal.
- `one_one_noscaffold` crushes the signal.

### 3. Same-format direct `1:1` schematic

Task: `fig5-direct-sweepga-1to1-schematic`

Commit on `main`:

```text
5bbcff1 feat: fig5-direct-sweepga-1to1-schematic (agent-2610)
```

Artifacts:

```text
paper_prep/_brainstorming/fig5_synteny_recombination_sweepga_1to1/
paper_prep/_brainstorming/fig5_synteny_recombination_sweepga_1to1/fig5_synteny_recombination_sweepga_1to1_full.pdf
paper_prep/_brainstorming/fig5_synteny_recombination_sweepga_1to1/selected_segments.sweepga_1to1.tsv
```

This was made in the same schematic style as the original, but it should be
treated as a failed diagnostic control. It plots direct `1:1` PAF rows and shows
why strict `1:1` is not the right representation for PHR recombination.

Observed failure:

- PAN027 chr9q/chr3q loses the chr3q donor structure.
- A misleading chr1p block appears in the strict direct result.
- PAR1 looks confusing because local PAF rows are not a single inherited path.

### 4. Joint-parent window-level direct sweepGA package

Task: `fig5-joint-parent-sweepga`

Commit on `main`:

```text
906b00c feat: fig5-joint-parent-sweepga (agent-2613)
```

Artifacts:

```text
paper_prep/_brainstorming/pedigree_direct_sweepga_joint_parent/
paper_prep/_brainstorming/fig5_synteny_recombination_joint_parent/
```

This fixed the parent-choice setup at the 500 kb-window level:

- One child query FASTA.
- One combined target FASTA containing both parental haplotypes.
- Joint filtering across both parent haplotypes.

Successful Slurm jobs:

```text
1704274 PAN027pat_vs_PAN011_joint
1704275 PAN027mat_vs_PAN010_joint
1704276 PAN028mat_vs_PAN027_joint
```

However, this is still window-limited. It is not the final answer.

Key window-level result:

- Joint `1:1` still crushes PHR donor signal.
- Joint `4:many` is the better visible window-level comparison.
- Raw `many:many` remains the audit/provenance layer.
- For selected events, joint `4:many` and raw `many:many` give the same rendered
  segment set.
- The joint window-level run still does not recover the PAN027 chr3q donor seen
  in the graph/untangle schematic.

## PAR1 literature status

Task: `research-par1-positive`

Artifact:

```text
paper_prep/_brainstorming/par1_positive_control_literature/REPORT.md
```

Conclusion:

- A paternal chrX/chrY PAR1 event is a legitimate internal positive-control
  class for detecting interchromosomal subtelomeric recombination.
- Say male PAR1 specifically. Do not say generic PAR recombination is obligate.
- PAR1 validates detection in a known context; it does not prove autosomal PHR
  events use the same mechanism.

Good wording:

```text
Positive-control paternal chrX/chrY PAR1 exchange, shown alongside candidate
autosomal PHR recombination events.
```

## What went wrong conceptually

1. The direct sweepGA runs were assumed to be whole-genome-ish, but they were
   actually 500 kb telomeric-window alignments.
2. Separate parent-haplotype filtering was not a joint parent-choice problem.
3. Strict `1:1` is not a path caller. It is a local orthogonal mapping filter.
   In PHR sequence it removes the multimap signal we need to inspect.
4. The same-format direct `1:1` schematic is visually useful only as evidence
   that `1:1` is the wrong filter. It should not be promoted into the paper.

## Next correct task

Create a new WG task for whole-genome direct sweepGA:

```text
Task: Whole-genome joint-parent sweepGA for Fig5 pedigree events
```

Required inputs:

```text
/moosefs/pangenomes/washu_pedigree/PAN010.fa.gz
/moosefs/pangenomes/washu_pedigree/PAN011.fa.gz
/moosefs/pangenomes/washu_pedigree/PAN027.fa.gz
/moosefs/pangenomes/washu_pedigree/PAN028.fa.gz
```

Required comparisons:

```text
PAN027 paternal hap2 query vs PAN011 both parental haplotypes
PAN027 maternal hap1 query vs PAN010 both parental haplotypes
PAN028 maternal hap1 query vs PAN027 both parental haplotypes
```

Important implementation requirements:

- Inspect `.fai` indexes to determine exact sequence-name patterns.
- Extract the transmitting child haplotype as query and both parent haplotypes
  as one combined target. Do this from full assemblies, not the 500 kb telomere
  FASTA.
- Run all heavy alignment on Slurm, not on the head node.
- Use node-local scratch (`$SLURM_TMPDIR` or `/tmp`) for sweepGA/FastGA temp
  files. Earlier `/dev/shm` failed in `FAtoGDB` cleanup.
- Preserve raw `many:many` PAFs as first-class artifacts.
- Apply filters jointly across the combined parental target:
  `1:1`, `1:many`, `2:many`, `4:many`, and raw/`many:many`.
- Evaluate whether `4:many` or raw `many:many` best preserves the old
  untangle-visible PAR1 and chr9q/chr3q structures.
- Only after inspecting the whole-genome results should a new same-format Fig5
  schematic be generated.

Strong recommendation:

- Treat `1:1` as a diagnostic control only.
- Treat raw `many:many` and `4:many` as the likely useful evidence layers.
- Do not use 500 kb local-window offsets in any figure intended to answer the
  whole-genome recombination-path question.

## Suggested WG prompt

```text
Run the corrected whole-genome direct sweepGA experiment for the WashU pedigree
Fig5 events. The previous direct and joint-parent packages used 500 kb
telomeric-window FASTAs and are not adequate for the full recombination-path
question. Use the full assembly FASTAs in /moosefs/pangenomes/washu_pedigree:
PAN010.fa.gz, PAN011.fa.gz, PAN027.fa.gz, and PAN028.fa.gz. Inspect the .fai
indexes for exact sequence names. Build one query FASTA for each transmitting
child haplotype and one combined target FASTA containing both parental
haplotypes for the relevant parent. Required comparisons: PAN027 paternal hap2
vs PAN011 both haplotypes; PAN027 maternal hap1 vs PAN010 both haplotypes;
PAN028 maternal hap1 vs PAN027 both haplotypes. Submit heavy sweepGA/FastGA
jobs through Slurm only, using node-local scratch. Preserve raw many:many PAFs.
Apply joint filters 1:1, 1:many, 2:many, 4:many, and keep raw many:many. Report
whether whole-genome raw/4:many recovers the PAR1 positive control and the
chr9q/chr3q PHR candidate structures seen in the graph/untangle schematic.
Generate same-format schematic comparisons only after validating the whole-genome
PAF support. Do not overwrite existing Fig5 schematic directories. Commit and
push all scripts, summaries, and review PDFs with WG provenance.
```

## Files to compare visually right now

Original graph/untangle view:

```text
paper_prep/_brainstorming/fig5_synteny_recombination_schematic/fig5_synteny_recombination_full.pdf
```

Failed separate-hap direct `1:1` diagnostic:

```text
paper_prep/_brainstorming/fig5_synteny_recombination_sweepga_1to1/fig5_synteny_recombination_sweepga_1to1_full.pdf
```

Window-level joint-parent diagnostics:

```text
paper_prep/_brainstorming/fig5_synteny_recombination_joint_parent/fig5_synteny_recombination_joint_parent_1to1_full.pdf
paper_prep/_brainstorming/fig5_synteny_recombination_joint_parent/fig5_synteny_recombination_joint_parent_4many_full.pdf
```

Use these only for debugging/filter behavior. The next real analysis should be
whole-genome or whole-chromosome.

