# Fig5 PAR1/PHR Candidate Panels

This directory contains an experimental candidate asset pack for reviewing a
compact Figure 5-style pedigree recombination panel. It does not modify
`submission/paper.tex`, `submission/fig/MainFigures/Fig5_pedigree_untangle.pdf`,
or any bibliography file.

## Outputs

- `fig5_par1_phr_candidate_panels.pdf`
- `fig5_par1_phr_candidate_panels.svg`
- `panel_event_summary.tsv`
- `plot_fig5_par1_phr_candidate_panels.py`

## Provenance

Drawing input:

`paper_prep/_brainstorming/fig5_sweepga_1to1_redraw/conservative_segments.tsv`

The plotted rectangles come from `conservative_segments.tsv`, the strict
`nb=1` sweepGA 1:1 no-scaffold primary-path table produced by the companion
redraw task. This correction does not plot nth-best, multimap, permissive
secondary, or non-primary PAF rows.

Annotation/provenance input:

`/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv`

This table is the existing WashU `odgi untangle` recombination patch output.
The script uses it only to recover community labels, score summaries, PHR/PAR
status, and interpretive labels when a conservative segment has the same query,
interval, donor arm, and donor haplotype.

The script filters three pre-selected, review-facing examples:

- PAN027 paternal hap2 vs PAN011 father, chrX p PAR1 positive control.
- PAN027 paternal hap2 vs PAN011 father, chr9 q terminal autosomal PHR candidate.
- PAN028 maternal hap1 vs PAN027 mother, chr3 q independent autosomal PHR candidate.

## Coordinate Convention

Axes, callouts, and `panel_event_summary.tsv` use 0-based half-open native
assembly coordinates parsed from source names such as
`PAN027#2#chr9.paternal:135704825-136204824_chr9_qarm`. A local segment
`[a,b)` is displayed as `chr:(start+a)-(start+b)`. For example, local
`[446944,472441)` in the PAN027 paternal chr9q source window is displayed as
`chr9:136,151,769-136,177,266`.

No CHM13 projection or liftover table is used here. The source names are native
sample assembly windows, so the figure deliberately labels the coordinates as
native assembly coordinates, not CHM13 coordinates.

Donor/source target names are also native assembly windows. The conservative
summary table does not carry exact target-side segment start/end fields, so
`panel_event_summary.tsv` records each donor source window but marks exact donor
segment intervals as not recovered in this lightweight presentation correction.

## Regeneration

From the repository root:

```bash
python3 paper_prep/_brainstorming/fig5_par1_phr_candidate_panels/plot_fig5_par1_phr_candidate_panels.py
```

Optional explicit input/output:

```bash
python3 paper_prep/_brainstorming/fig5_par1_phr_candidate_panels/plot_fig5_par1_phr_candidate_panels.py \
  --segments paper_prep/_brainstorming/fig5_sweepga_1to1_redraw/conservative_segments.tsv \
  --patches /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv \
  --out-dir paper_prep/_brainstorming/fig5_par1_phr_candidate_panels
```

## Interpretation Boundaries

Panel A is a known male PAR1 X/Y positive control and should be kept separate
from the autosomal PHR interpretation.

Panels B and C are candidate event-level examples compatible with chr3q/chr9q
C3 PHR exchange. They are not presented as clean full crossovers. In Panel B,
the terminal tract is mostly chr3q, while the chr15q segment is marked as a
smaller side fragment within the single selected 1:1 path, and the tiny chr20q
tail is treated as low-confidence. In Panel C, strict primary-path drawing
shows chr7p, chr16q, and chr20q side fragments in the selected chr3q window;
previous permissive chr9q patch calls are not drawn as alternate alignments.

The optional acrocentric/known-system panel was intentionally omitted here:
the finite four-panel asset keeps PAR1 plus the two autosomal chr3q/chr9q
examples legible and avoids a repetitive p-arm-dominated panel.

The asset is for review and presentation triage before any manuscript
integration.
