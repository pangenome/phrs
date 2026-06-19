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

Primary input:

`/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv`

This table is the existing WashU `odgi untangle` recombination patch output.
The script filters three pre-selected, review-facing examples from that table:

- PAN027 paternal hap2 vs PAN011 father, chrX p PAR1 positive control.
- PAN027 paternal hap2 vs PAN011 father, chr9 q terminal autosomal PHR candidate.
- PAN028 maternal hap1 vs PAN027 mother, chr3 q independent autosomal PHR candidate.

Existing sweepGA/native outputs were inspected during task execution, but this
compact figure is drawn from the recombination patch table so the segment
intervals, scores, and community status match the curated patch calls.

## Regeneration

From the repository root:

```bash
python3 paper_prep/_brainstorming/fig5_par1_phr_candidate_panels/plot_fig5_par1_phr_candidate_panels.py
```

Optional explicit input/output:

```bash
python3 paper_prep/_brainstorming/fig5_par1_phr_candidate_panels/plot_fig5_par1_phr_candidate_panels.py \
  --patches /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv \
  --out-dir paper_prep/_brainstorming/fig5_par1_phr_candidate_panels
```

## Interpretation Boundaries

Panel A is a known male PAR1 X/Y positive control and should be kept separate
from the autosomal PHR interpretation.

Panels B and C are candidate event-level examples compatible with chr3q/chr9q
C3 PHR exchange. They are not presented as clean full crossovers. In Panel B,
the terminal tract is mostly chr3q, while the chr15q segment is marked as a
secondary cross-community fragment and the tiny chr20q tail is treated as
low-confidence. In Panel C, the strongest support is the pair of chr9q h2
intervals near 390-470 kb; chr16q and chr20q fragments are secondary context.

The optional acrocentric/known-system panel was intentionally omitted here:
the finite four-panel asset keeps PAR1 plus the two autosomal chr3q/chr9q
examples legible and avoids a repetitive p-arm-dominated panel.

The asset is for review and presentation triage before any manuscript
integration.
