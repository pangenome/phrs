## Title

Interchromosomal similarities — n-chromosomes per region (HPRCv2)

## Bullets

- Orange traces: number of unique chromosomes matching each 100 kb window across CHM13; CEN / PAR / PHR / XTR painted as background reference.
- Spikes are not confined to PARs — every subtelomere lights up, plus centromeres and the acrocentric short arms.
- **PHR scale: median 105 kb, mean 144 kb, range 5–500 kb** (15,668 PHR sequences across 41/48 chromosome arms).
- That places typical subtelomeric pseudohomology on the same length scale as PAR2 (~334 kb) — and present at nearly every chromosome end, not just X/Y.
- Read the plot as: a PAR2-class exchange landscape replicated 41 times across the human genome.

## Primary figure

`p_num_chromosomes_wide.pdf` (worktree root) — produced by `plot-impg-coverage.R` (`p_num_chromosomes_wide`, line ~321; saved line ~552). This is the v1 slide 5 figure preserved verbatim.

Cross-reference for synthesizer: `paper_prep/figures/fig1/figure_fig1.pdf` panel **1b** is the publication-quality version of the same view (genome-wide num-chromosome heatmap, 100 kb windows). v1 traces are easier to read at conference distance, so keep the v1 plot here; reserve fig1b for the manuscript.

Optional annotation overlay (if synthesizer wants the PHR scale rendered into the figure rather than spoken — drop into the existing `plot-impg-coverage.R` after `p_num_chromosomes_wide` is built, ~line 334):

```r
# PHR scale annotation — median 105 kb, mean 144 kb (CROSSWALK C4 / Andrea ch. 01)
p_num_chromosomes_wide_annot <- p_num_chromosomes_wide +
  annotate("segment", x = 0, xend = 144000,
           y = -2, yend = -2,
           colour = "#1b7a3a", linewidth = 1.2) +
  annotate("text", x = 144000, y = -2,
           label = "PHR scale: median 105 kb, mean 144 kb (n=15,668)",
           hjust = -0.05, vjust = 0.5, size = 3, colour = "#1b7a3a")
ggsave("p_num_chromosomes_wide_phr_annot.pdf",
       plot = p_num_chromosomes_wide_annot,
       width = 16, height = 9)
```

## Speaker notes

This is the same view I showed last time — number of unique chromosomes per 100 kb window across CHM13, with CEN, PAR, PHR, and XTR painted as background. Last time I asked you to notice that the orange spikes pile up at the chromosome ends; this time I want to put numbers on it. Across our 466 HPRCv2 haplotypes we recovered 15,668 pseudohomologous regions — PHRs — spanning 41 of 48 chromosome arms. Median length 105 kb, mean 144 kb, range 5 to 500 kb. To anchor that scale: PAR2 is about 334 kb. So a typical subtelomeric PHR is on the same order as PAR2 — and unlike PAR2, this is happening at nearly every chromosome end. The takeaway is one sentence: **extended pseudohomology at nearly all subtelomeres, comparable in scale to canonical pseudoautosomal regions, but replicated dozens of times across the genome.** That is the central observation the rest of the talk explains — what these communities look like, who is in them, and why they exist.

## Time budget

60 seconds.

## Notes for synthesizer

- **Callbacks/setup.** Slide 04 is the genome-wide identity heatmap (avg identity per matching chromosome, 100 kb windows) with the chr18q-arm inset. This slide flips the same data from "how identical?" to "how many chromosomes share it?" — keep the visual continuity (same 100 kb window grid, same chromosome ordering on the y-axis, same Mbp x-axis). The phrase "subtelomeric patterns are known qualitatively but need more precise quantification" is set up on slide 04; slide 05 delivers the first quantification (PHR scale).
- **Forward setup.** This slide hands off to slide 06 (length distributions of inter-chromosomal matches) and slide 07 (all-vs-all heatmap). The "median 105 kb, mean 144 kb" number stated here should match the histogram summary on slide 06 — please verify the synthesis pass keeps these consistent (single source of truth: CROSSWALK §C4, citing Andrea end-to-end-report ch. 01).
- **PAR2 framing is a deliberate abstract callback.** The ABSTRACT explicitly compares PHR scale to PAR2; CROSSWALK flags this as a framing gap (REWRITE_PLAN TASK-03 — a PAR2-vs-typical-subtelomere length panel for Fig 1 is planned). The talk version of that comparison lives here as a one-liner ("PAR2 ~334 kb"). If the synthesizer ends up adding a dedicated PAR2-comparison slide, this bullet can be slimmed.
- **Figure provenance.** Primary figure is the v1 deck figure (`p_num_chromosomes_wide.pdf`), not the manuscript Fig 1b. Per task: keep v1. Fig 1b is noted only as the manuscript companion so the synthesizer doesn't accidentally swap them.
- **Annotation overlay is optional.** I provided the R snippet so the scale bar can be rendered into the figure if the synthesizer prefers a visual annotation over the bullet/speech version. Default: keep the figure clean and let the speaker say the numbers.
- **Numbers to lock down.** 15,668 PHRs / 41 of 48 arms / median 105 kb / mean 144 kb / range 5–500 kb / PAR2 ≈ 334 kb. All five appear consistently in ABSTRACT, CROSSWALK §C4, and `framing_synthesis.md` (table comparing PARs vs PHR communities, line ~50).
