## Title

Concerted evolution of human subtelomeres — what we saw, predicted, and recovered

## Bullets

- **Method (slide 3).** Implicit pangenome graph: wfmash all-vs-all over 18,827 telomere-anchored flanks, ~12% of all pairs evaluated — **230× above the Erdős-Rényi connectivity threshold** (`p* = log(n)/n ≈ 5.21×10⁻⁴`). No chromosome partitioning, no GFA: every haplotype is its own reference.
- **Empirical (slides 4–9).** **15,668 PHRs across 41/48 arms**, median 105 kb / mean 144 kb — PAR2-scale pseudohomology at nearly every chromosome end. Named clades: Xp/Yp & Xq/Yq via PARs, **acrocentric short arms** (C7, near-interchangeable), **10p–18p**, the big q-arm clade (**22q–21q–19q–1q–13q–17q**), and the **4q–10q DUX4** clade.
- **Mechanism (slides 10–12).** Hi-C/Pore-C/CiFi/Dip-C all independently recover community-structured 3D contacts (B/W 0.027–0.074; Mantel ρ=0.296, p=0.002; per-pair ρ=0.674 in CHM13 Hi-C). Median PHR (105 kb) sits at the base of a single meiotic loop — the **bouquet** is the predicted exchange venue.
- **Proof (slide 13).** **WashU T2T pedigree: 494/538 (92%) inter-chromosomal patches fall inside Leiden communities** — 133 gene-conversion-like, 16 crossover-like. The graph-derived community structure *predicts* where exchanges show up in real families.
- **Biology (slide 14).** D4Z4 / DUX4 (4q↔10q) is the disease-revealed instance of the same process — FSHD-modifying translocations are concerted evolution caught in the act.

> **Thesis: subtelomeres concertedly evolve through ongoing inter-chromosomal exchange — observable in pedigrees, predicted by 3D, recovered by an implicit pangenome graph across HPRC v2.**

## Primary figure

**Recommended:** `paper_prep/figures/ed8/figure_ed8.pdf` panel **(a)** — the four-link causal feedback loop (sequence sharing → 3D proximity → ectopic exchange → new shared segments → propagation), with edges color-coded by evidence type (solid blue = direct measurement, solid green = literature, dashed olive = inferred). This is the closer-friendly schematic; Andrea built it as the synthesis figure of the discussion (`end-to-end-report/report/07_integrated.md §1.6` → `paper_prep/figures/ed8/`).

**Why ed8(a) and not a new figure:** the closer should *recapitulate*, not introduce. Ed8(a) already encodes every pillar of the talk — sequence (slide 4–9), 3D (10–12), exchange (13–14) — as a single loop. No new compute, no SBATCH, no R needed.

If the synthesizer wants a thin annotation overlay that names each edge with its slide number (so the audience can map the loop back to what they saw), the following snippet wraps the existing PNG with four labels — runs locally in seconds, no data dependency:

```r
# slide_15_loop_callouts.R — overlay slide-number tags on the ed8(a) feedback loop
# Output: slide_15_loop_with_callouts.pdf  (drop-in for the slide background)
library(ggplot2); library(png); library(grid); library(cowplot)

panel_a <- readPNG("paper_prep/figures/ed8/figure_ed8.png")
bg      <- ggdraw() + draw_image(panel_a)

tag <- function(text, x, y) {
  ggdraw() + draw_label(text, x = 0.5, y = 0.5, hjust = 0.5, size = 9,
                        fontface = "bold", colour = "#1b3a6f") +
    theme(plot.background = element_rect(fill = "#FFF7E6",
                                         colour = "#1b3a6f", linewidth = 0.4))
}

# Four edge tags — coordinates are nominal; tune to the actual ed8(a) layout
ggdraw(bg) +
  draw_plot(tag("sequence sharing\n(slides 4–9)"),    x = 0.04, y = 0.78, width = 0.20, height = 0.07) +
  draw_plot(tag("3D proximity\n(slides 10–12)"),      x = 0.78, y = 0.78, width = 0.20, height = 0.07) +
  draw_plot(tag("ectopic exchange\n(slide 13 pedigree)"), x = 0.78, y = 0.10, width = 0.22, height = 0.07) +
  draw_plot(tag("propagation\n(slide 14 DUX4)"),      x = 0.04, y = 0.10, width = 0.20, height = 0.07)

ggsave("slides/v2/slide_15_loop_with_callouts.pdf", width = 12, height = 7)
```

(Synthesizer: ed8(a) is the recommended *primary* — the snippet above is optional polish, not a blocker. If the synthesizer skips the overlay, the bullets already do the slide-number tagging in words.)

## Speaker notes

This is the close. One slide, one breath, one thesis.

We started with a methodological commitment: the implicit pangenome graph. Wfmash over 18,827 telomere-anchored flanks, twelve percent of all pairs evaluated — and that twelve percent is two hundred and thirty times above the Erdős-Rényi threshold for graph connectivity. That is what licenses everything that followed: no chromosomal partitioning, every haplotype its own reference.

What that method showed us is that pseudohomology at PAR2 scale is replicated at nearly every chromosome end — fifteen thousand six hundred and sixty-eight PHRs across forty-one of forty-eight arms, median one hundred and five kilobases. We named the clades: PARs at Xp/Yp and Xq/Yq, the acrocentric short arms, ten-p with eighteen-p, the big q-arm clade — twenty-two q, twenty-one q, nineteen q, one q, thirteen q, seventeen q — and four q with ten q carrying DUX4.

We then asked: why are these in clades? Hi-C, Pore-C, CiFi, and single-cell Dip-C all independently recover community-structured three-dimensional contacts, with per-pair correlations as high as point-six-seven-four in CHM13. The median PHR fits at the base of a single meiotic loop — the bouquet is the predicted exchange venue.

And then the proof. The WashU T2T pedigree gives us four-hundred-ninety-four out of five-hundred-thirty-eight inter-chromosomal patches — ninety-two percent — falling *inside* the Leiden communities the graph predicted. That is the loop closing: sequence similarity, 3D proximity, and observed family-level exchange all agreeing.

D4Z4 and DUX4 are the disease-revealed instance — FSHD translocations are concerted evolution caught in the act on a clinically named locus.

The thesis: **subtelomeres concertedly evolve through ongoing inter-chromosomal exchange — observable in pedigrees, predicted by 3D, recovered by an implicit pangenome graph across HPRC v2.** Thank you.

## Time budget

**Target: 70 seconds.** Roughly 10 s on the method recap (Erdős-Rényi number stays in), 15 s on the empirical pillar (PHR scale + named clades — name 3–4 of them, don't list all), 15 s on the mechanism (3D convergence + bouquet, one Mantel number), 15 s on the pedigree proof (the 92 % number is the headline, do not skip it), 5 s on DUX4 as the disease-revealed instance, 10 s on the one-line thesis read aloud verbatim. The thesis sentence is the last thing the audience hears — do not paraphrase it on the fly.

## Notes for synthesizer

- **What slide 14 sets up.** Slide 14 lands DUX4/FSHD as the case study in concerted evolution caught in the act on a disease locus. Slide 15 must *not* re-explain DUX4; it cites it as the mechanistic exemplar in one bullet and moves on. The transition from 14 → 15 is "and that's not just one locus — here's the whole picture."
- **Callback discipline.** Every bullet on this slide is tagged to a prior slide number. **Do not drop the tags** in compression. The whole point of a closer is recapitulation; the audience needs the visual cue ("slide 13," "slides 10–12") to remember they already saw the evidence and trust the synthesis.
- **The thesis sentence is locked from the task spec.** "Subtelomeres concertedly evolve through ongoing inter-chromosomal exchange — observable in pedigrees, predicted by 3D, recovered by an implicit pangenome graph across HPRC v2." That exact wording is what the task asked for and what should appear on the slide. Synthesizer: this is a hard constraint, not a suggestion. Render as a blockquote / pull-quote, large type, bottom of slide.
- **Title-callback.** The manuscript title (slide 01) is "Concerted evolution and unorthodox recombination of human subtelomeres" (`paper_prep/synthesis/ABSTRACT.md`). Slide 15's headline echoes the "concerted evolution" half deliberately so the deck closes by closing the title — first slide and last slide bracket the same phrase.
- **Numbers to lock down (single source of truth: CROSSWALK + ABSTRACT).**
  - 18,827 telomere-anchored flanks; ~12% sampling; ER `p* ≈ 5.21×10⁻⁴`; ~230× threshold (slide 3)
  - 15,668 PHRs / 41 of 48 arms / median 105 kb / mean 144 kb (slide 5; matches `framing_synthesis.md`)
  - PAR2 ≈ 334 kb (slide 4 anchor)
  - Hi-C B/W 0.027–0.074, p 6.0e-18 to 9.4e-03; Mantel ρ=0.296, p=0.002; CHM13 per-pair ρ=0.674 (slides 10–12; from `07_integrated.md`)
  - WashU pedigree 494/538 = 92% within-community; 133 gene_conversion_like, 16 crossover_like (slide 13; from `14_pedigree_recombination.md`)
  - C1 = chr4_q/chr10_q, median 22 DUX4L, Mann-Whitney p=5.3e-6 (slide 14; from `07_integrated.md` D4Z4-CTCF section)
  - **All six numeric clusters above must agree with the slides they call back to.** If the synthesizer notices a drift between, e.g., slide 13's pedigree number and the number cited here, the slide-13 author is the source of truth — file a `wg msg` rather than silently editing this slide.
- **Figure choice rationale.** Ed8(a) is preferred because (i) it already exists, (ii) Andrea built it specifically as the discussion synthesis figure (`paper_prep/figures/ed8/caption.md`), (iii) it encodes the talk's argument as a single loop. The optional R overlay in *Primary figure* is a "nice to have" — if the synthesizer is short on time, ship ed8(a) bare.
- **Forward setup.** This is the last slide before Q&A. Do **not** add a "future work" / "thanks" slot inside this slide — the task spec is explicit (15 slides total, this is slide 15). Acknowledgements, if any, belong in a separate appendix slide that is not part of the 15-slide budget.
- **If compression is forced** (deck overruns and the synthesizer needs to cut): the order to drop is (1) the optional R overlay, (2) the named-clades enumeration in bullet 2 (drop 10p–18p and the big q-arm clade list, keep PAR + acro + 4q/10q), (3) the per-technology B/W numbers in bullet 3 (keep Mantel ρ=0.296 as the headline). **Do not drop:** the 230× ER number, the 92 % pedigree number, or the thesis pull-quote. Those three are the slide.
- **No SBATCH needed.** Ed8(a) PDF/PNG already exists at `paper_prep/figures/ed8/figure_ed8.{pdf,png}`. The optional overlay R script runs in seconds in the agent worktree and writes a single PDF.
- **No edits to other v2 slides.** Task spec: single output file, single commit. Cross-slide concerns are reported here for the synthesizer; this slide does not modify any neighbor.
