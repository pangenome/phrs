## Title
Genome-wide identity heatmap — interchromosomal homology at PAR2 scale

## Bullets
- 466 HPRCv2 haplotypes vs CHM13, per-position max identity to any matching chromosome (100 kb windows).
- Most of each chromosome is silent; dense red bands appear where assemblies reach the telomeres.
- chr18 q-arm inset: tight subtelomeric block of >98% identity matches from many chromosomes.
- These blocks span **10s–100s of kb — comparable in scale to PAR2 (~334 kb on Xq/Yq)**.
- PAR2-scale pseudo-homology at nearly every subtelomere → motivates the all-vs-all view that follows.

## Primary figure
**Recommended:** `paper_prep/figures/fig1/figure_fig1.pdf`, **panel (a)** — genome-wide stacked identity heatmap with the chr18 q-arm subtelomeric inset already embedded. This is the polished version of the v1 slide-4 visual concept (1:1 substitute, no rebuild needed).

**Alternative if a single-chromosome zoom is preferred:** `identity_heatmap_chr18.zoom_last1mb.pdf` (the chr18q subtelomeric panel alone, larger and more legible than the inset).

**PAR2 callout overlay (new; for synthesizer):** add a small callout on top of panel (a) — preferred placement is over the chrX/chrY rows of the main heatmap (or in white space at the right edge of the inset). Suggested text:

> Extended interchrom homology spans 10s–100s kb at nearly all subtelomeres — *comparable in scale to PAR2* (~334 kb, chrXq/chrYq).

If the synthesizer wants to render the overlay programmatically rather than annotate in Keynote/PowerPoint, the following ggplot2 snippet wraps any imported panel-1a PNG with the callout. It needs no data and no SBATCH:

```r
# slide_04_par2_callout.R — overlay a PAR2-scale callout on panel 1a
library(ggplot2); library(png); library(grid); library(cowplot)

panel_a <- readPNG("paper_prep/figures/fig1/figure_fig1.png")  # or a cropped 1a export
bg      <- ggdraw() + draw_image(panel_a)

callout <- ggdraw() +
  draw_label(
    "Extended interchrom homology: 10s–100s kb\n(comparable to PAR2, ~334 kb on Xq/Yq)",
    x = 0.5, y = 0.5, hjust = 0.5, size = 11, fontface = "bold"
  ) +
  theme(plot.background = element_rect(fill = "#FFF7E6", colour = "black", linewidth = 0.6))

# place callout in upper-right of the panel; tune x/y/width/height to taste
ggdraw(bg) +
  draw_plot(callout, x = 0.62, y = 0.82, width = 0.36, height = 0.10)

ggsave("slides/v2/slide_04_panel_with_par2_callout.pdf", width = 12, height = 7)
```

## Speaker notes
This is the empirical foundation of the talk. For every 100 kb window across 466 HPRCv2 haplotypes we plot the maximum alignment identity to any *other* chromosome. Most of each chromosome is silent — alignments only hit the same chromosome, as expected. But at the telomeres, dense red bands appear wherever assemblies reach the chromosome end. These are the inter-chromosomal exchange blocks. The chr18 q-arm inset zooms in on one: a tight band at the very end where many other chromosomes match at over 98% identity, extending tens to hundreds of kilobases inward from the telomere. Here is the reframe — this scale is what PAR2 looks like on the sex chromosomes. PAR2 is about 334 kb. We are seeing PAR2-scale pseudo-homology at nearly every subtelomere. That dramatically expands the known scope of pseudohomologous regions in the human genome, which is the thesis the next slides will quantify.

## Time budget
70 seconds.

## Notes for synthesizer
- **Concept preserved from v1 slide 4** (genome-wide 100kb identity heatmap + chr18q inset). Visual is a clean 1:1 substitute with the publication panel.
- **One new asset needed:** the PAR2 callout text. It is the *only* new content vs v1; it lands the abstract's central reframing ("comparable in scale to PAR2"). If you skip the overlay, at minimum say it on screen as a subtitle bullet.
- **Continuity inbound (slide 03):** the previous slide ends on IMPG / all-vs-all alignment as the *method*. This slide is the first *result* — open with "what does that look like?" The dense bands at telomeres are the answer.
- **Continuity outbound (slide 05):** sets up the unique-chromosomes-per-region view (v1 slide 5 / Fig 1b). The natural pivot at the end is "*how many* chromosomes are mixing here?" — which is what slide 05 quantifies.
- **PAR2 number provenance:** 334 kb is the canonical PAR2 length on chrXq/chrYq (T2T-CHM13). It is referenced repeatedly in `paper_prep/synthesis/AUDIT_REPORT.md` (C4) and `paper_prep/synthesis/CROSSWALK.md` (C4 entry, "PAR2 anchor 334 kb"). Safe to cite verbally.
- **Do not relitigate the inset.** v1 used chr18q, abstract / Fig 1a both keep chr18q — keep chr18q.
