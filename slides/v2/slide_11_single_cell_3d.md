## Title
Single-cell 3D — and it works in haploid sperm

## Bullets
- **GM12878 Dip-C, 16 cells (Tan 2018, remapped to T2T-CHM13v2.0):** community arms 6.9% closer in 3D within community than between (W/B = 0.931, Wilcoxon p = 3.8 × 10⁻⁴, Mantel ρ = 0.296).
- **Sperm scHi-C, 20 cells (Xu et al. 2025):** **60% closer** within community (W/B = 0.401, Fisher p = 3.9 × 10⁻⁵¹) — in the haploid, hyper-condensed sperm nucleus.
- **Negative control — pseudo-community "S" of 7 zero-sharing arms** (chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q): 11% *farther* in GM12878, 40% *farther* in sperm. Sequence sharing is necessary; community label alone is not enough.
- Same pattern across diploid soma and haploid germline → 3D clustering of subtelomeres is **not** a Hi-C population artefact and **not** restricted to interphase chromatin.
- Sperm is the bridge: it puts the signal in the germline cell that actually transmits the recombination — setting up mouse meiotic Hi-C next.

## Primary figure
**Recommended:** `paper_prep/figures/fig3/figure_fig3.pdf`, **panel (c)** — per-cell C-community W/B vs S_all (negative control), already plotted for both GM12878 (16/16 cells with W/B < 1 inside C; 0/16 inside S_all) and sperm (20/20 vs 1/20). One panel, both datasets, both directions of effect.

**Alternative if a sperm-only zoom is preferred:** `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_mantel_scatter.pdf` paired with `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_16cells_phr_mantel_scatter.pdf` side-by-side. The two-panel composition below stitches them with matching axes and a shared S_all callout — pure ggplot2/cowplot, no SBATCH:

```r
# slide_11_dipc_sperm_pair.R — side-by-side Mantel panels for GM12878 + sperm
library(ggplot2); library(magick); library(grid); library(cowplot)

gm  <- ggdraw() + draw_image(magick::image_read_pdf(
  "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_16cells_phr_mantel_scatter.pdf",
  density = 200))
spm <- ggdraw() + draw_image(magick::image_read_pdf(
  "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_mantel_scatter.pdf",
  density = 200))

label <- function(txt) ggdraw() + draw_label(txt, fontface = "bold", size = 12, hjust = 0.5)

panel <- plot_grid(
  plot_grid(label("GM12878 Dip-C (n = 16)\nW/B = 0.931 · ρ = 0.296 · p = 3.8e-04"), gm,
            ncol = 1, rel_heights = c(0.12, 1)),
  plot_grid(label("Sperm scHi-C (n = 20, Xu 2025)\nW/B = 0.401 · 60% closer · p = 3.9e-51"), spm,
            ncol = 1, rel_heights = c(0.12, 1)),
  ncol = 2
)

neg <- ggdraw() + draw_label(
  "Negative control: 7 zero-sharing arms (S_all)\nGM12878: 11% FARTHER  ·  Sperm: 40% FARTHER",
  size = 11, fontface = "italic"
) + theme(plot.background = element_rect(fill = "#FFF7E6", colour = "black", linewidth = 0.5))

plot_grid(panel, neg, ncol = 1, rel_heights = c(1, 0.12))
ggsave("slides/v2/slide_11_dipc_sperm_pair.pdf", width = 12, height = 6)
```

## Speaker notes
The bulk Hi-C signal you just saw could in principle be a population average artefact — many cells with weak or even absent contacts averaging into a coherent block. Single-cell 3D rules that out. Tan and colleagues' Dip-C reconstructs explicit 3D coordinates per allele in individual GM12878 cells; remapped to T2T-CHM13 with MAPQ-zero retention so we don't lose subtelomeric reads, 16 cells show community-member arms about 7% closer to each other than to non-members — Wilcoxon p of 3.8 × 10⁻⁴, Mantel ρ of 0.30. That's a small effect per cell, but it's there in essentially every cell. Now the punch line: we ran the same pipeline on twenty sperm cells from Xu and colleagues. Sperm is haploid, the chromatin is hyper-condensed, and the nuclear architecture is nothing like an interphase lymphoblast. The within-versus-between ratio is 0.40 — community arms are sixty percent closer in 3D — Fisher p of 3.9 × 10⁻⁵¹. Same direction, much stronger. The negative control closes the loop: a pseudo-community of the seven arms that share *no* subtelomeric sequence with anything else moves the *opposite* way — 11% farther in GM12878, 40% farther in sperm. Sequence sharing is necessary for clustering; a community label alone is not. So the 3D signal survives the bulk-to-single-cell test and survives the haploid-germline test. That last point is the bridge — sperm is the gamete that actually carries the recombination forward, and the next slide takes us into the meiotic cells where that recombination happens.

## Time budget
60 seconds.

## Notes for synthesizer
- **NEW slide vs v1:** v1 has no single-cell 3D content; this is a fresh contribution that lands two abstract claims at once — that the 3D signal is per-cell and that it generalizes to germline architecture. Both are needed for the abstract's "we hypothesize that these patterns are maintained by recombination facilitated by the physical proximity of subtelomeres" line.
- **Continuity inbound (slide 10):** previous slide is the bulk Hi-C / Pore-C 3D signal (community-level B/W < 1, Mantel rho on bulk maps). Open this slide with a single transition sentence: "Bulk Hi-C is a population average — does the signal survive at single-cell resolution?" Then dump straight into GM12878 → sperm.
- **Continuity outbound (slide 12):** next slide is mouse meiotic Hi-C (zygotene peak, lepto→zygo→pachy→diplo trajectory). The natural pivot is sperm → meiosis: "if the haploid product still shows it, what about the cells where the recombination happens?" Land the last bullet on that pivot — do **not** redo the meiosis story here.
- **Negative control is load-bearing.** The S_all pseudo-community result (11% / 40% *farther*) is the single best rebuttal to the "you're just measuring chromosome-territory crowding" objection. If the visual gets cropped, keep S_all in the bullet text at minimum.
- **Citations:** Tan et al. 2018 (Dip-C, GM12878); Xu et al. 2025 (sperm scHi-C). Both are in `paper_prep/synthesis/REFERENCES.bib`. Safe to cite verbally.
- **Provenance for numbers:** `end-to-end-report/report/06_dipc_validation.md` §"Community 3D enrichment (T2T)" and §"3D genome validation: sperm single-cell". Source TSVs live under `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/` (GM12878) and `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/` (sperm). Per-cell PDFs are listed in **Primary figure**.
- **Scope reminder:** AUDIT_REPORT / CROSSWALK currently slot Dip-C and sperm into "SI" or "do not include" for the manuscript itself. For the *talk* they are central — they are what convinces the room that the 3D claim is robust. Do not let manuscript scope decisions shrink the talk version.
- **Don't oversell Mantel ρ on sperm.** Sperm Mantel ρ = 0.202 (p = 0.023, significant but modest); the headline number is W/B = 0.401. Keep emphasis on W/B + the 60% framing.
- **Time discipline:** 60 s is tight. Single figure, three numbers (GM12878 6.9%, sperm 60%, S_all neg-control), one bridge sentence to slide 12. If forced to cut, drop the Wilcoxon-vs-Mantel detail; never drop the negative control.
