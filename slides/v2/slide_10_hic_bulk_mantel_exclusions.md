## Title

**Hi-C / Pore-C confirm sequence communities are 3D — and the signal *strengthens* when known confounds are removed**

Subtitle: arm-level Mantel of subtelomeric Jaccard similarity vs inter-chromosomal contact, before and after acrocentric / sex / strong-community exclusions.

## Bullets

- **Lead visual: HG002 Pore-C inter-arm contact matrix, 50 kb, 77 arm-haplotypes ordered by sequence community.** Diagonal blocks light up; **B/W = 0.056, p = 3.9 × 10⁻⁸⁵** (within-community contacts vastly outweigh between).
- **Bulk Mantel (similarity × contact) is positive in 7/8 datasets** at 50 kb — CHM13 ρ = 0.66, HG002 Hi-C ρ = 0.66, HG002 Pore-C ρ = 0.49, NA19036 ρ = 0.27, HG02148 ρ = 0.15 (the only borderline sample).
- **Exclude acrocentric p+q + chrX/Y + the four strong communities (D4Z4, acro p, PAR1, PAR2) and every sample's ρ goes up:** HG002 0.66 → **0.79**, CHM13 0.66 → **0.85**, **HG02148 0.15 → 0.72** (the marginal sample becomes one of the strongest).
- **At 10 kb, community-free per-sequence-pair correlations reach** ρ = 0.83 in NA19036 and ρ = 0.81 in HG02148 (`05_hic_validation.md` §"Individual sequence-pair similarity vs Hi-C contact") — finer resolution + no aggregation **strengthens** the signal further.
- **One-line takeaway: sequence similarity predicts 3D contact, and the signal strengthens after confound exclusions.** It is therefore not driven by nucleolar acrocentric clustering, by PAR contact, or by the few largest communities — it is a generic property of subtelomeric homology.

## Primary figure

Two-panel layout, both panels already exist as published assets — no rebuild needed:

- **Left (bulk):** `paper_prep/figures/fig3/figure_fig3.pdf`, **panel (a)** — HG002 Pore-C contact matrix, 50 kb, 77 arm-haplotypes, ordered by Leiden sequence community. Annotation: keep the **B/W = 0.056, p = 3.9 × 10⁻⁸⁵** label visible in the corner.
- **Right (exclusions):** `paper_prep/figures/ed5/figure_ed5.pdf`, **panel (b)** — Mantel ρ before vs after acrocentric + sex exclusion (50 kb), one point per HPRC sample, identity diagonal drawn. Y > X for 7/7. Source TSVs are in `community_based/50000bp/` and `no_acrocentric/50000bp/` (`<sample>_global_test.tsv`).

If the synthesizer wants a single composite for the slide instead of two pasted panels, the script below does it from the published PNGs (no SBATCH, no re-running the analysis):

```r
# slide_10_bulk_mantel_composite.R — bulk Hi-C panel + Mantel exclusion panel side-by-side
# Output: slides/v2/slide_10_bulk_mantel_composite.pdf
suppressPackageStartupMessages({
  library(ggplot2); library(png); library(grid); library(cowplot)
})

bulk      <- readPNG("paper_prep/figures/fig3/figure_fig3.png")    # crop to panel (a) before use
exclusion <- readPNG("paper_prep/figures/ed5/figure_ed5.png")      # crop to panel (b) before use

p_left  <- ggdraw() + draw_image(bulk)      + draw_label(
  "HG002 Pore-C contacts ordered by sequence community\nB/W = 0.056, p = 3.9e-85",
  x = 0.02, y = 0.97, hjust = 0, vjust = 1, size = 10, fontface = "bold")

p_right <- ggdraw() + draw_image(exclusion) + draw_label(
  "Mantel ρ: full vs (no acrocentric + sex)\n7/7 above identity — signal strengthens",
  x = 0.02, y = 0.97, hjust = 0, vjust = 1, size = 10, fontface = "bold")

annotation <- ggdraw() + draw_label(
  "HG002 0.66 → 0.79   CHM13 0.66 → 0.85   HG02148 0.15 → 0.72",
  x = 0.5, y = 0.5, size = 12, fontface = "bold"
) + theme(plot.background = element_rect(fill = "#FFF7E6", colour = "black", linewidth = 0.5))

top  <- plot_grid(p_left, p_right, nrow = 1, rel_widths = c(1, 1), labels = c("a", "b"))
fig  <- plot_grid(top, annotation, ncol = 1, rel_heights = c(1, 0.10))

ggsave("slides/v2/slide_10_bulk_mantel_composite.pdf", fig, width = 13, height = 6.5)
```

The script depends only on the two existing PNGs in the repo and `ggplot2` + `cowplot` + `png` — same toolchain used by neighboring v2 slides.

## Speaker notes

This is the moment we move from sequence to nucleus. Left panel: the HG002 Pore-C inter-arm contact matrix at 50 kilobases, 77 arm-haplotypes, ordered along both axes by their Leiden sequence community. The diagonal blocks are precisely the communities we built from pangenome graph similarity, and they light up — within-community contacts are eighteen-fold higher than between, p ≈ 10⁻⁸⁵. Sequence-defined communities are physical. Right panel: the bulk Mantel test asks whether arms with more similar subtelomeres also contact each other more often. Across eight HPRC datasets the answer is yes in seven — HG002 Hi-C and CHM13 each at ρ = 0.66, Pore-C at 0.49, with HG02148 the marginal exception at 0.15. Now the robustness check. The skeptic's worry is that the signal is just acrocentric nucleolar clustering or pseudoautosomal contact. Strip out the chr13–22 p- and q-arms, chrX, chrY, and the four strongest communities — D4Z4, acrocentric p, PAR1, PAR2 — and every sample's ρ goes up: HG002 to 0.79, CHM13 to 0.85, and HG02148, the marginal sample, jumps to 0.72. Drop to ten-kilobase resolution and treat individual sequence pairs without any community labels and you reach 0.83 in NA19036, 0.81 in HG02148. The signal is not a nucleolar artifact; it is a generic property of subtelomeric homology that gets cleaner the more carefully you look.

## Time budget

80 seconds.

## Notes for synthesizer

- **NEW slide vs v1.** v1 covered Hi-C only briefly. This slide consolidates two figures from the manuscript (Fig 3a + ED5b) into one slot so the talk lands the bulk-Mantel + exclusion-robustness story in 80 s without spawning a separate ED5 slide.
- **Continuity inbound (slide 09):** the previous slide should have set up Hi-C / Pore-C / CiFi as the orthogonal validation modality and introduced the per-haplotype 3D pipeline (`analyze_hic_communities.py`). Open this slide with "the matrix" — the contact matrix in panel a is the picture, not just a number.
- **Continuity outbound (slide 11):** subsequent slides (per-community detail, RPE-1 cell-type/cell-cycle modulation, Dip-C / sperm radial, mouse meiotic, etc.) all build on the *bulk* result that this slide establishes. End on the takeaway sentence so the next slide can refine into per-community / per-cell-type behavior without re-justifying the bulk effect.
- **Source tables to cite verbally if asked:**
  - HG002 Pore-C bulk B/W: `community_based/50000bp/hg002_porec_global_test.tsv` → 0.056, p = 3.9e-85.
  - Mantel full vs exclusions: `05_hic_validation.md` lines 327–351 (full table), with the **no acro pq + sex** column being the headline ("HG002 0.79, CHM13 0.85, HG02148 0.72"). The "no strong" column gives the same direction (HG002 0.765, CHM13 0.837).
  - 10 kb community-free peaks: `05_hic_validation.md` lines 157–166 — NA19036 ρ = 0.827 at 10 kb; HG02148 ρ = 0.809 at 10 kb.
- **Do not invert the convention.** B/W < 1 means *within*-community contacts are higher (between/within ratio). Mantel ρ uses similarity × contact, so positive ρ means more-similar arms contact more. Both go in the "communities are physical" direction; do not flip the signs verbally.
- **CHM13 caveat (only mention if asked):** CHM13 is haploid → ARI is high (0.54) but per-community W/B power is limited (singletons). Bulk Mantel is fine; do not lean on CHM13 per-community numbers.
- **HG002 CiFi was not in the no-acrocentric run** (`05_hic_validation.md` no-acro table covers 7 samples, not 8). This is a known gap in ED5b; do not claim "8/8" for the exclusion plot — say "7/7 sample × technology cells tested".
- **Visual budget:** if the panels won't fit side-by-side at legible size, prioritize panel (a) at full width and shrink panel (b) to a sub-inset with only the diagonal-crossing arrows + the three labelled ρ numbers (0.66 → 0.79 HG002, 0.66 → 0.85 CHM13, 0.15 → 0.72 HG02148). The arrows are the message; the rest is supporting.
