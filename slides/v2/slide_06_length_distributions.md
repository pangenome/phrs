## Title

Length distributions of inter-chromosomal matches per arm — the outliers are named clades

## Bullets

- Faceted histograms (kept verbatim from v1 slide 6): one panel per chromosome arm; **p-arm = blue, q-arm = orange**; pink fill marks the five introvert arms with no inter-chromosomal hits (2p, 3p, 5p, 8q, 11q, 14q).
- The bulk of arms cluster around the population scale stated on slide 05 — **median 105 kb, mean 144 kb, range 5–500 kb** (CROSSWALK §C4; n=15,668 PHR sequences across 41/48 arms).
- **The fat-right-tail outliers are not noise — they are the abstract's clades.** Acrocentric short arms (13p, 14p, 15p, 21p, 22p = **C7**) carry the heaviest, most uniform long-length distributions: fully homogenized rDNA-adjacent arms. PAR2 (Xq, Yq = **C14**) and PAR1 (Xp, Yp = **C15**) sit at the canonical pseudoautosomal scale (~334 kb on Xq/Yq). 4q and 10q (= **C1**, D4Z4 / DUX4) and the 10p–18p Linardopoulou pair (= **C2**) round out the six abstract-anchored clades.
- Pink panels are a biological signal, not missing data: arms without enough cross-chromosomal sharing to enter the 41×41 matrix.
- One-sentence reframe: **the shape of the per-arm length distribution recapitulates the community partition you are about to see on slide 09**.

## Primary figure

**Reuse v1 verbatim:** `/moosefs/guarracino/HPRCv2/PHR_III/plots/all-vs-all.1Mb.p95.id95.len_length_dist_by_chr_arm.pdf` (the exact figure on v1 slide 6 — faceted histograms by chromosome arm, p blue / q orange, pink for missing introvert arms). The underlying data is `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv` (18,827 sequences; columns `seq, arm, self_chr, region_start, region_end, chrs_involved, arms_involved`; PHR length per row = `region_end - region_start`).

**One new asset (this slide's deliverable):** a clade-callout overlay that names the four outlier groups directly on the facet grid. Renders to PDF / PNG locally with vanilla ggplot2 — the synthesizer can either drop the rendered overlay box onto the v1 figure in Keynote, or rebuild the histogram panel and let the overlay annotate live (data file is small enough to read directly with `readr::read_tsv`).

```r
# slide_06_clade_callouts.R
# Renders a callout legend that names the outlier facets on the v1 length-
# distribution facet grid. Pure ggplot2; no SBATCH; runs in seconds.
library(ggplot2)

cal <- data.frame(
  row   = 5:1,                                        # C1 lands at top of table
  C     = c("C7", "C14", "C15", "C1", "C2"),
  arms  = c("13p, 14p, 15p, 21p, 22p",
            "Xq, Yq",
            "Xp, Yp  (+18q, n=1)",
            "4q, 10q",
            "10p, 18p"),
  clade = c("acrocentric short arms — fully homogenized (rDNA-adjacent)",
            "PAR2 — pseudoautosomal q-end (~334 kb scale)",
            "PAR1 — pseudoautosomal p-end",
            "4q–10q DUX4 / D4Z4 — long-tail, copy-number diverse",
            "10p–18p — Linardopoulou 2005 pair"),
  fill  = c("#E5D8FA",   # C7 acrocentric (matches slide_09 palette)
            "#CDEAD3",   # C14 PAR2
            "#CDEAD3",   # C15 PAR1
            "#FDE2C8",   # C1 DUX4
            "#FFF1AA"),  # C2 Linardopoulou
  stringsAsFactors = FALSE
)

cols <- data.frame(
  x     = c(0.05, 0.22, 0.55),
  width = c(0.14, 0.30, 1.20),
  field = c("C", "arms", "clade"),
  stringsAsFactors = FALSE
)

cells <- do.call(rbind, lapply(seq_len(nrow(cols)), function(j) data.frame(
  x        = cols$x[j] + cols$width[j] / 2,
  width    = cols$width[j],
  y        = cal$row,
  text     = cal[[ cols$field[j] ]],
  fill     = cal$fill,
  fontface = ifelse(j == 1, "bold", "plain"),
  stringsAsFactors = FALSE
)))

hdr <- data.frame(
  x        = cols$x + cols$width / 2,
  width    = cols$width,
  y        = 6,
  text     = c("C", "Arms (outlier facets)", "Named clade — why the tail is long"),
  fill     = "#EEEEEE",
  fontface = "bold",
  stringsAsFactors = FALSE
)

title <- data.frame(
  x        = 0.85,
  width    = 1.6,
  y        = 7,
  text     = "The fat-right-tail facets are the abstract's clades",
  fill     = "#FFFFFF",
  fontface = "bold",
  stringsAsFactors = FALSE
)

note <- data.frame(
  x        = 0.85,
  width    = 1.6,
  y        = 0,
  text     = "Pink facets (2p, 3p, 5p, 8q, 11q, 14q) = introvert arms, no cross-chrom hits.",
  fill     = "#FFFFFF",
  fontface = "italic",
  stringsAsFactors = FALSE
)

p <- ggplot() +
  geom_tile(data = rbind(cells, hdr),
            aes(x = x, y = y, width = width, height = 0.95, fill = I(fill)),
            colour = "grey75") +
  geom_text(data = cells, aes(x = x - width/2 + 0.01, y = y, label = text,
                              fontface = fontface),
            hjust = 0, size = 3.0) +
  geom_text(data = hdr, aes(x = x - width/2 + 0.01, y = y, label = text,
                            fontface = fontface),
            hjust = 0, size = 3.2) +
  geom_text(data = title, aes(x = x - width/2 + 0.01, y = y, label = text,
                              fontface = fontface),
            hjust = 0, size = 3.6) +
  geom_text(data = note, aes(x = x - width/2 + 0.01, y = y, label = text,
                             fontface = fontface),
            hjust = 0, size = 2.6, colour = "grey35") +
  coord_cartesian(xlim = c(0, 1.7), ylim = c(-0.5, 7.5), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(6, 6, 6, 6))

ggsave("slide_06_clade_callouts.pdf", p, width = 8.0, height = 3.2)
ggsave("slide_06_clade_callouts.png", p, width = 8.0, height = 3.2, dpi = 200)
```

The five colored rows are exactly the abstract-anchored outlier clades (C1, C2, C7, C14, C15); palette matches slide_09 so the eye carries the same color-coding from histograms → community PCA. Place the rendered callout in the upper-right white space of the v1 facet grid (or beneath it as a banner) — this is the only on-slide change vs v1.

## Speaker notes

This is the same panel I showed last time — one histogram per chromosome arm, p-arms in blue, q-arms in orange, pink for the introvert arms with no cross-chromosomal hits at all. Most arms cluster around the population scale you saw on the previous slide: median around 105 kilobases, mean 144, range 5 to 500. What I want you to notice now is the **shape** of the outlier facets — the ones with the heaviest right tails, the most uniformly long distributions. Those are not noise, those are the clades the abstract names. The five acrocentric short arms — 13p, 14p, 15p, 21p, 22p — community 7 — fully homogenized, rDNA-adjacent. Xq with Yq — community 14 — PAR2, around 334 kb, the canonical pseudoautosomal scale. Xp with Yp — community 15 — PAR1. 4q with 10q — community 1 — the DUX4 / D4Z4 pair, long-tail and copy-number diverse. 10p with 18p — community 2 — Linardopoulou's 2005 pair. The pink facets — 2p, 3p, 5p, 8q, 11q, 14q — that absence is also a signal: those arms simply do not participate in the inter-chromosomal exchange landscape. So a single sentence carries this slide: **the shape of these per-arm distributions already recapitulates the community partition I will show you on slide 9.** The histograms know the clades before we even cluster.

## Time budget

50 seconds. Roughly: 10 s recap of the v1 visual ("same panel as last time, p blue / q orange, pink = introvert arms, scale anchored to the previous slide"); 30 s naming the five outlier clades on the callout (C7 acrocentric → C14 PAR2 → C15 PAR1 → C1 D4Z4 → C2 Linardopoulou, ~6 s each); 10 s landing the reframe ("the shape already recapitulates the community partition you'll see on slide 9").

## Notes for synthesizer

- **Concept preserved from v1 slide 6.** Faceted histograms by chromosome arm (p blue / q orange / pink for missing introvert arms). Do not rebuild the panel — the v1 PDF at `/moosefs/guarracino/HPRCv2/PHR_III/plots/all-vs-all.1Mb.p95.id95.len_length_dist_by_chr_arm.pdf` is what goes on the slide.
- **One new asset:** the clade-callout legend (R script above). It is the only addition vs v1. The script needs no data and no SBATCH; the rendered PDF / PNG sits over white space in the v1 facet grid (upper-right) or below it as a banner.
- **Continuity inbound (slide 05).** Slide 05 lands the population PHR scale (median 105 kb, mean 144 kb, range 5–500 kb, comparable to PAR2 ~334 kb) and the 15,668-PHRs-across-41-arms count. Slide 06 takes that same scale and shows it varies by arm in a *biologically structured* way — the outlier arms are the named clades. Lock these numbers to a single source of truth: CROSSWALK §C4, citing Andrea ch. 01.
- **Continuity outbound (slide 07 / 08 / 09).** This slide's reframe ("the shape already encodes the community partition") sets up the all-vs-all heatmap (slide 07) and the keystone PCA-by-community slide (slide 09). The five clade names spoken here (C1, C2, C7, C14, C15) MUST be the same five spoken on slide 09 — same colors, same arm lists. Synthesizer: keep the slide_06 and slide_09 callout palettes consistent (the R scripts already share the C1/C2/C7/C14/C15 fill colors).
- **Abstract-clade ↔ outlier-facet map (load-bearing — copy verbatim if compressing):**
  - C1 (4q, 10q) — DUX4 / D4Z4 long-tail; chapter 03 (28 DUX4L genes), chapter 04 (43.4 % type-discordance).
  - C2 (10p, 18p) — Linardopoulou 2005 Fig. 5 pair.
  - C7 (13p, 14p, 15p, 21p, 22p) — acrocentric short arms, fully homogenized (rDNA-adjacent). Heaviest, most uniform long-length distributions on the panel.
  - C14 (Xq, Yq) — PAR2, ~334 kb canonical scale.
  - C15 (Xp, Yp; 18q n=1) — PAR1.
  - Source: `paper_prep/synthesis/CROSSWALK.md` §C5; `end-to-end-report/report/01_pipeline.md` §"Arm-level community detection"; v1 slide 10 community list.
- **Missing-introvert callout (2p, 3p, 5p, 8q, 11q, 14q).** Same six arms appear on v1 slides 7, 8, 9, 10 in the upper-right corner. Keep the v1 phrasing here too — "introvert arms" is Erik's term and is consistent with "no inter-chromosomal hits" (their PHR length distribution is empty, hence the pink fill in v1).
- **Numbers to lock down across slides 05 / 06.** 15,668 PHRs / 41 of 48 arms / median 105 kb / mean 144 kb / range 5–500 kb / PAR2 ≈ 334 kb. All sourced from CROSSWALK §C4 (citing Andrea ch. 01) and `paper_prep/synthesis/ABSTRACT.md`. If a synthesis pass adjusts any of these, update slide 05 and slide 06 together.
- **Why the outliers are outliers (one-line each, for Q&A).** C7: rDNA recombination homogenizes the entire short arm; C14/C15: pseudoautosomal obligate crossover homogenizes both arms across X/Y; C1: D4Z4 / DUX4 macrosatellite copy-number polymorphism makes the right tail extend; C2: 10p–18p is the historical Linardopoulou exchange pair (long-known sequence-level homology block).
- **Do NOT touch the panel rendering.** The v1 PDF is on `/moosefs/guarracino/...` — out-of-worktree, but accessible read-only. Per task: keep v1. The slide deliverable is the *callout* on top of v1, not a rebuild.
- **If 50 s is cut to 35 s** (deck pacing): drop C2 (Linardopoulou) from the spoken walk — keep C7 acrocentric, C14 PAR2, C15 PAR1, C1 D4Z4. C2 is in the abstract clade list but is the smallest-effect outlier on the histograms and is recoverable on slide 09 if needed.
- **Cross-slide concern.** The v1 figure has small per-facet `n=` annotations (e.g., n=428, n=446); per the v1 summary these are sample counts per arm. Speaker should NOT verbalize per-facet counts at 50 s — they are visible on the figure for anyone who wants to read them. Save the per-arm n discussion for Q&A.
