## Title

Gene biology aside — DUX4, OR4F, TAR1 (the biology is interesting too)

## Bullets

- **DUX4 (FSHD locus).** Annotated on 18 q-arms across the pangenome, but only chr4q and chr10q (community **C1**) carry the full **D4Z4 macrosatellite** array (median **22** DUX4L copies; the FSHD-permissive 4qA haplotype lives here). All other arms: 0–2 copies — DUX4 has scattered, but the disease-relevant repeat unit has not.
- **OR4F (olfactory receptors).** 4 OR4F paralogs span 16 arms, **5,023** gene-copy entries total. Pseudogenisation runs as a clean per-arm gradient: **11.1% pseudogene at chr7p → 99.8% at chr15q** (population mean 62.1%). Same gene, same neighborhood — but the decay clock has run longer at one end.
- **TAR1 (telomere-associated repeat).** 21,544 entries across **94.6%** of all 15,668 PHR sequences and all 41 arms — universal except at PAR1 (chrX_p / chrY_p, **0.5%**). PAR1 has obligate meiotic crossover; satellite-mediated exchange anchors are evidently not required there.
- One sentence: **distinct biological histories — disease, decay, and exchange machinery — write themselves into the same subtelomeric architecture.**

## Primary figure

`slides/v2/slide_14_gene_biology.pdf` (this worktree) — produced by `slides/v2/slide_14_gene_biology.R`. Three panels, 16 × 5.2 in:

- **Panel a (DUX4):** boxplot of DUX4L copies per haplotype on chr4q and chr10q (the only arms with full D4Z4 arrays); red dashed line at C1 median = 22. Annotation reminds the audience that DUX4 is annotated on 18 q-arms but only C1 carries the macrosatellite.
- **Panel b (OR4F):** per-arm pseudogene fraction sorted ascending; chr7p (11.1%) and chr15q (99.8%) extremes labelled in red; population mean 62.1% as red dashed reference.
- **Panel c (TAR1):** per-arm TAR1 prevalence (% sequences carrying TAR1) sorted ascending; PAR1 arms (chrXp 0.3%, chrYp 1.1%) highlighted in blue against grey for autosomal arms (≥73%); all-PHR mean 94.6% as grey dashed reference.

Data sources (all already on disk; the R script reads them directly):

- `/moosefs/guarracino/HPRCv2/PHR_III/plots/d4z4_dux4l_by_community.tsv` (1,253 chr4q/chr10q haplotypes)
- `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv` (16 arms, 5,023 OR4F entries)
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_tar1_by_arm.tsv` (41 arms)

Cross-reference: `paper_prep/figures/ed4/figure_ed4.pdf` panels c (high-copy gene families) and d (OR4F pseudogenisation gradient) are the manuscript versions of the OR4F + DUX4 content; this slide's figure is the talk version (one row, three panels, larger labels) and adds the TAR1 panel that ED4 omits. Do **not** swap them — keep ED4 for the manuscript, slide_14 figure for the talk.

## Speaker notes

A quick aside before I move to the punchline. The talk has been about geometry — communities, exchange, scale — but the biology inside these regions is genuinely interesting in its own right, and worth flagging. Three vignettes, ten seconds each.

First — DUX4. We see DUX4 annotated across 18 q-arms in the pangenome. But the medically relevant biology — the D4Z4 macrosatellite array whose contraction causes facioscapulohumeral muscular dystrophy — lives only at chr4q and chr10q, the C1 community. Median 22 copies of DUX4L per haplotype on those two arms; everywhere else, just 0–2. So DUX4 has spread, but FSHD is geometrically constrained to one community out of fifteen.

Second — OR4F olfactory receptors. Four paralogs, sixteen arms, five thousand gene-copy entries. The pseudogenisation rate sweeps from eleven percent at chr7p to ninety-nine point eight percent at chr15q. Same gene family, distributed across the same subtelomeric exchange network — but the decay clock has been running for very different lengths of time at different ends.

Third — TAR1, the telomere-associated repeat. Ninety-four point six percent of all our subtelomeric sequences carry TAR1, across all forty-one arms. The one exception is PAR1 — chrX_p and chrY_p — at half a percent. PAR1 has obligate meiotic crossover. So the satellite seems to mark places that *need* a sequence-anchor for exchange; PAR1 doesn't, and TAR1 isn't there.

Three different biological readouts — disease, decay, exchange machinery — and all three write themselves into the architecture we just spent the talk mapping. **This is a digression from the core argument**, and if I'm running short I'll skip it. But if you remember one thing: this is the substrate the rest of human biology cares about — FSHD lives here, the olfactory repertoire ages here, and the exchange machinery itself leaves footprints here.

## Time budget

50 seconds. **Compressible to 0** — explicitly skippable if the previous slides have run long. If kept, run all three vignettes at ~15 s each plus 5 s framing.

## Notes for synthesizer

- **Framing — this is an ASIDE, not the core argument.** The slide opens "and the biology is interesting too" and the speaker note ends with "this is a digression from the core argument; if I'm running short I'll skip it." The synthesizer should preserve this skippable framing. The talk's main argument is geometry / exchange / population-scale (slides 1–13 and 15); slide 14 is a "look how rich this is" detour. Do not promote any of the three vignettes to a load-bearing claim of the talk.
- **Compressible to 0.** This is the cleanest slide to drop if the talk runs long. The 15-min slot is tight (15 slides, ~60 s each on average, with anchor slides — title, central PHR scale, PAR/PHR comparison, conclusions — needing more than 60 s). If timing is tight, slide 14 is the cut. The synthesizer should flag this explicitly to the speaker (Erik Garrison).
- **Numbers to lock down (single source of truth: ABSTRACT.md, CROSSWALK §C12, end-to-end-report 02_annotation.md and 03_gene_enrichment.md):**
  - DUX4: 18 q-arms (annotation), C1 = chr4q + chr10q (the D4Z4 community), median 22 DUX4L per haplotype on C1, 0–2 elsewhere. Andrea ch. 03 cites "DUX4L1–DUX4L44 (28 genes, 2 arms)" — the 28 distinct DUX4L paralog labels are the canonical CHM13 annotation; the 22-median is the per-haplotype copy count from `d4z4_dux4l_by_community.tsv`. These are not in conflict — different denominators.
  - OR4F gradient: 11.1% (chr7p) → 99.8% (chr15q), 16 arms, 5,023 entries, mean 62.1%. Verbatim from CROSSWALK §C12 / SURVEY_10_11_12 line 344.
  - TAR1: 94.6% of 15,668 sequences carry TAR1 across 41 arms; PAR1 = 0.3% (chrXp) and 1.1% (chrYp). 21,544 total TAR1 entries. From end-to-end-report 02_annotation.md ("TAR1 prevalence" subsection).
- **Figure provenance.** The slide figure is **slide-specific** (newly generated by `slide_14_gene_biology.R` in this worktree, using existing cluster TSVs — no SBATCH, no new data). It overlaps panels c+d of `paper_prep/figures/ed4/figure_ed4.pdf` (manuscript ED4) but adds the TAR1 panel and reformats for talk display. The synthesizer should not substitute ED4 verbatim — ED4 has 4 panels including a GO:BP plot and a copy-weighted comparison panel that are off-topic for this slide.
- **R script renders cleanly.** `Rscript slides/v2/slide_14_gene_biology.R` regenerates `.pdf` and `.png` from the cluster TSVs in ~3 s; no external dependencies beyond base R + readr/dplyr (both already in the project's R environment). Do not require patchwork / cowplot / gridExtra — they are NOT installed in the worktree's R environment.
- **Callbacks/setup.** The previous slide(s) (12, 13) cover within-community heterogeneity / cross-arm exchange status; slide 14 picks up "what's actually IN these regions?" framing without claiming any of it as a finding the speaker needs to defend. Slide 15 (closing / conclusions) should NOT cite slide 14 numbers as load-bearing — slide 14 is decorative.
- **Forward setup.** Slide 15 (concluding slide) returns to the geometry argument. Slide 14's three vignettes don't feed into it; slide 14 is a parallel "by the way" detour and the speaker explicitly returns to the main thread on slide 15.
- **PAR1 callback.** The TAR1 panel says "PAR1 is the only place that doesn't need a satellite-mediated exchange anchor." This implicitly connects to the PAR-vs-PHR framing established earlier in the deck (slide 05's "PAR2-class exchange landscape replicated 41 times"). The synthesizer can choose to make this callback more explicit if desired — but the bullet/note already gestures at it.
- **No new R packages required.** Verified the worktree has base R + ggplot2 + dplyr + readr but NOT patchwork / cowplot / gridExtra. Any future re-cut should stay within base R `layout()` / `par(mfrow=...)`.
