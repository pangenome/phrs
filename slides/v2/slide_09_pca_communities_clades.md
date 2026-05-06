## Title
All-vs-all PCA — 15 arm-level communities, named for the abstract's clades (the keystone slide)

## Bullets
- Same PCA layout as v1 slide 10: each point is one chromosome arm in PHR-similarity space; clusters = Leiden k=15 communities (C1–C15) on the 41×41 arm-level Jaccard distance matrix (chapter 01 §"Community detection"; arms `hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`).
- **Every named clade in the abstract is a community on this plot.** PAR1 = **C15** (Xp/Yp + 18q outlier). PAR2 = **C14** (Xq/Yq). Acrocentric short arms = **C7** (13p,14p,15p,21p,22p). The abstract's 22q/21q/19q/1q/13q/17q clade = **C6** (exact match). 10p–18p Linardopoulou pair = **C2**. 4q–10q DUX4 = **C1**.
- Three interpretive zones replace v1's question-marked overlays: **PAR-driven** (lower-left, C14+C15), **concerted-exchange PHR core** (center, C6 + C7 + C3), **DUX4 / D4Z4** (right, C1 — confidence upgraded; CROSSWALK §C5 ties it to chapter 03's 28 DUX4L genes and chapter 04's 43.4 % type-discordance).
- Same PC1 / PC2 axes as v1 (16.05 % / 11.2 %); same five missing introvert arms (2p, 3p, 5p, 8q, 11q, 14q). UPGMA k=14 recovers 12/15 of these communities — independent confirmation that the partition is not a Leiden artefact.
- **One slide that ties methods to abstract.** Every clade word the audience just read in the abstract has a colored cluster here, with an arm list, with a recombinational interpretation. After this slide the rest of the talk is *evidence for* what is grouped on this plot.

## Primary figure
**Reuse:** `slides/20260204_Subtelomics_overview_EG.pdf` page 10 — keep the PCA scatter, the chromosome-arm point labels, the PC1 / PC2 axes, and the three interpretive arrows verbatim. Do **not** re-render the scatter (the data lives outside the worktree at `/moosefs/guarracino/HPRCv2/PHR_III/similarity/…`).

**One new asset (this slide's deliverable):** replace v1's right-margin community list with an abstract-anchored legend table. Renders to PDF / PNG locally with vanilla ggplot2 (no `gridExtra`, no `cowplot` — confirmed working in this worktree's R install). Place to the right of (or beneath) the imported v1 PCA.

```r
# slide_09_clade_legend.R
# Renders the abstract-anchored legend that replaces the v1 community list on
# the right margin of the PCA scatter. Runs in seconds, no data dependencies,
# pure ggplot2.
library(ggplot2)

leg <- data.frame(
  row   = 15:1,                       # so C1 lands at top of the rendered table
  C     = sprintf("C%d", 1:15),
  arms  = c("4q, 10q",
            "10p, 18p",
            "3q, 7p, 9q, 11p, 16q, 19p",
            "7q, 12q",
            "6p, 9p, 12p, 20q",
            "1q, 13q, 17q, 19q, 21q, 22q",
            "13p, 14p, 15p, 21p, 22p",
            "15q",
            "16p",
            "17p",
            "1p, 5q, 6q, 8p",
            "2q, 20p",
            "4p",
            "Xq, Yq",
            "18q (n=1), Xp, Yp"),
  clade = c("4q–10q DUX4 / D4Z4",
            "10p–18p (Linardopoulou 2005)",
            "f7501 duplicons (fixed + AFR-enriched)",
            "private 7q/12q pair",
            "RPL23A / WASH duplicons",
            "concerted q-arm clade (22/21/19/1/13/17q)",
            "acrocentric short arms (rDNA-adjacent)",
            "chr15_q (single arm)",
            "chr16_p (single arm)",
            "chr17_p (single arm)",
            "OR4F21 sharing (Linardopoulou block 5)",
            "2q/20p pair",
            "chr4_p (single arm)",
            "PAR2 (Xq/Yq)",
            "PAR1 (Xp/Yp)"),
  fill  = c("#FDE2C8","#FFF1AA","#FFFFFF","#FFFFFF","#FFFFFF",   # C1 DUX4, C2 10p18p
            "#D6E8FF","#E5D8FA","#FFFFFF","#FFFFFF","#FFFFFF",   # C6 q-arm, C7 acro
            "#FFFFFF","#FFFFFF","#FFFFFF","#CDEAD3","#CDEAD3"),  # C14 PAR2, C15 PAR1
  stringsAsFactors = FALSE
)

cols <- data.frame(
  x      = c(0.05, 0.30, 0.70),
  width  = c(0.22, 0.38, 1.05),
  field  = c("C","arms","clade"),
  stringsAsFactors = FALSE
)

cells <- do.call(rbind, lapply(seq_len(nrow(cols)), function(j) data.frame(
  x        = cols$x[j] + cols$width[j]/2,
  width    = cols$width[j],
  y        = leg$row,
  text     = leg[[ cols$field[j] ]],
  fill     = leg$fill,
  fontface = ifelse(j == 1, "bold", "plain"),
  stringsAsFactors = FALSE
)))

hdr <- data.frame(
  x        = cols$x + cols$width/2,
  width    = cols$width,
  y        = 16,
  text     = c("Community","Arms","Abstract clade / interpretation"),
  fill     = "#EEEEEE",
  fontface = "bold",
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
            hjust = 0, size = 3.3) +
  coord_cartesian(xlim = c(0, 1.7), ylim = c(0.5, 16.5), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(6, 6, 6, 6))

ggsave("slide_09_clade_legend.pdf", p, width = 8.0, height = 5.5)
ggsave("slide_09_clade_legend.png", p, width = 8.0, height = 5.5, dpi = 200)
```

The colored rows are exactly the abstract-anchored six clades (C1 / C2 / C6 / C7 / C14 / C15); everything else stays neutral so the eye lands on the rows that match the abstract's wording. The synthesizer can dock this PNG to the right of the v1 PCA scatter (the shape — 8.0 × 5.5 inches — matches the v1 right margin).

## Speaker notes
This is the keystone of the talk — every clade word in the abstract is a colored cluster on this plot. The PCA is on the 41×41 arm-level Jaccard distance matrix; clusters are Leiden communities, k=15 chosen by silhouette. UPGMA on the same matrix recovers twelve of these fifteen — the partition is not an algorithm artefact.

Walk the audience around the plot in the order of the abstract. The pseudoautosomals — PAR2 is **C14**, Xq with Yq; PAR1 is **C15**, Xp with Yp. Acrocentric short arms — **C7** — thirteen-p, fourteen-p, fifteen-p, twenty-one-p, twenty-two-p, all five of them, the rDNA-adjacent homogenization clade. The novel q-arm clade the abstract calls "twenty-two-q, twenty-one-q, nineteen-q, one-q, thirteen-q, seventeen-q" is **C6** — an exact six-arm match, the largest non-acrocentric inter-chromosomal community we see. Ten-p with eighteen-p — Linardopoulou's 2005 pair — is **C2**. And four-q with ten-q, the DUX4/D4Z4 pair, is **C1** — what v1 marked with a question mark is now an established clade in the abstract, twenty-eight DUX4L genes, copy-number diversity, type-discordance forty-three percent.

So this slide says: the abstract's vocabulary maps onto a single empirical structure. Everything that follows is evidence *about* clusters you can already see here.

## Time budget
80 seconds. Roughly: 10 s on "every point is an arm, clusters are communities, k=15 by silhouette, UPGMA agrees on 12/15"; 50 s walking the six abstract-anchored communities (PAR2, PAR1, acrocentric, q-arm clade C6, 10p–18p, 4q–10q DUX4) — about 8 s each, name the clade, point at the cluster, name the arms; 20 s closing line "everything that follows is evidence about clusters you can already see here" + segue to the next slide.

## Notes for synthesizer
- **Layout is preserved from v1 slide 10.** Do not re-render the PCA scatter — keep the v1 panel as-is (PC axes, point shapes, point labels, the three interpretive arrows). The only on-slide change is the right-margin legend, which the R script above generates.
- **Update the three interpretive arrows from v1**: "PARs-driven" → "PAR1 / PAR2 (C15 / C14)"; "PHRs-driven" → "concerted-exchange core (C6 + C7 + C3)"; "DUX4-driven?" → "4q–10q DUX4 (C1)" — *drop the question mark*; the abstract treats it as established.
- **Cross-slide callbacks:**
  - **Slide 07 (v1 all-vs-all heatmap, v2 slide 07/08)** sets up the same arm-level distance matrix; this slide is the projection of that matrix into 2-D.
  - **Slide 08 (PCA by superpop, v1 slide 9)** sets up the projection coordinates without color-coding by community; this slide reuses those coordinates with the community color-coding.
  - **Outbound to slide 10 (within-community heterogeneity / Fig 2a)** — the speaker should land "C7 is the only community where paralog distance is less than allelic distance — the acrocentric p-arms are *fully* homogenized; we will see that on the next slide." That is the cleanest hand-off.
  - **Outbound to slide 11 (population structure / Fig 2c-d)** — communities here become the units for cross-arm Fst and the out-of-Africa tree.
- **Abstract-clade ↔ community map (the load-bearing block — copy verbatim if compressing):**
  - C1 (4q, 10q) ↔ "4q–10q DUX4 with copy-number diversity"
  - C2 (10p, 18p) ↔ "10p–18p" (Linardopoulou 2005 Fig. 5)
  - C6 (1q, 13q, 17q, 19q, 21q, 22q) ↔ "tightly linked clade involving 22q, 21q, 19q, 1q, 13q, and 17q" (exact six-arm match)
  - C7 (13p, 14p, 15p, 21p, 22p) ↔ "acrocentric short arms"
  - C14 (Xq, Yq) ↔ "Xq/Yq via the pseudoautosomal regions" → PAR2
  - C15 (Xp, Yp; 18q n=1) ↔ "Xp/Yp via the pseudoautosomal regions" → PAR1
  - Source: `paper_prep/synthesis/CROSSWALK.md` §C5; `end-to-end-report/report/01_pipeline.md` §"Arm-level community detection".
- **MDS-vs-PCA framing caveat (CROSSWALK §C6 / §3 row):** Andrea's report uses MDS / PCoA on the Jaccard distance, not strict PCA on a feature matrix. The v1 slide called it "PCA" and the abstract uses "principal component … analyses". For the BoG talk this is fine — MDS on a Jaccard distance is metric-equivalent to PCoA, and the audience reads "PCA" as "the 2-D projection of a similarity structure". Speaker should use "PCA" verbally; the methods writer (REWRITE_PLAN TASK-16) is responsible for resolving this in the manuscript.
- **Missing introvert arms (2p, 3p, 5p, 8q, 11q, 14q)** — keep the upper-right callout from v1; do NOT reframe. Their absence is not a quality issue, it is a biological signal (these arms have ~no inter-chromosomal hits, so they have no row in the 41×41 matrix). One sentence in the speaker notes if a Q&A pushes on it.
- **Do not** re-color the points by superpopulation here — that is the previous slide's job. Mixing collapses both slides into one and erases the keystone framing.
- **If the time budget is cut** (e.g., 60 s total for the deck): drop C3 and C5 from the spoken walk (they are duplicon-sharing groups, not in the abstract clade list); never drop PAR2 / PAR1 / C7 / C6 / C2 / C1 — those are the abstract-anchored six.
- **R script is light** (renders a 7.5×5 inch table grob in seconds; no data files, no SBATCH). If the script fails for missing `gridExtra`, fall back to a hand-typed table block in Keynote — content is what matters, rendering pathway is interchangeable.
