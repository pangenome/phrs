## Title

**All-vs-all in 2D — colored by chromosome (left) vs superpopulation (right)**

Subtitle: same projection of the 15,668-sequence Jaccard distance matrix; left panel asks *what arm is it from?*, right panel asks *what population?* — only the left answer holds up.

## Bullets

- Two-panel layout, **identical points**, only the coloring changes. Left: each point is one of 15,668 subtelomeric flanks colored by source chromosome (chr1–22, X, Y); shape encodes p / q arm. Right: same scatter recolored by 1KGP superpopulation (AFR / AMR / EAS / EUR / SAS).
- **Left panel: clusters are arm-community shaped.** Points group by Leiden community (15 communities over 41 arms — D4Z4 chr4_q+chr10_q, acrocentric p-arms chr13/14/15/21/22 p, PAR1 chrX_p+chrY_p, PAR2 chrX_q+chrY_q, the 6-arm q-arm community chr1_q/13_q/17_q/19_q/21_q/22_q, etc.). Visible structure ≈ arm-community structure.
- **Right panel: superpopulations are mixed across all clusters.** The arm-community clusters do not split by AFR/AMR/EAS/EUR/SAS — population structure is **real but secondary** to arm-community structure.
- **Population structure *is* there, just smaller.** Hudson Fst on cross-arm affinity (chapter 04): mean **Fst = 0.044**; AFR vs non-AFR pairs **0.10–0.15**, non-AFR–non-AFR pairs **−0.05 to 0.02**. The AFR-deepest split mirrors the human out-of-Africa tree (chapter 12 novel contribution #19) — a population signal exists, but it lives at finer scale than the dominant arm-community signal driving the global 2D layout.
- **Missing introvert arms** (no inter-chromosomal PHR detected, excluded from this projection): **2p, 3p, 5p, 8q, 11q, 14q** — same six arms across both panels.

## Primary figure

Two-panel side-by-side. **Both panels are existing artifacts in Andrea's pipeline output:**

- Left panel — by chromosome: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-chromosome.pdf` (PNG companion: `hprcv2.1Mb.subtelo.mds.color-by-chromosome.png`).
- Right panel — by superpopulation: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-superpopulation.pdf` (PNG companion: `hprcv2.1Mb.subtelo.mds.color-by-superpopulation.png`).

**Important — what these are.** Both produced by `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R` line 556 — `cmdscale(as.dist(jaccard_dist_df), eig = TRUE, k = 5)` on the 15,668 × 15,668 Jaccard distance matrix from `pggb -p 95` + `odgi similarity --all -P`. Variance-explained labels on the axes come from `fit_full$eig / sum(abs(fit_full$eig)) * 100` (script line 579). This is **classical MDS (= PCoA on a Jaccard distance matrix)**, not PCA on a feature matrix — see Notes for synthesizer.

If a side-by-side composite is needed for the slide deck (no SBATCH; runs in agent worktree):

```r
# slide_08_pca_combined.R — side-by-side composite of color-by-chromosome
# and color-by-superpopulation MDS panels. Inputs are Andrea's existing
# panel PDFs; this script just stitches them.
library(magick)
library(grid)
library(gridExtra)

base <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity"
left  <- image_read_pdf(file.path(base, "hprcv2.1Mb.subtelo.mds.color-by-chromosome.pdf"),       density = 200)
right <- image_read_pdf(file.path(base, "hprcv2.1Mb.subtelo.mds.color-by-superpopulation.pdf"),  density = 200)

# pad to equal height, then concatenate horizontally
h     <- max(image_info(left)$height, image_info(right)$height)
left  <- image_extent(left,  geometry_size_pixels(width = image_info(left)$width,  height = h),
                     gravity = "center", color = "white")
right <- image_extent(right, geometry_size_pixels(width = image_info(right)$width, height = h),
                     gravity = "center", color = "white")
combo <- image_append(c(left, right))
image_write(combo, "slide_08_pca_combined.png", format = "png")
```

**Annotation overlay (recommended)** — drop a single shared caption beneath the two panels: *"Same MDS / PCoA projection (cmdscale on 1 − Jaccard, k = 5); n = 15,668 flanks across 41 arms; missing introvert arms: 2p, 3p, 5p, 8q, 11q, 14q."* Speaker delivers the "left = arm-community, right = mixed → population structure is secondary" beat verbally.

## Speaker notes

This is the same all-vs-all object as the heatmap two slides back, just projected into two dimensions so we can ask two questions of one picture. Each point is one of 15,668 subtelomeric flanks. The projection is classical MDS — PCoA on the Jaccard distance matrix from the pangenome graph. **One scatter, two colorings.**

On the left I've colored by source chromosome. The clusters you see are arm-community-shaped: D4Z4 pulls chr4q and chr10q together; the acrocentric p-arms — chr13p, 14p, 15p, 21p, 22p — collapse into one cluster; PAR1 puts Xp on top of Yp; the six-arm q-arm clade with chr1q, 13q, 17q, 19q, 21q, 22q sits as its own cloud. The structure on the left panel is the arm-community structure from slide 07's heatmap, replotted.

On the right I've taken the *same* points and recolored them by 1000 Genomes superpopulation — AFR, AMR, EAS, EUR, SAS. The colors mix across every cluster. The clusters do **not** resolve by population. So the headline is: in this 2D view, what dominates is *which chromosome arm a sequence comes from*, not *which population the haplotype is from*.

That doesn't mean population structure is absent — it means it's secondary. Chapter 4 of our analysis tests it directly with Hudson's Fst on cross-arm affinity: mean Fst is 0.044, AFR-versus-non-AFR pairs sit at 0.10 to 0.15, and non-AFR pairs sit at zero. Same pattern you'd recognize from any human-population genetics study — AFR is the deepest split, mirrors out-of-Africa. So the population signal is present and quantifiable; it just doesn't dominate the global 2D layout, because the arm-community signal is much stronger. **Two-tier hierarchy: arm-community first, population structure within.**

## Time budget

**Target: 70 seconds.** Roughly 15 s setting up the two-panel framing ("same plot, two colorings, different questions"); 25 s on the left panel — call out 3–4 specific clusters by name (D4Z4 4q/10q, acrocentric p, PAR1, q-arm clade); 20 s on the right panel mixing → "population structure is real but secondary"; 10 s landing the Fst numbers (0.044 mean, 0.10–0.15 AFR vs non-AFR, out-of-Africa shape).

## Notes for synthesizer

- **CRITICAL — PCA vs MDS labeling (CROSSWALK §3 C6).** The v1 deck slides 8 and 9 both said "PCA". The artifact is **MDS / PCoA**, computed by `cmdscale(as.dist(jaccard_dist_df), eig = TRUE, k = 5)` at `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R` line 556, output cached to `hprcv2.1Mb.subtelo.full_mds.rds`. There is **no `*pca*.rds` in the HPRCv2 tree** (verified by `find /moosefs/guarracino/HPRCv2 -name "*pca*"` — only `HG002.GRCh38_no_alt.dipcall.bed` matches, unrelated). I have re-titled this slide accordingly. Two paths going forward; the synthesizer must pick one and propagate:
    - **(a) Relabel only.** Change abstract wording from "Principal component and community detection analyses" → "MDS / PCoA and community detection analyses on the Jaccard distance matrix". Erik's clarification per CROSSWALK is "abstract is not locked"; this is the lower-friction option and it is what the existing artifact actually shows.
    - **(b) Run a proper PCA** on a haplotype × PHR-presence-binary feature matrix (REWRITE_PLAN TASK-16). This produces a strict PCA but adds work and a new artifact; the *biological* conclusion (arm-community first, population secondary) does not change.
    - **My recommendation: (a)** for the talk this week, with TASK-16 as the manuscript follow-up. The talk slide says "MDS / PCoA"; if the abstract is later rewritten to PCA, the slide can be relabeled.
- **The v1 axis labels were "Dimension 1 (16.05%)" and "Dimension 2 (11.2%)".** Those numbers come from `var_explained_full[1]` / `var_explained_full[2]` in the script (line 579–581) — i.e., from `fit_full$eig / sum(abs(fit_full$eig)) * 100`. They are the **MDS** variance-explained, not PCA variance-explained. They are real numbers from the actual run; keep the values, just relabel the projection.
- **Combining v1 slides 8 + 9 into one combined slide is the right move** — they share axes, share points, share the missing-introvert-arms callout in the corner, and the comparison is the *content*: arm coloring "explains" the layout while superpop coloring does not. Showing them as separate slides loses 30+ seconds of redundant setup.
- **Callbacks.** Slide 07 (all-vs-all heatmap with arm-level dendrogram) sets up the same Jaccard distance matrix at the arm × arm level; this slide unfolds the underlying *sequence × sequence* matrix into 2D. Pre-cue verbally: "same matrix, zoomed out from 41 × 41 arm averages to 15,668 individual flanks."
- **Forward setup for slide 09 (community coloring) and slide 10 (DUX4 / acrocentric / PAR mechanism slides).** This slide deliberately uses the **chromosome** and **superpopulation** colorings, not the **community** coloring, so that slide 09 can introduce the Leiden communities as the *labeled* version of what the audience just saw cluster naturally. Do NOT pre-empt the community coloring here.
- **Hudson Fst footnote.** The numbers (mean 0.044, AFR vs non-AFR 0.10–0.15, non-AFR pairs −0.05 to 0.02, AMR-EUR closest at Fst 0.22 in the out-of-Africa tree) are CROSSWALK §1 ch.04 + ch.12 / end-to-end-report 04_heterogeneity.md lines 105–113. They support the abstract's "resolve subtelomeric clustering across human populations" claim *independent* of whether the projection is PCA or MDS. If the synthesizer rewrites the abstract per (a) above, the C6 Results paragraph should still cite these Fst numbers — they are the substance of the population-clustering claim.
- **Missing introvert arms (2p, 3p, 5p, 8q, 11q, 14q).** Six arms with no detected inter-chromosomal PHR. These are absent from the projection by construction (no Jaccard rows). Mention as a one-line corner annotation; do NOT spend speaker time explaining unless asked in Q&A — the explanation belongs in slide 11/12 (community detection caveats).
- **Figure provenance.** Both panel PDFs already exist in Andrea's PHR_III/similarity/ directory; no new alignment, no SBATCH, no graph rebuild needed. The optional R stitching script is light I/O only and runs in seconds inside the agent worktree.
- **Numbers to lock down (single source of truth).** 15,668 flanks / 41 of 48 arms / 6 missing introvert arms / Dim1 16.05% / Dim2 11.2% / Hudson Fst mean 0.044, AFR vs non-AFR 0.10–0.15. All consistent with CROSSWALK §3 C6 + end-to-end-report 01_pipeline.md and 04_heterogeneity.md.
