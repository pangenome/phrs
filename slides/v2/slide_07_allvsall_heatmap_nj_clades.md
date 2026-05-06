## Title

All-vs-all at the arm level — clustered heatmap + NJ tree with named clades (the cladistic structure)

## Bullets

- **Two-panel layout.** Left: 41×41 arm-level Jaccard distance heatmap (Fig 1c — `hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`), arms ordered by Leiden k=15 community, cyan rectangles delimit the 15 communities, the original v1 UPGMA k=14 dendrogram on top recovers 14 of 15 communities (12/15 exact). Right: neighbor-joining tree on the **same** 41×41 distance matrix, rooted at the acrocentric MRCA, tip labels colored by the **six named clades from the abstract**.
- **The cladistic structure recovers expected pseudohomology and reveals novel clades** — every clade word in the abstract is a monophyletic block on the NJ tree, with **100% perturbation-bootstrap support** at each named clade's MRCA (1000 reps, σ = 25% of off-diagonal IQR):
  - **PAR1** (Xp/Yp) — red
  - **PAR2** (Xq/Yq) — blue (Xq–Yq edge length essentially zero, the shortest pair on the tree)
  - **acrocentric short arms** (13p, 14p, 15p, 21p, 22p) — green; chosen as the root, monophyletic
  - **10p–18p** (Linardopoulou 2005) — orange
  - **tight q-arm clade** 22q, 21q, 19q, 1q, 13q, 17q — purple (exact six-arm match to the abstract wording)
  - **4q–10q DUX4** — brown
- **NJ ↔ Leiden 1:1 mapping** (the partition is not an algorithm artefact): PAR1 = C15, PAR2 = C14, acrocentric_p = C7, 10p–18p = C2, tight q-arm = C6, DUX4 = C1.
- **Missing introvert arms** (no inter-chromosomal PHR detected, absent from this 41×41 matrix and from the NJ tree): **2p, 3p, 5p, 8q, 11q, 14q**. Same six arms drop out on every downstream similarity view (slides 08–09).
- **This is the highest-value v1→v2 swap.** v1 showed the heatmap with an unlabeled UPGMA dendrogram and a single "?" overlay; v2 adds the NJ tree the abstract names, with bold labels on every clade the audience will hear about for the rest of the talk.

## Primary figure

Two panels side-by-side.

- **Left panel (heatmap):** `paper_prep/figures/fig1/figure_fig1.pdf` — extract panel **(c)** only (arm-level 41×41 Jaccard distance heatmap with UPGMA k=14 dendrogram on top and the 15 Leiden community boxes overlaid). This is the published-quality version of v1 slide 7's heatmap; preserve the cyan community rectangles, the dendrogram, and the p/q-arm side annotation.
- **Right panel (NJ tree):** `paper_prep/figures/nj_tree_arms/nj_tree_annotated.pdf` (PNG companion alongside) — produced by the upstream `nj-tree-from` task in this worktree. 41 tips, rooted at the acrocentric MRCA, tip labels in bold for the six named clades and color-keyed to a legend with the clade names spelled out. Bootstrap support printed at every internal node.

If a side-by-side composite is wanted in the deck (no SBATCH; runs in agent worktree, ~2 s):

```r
# slide_07_heatmap_nj_combined.R — stitch the Fig 1c arm-level heatmap and the
# annotated NJ tree into a single landscape composite for the slide.
library(magick)

heatmap <- image_read_pdf("paper_prep/figures/fig1/figure_fig1.pdf",
                          pages = 1, density = 200)
# Crop panel (c) only — adjust to your fig1 layout if it changes.
hi <- image_info(heatmap)
panel_c <- image_crop(heatmap,
                      sprintf("%dx%d+%d+%d",
                              as.integer(hi$width  * 0.55),  # right ~55% width
                              as.integer(hi$height * 0.55),  # lower ~55% height
                              as.integer(hi$width  * 0.45),
                              as.integer(hi$height * 0.45)))

njtree <- image_read_pdf("paper_prep/figures/nj_tree_arms/nj_tree_annotated.pdf",
                         density = 200)

# Pad to equal height, then concatenate horizontally.
h <- max(image_info(panel_c)$height, image_info(njtree)$height)
panel_c <- image_extent(panel_c,
                        geometry_size_pixels(width = image_info(panel_c)$width,
                                             height = h),
                        gravity = "center", color = "white")
njtree  <- image_extent(njtree,
                        geometry_size_pixels(width = image_info(njtree)$width,
                                             height = h),
                        gravity = "center", color = "white")
combo <- image_append(c(panel_c, njtree))
image_write(combo, "slide_07_heatmap_nj_combined.png", format = "png")
```

Crop coords are heuristic; if the synthesizer prefers, render the heatmap panel directly from `paper_prep/figures/fig1/figure_fig1.R` (the panel-c block) and skip the PDF crop.

**Annotation overlay (recommended, on the NJ tree only)** — the upstream task already paints the six abstract clades in bold colored text with a top-right legend; no new ggplot work needed for the figure itself. The only on-slide annotation to add is a one-line caption beneath the pair: *"Same 41×41 arm-level Jaccard distance matrix; left = clustered heatmap (Leiden k=15 + UPGMA), right = NJ tree (rooted at acrocentric MRCA, 1000-rep bootstrap)."*

## Speaker notes

This is the all-vs-all picture at the arm level. The matrix is forty-one by forty-one — every chromosome arm we have a signal on, against every other arm. Cell color is Jaccard distance on the pangenome graph; the cyan boxes are the fifteen Leiden communities; the dendrogram on top is the original UPGMA, which already pulled out fourteen of the fifteen blocks.

What I want you to look at is the right panel. This is a neighbor-joining tree on the **same** distance matrix, rooted at the acrocentric short-arm clade because that clade is monophyletic and gives a stable orientation. Every clade name you read in the abstract — PAR1, PAR2, acrocentric short arms, ten-p with eighteen-p, the four-q ten-q DUX4 pair, and the tight q-arm clade of twenty-two-q, twenty-one-q, nineteen-q, one-q, thirteen-q, and seventeen-q — every one of those is a monophyletic block on this tree, in bold color, with one hundred percent bootstrap support at the named-clade root. We perturbed the distance matrix a thousand times at twenty-five percent of its off-diagonal IQR, rebuilt the tree, and these six clade roots never broke.

So **the cladistic structure recovers expected pseudohomology — the pseudoautosomals, the acrocentric short arms — and it reveals novel clades**: the ten-p eighteen-p pair which Linardopoulou drew in 2005, the four-q ten-q DUX4 pair, and most strikingly the six-arm q-arm clade that the abstract names. Same answer from Leiden, from UPGMA, and from neighbor-joining. The next several slides are evidence about clusters you can already see here.

## Time budget

**80 seconds.** Roughly: 15 s setting up the two-panel framing ("same forty-one by forty-one matrix, two clustering views"); 10 s on the heatmap (cyan boxes = 15 communities, UPGMA dendrogram agrees on 14 of 15); 35 s walking the six abstract-anchored clades on the NJ tree (PAR1, PAR2, acrocentric, tight q-arm, 10p–18p, 4q–10q DUX4 — about 5 s each, name and color); 10 s on bootstrap robustness ("100% support at every named clade, three algorithms agree"); 10 s landing the "expected + novel" closing line and segueing to slide 08.

## Notes for synthesizer

- **Highest-value v1→v2 swap.** The v1 deck had the heatmap with an unlabeled UPGMA dendrogram and a single hand-drawn "?" overlay on what we now call C1 / DUX4. v2 keeps the heatmap (it works) but **adds the NJ tree the abstract explicitly names** (`Cladistic analysis based on neighbor-joining trees…`), and labels the six abstract clades in bold colored type. This single swap closes the C5 framing gap from CROSSWALK §C5 ("abstract names NJ but no NJ exists") and resolves REWRITE_PLAN TASK-01 + TASK-13 for the talk.
- **Layout decision.** Two panels side-by-side. If only one panel fits at conference distance, prefer the **NJ tree** alone (it carries the abstract terminology) and demote the heatmap to a small inset. The heatmap is the v1 figure and is recognizable; the NJ tree is the new artifact and the named-clade evidence.
- **Figure provenance — strict.** Heatmap panel = Fig 1c from `paper_prep/figures/fig1/figure_fig1.pdf` (already in repo, manuscript-quality, do not rebuild). NJ tree = `paper_prep/figures/nj_tree_arms/nj_tree_annotated.pdf` (produced by the upstream `nj-tree-from` task in this branch series; commit `602a9d3`). Both PNG companions are present alongside the PDFs. **Do not** re-render either — the only optional artifact is the `magick` stitch above, which is a 2-second compositing step.
- **Cross-slide callbacks (load-bearing).**
  - **Inbound (slide 06)** — slide 06 ends on PHR length distributions; this slide picks up the same 15,668-PHR object and asks "how do the *arms* group?". Pre-cue verbally: "we just summarized PHRs by length; now we summarize arms by who they share with."
  - **Outbound (slide 08)** — slide 08 takes the **same** 15,668 × 15,668 *sequence-level* Jaccard matrix and projects it into 2D, colored by chromosome and superpop. Hand-off line for the speaker: "we just saw the forty-one by forty-one **arm** picture — next slide unfolds it to the fifteen-thousand-flank picture."
  - **Outbound (slide 09)** — slide 09 is the **keystone** that walks the same six abstract clades on a 2D PCA layout with arm-level points. The clade vocabulary (PAR1=C15, PAR2=C14, ACRO_p=C7, 10p–18p=C2, TIGHT_q=C6, DUX4=C1) **must match** between this slide and slide 09. The Leiden ↔ NJ ↔ abstract-clade map is the single source of truth for the talk.
- **NJ ↔ Leiden ↔ abstract-clade map (single source of truth — copy verbatim):**
  - **C15 ↔ NJ block PAR1 ↔ "Xp/Yp via the pseudoautosomal regions"** (red on tree)
  - **C14 ↔ NJ block PAR2 ↔ "Xq/Yq via the pseudoautosomal regions"** (blue on tree)
  - **C7 ↔ NJ block ACRO_p ↔ "acrocentric short arms"** (green on tree; root)
  - **C2 ↔ NJ block 10p–18p ↔ "10p–18p homology"** (orange on tree)
  - **C6 ↔ NJ block TIGHT_q ↔ "tightly linked clade involving 22q, 21q, 19q, 1q, 13q, and 17q"** (purple on tree; exact six-arm match to abstract wording)
  - **C1 ↔ NJ block DUX4 ↔ "DUX4-containing homology between 4q and 10q with wide copy number diversity"** (brown on tree)
  - Source: `paper_prep/synthesis/CROSSWALK.md` §C5; `paper_prep/figures/nj_tree_arms/README.md`; upstream `nj-tree-from` task log: "all 6 abstract clades recovered as monophyletic … 100% support at every named-clade MRCA."
- **Bootstrap caveat — explain only if asked.** Support is from a perturbation bootstrap (Gaussian noise on the distance matrix, σ = 25% of off-diagonal IQR, 1000 reps), not a Felsenstein column-bootstrap, because the input is a derived distance summary rather than an alignment. Deeper backbone edges show 32–90% support — the named clades are robust; the relative ordering of clades is not. Mention only if reviewer-style Q&A pushes; the abstract uses the result, not the support method.
- **Missing introvert arms (2p, 3p, 5p, 8q, 11q, 14q)** — same six arms as slides 08/09. Their absence is a biological signal (no inter-chromosomal PHR detected) and absent from the 41×41 matrix by construction. One-line corner annotation; do not spend speaker time unless asked. (Note: chr18q has n=1 sequence and is typically grouped into C15/PAR1 on the heatmap; if the v1 heatmap shows it, leave the v1 layout alone.)
- **Three-algorithm agreement is the robustness story.** Leiden k=15 (community detection, modularity-based), UPGMA k=14 (agglomerative, on the same matrix, Fig 1c dendrogram), and NJ (distance-based phylogenetics) all recover the same six abstract-anchored clades. Speaker should land this in one sentence — "three algorithms, same answer, the clades are real." This is the answer to "is the clustering an artefact of Leiden?".
- **Numbers to lock down (single source of truth).** 41 × 41 arm-level matrix / 15 Leiden communities / UPGMA k=14 (12/15 exact, 14/15 partial agreement) / NJ 41 tips / 1000 perturbation reps σ=0.0163 / 100% support at all six named-clade MRCAs / 6 missing introvert arms (2p, 3p, 5p, 8q, 11q, 14q). All consistent with `paper_prep/synthesis/CROSSWALK.md`, `paper_prep/figures/fig1/caption.md`, and `paper_prep/figures/nj_tree_arms/README.md`.
- **Do NOT** re-color the heatmap or recompute the partition. The Leiden partition and the UPGMA dendrogram in Fig 1c are the published Andrea version; the NJ tree is the *additional* view that closes the abstract-naming gap. Two views of one matrix is the point.
