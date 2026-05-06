# Slide 07a Heatmap Ordering and p/q Label Review

## Bottom line

Current review-zoom slide 07a is **not a pure tree-ordered heatmap** and it is
**not ordered by the NJ tree**. The current `s07_fig1_panel_c.png` crop comes
from Fig 1 panel c. Its heatmap rows and columns are ordered primarily by the
15 arm-level Leiden communities (`C1` through `C15`), with UPGMA used only as a
secondary within-community ordering signal. The top dendrogram is UPGMA on the
same 41 x 41 arm-level Jaccard distance matrix, but the block-diagonal order is
community-first.

For a reviewer who wants "tree on the left", a cleaner alternative is feasible:
use the same 41 x 41 matrix, order both rows and columns by the UPGMA
average-linkage tree, place the dendrogram on the left, and color the row and
column labels by arm side: **p = red**, **q = blue**. A candidate render is in
this folder.

No deck source was edited.

## Files inspected

- Current zoom deck source: `slides/v2-review-zoom/_typst/zoom_review_deck.typ`
  uses `assets/s07_fig1_panel_c.png` for slide `07a`.
- Current 07a asset:
  `slides/v2-review-zoom/_typst/assets/s07_fig1_panel_c.png`, 1650 x 1650 PNG.
- Current slide 07 narrative/source:
  `slides/v2/slide_07_allvsall_heatmap_nj_clades.md`.
- Current Fig 1 source and caption:
  `paper_prep/figures/fig1/figure_fig1.R`,
  `paper_prep/figures/fig1/caption.md`, and
  `paper_prep/figures/fig1/sources.tsv`.
- Old overview deck: `slides/20260204_Subtelomics_overview_EG.pdf`
  (PDF 1.4, 3,280,684 bytes) and companion summary
  `slides/20260204_Subtelomics_overview_EG.summary.md`.
- Upstream/off-tree matrix plot script and logs:
  `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R`
  and `/moosefs/guarracino/HPRCv2/PHR_III/similarity/plot-sim-1691823.err`.
- Matrix and assignments:
  `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`,
  `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`,
  and
  `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-upgma-k14.assignments.tsv`.

## Current 07a ordering

`paper_prep/figures/fig1/figure_fig1.R` is explicit:

- It reads `hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`.
- It reads Leiden k=15 assignments from
  `hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`.
- It computes `hc_upgma <- hclust(as.dist(D), method = "average")`.
- It then orders arms as:
  `arrange(Community, upgma_pos)`.

That means:

- **Primary ordering:** Leiden community, canonical `C1` to `C15`.
- **Secondary ordering:** the UPGMA leaf position within each Leiden community.
- **Tree shown:** UPGMA average-linkage dendrogram, computed from the same
  41 x 41 matrix.
- **Not used for ordering:** NJ tree.
- **Manual component:** only the canonical community order (`C1` to `C15`);
  the arm order itself is data-derived from Leiden plus UPGMA position.

Current Fig 1 / slide 07a order:

```text
chr10_q chr4_q chr10_p chr18_p chr11_p chr16_q chr19_p chr3_q chr7_p chr9_q chr12_q chr7_q chr12_p chr20_q chr6_p chr9_p chr1_q chr13_q chr17_q chr19_q chr21_q chr22_q chr13_p chr14_p chr15_p chr21_p chr22_p chr15_q chr16_p chr17_p chr1_p chr5_q chr6_q chr8_p chr2_q chr20_p chr4_p chrX_q chrY_q chrX_p chrY_p
```

The 41-arm matrix itself is missing seven of the possible 48 chromosome arms:

```text
chr2_p chr3_p chr5_p chr8_q chr11_q chr14_q chr18_q
```

The existing slide callout only names the six "missing introvert arms"
`2p, 3p, 5p, 8q, 11q, 14q`. Do not silently add `18q` to that phrase unless
the biological wording is reviewed; the matrix excludes `chr18_q`, but the
current deck text treats the six-arm "introvert" list as a specific callout.

## Comparison to old overview slide 7

The old deck summary describes slide 7 as:

- Title: `All-vs-all - Heatmap`.
- Dendrograms on both the top and left.
- Red/blue side annotation for p versus q arms.
- Color scale around `0.2, 0.6, 0.8, 1`.
- Missing introvert arm callout: `2p, 3p, 5p, 8q, 11q, 14q`.

That matches the upstream `plot-similarity-subtelo.R` output pattern for
`hprcv2.1Mb.subtelo.heatmap.arm.clustered.pdf`:

- It builds an arm-level **similarity** matrix from the sequence-level Jaccard
  distance object.
- It plots with `pheatmap`.
- It uses `clustering_distance_rows = as.dist(1 - chrom_arm_sim_mat)`.
- It uses `clustering_distance_cols = as.dist(1 - chrom_arm_sim_mat)`.
- It uses `clustering_method = "average"`.
- It includes `annotation_row` and `annotation_col` for p/q arm side colors.

So the likely old slide 7 view was a **UPGMA-clustered Jaccard similarity
heatmap**, with rows and columns independently tree ordered by pheatmap. The
current 07a view is a **Leiden-ordered Jaccard distance heatmap**, with cyan
Leiden boxes and only a top UPGMA dendrogram. The underlying arm-level signal is
the same data family, but the display and ordering algorithm changed:

| Feature | Old overview slide 7 | Current review-zoom 07a |
| --- | --- | --- |
| Matrix scale | Jaccard similarity, visually red = more similar | Jaccard distance, Fig 1c magma scale |
| Ordering | UPGMA tree clustering in pheatmap | Leiden community first, UPGMA position second |
| Dendrogram placement | Top and left | Top only |
| Community annotation | p/q side colors; no Leiden boxes in summary | Cyan Leiden k=15 boxes and labels |
| Main claim | Clustered all-vs-all arm similarity | 15 Leiden communities, with UPGMA agreement |

The current pipeline also made the partition explicit: the March 28 run log
reports UPGMA k=14 with silhouette 0.342 and loaded Leiden k=15 at resolution
1.16 with silhouette 0.347. UPGMA and Leiden agree on the major named clades,
but they are not identical partitions. In particular, UPGMA merges the
f7501-adjacent Leiden region around C3/C8/20p and separates `2q`, while Leiden
keeps 15 communities.

## Candidate tree-left prototype

Candidate assets in this folder:

- `candidate_heatmap_upgma_tree_left_pq.png`
- `candidate_heatmap_upgma_tree_left_pq.pdf`
- `candidate_upgma_tree_order.tsv`
- `make_candidate_heatmap.R`

The candidate uses:

- Input matrix:
  `hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`.
- Ordering:
  UPGMA average linkage (`hclust(as.dist(D), method = "average")`).
- Display value:
  Jaccard similarity (`1 - distance`), so red/orange again means more similar.
- Dendrogram:
  one left-side UPGMA tree that defines both row and column order.
- p/q labels:
  p-arm labels red (`#D13F3F`), q-arm labels blue (`#245FAD`), with a small
  legend.
- Leiden status:
  thin right-side ticks mark Leiden community membership, so tree order and
  community membership are no longer conflated.

Candidate UPGMA tree order:

```text
chrX_q chrY_q chrX_p chrY_p chr4_p chr13_p chr21_p chr14_p chr15_p chr22_p chr12_q chr7_q chr10_p chr18_p chr17_p chr16_p chr9_p chr6_p chr12_p chr20_q chr2_q chr1_p chr8_p chr5_q chr6_q chr15_q chr11_p chr19_p chr3_q chr20_p chr9_q chr16_q chr7_p chr10_q chr4_q chr17_q chr13_q chr1_q chr19_q chr21_q chr22_q
```

## Recommendation for final integration

Use one of these two approaches, but label it unambiguously:

1. **Keep current Fig 1c crop / community-first story.**
   Add or preserve a caption along these lines:
   `Rows/columns ordered by Leiden k=15 community; UPGMA average-linkage dendrogram computed on the same 41 x 41 arm-level Jaccard distance matrix.`
   This keeps the Fig 1 community block story intact. If possible in final
   rendering, recolor row/column labels by p/q side or add a small p/q legend.

2. **Switch 07a to the candidate tree-left view.**
   Caption it as:
   `Rows/columns ordered by UPGMA average-linkage tree on the 41 x 41 arm-level Jaccard distance matrix; label color: p = red, q = blue; right ticks show Leiden k=15 membership.`
   This directly addresses the reviewer confusion about whether the heatmap is
   tree sorted. The cost is that the block-diagonal Leiden ordering is less
   visually strict than in Fig 1c.

Do **not** describe the heatmap order as NJ. The NJ tree is slide 07b's separate
view of the same 41-arm matrix.
