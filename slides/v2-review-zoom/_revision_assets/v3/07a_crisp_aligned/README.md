# Slide 07a Crisp Aligned Candidate

## Summary

This directory contains a regenerated slide 07a candidate built from the source
41x41 arm-level Jaccard distance matrix. The candidate uses one method only:
**UPGMA average linkage** via `hclust(as.dist(D), method = "average")`.

Rows, columns, and the left tree leaves are all drawn from the same ordered arm
list. The PDF is vector-first: heatmap cells, grid lines, labels, and tree
segments are drawn directly, not copied from a raster crop. The PNG is rendered
directly from the same scene at 4800 x 2700 pixels.

No deck source was edited.

## Files

| File | Purpose |
|---|---|
| `candidate_07a_upgma_crisp_aligned.pdf` | Vector-first candidate for deck integration. |
| `candidate_07a_upgma_crisp_aligned.png` | Direct high-resolution PNG render, 4800 x 2700. |
| `arm_order_upgma.tsv` | Top-to-bottom row order and left-to-right column order used by the heatmap and UPGMA tree. |
| `order_validation.tsv` | Machine-readable validation that the tree tip order equals the heatmap row and column order. |
| `make_07a_crisp_aligned.R` | Reproducible renderer. |

## Inputs Inspected

- Existing v2 review-zoom slide 07a deck reference:
  `slides/v2-review-zoom/_typst/zoom_review_deck.typ` references
  `../_revision_assets/07a_heatmap_tree_pq/candidate_heatmap_upgma_tree_left_pq.png`.
- Existing v2 07a asset and generator:
  `slides/v2-review-zoom/_revision_assets/07a_heatmap_tree_pq/candidate_heatmap_upgma_tree_left_pq.png`,
  `slides/v2-review-zoom/_revision_assets/07a_heatmap_tree_pq/candidate_heatmap_upgma_tree_left_pq.pdf`,
  `slides/v2-review-zoom/_revision_assets/07a_heatmap_tree_pq/candidate_upgma_tree_order.tsv`,
  and `slides/v2-review-zoom/_revision_assets/07a_heatmap_tree_pq/make_candidate_heatmap.R`.
- Earlier slide 07 fanout context:
  `slides/v2/slide_07_allvsall_heatmap_nj_clades.md` and
  `slides/v2-review-zoom/_revision_assets/07b_tree_options/README.md`.
- Source matrix:
  `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`.
- Source Leiden assignments used only as TSV metadata:
  `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`.

## Diagnosis

The v2 07a issue is not an actual row-order mismatch.

The old candidate order table and this v3 order table have the same `ChromArm`
sequence. A direct check of the two arm lists returned `TRUE`.

The visible tree/heatmap mismatch comes from plotting geometry in the previous
generator:

- The left tree panel used a y-window of `0.5..n+0.5`.
- The heatmap panel used `-4.4..n+0.5` so it could also hold rotated column
  labels and a caption below the matrix.
- Because those panels had different native y-ranges, one row in the tree was
  physically taller than one row in the heatmap. With `n = 41`, the heatmap
  matrix was compressed by the extra 4.9 y-units below the cells, making the
  tree appear about 12% taller.

The fuzziness risk is raster-related. The review deck ultimately consumes a PNG,
and the older Fig 1c crop available in the deck asset tree is only 1650 x 1650
with unknown crop geometry. This v3 candidate avoids crop scaling by drawing a
native PDF and rendering the PNG directly at slide scale.

## Rendering Method

The renderer computes the UPGMA order once:

```r
hc <- hclust(as.dist(D), method = "average")
tree_order <- hc$labels[hc$order]
```

It then uses `tree_order` for all three visual components:

- left UPGMA tree leaf positions,
- heatmap rows,
- heatmap columns.

The tree and heatmap share the same physical row pitch:

```text
row_center_y = heatmap_bottom + (n - row_index + 0.5) * cell_height
tree_tip_y   = heatmap_bottom + (leaf_y - 0.5) * cell_height
```

The source matrix diagonal is retained in the displayed heatmap because the
source diagonal stores within-arm mean distance rather than self-distance.
`as.dist(D)` ignores the diagonal for UPGMA ordering.

## Labeling

The visual is labeled **UPGMA** and does not use an NJ tree. It includes crisp
vector p/q label coloring:

- p-arm labels: red
- q-arm labels: blue

Leiden k=15 is retained in `arm_order_upgma.tsv` as metadata but is not used as
the tree or heatmap ordering method.

## Order Validation

`order_validation.tsv` contains the required validation check:

```text
tree_tip_order_equals_heatmap_row_order        TRUE
tree_tip_order_equals_heatmap_column_order     TRUE
heatmap_rows_equal_heatmap_columns             TRUE
matrix_rows_equal_tree_labels                  TRUE
ordering_method                                UPGMA average linkage via hclust(method = 'average')
tree_method_label                              UPGMA
rendering_method                               PDF vector cells/tree; PNG rendered directly from the same vector scene
source_matrix_rows                             41
source_matrix_columns                          41
max_source_asymmetry                           0
raw_diagonal_min                               0.009725
raw_diagonal_max                               0.597977
```

The renderer stops with an error if any of the first four checks are not `TRUE`.

## Validation Run

Commands run from the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v3/07a_crisp_aligned/make_07a_crisp_aligned.R
file slides/v2-review-zoom/_revision_assets/v3/07a_crisp_aligned/candidate_07a_upgma_crisp_aligned.pdf slides/v2-review-zoom/_revision_assets/v3/07a_crisp_aligned/candidate_07a_upgma_crisp_aligned.png
Rscript -e 'x<-read.table("slides/v2-review-zoom/_revision_assets/07a_heatmap_tree_pq/candidate_upgma_tree_order.tsv", header=TRUE, sep="\t", check.names=FALSE); y<-read.table("slides/v2-review-zoom/_revision_assets/v3/07a_crisp_aligned/arm_order_upgma.tsv", header=TRUE, sep="\t", check.names=FALSE); cat(identical(x$ChromArm, y$ChromArm), "\n")'
```

Observed results:

```text
candidate_07a_upgma_crisp_aligned.pdf: PDF document, version 1.4
candidate_07a_upgma_crisp_aligned.png: PNG image data, 4800 x 2700, 8-bit/color RGB, non-interlaced
old 07a UPGMA arm order equals v3 arm order: TRUE
```

The candidate does not use a stretched raster crop.
