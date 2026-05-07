# Review-Zoom v5 07a Tree-Then-Community Heatmaps

## Purpose

This directory contains a two-step heatmap pair for the review-zoom v5 slide 07a
story:

1. Show `07a_tree_ordered_heatmap.png` first: the 41-arm Jaccard similarity
   matrix is ordered by the current tree order, and the UPGMA tree is drawn on
   the side with rows and columns aligned.
2. Show `07a_community_ordered_heatmap.png` next: the same matrix, similarity
   scale, and heatmap palette are sorted by Leiden C1..C15 community blocks,
   with no side tree and with community bands/boxes along the axes.

This is intended to complement the existing 07a/07b story, not redesign the
deck. The first slide lets the tree reveal the structure; the second slide
summarizes that structure as the C1..C15 community grouping.

## Rendered Assets

| File | Role |
| --- | --- |
| `07a_tree_ordered_heatmap.png` | Slide-ready 4800 x 2700 raster render of the tree-ordered heatmap. |
| `07a_tree_ordered_heatmap.pdf` | Matching vector-first PDF render. |
| `07a_community_ordered_heatmap.png` | Slide-ready 4800 x 2700 raster render of the community-ordered heatmap. |
| `07a_community_ordered_heatmap.pdf` | Matching vector-first PDF render. |
| `make_07a_tree_then_community_heatmap.R` | Reproducible renderer. |
| `arm_order_tree.tsv` | Top-to-bottom and left-to-right order for the tree-ordered asset. |
| `arm_order_community.tsv` | Top-to-bottom and left-to-right order for the community-ordered asset. |
| `community_blocks.tsv` | C1..C15 community block extents in displayed coordinates. |
| `arm_inclusion_audit.tsv` | Per-arm audit over all 48 possible chr1..chr22/chrX/chrY p/q arms. |
| `source_audit.tsv` | Matrix, order, community, interval, label-size, and color-scale provenance. |
| `asset_manifest.tsv` | Rendered asset existence, expected dimensions, and byte counts. |
| `render_validation.tsv` | Machine-readable render/order checks. |

## Reproduction

Run from the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v5/07a_tree_then_community_heatmap/make_07a_tree_then_community_heatmap.R
```

The script takes optional positional arguments:

```text
1. arm-level distance matrix TSV
2. arm-level Leiden k=15 assignment TSV
3. paper_prep/figures/fig1/architecture_per_arm.tsv
4. CHM13 PHR BED
5. projected all-vs-all interval TSV
6. CHM13 chrom.sizes
7. output directory
```

Defaults are the current HPRCv2 matrix/community paths and this directory as the
output directory.

## Data And Ordering

The matrix source is:

```text
/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv
```

The cell color encodes Jaccard similarity:

```text
similarity = 1 - distance
```

Values are clipped to `[0, 1]`. Both assets use the same fixed palette and
legend limits:

```text
colorRampPalette(c("#F8FBFF", "#FEE08B", "#F46D43", "#A50026"))
```

The tree-ordered view computes the same current v3-style UPGMA order:

```r
hclust(as.dist(D), method = "average")
```

The tree leaf order is used for the side tree, heatmap rows, and heatmap
columns. This render is therefore UPGMA/tree ordered; it is not an NJ render.

The community-ordered view sorts by observed arm-level Leiden communities in
numeric order:

```text
C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15
```

Within each community block, arms retain the UPGMA leaf order. The community
view has no side tree. Community grouping is shown with a colored top band, a
matching left band, and colored diagonal block boxes.

## Label And Color Policy

Per the 2026-05-07 correction, chromosome arm labels on both x and y axes are
enlarged roughly 1.5x compared with the v3/current label size:

```text
tree view row/column label fonts:      8.9 / 8.6 pt
community view row/column label fonts: 8.6 / 8.3 pt
```

p/q label coloring is preserved:

```text
p-arm labels: #CC3B38
q-arm labels: #1F5EA8
```

The heatmap color scale was not changed while making labels larger. Community
bands use separate categorical colors and do not alter the cell palette or
legend.

## Audit Notes

`arm_inclusion_audit.tsv` covers all 48 possible chromosome arms, not only the
41 arms that enter the matrix. This makes the rendered set and missing arms
explicit:

```text
Matrix-excluded arms: chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q
```

The audit intentionally separates arm-level community assignment from called
CHM13 PHR interval status. Do not imply that every community-assigned arm has a
called interval in the repo-root CHM13 PHR BED.

Community-assigned arms that are included in the 41-arm heatmaps but have no
called interval row in `chm13.phrs.bed`:

```text
chr6_p, chr13_p, chrY_p, chrY_q
```

`chr13_p` has a non-empty CHM13#0 projected interval row in
`/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv`, but it is
not present in the repo-root `chm13.phrs.bed`. `chrY_p` and `chrY_q` have no
CHM13#0 projected all-vs-all rows in the audited table. The audit keeps those
statuses separate rather than collapsing them into a single "missing" claim.

## Validation

The renderer writes `render_validation.tsv` and stops if the order checks fail.
The current validation records:

```text
tree tip order == tree heatmap rows/columns: TRUE
community order rows == columns: TRUE
community order is C1..C15 numeric order: TRUE
matrix rows equal tree labels: TRUE
rendered PNG/PDF files exist and have nonzero byte counts: TRUE
```

`asset_manifest.tsv` records the expected dimensions: 4800 x 2700 pixels for
both PNG assets and 13.333 x 7.5 inches for both PDF assets.

No deck source was edited by this asset generation step.
