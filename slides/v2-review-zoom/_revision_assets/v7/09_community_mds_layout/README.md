# Slide 09 v7: community-labeled MDS layout

Task: `review-zoom-v7-slide09-community-mds-layout`

This directory contains the v7 replacement candidate for review-zoom slide 09.
It rebuilds the current all-community slide as an MDS-first plot with the same
visual grammar as slide 08a: a 12 x 10 inch ggplot frame, right-side legend,
large readable axes, `theme_bw()`, and a fixed 1:1 coordinate scale.

## Outputs

Run the renderer from the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v7/09_community_mds_layout/make_community_mds_layout.R
```

The script writes:

- `mds_community_layout.png` - 3600 x 3000 px PNG, slide-ready replacement image.
- `mds_community_layout.pdf` - vector companion.
- `label_positions.tsv` - emitted C1-C15 anchor, callout, arm, and point-count table.
- `validation_summary.tsv` - machine-readable checks for MDS source, 1:1 axes, all labels, and non-PCA status.

## Source Data

- MDS coordinates: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.full_mds.rds`
- Community assignments: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
- Original MDS generator: `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R`
- Layout reference: `slides/v2-review-zoom/_typst/assets/s08a_mds_chrom.png`
- Previous all-community label artifact: `slides/v2-review-zoom/_revision_assets/v3/09_all_communities_1to1/label_positions.tsv`

The task input mentioned `slides/v2-review-zoom/_revision_assets/v3/09_mds_community_labels`,
but that path is not present in this worktree. The existing v3 all-community
label output is `v3/09_all_communities_1to1`, and this v7 renderer reuses its
fixed callout centers while recomputing anchors and counts from the current MDS
and community assignment sources.

## MDS And Community Meaning

The coordinate source is a cached classical MDS result from:

```r
cmdscale(as.dist(jaccard_dist_df), eig = TRUE, k = 5)
```

The displayed dimensions are:

- Dimension 1: 15.55%
- Dimension 2: 10.80%

Each plotted point is one sequence-level HPRCv2 subtelomeric flank from the
cached MDS object. Points inherit C1-C15 from their chromosome arm through the
arm-level Leiden assignment table. These communities were assigned from
arm-level graph-path Jaccard distances. They were not assigned from 3D contact
data or gene labels.

## Label Strategy

All 15 communities are labeled directly on the plot.

- The point cloud is colored by stable C1-C15 community.
- Arm shape follows slide 08a grammar: p-arm points are circles and q-arm points are triangles.
- One median MDS anchor is computed per community.
- One fixed callout is placed per community using the v3 all-community label positions.
- Colored leader lines connect community anchors to labels.
- Label text uses the stable community ID plus the arm list from the assignment TSV.

This intentionally labels communities rather than individual sequence-level
flanks so the dense clouds remain readable at slide scale.

## Color Strategy

The six named-clade colors used by earlier review slides are preserved:

| Community | Meaning | Color |
|---|---|---|
| C1 | DUX4 / D4Z4, 4q/10q | `#A65628` |
| C2 | 10p-18p | `#FF7F00` |
| C6 | tight q-arm clade | `#984EA3` |
| C7 | acrocentric p-arms | `#4DAF4A` |
| C14 | PAR2, Xq/Yq | `#377EB8` |
| C15 | PAR1, Xp/Yp | `#E41A1C` |

The remaining C labels use the stable muted colors from the v3 all-community
asset so every community has a distinct point, line, label, and legend color.

## Validation

The renderer enforces and records the main task checks:

- MDS is used from `hprcv2.1Mb.subtelo.full_mds.rds`; no PCA coordinates are used.
- Axes are fixed at a 1:1 coordinate scale through `coord_fixed(ratio = 1)`.
- X and Y limits are identical: `-0.68,0.68`.
- All C1-C15 communities are present in `label_positions.tsv`.
- The plot includes the support text that C1-C15 communities come from arm-level graph-path Jaccard distances, not 3D/gene labels.

Validation command:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v7/09_community_mds_layout/make_community_mds_layout.R
git diff --check
```
