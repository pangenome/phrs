# Slide 09 v3: all communities labeled, 1:1 MDS / PCoA

Task: `review-zoom-v3-slide09-label-all-communities`

## Outputs

- `mds_pcoa_all_communities_1to1.png` - square 300 dpi candidate for slide 09.
- `mds_pcoa_all_communities_1to1.pdf` - vector companion.
- `make_all_communities_1to1.R` - reproducible renderer.
- `label_positions.tsv` - exact community label positions and leader-line anchors.
- `validation_summary.tsv` - renderer-emitted checks for labels, terminology, sources, and equal-scale limits.

## Recommendation

Use `mds_pcoa_all_communities_1to1.png` as the v3 slide-09 replacement candidate.
The slide title and axis labels should say `MDS / PCoA`, not `PCA`.

## Data Provenance

The renderer uses the canonical HPRCv2 similarity outputs already cited by the
v2 review deck:

- Coordinates: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.full_mds.rds`
- Arm/community assignments: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
- Original generator: `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R`

Terminology check: this is classical MDS, equivalently PCoA for the displayed
distance projection. The off-tree generator computes:

```r
cmdscale(as.dist(jaccard_dist_df), eig = TRUE, k = 5)
```

on the Jaccard distance matrix, then stores the result as
`hprcv2.1Mb.subtelo.full_mds.rds`. The v3 candidate therefore avoids `PCA`
wording unless a true feature-matrix PCA is generated later.

Axis percentages come from the cached MDS eigenvalues:

- MDS dimension 1: 15.55%
- MDS dimension 2: 10.80%

The plotted points are the 15,668 sequence-level PHR flanks from the cached MDS
object. Each point inherits its arm-level Leiden community from the assignment
TSV through the `ChromArm` field. The label text uses the `Arms` column from the
same assignment TSV.

## Label Placement Strategy

The plot deliberately labels communities, not individual points. This keeps the
dense lower-center cluster readable while satisfying the requirement that every
community C1-C15 be identifiable.

Exact strategy:

- Compute one community anchor per Leiden community as the median x/y coordinate
  of all sequence-level PHR flanks assigned to that community.
- Place one fixed callout box for each community C1-C15 around the perimeter of
  the square panel.
- Draw a colored leader line from each median community anchor to its callout.
- Mark each median anchor with a small filled point so the leader-line origin is
  visible.
- Use the community label plus arm list in each callout, for example
  `C6 / 1q 13q 17q 19q 21q 22q`.
- Do not label individual PHR flank points.

The dense lower-center communities are explicitly fanned out across the lower
and right margins:

| Community | Arms | Callout x | Callout y | Placement note |
|---|---|---:|---:|---|
| C1 | 4q, 10q | 0.52 | 0.18 | Right margin, above the crowded lower-center labels. |
| C2 | 10p, 18p | -0.18 | -0.55 | Bottom margin, separated from C4/C7. |
| C3 | 3q, 7p, 9q, 11p, 16q, 19p | -0.50 | 0.27 | Left margin, aligned with the left MDS cloud. |
| C4 | 7q, 12q | 0.05 | -0.55 | Bottom margin for the central-lower cluster. |
| C5 | 6p, 9p, 12p, 20q | -0.42 | -0.55 | Bottom-left margin near the low-y cloud. |
| C6 | 1q, 13q, 17q, 19q, 21q, 22q | 0.46 | 0.50 | Top-right margin near the high-x/high-y cloud. |
| C7 | 13p, 14p, 15p, 21p, 22p | 0.30 | -0.55 | Bottom-right margin, separated from C2/C4. |
| C8 | 15q | -0.50 | -0.08 | Left margin near the left-lower singleton cloud. |
| C9 | 16p | 0.52 | -0.30 | Right margin for the lower-center singleton. |
| C10 | 17p | 0.52 | 0.06 | Right margin above C14/C15. |
| C11 | 1p, 5q, 6q, 8p | -0.50 | 0.47 | Top-left margin near the high-y cloud. |
| C12 | 2q, 20p | -0.50 | 0.09 | Left margin between C3 and C8. |
| C13 | 4p | 0.52 | -0.42 | Lower-right margin for the central singleton. |
| C14 | Xq, Yq | 0.52 | -0.06 | Right margin, above C15 to separate PAR2 from PAR1. |
| C15 | Xp, Yp | 0.52 | -0.18 | Right margin, below C14 to separate PAR1 from PAR2. |

The same coordinates are emitted to `label_positions.tsv` during rendering.

## Color Strategy

The six v2/slide-07 named-clade colors are preserved:

- C1 DUX4 / D4Z4: brown `#A65628`
- C2 10p-18p: orange `#FF7F00`
- C6 tight q-arm clade: purple `#984EA3`
- C7 acrocentric p-arms: green `#4DAF4A`
- C14 PAR2: blue `#377EB8`
- C15 PAR1: red `#E41A1C`

The remaining communities use muted distinct colors so that all C1-C15 labels
remain visible without changing the established v2 color semantics for the
named clades.

## 1:1 Aspect Validation

The renderer enforces equal x/y scale in three ways:

- `coord_fixed(ratio = 1)`
- identical x and y limits: `-0.68` to `0.68`
- square output dimensions: `12 x 12` inches for both PNG and PDF

This means one unit on MDS dimension 1 is the same screen distance as one unit
on MDS dimension 2. The plot is not a stretched rectangular view.

## Validation

Commands run from the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v3/09_all_communities_1to1/make_all_communities_1to1.R
git diff --check
```

Task-specific checks:

- All 15 labels C1-C15 are present in `label_positions.tsv` and in the plotted
  callout layer.
- The README documents the exact label placement strategy and the data
  provenance.
- The candidate uses `MDS / PCoA` terminology and does not label the plot as
  PCA.
- The candidate uses equal-scale `coord_fixed(ratio = 1)` with identical x/y
  limits and square output dimensions.
