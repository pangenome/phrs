# Validation Note

Task: `review-zoom-v7-slide09-community-mds-layout`

## Checks

- **MDS is used:** `make_community_mds_layout.R` reads
  `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.full_mds.rds`,
  the cached classical MDS / cmdscale result from graph-path Jaccard distances.
  The rendered slide asset does not use PCA coordinates.
- **Axes are 1:1:** the plot uses `coord_fixed(ratio = 1)` with identical x
  and y limits of `-0.68,0.68`, and this is recorded in
  `validation_summary.tsv`.
- **All communities are labeled:** `label_positions.tsv` contains exactly
  `C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15`, and the plot draws one
  direct callout label for each community.
- **Community source is explicit:** the subtitle and patch note state that
  C1-C15 were assigned from arm-level graph-path Jaccard distances, not from
  3D contact maps or gene labels.
- **Slide 08a layout grammar is matched:** the PNG is 3600 x 3000 px from a
  12 x 10 inch ggplot render with a right-side legend, readable axes, and
  comparable point shapes for p/q arms.

## Commands Run

```bash
Rscript slides/v2-review-zoom/_revision_assets/v7/09_community_mds_layout/make_community_mds_layout.R
file slides/v2-review-zoom/_revision_assets/v7/09_community_mds_layout/mds_community_layout.png
Rscript -e 'v <- read.delim("slides/v2-review-zoom/_revision_assets/v7/09_community_mds_layout/validation_summary.tsv"); print(v, row.names = FALSE)'
git diff --check
```
