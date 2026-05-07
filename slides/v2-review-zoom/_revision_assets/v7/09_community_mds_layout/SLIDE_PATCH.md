# Slide 09 Replacement Patch

Task: `review-zoom-v7-slide09-community-mds-layout`

Do not apply this patch in this task. It is for
`review-zoom-v7-fanin-render` when assembling the v7 deck.

## Replace Slide 09 Image

In `slides/v2-review-zoom/_typst/zoom_review_deck.typ`, replace the current
slide 09 asset reference:

```typst
#figure-slide(
  "09",
  "MDS / PCoA: all Leiden communities C1-C15 labeled",
  "../_revision_assets/v3/09_all_communities_1to1/mds_pcoa_all_communities_1to1.png",
  source: "v3/09_all_communities_1to1/make_all_communities_1to1.R; validation_summary.tsv confirms all C1-C15 labels",
)
```

with:

```typst
#figure-slide(
  "09",
  "MDS: all Leiden communities C1-C15 labeled",
  "../_revision_assets/v7/09_community_mds_layout/mds_community_layout.png",
  source: "v7/09_community_mds_layout/make_community_mds_layout.R; MDS coordinates from hprcv2.1Mb.subtelo.full_mds.rds; C1-C15 from arm-level graph-path Jaccard Leiden assignments; validation_summary.tsv confirms 1:1 axes and all labels",
)
```

## Caption Or Speaker Note

Recommended concise support note:

> C1-C15 are arm-level Leiden communities assigned from graph-path Jaccard
> distances. They are not assigned from 3D contact maps or gene labels.

## Integration Notes

- The replacement should lead with MDS terminology.
- The asset uses the same 12 x 10 inch / 3600 x 3000 plot frame as slide 08a.
- The plot itself uses `coord_fixed(ratio = 1)` with identical x and y limits.
- All communities C1-C15 are directly labeled on the plot and colored by the
  stable community palette.
- Do not use the old `assets/s09_pca_communities.png` image for slide 09.
