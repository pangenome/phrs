# Slide Patch Recommendation

Task: `review-zoom-v8-mds-superpop-community-polish`.

Do not edit the Typst deck in this task. These are the recommended replacements for `review-zoom-v8-fanin-render`.

## Slide 08b

Replace the current slide 08b figure asset with:

```typst
#figure-slide(
  "08b",
  "MDS colored by continental superpopulation",
  "../_revision_assets/v8/mds_superpop_community_polish/superpopulation_mds_08a_matched.png",
  source: "v8/mds_superpop_community_polish/make_mds_superpop_community_polish.R; MDS D1-D2 from hprcv2.1Mb.subtelo.full_mds.rds; 08a-matched 12x10 MDS render style; color is continental superpopulation",
)
```

Exact caption/source line to use if the slide is captioned outside `figure-slide`:

`MDS D1-D2 from hprcv2.1Mb.subtelo.full_mds.rds; rendered in the same 12 x 10 theme_bw MDS format as slide 08a; color is continental superpopulation and shape is p/q arm.`

## Slide 08b.1

Replace the current slide 08b.1 figure asset with:

```typst
#figure-slide(
  "08b.1",
  "Nearest same-superpopulation neighbor in MDS space",
  "../_revision_assets/v8/mds_superpop_community_polish/nearest_same_superpop_distance_boxplot.png",
  source: "v8/mds_superpop_community_polish/nearest_same_superpop_mds_distances.tsv; metric is nearest other same-superpopulation Euclidean distance in displayed D1-D2 MDS, self excluded; Wilcoxon BH and Cliff delta TSVs",
)
```

Exact caption/source line:

`For each subtelomeric MDS point, distance is to the nearest other point from the same continental superpopulation in displayed D1-D2 MDS space. Self is excluded; this is not a centroid metric and not an all-pairwise within-population distance. Pairwise Wilcoxon tests use BH correction; effect sizes are median differences and Cliff's delta.`

Traceability files for slide 08b.1:

- `nearest_same_superpop_mds_distances.tsv`
- `nearest_same_superpop_mds_summary.tsv`
- `nearest_same_superpop_pairwise_wilcoxon.tsv`
- `nearest_same_superpop_effect_sizes.tsv`
- `nearest_same_superpop_global_tests.tsv`

## Slide 09

Replace the current slide 09 figure asset with:

```typst
#figure-slide(
  "09",
  "MDS: all Leiden communities C1-C15 labeled",
  "../_revision_assets/v8/mds_superpop_community_polish/community_mds_labeled.png",
  source: "v8/mds_superpop_community_polish/make_mds_superpop_community_polish.R; MDS from hprcv2.1Mb.subtelo.full_mds.rds; C1-C15 from arm-level graph-path Jaccard Leiden assignments; direct labels from community_mds_label_positions.tsv; 1:1 axes; not PCA",
)
```

Exact caption/source line:

`Classical MDS on 1 - graph-path Jaccard distances; C1-C15 are arm-level Leiden communities. Direct community labels and leader lines use positions recorded in community_mds_label_positions.tsv. This is MDS, not PCA.`

## Required Caveats

- Slide 08b.1 must retain the nearest-neighbor wording: nearest other same-superpopulation point, self excluded.
- Do not describe slide 08b.1 as centroid distance, RMS radius, all-pairwise distance, or average in-group distance.
- Slide 09 must be labeled as MDS in title/footer/caption; do not call it PCA.
