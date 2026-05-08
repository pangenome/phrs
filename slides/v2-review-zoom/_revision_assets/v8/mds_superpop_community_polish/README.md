# V8 MDS Superpopulation and Community Polish

Task: `review-zoom-v8-mds-superpop-community-polish`.

## Deliverables

- `superpopulation_mds_08a_matched.png` / `.pdf`: slide 08b replacement, colored by continental superpopulation but rendered in the 12 x 10 slide 08a MDS grammar.
- `nearest_same_superpop_distance_boxplot.png` / `.pdf`: slide 08b.1 replacement using boxplots plus jitter and concise global/pairwise statistical annotations.
- `nearest_same_superpop_mds_distances.tsv`: raw per-point nearest same-superpopulation D1-D2 MDS distance table.
- `nearest_same_superpop_mds_summary.tsv`: robust per-superpopulation summaries.
- `nearest_same_superpop_pairwise_wilcoxon.tsv`: pairwise Wilcoxon rank-sum tests with BH correction.
- `nearest_same_superpop_effect_sizes.tsv`: pairwise median differences and Cliff's delta effect sizes.
- `nearest_same_superpop_global_tests.tsv`: Kruskal-Wallis global test.
- `community_mds_labeled.png` / `.pdf`: slide 09 replacement, colored by C1-C15 arm-level Leiden community with non-overlapping direct labels and leader lines.
- `community_mds_label_positions.tsv`: explicit community label anchors and label positions.
- `source_manifest.tsv`, `validation_summary.tsv`, `VALIDATION.md`, and `SLIDE_PATCH.md`: provenance, checks, and fan-in instructions.

## Metric

The slide 08b.1 metric is exactly the nearest same-superpopulation neighbor distance in the displayed MDS panel. For each sequence-level subtelomeric point, the script computes Euclidean distances in dimensions D1 and D2 to all other points from the same continental superpopulation, excludes self by setting the diagonal to `Inf`, and retains the minimum.

This is not a centroid distance, not a within-superpopulation all-pairwise distribution, and not an average against all in-group points.

## Main 2D Summary

| Superpop | MDS points | Samples | Median nearest distance | Mean nearest distance | Q1-Q3 | Q95 | Max |
|---|---:|---:|---:|---:|---:|---:|---:|
| AFR | 4,549 | 67 | 6.53e-05 | 4.11e-04 | 1.57e-05-2.10e-04 | 1.27e-03 | 6.07e-02 |
| AMR | 2,974 | 44 | 4.60e-05 | 5.02e-04 | 1.10e-05-1.61e-04 | 1.35e-03 | 7.49e-02 |
| EAS | 3,438 | 52 | 3.87e-05 | 4.33e-04 | 9.08e-06-1.52e-04 | 1.27e-03 | 4.80e-02 |
| EUR | 2,206 | 33 | 5.34e-05 | 5.91e-04 | 1.45e-05-1.82e-04 | 1.42e-03 | 8.18e-02 |
| SAS | 2,501 | 37 | 5.89e-05 | 5.88e-04 | 1.70e-05-1.94e-04 | 1.96e-03 | 1.10e-01 |

Global Kruskal-Wallis p-value: `3.51e-26`.

Largest BH-significant pairwise effects:

| Contrast | Median difference | Cliff's delta | BH q |
|---|---:|---:|---:|
| AFR vs EAS | +2.66e-05 | +0.125 | 7.21e-21 |
| EAS vs SAS | -2.02e-05 | -0.119 | 1.80e-14 |
| AFR vs AMR | +1.93e-05 | +0.093 | 2.96e-11 |
| EAS vs EUR | -1.47e-05 | -0.091 | 1.62e-08 |
| AMR vs SAS | -1.29e-05 | -0.085 | 1.03e-07 |

## Reproduction

Run from the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v8/mds_superpop_community_polish/make_mds_superpop_community_polish.R
```
