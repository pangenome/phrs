# Slide 08b Nearest Same-Superpopulation MDS Distance

Task: `review-zoom-v7-slide08b-nearest-superpop-mds`.

## Deliverables

- `superpopulation_mds_original_style.png` / `.pdf`: replacement source-style MDS scatter plot using the original D1-D2 coordinates, original continental superpopulation colors, p/q arm shapes, and 1:1 MDS-axis scaling.
- `nearest_same_superpop_distance_distribution.png` / `.pdf`: violin/box/jitter plot of each point's nearest same-superpopulation neighbor distance in displayed D1-D2 MDS space.
- `nearest_same_superpop_mds_distances.tsv`: row-level nearest-neighbor table, including each point's 2D nearest same-superpopulation neighbor and a secondary 5D sensitivity nearest neighbor.
- `nearest_same_superpop_mds_summary.tsv`: per-superpopulation summaries for the 2D slide metric plus the optional 5D sensitivity check.
- `source_manifest.tsv`: input paths, existence checks, and SHA-256 hashes where available.
- `VALIDATION.md`: explicit validation note for metric choice and source checks.
- `SLIDE_PATCH.md`: recommended deck insertion guidance for the v7 fan-in renderer.
- `make_nearest_same_superpop_mds.R`: reproducible generator for all files in this directory.

## Metric

The slide metric is nearest same-superpopulation neighbor distance in the displayed MDS panel.

For each sequence-level subtelomeric MDS point with a continental superpopulation label, the script computes Euclidean distance in dimensions D1 and D2 to every other point from the same superpopulation, excludes the point itself by setting the diagonal distance to `Inf`, and keeps the minimum distance and nearest-neighbor identity.

This metric is deliberately not a centroid distance, not a within-population all-pairwise distribution, and not an average against all in-group points. Pairwise distances are only used internally to identify each point's nearest eligible neighbor.

Lower values mean each subtelomere has a closer same-superpopulation neighbor in sequence-similarity MDS space.

## Source Confirmation

- Coordinate source: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.full_mds.rds` (`15668` rows x `5` exported dimensions used here; cached R object stores `5` MDS dimensions).
- Direct row-label source: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv`.
- Original sample metadata source used by the upstream plot script: `/moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv`.
- Original plotting script: `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R`; the MDS call is `cmdscale(as.dist(jaccard_dist_df), eig = TRUE, k = 5)`.
- Current slide 08b asset: `slides/v2-review-zoom/_typst/assets/s08b_mds_superpop.png`.
- Exact upstream PNG: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-superpopulation.png`.
- Upstream PNG SHA-256: `3999f3eb17ec07c231e936258628e3dbbd78ebf66c369b3bb3db8c1f8ad66fb8`.
- Current slide PNG SHA-256: `3999f3eb17ec07c231e936258628e3dbbd78ebf66c369b3bb3db8c1f8ad66fb8`; the script asserts identity when both files are present.

Terminology: this is classical MDS / PCoA on a Jaccard distance matrix, not PCA on a feature matrix. The displayed slide uses D1 and D2. Axis percentages are computed from `fit_full$eig / sum(abs(fit_full$eig)) * 100`, giving D1 = 15.55% and D2 = 10.80%.

## Main 2D Results

| Superpop | MDS points | Samples | Median nearest distance | Mean nearest distance | Q1-Q3 | Q95 | Q99 | Max | Exact zero distances |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| AFR | 4,549 | 67 | 6.53e-05 | 4.11e-04 | 1.57e-05-2.10e-04 | 1.27e-03 | 7.53e-03 | 6.07e-02 | 12 |
| AMR | 2,974 | 44 | 4.60e-05 | 5.02e-04 | 1.10e-05-1.61e-04 | 1.35e-03 | 9.91e-03 | 7.49e-02 | 12 |
| EAS | 3,438 | 52 | 3.87e-05 | 4.33e-04 | 9.08e-06-1.52e-04 | 1.27e-03 | 8.09e-03 | 4.80e-02 | 10 |
| EUR | 2,206 | 33 | 5.34e-05 | 5.91e-04 | 1.45e-05-1.82e-04 | 1.42e-03 | 1.03e-02 | 8.18e-02 | 12 |
| SAS | 2,501 | 37 | 5.89e-05 | 5.88e-04 | 1.70e-05-1.94e-04 | 1.96e-03 | 1.04e-02 | 1.10e-01 | 6 |

The nearest-neighbor values are small because many subtelomeric points have an extremely close same-superpopulation point in the displayed D1-D2 MDS panel. The violin/box/jitter figure uses a pseudo-log distance axis to show both near-zero distances and sparse outliers without changing the raw TSV values.

## 5D Sensitivity

The slide figure should use the displayed 2D MDS metric. A secondary nearest-neighbor check in all five cached MDS dimensions is included in the TSV for auditability only.

| Superpop | Median nearest distance | Mean nearest distance | Q1-Q3 | Q95 | Max |
|---|---:|---:|---:|---:|---:|
| AFR | 2.03e-04 | 1.30e-03 | 4.72e-05-6.51e-04 | 3.04e-03 | 2.39e-01 |
| AMR | 1.31e-04 | 1.34e-03 | 2.85e-05-4.57e-04 | 3.21e-03 | 1.93e-01 |
| EAS | 1.11e-04 | 1.27e-03 | 2.20e-05-4.60e-04 | 2.93e-03 | 1.56e-01 |
| EUR | 1.59e-04 | 1.46e-03 | 3.70e-05-5.42e-04 | 3.06e-03 | 2.61e-01 |
| SAS | 1.86e-04 | 1.53e-03 | 4.56e-05-6.33e-04 | 3.69e-03 | 2.24e-01 |

## Sample Sizes

- Sequence-level MDS point counts: `AFR=4549, AMR=2974, EAS=3438, EUR=2206, SAS=2501`.
- Distinct sample counts: `AFR=67, AMR=44, EAS=52, EUR=33, SAS=37`.

## Reproduction

Run from the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v7/08b_nearest_same_superpop_mds/make_nearest_same_superpop_mds.R
```

The script regenerates the TSVs, PNG/PDF figures, source manifest, README, validation note, and slide patch note in this directory.

## Validation Note

No centroid, RMS-radius, covariance-area, hull-area, all-pairwise summary, or in-group average metric is used as the slide metric. The code path for the slide distribution is the pointwise nearest-neighbor minimum in displayed D1-D2 MDS space.
