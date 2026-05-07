# Slide 08b Within-Population Pairwise Variation

## Metric

The metric is the average Euclidean distance between every unordered pair of sequence-level PHR-flank points from the same superpopulation in the displayed MDS dimensions 1 and 2.

This replaces the v2 dispersion framing: the candidate recommendation is based on same-superpopulation point-to-point distances only.

## Source Confirmation

- Current slide 08b asset: `slides/v2-review-zoom/_typst/assets/s08b_mds_superpop.png`.
- Exact upstream PNG: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-superpopulation.png`.
- Upstream PNG SHA-256: `3999f3eb17ec07c231e936258628e3dbbd78ebf66c369b3bb3db8c1f8ad66fb8`.
- Current slide PNG SHA-256: `3999f3eb17ec07c231e936258628e3dbbd78ebf66c369b3bb3db8c1f8ad66fb8`; the current slide asset is byte-identical to the HPRCv2 pipeline output.
- Coordinate source: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.full_mds.rds` (`15668` rows x `5` dimensions).
- Superpopulation label source: `/moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv` plus hard-coded missing-sample overrides in `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R`.
- Row-level label export cross-check: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv`; matches recomputed superpopulation labels for 15668 rows.
- Generation script: `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R`; the relevant call is `cmdscale(as.dist(jaccard_dist_df), eig = TRUE, k = 5)`.

Terminology: this is **classical MDS / PCoA** on a Jaccard distance matrix, not PCA on a feature matrix. The displayed slide uses dimensions 1 and 2; the cached R object contains five dimensions. Axis percentages in the source plot are computed from `fit_full$eig / sum(abs(fit_full$eig)) * 100`, giving D1 = 15.55% and D2 = 10.80%.

## Main 2D Results

| Superpop | MDS points | Samples | Same-pop pairs | Mean distance | Sample-bootstrap 95% CI | Median distance | Q1-Q3 | IQR |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| AFR | 4,549 | 67 | 10,344,426 | 0.422 | 0.417-0.425 | 0.445 | 0.241-0.583 | 0.342 |
| AMR | 2,974 | 44 | 4,420,851 | 0.420 | 0.416-0.423 | 0.447 | 0.235-0.583 | 0.348 |
| EAS | 3,438 | 52 | 5,908,203 | 0.417 | 0.414-0.421 | 0.442 | 0.231-0.581 | 0.351 |
| EUR | 2,206 | 33 | 2,432,115 | 0.422 | 0.417-0.426 | 0.450 | 0.235-0.585 | 0.350 |
| SAS | 2,501 | 37 | 3,126,250 | 0.418 | 0.414-0.422 | 0.442 | 0.232-0.583 | 0.350 |

The 2D point-to-point summaries are close across superpopulations, so the displayed panel does not support a claim that one superpopulation has uniquely larger within-population spread.

## 5D Sensitivity

| Superpop | Mean distance | Sample-bootstrap 95% CI | Median distance | Q1-Q3 | IQR |
|---|---:|---:|---:|---:|---:|
| AFR | 0.600 | 0.596-0.603 | 0.652 | 0.490-0.755 | 0.265 |
| AMR | 0.600 | 0.596-0.604 | 0.657 | 0.493-0.757 | 0.263 |
| EAS | 0.596 | 0.592-0.600 | 0.655 | 0.485-0.752 | 0.266 |
| EUR | 0.600 | 0.595-0.605 | 0.660 | 0.497-0.757 | 0.260 |
| SAS | 0.598 | 0.594-0.603 | 0.656 | 0.489-0.755 | 0.265 |

The five-dimensional check uses all cached MDS dimensions and preserves the same qualitative read as the displayed two-dimensional panel.

## Outputs

- `within_pop_pairwise_2d_distribution.png` / `.pdf`: main candidate plot for slide 08b; violin/box shows sampled same-superpopulation pair distances, diamond shows the exact mean, and the line shows the sample-bootstrap CI.
- `within_pop_pairwise_2d_vs_5d_sensitivity.png` / `.pdf`: sensitivity check comparing the same metric in displayed 2D coordinates and cached 5D coordinates.
- `within_pop_pairwise_summary.tsv`: exact per-superpopulation summaries for 2D and 5D distances, including n points, n samples, number of pairs, mean, median, quartiles, IQR, and CI columns.
- `within_pop_pairwise_distance_sample.tsv`: sampled pair distances used for the violin/box display; this avoids writing tens of millions of pair rows while preserving an auditable plotted distribution sample.
- `make_within_pop_pairwise.R`: reproducible generator for the outputs in this directory.

## Limitations

- Unequal sample sizes: sample labels in this source are `AFR=67, AMR=44, EAS=52, EUR=33, SAS=37`; sequence-level point counts are `AFR=4549, AMR=2974, EAS=3438, EUR=2206, SAS=2501`. The number of same-population pairs scales with point count, so the summary table reports both n points and n pairs.
- Non-independent PHR flanks: the plotted and quantified unit is a sequence-level subtelomeric flank, not one independent individual. Each sample contributes many arm/haplotype flanks, and pair distances share points.
- Pairwise distances are also non-independent because each point appears in many pairs; the CI is a descriptive sample-bootstrap interval and should not be read as an independent-pair hypothesis test.
- MDS dimensionality: the main candidate uses only the displayed D1-D2 panel. The 5D sensitivity check is included because the cached MDS object stores five coordinates, but both are low-dimensional summaries of the full Jaccard distance structure.
- Coordinate method: because this is classical MDS / PCoA on `1 - Jaccard`, the reported values are distances in MDS coordinate space, not direct feature-space PCA distances.
- CHM13: the source script hard-codes CHM13 as EUR for superpopulation coloring; this analysis follows the slide source exactly.

## Speaker-Ready Sentence

"Using same-superpopulation point-to-point distances in the displayed MDS panel, the populations have similar within-population spread, so the visual structure is better read as shared arm/community geometry than as one population spreading more than the rest."
