# Slide 08b Superpopulation Dispersion

## Source Confirmation

- Current slide 08b asset: `slides/v2-review-zoom/_typst/assets/s08b_mds_superpop.png`.
- Exact upstream PNG: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-superpopulation.png`.
- SHA-256 for both PNGs: `3999f3eb17ec07c231e936258628e3dbbd78ebf66c369b3bb3db8c1f8ad66fb8`; the current slide asset is byte-identical to the HPRCv2 pipeline output.
- Coordinate source: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.full_mds.rds` (`15668` rows x `5` dimensions).
- Direct superpopulation label source: `/moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv` plus the hard-coded missing-sample overrides in `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R`.
- Row-level label export used as a cross-check: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv`; matches recomputed superpopulation labels for 15668 rows.
- Generation script: `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R`; the relevant call is `cmdscale(as.dist(jaccard_dist_df), eig = TRUE, k = 5)`.

Terminology: this is **classical MDS / PCoA** on a Jaccard distance matrix, not PCA on a feature matrix. The displayed slide uses MDS dimensions 1 and 2; the cached R object contains five dimensions. Axis percentages in the source plot are computed from `fit_full$eig / sum(abs(fit_full$eig)) * 100`, giving D1 = 15.55% and D2 = 10.80%.

## Metric Choice

Chosen talk metric: **RMS radius in the displayed 2D MDS panel**.

`RMS radius = sqrt(mean(||x_i - c_g||^2))`, where `x_i` is a sequence-level point in MDS dimensions 1 and 2, and `c_g` is the centroid for superpopulation `g`.

I chose RMS radius over covariance-ellipse area and convex-hull area because it is easy to say on a slide, stable under unequal sample sizes, and directly describes how far a typical point lies from its own superpopulation centre. Convex hull area was evaluated but is more sensitive to sample size and outliers; covariance ellipse area is concise statistically but harder to explain in a talk.

## Results

| Superpop | MDS points | Samples | RMS radius | Sample-bootstrap 95% CI | Mean squared radius | 68% ellipse area | Convex hull area | RMS vs non-AFR mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| AFR | 4549 | 67 | 0.342 | 0.341-0.344 | 0.117 | 0.410 | 0.423 | 1.00x |
| AMR | 2974 | 44 | 0.341 | 0.340-0.343 | 0.116 | 0.411 | 0.424 | 1.00x |
| EAS | 3438 | 52 | 0.339 | 0.338-0.341 | 0.115 | 0.407 | 0.418 | 1.00x |
| EUR | 2206 | 33 | 0.343 | 0.340-0.345 | 0.118 | 0.415 | 0.422 | 1.01x |
| SAS | 2501 | 37 | 0.340 | 0.339-0.342 | 0.116 | 0.409 | 0.420 | 1.00x |

AFR RMS radius is `0.342` versus a non-AFR mean of `0.341` (`1.00x`).
The candidate metrics agree: the displayed coordinates do **not** support a claim that AFR is unusually dispersed in the 2D MDS space. The visual impression is more likely due to point density, color salience, and every superpopulation sampling the same arm-defined clusters.

## Candidate Plots

- `superpop_dispersion_rms_radius.png` / `.pdf`: recommended slide candidate. Left panel shows the displayed MDS coordinates with AFR as the worked example for the RMS-radius definition; right panel shows the chosen metric for all superpopulations with sample-bootstrap CIs.
- `superpop_dispersion_metric_sensitivity.png` / `.pdf`: sensitivity plot comparing RMS radius, covariance ellipse area, convex hull area, and distance to the global centroid as ratios to the non-AFR mean.
- `superpop_dispersion_metrics.tsv`: full metric table used by the plots and this note.

## Interpretation

The defensible read is not "AFR spans more of the MDS space." It is: **all five superpopulations span nearly the same displayed MDS space**, consistent with slide 08b's main claim that the dominant 2D structure is arm/community structure rather than continental superpopulation.

## Limitations

- Unequal sample sizes: sample labels in this source are `AFR=67, AMR=44, EAS=52, EUR=33, SAS=37`; sequence-level point counts are `AFR=4549, AMR=2974, EAS=3438, EUR=2206, SAS=2501`. The chosen metric is less sample-size-sensitive than a hull, and the plotted CI bootstraps samples within each superpopulation, but the points are still not independent.
- Unit of analysis: the plotted and quantified unit is the **sequence-level subtelomeric flank** used in slide 08b, not an arm-level aggregate and not one point per assembly. Each sample contributes many arm/haplotype flanks, so arm composition dominates the geometry.
- MDS dimensionality: the source RDS stores five MDS dimensions, but slide 08b displays only dimensions 1 and 2. The 5D centroid metric was also computed in `superpop_dispersion_metrics.tsv` and gives the same conclusion.
- Coordinate method: because this is classical MDS / PCoA on `1 - Jaccard`, distances in the first two displayed dimensions are a low-dimensional approximation to the full Jaccard distance structure, not direct feature-space PCA distances.
- CHM13: the source script hard-codes CHM13 as EUR for superpopulation coloring; the metric follows the slide source exactly.

## Speaker-Ready Sentence

"Quantifying the same MDS panel, AFR is not more dispersed: its RMS radius is 0.342 versus 0.341 for the other superpopulations, so the 2D spread is shared across populations and mainly reflects arm-community structure."
