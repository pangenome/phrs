# V9 Labels and Superpopulation Stats Polish

Task: `review-zoom-v9-labels-superpop-stats-polish`.

## Deliverables

- `community_assignment_method_schematic_v9_readable.svg`: slide 07j.2 schematic with larger embedded labels and assignment-input wording.
- `slide07j2_typst_patch.typ`: drop-in Typst guidance for the slide 07j.2 macro; text sizes are about 25% larger than the current deck macro.
- `nearest_same_superpop_distance_boxplot_bracketed.png` / `.pdf`: slide 08b.1 boxplot with explicit bracket lines, stars, and BH-adjusted p-values.
- `nearest_same_superpop_mds_distances.tsv`: raw per-point nearest same-superpopulation D1-D2 MDS distance table.
- `nearest_same_superpop_mds_summary.tsv`: robust per-superpopulation summaries.
- `nearest_same_superpop_pairwise_wilcoxon.tsv`: full pairwise Wilcoxon rank-sum table with BH correction, p-value display strings, stars, and an on-slide flag.
- `nearest_same_superpop_brackets_shown.tsv`: the subset of pairwise contrasts drawn as brackets.
- `nearest_same_superpop_effect_sizes.tsv`: pairwise median differences and Cliff's delta effect sizes.
- `nearest_same_superpop_global_tests.tsv`: Kruskal-Wallis global test.
- `source_manifest.tsv`, `validation_summary.tsv`, `VALIDATION.md`, and `SLIDE_PATCH.md`: provenance, checks, and fan-in instructions.

The script also regenerates the v8-derived slide 08b and slide 09 support assets in this v9 folder for traceability, but `SLIDE_PATCH.md` only asks fan-in to change slides 07j.2 and 08b.1.

## Slide 07j.2 Label Polish

The 07j.2 guidance keeps the same provenance line while moving the visible slide wording toward short, readable statements. The replacement macro enlarges the main title, stat-card values, card labels, callout, robustness note, and bottom caption by about 25% relative to the current deck macro.

The schematic note now says that biological annotations are added after clustering. It avoids ambiguous caveat wording while preserving the methodological claim that graph similarity was the input to community assignment.

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

Bracketed on-slide pairwise contrasts:

| Contrast | Bracket label | Median difference | Cliff's delta |
|---|---:|---:|---:|
| AFR vs EAS | ***  BH p=7.2e-21 | +2.66e-05 | +0.125 |
| EAS vs SAS | ***  BH p=1.8e-14 | -2.02e-05 | -0.119 |
| AFR vs AMR | ***  BH p=3.0e-11 | +1.93e-05 | +0.093 |
| EAS vs EUR | ***  BH p=1.6e-08 | -1.47e-05 | -0.091 |
| AMR vs SAS | ***  BH p=1.0e-07 | -1.29e-05 | -0.085 |

Largest BH-significant pairwise effects in the full table:

| Contrast | Median difference | Cliff's delta | BH q |
|---|---:|---:|---:|
| AFR vs EAS | +2.66e-05 | +0.125 | 7.21e-21 |
| EAS vs SAS | -2.02e-05 | -0.119 | 1.80e-14 |
| AFR vs AMR | +1.93e-05 | +0.093 | 2.96e-11 |
| EAS vs EUR | -1.47e-05 | -0.091 | 1.62e-08 |
| AMR vs SAS | -1.29e-05 | -0.085 | 1.03e-07 |

## Statistical Wording

- KW = Kruskal-Wallis global non-parametric test across all five superpopulation groups.
- Pairwise Wilcoxon = rank-sum comparisons between two groups at a time.
- BH = Benjamini-Hochberg FDR correction over the ten pairwise tests.
- Stars use the BH-adjusted p-value: `***` for <0.001, `**` for <0.01, `*` for <0.05.

## Reproduction

Run from the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v9/labels_superpop_stats_polish/make_labels_superpop_stats_polish.R
```
