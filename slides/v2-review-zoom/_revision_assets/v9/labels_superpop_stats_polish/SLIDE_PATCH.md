# Slide Patch Recommendation

Task: `review-zoom-v9-labels-superpop-stats-polish`.

Do not edit the Typst deck in this task. These are the recommended replacements for `review-zoom-v9-fanin-render`.

## Slide 07j.2

Replace the current `community-method-stat`, `community-method-card`, and `community-assignment-method-slide` definitions with the drop-in patch file:

`slides/v2-review-zoom/_revision_assets/v9/labels_superpop_stats_polish/slide07j2_typst_patch.typ`

The patch uses this updated schematic inside the method slide:

`../_revision_assets/v9/labels_superpop_stats_polish/community_assignment_method_schematic_v9_readable.svg`

Keep the existing slide call site and source/provenance, or update only the source prefix to mention the v9 patch:

```typst
#community-assignment-method-slide(
  "07j.2",
  "Community assignment method",
  source: "v9/labels_superpop_stats_polish/slide07j2_typst_patch.typ; subtelomeric_analysis_report.md sections 5 and 6.1; HPRCv2 plot-similarity-subtelo.R; arm distance matrix and Leiden k15 assignment TSVs",
)
```

Visible wording notes:

- Use `Community definitions use graph similarity only.` for the callout title.
- Use `Biological names and 3D contact evidence are annotations interpreted after clustering.` for the callout body.
- Keep the bottom caption at 10.4 pt: `Arm-level C1-C15 partition across 41 detected-signal arms; the sequence-level 50-community partition is separate.`

## Slide 08b.1

Replace the current slide 08b.1 figure asset with:

```typst
#captioned-figure-slide(
  "08b.1",
  "Nearest same-superpopulation neighbor in MDS space",
  "Nearest same-superpopulation neighbor in MDS space",
  "../_revision_assets/v9/labels_superpop_stats_polish/nearest_same_superpop_distance_boxplot_bracketed.png",
  [
    For each subtelomeric MDS point, distance is to the nearest other point from the same continental superpopulation in displayed D1-D2 MDS space. Self is excluded. Boxes show distributions; printed values are means. KW is the Kruskal-Wallis global non-parametric test across groups; pairwise Wilcoxon is rank-sum group comparison; BH is Benjamini-Hochberg FDR correction over pairwise tests. Brackets show the five strongest BH-significant contrasts.
  ],
  source: "v9/labels_superpop_stats_polish/nearest_same_superpop_mds_distances.tsv; bracketed boxplot summary; Wilcoxon BH and Cliff delta TSVs",
)
```

The plot caption already defines KW, pairwise Wilcoxon, and BH. Keep the nearest-neighbor wording exactly: nearest other same-superpopulation point in displayed D1-D2 MDS space, self excluded. The v9.1 tweak removes jittered points, displays only boxplots, prints group means, and fixes the y-axis display to `0`-`1e-3`.

Traceability files for slide 08b.1:

- `nearest_same_superpop_mds_distances.tsv`
- `nearest_same_superpop_mds_summary.tsv`
- `nearest_same_superpop_pairwise_wilcoxon.tsv`
- `nearest_same_superpop_brackets_shown.tsv`
- `nearest_same_superpop_effect_sizes.tsv`
- `nearest_same_superpop_global_tests.tsv`

## Required Caveats

- Slide 08b.1 must retain the nearest-neighbor wording: nearest other same-superpopulation point, self excluded.
- Do not describe slide 08b.1 as centroid distance, RMS radius, all-pairwise distance, or average in-group distance.
- If the fan-in render finds the five brackets too dense at final slide scale, keep the first three contrasts and leave the full table path in the caption/source.
