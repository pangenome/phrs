# Slide Patch Recommendation

Task: `review-zoom-v9-broom-jaccard-method-slides`

Do not edit the Typst deck in this task. These are the recommended insertion
instructions for `review-zoom-v9-fanin-render`.

## Insertion Point

File: `slides/v2-review-zoom/_typst/zoom_review_deck.typ`

Insert two slides immediately after the existing `#jaccard-workflow-slide("07j",
...)` block and the `#pagebreak()` that follows it. This places the cartoons
between:

```typst
#jaccard-workflow-slide(
  "07j",
  "PHR Jaccard workflow",
  source: "HPRCv2 PHR similarity workflow: IMPG interval calls, PGGB graph over all paths, ODGI path Jaccard, arm-bundle averaging",
)

#pagebreak()
```

and the current PGGB graph slide:

```typst
#figure-slide(
  "07j.1",
  "PGGB graph main component: ODGI 2D layout",
  "../_revision_assets/v6/pggb_graph_black/pggb_graph_2d_black.png",
  source: "v6/pggb_graph_black; derived from v5/pggb_graph_odgi; ODGI layout TSV component 8 main component; 727,156 nodes; rotated 16:9; charcoal-on-white recolor",
)
```

Recommended slide IDs are `07j.a` and `07j.b` to avoid renumbering existing
`07j.1` and `07j.2`. If the final deck requires strictly numeric IDs, renumber
the later method slides consistently during fan-in.

## Typst Snippet

The same snippet is saved as `broom_jaccard_slides.typ`.

```typst
#captioned-figure-slide(
  "07j.a",
  "Chromosome-end bundle",
  "A chromosome-end bundle is a population collection",
  "../_revision_assets/v9/broom_jaccard_method/broom_1.png",
  [All HPRCv2 haplotype PHR paths assigned to one chromosome end are treated as one subtelomere bundle. The unit is the population collection for that end, not a single reference interval.],
  source: "Erik Garrison broom cartoon copied from repo-root broom_1.png; method wording: subtelomeric_analysis_report.md sections 5 and 6.1",
)

#pagebreak()

#captioned-figure-slide(
  "07j.b",
  "Bundle Jaccard",
  "Compare two chromosome-end bundles by graph-node overlap",
  "../_revision_assets/v9/broom_jaccard_method/brooms_compare.png",
  [For paths drawn from two chromosome-end bundles, Jaccard is intersection / union of the PGGB variation graph nodes they traverse. Averaged path-pair overlaps form the arm-level matrix for heatmaps, Leiden, and MDS; A x A can be below 1 when a bundle has heterogeneous haplotypes/paths.],
  source: "Erik Garrison broom cartoon copied from repo-root brooms_compare.png; method wording: subtelomeric_analysis_report.md sections 5 and 6.1 and current slide 07j",
)

#pagebreak()
```

## Slide 1

Title/label:

- Header label: `Chromosome-end bundle`
- Visible title: `A chromosome-end bundle is a population collection`

Asset:

- `../_revision_assets/v9/broom_jaccard_method/broom_1.png`

Caption:

All HPRCv2 haplotype PHR paths assigned to one chromosome end are treated as one
subtelomere bundle. The unit is the population collection for that end, not a
single reference interval.

Source line:

Erik Garrison broom cartoon copied from repo-root broom_1.png; method wording:
subtelomeric_analysis_report.md sections 5 and 6.1

## Slide 2

Title/label:

- Header label: `Bundle Jaccard`
- Visible title: `Compare two chromosome-end bundles by graph-node overlap`

Asset:

- `../_revision_assets/v9/broom_jaccard_method/brooms_compare.png`

Caption:

For paths drawn from two chromosome-end bundles, Jaccard is intersection / union
of the PGGB variation graph nodes they traverse. Averaged path-pair overlaps
form the arm-level matrix for heatmaps, Leiden, and MDS; A x A can be below 1
when a bundle has heterogeneous haplotypes/paths.

Source line:

Erik Garrison broom cartoon copied from repo-root brooms_compare.png; method
wording: subtelomeric_analysis_report.md sections 5 and 6.1 and current slide
07j

## Required Meaning

- Use two slides, not one crowded slide.
- Slide 1 introduces a chromosome-end subtelomere bundle as the whole HPRCv2
  haplotype PHR path collection assigned to one chromosome end.
- Slide 1 must not imply the bundle is one reference interval.
- Slide 2 explains that the comparison is graph-node Jaccard on PGGB/ODGI
  variation graph nodes, summarized into the arm-level similarity/distance
  matrix.
- Slide 2 must preserve the compact self-comparison caveat: A x A can be below
  1 when a chromosome-end bundle contains heterogeneous haplotypes/paths.
- The two slides should remain upstream of the heatmaps, Leiden community views,
  and MDS, because those views consume the same arm-level Jaccard matrix.

## Avoid

- Do not merge the two cartoons into a single crowded slide.
- Do not describe the cartoon as one reference interval or one reference
  haplotype.
- Do not say self-comparisons are always 1.
- Do not add a large equation block; the required formula can stay as
  "intersection / union" in the caption.
