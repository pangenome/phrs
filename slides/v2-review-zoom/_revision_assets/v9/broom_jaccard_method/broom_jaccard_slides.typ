// Paste this snippet into slides/v2-review-zoom/_typst/zoom_review_deck.typ
// after the existing #jaccard-workflow-slide("07j", ...) block and its
// pagebreak, before the current #figure-slide("07j.1", ... PGGB graph ...).
// It assumes the deck's existing captioned-figure-slide helper is in scope.

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
