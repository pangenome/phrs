# Slide 11 explicit distance-label candidate

Task: `review-zoom-v3-slide11-wb-labels`

## Decision

Slide 11 is a single-cell 3D distance ratio, not a bulk contact ratio.
The historical shorthand `W/B` on this panel means:

`mean within-community 3D distance / mean between-community 3D distance`

That is different from the bulk-contact shorthand `B/W` used for slide
10a, where the numerator is between-community contact and the denominator
is within-community contact. Do not mix those directions when revising
slide 11.

## Candidate output

Generated files in this directory:

| File | Purpose |
|---|---|
| `slide11_explicit_distance_labels_candidate.png` | Preferred slide-11 replacement/candidate visual with the numerator and denominator written out on the plot. |
| `slide11_explicit_distance_labels_candidate.pdf` | Vector copy of the same candidate for inspection or downstream conversion. |
| `slide11_explicit_distance_summary.tsv` | Small audit table of the plotted per-cell groups and counts below 1. |
| `make_slide11_explicit_distance_labels.R` | Reproducible generator using the same Dip-C and sperm per-cell TSV inputs as the v2 redesign candidate. |

The candidate visual deliberately avoids the shorthand label and instead
uses the full y-axis text:

`Within-community 3D distance / between-community 3D distance`

The on-plot direction text says that lower than 1 means same-community
arms are closer. The `S_all` group is labeled as a zero-sharing control:
seven arms with zero subtelomeric sequence sharing that should not cluster
if sequence sharing is the driver.

## Suggested deck replacement

In `slides/v2-review-zoom/_typst/zoom_review_deck.typ`, replace the slide
11 visual path with:

```typst
#figure-slide(
  "11",
  "Single-cell 3D distance: same-community arms are closer per cell",
  "../_revision_assets/v3/11_wb_labels/slide11_explicit_distance_labels_candidate.png",
  source: "v3/11_wb_labels/make_slide11_explicit_distance_labels.R; Dip-C and sperm per-cell TSVs",
)
```

Suggested caption:

`Per-cell ratio is mean within-community 3D distance divided by mean between-community 3D distance. Lower than 1 means same-community arms are closer; S_all is a zero-sharing negative control that moves the opposite way.`

One-sentence speaker line:

`This slide is a distance test: lower than 1 means same-community arms are closer in individual GM12878 and sperm nuclei, while the zero-sharing S_all control moves the opposite way.`

## Inspected sources

- `slides/v2-review-zoom/_typst/zoom_review_deck.typ`: slide 11 currently uses the v2 single-cell purpose candidate.
- `slides/v2-review-zoom/_revision_assets/hic_methods/README.md`: confirms slide 11 is single-cell 3D distance, not contact frequency, and records that values below 1 mean within-community arms are closer.
- `slides/v2-review-zoom/_revision_assets/hic_visual_redesign/README.md`: recommends keeping the Fig. 3c concept but retitling around its purpose and retaining the `S_all` negative-control result.
- `slides/v2-review-zoom/_revision_assets/hic_visual_redesign/make_hic_visual_redesign.R`: v2 candidate source for the same Dip-C and sperm per-cell data.
- `paper_prep/figures/fig3/figure_fig3.R`: original Fig. 3 panel C source logic and per-cell counts below 1.

## Validation

- Candidate visual contains no unexplained distance-ratio shorthand.
- The distance metric direction is explicit: lower than 1 means
  same-community arms are closer.
- The `S_all` group is explained as a zero-sharing negative control.
- A one-sentence speaker line is included above.
