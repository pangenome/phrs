# Slide Patch: Hi-C MAPQ0 Methods

Task: `review-zoom-v7-hic-mapq0-method-slide`

Do not integrate this task directly into the final deck here. This file is for
`review-zoom-v7-fanin-render`.

## Insertion Point

File: `slides/v2-review-zoom/_typst/zoom_review_deck.typ`

Insert the new slide after the existing `#pagebreak()` that follows
`#method-slide("10m", ...)` and before the
`#figure-slide("10a", "Sequence communities co-localize in 3D", ...)` block. In
the current v6 deck this is between the pagebreak at line 663 and the first
Hi-C/3D validation figure at lines 665-670.

Recommended sequence:

1. Keep slide `10m` as the high-level validation transition.
2. Add this MAPQ0 methods/background slide as `10m.1`.
3. Then continue to slide `10a` and the rest of the Hi-C/3D validation block.

## Typst Snippet

```typst
#figure-slide(
  "10m.1",
  "Making Hi-C work at subtelomeric repeats",
  "../_revision_assets/v7/hic_mapq0_methods/hic_mapq0_methods_flow.svg",
  source: "v7/hic_mapq0_methods/README.md; report anchors 05:5-11, 06:5-7, 07:34-42, 10:35-37",
)

#pagebreak()
```

## Required Slide Meaning

The slide must retain these claims:

- The original standard filtering was biased by ignoring MAPQ0 reads.
- MAPQ >=20 / removing multimappers deletes most subtelomeric signal.
- MAPQ=0 reanalysis is required for Hi-C, Pore-C, CiFi, and Dip-C validation in
  subtelomeric repeats.
- Multimappers are not duplicated across all possible placements; each
  multimapped read/segment keeps one primary/randomly chosen alignment.
- MAPQ0 adds symmetric noise, so the result supports aggregate community
  enrichment and not precise pair-level claims inside repetitive PHRs.
- Interpretation should point to flanking unique-sequence controls and
  O/E/contact normalization.

## Source Anchors For Footer Or Notes

- `end-to-end-report/report/05_hic_validation.md:5-11`
- `end-to-end-report/report/06_dipc_validation.md:5-7`
- `end-to-end-report/report/07_integrated.md:34-42`
- `end-to-end-report/report/10_limitations.md:35-37`
- `subtelomeric_analysis_report.md:1647-1653`
