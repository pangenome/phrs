# Broom Jaccard Method Slides

Task: `review-zoom-v9-broom-jaccard-method-slides`

This directory packages Erik's broom cartoons and a slide-ready handoff for a
two-slide explanation of the PHR bundle comparison / Jaccard method. The final
deck was not edited here; integration is left to `review-zoom-v9-fanin-render`.

## Files

- `broom_1.png`: copied local asset for slide 1, showing one chromosome-end
  subtelomere bundle.
- `brooms_compare.png`: copied local asset for slide 2, showing overlap between
  two chromosome-end bundles.
- `broom_jaccard_slides.typ`: paste-ready Typst snippet using the current
  deck's `captioned-figure-slide` helper.
- `SLIDE_PATCH.md`: exact fan-in insertion instructions and slide wording.
- `source_manifest.tsv`: source paths, dimensions, and checksums for copied
  bitmap assets.

## Source Assets

The task described the source assets as repo-root inputs:

- `broom_1.png`
- `brooms_compare.png`

In this isolated worktree those files were not present at the worktree root, so
the copies in this directory were taken from the parent checkout root:

- `/moosefs/erikg/phrs/broom_1.png`
- `/moosefs/erikg/phrs/brooms_compare.png`

Checksums in `source_manifest.tsv` confirm the local copied assets are identical
to those source files.

## Method Wording

Use these two concepts in the deck.

Slide 1, bundle concept:

> All HPRCv2 haplotype PHR paths assigned to one chromosome end are treated as
> one subtelomere bundle. The unit is the population collection for that end,
> not a single reference interval.

Slide 2, comparison metric:

> For paths drawn from two chromosome-end bundles, Jaccard is intersection /
> union of the PGGB variation graph nodes they traverse. Averaged path-pair
> overlaps form the arm-level matrix for heatmaps, Leiden, and MDS; A x A can be
> below 1 when a bundle has heterogeneous haplotypes/paths.

The second sentence is intentionally compact. The deck already has the broader
workflow slide, so the cartoon should carry the explanation rather than adding a
large equation block.

## Placement Rationale

The current deck has the relevant sequence:

1. `07m`: method transition, "How sequence sharing becomes communities"
2. `07j`: "PHR Jaccard workflow"
3. `07j.1`: "PGGB graph main component: ODGI 2D layout"
4. `07j.2`: "Community assignment method"
5. `07a.1` onward: heatmaps, then MDS downstream

Insert these two broom slides immediately after `07j` and before `07j.1`. That
keeps them near the Jaccard workflow, before the PGGB graph and community
assignment method, and upstream of heatmaps, Leiden community views, and MDS.

## Accuracy Notes

- Say "chromosome end" or "subtelomere bundle"; do not reduce the bundle to one
  reference interval.
- Say "variation graph nodes" or "PGGB graph nodes".
- Keep "Jaccard = intersection / union" tied to graph nodes traversed by PHR
  paths drawn from the two chromosome-end bundles.
- Preserve the self-comparison caveat: A x A need not be 1 because the bundle
  may contain heterogeneous haplotypes/paths.
