# Slide Patch Recommendation

Task: `review-zoom-v9-leiden-figure-slide`

Do not edit the Typst deck in this task. These are the recommended insertion instructions for `review-zoom-v9-fanin-render`.

## Insertion Point

File: `slides/v2-review-zoom/_typst/zoom_review_deck.typ`

Insert this slide before slide `07a.2`, between the two heatmap/community figures:

1. Keep the current `#figure-slide("07a.1", "Tree-ordered arm similarity heatmap", ...)` block.
2. Keep the `#pagebreak()` that follows slide `07a.1`.
3. Insert the Leiden explanation slide.
4. Add a `#pagebreak()`.
5. Continue with the current `#figure-slide("07a.2", "Same matrix grouped by Leiden community", ...)` block.

In the current deck, this is the gap after the `07a.1` block and before:

```typst
#figure-slide(
  "07a.2",
  "Same matrix grouped by Leiden community",
  "../_revision_assets/v5/07a_tree_then_community_heatmap/07a_community_ordered_heatmap.png",
  source: "v5/07a_tree_then_community_heatmap; same Jaccard similarity palette/scale, no side tree, C1-C15 boxes and bands",
)
```

## Typst Snippet

```typst
#captioned-figure-slide(
  "07a.1b",
  "Why Leiden for community assignment?",
  "Leiden adds refinement before aggregation",
  "../_revision_assets/v9/leiden_figure_slide/leiden_algorithm_fig3_official.png",
  [Leiden adds a refinement step so communities remain well connected; we use it on the arm-level graph-path Jaccard similarity network.],
  source: "Adapted from Traag, V. A., Waltman, L. & van Eck, N. J. From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reports 9, 5233 (2019), CC BY 4.0.",
)

#pagebreak()
```

## Required Slide Meaning

The slide should preserve these points:

- Leiden is introduced here because the next heatmap is sorted by Leiden arm-level communities.
- The visual should remain the official Figure 3 from Traag, Waltman, and van Eck, not a mirrored copy.
- The key difference to communicate is the refinement step: singleton partition -> node moves -> refinement -> aggregation -> repeat.
- The local application is our arm-level graph-path Jaccard similarity network.

## Required Caption

Use this compact caption verbatim unless the fan-in renderer must make a small typographic adjustment:

Leiden adds a refinement step so communities remain well connected; we use it on the arm-level graph-path Jaccard similarity network.

## Credit And License

Use this credit line:

Adapted from Traag, V. A., Waltman, L. & van Eck, N. J. From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reports 9, 5233 (2019), CC BY 4.0.

Traceability files:

- `asset_manifest.tsv`
- `LICENSE.md`
- `CREDIT.md`

## Avoid

- Do not use the CSDN mirror or any other non-official mirror.
- Do not present Leiden as a separate measurement; it is the community assignment step applied to the arm-level graph-path Jaccard similarity network.
- Do not insert this slide after `07a.2`; its purpose is to explain the community sorting before the community-ordered heatmap appears.
