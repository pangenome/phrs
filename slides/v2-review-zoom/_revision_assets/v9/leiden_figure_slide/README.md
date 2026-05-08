# Review Zoom v9 Leiden Figure Slide

Task: `review-zoom-v9-leiden-figure-slide`

This bundle prepares a Leiden algorithm explanation slide for v9. The slide is intended to be inserted before slide `07a.2`, between the tree-ordered heatmap (`07a.1`) and the community-ordered heatmap (`07a.2`).

## Purpose

The preceding slide shows the arm-level graph-path Jaccard matrix in tree order. The following slide sorts the same matrix by Leiden community. This inserted slide explains why Leiden appears in that ordering: Leiden adds a refinement step between local node moves and aggregation, reducing the risk of poorly connected communities before the algorithm repeats at the aggregate level.

Use the official paper Figure 3 as the main visual because it shows the method sequence at a glance:

1. singleton partition,
2. node moves,
3. refinement,
4. aggregation,
5. repeated node moves and refinement on the aggregate network.

## Contents

| File | Purpose |
|---|---|
| `leiden_algorithm_fig3_official.png` | Official Springer Nature Figure 3 PNG, downloaded from the provided direct image URL. |
| `asset_manifest.tsv` | Machine-readable source, license, credit, hash, size, and intended-use manifest. |
| `LICENSE.md` | CC BY 4.0 license summary and license URL. |
| `CREDIT.md` | Required credit line and source URLs. |
| `SLIDE_PATCH.md` | Fan-in integration guidance and Typst snippet. |

## Slide Text

Recommended slide number: `07a.1b`

Recommended title: `Why Leiden for community assignment?`

Compact caption:

Leiden adds a refinement step so communities remain well connected; we use it on the arm-level graph-path Jaccard similarity network.

Recommended source/credit line:

Adapted from Traag, V. A., Waltman, L. & van Eck, N. J. From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reports 9, 5233 (2019), CC BY 4.0.

## Source And Provenance

- Official Figure 3 page: https://www.nature.com/articles/s41598-019-41695-z/figures/3
- Official image URL: https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41598-019-41695-z/MediaObjects/41598_2019_41695_Fig3_HTML.png
- Paper DOI: https://doi.org/10.1038/s41598-019-41695-z
- Local SHA-256: `e724629a09c83f753951c137648feba8b4efd962ca393bd613d70d675e72711b`
- Local dimensions: `1650x1249`
- Local size: `287945` bytes

The asset was downloaded from the official Springer Nature media URL, not from the CSDN mirror.

## Integration Notes

Do not edit the final Typst deck in this task. `review-zoom-v9-fanin-render` should apply the insertion described in `SLIDE_PATCH.md`.

The insertion point in the current deck is:

1. Keep `#figure-slide("07a.1", "Tree-ordered arm similarity heatmap", ...)`.
2. Insert the Leiden explanation slide after the `#pagebreak()` that follows `07a.1`.
3. Continue with the existing `#figure-slide("07a.2", "Same matrix grouped by Leiden community", ...)`.

## Validation

- Confirmed the deck was not edited directly.
- Confirmed the downloaded image is a PNG with dimensions `1650x1249`.
- Confirmed local SHA-256 hash: `e724629a09c83f753951c137648feba8b4efd962ca393bd613d70d675e72711b`.
- Confirmed license and credit are recorded separately for fan-in reuse.
