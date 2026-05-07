# Slide Patch Recommendation

Task: `review-zoom-v7-slide08b-nearest-superpop-mds`.

Do not integrate this directly in this task. This note is for `review-zoom-v7-fanin-render`.

## Recommended Deck Insertion

Split the current slide 08b content into two slides.

### Slide 08b.1 - Source MDS View

- Use asset: `slides/v2-review-zoom/_revision_assets/v7/08b_nearest_same_superpop_mds/superpopulation_mds_original_style.png`.
- Suggested title: `MDS colored by continental superpopulation`.
- Suggested note: `Same D1-D2 MDS coordinates as the original slide 08b view; color is continental superpopulation and shape is p/q arm. Axes are rendered 1:1.`
- Keep this as the familiar orientation/context slide before the metric slide.

### Slide 08b.2 - Nearest Same-Superpopulation Distance

- Use asset: `slides/v2-review-zoom/_revision_assets/v7/08b_nearest_same_superpop_mds/nearest_same_superpop_distance_distribution.png`.
- Suggested title: `Nearest same-superpopulation neighbor in MDS space`.
- Suggested note: `For each subtelomeric MDS point, distance is to the nearest other point from the same continental superpopulation in displayed D1-D2 MDS space. Lower values mean each subtelomere has a closer same-superpopulation neighbor in sequence-similarity MDS space.`
- Include the sample-size labels already shown on the y-axis.
- If space is tight, keep only the title, the figure, and the lower-values interpretation note.

## Required Caveat

State or preserve this caveat in speaker notes: `Nearest same-superpopulation neighbor only; self is excluded. This is not a centroid metric, not an all-pairwise within-population distance plot, and not an average against all in-group points.`

## Assets

- `superpopulation_mds_original_style.png` / `.pdf`
- `nearest_same_superpop_distance_distribution.png` / `.pdf`
- `nearest_same_superpop_mds_distances.tsv`
- `nearest_same_superpop_mds_summary.tsv`
- `README.md`
- `VALIDATION.md`
