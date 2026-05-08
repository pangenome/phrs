# BoG 2026 Review Zoom v9 Revision Notes

Task: `review-zoom-v9-fanin-render`

Rendered draft:

- `slides/v2-review-zoom/BoG_2026_review_zoom_v9.pdf`
- proof PNGs: `slides/v2-review-zoom/_typst/page-01.png` through `page-52.png`
- Typst source: `slides/v2-review-zoom/_typst/zoom_review_deck.typ`

## Fan-In Summary

v9 integrates the final polish handoffs under
`slides/v2-review-zoom/_revision_assets/v9/`.

| Patch handoff | Integrated deck changes |
| --- | --- |
| `06a_q_axis_kbp/SLIDE_PATCH.md` | Slide 06a now uses the 10 kbp per-arm chromosome-end heatstrip with q arms flipped so the panels read p telomere to q telomere. The source/footer wording uses `>500 kbp not measured` and all slide-06a length units are kbp. |
| `broom_jaccard_method/SLIDE_PATCH.md` | Added two method cartoon slides after the PHR Jaccard workflow: one defining a chromosome-end bundle as the full HPRCv2 haplotype PHR path collection for that end, and one explaining graph-node Jaccard as intersection / union over PGGB nodes with the same-arm self-bundle caveat. |
| `labels_superpop_stats_polish/SLIDE_PATCH.md` | Slide 07j.2 uses the v9 readable community assignment schematic, removes the old label-exclusion wording, and increases method text size. Slide 08b.1 uses the bracketed nearest same-superpopulation MDS distance boxplot with stars and BH p-values, and the deck caption defines the metric, KW, pairwise Wilcoxon, and BH in plain language. |
| `leiden_figure_slide/SLIDE_PATCH.md` | Inserted the official Springer Nature Leiden Figure 3 slide between the tree-ordered heatmap and the Leiden community-ordered heatmap. The slide credits Traag, Waltman, and van Eck, Scientific Reports 9, 5233 (2019), CC BY 4.0. |
| `slide10m2_better_3d_viz/SLIDE_PATCH.md` | Slide 10m.2 now uses the clearer CHM13 PHR/subtelomeric Hi-C contact-space MDS visualization. The embedded caption states that it is bulk Hi-C contact-space MDS, not a physical single-cell genome reconstruction. |

## Changed Slide Map

| Proof page | Slide | v9 change |
| --- | --- | --- |
| `page-09.png` | 06a | q-arm orientation flipped; 10 kbp bins and `>500 kbp not measured` wording. |
| `page-13.png` | 07j.a | Added chromosome-end bundle broom cartoon slide. |
| `page-14.png` | 07j.b | Added bundle Jaccard broom cartoon slide with intersection / union and A x A caveat. |
| `page-16.png` | 07j.2 | Larger community assignment method slide with v9 schematic and revised annotation wording. |
| `page-18.png` | 07a.1b | Added official Leiden Figure 3 slide before slide 07a.2. |
| `page-25.png` | 08b.1 | Nearest same-superpopulation MDS distance boxplot now includes KW, pairwise Wilcoxon, BH correction, stars, p-values, and brackets. |
| `page-29.png` | 10m.2 | Replaced whole-arm CHM13 MDS with CHM13 PHR/subtelomeric contact-space MDS and explicit non-physical-reconstruction caption. |

Slides after the new method/Leiden insertions are shifted relative to v8, so
many proof PNG filenames from `page-13.png` onward changed even when their
visible slide content did not.

## Provenance Notes

All newly introduced data figures have an on-slide source/provenance line.

The Leiden figure comes from the official Springer Nature article image bundle,
with local provenance in `v9/leiden_figure_slide/asset_manifest.tsv`,
`LICENSE.md`, and `CREDIT.md`. The deck footer keeps the CC BY 4.0 credit on
the slide.

Slide 08b.1 keeps the required metric wording: nearest other point from the
same continental superpopulation in displayed D1-D2 MDS space, with self
excluded. It is not a centroid distance, RMS radius, all-pairwise distance, or
average in-group distance.

Slide 10m.2 uses the v9 contact-space replacement figure generated from
`chm13_hic.dist_matrix.tsv` and `chm13_subtelomeric_regions.bed`. The figure
caption distinguishes bulk Hi-C contact-space MDS from a physical single-cell
3DG reconstruction.

## Validation

Commands run from `slides/v2-review-zoom/_typst`:

```bash
typst compile --root .. zoom_review_deck.typ ../BoG_2026_review_zoom_v9.pdf
typst compile --root .. --ppi 144 zoom_review_deck.typ page-{0p}.png
```

Render results:

- PDF compile completed without Typst errors.
- Proof export produced 52 PNG pages at 1920 x 1080.
- All 52 proof PNG pages were decoded and checked for 1920 x 1080 geometry and
  nonzero RGB range. The changed content pages were also checked for nonzero
  pixel variance: `page-09.png`, `page-13.png`, `page-14.png`, `page-16.png`,
  `page-18.png`, `page-25.png`, and `page-29.png` all returned `OK`.
- Visual inspection covered the requested changed proof pages: 06a, 07j.a,
  07j.b, 07j.2, 07a.1b, 08b.1, and 10m.2.
- Source grep found no remaining legacy label-exclusion wording or old
  500-kbp-limit wording in `zoom_review_deck.typ`.
- `git diff --check` passed.

No Cargo validation was run because this is a documentation/render artifact
task and no Rust code was changed.
