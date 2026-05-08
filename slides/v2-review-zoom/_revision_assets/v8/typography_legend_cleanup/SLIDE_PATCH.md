# v8 Typography and Legend Cleanup Patch

Task: `review-zoom-v8-typography-legend-cleanup`

Scope: provide slide-ready assets and exact fan-in instructions only. The final
Typst deck is intentionally not edited in this task.

## Source Inspection

V7 proof PNGs inspected:

| Slide | V7 proof page | Current source in deck | Finding |
| --- | --- | --- | --- |
| 12 | `_typst/page-33.png` | `_revision_assets/hic_visual_redesign/slide_12_mouse_zygotene_trajectory_pairing.png` | Useful scatter plus stage trajectory, but right-side framing text and inset labels are too small for talk-speed reading. |
| 13b | `_typst/page-36.png` | `_typst/assets/s13_pedigree_bottom.png` | Bottom legend consumes height and is not needed because event labels are printed directly in the panel. |
| 14b | `_typst/page-39.png` | `_revision_assets/v5/gene_enrichment_figures/ranked_copy_aware_gene_signals.png` | Long row text, five-color legend, and small caption make the support-count message slow to read. |
| 14c | `_typst/page-40.png` | `_revision_assets/v5/gene_enrichment_figures/community_family_signal_map.png` | 9 x 9 tile map, tiny axis labels, tiny in-tile labels, and a small color legend are not quickly legible. |

Generated assets live in:

`slides/v2-review-zoom/_revision_assets/v8/typography_legend_cleanup/`

## Slide 12

Problem: font is about 50% too small, especially the right callout and stage
trajectory labels.

Fix: replace the v7 1800 x 1200 asset with the regenerated 3200 x 1800 asset:

- `slide12_mouse_zygotene_large_text.png`
- `slide12_mouse_zygotene_large_text.pdf`
- source values: `slide12_stage_mantel_rho.tsv`

What changed:

- Preserved the important figure content: zygotene scatter, fitted trend,
  Spearman statistic, and four-stage Mantel-rho trajectory.
- Enlarged plot title, axes, statistics, stage labels, and bouquet callout.
- Replaced the small "Readable framing" panel with a larger talk-speed
  sentence that explains why zygotene is the relevant stage.

Fan-in Typst replacement:

```typst
#figure-slide(
  "12",
  "Mouse zygotene: the bouquet-stage 3D signal",
  "../_revision_assets/v8/typography_legend_cleanup/slide12_mouse_zygotene_large_text.png",
  source: "v8/typography_legend_cleanup/make_typography_legend_cleanup.R; mouse Zuo 2021 stage tables; zygotene pair correlation and stage Mantel rho",
)
```

## Slide 13b

Problem: the bottom legend is visually detached from the detailed pedigree
panels. It is not needed for the talk version because the labels are printed
inside the panels.

Fix: delete the bottom legend using the materialized top crop:

- `slide13b_pedigree_bottom_no_unused_legend.png`
- crop helper: `crop_png_top.py`

The crop keeps the position axes and panel content, but removes the detached
legend band. Do not replace it with an explanatory legend. The source note can
say that event labels are direct in-panel annotations.

Patch snippet provided:

- `slide13b_remove_unused_legend.typ`

Replace the existing `figure-slide("13b", ...)` with:

```typst
#figure-slide(
  "13b",
  "Backup: detailed pedigree exchange events",
  "../_revision_assets/v8/typography_legend_cleanup/slide13b_pedigree_bottom_no_unused_legend.png",
  source: "v8/typography_legend_cleanup/crop_png_top.py; materialized crop of s13_pedigree_bottom.png with unused bottom legend removed; event labels are direct in-panel annotations",
)
```

## Slide 14b

Problem: font could be doubled. The v7 bar plot spends too much space on a
low-information bar display, long metric strings, and a five-color legend.

Fix: replace with:

- `slide14b_candidate_signals_talk_ready.png`
- `slide14b_candidate_signals_talk_ready.pdf`
- copied source table: `slide14b_ranked_signal_support_source.tsv`

What changed:

- Retained all eight candidate signals from the v5 support table.
- Sorted by support count for faster reading.
- Removed the color legend entirely; all bars use one direct encoding.
- Put support labels at bar ends in large text.
- Kept the statistical caveat on-slide: 116 family-community rows, 0
  BH-significant rows, best q = 0.071.

Fan-in Typst replacement:

```typst
#figure-slide(
  "14b",
  "Copy-aware candidate signals ranked by community-arm support",
  "../_revision_assets/v8/typography_legend_cleanup/slide14b_candidate_signals_talk_ready.png",
  source: "v8/typography_legend_cleanup/make_typography_legend_cleanup.R; v5 ranked_signal_support.tsv; bars are support counts, not q-values or BH-significant effects",
)
```

## Slide 14c

Problem: font is roughly 70% too small/illegible. The full 9 x 9 map and
legend require close reading, which is not suitable for a live talk.

Fix: replace with:

- `slide14c_community_family_map_talk_ready.png`
- `slide14c_community_family_map_talk_ready.pdf`
- support table: `slide14c_simplified_map_support.tsv`

What changed:

- Condensed the map to seven communities and six readable signal rows:
  `OR`, `RPL`, `SEPTIN`, `DDX11L`, `WASH/FAM138`, and `specific anchor`.
- Combined WASH and FAM138 into one row because they tell the same
  duplicon-module story at talk speed.
- Replaced the color legend with direct in-tile support labels; fill intensity
  is only secondary visual emphasis.
- Preserved the key biology: C3 broad OR/duplicon support, C5/C11/C12 recurring
  duplicated-family support, C1 DUX4L, C7 MTCO, and C15 PAR1 anchors.
- Kept the caveat: presence support only; not BH-significant enrichment.

Fan-in Typst replacement:

```typst
#figure-slide(
  "14c",
  "Community/family map separates support from statistical proof",
  "../_revision_assets/v8/typography_legend_cleanup/slide14c_community_family_map_talk_ready.png",
  source: "v8/typography_legend_cleanup/make_typography_legend_cleanup.R; v5 community_family_map_support.tsv simplified for talk legibility; direct labels replace legend",
)
```

## Asset Generation

Run from repo root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v8/typography_legend_cleanup/make_typography_legend_cleanup.R
```

The script uses copied v5 support TSVs for slides 14b/14c, the mouse
Zuo 2021 stage tables under the existing `/moosefs/guarracino/HPRCv2/...`
analysis tree for slide 12, and `crop_png_top.py` plus the current
`_typst/assets/s13_pedigree_bottom.png` source image for slide 13b.
Environment overrides:

- `MOUSE_T2T_50000BP_DIR`
- `HPRCV2_FISHER_TSV`
