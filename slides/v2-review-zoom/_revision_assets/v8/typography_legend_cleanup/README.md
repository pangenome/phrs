# v8 Typography and Legend Cleanup Assets

This folder contains slide-ready assets and fan-in patch guidance for
`review-zoom-v8-typography-legend-cleanup`.

The final deck is not edited here. Integrate through
`review-zoom-v8-fanin-render`.

## Contents

| File | Purpose |
| --- | --- |
| `SLIDE_PATCH.md` | Exact before/after notes and Typst replacements for slides 12, 13b, 14b, and 14c. |
| `LEGIBILITY_CHECKLIST.md` | Talk-speed legibility checklist. |
| `make_typography_legend_cleanup.R` | Regenerates the slide 12, 13b, 14b, and 14c assets from source tables/images. |
| `crop_png_top.py` | Minimal deterministic PNG cropper used for slide 13b because no pedigree regeneration recipe is available. |
| `slide13b_remove_unused_legend.typ` | Typst replacement snippet to use the no-legend slide 13b crop. |
| `asset_manifest.tsv` | Generated asset manifest. |
| `slide12_mouse_zygotene_large_text.png` / `.pdf` | Slide 12 replacement with larger labels and preserved zygotene scatter plus stage trajectory. |
| `slide12_stage_mantel_rho.tsv` | Stage trajectory values used in the slide 12 asset. |
| `slide13b_pedigree_bottom_no_unused_legend.png` | Slide 13b materialized crop with the unused bottom legend removed. |
| `slide14b_candidate_signals_talk_ready.png` / `.pdf` | Slide 14b replacement with direct labels and no legend. |
| `slide14b_ranked_signal_support_source.tsv` | Copied v5 support source for slide 14b. |
| `slide14c_community_family_map_talk_ready.png` / `.pdf` | Slide 14c condensed map with direct labels and no legend. |
| `slide14c_simplified_map_support.tsv` | Condensed support table used for slide 14c. |

## Regeneration

From repo root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v8/typography_legend_cleanup/make_typography_legend_cleanup.R
```

The script reads copied v5 gene-enrichment support tables from
`slides/v2-review-zoom/_revision_assets/v5/gene_enrichment_figures/` and the
existing mouse Zuo 2021 tables under `/moosefs/guarracino/HPRCv2/...`. It also
uses `crop_png_top.py` to materialize the slide 13b no-legend crop from the
current deck asset.

Optional environment overrides:

- `MOUSE_T2T_50000BP_DIR`
- `HPRCV2_FISHER_TSV`
