# BoG 2026 Review Zoom Deck v4 Revision Notes

Task: `review-zoom-v4-slide10a-xaxis-orientation-fix`

Scope: `slides/v2-review-zoom`

This is a focused v4 correction for slide `10a` only. It renders
`BoG_2026_review_zoom_v4.pdf` from `_typst/zoom_review_deck.typ` while
preserving the v1/v2/v3 PDFs:

- `BoG_2026_review_zoom.pdf`
- `BoG_2026_review_zoom_v2.pdf`
- `BoG_2026_review_zoom_v3.pdf`
- `BoG_2026_review_zoom_v4.pdf`

## What Changed

- Replaced slide `10a`'s deck asset with
  `_revision_assets/v4/10a_xaxis_orientation/candidate_10a_xaxis_orientation.png`.
- Added `_revision_assets/v4/10a_xaxis_orientation/make_10a_xaxis_orientation.R`,
  which regenerates the HG002 Pore-C matrix from the same source files used by
  v3.
- Added `_revision_assets/v4/10a_xaxis_orientation/orientation_audit.tsv`.
  This TSV records the displayed row order, the corrected displayed column order,
  the v3 implicit displayed column order, `v3_x_axis_mirrored = TRUE`, and the
  corrected X-axis policy.
- Added explicit first/last X-axis and Y-axis labels to the slide 10a asset:
  X runs left-to-right from `chr4_MATERNAL_q` in C1 to `chrY_PATERNAL_p` in C15;
  Y runs top-to-bottom from the same first ordered row to the same last ordered
  row.
- Updated the Typst deck metadata to v4 and changed the slide `10a` source path
  and footer provenance to the new v4 asset/audit.

No other slide visual was redesigned for v4.

## Why v3 Missed It

The v3 slide 10a generator correctly proved that the analytical row names and
column names were identical after ordering. That was necessary but not
sufficient.

The HG002 Pore-C matrix is symmetric, so a left/right visual inversion can still
look plausible: diagonal community blocks remain visually coherent, and row-name
equals column-name checks still pass.

The original manuscript Fig. 3 panel A uses base R:

```r
image(seq_len(n), seq_len(n), t(vals_norm)[, n:1], ...)
```

That transform is appropriate for `image()`, where matrix dimensions map through
the image coordinate convention. The v3 slide asset instead rendered with
`rasterImage(as.raster(t(color_matrix)[, n:1]))`. For `rasterImage()`, raster
rows are already drawn top-to-bottom and raster columns left-to-right. Reusing
the base-`image()` transform made the displayed X axis map left-to-right to
`reverse(ordered_arms)`.

The v4 generator fixes the visual policy directly:

```r
rasterImage(as.raster(color_matrix), 0.5, 0.5, n + 0.5, n + 0.5,
            interpolate = FALSE)
```

Rows are displayed top-to-bottom as `ordered_arms[1:n]`, columns are displayed
left-to-right as `ordered_arms[1:n]`, and community boxes are computed in those
displayed coordinates.

## Preserved Source Interpretation

The v4 correction preserves the slide 10a data interpretation:

- Contact matrix:
  `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_contact_matrix.tsv`
- Sequence community table:
  `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`
- Source statistic TSV:
  `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_global_test.tsv`

The generator recomputes and asserts:

- source matrix shape: 77 x 77
- source row names equal source column names
- source matrix symmetric
- all matrix arm-haplotypes map to sequence communities
- within mean = `0.0274578380817`
- between mean = `0.00152475029036`
- B/W = `0.055530602`, displayed as `0.056`
- p-value = `3.856e-85`, displayed as `3.9e-85`

## New v4 Asset Manifest

All paths below are relative to `slides/v2-review-zoom`.

| Path | Purpose |
| --- | --- |
| `_revision_assets/v4/10a_xaxis_orientation/candidate_10a_xaxis_orientation.png` | Corrected slide 10a PNG used by the deck. |
| `_revision_assets/v4/10a_xaxis_orientation/candidate_10a_xaxis_orientation.pdf` | Same corrected figure as a standalone PDF. |
| `_revision_assets/v4/10a_xaxis_orientation/make_10a_xaxis_orientation.R` | Generator with source assertions and visual orientation policy. |
| `_revision_assets/v4/10a_xaxis_orientation/orientation_audit.tsv` | Per-index displayed row order, corrected displayed X order, v3 implicit mirrored X order, and corrected policy. |
| `_revision_assets/v4/10a_xaxis_orientation/ordered_arm_haplotypes.tsv` | Display order table with source row/column indices. |
| `_revision_assets/v4/10a_xaxis_orientation/sequence_community_boxes.tsv` | Displayed-coordinate C1-C15 box extents. |
| `_revision_assets/v4/10a_xaxis_orientation/matrix_order_audit.tsv` | Source, statistic, and orientation summary checks. |
| `_revision_assets/v4/10a_xaxis_orientation/README.md` | Human-readable explanation of the v4 asset. |

## Build

Rendered from `slides/v2-review-zoom/_typst` with Typst 0.13.1:

```bash
typst compile --root .. zoom_review_deck.typ ../BoG_2026_review_zoom_v4.pdf
typst compile --root .. --ppi 144 zoom_review_deck.typ page-{0p}.png
```

The render log is `slides/v2-review-zoom/_typst/render.log`.

## Page Map Delta

The page map remains 33 pages. Page 20 / slide `10a` is the only visual asset
changed in v4:

| Page | Label | v4 visual focus |
| ---: | --- | --- |
| 20 | 10a | Corrected X-axis HG002 Pore-C community matrix with explicit orientation labels. |

## Validation

- Typst compiled `BoG_2026_review_zoom_v4.pdf` successfully.
- `pdfinfo BoG_2026_review_zoom_v4.pdf` reports 33 pages and a 959.76 x 540 pt page size.
- PNG export produced 33 page PNGs at 1920 x 1080.
- Page 20 / slide `10a` is nonblank: ImageMagick reported
  `page-20.png 1920x1080 mean=63491.8 extrema=0-65535`.
- The corrected slide 10a source asset is nonblank:
  `candidate_10a_xaxis_orientation.png 1800x1800 mean=59990.9 extrema=0-65535`.
- The v4 PNG is not a copy of the v3 PNG:
  `cmp -s v3/candidate_10a_axis_box_fix.png v4/candidate_10a_xaxis_orientation.png`
  returned exit code `1`.
- `orientation_audit.tsv` exists and states `v3_x_axis_mirrored = TRUE` plus the
  corrected X-axis policy.
- `git diff --check` passes.
- A text scan for `.wg-worktrees` under `slides/v2-review-zoom` returns no
  matches.
