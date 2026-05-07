# BoG 2026 Review Zoom Deck v6 Revision Notes

Tasks: `review-zoom-v6-dipc-validation-slides`,
`review-zoom-v6-pggb-graph-black`, and
`review-zoom-v6-community-assignment-method-slide`

Scope: `slides/v2-review-zoom`

This v6 render starts from the v5 Typst deck, adds a concise single-cell 3D
validation section after the existing slide `11`, recolors the PGGB graph
projection for readability, and inserts one arm-level community assignment
methods slide between the PGGB graph and the tree/community heatmap sequence.
The Dip-C/sperm section reuses already-rendered PDF plots where possible,
converts them to slide PNGs, and adds two lightweight summary plots built only
from existing TSV outputs. No Dip-C, hickit, sperm scHi-C, 3D modelling, PGGB,
ODGI similarity, or community detection analysis was re-run for these slide
edits.

The rendered deck is:

- `BoG_2026_review_zoom_v6.pdf`

Prior reference renders remain in place:

- `BoG_2026_review_zoom.pdf`
- `BoG_2026_review_zoom_v2.pdf`
- `BoG_2026_review_zoom_v3.pdf`
- `BoG_2026_review_zoom_v4.pdf`
- `BoG_2026_review_zoom_v5.pdf`

## What Changed From v5

- Re-rendered slide `07j.1`, the PGGB ODGI 2D layout main-component view, with
  dark charcoal graph marks on a white background. This keeps the v5
  component-8/main-component decision and source provenance while removing the
  low-contrast blue graph styling for projection.
- Inserted slide `07j.2`, "How we assigned PHR communities", immediately after
  the PGGB graph slide and before the heatmap/tree sequence. It summarizes the
  arm-level `C1`-`C15` assignment pipeline: 15,668 paths, one `pggb -p 95`
  graph, `odgi similarity --all -P` graph-path Jaccard, the 41 x 41 arm-level
  distance matrix, Leiden clustering with silhouette-selected resolution 1.16,
  silhouette 0.347, and the UPGMA comparison.
- Inserted slide `11a`, a 2x2 panel of existing rendered GM12878 Dip-C and
  sperm scHi-C Mantel/radial plots. The source PDFs were copied into
  `_revision_assets/v6/dipc_validation/source_pdfs/` and converted with Poppler
  into `_revision_assets/v6/dipc_validation/pdf_pngs/`.
- Inserted slide `11b`, a W/B negative-control summary contrasting
  sequence-sharing C communities with non-sharing `S_all` arms in GM12878 and
  sperm. This plot is generated from existing summary/per-cell/per-community
  TSVs only.
- Inserted slide `11c`, a community-free per-cell rho distribution for GM12878
  and sperm. This plot is generated from existing community-free per-cell and
  arm-level TSVs only.
- Updated deck metadata and default footer text from v5 to v6.
- Did not add PBMC to the foreground slide section. The source report says PBMC
  community-free analysis is unavailable, so the v6 slides avoid implying that
  it exists.

## Community Assignment Method Slide

Slide `07j.2` uses:

- `_revision_assets/v6/community_assignment_method/community_assignment_method_schematic.svg`

The schematic is hand-authored and introduces no new analysis. It visualizes the
same source workflow documented in `subtelomeric_analysis_report.md` sections 5
and 6.1 and in:

- `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R`
- `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`

Slide-specific caveats are kept on the slide and in the asset README:

- No gene labels or 3D data were used to define communities.
- Biological labels and interpretations such as D4Z4, acrocentric p, PAR1/PAR2,
  f7501, and OR4F were assigned after clustering.
- The slide describes only the arm-level `C1`-`C15` partition. The separate
  sequence-level 50-community partition is a finer-grained analysis.
- The 41 x 41 matrix contains only arms with detected inter-chromosomal PHR
  signal; seven zero-signal arms were excluded.

## Slide Metrics

The v6 slide text follows `end-to-end-report/report/06_dipc_validation.md` and
`paper_prep/surveys/SURVEY_06_dipc_validation.md`.

| Slide | Metric | Value shown |
| --- | --- | --- |
| `11a` | GM12878 Dip-C W/B | 0.931, 6.9% closer |
| `11a` | GM12878 Fisher combined p | 2.4e-05 |
| `11a` | GM12878 Mantel | rho = 0.296, p = 0.002 |
| `11a` | Sperm scHi-C W/B | 0.401, 60% closer |
| `11a` | Sperm Fisher combined p | 3.9e-51 |
| `11a` | Sperm Mantel | rho = 0.202, p = 0.023 |
| `11b` | GM12878 `S_all` W/B | 1.106, 11% farther; 0/16 cells below 1 |
| `11b` | Sperm `S_all` W/B | 1.397, 40% farther; 1/20 cells below 1 |
| `11c` | GM12878 community-free per-cell rho | median rho = 0.093; 15/16 positive |
| `11c` | Sperm community-free per-cell rho | median rho = 0.029; 15/20 positive |
| `11c` | Sperm arm-level pooled caveat | rho = -0.048, p = 0.197 |

## Existing Rendered Assets Reused

All paths are relative to `slides/v2-review-zoom`.

| v6 asset | Original source | Slide use |
| --- | --- | --- |
| `_revision_assets/v6/dipc_validation/source_pdfs/gm12878_mantel_scatter.pdf` and `pdf_pngs/gm12878_mantel_scatter.png` | `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_mantel_scatter.pdf` | Slide `11a`, top-left panel. |
| `_revision_assets/v6/dipc_validation/source_pdfs/gm12878_radial_community.pdf` and `pdf_pngs/gm12878_radial_community.png` | `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_radial_community.pdf` | Slide `11a`, top-right panel. |
| `_revision_assets/v6/dipc_validation/source_pdfs/sperm_all20_mantel_scatter.pdf` and `pdf_pngs/sperm_all20_mantel_scatter.png` | `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_mantel_scatter.pdf` | Slide `11a`, bottom-left panel. |
| `_revision_assets/v6/dipc_validation/source_pdfs/sperm_all20_radial_community.pdf` and `pdf_pngs/sperm_all20_radial_community.png` | `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_radial_community.pdf` | Slide `11a`, bottom-right panel. |
| `_revision_assets/v6/dipc_validation/source_pdfs/sperm_all20_by_arm_type_arm.pdf` and `pdf_pngs/sperm_all20_by_arm_type_arm.png` | `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/overlay_plots/sperm_all20.by_arm-type.arm.pdf` | Staged/converted for provenance; not placed in the concise v6 section. |
| `_revision_assets/v6/dipc_validation/source_pdfs/sperm_all20_by_arm_type_per_cell.pdf` and `pdf_pngs/sperm_all20_by_arm_type_per_cell.png` | `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/overlay_plots/sperm_all20.by_arm-type.per-cell.pdf` | Staged/converted for provenance; not placed in the concise v6 section. |

Conversion command, run by
`_revision_assets/v6/dipc_validation/prepare_dipc_validation_assets.sh`:

```bash
guix shell ghostscript poppler -- pdftoppm -r 220 -png -singlefile <source.pdf> <asset-prefix>
```

The conversion hashes and sizes are recorded in:

- `_revision_assets/v6/dipc_validation/source_manifest.tsv`
- `_revision_assets/v6/dipc_validation/conversion_log.tsv`

## PGGB Graph Readability Recolor

Slide `07j.1` now uses:

- `_revision_assets/v6/pggb_graph_black/pggb_graph_2d_black.png`

This asset was re-rendered from the existing ODGI layout TSV rather than by
running a new graph extraction, layout, ODGI draw, or gfalook job. It preserves
the v5 provenance:

- ODGI layout TSV component `8`
- Main graph component, not the full graph
- 727,156 layout nodes
- Existing ODGI X/Y coordinates plotted as Y/X for the 16:9 slide frame

The reproducible render pipeline is:

- `_revision_assets/v6/pggb_graph_black/render_pggb_graph_black.sh`
- `_revision_assets/v6/pggb_graph_black/render_pggb_layout_component8_black.R`
- `_revision_assets/v6/pggb_graph_black/render_log.tsv`

Palette parameters are recorded in `render_log.tsv`: charcoal graph marks
(`#111111`, alpha `0.30`) on a white background with a neutral border
(`#b8c0cc`). No SLURM job was used because this lightweight render reads only
the existing 41 MB layout TSV.

## Generated Summary Assets

The generated summary plots are reproducible with:

```bash
cd slides/v2-review-zoom/_revision_assets/v6/dipc_validation
./prepare_dipc_validation_assets.sh
```

Generated outputs:

| Asset | Source TSVs |
| --- | --- |
| `_revision_assets/v6/dipc_validation/plots/wb_negative_control.png` | `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_summary.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_per_cell.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_per_community_per_cell.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_summary.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_per_cell.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_per_community_per_cell.tsv`. |
| `_revision_assets/v6/dipc_validation/plots/community_free_rho_distribution.png` | `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_community_free_per_cell.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_community_free_arm.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_community_free_per_cell.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_community_free_arm.tsv`. |

Summary TSVs written by the plot script:

- `_revision_assets/v6/dipc_validation/plots/wb_negative_control_summary.tsv`
- `_revision_assets/v6/dipc_validation/plots/community_free_rho_summary.tsv`

## Render Validation

Commands run from `slides/v2-review-zoom/_typst`:

```bash
typst compile --root .. zoom_review_deck.typ ../BoG_2026_review_zoom_v6.pdf
typst compile --root .. --ppi 144 zoom_review_deck.typ page-{0p}.png
```

Observed checks:

- PDF compiles successfully with Typst 0.13.1.
- PDF compiles successfully with Typst 0.13.1.
- Page size scan reports `/MediaBox [0 0 959.76 540]`.
- Page PNG export produced `page-01.png` through `page-43.png`, each
  1920 x 1080 RGBA.
- Page `13` / slide `07j.1` uses the v6 black PGGB graph asset and is nonblank.
- Page `14` / slide `07j.2` is the new community assignment methods slide and
  is nonblank/readable at 1920 x 1080.
- Page `15` and page `16` are the shifted tree/community heatmap context pages.
- Dip-C/sperm validation slides are shifted to `page-27.png` through
  `page-29.png`; context exports include `page-26.png` and `page-30.png`.
- `git diff --check` passes.
- No stale prior-agent worktree absolute paths are embedded under
  `slides/v2-review-zoom`.

Nonblank PNG scan for the new slides and nearby context:

| Page PNG | Slide | Mean RGBA value | Extrema |
| --- | --- | ---: | --- |
| `page-13.png` | `07j.1` | 252.46 | 0-255 |
| `page-14.png` | `07j.2` | 247.17 | 10-255 |
| `page-15.png` | `07a.1` | 248.70 | 0-255 |
| `page-16.png` | `07a.2` | 248.18 | 0-255 |
| `page-26.png` | `11` | 250.54 | 0-255 |
| `page-27.png` | `11a` | 246.51 | 0-255 |
| `page-28.png` | `11b` | 246.67 | 0-255 |
| `page-29.png` | `11c` | 247.18 | 0-255 |
| `page-30.png` | `12` | 250.62 | 0-255 |
