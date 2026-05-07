# BoG 2026 Review Zoom Deck v6 Manifest

Generated on 2026-05-07 UTC from `zoom_review_deck.typ`.

## Artifacts

- `../BoG_2026_review_zoom_v6.pdf`
- `zoom_review_deck.typ`
- `render.log`
- `page-01.png` through `page-43.png`
- `../REVISION_NOTES_V6.md`
- `../BoG_2026_review_zoom.pdf`, `../BoG_2026_review_zoom_v2.pdf`,
  `../BoG_2026_review_zoom_v3.pdf`, `../BoG_2026_review_zoom_v4.pdf`, and
  `../BoG_2026_review_zoom_v5.pdf` remain in place as prior reference renders.

## Page Map

| Page | Label | Visual focus |
| ---: | --- | --- |
| 01 | 01 | Title text-only focus page |
| 02 | 02m | Method transition: implicit graph over chromosome ends |
| 03 | 02 | Backup interval-tree schematic |
| 04 | 03a | IMPG workflow |
| 05 | 03b | Erdos-Renyi connectivity replacement plot |
| 06 | 04 | HPRCv2 interchrom karyogram |
| 07 | 04b | Backup manuscript Fig 1a genome-wide identity heatmap |
| 08 | 05 | Backup genome-wide count view |
| 09 | 06a | Violin/distribution PHR length view with 500 kb cap |
| 10 | 06b | Clade story matrix |
| 11 | 07m | Method transition: sequence sharing to communities |
| 12 | 07j | PHR path-to-Jaccard similarity workflow |
| 13 | 07j.1 | PGGB graph main component from ODGI 2D layout |
| 14 | 07j.2 | Community assignment method: graph-path Jaccard to arm-level Leiden calls |
| 15 | 07a.1 | Tree/UPGMA-ordered heatmap with side tree |
| 16 | 07a.2 | Same matrix ordered by Leiden C1-C15 community, no side tree |
| 17 | 07b | Rooted readable NJ tree backup |
| 18 | 07c | Unrooted NJ with acrocentric p-arm audit |
| 19 | 08m | Method transition: MDS / PCoA and pairwise population variation |
| 20 | 08a | Backup chromosome-colored MDS / PCoA |
| 21 | 08b | Within-population pairwise variation metric |
| 22 | 09 | 1:1 all-community MDS / PCoA labels |
| 23 | 10m | Method transition: explicit 3D contact validation language |
| 24 | 10a | X-axis-orientation-corrected Pore-C community matrix |
| 25 | 10b | Mantel exclusion clarity plot |
| 26 | 11 | Explicit within-community vs between-community single-cell distance plot |
| 27 | 11a | Dip-C/sperm Mantel and radial panels from existing rendered PDFs |
| 28 | 11b | S_all negative-control W/B summary |
| 29 | 11c | Community-free per-cell rho distribution |
| 30 | 12 | Mouse zygotene plus stage trajectory |
| 31 | 12b | Human sequence-similarity vs 3D-contact arm-pair analog |
| 32 | 13a | Pedigree proof, top readability crop |
| 33 | 13b | Pedigree details, bottom readability crop |
| 34 | 14m | Copy-number-aware enrichment method transition |
| 35 | 14a | Copy-number-aware method boundary figure |
| 36 | 14b | Ranked copy-aware candidate support signals |
| 37 | 14c | Community/family support map with interval-scope caveat |
| 38 | 14d | Backup DUX4/D4Z4 C1 genome-browser panel |
| 39 | 14e | Backup OR4F-rich C3 genome-browser panel |
| 40 | 14f | Backup OR4F pseudogene endpoint C8 genome-browser panel |
| 41 | 14g | Backup TAR1-rich C2 genome-browser panel |
| 42 | 14h | Backup C7 acrocentric p-arm genome-browser panel |
| 43 | 15 | Closing text-only focus page |

## Validation

- PDF compiles with Typst 0.13.1 using `--root ..`.
- PDF page size scan reports `/MediaBox [0 0 959.76 540]`, matching 13.33 in
  x 7.5 in and 16:9.
- PNG export produced 43 pages at 1920 x 1080 RGBA.
- Page 12 / slide `07j` remains the PHR path-to-Jaccard workflow slide.
- Page 13 / slide `07j.1` remains the PGGB graph main-component ODGI layout
  view and now uses
  `_revision_assets/v6/pggb_graph_black/pggb_graph_2d_black.png`.
- Page 14 / slide `07j.2` is the new arm-level community assignment method
  slide. It explicitly covers 15,668 paths, the 41 x 41 arm matrix, Leiden
  weighting, resolution 1.16, silhouette 0.347, the UPGMA comparison, and the
  no-gene-labels/no-3D-data caveat.
- Page 15 / slide `07a.1` and page 16 / slide `07a.2` remain the v5 heatmap
  pair.
- Page 24 / slide `10a` preserves the v4 corrected X-axis orientation asset:
  `_revision_assets/v4/10a_xaxis_orientation/candidate_10a_xaxis_orientation.png`.
- Page 27 / slide `11a` uses four existing rendered PDFs converted to PNG:
  GM12878 Mantel, GM12878 radial, sperm Mantel, and sperm radial.
- Page 28 / slide `11b` includes the required negative-control values:
  GM12878 `S_all` W/B = 1.106 with 0/16 cells below 1, and sperm `S_all`
  W/B = 1.397 with 1/20 cells below 1.
- Page 29 / slide `11c` includes the required community-free per-cell rho
  values: GM12878 median rho = 0.093 with 15/16 positive cells, and sperm
  median rho = 0.029 with 15/20 positive cells. The sperm arm-level pooled
  rho is explicitly labeled as a caveat.
- The v6 Dip-C section does not imply PBMC community-free analysis exists.
- Page 34 through page 37 provide the v5 copy-number-aware enrichment section,
  shifted by the three new Dip-C validation slides.
- The PGGB graph asset was re-rendered from the existing component-8 ODGI
  layout TSV with charcoal marks (`#111111`, alpha `0.30`) on white, preserving
  the v5 main-component provenance and avoiding blue graph strokes.
- RGB scan of `_revision_assets/v6/pggb_graph_black/pggb_graph_2d_black.png`
  reports 42,327 nonwhite pixels, 21,666 dark pixels, and 0 saturated-blue or
  blue-dominant pixels. A page-13 graph-region crop excluding the deck header
  and footer also reports 0 saturated-blue or blue-dominant pixels.
- Critical page PNGs are nonblank by PNG pixel scan:

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

- `_revision_assets/v6/dipc_validation/README.md` records exact sources and
  whether each asset came from an existing PDF or from a TSV-generated summary
  plot.
- `_revision_assets/v6/community_assignment_method/README.md` records the source
  anchors, matrix/assignment audits, arm-level caveats, and the hand-authored
  schematic used by slide `07j.2`.
- `_revision_assets/v6/pggb_graph_black/README.md` and `render_log.tsv` record
  the PGGB black graph source paths, retained component decision, render command,
  palette parameters, and no-SLURM provenance.
- `_revision_assets/v6/dipc_validation/source_manifest.tsv`,
  `_revision_assets/v6/dipc_validation/conversion_log.tsv`,
  `_revision_assets/v6/dipc_validation/plots/wb_negative_control_summary.tsv`,
  and `_revision_assets/v6/dipc_validation/plots/community_free_rho_summary.tsv`
  record the reproducible asset pipeline outputs.
- The slide directory contains no stale prior-agent worktree absolute paths.
- Detailed provenance for the v6 additions is recorded in
  `../REVISION_NOTES_V6.md`.
