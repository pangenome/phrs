# BoG 2026 Review Zoom Deck v5 Manifest

Generated on 2026-05-07 UTC from `zoom_review_deck.typ`.

## Artifacts

- `../BoG_2026_review_zoom_v5.pdf`
- `zoom_review_deck.typ`
- `render.log`
- `page-01.png` through `page-38.png`
- `../REVISION_NOTES_V5.md`
- `../BoG_2026_review_zoom.pdf`, `../BoG_2026_review_zoom_v2.pdf`,
  `../BoG_2026_review_zoom_v3.pdf`, and `../BoG_2026_review_zoom_v4.pdf`
  remain in place as prior reference renders.

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
| 13 | 07a.1 | Tree/UPGMA-ordered heatmap with side tree |
| 14 | 07a.2 | Same matrix ordered by Leiden C1-C15 community, no side tree |
| 15 | 07b | Rooted readable NJ tree backup |
| 16 | 07c | Unrooted NJ with acrocentric p-arm audit |
| 17 | 08m | Method transition: MDS / PCoA and pairwise population variation |
| 18 | 08a | Backup chromosome-colored MDS / PCoA |
| 19 | 08b | Within-population pairwise variation metric |
| 20 | 09 | 1:1 all-community MDS / PCoA labels |
| 21 | 10m | Method transition: explicit 3D contact validation language |
| 22 | 10a | X-axis-orientation-corrected Pore-C community matrix |
| 23 | 10b | Mantel exclusion clarity plot |
| 24 | 11 | Explicit within-community vs between-community single-cell distance plot |
| 25 | 12 | Mouse zygotene plus stage trajectory |
| 26 | 12b | Human sequence-similarity vs 3D-contact arm-pair analog |
| 27 | 13a | Pedigree proof, top readability crop |
| 28 | 13b | Pedigree details, bottom readability crop |
| 29 | 14m | Copy-number-aware enrichment method transition |
| 30 | 14a | Copy-number-aware method boundary figure |
| 31 | 14b | Ranked copy-aware candidate support signals |
| 32 | 14c | Community/family support map with interval-scope caveat |
| 33 | 14d | Backup DUX4/D4Z4 C1 genome-browser panel |
| 34 | 14e | Backup OR4F-rich C3 genome-browser panel |
| 35 | 14f | Backup OR4F pseudogene endpoint C8 genome-browser panel |
| 36 | 14g | Backup TAR1-rich C2 genome-browser panel |
| 37 | 14h | Backup C7 acrocentric p-arm genome-browser panel |
| 38 | 15 | Closing text-only focus page |

## Validation

- PDF compiles with Typst 0.13.1 using `--root ..`.
- PDF page scan reports 38 `/Type /Page` entries.
- PDF page size scan reports `/MediaBox [0 0 959.76 540]`, matching 13.33 in
  x 7.5 in and 16:9.
- PNG export produced 38 pages at 1920 x 1080 RGBA.
- Page 12 / slide `07j` is the PHR path-to-Jaccard workflow slide. It names
  IMPG, PGGB, ODGI/Jaccard, arm/haplotype bundles, and same-arm self-bundle
  averaging, and states why `A x A` can be below 1.
- Page 13 / slide `07a.1` is the tree/UPGMA-ordered heatmap with side tree.
- Page 14 / slide `07a.2` is the same matrix in community order with no side
  tree and with community bands/boxes.
- The heatmap assets are the v5 4800 x 2700 renders with larger arm labels and
  shared Jaccard similarity palette/scale.
- Page 22 / slide `10a` preserves the v4 corrected X-axis orientation asset:
  `_revision_assets/v4/10a_xaxis_orientation/candidate_10a_xaxis_orientation.png`.
- Page 29 through page 32 provide a copy-number-aware enrichment section with
  method/caveat wording and three support-focused figure slides.
- Critical page PNGs are nonblank by PNG pixel scan:

| Page PNG | Slide | Mean RGBA value | Extrema |
| --- | --- | ---: | --- |
| `page-12.png` | `07j` | 248.88 | 10-255 |
| `page-13.png` | `07a.1` | 248.70 | 0-255 |
| `page-14.png` | `07a.2` | 248.18 | 0-255 |
| `page-22.png` | `10a` | 247.05 | 0-255 |
| `page-29.png` | `14m` | 246.53 | 10-255 |
| `page-30.png` | `14a` | 247.68 | 0-255 |
| `page-31.png` | `14b` | 245.12 | 0-255 |
| `page-32.png` | `14c` | 235.25 | 0-255 |

- `_revision_assets/v4/10a_xaxis_orientation/orientation_audit.tsv` records
  `v3_x_axis_mirrored = TRUE` and the corrected X-axis policy.
- The slide directory contains no stale prior-agent absolute worktree paths.
- Detailed provenance for the v5 additions is recorded in
  `../REVISION_NOTES_V5.md`.
