# BoG 2026 Review Zoom Deck v7 Manifest

Generated on 2026-05-07 UTC from `zoom_review_deck.typ`.

## Artifacts

- `../BoG_2026_review_zoom_v7.pdf`
- `zoom_review_deck.typ`
- `render.log`
- `page-01.png` through `page-46.png`
- `../REVISION_NOTES_V7.md`
- Prior reference PDFs `../BoG_2026_review_zoom.pdf`,
  `../BoG_2026_review_zoom_v2.pdf`,
  `../BoG_2026_review_zoom_v3.pdf`,
  `../BoG_2026_review_zoom_v4.pdf`,
  `../BoG_2026_review_zoom_v5.pdf`, and
  `../BoG_2026_review_zoom_v6.pdf` remain in place.

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
| 09 | 06a | v7 PHR length histogram with 500 kb analysis-window ceiling |
| 10 | 06b | Clade story matrix |
| 11 | 07m | Method transition: sequence sharing to communities |
| 12 | 07j | PHR path-to-Jaccard similarity workflow |
| 13 | 07j.1 | PGGB graph main component from ODGI 2D layout |
| 14 | 07j.2 | Community assignment method: graph-path Jaccard to arm-level Leiden calls |
| 15 | 07a.1 | Tree/UPGMA-ordered heatmap with side tree |
| 16 | 07a.2 | Same matrix ordered by Leiden C1-C15 community, no side tree |
| 17 | 07b | Rooted readable NJ tree backup |
| 18 | 07c | Unrooted NJ with acrocentric p-arm audit |
| 19 | 08m | Method transition: MDS and nearest same-superpopulation neighbor spread |
| 20 | 08a | Backup chromosome-colored MDS / PCoA |
| 21 | 08b | Original-style superpopulation-labeled MDS |
| 22 | 08b.1 | Nearest same-superpopulation MDS-neighbor distance distribution |
| 23 | 09 | v7 all-community MDS labels, 1:1 axes, not PCA |
| 24 | 10m | Method transition: explicit 3D contact validation language |
| 25 | 10m.1 | MAPQ0/multimapper methods for repetitive subtelomeric contacts |
| 26 | 10m.2 | CHM13 Hi-C 3D MDS contact-space view |
| 27 | 10a | X-axis-orientation-corrected Pore-C community matrix |
| 28 | 10b | Mantel exclusion clarity plot |
| 29 | 11 | Explicit within-community vs between-community single-cell distance plot |
| 30 | 11a | Dip-C/sperm Mantel and radial panels from existing rendered PDFs |
| 31 | 11b | S_all negative-control W/B summary |
| 32 | 11c | Community-free per-cell rho distribution |
| 33 | 12 | Mouse zygotene plus stage trajectory |
| 34 | 12b | Human sequence-similarity vs 3D-contact arm-pair analog |
| 35 | 13a | Pedigree proof, top readability crop |
| 36 | 13b | Pedigree details, bottom readability crop |
| 37 | 14m | Report-backed gene enrichment transition with conservative statistics |
| 38 | 14a | Report-backed gene-family architecture summary SVG |
| 39 | 14b | Candidate signals ranked by community-arm support counts |
| 40 | 14c | Community/family support map with statistical-proof caveat |
| 41 | 14d | Backup DUX4/D4Z4 C1 genome-browser panel |
| 42 | 14e | Backup OR4F-rich C3 genome-browser panel |
| 43 | 14f | Backup OR4F pseudogene endpoint C8 genome-browser panel |
| 44 | 14g | Backup TAR1-rich C2 genome-browser panel |
| 45 | 14h | Backup C7 acrocentric p-arm genome-browser panel |
| 46 | 15 | Closing text-only focus page |

## Validation

- PDF compiles with Typst 0.13.1 using `--root ..`.
- Typst emitted no warnings or errors during PDF or PNG export.
- PDF `/MediaBox` scan reports `[0 0 959.76 540]`, matching 13.33 in x
  7.5 in and 16:9.
- PNG export produced 46 pages at 1920 x 1080 RGBA.
- Page 09 / slide `06a` uses the v7 histogram and visibly labels the 500 kb
  analysis-window ceiling.
- Page 21 / slide `08b` uses the original-style superpopulation MDS.
- Page 22 / slide `08b.1` uses nearest same-superpopulation neighbor distance
  in displayed D1-D2 MDS space. It is explicitly not centroid/all-pairwise.
- Page 23 / slide `09` is MDS, not PCA; all C1-C15 communities are labeled and
  the footer records 1:1 axes and graph-path Jaccard Leiden assignments.
- Page 25 / slide `10m.1` records the MAPQ0/multimapper handling caveat.
- Page 26 / slide `10m.2` is the CHM13 Hi-C 3D MDS contact-space view and is
  labeled as not a physical single-cell reconstruction.
- Page 37 / slide `14m` and page 38 / slide `14a` cite the report-backed
  gene-enrichment interpretation and caveat the canonical Fisher screen:
  116 tested rows, 0 BH-significant rows, with C3 OR and C7 MTCO as candidate
  presence patterns.
- Page 39 / slide `14b` states that bars are support counts, not q-values or
  BH-significant effects.
- Page 40 / slide `14c` frames the map as report-backed presence patterns, not
  definitive enriched classes.
- Changed pages are nonblank by pixel scan and readable by visual inspection:

| Page PNG | Slide | Mean RGBA | Nonwhite pixels | Darkish pixels |
| --- | --- | ---: | ---: | ---: |
| `page-09.png` | `06a` | 238.76 | 21.74% | 18.91% |
| `page-21.png` | `08b` | 250.88 | 8.84% | 6.09% |
| `page-22.png` | `08b.1` | 248.92 | 11.86% | 9.16% |
| `page-23.png` | `09` | 250.38 | 10.04% | 7.38% |
| `page-25.png` | `10m.1` | 243.93 | 37.54% | 26.68% |
| `page-26.png` | `10m.2` | 250.93 | 19.93% | 18.86% |
| `page-37.png` | `14m` | 245.24 | 46.27% | 14.66% |
| `page-38.png` | `14a` | 234.09 | 53.55% | 39.94% |
| `page-39.png` | `14b` | 245.07 | 15.07% | 12.73% |
| `page-40.png` | `14c` | 235.22 | 55.12% | 27.36% |

- Detailed provenance and caveats for v7 are recorded in
  `../REVISION_NOTES_V7.md`.
