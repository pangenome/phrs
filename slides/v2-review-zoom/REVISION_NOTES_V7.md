# BoG 2026 Review Zoom Deck v7 Revision Notes

Task: `review-zoom-v7-fanin-render`

Scope: `slides/v2-review-zoom`

This v7 render starts from the v6 Typst deck and fans in the completed v7 slide
feedback assets. It keeps the successful v6 additions: the charcoal PGGB graph,
the community-assignment method slide, the Dip-C/sperm validation panel, the
negative-control and community-free 3D validation summaries, and the v6
one-visual-per-page review-zoom layout.

The rendered deck is:

- `BoG_2026_review_zoom_v7.pdf`

The proof PNG export is:

- `_typst/page-01.png` through `_typst/page-46.png`

Prior reference renders remain in place:

- `BoG_2026_review_zoom.pdf`
- `BoG_2026_review_zoom_v2.pdf`
- `BoG_2026_review_zoom_v3.pdf`
- `BoG_2026_review_zoom_v4.pdf`
- `BoG_2026_review_zoom_v5.pdf`
- `BoG_2026_review_zoom_v6.pdf`

## What Changed From v6

- Replaced slide `06a` with the v7 restored PHR length histogram:
  `_revision_assets/v7/06a_length_histogram_restore/phr_length_histogram_restore.png`.
  The slide now shows 25 kb histogram bins and calls out the 500 kb analysis
  ceiling, rather than the rejected named-clade violin/grouped view.
- Updated the MDS/population transition slide `08m` so the population metric is
  explicitly nearest same-superpopulation MDS neighbor distance, not centroid
  distance, all-pairwise distance, or an average against all in-group points.
- Replaced v6 slide `08b` with an original-style superpopulation-labeled MDS
  view from
  `_revision_assets/v7/08b_nearest_same_superpop_mds/superpopulation_mds_original_style.png`.
- Added slide `08b.1` with the nearest same-superpopulation neighbor distance
  distribution from
  `_revision_assets/v7/08b_nearest_same_superpop_mds/nearest_same_superpop_distance_distribution.png`.
- Replaced slide `09` with
  `_revision_assets/v7/09_community_mds_layout/mds_community_layout.png`.
  The slide title and footer now lead with MDS, confirm the 1:1 axis rendering,
  and state that it is not PCA.
- Added slide `10m.1`, the MAPQ0/multimapper methods slide from
  `_revision_assets/v7/hic_mapq0_methods/hic_mapq0_methods_flow.svg`.
- Added slide `10m.2`, the best available direct Hi-C 3D visual:
  `_revision_assets/v7/hic_3d_plots/pngs/chm13_hic_mds_3d_coords.png`.
  The footer distinguishes this as a 3D MDS embedding of contact frequencies,
  not a physical single-cell genome reconstruction.
- Tightened the gene-enrichment transition slide `14m` with report-backed
  language from `subtelomeric_analysis_report.md` section 9 and
  `end-to-end-report/report/03_gene_enrichment.md`.
- Replaced slide `14a` with the report-backed summary SVG from
  `_revision_assets/v7/gene_enrichment_report_backed/gene_enrichment_report_backed_summary.svg`.
- Updated slide `14b` and `14c` titles/footers so the existing support-count
  figures are described as candidate/community-arm support, not q-values,
  definitive enriched classes, or BH-significant effects.
- Updated deck metadata and default footer text from v6 to v7.

## Key Caveats Preserved In The Deck

- Slide `06a`: the right edge is a measurement ceiling. The analysis searched
  terminal 500 kb windows, so the 500 kb pile-up is right-censoring, not
  evidence that longer shared sequence is absent.
- Slide `08b.1`: the distance metric is each point's nearest other
  same-superpopulation neighbor in displayed D1-D2 MDS space. Self is excluded.
  It is not a centroid, all-pairwise, hull, covariance, or average in-group
  metric.
- Slide `09`: the community view uses classical MDS coordinates from
  `hprcv2.1Mb.subtelo.full_mds.rds`, not PCA. C1-C15 labels come from arm-level
  graph-path Jaccard Leiden assignments, not 3D contact maps or gene labels.
- Slide `10m.1`: MAPQ0 reanalysis keeps one primary/randomly chosen alignment
  per multimapped read/segment. The signal supports aggregate community
  enrichment with controls and normalization, not precise pair-level claims
  inside repetitive PHRs.
- Slide `10m.2`: the CHM13 plot is a 3D MDS contact-frequency embedding. It is
  not a physical single-cell genome reconstruction.
- Slides `14m`-`14c`: the canonical HPRCv2 community-family Fisher screen has
  116 tested rows and 0 BH-significant rows. C3 OR and C7 MTCO are near-miss
  candidate presence patterns with BH q = 0.07118, not definitive enriched
  classes. GO/ORA outputs are historical context or method contrast unless
  revalidated with calibrated/permutation methods.

## Render Validation

Commands run from `slides/v2-review-zoom/_typst`:

```bash
typst compile --root .. zoom_review_deck.typ ../BoG_2026_review_zoom_v7.pdf
rm -f page-*.png
typst compile --root .. --ppi 144 zoom_review_deck.typ page-{0p}.png
```

Observed checks:

- PDF compile succeeded with Typst 0.13.1 and emitted no warnings or errors.
- PNG export produced `page-01.png` through `page-46.png`, each 1920 x 1080
  RGBA.
- PDF `/MediaBox` scan reports `[0 0 959.76 540]`, matching the 13.33 in x
  7.5 in 16:9 slide size.
- The changed v7 proof PNGs are nonblank and readable at slide scale by visual
  inspection and pixel scan.
- `git diff --check` passes.

Changed-slide PNG pixel scan:

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
