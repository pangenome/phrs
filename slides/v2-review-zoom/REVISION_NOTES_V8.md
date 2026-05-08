# BoG 2026 Review Zoom v8 Revision Notes

Task: `review-zoom-v8-fanin-render`

Rendered draft:

- `slides/v2-review-zoom/BoG_2026_review_zoom_v8.pdf`
- proof PNGs: `slides/v2-review-zoom/_typst/page-01.png` through `page-49.png`
- Typst source: `slides/v2-review-zoom/_typst/zoom_review_deck.typ`

## Fan-In Summary

v8 integrates all slide patch handoffs under
`slides/v2-review-zoom/_revision_assets/v8/`.

| Patch handoff | Integrated deck changes |
| --- | --- |
| `06a_length_alternatives/SLIDE_PATCH.md` | Slide 06a now uses the 10 kb per-arm chromosome-end heatstrip. The slide title and figure text keep the terminal 500 kb discovery ceiling explicit, avoiding the v7 25 kb mega-bin look. |
| `mds_superpop_community_polish/SLIDE_PATCH.md` | Slide 08b now uses the 08a-matched superpopulation MDS render. Slide 08b.1 now uses the nearest same-superpopulation boxplot with self excluded, robust summaries, BH-corrected Wilcoxon notes, and Cliff delta effect-size notes. Slide 09 now uses the labeled community MDS with non-overlapping manual label positions and is explicitly called MDS, not PCA. |
| `hic_dipc_clarity_split/SLIDE_PATCH.md` | Added a captioned figure helper. Slide 10m.2 now states that the plotted CHM13 regions are whole p-arm/centromere/q-arm intervals, not PHRs, not terminal 500 kb windows, and not 1 Mb flanks/windows; it is contact-space MDS, not physical reconstruction. The old dense slide 11a was replaced by four readable slides: GM12878 Dip-C proximity, GM12878 Dip-C radial, sperm scHi-C proximity, and sperm scHi-C radial. Slides 11b and 11c use reduced-text v8 assets and explanatory captions. |
| `typography_legend_cleanup/SLIDE_PATCH.md` | Slide 12 uses the enlarged mouse zygotene asset. Slide 13b uses the no-unused-legend crop. Slides 14b and 14c use the talk-ready candidate-signal and support-map assets with larger type and direct labels. |
| `chm13_ucsc_examples/SLIDE_PATCH.md` | The old v7 slides 14d through the closing slide were removed and replaced by six selected CHM13/hs1 UCSC browser examples. The six raster slides include title, rationale, browser panel, coordinate range, manifest row, and source footer. |

## Changed Slide Map

| Proof page | Slide | v8 change |
| --- | --- | --- |
| `page-09.png` | 06a | 10 kb per-arm/end PHR length heatstrip; explicit 500 kb ceiling. |
| `page-21.png` | 08b | 08a-matched superpopulation MDS scale/style. |
| `page-22.png` | 08b.1 | Nearest same-superpopulation D1-D2 MDS distance boxplots; self excluded; robust summary/stat notes. |
| `page-23.png` | 09 | Community MDS labels with leader lines and reduced overlap. |
| `page-26.png` | 10m.2 | CHM13 whole-arm Hi-C contact-space MDS captioned with region definition and not-physical-reconstruction caveat. |
| `page-30.png` | 11a.1 | GM12878 Dip-C proximity panel; W/B and proximity direction defined. |
| `page-31.png` | 11a.2 | GM12878 Dip-C radial panel with source/provenance line. |
| `page-32.png` | 11a.3 | Sperm scHi-C proximity panel; W/B and proximity direction defined. |
| `page-33.png` | 11a.4 | Sperm scHi-C radial panel with source/provenance line. |
| `page-34.png` | 11b | Reduced-text W/B negative-control plot; W/B direction explained. |
| `page-35.png` | 11c | Reduced-text community-free rho distribution; positive rho direction explained. |
| `page-36.png` | 12 | Enlarged mouse zygotene/stage-trajectory text. |
| `page-39.png` | 13b | Bottom unused legend removed; event labels remain direct in-panel annotations. |
| `page-42.png` | 14b | Larger copy-aware candidate-signal support bars; no unused legend. |
| `page-43.png` | 14c | Talk-ready community/family support map with larger direct labels. |
| `page-44.png` | CHM13 example 01 | chr4q C1 D4Z4/DUX4L UCSC example. |
| `page-45.png` | CHM13 example 02 | chr3q C3 OR4F/f7501 UCSC example. |
| `page-46.png` | CHM13 example 03 | chr15q C8 OR4F endpoint UCSC example. |
| `page-47.png` | CHM13 example 04 | chr18p C2 repeat/gene-context UCSC example. |
| `page-48.png` | CHM13 example 05 | chr21p C7 acrocentric p-arm UCSC example. |
| `page-49.png` | CHM13 example 06 | chrXp C15 PAR1/SHOX UCSC example. |

## Provenance Notes

All newly introduced data figures have an on-slide source/provenance line.
For the split 10m.2 and 11a-11c slides, the deck footers point to
`v8/hic_dipc_clarity_split/Source_manifest.tsv`, which records dataset,
technology, sample/cell counts, source TSV/PDF/PNG files, region definitions,
and the distance/proximity convention.

For the MDS slides, the deck footers point to
`v8/mds_superpop_community_polish/make_mds_superpop_community_polish.R` and
the generated nearest-neighbor/statistics TSVs. Slide 08b.1 keeps the required
wording: nearest other same-superpopulation point in displayed D1-D2 MDS space,
self excluded; it is not a centroid, radius, or all-pairwise metric.

For the CHM13 browser examples, each slide raster contains its own UCSC hs1 /
CHM13 footer with PHR BED coordinates, browser window, selection manifest row,
and copied panel filename.

## Validation

Commands run from this worktree:

```bash
typst compile --root .. zoom_review_deck.typ ../BoG_2026_review_zoom_v8.pdf
typst compile --root .. --ppi 144 zoom_review_deck.typ page-{0p}.png
```

Render results:

- PDF compile completed without Typst errors.
- Proof export produced 49 PNG pages at 1920 x 1080.
- Changed proof PNGs were checked with a PNG decoder/stat pass; every changed
  proof page had nonzero pixel variance and was marked `OK`.
- Targeted visual inspection covered the requested changed proof pages:
  06a, 08b, 08b.1, 09, 10m.2, 11a.1-11a.4, 11b, 11c, 12, 13b, 14b, 14c, and
  the selected CHM13 UCSC examples.
- `git diff --check` passed.

No Cargo validation was run because this is a documentation/render artifact
task and no `Cargo.toml` was found at repository depth 3.
