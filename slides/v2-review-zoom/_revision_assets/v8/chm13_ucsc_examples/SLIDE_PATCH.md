# Slide Patch: CHM13 UCSC Browser Examples

Task: `review-zoom-v8-chm13-ucsc-example-selection`

## Scope

Do not edit `slides/v2-review-zoom/_typst/zoom_review_deck.typ` in this task.
This patch is for the v8 fan-in renderer.

Keep the v7 report-backed gene section through slide `14c`:

- Page 37 / slide `14m`: report-backed gene enrichment transition.
- Page 38 / slide `14a`: report-backed gene-family architecture summary.
- Page 39 / slide `14b`: copy-aware candidate signal support counts.
- Page 40 / slide `14c`: community/family support map with caveat.

Remove v7 slides `14d` through the end of the deck:

- Page 41 / slide `14d`: `Backup: DUX4/D4Z4 paired C1 chr4q/chr10q PHR view`
- Page 42 / slide `14e`: `Backup: OR4F-rich C3 chr3q PHR`
- Page 43 / slide `14f`: `Backup: OR4F pseudogene endpoint in C8 chr15q`
- Page 44 / slide `14g`: `Backup: TAR1-rich C2 chr18p PHR`
- Page 45 / slide `14h`: `Backup: C7 acrocentric p-arm panel uses one track grammar`
- Page 46 / slide `15`: `Closing: sequence sharing, 3D proximity, exchange`

Insert these six CHM13/hs1 UCSC examples after slide `14c`, in this order:

| New slide | Browser example | Community | Coordinate range | Manifest row | Asset |
| --- | --- | --- | --- | ---: | --- |
| `14d` | chr4q D4Z4/DUX4L | C1 | `chr4:193,304,946-193,574,945` | 28 | `selected_example-01.png` |
| `14e` | chr3q OR4F/f7501 cluster | C3 | `chr3:200,723,449-201,105,948` | 26 | `selected_example-02.png` |
| `14f` | chr15q OR4F endpoint | C8 | `chr15:99,565,696-99,753,195` | 11 | `selected_example-03.png` |
| `14g` | chr18p repeat/gene context | C2 | `chr18:1-397,502` | 16 | `selected_example-04.png` |
| `14h` | chr21p acrocentric p-arm | C7 | `chr21:1-727,502` | 22 | `selected_example-05.png` |
| `14i` | chrXp PAR1/SHOX | C15 | `chrX:1-750,002` | 36 | `selected_example-06.png` |

## Direct Use

The easiest v8 integration path is to insert each `selected_example-0N.png` as
a full-page figure after v7 slide `14c`. These PNGs are already 16:9 slide
rasters with:

- visible title containing chromosome end and community;
- one-line reason for inclusion;
- large real UCSC browser panel from
  `slides/chm13-phr-ucsc-browser/_assets/ucsc/panels`;
- source/data footer: UCSC hs1/CHM13 browser, PHR BED track, exact coordinate
  range, manifest row, and copied panel filename.

The paired PDF `selected_chm13_ucsc_examples.pdf` contains the same six pages.
If the fan-in renderer prefers live Typst layout rather than full-page rasters,
use `selected_chm13_ucsc_examples.typ` as the layout source and the copied panel
PNGs under `source_panels/`.

## Notes For Speaker

- The selected set is intentionally six examples, not the full 37-panel suite.
- The set covers the requested biological categories: D4Z4/DUX4L, OR4F-rich
  ends, TAR1-rich C2 repeat/gene context, acrocentric p-arm biology, PAR1, and
  the distinction between terminal telomere and internal subtelomeric PHR.
- The chr18p panel is the best real UCSC substitute for the older v3 TAR1 slide:
  the cached UCSC suite does not include a dedicated TAR1 RepeatMasker lane, so
  this slide should be described as C2 repeat/gene context with report-backed
  TAR1 richness rather than as a direct TAR1 track image.
- Keep the browser screenshots large. Do not place them inside multi-panel cards
  or shrink them below full slide width.

## Source Files

- Selection manifest: `selection_manifest.tsv`
- Slide preview source: `selected_chm13_ucsc_examples.typ`
- Slide-ready PDF: `selected_chm13_ucsc_examples.pdf`
- Slide-ready PNGs: `selected_example-01.png` through `selected_example-06.png`
- Copied browser panels: `source_panels/*.png`
