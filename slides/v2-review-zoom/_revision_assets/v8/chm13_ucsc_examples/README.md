# CHM13 UCSC Browser Example Selection

Task: `review-zoom-v8-chm13-ucsc-example-selection`

This directory packages six selected CHM13/hs1 UCSC Browser examples for the
v8 review zoom fan-in. The goal is to replace v7's late backup gene-browser
slides with a small set of interpretable real UCSC panels.

## Files

- `selection_manifest.tsv`: selected examples, communities, exact coordinate
  ranges, rationales, source panel paths, copied panel paths, and slide assets.
- `SLIDE_PATCH.md`: exact v7 slide removal/insertion instructions for the v8
  fan-in renderer.
- `selected_chm13_ucsc_examples.typ`: self-contained Typst source for the six
  slide-ready pages.
- `selected_chm13_ucsc_examples.pdf`: six-page preview PDF.
- `selected_example-01.png` through `selected_example-06.png`: direct-use 16:9
  slide rasters.
- `source_panels/*.png`: local copies of the selected real UCSC browser panels.

## Selected Examples

| Slide | End | Community | Why it is selected |
| --- | --- | --- | --- |
| `14d` | chr4q | C1 | Best D4Z4/DUX4L real-UCSC example. |
| `14e` | chr3q | C3 | OR4F/f7501 multi-arm cluster context. |
| `14f` | chr15q | C8 | Clean OR4F-rich singleton endpoint. |
| `14g` | chr18p | C2 | Repeat/gene context for report-backed TAR1-rich C2 biology. |
| `14h` | chr21p | C7 | Acrocentric p-arm PHR example with MTCO/SNX18 pseudogene labels. |
| `14i` | chrXp | C15 | PAR1/SHOX example that also shows PHR is subtelomeric, not the telomere. |

## Rendering

From the repository root:

```bash
typst compile --root . --ppi 144 \
  slides/v2-review-zoom/_revision_assets/v8/chm13_ucsc_examples/selected_chm13_ucsc_examples.typ \
  slides/v2-review-zoom/_revision_assets/v8/chm13_ucsc_examples/selected_chm13_ucsc_examples.pdf

typst compile --root . --ppi 144 \
  slides/v2-review-zoom/_revision_assets/v8/chm13_ucsc_examples/selected_chm13_ucsc_examples.typ \
  slides/v2-review-zoom/_revision_assets/v8/chm13_ucsc_examples/selected_example-{0p}.png

for i in 1 2 3 4 5 6; do
  mv "slides/v2-review-zoom/_revision_assets/v8/chm13_ucsc_examples/selected_example-${i}.png" \
     "slides/v2-review-zoom/_revision_assets/v8/chm13_ucsc_examples/selected_example-0${i}.png"
done
```

## Validation

- All selected source panels are from
  `slides/chm13-phr-ucsc-browser/_assets/ucsc/panels`.
- Each slide page includes the chromosome end, community, one-line rationale,
  UCSC hs1/CHM13 source, PHR BED interval, exact browser coordinate range,
  manifest row, and copied source panel filename.
- The selected set covers D4Z4/DUX4L, OR4F-rich endpoints, C2 repeat/gene
  context, acrocentric p-arm biology, PAR1 coding genes, and a clear
  subtelomere-not-telomere example.
