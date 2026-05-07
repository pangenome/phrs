# Validation

Task: `review-zoom-v7-slide06a-length-histogram-restore`

## Checks Performed

- Ran `Rscript slides/v2-review-zoom/_revision_assets/v7/06a_length_histogram_restore/make_06a_length_histogram_restore.R`.
- Confirmed the generator wrote both `phr_length_histogram_restore.png` and `phr_length_histogram_restore.pdf`.
- Confirmed the generator wrote `length_distribution_summary.tsv` and `histogram_bins_25kb.tsv`.
- Confirmed the output PNG is nonblank with PNG pixel statistics: `3072x1728`, `1,085,577` non-white pixels out of `5,308,416` total pixels (`20.45%`).
- Confirmed the output PNG contains the visible ceiling wording: `analysis window ends at 500 kb; longer shared sequence is not measured`.
- Confirmed the final Typst deck was not edited.

## Visual Validation Note

The PNG is nonblank and the 500 kb ceiling is visible. The rendered plot includes an orange cap line at 500 kb, the cap-bin highlight, and a large callout using the required wording:

`analysis window ends at 500 kb; longer shared sequence is not measured`

## Reproducibility

Regenerate from the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v7/06a_length_histogram_restore/make_06a_length_histogram_restore.R
```
