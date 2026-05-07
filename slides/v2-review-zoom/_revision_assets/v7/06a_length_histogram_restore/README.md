# Slide 06a Length Histogram Restore

Task: `review-zoom-v7-slide06a-length-histogram-restore`

## Purpose

This directory replaces the disliked slide 06a named-clade violin/grouped view with a simpler histogram-style distribution. The asset answers one question directly:

How long are called inter-chromosomal PHR intervals, given that discovery was limited to terminal 500 kb windows?

The restored plot intentionally avoids the artificial named-clade grouping used in the current v6 slide 06a. It plots all non-empty called intervals in one 25 kb histogram, marks the median, and makes the 500 kb measurement ceiling visible.

## Data Source

The generator reads the canonical length TSV directly:

`/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv`

Override with `PHR_LENGTH_TSV=/path/to/file.tsv` if the source moves.

The plotted length is computed as:

`(region_end - region_start) / 1000`

Rows with `region_start = .` are zero-signal rows and are excluded from the interval-length histogram because they do not have a called interval length.

## Analysis Ceiling

This analysis searched/measured only within terminal 500 kb windows. The exact ceiling wording is printed on the plot:

`analysis window ends at 500 kb; longer shared sequence is not measured`

Values beyond 500 kb are not measured. They should not be interpreted as absent.

## Generated Files

| File | Purpose |
| --- | --- |
| `make_06a_length_histogram_restore.R` | Reproducible generator using base R plus `ggplot2`. |
| `phr_length_histogram_restore.png` | Talk-ready 16:9 PNG replacement asset for slide 06a. |
| `phr_length_histogram_restore.pdf` | PDF version of the same asset. |
| `length_distribution_summary.tsv` | Source, counts, quantiles, cap-hit count, and ceiling note. |
| `histogram_bins_25kb.tsv` | Exact 25 kb bin counts used by the plot. |
| `SLIDE_PATCH.md` | Fan-in instructions for replacing slide 06a. |
| `VALIDATION.md` | Render and visibility validation notes. |

## Regenerate

From the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v7/06a_length_histogram_restore/make_06a_length_histogram_restore.R
```

Optional source override:

```bash
PHR_LENGTH_TSV=/path/to/all-vs-all.1Mb.p95.id95.len.tsv \
  Rscript slides/v2-review-zoom/_revision_assets/v7/06a_length_histogram_restore/make_06a_length_histogram_restore.R
```

## Summary

The rendered asset reports:

- `15,668` non-empty called PHR intervals.
- `41/48` arms with called intervals.
- Median called interval length: `105 kb`.
- 90th percentile: `330 kb`.
- 95th percentile: `500 kb`.
- `1,087` intervals reported at the 500 kb ceiling.
- Zero-signal arms in the TSV: `2p, 3p, 5p, 8q, 11q, 14q, 18q`.

## Prior Assets Inspected

- Older histogram assets: `slides/v2-review-zoom/_typst/assets/s06_length_dist*.png`.
- v2/v1 length redesign: `slides/v2-review-zoom/_revision_assets/06_length_redesign/`.
- v3 violin/censor candidate: `slides/v2-review-zoom/_revision_assets/v3/06_violin_censor/`.
- Current v6 deck source for slide conventions: `slides/v2-review-zoom/_typst/zoom_review_deck.typ`.

No final Typst deck integration was performed here; `review-zoom-v7-fanin-render` owns that step.
