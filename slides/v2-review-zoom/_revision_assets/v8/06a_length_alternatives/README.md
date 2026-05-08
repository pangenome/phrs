# Slide 06a Length Alternatives

Task: `review-zoom-v8-slide06a-length-alternatives`

## Purpose

This directory makes slide 06a talk-ready for v8 by replacing the v7 wide-bin view with finer-bin candidates and a readable per-chromosome-end representation.

The goal is to keep the part of v7 that works - an immediately understandable length distribution - while recovering the older all-chromosome/end information without a 41-entry color legend.

## Bin-Size Confirmation

Current v7 is a 25 kb-bin histogram, not 10 kb.

The v7 source confirms this in three places:

- `slides/v2-review-zoom/_revision_assets/v7/06a_length_histogram_restore/make_06a_length_histogram_restore.R` sets `binwidth_kb <- 25`.
- `slides/v2-review-zoom/_revision_assets/v7/06a_length_histogram_restore/histogram_bins_25kb.tsv` is the exact binned output.
- `slides/v2-review-zoom/_revision_assets/v7/06a_length_histogram_restore/README.md` describes "one 25 kb histogram" and lists `histogram_bins_25kb.tsv`.

The proposed v8 main bin width is 10 kb. A 5 kb sensitivity render is included because the underlying detection window uses a 5 kb step.

## Data Source

The generator reads the canonical length TSV directly:

`/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv`

Override with `PHR_LENGTH_TSV=/path/to/file.tsv` if the source moves.

The plotted length is computed as:

`(region_end - region_start) / 1000`

Rows with `region_start = .` are zero-signal rows and are excluded from the interval-length histograms because they do not have a called interval length. They are still represented in the heatstrip as no-call chromosome ends.

## Analysis Ceiling

This analysis searched/measured only terminal 500 kb windows. Values at or near 500 kb are right-censored by that analysis ceiling.

The figure text therefore avoids implying that sharing stops at 500 kb. The exact summary note written by the generator is:

`analysis searched/measured terminal 500 kb windows; values at/near 500 kb are right-censored`

## Recommendation

Primary recommendation when slide 06a should preserve the older all-chromosome/end information:

`phr_length_arm_heatstrip_10kb.png`

This heatstrip uses one row per chromosome end, split into p-arm and q-arm panels. Color encodes the within-end share of calls in each 10 kb length bin, so rows remain comparable without requiring a separate color for each of 41 signaled ends. The 500 kb edge is outlined and labeled as the ceiling.

Fastest single-distribution alternative:

`phr_length_histogram_10kb.png`

This keeps the v7 histogram idea, documents the 10 kb bin width on the axis/footer, and removes the wide 25 kb-bin look. Use it if the fan-in slide needs one quickly legible global distribution rather than per-end detail.

Sensitivity only:

`phr_length_histogram_5kb_sensitivity.png`

This is readable and shows the discrete 5 kb-step structure, but it is visually busier. Do not use it as the main slide unless Erik explicitly wants the finest-bin view.

## Generated Files

| File | Purpose |
| --- | --- |
| `make_06a_length_alternatives.R` | Reproducible generator using base R plus `ggplot2`. |
| `phr_length_histogram_10kb.png` | 16:9 PNG, recommended simple 10 kb histogram candidate. |
| `phr_length_histogram_10kb.pdf` | PDF companion for the 10 kb histogram. |
| `phr_length_histogram_5kb_sensitivity.png` | 16:9 PNG, 5 kb sensitivity candidate. |
| `phr_length_histogram_5kb_sensitivity.pdf` | PDF companion for the 5 kb sensitivity candidate. |
| `phr_length_arm_heatstrip_10kb.png` | 16:9 PNG, recommended per-end old-style replacement. |
| `phr_length_arm_heatstrip_10kb.pdf` | PDF companion for the per-end heatstrip. |
| `length_distribution_summary.tsv` | Source, v7/proposed bin sizes, counts, quantiles, and ceiling note. |
| `histogram_bins_10kb.tsv` | Exact 10 kb global histogram bins. |
| `histogram_bins_5kb.tsv` | Exact 5 kb global histogram bins. |
| `arm_length_summary.tsv` | Per-end counts and quantiles. |
| `arm_length_bins_10kb.tsv` | Per-end 10 kb bins used by the heatstrip. |
| `asset_manifest.tsv` | Machine-readable asset inventory and intended use. |
| `SLIDE_PATCH.md` | Fan-in instructions for v8 deck integration. |
| `VALIDATION.md` | Render and legibility validation notes. |

## Regenerate

From the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v8/06a_length_alternatives/make_06a_length_alternatives.R
```

Optional source override:

```bash
PHR_LENGTH_TSV=/path/to/all-vs-all.1Mb.p95.id95.len.tsv \
  Rscript slides/v2-review-zoom/_revision_assets/v8/06a_length_alternatives/make_06a_length_alternatives.R
```

## Rendered Summary

The rendered assets report:

- `15,668` non-empty called PHR intervals.
- `41/48` chromosome ends with called intervals.
- Zero-signal ends: `2p, 3p, 5p, 8q, 11q, 14q, 18q`.
- Median called interval length: `105 kb`.
- 90th percentile: `330 kb`.
- 95th percentile: `500 kb`.
- `1,087` intervals reported exactly at the 500 kb ceiling.

No final Typst deck integration was performed here; `review-zoom-v8-fanin-render` owns that step.
