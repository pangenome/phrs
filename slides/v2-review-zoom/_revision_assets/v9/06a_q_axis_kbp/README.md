# Slide 06a Q-Axis and kbp Patch

Task: `review-zoom-v9-slide06a-q-axis-kbp`

## Purpose

This directory provides the v9 replacement asset for slide 06a. It keeps the v8 per-arm/end heatstrip concept, but changes the q-arm orientation and standardizes slide-visible units to kbp.

## Recommended Asset

Use this slide-ready heatstrip:

`phr_length_arm_heatstrip_10kbp.png`

Use this companion if the fan-in render workflow prefers PDF figures:

`phr_length_arm_heatstrip_10kbp.pdf`

## Changes from v8

The p-arm panel remains left-to-right from 0 to 500 kbp.

The q-arm panel is flipped so its displayed axis runs 500 to 0 kbp from left to right. With p arms on the left and q arms on the right, the figure now reads naturally from p telomere to q telomere.

Slide-visible wording now uses `>500 kbp not measured`. The figure title, subtitle, axis, callouts, and source line use kbp for kilobase pairs.

## Data Source

The generator reads:

`/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv`

Override with `PHR_LENGTH_TSV=/path/to/file.tsv` if the source moves.

The plotted length is computed as:

`(region_end - region_start) / 1000`

Rows with `region_start = .` are zero-signal rows. They are excluded from interval-length calculations and shown in the heatstrip as `no called interval`.

## Generated Files

| File | Purpose |
| --- | --- |
| `make_06a_q_axis_kbp.R` | Reproducible v9 generator using base R plus `ggplot2`. |
| `phr_length_arm_heatstrip_10kbp.png` | 16:9 PNG, recommended slide-ready heatstrip. |
| `phr_length_arm_heatstrip_10kbp.pdf` | PDF companion for the heatstrip. |
| `length_distribution_summary.tsv` | Source path, bin width, q-axis orientation note, counts, quantiles, and terminal-window note. |
| `arm_length_summary.tsv` | Per-end counts and quantiles using kbp column names. |
| `arm_length_bins_10kbp.tsv` | Per-end 10 kbp bins used by the heatstrip, including the plotted x-position after the q-axis flip. |
| `asset_manifest.tsv` | Machine-readable asset inventory. |
| `SLIDE_PATCH.md` | Fan-in instructions for v9 deck integration. |
| `VALIDATION.md` | Render and wording validation notes. |

## Regenerate

From the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v9/06a_q_axis_kbp/make_06a_q_axis_kbp.R
```

Optional source override:

```bash
PHR_LENGTH_TSV=/path/to/all-vs-all.1Mb.p95.id95.len.tsv \
  Rscript slides/v2-review-zoom/_revision_assets/v9/06a_q_axis_kbp/make_06a_q_axis_kbp.R
```

## Rendered Summary

The rendered v9 heatstrip reports:

- `15,668` non-empty called PHR intervals.
- `41/48` chromosome ends with called intervals.
- Zero-signal ends: `2p, 3p, 5p, 8q, 11q, 14q, 18q`.
- Median called interval length: `105 kbp`.
- 90th percentile: `330 kbp`.
- 95th percentile: `500 kbp`.
- `1,087` intervals reported exactly at 500 kbp.

No final Typst deck integration was performed here; `review-zoom-v9-fanin-render` owns that step.
