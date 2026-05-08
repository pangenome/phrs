# Slide 06a Length Alternatives Validation

Task: `review-zoom-v8-slide06a-length-alternatives`

## Render Validation

Command run from the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v8/06a_length_alternatives/make_06a_length_alternatives.R
```

The generator wrote all expected candidate assets:

- `phr_length_histogram_10kb.png`
- `phr_length_histogram_10kb.pdf`
- `phr_length_histogram_5kb_sensitivity.png`
- `phr_length_histogram_5kb_sensitivity.pdf`
- `phr_length_arm_heatstrip_10kb.png`
- `phr_length_arm_heatstrip_10kb.pdf`

It also wrote the supporting summary tables:

- `length_distribution_summary.tsv`
- `histogram_bins_10kb.tsv`
- `histogram_bins_5kb.tsv`
- `arm_length_summary.tsv`
- `arm_length_bins_10kb.tsv`
- `asset_manifest.tsv`

## Dimension Check

`file slides/v2-review-zoom/_revision_assets/v8/06a_length_alternatives/*.png` reports all three PNGs as 3072 x 1728 RGB images, matching a 16:9 slide at 240 dpi.

`identify` was not available in this worktree, so `file` was used for the image dimension check.

## Visual Legibility Check

Manually inspected:

- `phr_length_histogram_10kb.png`
- `phr_length_histogram_5kb_sensitivity.png`
- `phr_length_arm_heatstrip_10kb.png`

Confirmed:

- Axis labels are readable at 16:9 slide size.
- The 10 kb histogram is not dominated by wide, blocky 25 kb bins.
- The 5 kb sensitivity histogram is readable, but busier than the 10 kb render.
- The heatstrip directly represents every chromosome end as a row, split into p/q panels, without a 41-color legend.
- The heatstrip labels zero-signal ends directly: `2p, 3p, 5p, 8q, 11q, 14q, 18q`.
- Each candidate has visible 500 kb ceiling text and/or right-edge marking.
- The ceiling wording makes clear that values at or near 500 kb are right-censored by the terminal-window search, not evidence that longer shared sequence is absent.

## Numeric Checks

`length_distribution_summary.tsv` confirms:

- current v7 bin width: `25 kb`
- proposed v8 bin width: `10 kb`
- sensitivity bin width: `5 kb`
- non-empty intervals: `15,668`
- ends with called intervals: `41/48`
- median length: `105 kb`
- 90th percentile: `330 kb`
- 95th percentile: `500 kb`
- exact 500 kb calls: `1,087`

## Repository Checks

- `git diff --check` passed with no whitespace errors.
- `find . -maxdepth 3 -name Cargo.toml -print` returned no Cargo project at this repo depth, so `cargo build` and `cargo test` are not applicable for this documentation/render task.

No final Typst deck files were edited.
