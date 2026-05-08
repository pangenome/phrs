# Slide 06a Q-Axis and kbp Validation

Task: `review-zoom-v9-slide06a-q-axis-kbp`

## Render Command

From the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v9/06a_q_axis_kbp/make_06a_q_axis_kbp.R
```

This wrote the slide-ready PNG/PDF and source TSVs under:

`slides/v2-review-zoom/_revision_assets/v9/06a_q_axis_kbp/`

## Output Files

- `phr_length_arm_heatstrip_10kbp.png`
- `phr_length_arm_heatstrip_10kbp.pdf`
- `length_distribution_summary.tsv`
- `arm_length_summary.tsv`
- `arm_length_bins_10kbp.tsv`
- `asset_manifest.tsv`

## File Checks

`file slides/v2-review-zoom/_revision_assets/v9/06a_q_axis_kbp/phr_length_arm_heatstrip_10kbp.png` reports a 3072 x 1728 RGB PNG, matching a 16:9 slide at 240 dpi.

`file slides/v2-review-zoom/_revision_assets/v9/06a_q_axis_kbp/phr_length_arm_heatstrip_10kbp.pdf` reports a PDF document.

## Visual Checks

- The p-arm panel is on the left and reads 0 to 500 kbp from left to right.
- The q-arm panel is on the right and reads 500 to 0 kbp from left to right.
- The combined p/q layout reads p telomere to q telomere.
- The 10 kbp per-arm/end heatstrip concept from v8 is preserved.
- The terminal-window callout reads `>500 kbp not measured`.
- The visible slide text uses kbp units.

## Data Checks

`length_distribution_summary.tsv` confirms:

- q panel display is flipped: left edge is 500 kbp, right edge is 0 kbp.
- bin width: `10 kbp`.
- non-empty called PHR intervals: `15,668`.
- chromosome ends with called intervals: `41/48`.
- zero-signal ends: `2p, 3p, 5p, 8q, 11q, 14q, 18q`.
- median called interval length: `105 kbp`.
- 90th percentile: `330 kbp`.
- 95th percentile: `500 kbp`.
- intervals reported exactly at 500 kbp: `1,087`.
- terminal-window note: `analysis measured terminal 500 kbp windows; >500 kbp was not measured`.
