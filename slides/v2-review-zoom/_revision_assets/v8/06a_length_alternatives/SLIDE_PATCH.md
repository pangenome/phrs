# Slide 06a Patch Instructions

Task: `review-zoom-v8-slide06a-length-alternatives`

## Recommendation

Use this asset if slide 06a should recover the old all-chromosome/end information while staying legible:

`slides/v2-review-zoom/_revision_assets/v8/06a_length_alternatives/phr_length_arm_heatstrip_10kb.png`

Use the matching PDF if the fan-in deck workflow prefers PDF figures:

`slides/v2-review-zoom/_revision_assets/v8/06a_length_alternatives/phr_length_arm_heatstrip_10kb.pdf`

Rationale: this is the best compromise between the older per-end distribution plot and Erik's v7 preference for a clearer length-distribution slide. It represents every chromosome end as a row, splits p and q arms into separate panels, uses direct chromosome labels, and avoids an unreadable 41-end color legend. It uses 10 kb bins, not the v7 25 kb bins, and the 500 kb ceiling is labeled and outlined.

## Simple Histogram Alternative

If the v8 fan-in slide needs the fastest possible single global distribution, use:

`slides/v2-review-zoom/_revision_assets/v8/06a_length_alternatives/phr_length_histogram_10kb.png`

PDF companion:

`slides/v2-review-zoom/_revision_assets/v8/06a_length_alternatives/phr_length_histogram_10kb.pdf`

This keeps the v7 histogram concept but changes the bin width from 25 kb to 10 kb, making the distribution less blocky while preserving the explicit 500 kb right-censoring note.

## Sensitivity Asset

Do not use as the main slide unless Erik asks for the finest-bin version:

`slides/v2-review-zoom/_revision_assets/v8/06a_length_alternatives/phr_length_histogram_5kb_sensitivity.png`

PDF companion:

`slides/v2-review-zoom/_revision_assets/v8/06a_length_alternatives/phr_length_histogram_5kb_sensitivity.pdf`

The 5 kb version is readable and matches the detection step, but it is busier than the 10 kb candidates.

## Suggested Slide Label

For the heatstrip:

`PHR length distribution by chromosome end`

For the simple histogram fallback:

`Called PHR lengths in the 500 kb discovery window`

## Suggested Source Footer

For the heatstrip:

`v8/06a_length_alternatives/make_06a_length_alternatives.R; arm_length_bins_10kb.tsv; terminal 500 kb ceiling`

For the 10 kb histogram:

`v8/06a_length_alternatives/make_06a_length_alternatives.R; histogram_bins_10kb.tsv; terminal 500 kb ceiling`

## Bin-Size Note

Current v7 is 25 kb bins:

- v7 generator: `binwidth_kb <- 25`
- v7 binned output: `histogram_bins_25kb.tsv`

Proposed v8 main bin width is 10 kb. The 5 kb render is sensitivity only.

## Speaker Framing

This slide shows called interval lengths inside the discovery window. The analysis searched/measured terminal 500 kb windows, so values at or near 500 kb are right-censored by the window. The pile-up at the right edge should be read as "we reached the analysis ceiling," not as evidence that sharing stops there.

The heatstrip version restores the information content of the old all-chromosome/end plot without trying to color 41 chromosome ends separately. The audience can scan p and q panels by chromosome number and see which ends concentrate at short, intermediate, or ceiling-length PHR calls.

## Integration Constraint

Do not edit the final Typst deck directly in this task. Integration belongs to `review-zoom-v8-fanin-render`.
