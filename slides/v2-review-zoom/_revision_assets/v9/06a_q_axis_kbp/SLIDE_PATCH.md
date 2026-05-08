# Slide 06a Patch Instructions

Task: `review-zoom-v9-slide06a-q-axis-kbp`

## Recommendation

Use this asset for slide 06a:

`slides/v2-review-zoom/_revision_assets/v9/06a_q_axis_kbp/phr_length_arm_heatstrip_10kbp.png`

Use the matching PDF if the fan-in deck workflow prefers PDF figures:

`slides/v2-review-zoom/_revision_assets/v9/06a_q_axis_kbp/phr_length_arm_heatstrip_10kbp.pdf`

This preserves the v8 10 kbp per-arm/end heatstrip while fixing the q-arm orientation. The p panel reads 0 to 500 kbp from left to right. The q panel reads 500 to 0 kbp from left to right, so the combined layout follows the normal p telomere to q telomere chromosome convention.

## Suggested Slide Label

`PHR length distribution by chromosome end`

## Suggested Source Footer

`v9/06a_q_axis_kbp/make_06a_q_axis_kbp.R; arm_length_bins_10kbp.tsv; >500 kbp not measured`

## Suggested Caption or Speaker Framing

This slide shows called PHR interval lengths by chromosome end. Rows are chromosome ends; q arms are flipped so the heatstrip reads p telomere to q telomere. Color is the within-end share of calls in each 10 kbp bin. Terminal windows were 500 kbp, so `>500 kbp` was not measured.

## Integration Constraint

Do not edit the final Typst deck directly in this task. Integration belongs to `review-zoom-v9-fanin-render`.
