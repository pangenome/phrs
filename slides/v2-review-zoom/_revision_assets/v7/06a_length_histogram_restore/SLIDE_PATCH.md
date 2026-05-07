# Slide 06a Patch Instructions

Task: `review-zoom-v7-slide06a-length-histogram-restore`

## Asset To Use

Replace the current slide 06a visual with:

`slides/v2-review-zoom/_revision_assets/v7/06a_length_histogram_restore/phr_length_histogram_restore.png`

Use the matching PDF only if the fan-in deck workflow prefers PDF figures:

`slides/v2-review-zoom/_revision_assets/v7/06a_length_histogram_restore/phr_length_histogram_restore.pdf`

## Suggested Slide Label

`PHR lengths within the 500 kb discovery window`

## Suggested Source Footer

`v7/06a_length_histogram_restore/make_06a_length_histogram_restore.R; length_distribution_summary.tsv; 500 kb analysis ceiling`

## Suggested Caption Or Speaker Note

Histogram shows all non-empty called inter-chromosomal PHR intervals from the length TSV, with 25 kb bins and no named-clade grouping. The right edge is a measurement ceiling: analysis window ends at 500 kb; longer shared sequence is not measured.

## Speaker Framing

This slide is a measured-length distribution, not proof that sharing stops at 500 kb. Calls were searched only inside terminal 500 kb windows, so the 500 kb pile-up should be read as right-censoring at the analysis window.

## Integration Constraint

Do not keep the v6 named-clade violin/grouped asset on slide 06a. Erik rejected that grouping as hard to interpret. Use this restored histogram asset for slide 06a and leave detailed clade naming to the surrounding slide sequence.
