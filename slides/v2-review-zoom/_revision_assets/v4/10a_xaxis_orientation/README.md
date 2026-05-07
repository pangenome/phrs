# Review zoom v4 slide 10a X-axis orientation fix

Task: `review-zoom-v4-slide10a-xaxis-orientation-fix`.

## Outputs

- `candidate_10a_xaxis_orientation.png`: corrected square PNG candidate (1800x1800).
- `candidate_10a_xaxis_orientation.pdf`: corrected square PDF candidate (10x10 inches).
- `make_10a_xaxis_orientation.R`: generator with source/statistic assertions and display-axis audit.
- `orientation_audit.tsv`: per-index row order, corrected displayed X order, v3 implicit displayed X order, and corrected X-axis policy.
- `ordered_arm_haplotypes.tsv`: final displayed row/column order plus the v3 implicit X order.
- `sequence_community_boxes.tsv`: exact displayed-coordinate community box coordinates.
- `matrix_order_audit.tsv`: source, statistic, and orientation checks.

## Source Files Inspected

- v3 slide 10a generator: `slides/v2-review-zoom/_revision_assets/v3/10a_axis_box_fix/make_10a_axis_box_fix.R`.
- v3 slide 10a asset: `slides/v2-review-zoom/_revision_assets/v3/10a_axis_box_fix/candidate_10a_axis_box_fix.png` (1800x1800).
- v3 audit: `slides/v2-review-zoom/_revision_assets/v3/10a_axis_box_fix/matrix_order_audit.tsv`.
- Original manuscript Fig. 3 generator: `paper_prep/figures/fig3/figure_fig3.R:32-80`; panel A uses base `image(..., t(vals_norm)[, n:1])`.
- v2 redesign generator: `slides/v2-review-zoom/_revision_assets/hic_visual_redesign/make_hic_visual_redesign.R:28-93`; it also used the base-image transform shape with `rasterImage()`.
- Contact matrix: `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_contact_matrix.tsv`.
- Sequence community table for boxes: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`.
- Source B/W and p-value TSV: `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_global_test.tsv`.

## Orientation Finding

- The v3 validation asserted that analytical row and column names were identical after ordering. That is true, but it does not validate the displayed X axis.
- The v3 renderer used `rasterImage(as.raster(t(color_matrix)[, n:1]))`. That transform is appropriate for base `image()` semantics in the original Fig. 3 panel, but `rasterImage()` already displays raster rows top-to-bottom and columns left-to-right.
- Therefore v3's visual X axis was mirrored left/right: its implicit displayed column order was `reverse(ordered_arms)`. The symmetric contact matrix made the mistake easy to miss by visual block structure and by row/column-name equality checks.
- Corrected X-axis policy: display columns left-to-right as ordered_arms[1:n], matching the sequence-community order; draw rasterImage(as.raster(color_matrix)) directly and compute x boxes in this displayed column coordinate system.
- Corrected row-axis policy: display rows top-to-bottom as ordered_arms[1:n]; compute y boxes with y_min = n - row_end + 0.5 and y_max = n - row_start + 1.5.

## Statistic Check

- Source within mean = 0.0274578380817 over 256 within-community pairs.
- Source between mean = 0.00152475029036 over 2670 between-community pairs.
- B/W = between / within = 0.055530602, displayed as 0.056.
- p-value = 3.856e-85, displayed as 3.9e-85.

The generator recomputes the within and between means from the HG002 Pore-C contact matrix plus the expanded sequence-community table, then asserts that the recomputed values match `hg002_porec_global_test.tsv`. The B/W statistic, p-value, and sequence-community interpretation are preserved because those checks pass.

## Validation

- v3 X-axis mirrored stated in `orientation_audit.tsv`: TRUE.
- Corrected displayed X axis starts at `chr4_MATERNAL_q` and ends at `chrY_PATERNAL_p`.
- Corrected PNG square: TRUE.
- Community boxes are computed in displayed coordinates after the corrected raster policy.
- The corrected asset is generated from source matrix values, not copied from the v3 PNG.
