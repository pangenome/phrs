# Review zoom v3 slide 10a axis/box fix

Task: `review-zoom-v3-slide10a-axis-box-fix`.

## Outputs

- `candidate_10a_axis_box_fix.png`: square PNG candidate (1800x1800).
- `candidate_10a_axis_box_fix.pdf`: square PDF candidate (10x10 inches).
- `make_10a_axis_box_fix.R`: generator with assertions.
- `ordered_arm_haplotypes.tsv`: final row/column order.
- `sequence_community_boxes.tsv`: exact community box coordinates.
- `matrix_order_audit.tsv`: source and ordering checks.

## Source Files

- Contact matrix: `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_contact_matrix.tsv`.
- Sequence community table for boxes: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`.
- Source B/W and p-value TSV: `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_global_test.tsv`.
- v2 slide 10a asset inspected: `slides/v2-review-zoom/_revision_assets/hic_visual_redesign/slide_10a_square_matrix_candidate.png` (1800x1800).
- v2 redesign generator inspected: `slides/v2-review-zoom/_revision_assets/hic_visual_redesign/make_hic_visual_redesign.R:28-44` and `:67-80`.
- Manuscript Fig. 3 generator inspected: `paper_prep/figures/fig3/figure_fig3.R:32-45` and `:63-80`.

## Rows, Columns, Boxes

- rows = the 77 HG002 Pore-C arm-haplotypes from the first column of `hg002_porec_contact_matrix.tsv`, reordered by arm-level sequence community.
- columns = the same 77 HG002 Pore-C arm-haplotypes from the matrix header, in the identical post-order list.
- boxes = arm-level sequence communities from `hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`, expanded from base arms such as `chr4_q` to matrix labels such as `chr4_MATERNAL_q` and `chr4_PATERNAL_q`.

The final ordered row list and column list are asserted identical by name after ordering. The generator stops before plotting if they diverge.

## Audit Findings

- Source contact matrix shape: 77 x 77.
- Source row names equal source column names: TRUE.
- Source matrix symmetric: TRUE; max absolute symmetry delta = 0e+00.
- All 77 source matrix arm-haplotypes map to sequence communities: TRUE.
- Sequence communities represented in the plotted matrix: 15 (C1-C15).
- v2 box/order table communities represented: 42 contact-derived communities from `hg002_porec_hic.arm-leiden.communities.tsv`.
- v2 contact-community partition agreement with the sequence-community partition: 2712/2926 arm-pairs.
- Was the v2 plot transposed incorrectly? No. The source matrix rows and columns were identical, and the v2 renderer used the usual R display transform (`t(colors)[, n:1]`) to put the first row at the top. The real v2 defect for this task is that community boxes/order came from the contact-community table, not the sequence-community table behind the stated interpretation and B/W statistic.
- Candidate transpose policy: no analytical matrix transpose is applied. The raster is transformed only for R image coordinates, and the blue boxes are computed in that same displayed coordinate system after the visual y-axis reversal.

## Statistic Check

- Source within mean = 0.0274578380817 over 256 within-community pairs.
- Source between mean = 0.00152475029036 over 2670 between-community pairs.
- B/W = between / within = 0.055530602, displayed as 0.056.
- p-value = 3.856e-85, displayed as 3.9e-85.

The generator recomputes the within and between means from `hg002_porec_contact_matrix.tsv` plus the expanded sequence-community table, then asserts that the recomputed values match `hg002_porec_global_test.tsv`. The B/W and p-value are preserved in the candidate only because those checks pass.

## Validation

- README row/column/order checks included: TRUE.
- v2 transpose finding included: TRUE.
- Candidate PNG square: TRUE.
- Candidate PDF square: TRUE (10 x 10 inch device).
- Community boxes align on both axes by construction: `sequence_community_boxes.tsv` stores matching x and y extents for each contiguous sequence-community block, and the plot draws those extents after the display coordinate transform.
