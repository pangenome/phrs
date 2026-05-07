# Community Assignment Method Slide

Task: `review-zoom-v6-community-assignment-method-slide`

Deck placement: slide `07j.2`, inserted after the PGGB graph slide `07j.1`
and before the tree/community heatmap sequence.

## Purpose

This slide gives a concise, non-defensive explanation of how the arm-level
PHR communities `C1`-`C15` were assigned. It is meant to make the heatmaps and
community calls credible without turning the talk into a methods dump.

## Slide Claims

- Start with 15,668 PHR paths from HPRCv2 haplotypes/arms with detected
  inter-chromosomal PHR signal.
- Build one PGGB graph with `pggb -p 95`, then compute all-vs-all graph-path
  Jaccard with `odgi similarity --all -P`.
- Convert graph-path Jaccard to distance with `distance = 1 - Jaccard`.
- Collapse path-level distances to a 41 x 41 arm-level distance matrix by
  averaging all haplotype/path pair distances for each arm pair `A x B`.
- Run Leiden on a fully connected weighted graph of arms, using
  `w_ij = exp(-d_ij / median(d))`.
- Select the Leiden resolution by a silhouette scan from 0.1 to 3.0 in 0.01
  steps. The selected arm-level partition has 15 communities, resolution 1.16,
  and mean silhouette 0.347.
- Compare against UPGMA average-linkage clustering on the same distance matrix:
  14 communities, mean silhouette 0.342, and exact agreement on 12 of 15 Leiden
  communities. The differences are boundary cases around f7501-like arms.
- Biological labels and interpretations such as D4Z4, acrocentric p, PAR1,
  PAR2, f7501, and OR4F were added after clustering and were not inputs.

## Caveats Preserved

- The slide describes the arm-level `C1`-`C15` partition only. The separate
  sequence-level 50-community partition is a finer-grained analysis and is not
  used to define these arm-level calls.
- The 41 x 41 matrix includes only arms with detected inter-chromosomal PHR
  signal. Seven zero-signal arms were excluded rather than clustered.
- The slide does not imply that CHM13 called PHR intervals exist for every
  community-assigned arm.
- No gene labels, biological labels, or 3D data were used to define the
  communities.

## Source Anchors

- `subtelomeric_analysis_report.md`, sections 5 and 6.1.
- `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R`.
- `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`.
- `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`.

## Data Audits

Commands run from the repository worktree:

```bash
wc -l /moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv \
  /moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv

awk -F '\t' 'NR==1{print NF-1 " arms in matrix header"} NR>1{n++} END{print n " matrix data rows"}' \
  /moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv

awk -F '\t' 'NR>1{c[$2]++} END{for (k in c) print k, c[k]}' \
  /moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv | sort -V
```

Observed results:

- Both matrix and assignment files have 42 lines, consistent with a header plus
  41 arms.
- The distance matrix has 41 arm columns and 41 data rows.
- The Leiden assignments contain `C1` through `C15`, with community sizes
  `2, 2, 6, 2, 4, 6, 5, 1, 1, 1, 4, 2, 1, 2, 2`.

## Slide Asset

- `community_assignment_method_schematic.svg` is a hand-authored schematic of
  the assignment workflow. It does not introduce new analysis; it visualizes the
  already documented path-to-graph-to-arm-matrix-to-clustering procedure.
