# Figure 1 — Population-scale subtelomeric communities (the landscape)

**Caption.** Inter-chromosomal subtelomeric sharing partitions 41 chromosome
arms into 15 communities across 465 near-complete assemblies
(18,827 telomere-anchored 500 kb flanks; 15,668 PHRs).
**(a)** Genome-wide stacked identity heatmap, 465 near-complete assemblies × 24
chromosomes (per-position maximum identity to any matching chromosome,
100 kb windows). Telomere-anchored high-identity blocks mark the
inter-chromosomal exchange landscape.
**(b)** Genome-wide heatmap of the number of chromosomes sharing each
subtelomeric position with the called PHR-BED overlay.
**(c)** 41 × 41 arm-level Jaccard distance heatmap
(`hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`), arms ordered by Leiden
k = 15 community
(`hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`); cyan rectangles
delimit the 15 communities. The UPGMA k = 14 dendrogram (top, computed on
the same matrix) recovers 14 of 15 Leiden communities (agreement 14/15).
**(d)** Per-arm architecture category from Leiden membership:
*fully interchangeable* (red, 9/41 — C7 acrocentric p-arms + C14 PAR2 +
C15 PAR1, all negative-silhouette communities with allele–paralog distance
reversed or near equal; SURVEY_04 §1.1–§1.2);
*homogeneous* (blue, 4/41 — single-arm communities C8/C9/C10/C13);
*polymorphic* (teal, 28/41 — multi-arm members with arm identity preserved).
Bar height = per-arm cross-arm sequence rate (cross-arm sequences in
`cross_arm_affinity_sequences.tsv` ÷ total in
`hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv`); chrX_q 99.7 %,
chr21_p 94.0 %, chr11_p 74.1 % (highest autosomal). Counts refine the
preliminary 8/34/7 skeleton split into 4/28/9.
