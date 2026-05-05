# Extended Data Figure 2 — Sequence-level (50-community) detail

**(a)** UMAP layout (first 2 of 3 dims from
`hprcv2.1Mb.subtelo.umap.rds`) of all 15,668 PHRs, coloured by sequence-level
Leiden community (50; k = 75, res = 0.8, modularity 0.97). Top 8 labelled at
centroid: C4 (770, acro-p), C3 (712, D4Z4), C13 (508), C26 (490),
C25 (487), C18 (468), C9 (466), C8 (454).
**(b)** Within-community Jaccard distance (= 1 − jaccard.similarity)
distributions for the 8 arm-level communities flagged as bimodal in
`SURVEY_04 §1.10` (C1, C2, C3, C5, C6, C7, C11, C12). Up to 50,000 random
within-community pairs are shown per panel; pairs were extracted by streaming
`similarity.tsv.gz` (~12.5 GB) once with
`extract_within_community_jaccard.sh`. Allele/paralog separation is visible in
C5, C7, C11.
**(c)** Cross-arm affinity radial plot. 2,484 cross-arm sequences from
`cross_arm_affinity_sequences.tsv` (cross_arm_affinity > 1) → 48 (origin →
affinity) edges between 41 partition arms; chord width ∝ √(#sequences); chord
hue = arm-Leiden community of the destination. Top edges: chrY_q→chrX_q 329
(PAR2), chrY_p→chrX_p 287 (PAR1), chr15_p→chr14_p 280, chr10_q→chr4_q 220
(D4Z4).
**(d)** Confusion matrix between the arm-Leiden (15 rows) and seq-Leiden
(50 cols) partitions; cell colour = # PHRs (sqrt). ARI 0.35, NMI 0.76,
n = 15,668. C3 and C11 are the most fragmented; C1 (D4Z4) is the most
monolithic.

**Source data.** Listed in `sources.tsv`.
