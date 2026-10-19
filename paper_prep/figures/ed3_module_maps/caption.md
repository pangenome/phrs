# Extended Data Figure 3 - Subtelomeric module maps and module architecture

**Module maps.** Each chromosome end is cut into 25 kb bins; each bin is labelled with the
partner arm that contributes the most aligned base pairs to it, aggregated across 233 sampled
individuals in the HPRCv2 all-vs-all flank alignments. Fill colour is that partner's arm-level
Leiden community (C1-C15). Tiles at least three bins wide carry the partner arm label. 26 of the
46 arm ends receive enough inter-arm alignment to be tiled; the most frequent top partners are
the acrocentric p-arms (chr14p, chr22p).

**What the panel shows about architecture.** The shared modules are reassorted in order but
consistent in orientation. Across 764 arm pairs and 9.6 Gbp of inter-arm alignment, 98.6% of
matched base pairs lie in the same orientation relative to the telomere (P-P 97.8%, Q-Q 100.0%,
P-Q 98.4%, with homolog-versus-homolog pairs excluded); the per-pair median is 100% and only
7 of 764 pairs fall below 90%. Order is not conserved: concordance across shared partners is
negligible and does not improve within communities. What tracks community membership is the
module repertoire (Jaccard 0.83 within versus 0.57 between communities, permutation
p = 5e-10) and the length of shared blocks (median 13.3 versus 9.9 kb).

Sources: HPRCv2 all-vs-all telomeric flank alignments; `scripts/cladistics/arm_module_{tiling,
orientation,order}.tsv`; `data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`.
