# Neighbor-joining tree of HPRC v2 subtelomeric arm Jaccard distances

`nj_tree.R` consumes the precomputed 41x41 arm-level Jaccard distance matrix at
`/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`,
builds an unrooted neighbor-joining tree with `ape::nj`, roots it at the MRCA
of the acrocentric short-arm clade (which is monophyletic), and renders an
annotated tree in PDF and PNG. Every named clade highlighted in the abstract
is recovered as a monophyletic group on the NJ tree and corresponds one-to-one
to a community in the existing Leiden k=15 arm-level partition
(`hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`): PAR1/Xp,Yp ↔ C15;
PAR2/Xq,Yq ↔ C14; acrocentric short arms 13p/14p/15p/21p/22p ↔ C7;
10p-18p ↔ C2; the tight q-arm clade 22q/21q/19q/1q/13q/17q ↔ C6; and the
DUX4-containing 4q/10q ↔ C1. A 1000-replicate perturbation bootstrap
(Gaussian noise added to the precomputed distances at sigma = 25% of the
off-diagonal IQR; ~1 s total) puts the support at the MRCA of every named
clade at 100%, indicating that all six abstract clades sit on edges that are
robust to substantial perturbation of the distance matrix; deeper/internal
backbone edges show much lower support (32-90%, see PDF). True character-level
bootstrap was not possible because the input is a derived distance summary
rather than an alignment.
