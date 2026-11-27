**Figure 5 | Concerted evolution of human subtelomeres.**

**a**, For each Leiden community, the median divergence among member chromosome arms divided by the
median allelic divergence of those arms. A value of 1.0 (dotted line) means non-homologous copies are as
similar as alleles. The acrocentric short arms (C7) and PAR1 (C15) sit at 1.0, the obligate case, and
PAR2 (C14) at 1.4; the cohesive communities C1, C11 and C3 sit at 3.3-3.5; the remainder are
facultative, from 7 to 27. C4 is a residual pairing, not a cohesive group (its two arms share almost no
near-identical sequence, and chr12q is not among the closest non-homologous partners of chr7q).

**b**, Cumulative distribution of substitution divergence, weighted by aligned base pairs, between
homologous copies of an arm (allelic, median 0.15%) and copies on different chromosomes (paralogous,
median 0.89%).

**c**, Paralogous-to-allelic divergence ratio against distance from the telomere-subtelomere junction,
one curve per chromosome arm (grey) and the median across 48 arms (black). The ratio rises from about
1.3 within the first 8 kb to 6-7 beyond 300 kb. The hatched zone is excluded: it is the TERRA promoter
array, a copy-variable repeat family whose copy number and phase differ between haplotypes, so same-arm
alignments there are short and partly paralogous and the ratio is undefined rather than imprecise.

**d**, All candidate ends offered simultaneously to the same similarity computation, scored against the
daughter's paternal chr9q across the tract that was previously called as an exchange. The father's chr9q
rebuilt from his own high-fidelity reads alone matches at graph dice 1.0000 inside the tract and 0.9998 outside,
and his chr3q is second at 0.9985 inside the tract. The released assembly matches at 0.9918, the same value as his
chr1p in every window of its last 60 kb, because that tail is a copy of his chr1p. His own other chr9q haplotype
has no alignment in the tract (0.9072 over the 47 outside windows it covers), so a non-homologous end carries a
block that one of his homologous copies lacks. Values: candidate_block_identity_9q.tsv, from
scripts/pedigree/fig5_v131_check/fig5d_competition.sbatch on the final corrected v1.1 assemblies (2026-11-26).

Divergence in **a**-**c** is X/(=+X) from alignment CIGARs in 2 kb windows over the terminal 500 kb of
each arm, inter-chromosomal pairs only. Dice in **d** is IMPG node-overlap similarity in 2 kb windows.
