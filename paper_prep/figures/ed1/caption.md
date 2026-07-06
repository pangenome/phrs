# Extended Data Figure 1 — Pipeline and per-arm flank inventory

**(a)** Pipeline overview: 465 near-complete assemblies (232 HPRC v2
individuals and CHM13) → 18,827 telomere-anchored 500 kb subtelomeric flanks across 48
arms → wfmash all-vs-all `-p 95 -t 48` (18,827 PAFs) → impg sliding-window
scan (id ≥ 0.95, ≥ 2 chrs, ≥ 5 alns/chr, ≥ 3 kb output) → 15,668 PHRs
(83.2 %) → pggb + odgi similarity (15,668 × 15,668 Jaccard) → Leiden 15
arm-level (silhouette 0.347) and 50 sequence-level (modularity 0.97).
**(b)** Per-arm flank counts (n = 18,826 flanks across 48 arms; median
447/arm). Acrocentric p (chr13_p 76 to chr22_p 223) and sex chromosomes
(chrY_p 92 to chrX_q 330) are systematically under-represented.
*n** = mixed-strand-caveat contigs from `contig_classifications.tsv`
(zero on autosomes).
**(c)** PHR length distribution for the 15,668 signal-bearing flanks
(median 105 kb, mean 144 kb; saturates at 500 kb because input flanks are
clipped).
**(d)** Chr18_q chimera (NA18982#1, JBKABS010000018.1, 84.4 Mb): 966 kb of
chrX PAR1 fused to the chr18 backbone across a 100 bp NNN scaffold join at
~83.37 Mb; mapq 60 in both wfmash v0.23 and minimap2 v2.30; 2,826 bp
(TTAGGG)n × 471 precedes the gap. Removed from the PHR set (15,669 → 15,668).

**Source data.** `all-vs-all.1Mb.p95.id95.len.tsv` (b counts, c lengths) and
`pq-classification/contig_classifications.tsv` (b QC overlay). Full paths and
metrics in `sources.tsv`.
