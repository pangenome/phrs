## Within-arm heterogeneity analysis

**What it does.** Tests whether sequences from the same chromosome arm form a homogeneous cluster or show internal heterogeneity — specifically, whether some sequences are more similar to a foreign arm's sequences than to their own (cross-arm affinity).

**Mechanism caveat** (applies throughout the heterogeneity section): These analyses detect outcomes of sequence homogenization but cannot distinguish among mechanisms — gene conversion, unequal crossover, reciprocal exchange, or other processes may all contribute.

### Allele vs paralog distance

**What it does.** Tests whether alleles (maternal vs paternal copy at the same chromosome arm within the same individual) are more similar to each other than paralogs (the most similar sequence from a different arm in the same community). If inter-chromosomal exchange were so frequent that arm identity is erased, paralogs would be as close as or closer than alleles.

**How.** For each individual with two haplotypes in a multi-arm community, the allelic distance (Jaccard distance between the maternal and paternal copy at the same arm) is paired with the paralog distance (Jaccard distance to the closest sequence from a different arm in the same community). The **Wilcoxon signed-rank test** is applied to these paired differences — a non-parametric test for whether the median difference between paired observations is zero. It does not assume normality and is robust to outliers.

| Community | N pairs | Median allelic | Median paralog | % paralog closer | Wilcoxon p | Direction |
|-----------|---------|---------------|----------------|-----------------|------------|-----------|
| C2 | 424 | 0.023 | 0.365 | 0.0% | 3.3e-71 | Allele closer |
| C12 | 437 | 0.037 | 0.727 | 0.2% | 3.0e-73 | Allele closer |
| C11 | 854 | 0.098 | 0.364 | 14.4% | 2.7e-103 | Allele closer |
| C5 | 805 | 0.032 | 0.439 | 14.9% | 1.4e-90 | Allele closer |
| C6 | 1,307 | 0.073 | 0.250 | 16.1% | 3.1e-159 | Allele closer |
| C4 | 429 | 0.109 | 0.549 | 20.5% | 3.2e-55 | Allele closer |
| C3 | 1,269 | 0.170 | 0.368 | 34.9% | 7.0e-52 | Allele closer |
| C1 | 265 | 0.293 | 0.375 | 44.2% | 6.0e-5 | Allele closer |
| **C7** | **156** | **0.481** | **0.353** | **70.5%** | **2.0e-7** | **Paralog closer** |
| Overall | 5,946 | 0.069 | 0.387 | 20.4% | <1e-300 | Allele closer |

**Result.** Alleles are significantly more similar than paralogs in 8 of 9 multi-arm communities (overall Wilcoxon p < 1e-300). The sole exception is C7 (acrocentric p-arms), where paralogs are closer than alleles in 70.5% of individuals (p = 2.0e-7, now significant).

**Conclusion.** Despite inter-chromosomal exchange, arm-specific identity is maintained in most communities. C7's reversal is now statistically significant (p = 2.0e-7, 156 pairs), confirming complete inter-chromosomal homogenization of acrocentric p-arms (silhouette = −0.029). C1 (D4Z4) shows the weakest allele advantage (44.2% paralog closer), consistent with its intermediate homogenization (silhouette = 0.147).

### Arm separation (silhouette analysis)

**What it does.** Quantifies how well sequences from different arms within the same community can be distinguished.

| Community | Arms | N sequences | Silhouette | Separation ratio | Interpretation |
|-----------|------|-------------|------------|------------------|----------------|
| C2 | chr10_p, chr18_p | 889 | 0.888 | 8.38 | Excellent separation: arms share community membership but retain distinct sequence content |
| C12 | chr2_q, chr20_p | 903 | 0.820 | 5.45 | Excellent separation |
| C5 | chr6_p, chr9_p, chr12_p, chr20_q | 1,728 | 0.611 | 3.60 | Moderate separation |
| C4 | chr7_q, chr12_q | 895 | 0.561 | 2.20 | Moderate separation |
| C6 | chr1_q, chr13_q, chr17_q, chr19_q, chr21_q, chr22_q | 2,701 | 0.521 | 4.00 | Moderate |
| C11 | chr1_p, chr5_q, chr6_q, chr8_p | 1,783 | 0.517 | 2.50 | Moderate |
| C3 | chr3_q, chr7_p, chr9_q, chr11_p, chr16_q, chr19_p | 2,655 | 0.413 | 2.28 | Moderate |
| C1 | chr4_q, chr10_q | 714 | 0.147 | 1.19 | Poor separation: substantial arm mixing |
| C15 | chrX_p, chrY_p | 419 | −0.025 | 0.97 | No separation |
| C7 | chr13_p, chr14_p, chr15_p, chr21_p, chr22_p | 763 | −0.029 | 1.02 | Negative: arms are near-interchangeable |
| C14 | chrX_q, chrY_q | 431 | −0.163 | 0.79 | Negative: X and Y sequences more similar to each other than to same-chromosome sequences |

Note: C8 (chr15_q), C9 (chr16_p), C10 (chr17_p), and C13 (chr4_p) are single-arm communities — no within-community arm separation is computed.

**Result.** Two well-separated communities (C2 sil=0.888, C12 sil=0.820) contain arms that share duplicon content but retain distinct sequence identity. Three poorly separated communities reveal extensive homogenization: C14 (PAR2, sil=−0.163, X/Y more similar to each other than to same-chromosome sequences), C7 (acrocentric p-arms, sil=−0.029, near-interchangeable), C15 (PAR1, sil=−0.025).

**Conclusion.** Acrocentric p-arms (C7) are near-interchangeable, consistent with sequence homogenization at rDNA-adjacent subtelomeric regions. PAR sequences (C14, C15) are indistinguishable between X and Y, as expected for pseudoautosomal regions that undergo obligate meiotic recombination.

### Cross-arm affinity

**What it does.** Quantifies the fraction of sequences more similar to a foreign arm within their community than to their own arm's majority.

**Key metrics.** 2,484 sequences (15.9% of 15,668) show cross-arm affinity at the arm level.

**Result.** Notable cross-arm affinity rates:

| Community | Arm | Cross-arm | Self-arm | Rate | Affinity toward |
|-----------|-----|-----------|----------|------|-----------------|
| C14 | chrX_q | 329 | 1 | 99.7% | chrY_q |
| C3 | chr11_p | 332 | 116 | 74.1% | chr3_q, chr19_p, chr9_q |
| C15 | chrX_p | 305 | 22 | 93.3% | chrY_p |
| C5 | chr6_p | 245 | 172 | 58.8% | chr20_q |
| C6 | chr22_q | 231 | 219 | 51.3% | chr21_q |
| C7 | chr14_p | 190 | 39 | 83.0% | other acrocentric p-arms |
| C1 | chr4_q | 146 | 211 | 40.9% | chr10_q |
| C3 | chr9_q | 130 | 320 | 28.9% | other C3 arms |
| C7 | chr21_p | 109 | 7 | 94.0% | other acrocentric p-arms |
| C7 | chr22_p | 96 | 127 | 43.0% | other acrocentric p-arms |
| C7 | chr15_p | 91 | 28 | 76.5% | other acrocentric p-arms |

**Conclusion.** chrX_q (C14) has 99.7% cross-arm affinity (329/330): essentially all X sequences are closer to chrY_q than to other chrX_q sequences, reflecting PAR2 homogenization. chrX_p (C15) shows the same pattern at 93.3% (305/327) for PAR1. chr11_p (C3) has 74.1% cross-arm affinity (332/448), the highest autosomal rate, consistent with extensive duplicon sharing across the f7501 cluster. In C7, all five acrocentric p-arms show high cross-arm rates (chr21_p 94.0%, chr14_p 83.0%, chr15_p 76.5%, chr13_p 88.2%, chr22_p 43.0%), consistent with inter-chromosomal homogenization of rDNA-adjacent subtelomeric regions.

### Population structure in cross-arm affinity

**What it does.** Tests whether the frequency of cross-arm affinity (carrying a foreign chromosome's subtelomeric sequence) differs between human superpopulations (AFR, AMR, EAS, EUR, SAS). If subtelomeric exchange events occurred at different rates or times in different populations, the frequency of cross-arm haplotypes should differ between superpopulations.

**How.** For each arm/community pair, a 2x5 contingency table (cross-arm vs self-arm x 5 superpopulations) is tested with **Fisher's exact test** — a non-parametric test for association between two categorical variables that is exact (does not rely on asymptotic approximations) and is appropriate for small expected cell counts. P-values are BH-corrected across 19 arm/community pairs (from 11 multi-arm communities). All 232 individuals have superpopulation annotation.

**Key metrics.** 10 of 19 pairs show significant superpopulation bias (p_adj < 0.05):

| Community | Arm | Cross | Self | Cross-arm superpop distribution | Self-arm superpop distribution | p_adj |
|-----------|-----|-------|------|--------------------------------|-------------------------------|-------|
| C1 | chr4_q | 146 | 211 | AFR=52; AMR=13; EAS=23; EUR=17; SAS=20 | AFR=32; AMR=53; EAS=47; EUR=27; SAS=38 | 4.7e-04 |
| C3 | chr16_q | 86 | 363 | AFR=60; AMR=5; EAS=2; EUR=1; SAS=6 | AFR=54; AMR=76; EAS=76; EUR=61; SAS=64 | 4.7e-04 |
| C5 | chr6_p | 245 | 172 | AFR=25; AMR=44; EAS=57; EUR=43; SAS=49 | AFR=74; AMR=29; EAS=19; EUR=14; SAS=18 | 4.7e-04 |
| C15 | chrX_p | 305 | 22 | AFR=71; AMR=62; EAS=52; EUR=43; SAS=48 | AFR=18; AMR=0; EAS=0; EUR=0; SAS=0 | 4.7e-04 |
| C6 | chr19_q | 59 | 398 | AFR=20; AMR=12; EAS=0; EUR=8; SAS=11 | AFR=93; AMR=67; EAS=79; EUR=55; SAS=61 | 7.6e-04 |
| C3 | chr9_q | 130 | 320 | AFR=45; AMR=17; EAS=19; EUR=18; SAS=12 | AFR=66; AMR=62; EAS=58; EUR=43; SAS=59 | 0.015 |
| C6 | chr22_q | 231 | 219 | AFR=49; AMR=41; EAS=43; EUR=28; SAS=46 | AFR=66; AMR=38; EAS=34; EUR=34; SAS=21 | 0.028 |
| C15 | chrY_p | 10 | 82 | AFR=5; AMR=1; EAS=0; EUR=1; SAS=0 | AFR=13; AMR=15; EAS=20; EUR=9; SAS=18 | 0.028 |
| C1 | chr10_q | 26 | 331 | AFR=4; AMR=3; EAS=3; EUR=1; SAS=11 | AFR=79; AMR=63; EAS=60; EUR=45; SAS=49 | 0.030 |
| C11 | chr6_q | 12 | 428 | AFR=1; AMR=1; EAS=0; EUR=4; SAS=4 | AFR=107; AMR=78; EAS=79; EUR=55; SAS=64 | 0.032 |

Notable patterns: chr16_q (C3) cross-arm is 70% AFR (60/86), consistent with the f7501 AFR-enrichment in the arm-level community detection. chr4_q (C1) cross-arm is AFR-enriched (52/146 = 36% vs baseline 28.9%). chrX_p (C15) self-arm is entirely AFR (18/18), indicating that the rare non-cross-arm PAR1 haplotypes are exclusively African.

**Fst across subtelomeric types**: **Fst** (fixation index) measures genetic differentiation between populations. Fst = 0 means allele frequencies are identical across populations; Fst = 1 means populations are fixed for different alleles. Here, the "allele" at each arm is binary: self-arm (0) or cross-arm (1). Hudson's Fst estimator is computed for each pair of superpopulations across the 10 arm/community pairs with the strongest signal and averaged:

| | AFR | AMR | EAS | EUR | SAS |
|---|-----|-----|-----|-----|-----|
| AFR | — | 0.102 | 0.152 | 0.108 | 0.103 |
| AMR | | — | 0.007 | 0.007 | 0.004 |
| EAS | | | — | −0.047 | 0.005 |
| EUR | | | | — | −0.003 |
| SAS | | | | | — |

**Result.** Mean Fst = 0.044. AFR is strongly differentiated from all non-AFR superpopulations (Fst 0.10–0.15), while non-AFR populations are nearly indistinguishable (Fst −0.05 to 0.01). The chr4_q cross-arm AFR enrichment (AFR=52/146 = 36% vs AFR=32/211 self = 15%) and the chr6_p cross-arm non-AFR enrichment are consistent with population-specific exchange rates at these loci.

**Conclusion.** The AFR/non-AFR Fst pattern parallels the known out-of-Africa demographic history. The data cannot distinguish population-specific exchange from drift or incomplete lineage sorting without outgroup comparison.

### Subtelomeric type discordance

**What it does.** Measures how often individuals are heterozygous for subtelomeric structural type — carrying one haplotype that matches the arm's own subtelomeric sequence (self-arm) and one that resembles a foreign chromosome's subtelomere (cross-arm). High discordance indicates that subtelomeric exchange is a common segregating polymorphism rather than a rare event.

**How.** Each individual's two haplotypes at a given arm are classified independently (the sequence-level community detection) as self-arm or cross-arm, then paired: **concordant** = both self-arm or both cross-arm; **discordant** = one of each. Discordance rate = discordant / (concordant + discordant) for individuals where both haplotypes are classified.

**Key metrics.** 9 of 35 arm/community pairs show >20% discordance rate:

| Arm | Community | Discordant / Total | Rate |
|-----|-----------|---------------------|------|
| chr22_q | C6 | 103 / 217 | 47.5% |
| chr22_p | C7 | 25 / 54 | 46.3% |
| chr9_q | C3 | 98 / 217 | 45.2% |
| chr4_q | C1 | 62 / 143 | 43.4% |
| chr6_p | C5 | 78 / 187 | 41.7% |
| chr11_p | C3 | 81 / 216 | 37.5% |
| chr15_p | C7 | 6 / 22 | 27.3% |
| chr16_q | C3 | 57 / 219 | 26.0% |
| chr14_p | C7 | 13 / 55 | 23.6% |

**Result.** The chr22_q (C6) discordance of 47.5% means nearly half of individuals carry two structurally different subtelomeric types at this locus. The chr4_q rate (43.4%) is higher than the "~20%" reported in a Dutch population (Mefford & Trask 2002, citing van Deutekom et al. 1996). Three C3 arms (chr9_q 45.2%, chr11_p 37.5%, chr16_q 26.0%) show high discordance, reflecting the extensive polymorphism in the f7501 cluster.

**Conclusion.** High discordance rates indicate subtelomeric exchange is a segregating polymorphism. The chr4_q discrepancy (43.4% vs ~20%) reflects the different measurement: the present analysis detects subtelomeric sequence affinity across the broader PHR region, which is more inclusive than full D4Z4 array translocation as defined by Southern blot or PFGE.

### Region length differences

**What it does.** Tests whether cross-arm and self-arm sequences differ in PHR region length.

**How.** Wilcoxon rank-sum test, 18 arm/community pairs, BH-corrected.

**Result.** Cross-arm sequences often differ significantly in length from self-arm sequences (14 of 18 pairs significant at p_adj < 0.05):

| Arm | Community | Cross median | Self median | p-adjusted | Direction |
|-----|-----------|-------------|-------------|------------|-----------|
| chr22_q | C6 | 20 kb | 25 kb | 1.3e-84 | Cross shorter |
| chr6_p | C5 | 100 kb | 153 kb | 1.4e-69 | Cross shorter |
| chr16_q | C3 | 210 kb | 155 kb | 1.9e-46 | Cross longer |
| chr9_q | C3 | 305 kb | 225 kb | 5.1e-34 | Cross longer |
| chrX_p | C15 | 500 kb | 20 kb | 1.3e-26 | Cross longer |
| chr13_q | C6 | 105 kb | 10 kb | 3.7e-17 | Cross longer |
| chr11_p | C3 | 220 kb | 25 kb | 3.7e-11 | Cross longer |
| chr19_p | C3 | 120 kb | 205 kb | 9.3e-8 | Cross shorter |
| chr22_p | C7 | 500 kb | 445 kb | 4.6e-8 | Cross longer |
| chr6_q | C11 | 245 kb | 145 kb | 3.2e-7 | Cross longer |
| chrY_p | C15 | 493 kb | 500 kb | 1.6e-6 | Cross shorter |
| chr19_q | C6 | 20 kb | 20 kb | 4.1e-4 | Means differ |
| chr4_q | C1 | 150 kb | 160 kb | 0.014 | Cross shorter |
| chr10_q | C1 | 145 kb | 135 kb | 0.028 | Cross longer |

**Conclusion.** Length differences confirm structurally distinct subtelomeric haplotype classes. Cross-arm sequences at chr6_p are shorter (100 kb vs 153 kb), while cross-arm sequences at chr9_q are longer (305 kb vs 225 kb), consistent with gain of foreign duplicon content in exchange events. chrX_p (C15) shows the most extreme difference: cross-arm sequences (500 kb) are 25x longer than self (20 kb), because nearly all chrX_p sequences cluster with chrY_p (PAR1 sharing) and the rare "self" haplotypes are truncated. Mefford & Trask (2002, citing Wilkie et al. 1991) reported a "260 kb size difference between the largest and smallest chr16_p subtelomere alleles," with "the longer form accounting for ~30% of alleles."

### Gene repertoire replacement scoring

**What it does.** Quantifies the fraction of a cross-arm sequence's gene content that matches its affinity arm rather than its own arm. A score of 1.0 means complete gene repertoire overlap.

**How.** For each cross-arm sequence, the fraction of genes matching the affinity arm vs own arm is computed (see mechanism caveat, the heterogeneity section intro).

| Community | Own arm | Affinity arm | Score | N affinity-specific genes |
|-----------|---------|-------------|-------|--------------------------|
| C7 | chr13_p | chr22_p | 1.000 | 49 |
| C7 | chr13_p | chr14_p | 1.000 | 46 |
| C7 | chr13_p | chr21_p | 1.000 | 43 |
| C7 | chr13_p | chr15_p | 1.000 | 39 |
| C15 | chrX_p | chrY_p | 1.000 | 10 |
| C14 | chrX_q | chrY_q | 1.000 | 22 |
| C11 | chr6_q | chr1_p | 0.724 | 35 |
| C7 | chr21_p | chr14_p | 0.720 | 41 |
| C3 | chr11_p | chr19_p | 0.698 | 75 |
| C7 | chr21_p | chr15_p | 0.694 | 33 |
| C7 | chr21_p | chr22_p | 0.687 | 32 |
| C3 | chr11_p | chr3_q | 0.677 | 63 |
| C3 | chr11_p | chr7_p | 0.667 | 61 |
| C3 | chr16_q | chr9_q | 0.635 | 53 |
| C7 | chr15_p | chr21_p | 0.594 | 19 |
| C7 | chr15_p | chr14_p | 0.589 | 18 |
| C7 | chr14_p | chr22_p | 0.555 | 14 |
| C7 | chr22_p | chr14_p | 0.562 | 13 |
| C1 | chr4_q | chr10_q | 0.500 | 3 |
| C1 | chr10_q | chr4_q | 0.494 | 1 |
| C6 | chr22_q | chr21_q | 0.333 | 0 |

**Result.** chr13_p achieves perfect conversion scores (1.000) toward all four other acrocentric p-arms — its cross-arm sequences share all genes with their affinity arm. PAR regions also show complete gene repertoire overlap: chrX_p→chrY_p (1.000, 10 genes) and chrX_q→chrY_q (1.000, 22 genes). The f7501 cluster shows high conversion: chr11_p→chr19_p (0.698, 75 affinity-specific genes), chr11_p→chr3_q (0.677, 63 genes). D4Z4 (C1) shows balanced exchange: chr4_q→chr10_q (0.500) and chr10_q→chr4_q (0.494).

**Conclusion.** The conversion scores reflect gene content overlap between cross-arm sequences and their affinity arm, with the highest values in the acrocentric p-arms and the C11 complex (chr1_p/chr5_q/chr6_q/chr8_p). The balanced scores at D4Z4 are consistent with bidirectional exchange maintaining shared gene content between chr4_q and chr10_q.

### TAR1 in cross-arm vs self-arm sequences

**What it does.** Tests whether TAR1 prevalence differs between cross-arm and self-arm sequences.

**How.** Fisher's exact test for prevalence, Wilcoxon for length; 19 arm/community pairs, BH-corrected.

**Result.** TAR1 prevalence is uniformly high (>85%) in both categories for most communities. Three pairs show significant TAR1 prevalence difference (p_adj < 0.05):

| Community | Arm | Cross TAR1% | Self TAR1% | p_adj | Pattern |
|-----------|-----|-------------|------------|-------|---------|
| C3 | chr16_q | 100% (86/86) | 71.1% (258/363) | 2.4e-10 | Cross higher |
| C7 | chr13_p | 88.1% (59/67) | 0% (0/9) | 1.6e-6 | Cross higher |
| C7 | chr14_p | 72.1% (137/190) | 94.9% (37/39) | 0.010 | Cross lower |

C15 (PAR1) is essentially TAR1-free in both cross and self: chrX_p cross 0.3% (1/305), chrY_p self 1.2% (1/82).

**Conclusion.** TAR1 is maintained regardless of exchange history in most communities. The C3/chr16_q pattern (100% cross vs 71% self) is consistent with cross-arm chr16_q sequences — which cluster with f7501 arms — carrying TAR1, while many self-arm chr16_q sequences lack it. A mechanistic role for TAR1 in facilitating exchange cannot be established from prevalence data alone.

### Telomere-proximal vs telomere-distal island distribution by exchange status

**What it does.** Tests whether cross-arm exchange shifts the positional distribution of internal (TTAGGG)n islands.

**How.** Islands classified as proximal (centromere-ward half of PHR) or distal (telomere-ward half). Chi-squared test for cross-arm vs self-arm distribution.

**Key metrics.** Cross-arm: 45.9% proximal (1,494/3,256). Self-arm: 45.9% proximal (6,922/15,096). Chi-squared = 0.00, p = 0.99, OR = 1.00.

**Result.** No difference in island positional distribution between cross-arm and self-arm sequences. The proximal/distal ratio is identical (45.9%) regardless of exchange status.

**Conclusion.** Cross-arm exchange does not shift the positional distribution of internal (TTAGGG)n islands within the PHR. The island distribution reflects the underlying duplicon architecture rather than exchange history.

### Within-community Jaccard distance structure

**What it does.** Tests whether within-community pairwise distances show bimodal structure (allelic vs paralog peaks).

**How.** Pairwise Jaccard distances sampled (up to 4,950 pairs per community) and examined for multi-modal structure.

**Key metrics.** C2: peaks near 0.02 and 0.37 (median allelic 0.023, median paralog 0.365). C12: peaks near 0.04 and 0.73 (median allelic 0.037, median paralog 0.727). Low-distance peak (0.01–0.05) = allelic variation; high-distance peaks = inter-arm paralog distances.

**Result.** Communities with well-separated bimodal peaks (C2, C12) have distinct arm-specific sequence types. Communities with diffuse distributions (C1, C7) show blurring of arm identity consistent with active homogenization.

**Conclusion.** Qualitatively consistent with Ambrosini et al.'s (2007) bimodal identity distribution (peaks at 98% and 91%), though Jaccard distances differ from sequence identity due to the k-mer-based metric.

### Two-domain subtelomeric model test

**What it does.** Tests the Flint et al. (1997) / Mefford & Trask (2002) two-domain model: a distal domain (near telomere) with blocks shared across many chromosome ends, and a proximal domain (centromere-ward) with blocks shared with few ends, separated by internal (TTAGGG)n tracts.

**How.** Three increasingly stringent tests applied to per-window output (5 kb windows, 5 kb step, 500 kb max; windows across 18,827 haplotypes). Focus arms: chr4p, chr4q, chr16p, chr18p, chr20p, chr22q (Flint originals + Mefford candidates).

**Test 1 — Monotonic gradient (Spearman)**: For each arm, all per-window n_chrs values were pooled across haplotypes and correlated with dist_from_telomere. **39/48 arms** show significant negative correlation (p < 0.05), including all 6 focus arms (chr4p rho = −0.85, chr4q −0.68, chr16p −0.32, chr18p −0.31, chr20p −0.47, chr22q −0.95). Per-sequence analysis confirms the gradient is nearly universal: **99.7%** of individual haplotype sequences with ≥5 windows (13,728/13,763) have negative within-sequence rho (median rho = −0.79). The remaining 9 arms break down as: 5 arms with ≤2 windows per sequence and trivially zero rho (chr2_p, chr3_p, chr8_q, chr11_q, chr14_q — sharing confined to a single window, so no gradient can manifest); 2 arms with non-significant rho (chr5_p p=0.32, chrY_p p=0.50); and 2 arms with significant positive correlation (chr18_q rho=+0.56, chrX_p rho=+0.04). chr18_q's positive correlation is consistent with its known atypical subtelomeric structure.

**Test 2 — Discrete breakpoint (piecewise linear regression)**: To distinguish the two-domain prediction (discrete boundary) from smooth exponential decay, a two-segment piecewise linear model was compared to a single linear model for each arm's mean n_chrs vs dist_from_telomere curve, using F-test for significance. **39/41 testable arms** show significantly better fit with the two-segment model (p < 0.05). The two non-significant arms are chr12q (p = 0.48; already well-fit by a single line, R² = 0.65) and chr19p (p = 0.10). Seven arms were excluded for having < 5 unique distance bins. Focus arm breakpoints (single → two-segment R²):

| Arm | Breakpoint (kb) | R² single | R² two-segment | R² improvement |
|-----|-----------------|-----------|----------------|----------------|
| chr4p | 70 | 0.13 | 0.81 | +0.68 |
| chr4q | 50 | 0.14 | 0.92 | +0.77 |
| chr16p | 295 | 0.44 | 0.64 | +0.20 |
| chr18p | 120 | 0.05 | 0.30 | +0.24 |
| chr20p | 165 | 0.24 | 0.47 | +0.23 |
| chr22q | 15 | 0.81 | 1.00 | +0.19 |

Breakpoint positions vary across arms (10–445 kb), consistent with Mefford's observation that the boundary between domains is arm-specific rather than at a fixed distance.

**Test 3 — TTAGGG island co-localization**: The two-domain model predicts that internal (TTAGGG)n tracts mark the domain boundary. RepeatMasker annotations within subtelomeric flanks were extracted from 467 per-sample BED files (29,274 total annotations: 23,636 TAR1, 3,722 (ACCCTA)n, 1,647 (TTAGGG)n). TAR1 blocks cluster at the telomeric tip (median 0.0–2.0 kb from telomere) — these are subtelomeric satellite, not domain separators. The relevant features are internal (TTAGGG)n/(ACCCTA)n blocks > 5 kb from the telomere (3,426 blocks across 19 of 41 testable arms).

For the 19 arms with internal TTAGGG blocks, the closest block to the breakpoint is within 25 kb for 11/19 arms and within 50 kb for 16/19 arms:

| Arm | Breakpoint (kb) | Closest ITS (kb from breakpoint) | N blocks |
|-----|-----------------|----------------------------------|----------|
| chr8p | 145 | 0.1 | 220 |
| chr20p* | 165 | 0.6 | 443 |
| chr13q | 90 | 1.5 | 3 |
| chr7p | 120 | 3.5 | 839 |
| chr16p* | 295 | 17.7 | 121 |
| chr19p | 215 | 18.8 | 255 |
| chr18p* | 120 | 20.3 | 883 |
| chr9p | 45 | 22.9 | 19 |
| chr20q | 135 | 30.8 | 210 |

(*) Mefford/Flint focus arm.

The overall Spearman correlation between breakpoint position and median ITS position is marginal (rho = 0.42, p = 0.08) because the median is diluted by scattered ITS throughout the subtelomeric zone. Not all arms have detectable internal TTAGGG blocks (only 19/41), consistent with either boundary erosion or domain boundaries defined by other features. The co-localization is strongest at specific arms (chr20p: 0.6 kb; chr8p: 0.1 kb) and absent at others, suggesting that internal (TTAGGG)n tracts mark the domain boundary on some but not all chromosome ends.

**Conclusion.** The pangenome-scale data supports the Flint/Mefford two-domain model at all three levels: (1) inter-chromosomal sharing decreases with distance from the telomere on 39/48 arms; (2) this decrease follows a two-phase pattern on 39/41 testable arms; (3) internal (TTAGGG)n blocks co-localize with the domain boundary within 25 kb on 11/19 testable arms. The gradient is present in 99.7% of individual haplotype sequences. Breakpoint positions are arm-specific (10–445 kb), extending the model from the handful of arms originally characterized by Flint et al. (1997) to the entire human chromosome complement.

### Files and scripts

**Scripts:**
| Script | Description |
|--------|-------------|
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze-within-arm-heterogeneity.R` | Cross-arm affinity, population bias, type discordance, gene replacement |
| `/moosefs/guarracino/HPRCv2/scripts/community/plot-within-arm-heterogeneity.R` | Heterogeneity visualization |
| `/moosefs/guarracino/HPRCv2/scripts/community/allele_vs_paralog_distance.R` | Allelic vs paralogous distance comparison |
| `/moosefs/guarracino/HPRCv2/scripts/community/compute_fst_superpop.py` | Pairwise Fst between superpopulations |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/analyze_polymorphic_arms.py` | Polymorphic arm subgroup characterization |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/test_d4z4_causality.py` | D4Z4 causality tests for C1 |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/test_two_domain.py` | Two-domain gradient test (Spearman) |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/test_two_domain_changepoint.py` | Two-domain changepoint + TTAGGG co-localization |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/extract_rm_ttaggg_subtelomeric.py` | RepeatMasker telomeric annotations |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/compute_its_breakpoint_coloc.py` | ITS-breakpoint co-localization |

**Input data:**
| File | Description |
|------|-------------|
| `/moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv` | HPRCv2 sample metadata (superpopulations) |

**Output files:**
| File | Description |
|------|-------------|
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/allele_vs_paralog_distance.tsv` | Allele vs paralog distance per community |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_affinity_sequences.tsv` | Per-sequence cross-arm classification |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_type_discordance.tsv` | Per-individual type discordance |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/gene_conversion_scores.tsv` | Gene replacement scores |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/fst_superpop_matrix.tsv` | Pairwise Fst (5 superpopulations) |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/two_domain_test.tsv` | Per-arm Spearman gradient (48 arms) |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/two_domain_changepoint.tsv` | Per-arm breakpoint analysis (41 arms) |


---

