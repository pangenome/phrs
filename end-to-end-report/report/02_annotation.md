## Annotation preprocessing

**What it does.** Intersects gene annotations (Liftoff GFF3) and TAR1 repeat entries (RepeatMasker) with the subtelomeric similarity regions to characterize gene and repeat content per community.

**How.** 464 haplotype-specific GFF3 files + CHM13 reference processed (465 total). The remote HPRC annotation index covers 462 haplotypes; HG002 (2 haplotypes) was annotated separately using JHU Liftoff v0.6 on the HG002v1.1 assembly. TAR1 entries extracted from RepeatMasker BED files, with CHM13 converted to PanSN naming and HG002 converted from bigBed to 10-column BED.

**Key metrics.** 18,827 total sequences → 15,668 retained (3,158 with no inter-chromosomal signal removed; chimeric chr18_q already excluded in the inter-chromosomal detection section).

**Result.**
- **Gene annotations**: 173,881 gene annotations extracted (374 unique genes) across 233 samples and 39 arms. Two arms — chr7_q and chr12_q — have zero gene annotations because their similarity regions are confined to 5–25 kb at the telomeric tip, below the detection threshold of Liftoff annotation pipelines.
- **TAR1 repeats**: 21,544 TAR1 entries across 14,816 sequences (94.6%) and all 41 arms.

**Conclusion.** The gene and repeat annotations provide the basis for characterizing the biological content of each community (the gene enrichment section) and testing whether exchange status affects gene repertoire (the heterogeneity section).

### TAR1 prevalence

**What it does.** Quantifies TAR1 (Telomere-Associated Repeat 1; Brown et al. 1990) prevalence across chromosome arms.

**How.** RepeatMasker TAR1 entries intersected with PHR regions; prevalence and density computed per arm.

**Key metrics.** Present in all 41 arms. Most autosomal arms: >99%. Near-absent: chrX_p (1/327, 0.3%), chrY_p (1/92, 1.1%). Acrocentric p-arms intermediate: chr15_p 73.1%, chr14_p 76.0%, chr13_p 77.6%, chr21_p 78.4%, chr22_p 78.9%. Highest TAR1 density: chr18_p (mean 4.00 copies/sequence), chr9_q (mean 2.82 copies/sequence).

**Result.** TAR1 is near-universal outside PAR (94.6% of all sequences). The near-absence from PAR1 (chrX_p/chrY_p) is consistent with pseudoautosomal regions using obligate meiotic crossover (Rouyer et al. 1986) rather than repeat-mediated exchange.

**Conclusion.** TAR1's near-universal presence makes it difficult to assess whether it plays a functional role in inter-chromosomal exchange or is a co-located passenger (the testable predictions section, prediction 3).

### TAR1 positional distribution

**What it does.** Maps where TAR1 sits within the PHR region, measuring distance from the telomeric end.

**Key metrics.** 66.9% of 21,544 TAR1 entries within 10 kb of the telomere; 70.3% within 25 kb (median 0.3 kb, mean 44.1 kb). TAR1 is overwhelmingly telomere-proximal. Most arms (both p and q) have median TAR1 distance <1 kb from the telomere. Exceptions with deeper TAR1: acrocentric p-arms (chr13_p–chr22_p: median 179–196 kb, reflecting large PHR regions where TAR1 sits near the telomeric end of a long duplicated zone), chr9_q (median 223.3 kb), chr8_p (140.4 kb), chr18_p (87.3 kb), chr19_p (68.6 kb), chr6_p (56.7 kb).

**Result.** TAR1 is a telomere-proximal element. In arms with short PHR regions (<25 kb), TAR1 sits at the very tip. In arms with long PHR regions (>100 kb), TAR1 remains near the telomere while the duplicated zone extends centromere-ward beyond it.

**Conclusion.** TAR1 positional distribution is consistent with it being a telomere-associated satellite that marks the boundary between terminal telomeric sequence and the subtelomeric duplicated zone.

### TAR1 per community

**What it does.** Quantifies TAR1 prevalence per community.

| Community | Sequences | With TAR1 | % | Mean count | Mean length (bp) |
|-----------|-----------|-----------|---|------------|-------------------|
| C15 (PAR1) | 419 | 2 | 0.5 | 1.50 | 1,412 |
| C7 (acrocentric p) | 763 | 587 | 76.9 | 1.28 | 1,436 |
| C4 (chr7_q/chr12_q) | 895 | 820 | 91.6 | 1.00 | 798 |
| C10 (chr17_p) | 448 | 410 | 91.5 | 1.04 | 1,014 |
| C3 (f7501) | 2,655 | 2,539 | 95.6 | 1.66 | 2,142 |
| C9 (chr16_p) | 442 | 428 | 96.8 | 1.06 | 1,298 |
| C8 (chr15_q) | 449 | 448 | 99.8 | 1.00 | 1,339 |
| C2 (chr10_p/chr18_p) | 889 | 887 | 99.8 | 2.51 | 1,962 |
| All others | — | — | 99.1–100 | 1.06–1.86 | 1,401–2,806 |

**Result.** C15 (PAR1) is essentially TAR1-free (0.5%). C7 (acrocentric p-arms) has the lowest non-PAR prevalence (76.9%), reflecting lower TAR1 density in rDNA-adjacent regions. C2 has the highest TAR1 density (mean 2.51 copies), consistent with chr18_p having the highest per-arm density (4.0 copies/sequence).

### Internal (TTAGGG)n islands

**What it does.** Identifies and quantifies short degenerate telomeric repeat tracts inside subtelomeric regions, distinct from the terminal telomeric array and from TAR1 satellite blocks.

**How.** Canonical (TTAGGG) and variant (TGAGGG, TCAGGG, TTGGGG) telomeric motifs searched across all 15,668 PHR sequences using `seqkit locate`; overlapping hits merged (12 bp gap tolerance) and filtered to 50–1000 bp. Ambrosini et al. (2007) identified these elements (150–823 bp, mostly 150–200 bp) at "duplicon boundaries inside subtelomeric regions."

**Key metrics.** 18,352 islands across 8,321 sequences (53.1%) and all 41 arms. Island lengths: median 79 bp, mean 102 bp. Most sequences carry 1–4 islands; maximum 22 per sequence.

**Result.** Highest island counts per arm: chr20_q (1,765 islands), chr12_q (1,149), chr16_p (898), chr18_p (851). The per-arm positional distribution shows telomere-proximal bias for most p-arms (16/21 with median < 500 bp from telomere). Five p-arms with larger PHR regions have deeper islands (chr11_p 153 kb, chr6_p 99 kb, chr20_p 65 kb, chr18_p 54 kb, chr16_p 1.1 kb). This is consistent with Ambrosini et al.'s (2007) observation that internal (TTAGGG)n-like sequences "almost always co-localize to duplicon boundaries" — islands mark past duplication breakpoints, so arms with more layered duplication history have islands distributed deeper into the subtelomeric zone.

**Conclusion.** Island lengths (median 79 bp) are shorter than Ambrosini et al.'s (2007) observation that "most [islands are] in the 150-200 bp range" (max 823 bp). The difference is methodological: our detection uses degenerate motif search (canonical + 3 variant hexamers + mixed patterns) with a 50 bp minimum across 15,668 sequences from 465 haplotypes, capturing shorter and more degenerate tracts that were not detected in Ambrosini's single-reference analysis.

### TTAGGG island boundary enrichment test

**What it does.** Tests whether internal (TTAGGG)n islands specifically mark the PHR boundary (the transition between duplicated and unique sequence).

**How.** Island positions compared against the centromere-ward boundary of each PHR region. KS test for non-uniformity; binomial test for boundary-proximal enrichment.

**Key metrics.** KS test stat=0.37, p < 1e-300 (non-uniform). Mean fractional distance = 0.54 (islands sit on average 54% of the way through the PHR region). 42.5% of islands within 5 kb of boundary (7,806/18,352).

**Result.** Islands are distributed throughout the PHR region (mean fractional position 0.54), with a slight enrichment near the boundary (42.5% within 5 kb). The KS test rejects uniformity (p < 1e-300) but the distribution is not strongly polarized toward either end.

**Conclusion.** This test evaluated the PHR outer boundary (the single transition from duplicated to unique sequence). Ambrosini et al.'s (2007) claim concerns co-localization with internal duplicon-to-duplicon boundaries throughout the subtelomeric zone — a different feature not tested here. A direct test would require mapping individual duplicon block boundaries within each PHR region.

### TTAGGG island length distribution

**What it does.** Characterizes the size distribution of internal (TTAGGG)n islands at population scale.

**Key metrics.** Mode at 50–74 bp (8,433 islands, 46.0%); monotonically decreasing at longer lengths: 75–99 bp (4,109, 22.4%), 100–149 bp (3,047, 16.6%), 150–199 bp (1,455, 7.9%), 200–299 bp (838, 4.6%), 300–499 bp (374, 2.0%), 500–1000 bp (96, 0.5%).

**Result.** Ambrosini et al. (2007) reported a mode at 150–200 bp; this subrange contains only 7.9% of population-scale islands, with no secondary mode.

**Conclusion.** The 50 bp minimum filter captures shorter degenerate tracts below Ambrosini et al.'s detection threshold, accounting for the shorter median (79 bp vs 150–200 bp) — a methodological, not biological, difference.

### TTAGGG island motif composition

**What it does.** Quantifies the canonical vs variant telomeric hexamer content of internal islands.

**Key metrics.** From 18,352-island dataset (296,406 total hexamer instances). Canonical TTAGGG: 52.3% of hexamer instances (154,886/296,406). Variants: TGAGGG 19.0%, TTGGGG 16.0%, TCAGGG 12.7%. Only 32.2% of islands are "pure canonical" (≥80% TTAGGG+CCCTAA); 47.2% are variant-dominant (<50% canonical).

**Result.** The high variant content — consistent with Ambrosini et al.'s (2007) identification of the same three dominant variant motifs — indicates these islands are substantially degenerate relative to the terminal telomeric array. Linardopoulou et al. (2005) found degenerate telomeric repeats enriched at "4% of subtelomeric DSB sites" (vs 0.5% background), suggesting these islands may have been appended during DSB repair.

**Conclusion.** The degenerate composition supports interpretation of internal (TTAGGG)n islands as ancient remnants of telomeric sequence incorporated during subtelomeric duplication events.

### TTAGGG island count by cross-arm status

**What it does.** Tests whether cross-arm exchange status affects the number of internal telomeric islands.

**Key metrics.** 8,321 sequences with islands: 1,569 (18.9%) cross-arm, 6,752 (81.1%) self-arm. Mean island count: cross-arm 2.08, self-arm 2.24. Mann-Whitney U test z = −1.89, p = 0.045.

**Result.** Marginally significant difference (p=0.045): self-arm sequences carry slightly more islands on average (2.24 vs 2.08), though the effect size is small.

**Conclusion.** Cross-arm exchange status does not affect internal telomeric island count.

### Terminal telomere tract length by community

**What it does.** Tests whether terminal telomere length varies across communities.

**How.** Terminal telomere repeat tract lengths (from `.telo.tsv`) matched to arm-level community assignments.

**Key metrics.** Kruskal-Wallis H = 100.89, p = 3.2e-15. Medians range from 7,638 bp (C10, chr17_p) to 9,418 bp (C13, chr4_p). Correlation between terminal telomere tract length and TTAGGG island count (restricted to the 8,321 sequences with at least one island): Spearman rho = −0.056 (p = 2.7e-7). Telomere lengths: min 470 bp, max 33,826 bp.

**Result.** Telomere lengths vary significantly across communities. C13 (chr4_p) and C15 (PAR1) have the longest telomeres (mean 9,640 bp and 9,266 bp); C9 (chr16_p) and C10 (chr17_p) the shortest (mean 8,360 bp and 8,178 bp). The correlation between telomere length and TTAGGG island count is weakly negative (rho = −0.056, p = 2.7e-7) — sequences with longer telomeres tend to have slightly fewer internal (TTAGGG)n islands, but the effect size is minimal.

**Conclusion.** Community membership is associated with terminal telomere length variation, but the biological significance of this association is unclear.

### Files and scripts

**Scripts:**
| Script | Description |
|--------|-------------|
| `/moosefs/guarracino/HPRCv2/scripts/preprocessing/preprocess-subtelomeric-annotations.R` | Gene/TAR1 annotation intersection with PHR regions |
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze-tar1-positional.R` | TAR1 positional distribution analysis |
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze-ttaggg-islands.py` | TTAGGG island boundary enrichment, length, cross-arm count |
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze-ttaggg-motifs.py` | TTAGGG island hexamer motif composition |
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze-island-exchange-status.py` | TTAGGG island distribution by exchange status |
| `/moosefs/guarracino/HPRCv2/scripts/community/ttaggg_boundary_enrichment.py` | TTAGGG boundary enrichment statistical testing |
| `/moosefs/guarracino/HPRCv2/scripts/community/telomere_length_by_community.py` | Terminal telomere tract length per community |

**Input data:**
| File | Description |
|------|-------------|
| `/moosefs/guarracino/HPRCv2/PHR_III/hprc_annotations/*.gff3.gz` | 464 haplotype-specific gene annotations |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprc_repeatmasker/*.RepeatMasker.bed.gz` | 465 haplotype-specific repeat annotations |

**Output files:**
| File | Description |
|------|-------------|
| `/moosefs/guarracino/HPRCv2/PHR_III/annotations/subtelomeric_annotations.1Mb.rds` | Gene + TAR1 annotations intersected with PHR regions |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_tar1_summary.tsv` | TAR1 prevalence per community |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_tar1_comparison.tsv` | TAR1 in cross-arm vs self-arm sequences |


