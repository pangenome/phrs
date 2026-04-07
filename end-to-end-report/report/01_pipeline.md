## Contig classification and telomere filtering

**What it does.** Classifies each assembled contig as p-arm, q-arm, or spanning (pq) and filters for telomere presence.

**How.** `classify_contigs.py` classifies each contig based on its alignment to CHM13 p-arm and q-arm coordinates, using only alignments where the query contig's chromosome of origin (via chromosome alias) matches the target CHM13 chromosome:
- **p-contig**: aligns only to the p-arm of its corresponding chromosome.
- **q-contig**: aligns only to the q-arm.
- **pq-contig**: spans both arms (T2T or near-T2T contig).
- **unclassified**: no alignments to p-arms or q-arms.

**Key metrics.** Input: 465 PAF mapping files. Filters: minimum contig length ≥ 1 Mb; telomere presence required (p-contigs need a p-telomere, q-contigs a q-telomere, pq-contigs at least one); no more than one telomere of each type per contig.

**Result.** 12,649 classified contigs (9,557 pq-contigs spanning both arms, 1,598 q-contigs, 1,494 p-contigs) in `pq-classification/contig_classifications.tsv`. Of these, 12,635 pass full validation; 14 have mixed-strand alignments to CHM13, so their p/q telomere assignment is less reliable, but they are used in downstream steps.

**Conclusion.** These classifications define which contigs have a telomere-anchored arm assignment and contribute to downstream subtelomeric analysis.

---

## Subtelomeric flank extraction

**What it does.** Extracts 500 kb of subtelomeric sequence internal to each telomere from all classified contigs.

**How.** For each classified contig with a telomere, extract the telomere plus 500 kb of adjacent subtelomeric sequence (inward from the telomere), then trim the telomere repeat tract itself. The result is exactly 500 kb of subtelomeric sequence internal to each telomere (the MIN_LEN ≥ 1 Mb filter ensures every contig is long enough). Each p-arm and q-arm telomere generates one flank; pq-contigs (spanning both arms) contribute two. GRCh38 contigs (no telomeres) and CHM13#0#chrY (masked PAR1) are excluded.

**Key metrics.** FLANK = 500 kb (motivated by Ambrosini et al. 2007, who stated that subtelomeric segmental duplications comprise "about 25% of the most distal 500 kb"); MIN_LEN = 1 Mb.

**Result.** 18,827 subtelomeric flank sequences across all 48 chromosome arms (24 chromosomes × 2 arms) in `hprcv2.1Mb.telo_500kb_trimmed.fa.gz`. Per-arm counts, sorted by chromosome:

| Arm | Count | | Arm | Count | | Arm | Count | | Arm | Count |
|-----|------:|-|-----|------:|-|-----|------:|-|-----|------:|
| chr1_p | 439 | | chr7_p | 438 | | chr13_p | 76 | | chr19_p | 434 |
| chr1_q | 445 | | chr7_q | 447 | | chr13_q | 449 | | chr19_q | 457 |
| chr2_p | 438 | | chr8_p | 458 | | chr14_p | 230 | | chr20_p | 453 |
| chr2_q | 450 | | chr8_q | 452 | | chr14_q | 390 | | chr20_q | 441 |
| chr3_p | 453 | | chr9_p | 448 | | chr15_p | 119 | | chr21_p | 116 |
| chr3_q | 436 | | chr9_q | 450 | | chr15_q | 449 | | chr21_q | 449 |
| chr4_p | 448 | | chr10_p | 443 | | chr16_p | 442 | | chr22_p | 223 |
| chr4_q | 357 | | chr10_q | 358 | | chr16_q | 449 | | chr22_q | 450 |
| chr5_p | 451 | | chr11_p | 448 | | chr17_p | 448 | | chrX_p | 327 |
| chr5_q | 447 | | chr11_q | 451 | | chr17_q | 451 | | chrX_q | 330 |
| chr6_p | 447 | | chr12_p | 454 | | chr18_p | 446 | | chrY_p | 92 |
| chr6_q | 440 | | chr12_q | 450 | | chr18_q | 457 | | chrY_q | 101 |


**Conclusion.** The 500 kb window captures the full subtelomeric duplicated zone. Counts range from 76 (chr13_p) to 458 (chr8_p). Lower counts at acrocentric p-arms and sex chromosomes reflect assembly difficulty in these regions.

---

## All-vs-all subtelomeric alignment

**What it does.** Aligns every subtelomeric flank against every other flank to detect shared sequence at ≥95% identity.

**How.** wfmash v0.23.0-41-gb5f0ff1c (`-p 95 -t 48 --quiet`): each of the 18,827 sequences serves as target in turn, with all 18,827 sequences as queries. The 95% identity threshold was chosen to capture the high-identity peak of Ambrosini et al.'s (2007) bimodal duplicon identity distribution (peaks at 91% and 98%), accepting that older exchanges (≤91% identity) are missed (see the limitations section, limitation 1).

**Key metrics.** Identity threshold: 95%. Output: 18,827 PAF.gz files.

**Result.** 18,827 PAF.gz files in `all-vs-all.1Mb.p95/`, indexed with impg for fast downstream queries.

**Conclusion.** These alignments provide the raw inter-sequence similarity data for inter-chromosomal region detection (the inter-chromosomal detection section) and pangenome graph construction (the pangenome graph section).

---

## Inter-chromosomal region detection

**What it does.** Identifies, for each subtelomeric flank, the region of inter-chromosomal similarity, that is the portion of the flank that shares sequence with other chromosomes at ≥95% identity.

**How.** `find-multichr-regions-incremental.py` slides windows from the telomere inward along each flank. For each window, it queries impg for alignments and counts how many different chromosomes have matches. Scanning stops early when consecutive windows fail to meet the threshold.

**Key metrics.** Window = 5 kb, step = 5 kb, max distance from telomere = 500 kb, min alignment identity = 0.95, min different chromosomes = 2, min alignment count per chromosome = 5, min consecutive failing windows = 4, min output region length = 3 kb.

**Result.** 18,826 data rows in `all-vs-all.1Mb.p95.id95.len.tsv` (one per flank, excluding the chimeric chr18_q contig — see below). 15,668 sequences (83.2%) have inter-chromosomal matches; 3,158 (16.8%) have none. The 15,668 PHR sequences span 41 of 48 chromosome arms. Similarity region lengths: median 105 kb, mean 144 kb, range 5–500 kb.

**Conclusion.** The majority of subtelomeric flanks share detectable sequence with other chromosomes. Seven arms have no inter-chromosomal signal at these thresholds.

### Chimeric contig exclusion

**What it does.** Identifies and removes a chimeric contig whose inter-chromosomal signal was artifactual.

**How.** The chr18_q flank from NA18982#1 (JBKABS010000018.1, 84.4 Mb) was flagged because its inter-chromosomal signal (chrX, chrY) arose from a scaffolding artifact: the contig fuses chr18 sequence with 966 kb of chrX PAR1 across a 100 bp NNN scaffold join.

**Key metrics.** Junction at query ~83.37–83.38 Mb confirmed by both wfmash and minimap2 v2.30 (mapq 60). Flagger labels the junction as NNN and the chrX portion as Hap (haploid). A 2,826 bp terminal telomeric tract (~471 TTAGGG repeats) precedes the NNN gap. NA18982#1 has no separate chrX contig.

**Result.** One sequence removed; 15,668 PHR sequences retained.

**Conclusion.** The NNN gap indicates a scaffolding artifact, not a contiguous assembly error.

### Zero-signal arms and one-copy regions

**What it does.** Compares the 7 arms with no inter-chromosomal signal against Ambrosini et al.'s (2007) 6 "one-copy DNA (TTAGGG)n-adjacent regions", that is chromosome ends where the sequence immediately adjacent to the terminal telomeric tract "is unique in the genome."

**How.** The 7 zero-signal arms (chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q) were compared to Ambrosini's 6 one-copy regions (chr7_q, chr8_q, chr11_q, chr12_q, chr18_q, chrX_p/chrY_p). For arms with signal (chr7_q, chr12_q, chrX_p/chrY_p), per-arm chromosome matching was computed from the `all-vs-all.1Mb.p95.id95.len.tsv` file.

**Key metrics.**
- 3 one-copy regions (chr8_q, chr11_q, chr18_q) are among the 7 zero-signal arms.
- chrX_p/chrY_p: 416/419 sequences with signal (99.3%) match exclusively chrX/chrY and no other chromosome.
- chr12_q: 449/449 (100%) match exclusively chr7/chr12.
- chr7_q: 424/446 (95.1%) match exclusively chr7/chr12; the 22 outliers match many chromosomes.
- chr7_q shared regions: median 40 kb (range 5–45 kb). chr12_q: median 25 kb (range 15–45 kb).

**Result.** Four of 6 one-copy regions are confirmed: chr8_q, chr11_q, and chr18_q have zero signal; chrX_p/chrY_p shares sequence only within the PAR1 pair. The remaining 2 one-copy regions (chr7_q, chr12_q) form a previously undetected private pair, analogous to chrX_p/chrY_p — Ambrosini listed them as independent one-copy regions, but the population-scale data reveals they share small regions (5–25 kb) exclusively with each other.

**Conclusion.** The population-scale analysis confirms Ambrosini et al.'s one-copy classification for 4 of 6 regions and refines the classification for chr7_q/chr12_q by revealing a private pair invisible in the single reference genome. The chr7_q/chr12_q pair forms community C4 (the community detection section).

---

## Pangenome graph construction and pairwise similarity

**What it does.** Builds a pangenome graph of all PHR sequences and computes pairwise Jaccard similarity between every pair of sequences. The Jaccard similarity measures the fraction of graph nodes (sequence segments) shared between any two sequences.

**How.** `pggb -p 95 -D /scratch` constructs the pangenome graph from the 15,668 PHR sequences. `odgi similarity --all -P` then computes the Jaccard similarity for every pair of sequences. Values range from 0 (no shared nodes) to 1 (identical graph traversals).

**Key metrics.** 15,668 input sequences. Graph construction identity threshold: 95%. Output: 15,668 × 15,668 Jaccard similarity matrix (12 GB compressed, 245 million pairwise entries including both directions).

**Result.** Pangenome graph and full pairwise similarity matrix in `*.similarity.tsv.gz`.

**Conclusion.** The Jaccard similarity matrix provides the distance measure for all downstream community detection (the community detection section) and heterogeneity analyses (the heterogeneity section). Unlike alignment identity, Jaccard similarity captures structural variation (insertions, deletions, inversions) through graph topology.

---

## Community detection

### Arm-level community detection

**What it does.** Groups the 41 chromosome arms with detected PHRs into communities based on their average pairwise subtelomeric sequence similarity. Seven of the 48 human chromosome arms (chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q) are excluded because no inter-chromosomal PHR was detected at these arms (the inter-chromosomal detection section).

**How.** The pairwise Jaccard similarity (the pangenome graph section) is converted to distance (1 − Jaccard) and averaged across all sequence pairs per arm pair to produce a 41 × 41 arm-level distance matrix. Two clustering methods are applied:

- **Leiden** (primary): a graph-based community detection algorithm (Traag et al. 2019) that partitions a fully-connected weighted graph into communities. The Leiden algorithm internally optimizes modularity at each resolution; the optimal resolution is selected by maximizing mean silhouette score across a scan from 0.1 to 3.0 in steps of 0.01. Edge weights are derived from the arm-level distance matrix using exponential decay (w_ij = exp(-d_ij / median(d))). Used throughout the analysis.
- **UPGMA** (comparison): hierarchical agglomerative clustering (`hclust(..., method = "average")`). The dendrogram is cut at the k (2–20) that maximizes silhouette score. Run by `plot-similarity-subtelo.R`.

**Key metrics.** Leiden: 15 communities (optimal resolution = 1.16, silhouette = 0.347). UPGMA: 14 communities (silhouette = 0.342). **Silhouette score** measures how well each element fits its assigned cluster: for element i, silhouette(i) = (b − a) / max(a, b), where a = mean distance to elements in the same cluster and b = mean distance to elements in the nearest other cluster. Range [−1, 1]: positive = well-clustered, 0 = on boundary, negative = likely misassigned.

**Result.** The 15 Leiden communities group the 41 arms as follows:

| Leiden | N arms | Arms | Description |
|--------|-------:|------|-------------|
| C1 | 2 | chr4_q, chr10_q | D4Z4 macrosatellite sharing |
| C2 | 2 | chr10_p, chr18_p | Recurrent transfer pair (Linardopoulou 2005 Fig. 5) |
| C3 | 6 | chr3_q, chr7_p, chr9_q, chr11_p, chr16_q, chr19_p | f7501 sites: fixed (chr3_q, chr19_p) + AFR-enriched (chr7_p, chr16_q) + variable (chr9_q, chr11_p) |
| C4 | 2 | chr7_q, chr12_q | Minimal telomeric tips, no gene annotations |
| C5 | 4 | chr6_p, chr9_p, chr12_p, chr20_q | Shared duplicon modules (RPL23A, WASH, DDX11L families) |
| C6 | 6 | chr1_q, chr13_q, chr17_q, chr19_q, chr21_q, chr22_q | 3 acrocentric q-arms (chr13_q, chr21_q, chr22_q) + chr1_q/chr17_q/chr19_q |
| C7 | 5 | chr13_p, chr14_p, chr15_p, chr21_p, chr22_p | All 5 acrocentric p-arms; rDNA-adjacent homogenization |
| C8 | 1 | chr15_q | Single arm |
| C9 | 1 | chr16_p | Single arm |
| C10 | 1 | chr17_p | Single arm |
| C11 | 4 | chr1_p, chr5_q, chr6_q, chr8_p | OR4F21 sharing (Linardopoulou 2005 block 5); variable f7501 sites |
| C12 | 2 | chr2_q, chr20_p | Well-separated pair; chr2_q is a variable f7501 site |
| C13 | 1 | chr4_p | Single arm |
| C14 | 2 | chrX_q, chrY_q | PAR2 sharing |
| C15 | 2 | chrX_p, chrY_p | PAR1 sharing |

**UPGMA comparison.** UPGMA produces 14 communities (vs 15 for Leiden) with comparable silhouette (0.342 vs 0.347). The two methods agree exactly on 12 of 15 Leiden communities: C1(D4Z4), C2(chr10p/18p), C4(chr7q/12q), C5(RPL23A/WASH), C6(q-arms), C7(acrocentric p), C9(chr16p), C10(chr17p), C11(OR4F21), C13(chr4p), C14(PAR2), C15(PAR1). They differ on 3 Leiden communities around the f7501 cluster: UPGMA merges Leiden C3(6 arms) + C8(chr15_q) + chr20_p from C12 into a single 8-arm community, while isolating chr2_q (the other half of Leiden C12) as a singleton. This confirms that the f7501-carrying arms (C3) form the most cohesive group in the dataset, with chr15_q and chr20_p on its boundary.

**Conclusion.** Both Leiden and UPGMA produce arm-level communities consistent with known biology: D4Z4 sharing (C1), acrocentric p-arm homogenization (C7), PAR1/PAR2 (C14/C15), and f7501 duplicon exchange (C3, C11). The two methods agree exactly on 12 of 15 communities; the 3 differences involve boundary assignments around the f7501 cluster.

### f7501 (L78442) direct alignment validation

**What it does.** Tests whether direct alignment of the f7501 cosmid reproduces Mefford & Trask's (2002) Fig. 3 per-arm distribution at population scale.

**How.** The L78442.1 cosmid (36.3 kb, the original f7501 clone from Mefford & Trask 2002, with 3 regions of similarity to olfactory receptor genes including the expressed *OR-A* gene) was aligned against each of the 18,827 subtelomeric flanks individually (minimap2 `-x asm20`, one-vs-one to avoid multi-mapping). Haplotypes with ≥30 kb matching bases were counted as carrying f7501.

**Key metrics.** Per-arm distribution across 465 haplotypes (233 samples × 2 haplotypes, minus 1 for CHM13 which has a single haplotype). Population enrichment tested by Fisher's exact test (one-sided greater) for each superpopulation independently (AFR=134, AMR=88, EAS=104, EUR=65, SAS=74 haplotypes); the table reports the most significantly enriched superpopulation per arm. The per-arm distribution reproduces Mefford & Trask's Figure 3 at population scale:

| Arm | Haps | % of 465 | Mefford status | Leiden | AFR | AMR | EAS | EUR | SAS | Best enriched | % | OR | p-value |
|-----|-----:|--------:|----------------|--------|----:|----:|----:|----:|----:|:---:|---:|----:|--------:|
| chr3_q | 427 | 91.8% | FIXED | C3 | 128 | 84 | 88 | 58 | 69 | AFR | 30.0 | 2.3 | 4.3e-02 |
| chr19_p | 421 | 90.5% | FIXED | C3 | 126 | 80 | 87 | 57 | 71 | SAS | 16.9 | 2.8 | 5.6e-02 |
| **chr15_q** | **398** | **85.6%** | **FIXED** | **C8** | **87** | **81** | **98** | **64** | **68** | **EUR** | **16.1** | **12.6** | **2.5e-04** |
| chr11_p | 273 | 58.7% | Variable | C3 | 83 | 50 | 66 | 32 | 42 | EAS | 24.2 | 1.3 | 1.6e-01 |
| **chr9_q** | **138** | **29.7%** | **Variable** | **C3** | **59** | **22** | **25** | **19** | **13** | **AFR** | **42.8** | **2.5** | **1.9e-05** |
| **chr16_q** | **85** | **18.3%** | **AFR-enriched** | **C3** | **67** | **6** | **2** | **2** | **8** | **AFR** | **78.8** | **17.4** | **6.6e-27** |
| **chr16_p** | **25** | **5.4%** | **Variable** | **C9** | **19** | **4** | **0** | **1** | **1** | **AFR** | **76.0** | **8.9** | **6.7e-07** |
| **chr2_q** | **22** | **4.7%** | **Variable** | **C12** | **2** | **0** | **3** | **0** | **17** | **SAS** | **77.3** | **23.0** | **6.8e-11** |
| **chr7_p** | **17** | **3.7%** | **AFR-enriched** | **C3** | **10** | **4** | **0** | **2** | **1** | **AFR** | **58.8** | **3.7** | **8.2e-03** |
| **chr6_p** | **13** | **2.8%** | **Variable** | **C5** | **2** | **8** | **0** | **2** | **1** | **AMR** | **61.5** | **7.4** | **7.0e-04** |
| **chr8_p** | **10** | **2.2%** | **Variable** | **C11** | **9** | **0** | **1** | **0** | **0** | **AFR** | **90.0** | **23.8** | **8.5e-05** |
| chr6_q | 5 | 1.1% | Variable | C11 | 3 | 1 | 0 | 1 | 0 | AFR | 60.0 | 3.8 | 1.5e-01 |
| chr5_q | 3 | 0.6% | Variable | C11 | 0 | 1 | 0 | 2 | 0 | EUR | 66.7 | 12.7 | 5.3e-02 |

**Result.** The f7501 distribution confirms Mefford & Trask (2002) and reveals new population-enrichment findings:

- **Fixed sites confirmed**: chr3_q (91.8%), chr19_p (90.5%), chr15_q (85.6%) are present in 86–92% of haplotypes, consistent with colonization before the out-of-Africa dispersal. chr3_q and chr19_p map to Leiden C3 (6 arms, tied with C6 as the largest arm-level community); chr15_q maps to C8 (a singleton community — chr15_q is the only arm in C8).
- **AFR enrichment confirmed**: chr16_q (p=6.6e-27) and chr7_p (p=8.2e-03) are significantly AFR-enriched, confirming Mefford & Trask's observation that "the high frequency of the block on chromosomes 7p and 16q in members of the African pygmy group is notable" and that "the block is rarely seen on these chromosomes in populations that migrated out of Africa." Both map to Leiden C3.
- **Three additional AFR-enriched arms identified**: chr8_p (9/10 = 90% AFR, p=8.5e-05, C11), chr16_p (19/25 = 76% AFR, p=6.7e-07, C9), chr9_q (59/138 = 43% AFR, p=1.9e-05, C3).
- **chr2_q is SAS-enriched**: 17/22 = 77% SAS (p=6.8e-11, C12).
- **chr6_p is AMR-enriched**: 8/13 = 62% AMR (p=7.0e-04, C5).
- **chr15_q refinement**: Mefford classified chr15_q as "FIXED" based on 52 individuals by FISH. At 233 samples, the overall prevalence (85.6%) confirms near-fixation, but with significant EUR enrichment (p=2.5e-04): 98.5% of EUR haplotypes carry f7501 (64/65) vs only 64.9% of AFR (87/134). Non-AFR populations are uniformly high: AMR 92.0%, SAS 91.9%, EAS 94.2%, EUR 98.5%. This is consistent with ongoing f7501 loss in African populations or incomplete lineage sorting at this locus.
- **Novel locations not in Mefford**: chr1_p (5 haplotypes, C11), chr20_p (2, C12), and chrX_q (1, C14) carry f7501 at ≥30 kb but are absent from Mefford's Figure 3 (which used FISH on 52 individuals). These may represent rare or recently acquired f7501 copies below the detection threshold of the original FISH survey.
- **chr4_q and chr19_q absence**: Mefford's Figure 3 lists chr4_q and chr19_q as variable f7501 sites detected by FISH. The L78442 cosmid produces zero alignments ≥30 kb to any chr19_q flank (457 tested) and only ~954 bp partial matches to chr4_q. The f7501 copies at these arms are too divergent for sequence alignment at asm20 stringency, despite being detectable by FISH cross-hybridization.

**Conclusion.** The f7501 distribution across 465 haplotypes reproduces the FISH-based patterns observed by Mefford & Trask (2002) in 52 individuals and extends them with three newly identified AFR-enriched sites, population-specific enrichments at chr2_q (SAS) and chr6_p (AMR), and three novel f7501 locations (chr1_p, chr20_p, chrX_q) not reported by FISH.

---

### Sequence-level community detection

**What it does.** Clusters individual subtelomeric sequences (not arm averages) into communities, resolving within-arm polymorphism invisible at the arm level.

**How.** Leiden community detection on the full 15,668 × 15,668 pairwise Jaccard distance matrix (the pangenome graph section). A k-NN graph is constructed with exponential decay edge weights. A joint scan over k (10, 25, 50, 75, 100, 125) and resolution (0.1 to 3.0 in steps of 0.1) selects the (k, resolution) pair that maximizes modularity within a target community count range (5–50). Optimal: k=75, resolution=0.8. Run by `detect_communities.R`.

**Key metrics.** Leiden: k-NN = 75, 50 communities (modularity = 0.97).

**Result.**

*Community structure.* Of the 50 Leiden communities: 18 are pure (single arm), 23 are near-pure (dominant arm > 90%), and 9 are mixed (dominant arm ≤ 90%). From the arm perspective: chr6_q is the most polymorphic (sequences split across 8 communities), followed by chr19_p (7), chr3_q (6), and chr7_q, chr11_p, chr5_q, chr16_q, chr20_p (5 each).

*Inter-chromosomal sharing.* The mixed communities reveal inter-chromosomal sequence sharing consistent with past exchange:

| Community | Size | Composition | Biological pattern |
|-----------|------|-------------|-------------------|
| C4 | 770 | chr14_p (30%), chr22_p (29%), chr15_p (15%), chr21_p (15%), chr13_p (10%) + 3 minor | All 5 acrocentric p-arms + 3 stray seqs; rDNA-adjacent concerted evolution |
| C3 | 712 | chr4_q (50%), chr10_q (50%) + 2 stray seqs (chr21_p, chr14_p) | D4Z4 macrosatellite sharing |
| C32 | 432 | chrX_q (76%), chrY_q (23%) + 1 stray (chrY_p) | Sex chromosome PAR2 sharing |
| C33 | 416 | chrX_p (78%), chrY_p (22%) | Sex chromosome PAR1 sharing |
| C40 | 352 | chr9_q (39%), chr11_p (30%), chr16_q (17%), chr16_p (7%), chr7_p (5%) + 1 minor | f7501-associated arms; multi-arm sharing |
| C13 | 508 | chr17_p (86%), chr11_p (14%) | chr17_p/chr11_p duplicon sharing |
| C26 | 490 | chr7_p (86%), chr11_p (12%), chr17_p (2%) | chr7_p/chr11_p duplicon sharing |
| C27 | 298 | chr6_q (88%), chr15_q (12%) | chr6_q/chr15_q pair |

**Conclusion.** The sequence-level clustering resolves the 15 arm-level communities into 50 finer-grained communities that capture within-arm polymorphism and inter-chromosomal mixing invisible at the arm level. chr6_q is the most polymorphic arm (8 communities), reflecting diverse subtelomeric content across haplotypes.

### Arm-level vs sequence-level community nesting

**What it does.** Quantifies how arm-level communities fragment into sequence-level communities, distinguishing monolithic from heterogeneous arm communities.

**How.** For each arm-level community, the fraction of sequences assigned to the largest sequence-level community is computed.

**Result.** The 15 arm-level communities fragment unevenly into the 50 sequence-level communities (Leiden k-NN=75, resolution=0.8, modularity=0.97). The most monolithic:

| Arm community | Arms | Largest seq community | Fraction | N seq communities |
|--------------|------|----------------------|----------|-------------------|
| C1 (D4Z4) | chr4_q, chr10_q | C3 | 99.4% | 3 |
| C7 (acro p) | 5 acrocentric p-arms | C4 | 99.6% | 3 |
| C13 (chr4_p) | chr4_p | C24 | 99.8% | 2 |
| C14 (PAR2) | chrX_q, chrY_q | C32 | 99.8% | 2 |
| C15 (PAR1) | chrX_p, chrY_p | C33 | 99.3% | 3 |
| C10 (chr17_p) | chr17_p | C13 | 97.1% | 3 |
| C8 (chr15_q) | chr15_q | C10 | 91.5% | 3 |

Intermediate (neither monolithic nor highly fragmented):

| Arm community | Arms | Largest seq community | Fraction | N seq communities |
|--------------|------|----------------------|----------|-------------------|
| C9 (chr16_p) | chr16_p | C11 | 71.0% | 4 |
| C4 (chr7_q, chr12_q) | chr7_q, chr12_q | C7 | 50.1% | 6 |
| C2 (chr10_p, chr18_p) | chr10_p, chr18_p | C5 | 49.8% | 3 |
| C12 (chr2_q, chr20_p) | chr2_q, chr20_p | C18 | 49.5% | 8 |

The most fragmented:

| Arm community | Arms | Largest seq community | Fraction | N seq communities |
|--------------|------|----------------------|----------|-------------------|
| C3 (f7501) | chr3_q, chr7_p, chr9_q, chr11_p, chr16_q, chr19_p | C26 | 18.0% | 16 |
| C11 (OR4F21) | chr1_p, chr5_q, chr6_q, chr8_p | C25 | 24.7% | 15 |
| C6 (q-arms) | chr1_q, chr13_q, chr17_q, chr19_q, chr21_q, chr22_q | C14 | 16.7% | 10 |
| C5 (RPL23A/WASH) | chr6_p, chr9_p, chr12_p, chr20_q | C8 | 26.3% | 8 |

Summary: 15 arm-level communities, 50 sequence-level communities (32 multi-arm), mean clustering silhouette 0.347 (arm) vs 0.602 (sequence). Partition agreement: ARI=0.35, NMI=0.76.

**Conclusion.** The D4Z4-containing arms (C1: chr4_q, chr10_q) and the acrocentric p-arms (C7) are the most monolithic: >99% of sequences in a single seq-level community, reflecting homogeneous inter-chromosomal sharing. PAR1 (C15) and PAR2 (C14) are similarly cohesive (>99%). The most fragmented arm-level communities (C3/f7501 with 16 seq-communities, C11/OR4F21 with 15) contain multiple chromosome arms each with extensive within-arm polymorphism. The sequence-level clustering (silhouette 0.60) substantially improves over arm-level (0.35), capturing within-arm subtypes invisible at coarser resolution.

### Polymorphic arm subgroup characterization

**What it does.** Classifies each sequence as self-arm or cross-arm and characterizes three patterns of subtelomeric polymorphism across the 41 chromosome arms (34 polymorphic, 7 monomorphic).

**How.** A sequence has **cross-arm affinity** when its subtelomeric sequence content is more similar to sequences from a different chromosome arm than to its own arm — measured by Jaccard distance on shared pangenome graph nodes (the pangenome graph section). Operationally, the sequence clusters into a community dominated by a foreign arm rather than its own. Such cross-arm similarity is attributed to inter-chromosomal exchange (Linardopoulou et al. 2005; Mefford & Trask 2002), though the present data measures similarity, not exchange events directly.

**Key metrics.** 1,740 of 15,668 sequences (11.1%) show cross-arm affinity at the k50 partition. Seven arms are fully absorbed (100% cross-arm): chrY_p, chrY_q, chr15_p, chr21_p, chr13_p, chr22_p, and chr10_q. Eight arms have 0% cross-arm affinity (chr10_p, chr12_p, chr17_q, chr18_p, chr1_q, chr20_q, chr21_q, chr9_q).

**Result.** Three patterns of cross-arm affinity emerge:

**Pattern 1 — Full cross-arm absorption**: Seven arms have no self-community — all their sequences cluster with communities dominated by a foreign arm:
- **chrY_p** (100%): All 92 sequences cluster with chrX_p (PAR1 sharing).
- **chrY_q** (100%): All 101 sequences cluster with chrX_q (PAR2 sharing).
- **chr15_p** (100%), **chr21_p** (100%), **chr13_p** (100%), **chr22_p** (100%): All four acrocentric p-arms cluster entirely with chr14_p-dominated communities, reflecting pervasive rDNA-adjacent sequence homogenization.
- **chr10_q** (100%): All 357 sequences cluster with chr4_q-dominated communities, driven by D4Z4 repeat sharing.

**Pattern 2 — Substantial cross-arm affinity**:
- **chr11_p** (62.7% cross-arm): 281/448 sequences cluster with communities dominated by chr9_q (C40: 107 seqs), chr17_p (C13: 73 seqs), chr7_p (C26: 59 seqs), or chr5_q (C25: 42 seqs). The most polymorphic autosomal arm.
- **chr16_q** (19.8%): 89/449 cross-arm toward chr9_q (C40: 61 seqs) and chr3_q (C23: 26 seqs).
- **chr6_q** (8.9%): 39/440 cross-arm toward chr2_q.
- **chr15_q** (8.5%): 38/449 cross-arm toward chr6_q.

**Pattern 3 — Nearly monomorphic**: Most arms (>95% self):
- **8 arms at 0% cross-arm**: chr10_p, chr12_p, chr17_q, chr18_p, chr1_q, chr20_q, chr21_q, chr9_q.
- **10 arms at 0.2–0.7% cross-arm** (1–3 stray sequences each): chr13_q (0.7%), chr2_q (0.7%), chr14_p (0.4%), chrX_p (0.3%), chrX_q (0.3%), chr4_q (0.3%), chr1_p (0.2%), chr4_p (0.2%), chr12_q (0.2%), chr19_q (0.2%).

Cross-arm affinity rates across all 41 arms:

| Arm | Total | Self | Cross-arm | Cross-arm % | Top foreign arm |
|-----|-------|------|-----------|-------------|-----------------|
| chr10_q | 357 | 0 | 357 | 100.0% | chr4_q (354) |
| chr13_p | 76 | 0 | 76 | 100.0% | chr14_p (76) |
| chr15_p | 119 | 0 | 119 | 100.0% | chr14_p (119) |
| chr21_p | 116 | 0 | 116 | 100.0% | chr14_p (115) |
| chr22_p | 223 | 0 | 223 | 100.0% | chr14_p (222) |
| chrY_p | 92 | 0 | 92 | 100.0% | chrX_p (90) |
| chrY_q | 101 | 0 | 101 | 100.0% | chrX_q (101) |
| chr11_p | 448 | 167 | 281 | 62.7% | chr9_q (107) |
| chr16_q | 449 | 360 | 89 | 19.8% | chr9_q (61) |
| chr6_q | 440 | 401 | 39 | 8.9% | chr2_q (21) |
| chr15_q | 449 | 411 | 38 | 8.5% | chr6_q (38) |
| chr7_q | 446 | 418 | 28 | 6.3% | chr13_q (20) |
| chr16_p | 442 | 416 | 26 | 5.9% | chr9_q (25) |
| chr22_q | 450 | 429 | 21 | 4.7% | chr19_q (21) |
| chr9_p | 448 | 429 | 19 | 4.2% | chr20_p (19) |
| chr7_p | 438 | 420 | 18 | 4.1% | chr9_q (17) |
| chr6_p | 417 | 403 | 14 | 3.4% | chr15_q (13) |
| chr19_p | 434 | 420 | 14 | 3.2% | chr20_p (10) |
| chr17_p | 448 | 435 | 13 | 2.9% | chr7_p (11) |
| chr20_p | 453 | 441 | 12 | 2.6% | chr9_p (5) |
| chr3_q | 436 | 425 | 11 | 2.5% | chr16_q (5) |
| chr8_p | 458 | 447 | 11 | 2.4% | chr19_p (10) |
| chr5_q | 447 | 439 | 8 | 1.8% | chr19_p (3) |
| chr13_q | 449 | 446 | 3 | 0.7% | chr18_p (3) |
| chr2_q | 450 | 447 | 3 | 0.7% | chr9_p (1) |
| chr14_p | 229 | 228 | 1 | 0.4% | chr4_q (1) |
| chrX_p | 327 | 326 | 1 | 0.3% | chr9_p (1) |
| chrX_q | 330 | 329 | 1 | 0.3% | chr14_p (1) |
| chr4_q | 357 | 356 | 1 | 0.3% | chr19_q (1) |
| chr1_p | 438 | 437 | 1 | 0.2% | chr5_q (1) |
| chr4_p | 448 | 447 | 1 | 0.2% | chr6_p (1) |
| chr12_q | 449 | 448 | 1 | 0.2% | chr7_q (1) |
| chr19_q | 457 | 456 | 1 | 0.2% | chr1_q (1) |
| chr10_p | 443 | 443 | 0 | 0.0% | — |
| chr12_p | 454 | 454 | 0 | 0.0% | — |
| chr17_q | 451 | 451 | 0 | 0.0% | — |
| chr18_p | 446 | 446 | 0 | 0.0% | — |
| chr1_q | 445 | 445 | 0 | 0.0% | — |
| chr20_q | 409 | 409 | 0 | 0.0% | — |
| chr21_q | 449 | 449 | 0 | 0.0% | — |
| chr9_q | 450 | 450 | 0 | 0.0% | — |

**Conclusion.** 1,740 of 15,668 sequences (11.1%) show cross-arm affinity at the 50-community partition. The key biological signals are: (1) sex chromosome PAR sharing (chrY absorbed by chrX at both PAR1 and PAR2); (2) acrocentric p-arm homogenization (chr13_p, chr15_p, chr21_p, chr22_p all absorbed by chr14_p-dominated communities); (3) D4Z4-driven chr10_q/chr4_q co-clustering (chr10_q fully absorbed into chr4_q communities); (4) chr11_p as the most polymorphic autosomal arm (62.7% cross-arm, distributed across chr9_q, chr17_p, chr7_p, and chr5_q communities). The 8 arms with 0% cross-arm (chr10_p, chr12_p, chr17_q, chr18_p, chr1_q, chr20_q, chr21_q, chr9_q) represent truly arm-specific subtelomeric content not shared with any other chromosome end.

Having defined community structure at both arm and sequence levels, the next sections characterize the gene and repeat content of each community.

### Files and scripts

**Scripts:**
| Script | Description |
|--------|-------------|
| `/moosefs/guarracino/HPRCv2/scripts/classify_contigs.py` | Contig classification by p/q-arm alignment and telomere filtering |
| `/moosefs/guarracino/HPRCv2/scripts/trim-telomeres.sh` | Telomere trimming from flank sequences |
| `/moosefs/guarracino/HPRCv2/scripts/phr_wfmash_array.sh` | SLURM array job for wfmash all-vs-one alignments |
| `/moosefs/guarracino/HPRCv2/scripts/phr_post_wfmash.sh` | Post-wfmash processing (impg indexing, PHR detection) |
| `/moosefs/guarracino/HPRCv2/scripts/find-multichr-regions-incremental.py` | Sliding-window inter-chromosomal region detection with early stopping |
| `/moosefs/guarracino/HPRCv2/scripts/extract_flanking_sequences.py` | Extract 100kb flanking regions centromere-ward of PHR boundaries |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R` | Distance matrix, MDS, UMAP, community detection |
| `/moosefs/guarracino/HPRCv2/scripts/community/detect_communities.R` | Leiden community detection (arm-level + sequence-level) |
| `/moosefs/guarracino/HPRCv2/scripts/community/extract-seq-assignments.R` | Sequence-level community assignments with metadata |
| `/moosefs/guarracino/HPRCv2/scripts/community/community-utils.R` | Shared utility functions |
| `/moosefs/guarracino/HPRCv2/scripts/community/compare-community-levels.R` | Arm-level vs sequence-level comparison (ARI, NMI) |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-seq-community-structure.R` | Sequence-level community visualization |

**Input data:**
| File | Description |
|------|-------------|
| `/moosefs/pangenomes/HPRCv2/*.fa.gz` | 465 HPRCv2 assemblies |
| `/moosefs/guarracino/HPRCv2/PHR_III/pq-classification/chm13.centromeres.approximate.bed` | CHM13 centromere coordinates |

**Intermediate files:**
| File | Description |
|------|-------------|
| `/moosefs/guarracino/HPRCv2/PHR_III/pq-classification/contig_classifications.tsv` | 12,649 classified contigs |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz` | 18,827 subtelomeric flanks |
| `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95/*.paf.gz` | 18,827 all-vs-all PAF alignments (88 GB) |
| `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv` | Inter-chromosomal region calls |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz` | 15,668 PHR sequences |
| `/moosefs/guarracino/HPRCv2/PHR_III/pggb/.../similarity.tsv.gz` | Jaccard pairwise similarity (10.8 GB compressed) |

**Output files:**
| File | Description |
|------|-------------|
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.dist_matrix.rds` | 15,668 × 15,668 distance matrix |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv` | 15 arm-level communities |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv` | 50 sequence-level communities |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv` | 41 × 41 arm-level distance matrix |


