# Subtelomeric Analysis Report

*Last updated: 2026-03-30*

## Overview

**What it does.** Identifies and characterizes inter-chromosomal subtelomeric sequence sharing across 232 HPRCv2 individuals, validates the findings with 3D genome data from multiple technologies and organisms.

**How.** Extract the terminal 500 kb from each chromosome arm across 465 near-complete assemblies (contigs ≥ 1 Mb). Align all-vs-all at ≥95% identity to find which arms share sequence across different chromosomes. Build a pangenome graph from these shared regions, then cluster arms into communities by graph similarity. Annotate communities with genes and repeats. Validate by testing whether same-community arms sit closer in 3D nuclear space — using Hi-C (6 human samples), Pore-C and CiFi (HG002), Dip-C single-cell (16 GM12878 cells), sperm single-cell (20 cells), RPE-1 (3 datasets), and mouse meiotic Hi-C (4 stages).

**Key metrics.** 232 individuals (465 near-complete assemblies), 48 chromosome arms (41 with inter-chromosomal signal), 15,668 PHR sequences, 15 arm-level communities, 50 sequence-level communities. 3D validation across 4 technologies, 2 cell types (lymphoblastoid cell lines [LCL], retinal pigment epithelium [RPE-1]), sperm, and mouse meiosis.

**Result.** Chromosome arms that share subtelomeric sequence cluster into discrete communities. These communities are reflected in 3D nuclear proximity — arms in the same community are physically closer than arms in different communities. This holds in bulk Hi-C, single-cell Dip-C, haploid sperm, and mouse meiotic cells.

**Conclusion.** Subtelomeric regions form a structured system of inter-chromosomal sharing, shaped by recurrent ectopic exchange and reflected in 3D nuclear organization across cell types and species.

---

## 1. Contig classification and telomere filtering

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

## 2. Subtelomeric flank extraction

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

## 3. All-vs-all subtelomeric alignment

**What it does.** Aligns every subtelomeric flank against every other flank to detect shared sequence at ≥95% identity.

**How.** wfmash v0.23.0-41-gb5f0ff1c (`-p 95 -t 48 --quiet`): each of the 18,827 sequences serves as target in turn, with all 18,827 sequences as queries. The 95% identity threshold was chosen to capture the high-identity peak of Ambrosini et al.'s (2007) bimodal duplicon identity distribution (peaks at 91% and 98%), accepting that older exchanges (≤91% identity) are missed (see §16, limitation 1).

**Key metrics.** Identity threshold: 95%. Output: 18,827 PAF.gz files.

**Result.** 18,827 PAF.gz files in `all-vs-all.1Mb.p95/`, indexed with impg for fast downstream queries.

**Conclusion.** These alignments provide the raw inter-sequence similarity data for inter-chromosomal region detection (§4) and pangenome graph construction (§5).

---

## 4. Inter-chromosomal region detection

**What it does.** Identifies, for each subtelomeric flank, the region of inter-chromosomal similarity, that is the portion of the flank that shares sequence with other chromosomes at ≥95% identity.

**How.** `find-multichr-regions-incremental.py` slides windows from the telomere inward along each flank. For each window, it queries impg for alignments and counts how many different chromosomes have matches. Scanning stops early when consecutive windows fail to meet the threshold.

**Key metrics.** Window = 5 kb, step = 5 kb, max distance from telomere = 500 kb, min alignment identity = 0.95, min different chromosomes = 2, min alignment count per chromosome = 5, min consecutive failing windows = 4, min output region length = 3 kb.

**Result.** 18,826 data rows in `all-vs-all.1Mb.p95.id95.len.tsv` (one per flank, excluding the chimeric chr18_q contig — see below). 15,668 sequences (83.2%) have inter-chromosomal matches; 3,158 (16.8%) have none. The 15,668 PHR sequences span 41 of 48 chromosome arms. Similarity region lengths: median 105 kb, mean 144 kb, range 5–500 kb.

**Conclusion.** The majority of subtelomeric flanks share detectable sequence with other chromosomes. Seven arms have no inter-chromosomal signal at these thresholds.

### 4.1 Chimeric contig exclusion

**What it does.** Identifies and removes a chimeric contig whose inter-chromosomal signal was artifactual.

**How.** The chr18_q flank from NA18982#1 (JBKABS010000018.1, 84.4 Mb) was flagged because its inter-chromosomal signal (chrX, chrY) arose from a scaffolding artifact: the contig fuses chr18 sequence with 966 kb of chrX PAR1 across a 100 bp NNN scaffold join.

**Key metrics.** Junction at query ~83.37–83.38 Mb confirmed by both wfmash and minimap2 v2.30 (mapq 60). Flagger labels the junction as NNN and the chrX portion as Hap (haploid). A 2,826 bp terminal telomeric tract (~471 TTAGGG repeats) precedes the NNN gap. NA18982#1 has no separate chrX contig.

**Result.** One sequence removed; 15,668 PHR sequences retained.

**Conclusion.** The NNN gap indicates a scaffolding artifact, not a contiguous assembly error.

### 4.2 Zero-signal arms and one-copy regions

**What it does.** Compares the 7 arms with no inter-chromosomal signal against Ambrosini et al.'s (2007) 6 "one-copy DNA (TTAGGG)n-adjacent regions", that is chromosome ends where the sequence immediately adjacent to the terminal telomeric tract "is unique in the genome."

**How.** The 7 zero-signal arms (chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q) were compared to Ambrosini's 6 one-copy regions (chr7_q, chr8_q, chr11_q, chr12_q, chr18_q, chrX_p/chrY_p). For arms with signal (chr7_q, chr12_q, chrX_p/chrY_p), per-arm chromosome matching was computed from the `all-vs-all.1Mb.p95.id95.len.tsv` file.

**Key metrics.**
- 3 one-copy regions (chr8_q, chr11_q, chr18_q) are among the 7 zero-signal arms.
- chrX_p/chrY_p: 416/419 sequences with signal (99.3%) match exclusively chrX/chrY and no other chromosome.
- chr12_q: 449/449 (100%) match exclusively chr7/chr12.
- chr7_q: 424/446 (95.1%) match exclusively chr7/chr12; the 22 outliers match many chromosomes.
- chr7_q shared regions: median 40 kb (range 5–45 kb). chr12_q: median 25 kb (range 15–45 kb).

**Result.** Four of 6 one-copy regions are confirmed: chr8_q, chr11_q, and chr18_q have zero signal; chrX_p/chrY_p shares sequence only within the PAR1 pair. The remaining 2 one-copy regions (chr7_q, chr12_q) form a previously undetected private pair, analogous to chrX_p/chrY_p — Ambrosini listed them as independent one-copy regions, but the population-scale data reveals they share small regions (5–25 kb) exclusively with each other.

**Conclusion.** The population-scale analysis confirms Ambrosini et al.'s one-copy classification for 4 of 6 regions and refines the classification for chr7_q/chr12_q by revealing a private pair invisible in the single reference genome. The chr7_q/chr12_q pair forms community C4 (§6).

---

## 5. Pangenome graph construction and pairwise similarity

**What it does.** Builds a pangenome graph of all PHR sequences and computes pairwise Jaccard similarity between every pair of sequences. The Jaccard similarity measures the fraction of graph nodes (sequence segments) shared between any two sequences.

**How.** `pggb -p 95 -D /scratch` constructs the pangenome graph from the 15,668 PHR sequences. `odgi similarity --all -P` then computes the Jaccard similarity for every pair of sequences. Values range from 0 (no shared nodes) to 1 (identical graph traversals).

**Key metrics.** 15,668 input sequences. Graph construction identity threshold: 95%. Output: 15,668 × 15,668 Jaccard similarity matrix (12 GB compressed, 245 million pairwise entries including both directions).

**Result.** Pangenome graph and full pairwise similarity matrix in `*.similarity.tsv.gz`.

**Conclusion.** The Jaccard similarity matrix provides the distance measure for all downstream community detection (§6) and heterogeneity analyses (§10). Unlike alignment identity, Jaccard similarity captures structural variation (insertions, deletions, inversions) through graph topology.

---

## 6. Community detection

### 6.1 Arm-level community detection

**What it does.** Groups the 41 chromosome arms with detected PHRs into communities based on their average pairwise subtelomeric sequence similarity. Seven of the 48 human chromosome arms (chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q) are excluded because no inter-chromosomal PHR was detected at these arms (§4).

**How.** The pairwise Jaccard similarity (§5) is converted to distance (1 − Jaccard) and averaged across all sequence pairs per arm pair to produce a 41 × 41 arm-level distance matrix. Two clustering methods are applied:

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

### 6.1.1 f7501 (L78442) direct alignment validation

**What it does.** Tests whether direct alignment of the f7501 cosmid reproduces Mefford & Trask's (2002) Fig. 3 per-arm distribution at population scale.

**How.** The L78442.1 cosmid (36.3 kb, the original f7501 clone from Mefford & Trask 2002, with 3 regions of similarity to olfactory receptor genes including the expressed *OR-A* gene) was aligned against each of the 18,827 subtelomeric flanks individually (minimap2 `-x asm20`, one-vs-one to avoid multi-mapping). Haplotypes with ≥30 kb matching bases were counted as carrying f7501.

**Key metrics.** Per-arm distribution across 465 near-complete assemblies (two haplotypes for each of the 232 individuals plus CHM13). Population enrichment tested by Fisher's exact test (one-sided greater) for each superpopulation independently over these 465 assemblies, with the single haploid CHM13 grouped into EUR (AFR=134, AMR=88, EAS=104, EUR=65, SAS=74 haplotypes; these sum to 465). The cohort demographics proper are 232 individuals / 464 HPRC haplotypes (EUR = 32 individuals, 64 haplotypes); CHM13 is a reference anchor, not a study individual. The table reports the most significantly enriched superpopulation per arm. The per-arm distribution reproduces Mefford & Trask's Figure 3 at population scale:

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
- **chr15_q refinement**: Mefford classified chr15_q as "FIXED" based on 52 individuals by FISH. Across the 465 assemblies, the overall prevalence (85.6%) confirms near-fixation, but with significant EUR enrichment (p=2.5e-04): 98.5% of EUR haplotypes carry f7501 (64/65) vs only 64.9% of AFR (87/134). Non-AFR populations are uniformly high: AMR 92.0%, SAS 91.9%, EAS 94.2%, EUR 98.5%. This is consistent with ongoing f7501 loss in African populations or incomplete lineage sorting at this locus.
- **Novel locations not in Mefford**: chr1_p (5 haplotypes, C11), chr20_p (2, C12), and chrX_q (1, C14) carry f7501 at ≥30 kb but are absent from Mefford's Figure 3 (which used FISH on 52 individuals). These may represent rare or recently acquired f7501 copies below the detection threshold of the original FISH survey.
- **chr4_q and chr19_q absence**: Mefford's Figure 3 lists chr4_q and chr19_q as variable f7501 sites detected by FISH. The L78442 cosmid produces zero alignments ≥30 kb to any chr19_q flank (457 tested) and only ~954 bp partial matches to chr4_q. The f7501 copies at these arms are too divergent for sequence alignment at asm20 stringency, despite being detectable by FISH cross-hybridization.

**Conclusion.** The f7501 distribution across 465 near-complete assemblies reproduces the FISH-based patterns observed by Mefford & Trask (2002) in 52 individuals and extends them with three newly identified AFR-enriched sites, population-specific enrichments at chr2_q (SAS) and chr6_p (AMR), and three novel f7501 locations (chr1_p, chr20_p, chrX_q) not reported by FISH.

---

### 6.2 Sequence-level community detection

**What it does.** Clusters individual subtelomeric sequences (not arm averages) into communities, resolving within-arm polymorphism invisible at the arm level.

**How.** Leiden community detection on the full 15,668 × 15,668 pairwise Jaccard distance matrix (§5). A k-NN graph is constructed with exponential decay edge weights. A joint scan over k (10, 25, 50, 75, 100, 125) and resolution (0.1 to 3.0 in steps of 0.1) selects the (k, resolution) pair that maximizes modularity within a target community count range (5–50). Optimal: k=75, resolution=0.8. Run by `detect_communities.R`.

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

### 6.3 Arm-level vs sequence-level community nesting

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

### 6.4 Polymorphic arm subgroup characterization

**What it does.** Classifies each sequence as self-arm or cross-arm and characterizes three patterns of subtelomeric polymorphism across the 41 chromosome arms (34 polymorphic, 7 monomorphic).

**How.** A sequence has **cross-arm affinity** when its subtelomeric sequence content is more similar to sequences from a different chromosome arm than to its own arm — measured by Jaccard distance on shared pangenome graph nodes (§5). Operationally, the sequence clusters into a community dominated by a foreign arm rather than its own. Such cross-arm similarity is attributed to inter-chromosomal exchange (Linardopoulou et al. 2005; Mefford & Trask 2002), though the present data measures similarity, not exchange events directly.

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

## 8. Annotation preprocessing

**What it does.** Intersects gene annotations (Liftoff GFF3) and TAR1 repeat entries (RepeatMasker) with the subtelomeric similarity regions to characterize gene and repeat content per community.

**How.** 464 haplotype-specific GFF3 files + CHM13 reference processed (465 total). The remote HPRC annotation index covers 462 haplotypes; HG002 (2 haplotypes) was annotated separately using JHU Liftoff v0.6 on the HG002v1.1 assembly. TAR1 entries extracted from RepeatMasker BED files, with CHM13 converted to PanSN naming and HG002 converted from bigBed to 10-column BED.

**Key metrics.** 18,827 total sequences → 15,668 retained (3,159 with no inter-chromosomal signal removed; chimeric chr18_q already excluded in §4).

**Result.**
- **Gene annotations**: 173,881 gene annotations extracted (374 unique genes) across 232 individuals and 39 arms. Two arms — chr7_q and chr12_q — have zero gene annotations because their similarity regions are confined to 5–25 kb at the telomeric tip, below the detection threshold of Liftoff annotation pipelines.
- **TAR1 repeats**: 21,544 TAR1 entries across 14,816 sequences (94.6%) and all 41 arms.

**Conclusion.** The gene and repeat annotations provide the basis for characterizing the biological content of each community (§9) and testing whether exchange status affects gene repertoire (§10).

### 8.1 TAR1 prevalence

**What it does.** Quantifies TAR1 (Telomere-Associated Repeat 1; Brown et al. 1990) prevalence across chromosome arms.

**How.** RepeatMasker TAR1 entries intersected with PHR regions; prevalence and density computed per arm.

**Key metrics.** Present in all 41 arms. Most autosomal arms: >99%. Near-absent: chrX_p (1/327, 0.3%), chrY_p (1/92, 1.1%). Acrocentric p-arms intermediate: chr15_p 73.1%, chr14_p 76.0%, chr13_p 77.6%, chr21_p 78.4%, chr22_p 78.9%. Highest TAR1 density: chr18_p (mean 4.00 copies/sequence), chr9_q (mean 2.82 copies/sequence).

**Result.** TAR1 is near-universal outside PAR (94.6% of all sequences). The near-absence from PAR1 (chrX_p/chrY_p) is consistent with pseudoautosomal regions using obligate meiotic crossover (Rouyer et al. 1986) rather than repeat-mediated exchange.

**Conclusion.** TAR1's near-universal presence makes it difficult to assess whether it plays a functional role in inter-chromosomal exchange or is a co-located passenger (§17.3, prediction 3).

### 8.2 TAR1 positional distribution

**What it does.** Maps where TAR1 sits within the PHR region, measuring distance from the telomeric end.

**Key metrics.** 66.9% of 21,544 TAR1 entries within 10 kb of the telomere; 70.3% within 25 kb (median 0.3 kb, mean 44.1 kb). TAR1 is overwhelmingly telomere-proximal. Most arms (both p and q) have median TAR1 distance <1 kb from the telomere. Exceptions with deeper TAR1: acrocentric p-arms (chr13_p–chr22_p: median 179–196 kb, reflecting large PHR regions where TAR1 sits near the telomeric end of a long duplicated zone), chr9_q (median 223.3 kb), chr8_p (140.4 kb), chr18_p (87.3 kb), chr19_p (68.6 kb), chr6_p (56.7 kb).

**Result.** TAR1 is a telomere-proximal element. In arms with short PHR regions (<25 kb), TAR1 sits at the very tip. In arms with long PHR regions (>100 kb), TAR1 remains near the telomere while the duplicated zone extends centromere-ward beyond it.

**Conclusion.** TAR1 positional distribution is consistent with it being a telomere-associated satellite that marks the boundary between terminal telomeric sequence and the subtelomeric duplicated zone.

### 8.3 TAR1 per community

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

### 8.4 Internal (TTAGGG)n islands

**What it does.** Identifies and quantifies short degenerate telomeric repeat tracts inside subtelomeric regions, distinct from the terminal telomeric array and from TAR1 satellite blocks.

**How.** Canonical (TTAGGG) and variant (TGAGGG, TCAGGG, TTGGGG) telomeric motifs searched across all 15,668 PHR sequences using `seqkit locate`; overlapping hits merged (12 bp gap tolerance) and filtered to 50–1000 bp. Ambrosini et al. (2007) identified these elements (150–823 bp, mostly 150–200 bp) at "duplicon boundaries inside subtelomeric regions."

**Key metrics.** 18,352 islands across 8,321 sequences (53.1%) and all 41 arms. Island lengths: median 79 bp, mean 102 bp. Most sequences carry 1–4 islands; maximum 22 per sequence.

**Result.** Highest island counts per arm: chr20_q (1,765 islands), chr12_q (1,149), chr16_p (898), chr18_p (851). The per-arm positional distribution shows telomere-proximal bias for most p-arms (16/21 with median < 500 bp from telomere). Five p-arms with larger PHR regions have deeper islands (chr11_p 153 kb, chr6_p 99 kb, chr20_p 65 kb, chr18_p 54 kb, chr16_p 1.1 kb). This is consistent with Ambrosini et al.'s (2007) observation that internal (TTAGGG)n-like sequences "almost always co-localize to duplicon boundaries" — islands mark past duplication breakpoints, so arms with more layered duplication history have islands distributed deeper into the subtelomeric zone.

**Conclusion.** Island lengths (median 79 bp) are shorter than Ambrosini et al.'s (2007) observation that "most [islands are] in the 150-200 bp range" (max 823 bp). The difference is methodological: our detection uses degenerate motif search (canonical + 3 variant hexamers + mixed patterns) with a 50 bp minimum across 15,668 sequences from 465 near-complete assemblies, capturing shorter and more degenerate tracts that were not detected in Ambrosini's single-reference analysis.

### 8.5 TTAGGG island boundary enrichment test

**What it does.** Tests whether internal (TTAGGG)n islands specifically mark the PHR boundary (the transition between duplicated and unique sequence).

**How.** Island positions compared against the centromere-ward boundary of each PHR region. KS test for non-uniformity; binomial test for boundary-proximal enrichment.

**Key metrics.** KS test stat=0.37, p < 1e-300 (non-uniform). Mean fractional distance = 0.54 (islands sit on average 54% of the way through the PHR region). 42.5% of islands within 5 kb of boundary (7,806/18,352).

**Result.** Islands are distributed throughout the PHR region (mean fractional position 0.54), with a slight enrichment near the boundary (42.5% within 5 kb). The KS test rejects uniformity (p < 1e-300) but the distribution is not strongly polarized toward either end.

**Conclusion.** This test evaluated the PHR outer boundary (the single transition from duplicated to unique sequence). Ambrosini et al.'s (2007) claim concerns co-localization with internal duplicon-to-duplicon boundaries throughout the subtelomeric zone — a different feature not tested here. A direct test would require mapping individual duplicon block boundaries within each PHR region.

### 8.6 TTAGGG island length distribution

**What it does.** Characterizes the size distribution of internal (TTAGGG)n islands at population scale.

**Key metrics.** Mode at 50–74 bp (8,433 islands, 46.0%); monotonically decreasing at longer lengths: 75–99 bp (4,109, 22.4%), 100–149 bp (3,047, 16.6%), 150–199 bp (1,455, 7.9%), 200–299 bp (838, 4.6%), 300–499 bp (374, 2.0%), 500–1000 bp (96, 0.5%).

**Result.** Ambrosini et al. (2007) reported a mode at 150–200 bp; this subrange contains only 7.9% of population-scale islands, with no secondary mode.

**Conclusion.** The 50 bp minimum filter captures shorter degenerate tracts below Ambrosini et al.'s detection threshold, accounting for the shorter median (79 bp vs 150–200 bp) — a methodological, not biological, difference.

### 8.7 TTAGGG island motif composition

**What it does.** Quantifies the canonical vs variant telomeric hexamer content of internal islands.

**Key metrics.** From 18,352-island dataset (296,406 total hexamer instances). Canonical TTAGGG: 52.2% of hexamer instances (154,886/296,406). Variants: TGAGGG 19.0%, TTGGGG 16.0%, TCAGGG 12.7%. Only 32.2% of islands are "pure canonical" (≥80% TTAGGG+CCCTAA); 47.2% are variant-dominant (<50% canonical).

**Result.** The high variant content — consistent with Ambrosini et al.'s (2007) identification of the same three dominant variant motifs — indicates these islands are substantially degenerate relative to the terminal telomeric array. Linardopoulou et al. (2005) found degenerate telomeric repeats enriched at "4% of subtelomeric DSB sites" (vs 0.5% background), suggesting these islands may have been appended during DSB repair.

**Conclusion.** The degenerate composition supports interpretation of internal (TTAGGG)n islands as ancient remnants of telomeric sequence incorporated during subtelomeric duplication events.

### 8.8 TTAGGG island count by cross-arm status

**What it does.** Tests whether cross-arm exchange status affects the number of internal telomeric islands.

**Key metrics.** 8,321 sequences with islands: 1,569 (18.9%) cross-arm, 6,752 (81.1%) self-arm. Mean island count: cross-arm 2.08, self-arm 2.24. Mann-Whitney U test z = −1.89, p = 0.04.

**Result.** Marginally significant difference (p=0.04): self-arm sequences carry slightly more islands on average (2.24 vs 2.08), though the effect size is small.

**Conclusion.** Cross-arm exchange status does not affect internal telomeric island count.

### 8.9 Terminal telomere tract length by community

**What it does.** Tests whether terminal telomere length varies across communities.

**How.** Terminal telomere repeat tract lengths (from `.telo.tsv`) matched to arm-level community assignments.

**Key metrics.** Kruskal-Wallis H = 89.63, p = 4.5e-13. Medians range from 7,638 bp (C10, chr17_p) to 9,404 bp (C13, chr4_p). Correlation between terminal telomere tract length and TTAGGG island count: across all 15,666 sequences (including those with 0 islands), Spearman rho = −0.039 (p = 1.2e-6); restricted to the 8,320 sequences with at least one island, rho = −0.057 (p = 2.4e-7). Per-arm breakdown: p-arm all N = 7,396 rho = −0.048 (p = 3.6e-5), p-arm with islands N = 4,059 rho = −0.078 (p = 6.2e-7); q-arm all N = 8,270 rho = −0.031 (p = 4.3e-3), q-arm with islands N = 4,261 rho = −0.040 (p = 9.9e-3). Telomere lengths: min 470 bp, median 8,249 bp, mean 8,677 bp, max 33,826 bp.

**Result.** Telomere lengths vary significantly across communities. C13 (chr4_p) and C15 (PAR1) have the longest telomeres (mean 9,601 bp and 9,247 bp); C9 (chr16_p) and C10 (chr17_p) the shortest (mean ~8,180 bp). The correlation between telomere length and TTAGGG island count is weakly negative across all arms and subsets (rho from −0.031 to −0.078) — sequences with longer telomeres tend to have slightly fewer internal (TTAGGG)n islands, but the effect size is minimal. The correlation is strongest for p-arm sequences restricted to those with islands (rho = −0.078), consistent with the cleaner telomere-length extraction from p-arm coordinates.

**Conclusion.** Community membership is associated with terminal telomere length variation, but the biological significance of this association is unclear.

## 9. Community gene enrichment analysis

**What it does.** Tests which genes are shared across community member arms, which are community-specific, and which span multiple communities.

**How.** Gene annotations (§8) grouped by arm-level community. Shared gene instances counted across arms within each community. Fisher's exact tests for per-community enrichment (116 tests, BH-corrected).

**Key metrics.** 374 unique genes across 39 arms (chr7_q and chr12_q have zero gene annotations). 576 shared gene instances across arms within communities. 93 genes community-specific (found in only one community). 216 genes in 2 or more communities.

### 9.1 Biotype composition

**Result.** Subtelomeric genes are predominantly pseudogenes and ncRNA, with low protein-coding content:

| Community | Arms | Total genes | Protein-coding (%) | Pseudogene (%) | ncRNA (%) |
|-----------|------|-------------|---------------------|-----------------|-----------|
| C1 | chr4_q, chr10_q | 59 | 5.1 | 86.4 | 8.5 |
| C2 | chr10_p, chr18_p | 14 | 21.4 | 28.6 | 50.0 |
| C3 | chr3_q, chr7_p, chr9_q, chr11_p, chr16_q, chr19_p | 195 | 6.7 | 55.4 | 37.9 |
| C4 | chr7_q, chr12_q | 0 | — | — | — |
| C5 | chr6_p, chr9_p, chr12_p, chr20_q | 101 | 5.9 | 61.4 | 32.7 |
| C6 | chr1_q, chr13_q, chr17_q, chr19_q, chr21_q, chr22_q | 34 | 8.8 | 55.9 | 35.3 |
| C7 | chr13_p, chr14_p, chr15_p, chr21_p, chr22_p | 118 | 5.1 | 65.3 | 29.7 |
| C8 | chr15_q | 72 | 4.2 | 50.0 | 45.8 |
| C9 | chr16_p | 136 | 5.1 | 56.6 | 38.2 |
| C10 | chr17_p | 39 | 5.1 | 64.1 | 30.8 |
| C11 | chr1_p, chr5_q, chr6_q, chr8_p | 162 | 4.3 | 57.4 | 38.3 |
| C12 | chr2_q, chr20_p | 129 | 3.9 | 57.4 | 38.8 |
| C13 | chr4_p | 12 | 8.3 | 66.7 | 25.0 |
| C14 | chrX_q, chrY_q | 60 | 8.3 | 58.3 | 33.3 |
| C15 | chrX_p, chrY_p | 28 | 32.1 | 39.3 | 28.6 |

C15 (PAR1) has 32.1% protein-coding genes, reflecting the functional gene content of the pseudoautosomal region. C4 (chr7_q/chr12_q) has zero gene annotations — its similarity regions are confined to 5–25 kb at the telomeric tip. C2 (chr10_p/chr18_p) has 21.4% protein-coding and 50.0% ncRNA.

### 9.2 Olfactory receptor genes

**What it does.** Tests whether olfactory receptor gene families are present across communities, as predicted by Mefford & Trask (2002) and Ambrosini et al. (2007).

**Key metrics.** 10 OR4F family genes detected across 7 communities (C3, C5, C8, C9, C11, C12, C14). OR4F5 and OR4F8P most widespread (14 arms each). IQSEC3 detected in C5 (chr12_p, 453 samples).

**Result.** Confirms at population scale that "human subtelomeres can contain genes, such as members of the olfactory receptor gene family" (Mefford & Trask 2002) and extends Ambrosini et al.'s (2007) OR duplicon architecture (Table 1, Block 2) to 465 near-complete assemblies.

### 9.3 Ambrosini subtelomere-specific blocks → Leiden communities

**What it does.** Maps each of Ambrosini et al.'s (2007) 11 subtelomere-specific duplicon block entries (Table 1, numbered 1–3, 5–8, 10–12; block 6' is a variant of block 6) to the present Leiden communities.

**Result.** Their anchor telomeres map systematically to the present Leiden communities:

| Ambrosini block | Anchor | Size | Copies | Community | Diagnostic gene(s) | Confirmed |
|-----------------|--------|------|--------|-----------|---------------------|-----------|
| 1 | chr1_p | 25 kb | 4 | C11 (chr1_p,chr5_q,chr6_q,chr8_p) | — | by arm membership |
| 2 | chr15_q | 88 kb | 7 | C8 (chr15_q singleton) | OR4F17 (1 arm, 2 samples), OR4F4 (1 arm, 409 samples) | gene-level |
| 3 | chr1_p | 35 kb | 1 | C11 | — | by arm membership |
| 5 | chr2_p | 17 kb | 5 | No PHR signal | RPL23AP7 (chr2_p excluded; no inter-chromosomal signal at ≥95%) | — |
| 6 | chr3_q | 38 kb | 5 | C3 | RPL23AP7-related | by arm membership |
| 6' | chr11_p | 11 kb | 1 | C3 | RYD5 (not in current annotations) | by arm membership |
| 7 (D4Z4) | chr4_q | 28 kb | 1* | C1 (chr4_q,chr10_q) | DUX4L1–DUX4L44 (28 genes, 2 arms) | gene-level |
| 8 (TUBB4q) | chr4_q | 14 kb | 6 | C1 | TUBB4q (not in current annotations) | by arm membership |
| 10 | chr2_q | 49 kb | 1 | C12 (chr2_q,chr20_p) | FBXO25 (not in current annotations) | by arm membership |
| 11 (IL9R) | chr9_q | 36 kb | 6 | C3 (chr3_q,chr7_p,chr9_q,chr11_p,chr16_q,chr19_p) | IL9RP1 (1 arm) | gene-level |
| 12 | chr12_p | 15 kb | 1 | C5 (chr6_p,chr9_p,chr12_p,chr20_q) | IQSEC3 (1 arm, 453 samples) | gene-level |

IL9R pseudogenes also appear in C2 (chr10_p/chr18_p: IL9RP2, IL9RP4), C6 (IL9RP4), C9 (chr16_p: IL9RP3, IL9R), and C14/C15 (PAR: IL9R, SPRY3), consistent with Mefford & Trask (2002) who reported IL9R pseudogenes at chr9_q, chr10_p, chr16_p, and chr18_p.

**Conclusion.** The concordance between single-genome duplicon blocks (Ambrosini 2007) and population-scale Leiden communities confirms that the community structure captures the known subtelomeric duplicon architecture.

### 9.4 Ambrosini subterminal families → Leiden communities

**What it does.** Maps each of Ambrosini et al.'s (2007) 6 subterminal duplication families (Table 2, A–F) to the present Leiden communities.

**Result.** Their anchor telomeres map to the present Leiden communities:

| Family | Anchor | Size | Subterminal copies | Community | Key genes (Ambrosini) | Gene status |
|--------|--------|------|-------------------|-----------|----------------------|-------------|
| A | chr2_p | 7 kb | 6 | No PHR signal | RPL23AP7, FAM41C | chr2_p excluded (no inter-chromosomal signal at ≥95% identity); FAM41C detected in C11 (1 arm) |
| B | chr4_p | 17 kb | 10 | C13 (chr4_p) | RPL23AP7, FAM41C | C13 has DDX11L16 (1 arm); FAM41C not detected in C13 |
| C | chr9_p | 10 kb | 6 | C5 (chr6_p,chr9_p,chr12_p,chr20_q) | DDX11L-like, CXYorf1 | DDX11L family confirmed in C5 (7 members, 3–4 arms) |
| D | chr10_q | 22 kb | 10 | C1 (chr4_q,chr10_q) | RPL23AP7, FAM41C | RPL23A pseudogenes present in C1 (11 members, 1–2 arms); FAM41C not detected in C1 |
| E | chr17_p | 21 kb | 5 | C10 (chr17_p) | — | RPL23AP21/45/47/88 present (1 arm each) |
| F | chr18_p | 15 kb | 1 | C2 (chr10_p,chr18_p) | — | IL9RP2, IL9RP4 present (1 arm each) |

Family A (chr2_p) has no PHR signal — chr2_p is one of 6 arms excluded in §4 for lacking inter-chromosomal signal at ≥95% identity. In Ambrosini et al.'s Table 2, Family A has 6 subterminal + 12 subtelomeric + 1 non-subtelomeric = 19 total copies, so it is predominantly subtelomeric. The lack of PHR signal at ≥95% identity likely reflects sequence divergence below the 95% threshold rather than absence of duplicated content.

FAM41C (lncRNA) is detected only in C11 (chr1_p/chr5_q/chr6_q/chr8_p, 1 arm, 215 samples), not at the anchor communities of families A, B, or D — indicating that the FAM41C copies associated with these subterminal families either lie below the 95% identity threshold or reside outside the PHR regions.

The RPL23A pseudogene family is the most widespread subtelomeric duplicon marker (RPL23AP45 spans 10 communities, 21 arms), but the specific member RPL23AP7 cited by Ambrosini is not annotated in the current Liftoff gene models; the RPL23A pseudogene nomenclature has been revised since Ambrosini's analysis.

**Family C's DDX11L-like gene signature is confirmed:** DDX11L family members (7 members) are present across 3–4 arms of C5, consistent with this family's inter-chromosomal distribution. Family F (chr18_p) maps to C2, which has the highest TAR1 density (mean 2.51 copies/sequence; §8), consistent with chr18_p's known repeat-rich subtelomeric architecture.

**Conclusion.** Family A (chr2_p) has no PHR signal, consistent with sequence divergence below the 95% threshold. Family C's DDX11L-like gene signature is confirmed across 3–4 arms of C5.

### 9.5 Community-specific genes

**What it does.** Identifies genes found in exactly one community, marking distinct subtelomeric identities.

| Community | N specific | Key genes | Biological significance |
|-----------|-----------|-----------|------------------------|
| C7 (acrocentric p-arms) | 48 | MTCO1P34, MTCO3P26, MTCO3P33, MTCO3P34, SNX18P15, SOWAHCP1, ASNSP2 | Mitochondrial pseudogenes (MTCO) and rDNA-associated loci; MTCO enrichment in acrocentric p-arm subtelomeres is a novel observation enabled by T2T-quality assemblies of these previously unresolved regions |
| C1 (chr4_q/chr10_q) | 26 | DUX4L pseudogenes (22), AGGF1P1, CLUHP4, DBET, LOC100996375 | D4Z4 macrosatellite array: DUX4L pseudogenes are copies of the DUX4 gene embedded within each D4Z4 repeat unit; the terminal copy on a permissive 4qA haplotype can produce pathogenic DUX4 protein causing FSHD (Lemmers et al. 2010); CTCF binds within D4Z4 as an insulator (Ottaviani et al. 2009). Chr10_q copies lack the stabilizing polyadenylation signal |
| C15 (PAR1) | 10 | SHOX, PPP2R3B, PLCXD1, GTPBP6, P2RY8, LOC124905300, LOC102724521, LINC00685, FABP5P13, KRT18P53 | 5 protein-coding genes including SHOX (short stature homeobox), a key growth regulator whose haploinsufficiency causes Leri-Weill dyschondrosteosis |
| C11 | 4 | FAM87B, LINC00115, LINC01409, LOC124903817 | lncRNAs |
| C3 (f7501) | 3 | FAM41AY1, FAM41AY2, LOC105375112 | lncRNAs specific to the 6-arm f7501 cluster |
| C14 (PAR2) | 1 | LOC124905309 | Single pseudogene |
| C6 | 1 | LOC124907874 | Single protein-coding gene |

### 9.6 Hub genes

**What it does.** Identifies genes present in 3 or more communities — the common duplicon backbone spanning multiple community boundaries.

| Gene | N communities | N arms | Biotype |
|------|--------------|--------|---------|
| RPL23AP45 | 10 | 21 | pseudogene |
| SEPTIN14P22 | 9 | 22 | pseudogene |
| DDX11L16 | 9 | 20 | transcribed pseudogene |
| FAM138D | 9 | 17 | lncRNA |
| LOC101929828 | 9 | 21 | lncRNA |
| LOC102723681 | 9 | 23 | pseudogene |
| RPL23AP60 | 9 | 18 | pseudogene |
| RPL23AP87 | 9 | 19 | transcribed pseudogene |
| RPL23AP88 | 9 | 19 | pseudogene |

RPL23A pseudogenes (ribosomal protein L23a) and SEPTIN14 pseudogenes are the most widespread subtelomeric duplicon markers, consistent with their identification as core subtelomeric duplicon components by Ambrosini et al. (2007). The DDX11L (DEAD/H-box helicase), WASH (Wiskott-Aldrich syndrome protein homolog), and MIR6859 families also span 7–9 communities. These gene families are among the most widespread subtelomeric duplicon markers. The subtelomeric WASH copies are pseudogenes of WASHC1 (the catalytic subunit of the WASH complex).

**Conclusion.** The predominance of pseudogenes and ncRNA across all communities (28.6–86.4% pseudogene) is consistent with telomere position effect (TPE). Mefford & Trask (2002, citing Baur et al. 2001) noted that TPE operates in human cells, with "reporter genes near telomeres expressed at ten times lower levels," and suggested that subtelomeric regions might "buffer genes in chromosome-specific regions and in proximal subtelomeric domains from telomere-mediated repression." The exception — PAR1 (C15, 32.1% protein-coding) — involves a region with obligate crossover recombination. Fisher's exact tests for gene family enrichment per community (116 tests, BH-corrected) yield no significant results — the qualitative enrichments described above (MTCO in C7, DUX4L in C1, OR4F in C3) reflect presence patterns but do not survive multiple testing correction.

### 9.7 Sequence-level vs arm-level enrichment comparison

**What it does.** Tests whether gene enrichment patterns are preserved when using the finer 50-community partition.

**Result.** The sequence-level enrichment uses 50 communities that are subpartitions of the 15 arm-level communities. Gene annotations are assigned at the arm level (§9.1–§9.6), where each community groups 1–6 chromosome arms with shared gene content. At the 50-community level, most communities are pure single-arm (18/50) or near-pure (23/50), making gene enrichment largely redundant with the arm-level analysis. The arm-level enrichment (§9.1–§9.6) provides the biologically meaningful gene characterization; the sequence-level partition's value is in resolving within-arm polymorphism (§6.4, §10) rather than gene content.

**Conclusion.** The arm-level analysis captures the major duplicon block structure; the sequence-level analysis captures within-arm polymorphism.

### 9.8 Singleton/doubleton QC

**What it does.** Identifies sequences assigned to a community where their arm is not the dominant arm (1–2 sequences only).

**Key metrics.** 29 singletons and 7 doubletons across 20 communities involving 25 arms.

**Conclusion.** These represent rare structural variants or borderline cluster assignments, not systematic errors.

## 10. Within-arm heterogeneity analysis

**What it does.** Tests whether sequences from the same chromosome arm form a homogeneous cluster or show internal heterogeneity — specifically, whether some sequences are more similar to a foreign arm's sequences than to their own (cross-arm affinity).

**Mechanism caveat** (applies throughout §10): These analyses detect outcomes of sequence homogenization but cannot distinguish among mechanisms — gene conversion, unequal crossover, reciprocal exchange, or other processes may all contribute.

### 10.1 Allele vs paralog distance

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

### 10.2 Arm separation (silhouette analysis)

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

### 10.3 Cross-arm affinity

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

### 10.4 Population structure in cross-arm affinity

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

Notable patterns: chr16_q (C3) cross-arm is 70% AFR (60/86), consistent with the f7501 AFR-enrichment in §6.1.1. chr4_q (C1) cross-arm is AFR-enriched (52/146 = 36% vs baseline 28.9%). chrX_p (C15) self-arm is entirely AFR (18/18), indicating that the rare non-cross-arm PAR1 haplotypes are exclusively African.

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

### 10.5 Subtelomeric type discordance

**What it does.** Measures how often individuals are heterozygous for subtelomeric structural type — carrying one haplotype that matches the arm's own subtelomeric sequence (self-arm) and one that resembles a foreign chromosome's subtelomere (cross-arm). High discordance indicates that subtelomeric exchange is a common segregating polymorphism rather than a rare event.

**How.** Each individual's two haplotypes at a given arm are classified independently (§6.2) as self-arm or cross-arm, then paired: **concordant** = both self-arm or both cross-arm; **discordant** = one of each. Discordance rate = discordant / (concordant + discordant) for individuals where both haplotypes are classified.

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

### 10.6 Region length differences

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

### 10.7 Gene repertoire replacement scoring

**What it does.** Quantifies the fraction of a cross-arm sequence's gene content that matches its affinity arm rather than its own arm. A score of 1.0 means complete gene repertoire overlap.

**How.** For each cross-arm sequence, the fraction of genes matching the affinity arm vs own arm is computed (see mechanism caveat, §10 intro).

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

### 10.8 TAR1 in cross-arm vs self-arm sequences

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

### 10.9 Telomere-proximal vs telomere-distal island distribution by exchange status

**What it does.** Tests whether cross-arm exchange shifts the positional distribution of internal (TTAGGG)n islands.

**How.** Islands classified as proximal (centromere-ward half of PHR) or distal (telomere-ward half). Chi-squared test for cross-arm vs self-arm distribution.

**Key metrics.** Cross-arm: 45.9% proximal (1,494/3,256). Self-arm: 45.9% proximal (6,922/15,096). Chi-squared = 0.00, p = 0.99, OR = 1.00.

**Result.** No difference in island positional distribution between cross-arm and self-arm sequences. The proximal/distal ratio is identical (45.9%) regardless of exchange status.

**Conclusion.** Cross-arm exchange does not shift the positional distribution of internal (TTAGGG)n islands within the PHR. The island distribution reflects the underlying duplicon architecture rather than exchange history.

### 10.10 Within-community Jaccard distance structure

**What it does.** Tests whether within-community pairwise distances show bimodal structure (allelic vs paralog peaks).

**How.** Pairwise Jaccard distances sampled (up to 4,950 pairs per community) and examined for multi-modal structure.

**Key metrics.** C2: peaks near 0.02 and 0.37 (median allelic 0.023, median paralog 0.365). C12: peaks near 0.04 and 0.73 (median allelic 0.037, median paralog 0.727). Low-distance peak (0.01–0.05) = allelic variation; high-distance peaks = inter-arm paralog distances.

**Result.** Communities with well-separated bimodal peaks (C2, C12) have distinct arm-specific sequence types. Communities with diffuse distributions (C1, C7) show blurring of arm identity consistent with active homogenization.

**Conclusion.** Qualitatively consistent with Ambrosini et al.'s (2007) bimodal identity distribution (peaks at 98% and 91%), though Jaccard distances differ from sequence identity due to the k-mer-based metric.

### 10.11 Two-domain subtelomeric model test

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

---

## 11. Summary of key biological findings

1. **Subtelomeric regions form discrete communities of inter-chromosomal similarity.**
   - *What*: 41 chromosome arms cluster into 15 communities (arm level) and 50 (sequence level).
   - *Key metrics*: 232 individuals, 465 near-complete assemblies, 15,668 PHR sequences.
   - *Conclusion*: The community structure reflects shared duplicon content acquired through inter-chromosomal exchange and is consistent across all samples.

2. **Three categories of subtelomeric architecture.**
   - *What*: (a) Homogeneous arms with unique content (8/41 with 0% cross-arm). (b) Polymorphic arms with multiple subtelomeric types (34/41 spanning 2+ seq-level communities). (c) Fully interchangeable inter-chromosomal sharing (7/41 with 100% cross-arm).
   - *Conclusion*: This classification is novel and extends the qualitative patchwork model of Mefford & Trask (2002) with a quantitative framework.

3. **Recurrent inter-chromosomal exchange.**
   - *What*: At the arm level, 15.9% of sequences (2,484/15,668) show cross-arm affinity. At the sequence level (50 communities), 11.1% (1,740/15,668) remain cross-arm. Discordance rates up to 47.5%.
   - *Key metrics*: Linardopoulou et al. (2005) established that "human subtelomeres are hot spots of interchromosomal recombination and segmental duplication."
   - *Conclusion*: The data quantifies this at population scale. Discordance demonstrates this is a segregating polymorphism, though individual exchange events cannot be dated.

4. **Extensive gene repertoire overlap at acrocentric and sex chromosome subtelomeres.**
   - *What*: Acrocentric p-arms: chr13_p achieves complete gene overlap (score 1.000) toward all other acrocentric p-arms. PAR1/PAR2: complete gene overlap (score 1.000).
   - *Key metrics*: chr14_p (N=229) shows 83.0% cross-arm affinity with other acrocentric p-arms.
   - *Conclusion*: Consistent with recurrent inter-chromosomal exchange at rDNA-adjacent regions and obligate meiotic recombination at PAR.

5. **Population-specific exchange histories.**
   - *What*: 10 arm/community pairs show significant superpopulation bias (Fisher's exact, BH-corrected; 10 unique arms).
   - *Key metrics*: chr16_q cross-arm 70% AFR (60/86); chrX_p self-arm 82% AFR (18/22).
   - *Conclusion*: Consistent with population-specific exchange, but small sample sizes (chrY_p N=10, chr16_q N=86) and alternative explanations (drift, ILS) cannot be excluded without outgroup data.

6. **Subtelomeric gene content dominated by pseudogenes.**
   - *What*: Protein-coding 4–9% in most communities; PAR1 exception at 32.1%.
   - *Key metrics*: RPL23AP45 spans 10 communities/21 arms; SEPTIN14P22 spans 9 communities/22 arms. DUX4L specific to C1, MTCO to C7, SHOX to C15.
   - *Conclusion*: Consistent with telomere position effect (Mefford & Trask 2002, citing Baur et al. 2001).

7. **TAR1 as subtelomeric marker.**
   - *What*: Present in 94.6% of sequences and all 41 arms; virtually absent from PAR1.
   - *Key metrics*: TAR1 prevalence consistent across 1Mb analysis.
   - *Conclusion*: Consistent with different recombination mechanisms at PAR vs non-PAR subtelomeres. Whether TAR1 facilitates exchange or is co-transferred remains unclear.

8. **3D genome organization mirrors sequence communities.**
   - *What*: Four technologies, three cell types, and two organisms confirm community-structured 3D co-localization.
   - *Key metrics*: Human Hi-C B/W 0.027–0.074 (8 datasets); RPE-1 B/W 0.017–0.025 (3 datasets); Dip-C T2T 6.9% closer (16 cells, Fisher p=2.4e-05, Mantel rho=0.30); Sperm W/B=0.401 (60% closer, 20 cells); Mouse Hi-C per-arm-pair rho=0.63–0.86, Mantel rho=-0.39 to -0.49 (4 meiotic stages, 1Mb windows, per-haplotype, PHR coords); 2Mb: Mantel rho=-0.47 to -0.52; 4Mb: Mantel rho=-0.58 to -0.64. Negative control: 7 non-sharing arms 11% farther in GM12878, 40% farther in sperm. Community-free arm-level: Dip-C rho=−0.34, human Hi-C rho=0.66–0.83.
   - *Conclusion*: Sequence similarity correlates with 3D proximity across all tested samples, cell types, and organisms.

9. **Flanking unique-sequence regions show equal or stronger 3D signal.**
   - *What*: 100 kb flanking regions show stronger clustering than duplicated PHR (HG002: flanking (B/W=0.002) vs PHR (B/W=0.027) enrichment for HG002).
   - *Conclusion*: Rules out multi-mapping artifact. 3D clustering involves broader chromosomal domains. Multi-mapping suppresses PHR signal; flanking shows the true signal magnitude.

10. **C4 (chr7_q/chr12_q) as minimal-PHR positive control.**
    - *What*: Significant in 4/5 diploid Hi-C samples despite only 5–25 kb PHR and zero gene annotations.
    - *Conclusion*: Even minimal shared subtelomeric sequence coincides with 3D co-localization. The 3D signal is not driven by gene content.

11. **Community-specific 3D predictions confirmed.**
    - *What*: C7 co-localizes with nucleolar association. C1 peripheral, consistent with D4Z4-proximal lamin tethering (Masny et al. 2004). C14/C15 strongest in male samples.
    - *Key metrics*: Singleton communities (C8, C10, C13) show enrichment in highest-depth samples — this reflects homolog-homolog contacts, not inter-chromosomal co-localization.

12. **Two-domain subtelomeric model supported genome-wide.**
    - *What*: Flint/Mefford gradient confirmed on 39/48 arms, breakpoint structure on 39/41 testable arms, in 99.7% of individual haplotype sequences.
    - *Key metrics*: TTAGGG co-localization within 25 kb on 11/19 arms. Breakpoints 10–445 kb.
    - *Conclusion*: Extends the model from the handful of arms originally characterized by FISH to the entire human chromosome complement.

---

## 12. 3D genome validation: Hi-C and Pore-C

**What it does.** Tests whether chromosome arms in the same Leiden community show elevated inter-chromosomal Hi-C/Pore-C contacts at their subtelomeric PHR intervals — i.e., whether sequence-defined communities have a physical counterpart in 3D nuclear organization.

**How.** Per-haplotype analysis (maternal/paternal kept as separate arms, giving ~75 arms for diploid samples), 50 kb resolution. Three complementary statistical tests are applied:

1. **Within/Between ratio (W/B)** and **bootstrap permutation test**: For each community, the mean Hi-C contact between all arm pairs within the community is compared to the mean contact expected by chance. The null distribution is generated by randomly permuting community labels 10,000 times and recomputing within-community contact each time. The p-value is the fraction of permutations where the random within-community contact equals or exceeds the observed value. P-values are corrected for multiple testing using the **Benjamini-Hochberg (BH) procedure**, which controls the false discovery rate (FDR) — the expected proportion of false positives among significant results — at q < 0.05.
2. **Global Mann-Whitney U test**: All within-community arm-pair contacts are pooled and compared against all between-community arm-pair contacts using a non-parametric rank test. This tests whether within-community contacts are systematically higher than between-community contacts across all communities simultaneously. The **W/B ratio** is the ratio of mean within-community contact to mean between-community contact.
3. **Mantel test** and **per-PHR-pair correlation**: Described in §12.3 and §12.5.

**Key metrics.** 339 bootstrap tests across 7 samples, plus 7 global tests and 7 Mantel tests. Script: `analyze_hic_communities.py`.

### 12.1 Community enrichment at PHR intervals

**What it does.** For each sample, tests whether arm pairs sharing a community label have higher inter-chromosomal Hi-C contact at their subtelomeric PHR regions than arm pairs in different communities.

**How.** 6 Hi-C samples — HG002, HG02559, HG00658, HG02148, NA19036 (diploid HPRC) + CHM13 (haploid) —, HG002 Pore-C and HG002 CiFi. For each sample, the script loads the `.mcool` contact matrix at 50 kb resolution, extracts contact values at PHR intervals (from `all-vs-all.1Mb.p95.id95.len.tsv`), assigns arms to the 15 Leiden communities (§6), and runs the bootstrap and global tests described above.

**Result.**

| Sample | Technology | B/W ratio | p-value |
|--------|-----------|-----------|---------|
| CHM13 | Hi-C | 0.071 | 6.0e-18 |
| HG002 | Hi-C | 0.027 | 4.0e-66 |
| HG002 | Pore-C | 0.056 | 3.9e-85 |
| HG002 | CiFi | 0.036 | 2.0e-74 |
| HG00658 | Hi-C | 0.056 | 7.6e-12 |
| HG02148 | Hi-C | 0.050 | 9.1e-05 |
| HG02559 | Hi-C | 0.074 | 9.4e-03 |
| NA19036 | Hi-C | 0.049 | 1.9e-07 |

**Conclusion.** All 8 datasets show significant community enrichment (all p < 0.01). B/W ratios < 1 indicate within-community contacts are higher than between-community contacts (ratio of between/within). HG002 Hi-C shows the strongest enrichment (B/W = 0.027).

### 12.2 Flanking region analysis (100kb centromere-ward of PHR)

**What it does.** Tests whether unique-sequence regions flanking PHR boundaries also show community-structured 3D clustering, serving as a control for multi-mapping artifact.

**Key metrics.** Flanking extraction: 15,666 sequences (using `--fasta-mapping` with per-haplotype FASTAs, 0 truncated). Flanking community detection: 2 communities (C1 = most arms, C2 = acrocentric p-arms).

*Flanking community-free per-sample results (Spearman correlation between flanking Jaccard similarity and Hi-C contact):*

| Sample | Technology | rho | p-value |
|--------|-----------|-----|---------|
| CHM13 | Hi-C | 0.136 | 7e-4 |
| HG002 | Hi-C | 0.131 | 9e-12 |
| HG02559 | Hi-C | 0.039 | ns |
| HG00658 | Hi-C | 0.050 | ns |
| HG02148 | Hi-C | NaN | — |
| NA19036 | Hi-C | NaN | — |
| HG002 | Pore-C | 0.038 | 0.049 |
| HG002 | CiFi | 0.034 | ns |

*Note: HG02148 and NA19036 report NaN because their assemblies have fragmented chromosomes (chr1, chr4, chr6, etc. as random contigs instead of single scaffolds), and all nonzero-Jaccard flanking pairs happen to involve those fragmented chromosomes.*

*Dip-C flanking overlay (radial positioning):*

| Cell type | Flanking particles | Non-flanking terminal | Mean radial (flank) | Mean radial (non-flank) | p-value |
|-----------|-------------------|----------------------|---------------------|------------------------|---------|
| GM12878 | 2,966 | 14,712 | 0.634 | 0.658 | 3.8e-9 |
| PBMC | 3,960 | 18,071 | 0.602 | 0.613 | 5.0e-3 |

**Conclusion.** Flanking regions are unique sequence (no multi-mapping risk). Community-free analysis shows significant positive correlations in CHM13 (rho=0.136) and HG002 Hi-C (rho=0.131), confirming that the similarity-contact relationship extends into flanking unique sequence. The weaker effect (compared to PHR) is expected since flanking regions have less inter-chromosomal similarity. Two assemblies (HG02148, NA19036) are untestable due to chromosome fragmentation.

### 12.3 Mantel test: sequence similarity vs Hi-C contact

**What it does.** Tests whether there is a continuous, graded relationship between arm-level sequence similarity and Hi-C contact frequency — not just a binary "within vs between" enrichment, but a proportional correlation where more similar arms have more contact.

**How.** The **Mantel test** (Mantel 1967) measures the correlation between two distance/similarity matrices. Here, two N×N matrices (N = number of shared arms, typically 31–38) are compared: (1) the arm-level Jaccard distance matrix from §5 (sequence similarity) and (2) the arm-level Hi-C contact matrix (3D proximity). The test computes a **Spearman rank correlation** (rho) between the upper-triangle elements of the two matrices, then assesses significance by permuting row/column labels 10,000 times and recomputing rho each time. The permutation p-value is the fraction of permutations where |rho_perm| ≥ |rho_obs|. Expected direction: **negative rho** — arms with small Jaccard distance (high similarity) should have high Hi-C contact, producing a negative distance-contact correlation.

**Result.** PHR regions:

| Sample | Technology | Mantel rho | Mantel p | n shared arms |
|--------|-----------|------------|----------|---------------|
| CHM13 | Hi-C | -0.656 | <0.0001 | 38 |
| HG002 | Hi-C | -0.657 | <0.0001 | 41 |
| HG02559 | Hi-C | -0.397 | 4.0e-04 | 37 |
| HG00658 | Hi-C | -0.276 | 9.3e-03 | 37 |
| HG02148 | Hi-C | -0.153 | 0.085 | 37 |
| NA19036 | Hi-C | -0.266 | 0.018 | 34 |
| HG002 | Pore-C | -0.486 | <0.0001 | 41 |
| HG002 | CiFi | -0.308 | <0.0001 | 41 |

7 of 8 datasets show significant negative Mantel correlation (all except HG02148, p=0.085). CHM13 and HG002 show the strongest correlations (rho ≈ -0.66).

*100 kb flanking regions:*

| Sample | Technology | Mantel rho | Mantel p | n shared arms |
|--------|-----------|------------|----------|---------------|
| CHM13 | Hi-C | -0.522 | <0.0001 | 38 |
| HG002 | Hi-C | -0.520 | <0.0001 | 41 |
| HG02559 | Hi-C | -0.323 | 0.002 | 37 |
| HG00658 | Hi-C | -0.288 | 0.005 | 37 |
| HG02148 | Hi-C | -0.127 | 0.124 | 37 |
| NA19036 | Hi-C | -0.226 | 0.032 | 34 |
| HG002 | Pore-C | -0.314 | <0.0001 | 41 |

**Conclusion.** Flanking Mantel rho values are comparable to PHR Mantel (CHM13: -0.522 flanking vs -0.656 PHR; HG002: -0.520 vs -0.657), confirming that the 3D proximity signal extends into unique-sequence regions centromere-ward of the PHR boundary. The flanking regions themselves contain no inter-chromosomal duplicated content, so this signal reflects broader chromosomal domain proximity rather than multi-mapping artifacts.

### 12.4 Independent Hi-C community detection

**What it does.** Tests whether Hi-C contacts, analyzed independently from sequence data, recover the same community structure as sequence similarity. This is the strongest test of concordance: two completely independent data sources (pangenome graph vs proximity ligation) are used to partition arms into communities, then the two partitions are compared.

**How.** For each sample: (1) export the O/E-normalized (observed/expected) inter-chromosomal contact matrix at PHR regions (50 kb, per-haplotype); (2) run Leiden community detection on this contact matrix (same algorithm as §6, but using Hi-C contacts instead of Jaccard similarity); (3) compare the resulting Hi-C communities to the 15 sequence-based communities using the **Adjusted Rand Index (ARI)**. The ARI measures agreement between two partitions of the same set of elements, corrected for chance: ARI = 1 means perfect agreement (identical partitions), ARI = 0 means agreement no better than random, and ARI < 0 means agreement worse than random. Unlike raw Rand Index, ARI adjusts for the number and size of clusters, so it is not inflated by partitions with many small clusters.

**Result.**

| Sample | Technology | ARI vs seq |
|--------|-----------|------------|
| CHM13 | Hi-C | 0.539 |
| HG002 | Hi-C | 0.296 |
| HG02559 | Hi-C | 0.128 |
| HG00658 | Hi-C | 0.133 |
| HG02148 | Hi-C | 0.123 |
| NA19036 | Hi-C | 0.165 |
| HG002 | Pore-C | 0.264 |
| HG002 | CiFi | 0.056 |

**Conclusion.** All 8 samples show positive ARI (>0), confirming non-random agreement. CHM13 shows highest concordance (ARI=0.539), followed by HG002 Hi-C (0.296) and Pore-C (0.264).

### 12.5 Per-PHR-pair sequence similarity vs 3D proximity

**What it does.** Correlates per-arm-pair sequence similarity with Hi-C contact at finer resolution than the Mantel test (§12.3). While the Mantel test uses a single arm-level distance matrix averaged across all 232 individuals, this test uses per-arm-pair mean Jaccard similarity computed from the full sequence-level pairwise comparisons, preserving more of the underlying variance.

**How.** For each pair of arms on different chromosomes, the mean Jaccard similarity is computed from all sequence pairs (e.g., 440 × 450 = 198,000 pairwise comparisons for two arms with ~440 sequences each). This mean Jaccard is then correlated with the Hi-C contact value for that arm pair using **Spearman rank correlation** (a non-parametric measure of monotonic association; rho = 1 means perfect positive monotonic relationship). Expected direction: positive rho — arm pairs with higher mean Jaccard similarity should have higher Hi-C contact.

**Result.**

| Sample | Technology | Arm pairs | Spearman rho | p-value |
|--------|-----------|-----------|-------------|---------|
| CHM13 | Hi-C | 688 | 0.674 | 3.0e-92 |
| HG002 | Hi-C | 803 | 0.657 | 1.7e-100 |
| HG02559 | Hi-C | 652 | 0.401 | 1.6e-26 |
| HG00658 | Hi-C | 652 | 0.281 | 2.4e-13 |
| HG02148 | Hi-C | 652 | 0.160 | 4.0e-05 |
| NA19036 | Hi-C | 550 | 0.263 | 3.8e-10 |
| HG002 | Pore-C | 803 | 0.485 | 1.6e-48 |
| HG002 | CiFi | 803 | 0.304 | 1.3e-18 |

**Conclusion.** All 8 samples show significant positive correlation. CHM13 and HG002 Hi-C strongest (rho ≈ 0.66), consistent with highest sequencing depth. Direction: higher Jaccard similarity → more Hi-C contact.

### 12.5b Sequence-level similarity vs Hi-C contact (community-free)

**What it does.** Correlates individual PHR sequence pairs directly — no community labels, no arm-level aggregation. For each pair of PHR sequences on different chromosomes, the Jaccard similarity (from the pangenome graph) is compared with the Hi-C contact at their exact genomic coordinates. This is the most direct test: do two specific subtelomeric sequences that share pangenome graph structure also tend to be in 3D proximity?

**How.** The PHR sequence names encode subranges of the assembly chromosomes used by the cooler (e.g., `HG002#1#CM088242.1:82897078-82942077_chr17_qarm` → `chr17_PATERNAL:82897078-82942077`). For each sample's 3D dataset, only that sample's PHR sequences are used. The cooler is queried at the exact coordinates of each inter-chromosomal sequence pair to obtain the balanced Hi-C contact sum, normalized by total region size (bp). Spearman rank correlation is computed across all inter-chromosomal pairs. Script: `sequence_hic_correlation.py`.

**Result.** PHR regions, all samples at 50 kb and 10 kb:

| Sample | Technology | Seq pairs | ρ (50kb) | p (50kb) | ρ (10kb) | p (10kb) |
|---|---|---|---|---|---|---|
| CHM13 | Hi-C | 652 | 0.716 | 1.2e-103 | 0.730 | 1.5e-115 |
| HG002 | Hi-C | 2,544 | 0.662 | <1e-300 | 0.667 | 1.7e-275 |
| HG02559 | Hi-C | 825 | 0.638 | 1.3e-95 | 0.709 | 2.1e-93 |
| HG00658 | Hi-C | 742 | 0.701 | 7.1e-111 | 0.739 | 9.9e-149 |
| HG02148 | Hi-C | 627 | 0.760 | 7.3e-119 | 0.809 | 2.4e-146 |
| NA19036 | Hi-C | 901 | 0.773 | 1.4e-179 | 0.827 | 7.7e-237 |
| HG002 | Pore-C | 2,830 | 0.381 | 1.2e-98 | 0.379 | 1.5e-97 |
| HG002 | CiFi | 2,757 | 0.191 | 3.7e-24 | 0.191 | 3.8e-24 |

RPE-1 (self-discovered similarity on full flanks):

| Dataset | Technology | Seq pairs | ρ (50kb) | p (50kb) | ρ (10kb) | p (10kb) |
|---|---|---|---|---|---|---|
| Async CiFi | CiFi | 4,048 | 0.315 | 3.6e-94 | 0.302 | 3.7e-86 |
| Async Pore-C | Pore-C | 4,048 | 0.435 | 4.1e-186 | 0.424 | 1.5e-176 |
| Mitotic CiFi | CiFi | 4,048 | 0.181 | 2.9e-31 | 0.185 | 1.8e-32 |

Mouse: see §15 for cross-species validation with 4 meiotic stages (Zuo et al. 2021).

**Conclusion.** Significant in ALL 8 datasets (all p < 3.7e-24). This community-free test confirms the relationship holds at individual sequence-pair resolution: specific subtelomeric sequences that share more pangenome graph nodes are in closer 3D proximity, without relying on any discrete community assignment. Hi-C correlations are strong (ρ = 0.64–0.83 across 6 samples at 10kb), with the strongest in NA19036 and HG02148. Pore-C is intermediate (ρ = 0.38), CiFi weakest (ρ = 0.19). Results are consistent across 50kb and 10kb resolutions, with slight strengthening at 10kb. The RPE-1 correlations (ρ = 0.18–0.44) are comparable to HG002 CiFi/Pore-C, confirming the signal generalizes across cell types. Mitotic attenuation (~2x) mirrors the community-based results (§12.14.5).

### 12.6 Sequence-level community validation (Hi-C)

**What it does.** Tests community-structured 3D contacts using the finer 50-community partition (§6.2) instead of the 15 arm-level communities. This tests whether the within-arm polymorphism captured by sequence-level communities also has a 3D counterpart. For a given individual, each haplotype-arm maps to one of the 50 sequence-level communities — many of which will be singletons (only 1 arm in that individual), reducing statistical power but increasing specificity.

**Result.** PHR regions:

| Sample | Technology | Seq communities | Seq global W/B p | Seq significant (BH q<0.05) |
|--------|-----------|----------------|------------------|----------------------------|
| CHM13 | Hi-C | 31 (4 multi-arm) | 0.217 | 1/4 |
| HG002 | Hi-C | 34 (27 multi-arm) | 3.0e-38 | 20/27 |
| HG02559 | Hi-C | 33 (22 multi-arm) | 0.035 | 7/22 |
| HG00658 | Hi-C | 32 (26 multi-arm) | 9.4e-06 | 0/26 |
| HG02148 | Hi-C | 34 (19 multi-arm) | 1.4e-06 | 4/19 |
| NA19036 | Hi-C | 36 (19 multi-arm) | 1.2e-05 | 9/19 |
| HG002 | Pore-C | 34 (27 multi-arm) | 2.5e-51 | 25/27 |

100kb flanking regions:

| Sample | Technology | Seq communities | Seq global W/B p | Seq significant (BH q<0.05) |
|--------|-----------|----------------|------------------|----------------------------|
| CHM13 | Hi-C | 31 (4 multi-arm) | 0.675 | 1/4 |
| HG002 | Hi-C | 34 (27 multi-arm) | 2.1e-127 | 27/27 |
| HG02559 | Hi-C | 33 (22 multi-arm) | 3.9e-07 | 8/22 |
| HG00658 | Hi-C | 32 (26 multi-arm) | 1.1e-19 | 12/26 |
| HG02148 | Hi-C | 34 (19 multi-arm) | 6.3e-17 | 10/19 |
| NA19036 | Hi-C | 36 (19 multi-arm) | 5.3e-18 | 11/19 |
| HG002 | Pore-C | 34 (27 multi-arm) | 2.2e-40 | 26/27 |

### 12.7 Resolution sensitivity (10 kb)

**What it does.** Tests whether the 3D signal is robust to resolution choice (10 kb vs 50 kb).

**How.** Re-analysis at 10 kb (5x finer; PHR median 105 kb = ~10 bins at 10 kb vs ~2 bins at 50 kb).

**Result.** PHR regions at 10 kb:

| Sample | Technology | B/W (10kb) | B/W (50kb) | Global p (10kb) |
|--------|-----------|------------|------------|-----------------|
| CHM13 | Hi-C | 0.072 | 0.071 | 1.2e-22 |
| HG002 | Hi-C | 0.028 | 0.027 | 4.0e-111 |
| HG02559 | Hi-C | 0.081 | 0.074 | 6.3e-03 |
| HG00658 | Hi-C | 0.055 | 0.056 | 3.6e-12 |
| HG02148 | Hi-C | 0.053 | 0.050 | 8.1e-06 |
| NA19036 | Hi-C | 0.065 | 0.049 | 1.2e-08 |
| HG002 | Pore-C | 0.055 | 0.056 | 2.6e-92 |

B/W ratios are comparable at 10 kb and 50 kb resolutions, confirming the signal is robust to resolution choice. All 7 datasets significant at both resolutions.

**Full multi-resolution B/W ratios (all 5 mcool resolutions):**

| Sample | 5kb | 10kb | 20kb | 50kb | 100kb |
|--------|-----|------|------|------|-------|
| CHM13 | 0.075 | 0.072 | 0.074 | 0.071 | 0.072 |
| HG002 | 0.029 | 0.028 | 0.028 | 0.027 | 0.021 |
| HG002 Pore-C | 0.058 | 0.055 | 0.054 | 0.056 | 0.049 |
| HG002 CiFi | 0.048 | 0.045 | 0.041 | 0.036 | 0.042 |
| HG02559 | 0.078 | 0.081 | 0.086 | 0.074 | 0.054 |
| HG00658 | 0.056 | 0.055 | 0.053 | 0.056 | 0.032 |
| HG02148 | 0.052 | 0.053 | 0.054 | 0.050 | 0.038 |
| NA19036 | 0.063 | 0.065 | 0.060 | 0.049 | 0.040 |

All B/W ratios remain well below 1.0 at every resolution, demonstrating that community enrichment is resolution-invariant. The slight increase at 100kb for some samples reflects reduced power as PHR regions (median 105 kb) span only ~1 bin.

Comparison of significant multi-arm communities at 10 kb vs 50 kb:

| Sample | Technology | 50 kb | 10 kb | Shared | Only 50 kb | Only 10 kb |
|--------|-----------|-------|-------|--------|------------|------------|
| CHM13 | Hi-C | 0 | 0 | — | — | — |
| HG002 | Hi-C | 11 | 9 | C2,C4,C6,C7,C10,C14,C15 (7) | C3,C11,C12,C13 (4) | C1,C9 (2) |
| HG02559 | Hi-C | 6 | 9 | C1,C4,C5,C7,C10,C11 (6) | — | C2,C12,C13 (3) |
| HG00658 | Hi-C | 0 | 6 | — | — | C4,C5,C8,C9,C11,C12 (6) |
| HG02148 | Hi-C | 4 | 4 | C9,C14,C15 (3) | C4 (1) | C2 (1) |
| NA19036 | Hi-C | 4 | 8 | C1,C2,C4,C14 (4) | — | C3,C5,C8,C11 (4) |
| HG002 | Pore-C | 11 | 12 | C1,C3,C4,C6,C7,C10,C11,C12,C13,C14,C15 (11) | — | C5 (1) |

**Conclusion.** The two resolutions identify largely overlapping sets. W/B ratios increase for 5/7 samples at 10 kb. Net effect: more detections at 10 kb (38 vs 26 at 50 kb). Core communities (C4, C14, C15 in males; C1, C2 across populations) are significant at both resolutions, confirming robustness.

100 kb flanking regions at 10 kb:

| Sample | Technology | B/W (10kb) | B/W (50kb) | Global p (10kb) |
|--------|-----------|------------|------------|-----------------|
| CHM13 | Hi-C | 0.071 | 0.057 | 8.0e-07 |
| HG002 | Hi-C | 0.001 | 0.002 | 1.1e-48 |
| HG02559 | Hi-C | 0.005 | 0.008 | 2.0e-03 |
| HG00658 | Hi-C | 0.003 | 0.006 | 4.8e-13 |
| HG02148 | Hi-C | 0.003 | 0.005 | 1.2e-05 |
| NA19036 | Hi-C | 0.003 | 0.006 | 2.2e-05 |
| HG002 | Pore-C | 0.032 | 0.034 | 5.0e-27 |

**Full multi-resolution flanking B/W ratios (all 5 mcool resolutions):**

| Sample | 5kb | 10kb | 20kb | 50kb | 100kb |
|--------|-----|------|------|------|-------|
| CHM13 | 0.062 | 0.071 | 0.057 | 0.057 | 0.052 |
| HG002 | 0.001 | 0.001 | 0.001 | 0.002 | 0.003 |
| HG002 Pore-C | 0.030 | 0.032 | 0.033 | 0.034 | 0.039 |
| HG02559 | 0.004 | 0.005 | 0.006 | 0.008 | 0.014 |
| HG00658 | 0.003 | 0.003 | 0.004 | 0.006 | 0.010 |
| HG02148 | 0.003 | 0.003 | 0.004 | 0.005 | 0.008 |
| NA19036 | 0.003 | 0.003 | 0.004 | 0.006 | 0.010 |

Flanking B/W ratios are comparable or stronger at finer resolutions (HG002 Hi-C: 0.001 at 5-10kb vs 0.002-0.003 at 50-100kb). All datasets show significant flanking enrichment at all 5 resolutions, confirming that the community-structured 3D signal extends to unique-sequence flanking regions.

### 12.8 p-arm vs q-arm 3D enrichment

**What it does.** Tests whether the 3D signal differs between p-arm-dominated and q-arm-dominated communities. Human chromosomes have distinct p-arm (short) and q-arm (long) subtelomeric architectures. If one arm type had systematically stronger 3D clustering, it could suggest a mechanistic difference.

**How.** Communities classified by arm composition: q-dominated (≥75% q-arms), p-dominated (≥75% p-arms), or mixed. Per-sample W/B ratios compared across categories using **Mann-Whitney U test** — a non-parametric test comparing the distributions of two independent groups.

**Key metrics.** q-dominated median enrichment 12.0x; p-dominated 10.2x; mixed 7.9x. Mann-Whitney U (q vs p) = 812, p = 0.90.

**Conclusion.** No significant difference — the 3D signal does not differ systematically between p-arm and q-arm regions, consistent with subtelomeric exchange operating at both chromosome ends.

### 12.9 Acrocentric exclusion confound control

**What it does.** Tests whether the 3D signal is driven by nucleolar association (acrocentric chromosomes constitutively associate with nucleoli) rather than community structure.

**How.** All acrocentric chromosome arms (chr13, 14, 15, 21, 22 — both p and q) excluded; analysis re-run at 50 kb.

**Result.**

| Sample | Technology | B/W (all arms) | B/W (no acro) | Global p (no acro) |
|--------|-----------|----------------|---------------|---------------------|
| CHM13 | Hi-C | 0.071 | 0.079 | 2.0e-16 |
| HG002 | Hi-C | 0.027 | 0.027 | 5.3e-66 |
| HG02559 | Hi-C | 0.074 | 0.069 | 5.1e-05 |
| HG00658 | Hi-C | 0.056 | 0.052 | 2.5e-19 |
| HG02148 | Hi-C | 0.050 | 0.051 | 1.6e-06 |
| NA19036 | Hi-C | 0.049 | 0.051 | 2.2e-14 |
| HG002 | Pore-C | 0.056 | 0.086 | 4.0e-59 |

**Full multi-resolution no-acrocentric B/W ratios (all 5 mcool resolutions):**

| Sample | 5kb | 10kb | 20kb | 50kb | 100kb |
|--------|-----|------|------|------|-------|
| CHM13 | 0.082 | 0.081 | 0.081 | 0.079 | 0.078 |
| HG002 | 0.032 | 0.031 | 0.030 | 0.027 | 0.021 |
| HG002 Pore-C | 0.079 | 0.071 | 0.073 | 0.086 | 0.088 |
| HG02559 | 0.074 | 0.077 | 0.081 | 0.069 | 0.051 |
| HG00658 | 0.050 | 0.049 | 0.048 | 0.052 | 0.029 |
| HG02148 | 0.054 | 0.055 | 0.056 | 0.051 | 0.039 |
| NA19036 | 0.065 | 0.067 | 0.061 | 0.051 | 0.041 |

All B/W ratios remain well below 1.0 and significant at every resolution, confirming that the no-acrocentric signal is resolution-invariant.

**Conclusion.** The global 3D signal persists in all 7 datasets when acrocentric arms are excluded (all p < 1e-05). B/W ratios are comparable with and without acrocentrics (HG002 Hi-C: 0.027 → 0.027; most samples similar or slightly higher). The acrocentric exclusion demonstrates that community-structured 3D contacts are not an artifact of nucleolar co-localization. This holds at all 5 resolutions (5kb-100kb).

### 12.10 O/E inter-chromosomal normalization

**What it does.** Controls for chromosome-level contact biases. In Hi-C data, smaller chromosomes tend to intermingle more (they occupy less nuclear volume), and chromosomes with higher mappability generate more contacts. These biases could inflate the W/B ratio if community members happen to be small chromosomes. The **O/E (observed/expected) normalization** removes these biases.

**How.** For each pair of arms (i, j), the expected contact E_ij is computed from the arm's marginal contact frequencies: E_ij = (row_sum_i × col_sum_j) / total_inter, where row_sum_i is the total inter-chromosomal contact of arm i with all other arms, and total_inter is the sum of all inter-chromosomal contacts. The O/E ratio = observed / expected. An O/E > 1 means the arm pair has more contact than expected from their marginal frequencies. The within-community and between-community O/E values are then compared as in §12.1.

**Result.**

| Sample | Technology | O/E within | O/E between | O/E ratio | Raw B/W (§12.1) |
|--------|-----------|------------|-------------|-----------|-----------------|
| CHM13 | Hi-C | 0.176 | 0.014 | 12.9x | 0.071 |
| HG002 | Hi-C | 0.068 | 0.002 | 34.4x | 0.027 |
| HG02559 | Hi-C | 0.101 | 0.007 | 14.4x | 0.074 |
| HG00658 | Hi-C | 0.154 | 0.012 | 13.0x | 0.056 |
| HG02148 | Hi-C | 0.134 | 0.008 | 17.3x | 0.050 |
| NA19036 | Hi-C | 0.153 | 0.008 | 19.2x | 0.049 |
| HG002 | Pore-C | 0.033 | 0.004 | 8.6x | 0.056 |

**Conclusion.** Small chromosomes tend to intermingle more in the nucleus, which could inflate the raw enrichment if community members happen to be on small chromosomes. The O/E normalization removes this bias by dividing each contact by the expected value based on chromosome sizes. After normalization, within-community contacts are still 8.6–34.4x higher than between-community contacts, confirming that the enrichment reflects genuine community-specific proximity, not chromosome-size effects.

### 12.11 Per-individual cross-arm affinity vs 3D enrichment

**What it does.** Tests whether the 3D signal depends on intact duplicated sequence — specifically, whether communities with low discordance show stronger Hi-C signal.

**Key metrics.** 54 sample × community pairs (5 diploid samples × multi-arm communities): Spearman rho = −0.31 (p = 0.024). Significant communities have lower mean discordance (0.087) than non-significant (0.139; Wilcoxon p = 0.013).

**Conclusion.** Consistent with the 3D signal depending on intact duplicated sequence: concordant (both haplotypes self-arm) generates the contact signal; discordant (one exchanged) dilutes it.

### 12.12 Dip-C lymphocyte vs monocyte cell-type split

**What it does.** Tests whether the 3D signal differs between lymphocytes and monocytes in PBMC Dip-C data.

**How.** 18 PBMC Dip-C cells (Tan et al. 2018): 15 lymphocytes and 3 monocytes/neutrophils (cells 9, 14, 18; identified in Tan et al. Fig 4D).

**Result (all 18 cells combined).**

| Metric | Value |
|--------|-------|
| W/B ratio | 0.983 (1.7% closer within-community) |
| Wilcoxon p | 0.305 (not significant) |
| Fisher combined p | 0.217 (not significant) |

*Data stored at `dipc_t2t/pbmc_hg19/` (18 cells, hg19 coordinates, PHR boundaries projected to hg19 via impg). Monocytes: cells 09, 14, 18 (Tan et al. 2018 Fig 3); lymphocytes: remaining 15 cells.*

**Conclusion.** PBMC cells (all 18 combined) do not show significant community-structured clustering at PHR-specific coordinates (W/B = 0.983, p = 0.305). The non-significance reflects (1) hg19 coordinate noise from PHR projection, (2) mixed cell-type composition (monocytes dilute signal), (3) small N = 18 cells, and (4) PHR-specific regions are smaller than the default 300kb window used in earlier analysis. The GM12878 T2T result (6.9% closer, p = 2.4e-5) uses native T2T coordinates without projection noise and remains the primary Dip-C finding.

### 12.13 Nuclear lamina cross-reference

**What it does.** Tests whether community radial positioning (Dip-C) correlates with Hi-C enrichment.

**Key metrics.** C1 (D4Z4) radial = 0.717 (peripheral), consistent with D4Z4-proximal lamin A/C tethering (Masny et al. 2004; Ottaviani et al. 2009). C14 (PAR2) most peripheral (0.839). C10 (chr17_p) most interior (0.520). C6 interior (0.572).

**Conclusion.** Lamina association may contribute to but is not the primary correlate of 3D community enrichment.

### 12.14 Cell-type validation: RPE-1 (CiFi + Pore-C)

**What it does.** Tests whether the subtelomeric community signal generalizes beyond LCLs to a different cell type, and whether it is modulated by cell cycle state.

**How.** Proximity ligation data from RPE-1 (hTERT-RPE-1), a near-diploid, non-transformed, telomerase-immortalized retinal pigment epithelial cell line, aligned to the RPE-1 diploid assembly (RPE1v1.1). Three datasets, enabling three independent comparisons:

| Dataset | Technology | Cell state | Contacts | Comparison |
|---|---|---|---|---|
| Async CiFi | PacBio CiFi | Asynchronous (interphase) | 27.5M | Cell-type control vs LCL |
| Async Pore-C | ONT Pore-C | Asynchronous (interphase) | 48.7M | Cross-platform (PacBio vs ONT) |
| Mitotic CiFi | PacBio CiFi | Mitotic arrest | 16.9M | Cell-cycle modulation |

Same unified 5-step workflow as human samples: O/E matrix export, independent Hi-C community detection (Leiden), ARI comparison, bootstrap permutation (10,000), Mantel test, per-PHR-pair correlation. 50 kb resolution, per-haplotype (92 arms = 46 chromosomes × 2 haplotypes). Sequence communities reused from HPRC (15 arm-level Leiden communities, §6).

#### 12.14.1 Community enrichment

**Key metrics.** Within-community vs between-community contact ratio (W/B), Mann-Whitney U p-value, number of significant communities (BH q<0.05).

**Result.**

**HPRC-community B/W across all resolutions:**

| Dataset | 5kb | 10kb | 20kb | 50kb | 100kb |
|---------|-----|------|------|------|-------|
| Async CiFi | 0.052 | 0.033 | 0.032 | 0.024 | 0.015 |
| Async Pore-C | 0.048 | 0.048 | 0.038 | 0.031 | 0.019 |
| Mitotic CiFi | 0.043 | 0.030 | 0.022 | 0.008 | 0.005 |

**RPE-1 self-discovered community B/W (§15b pipeline):**

| Dataset | 5kb | 10kb | 20kb | 50kb | 100kb |
|---------|-----|------|------|------|-------|
| Async CiFi | 0.031 | 0.022 | 0.021 | 0.012 | 0.008 |
| Async Pore-C | 0.034 | 0.032 | 0.032 | 0.022 | 0.019 |
| Mitotic CiFi | 0.027 | 0.022 | 0.014 | 0.005 | 0.003 |

The self-discovered communities show stronger enrichment (B/W 0.003–0.034 vs 0.005–0.052) because they capture RPE-1-specific sharing patterns directly from the cell line's own genome. All analyses use RPE-1's own PHR boundaries (`rpe1.all-vs-all.p95.id95.len.tsv`) for coordinate precision.

**RPE-1 flanking (100kb centromere-ward of PHR) B/W across all resolutions:**

| Dataset | 5kb | 10kb | 20kb | 50kb | 100kb |
|---------|-----|------|------|------|-------|
| Async CiFi | 0.006 | 0.005 | 0.005 | 0.006 | 0.006 |
| Async Pore-C | 0.007 | 0.010 | 0.007 | 0.014 | 0.018 |
| Mitotic CiFi | 0.020 | 0.024 | 0.015 | 0.015 | 0.019 |

All datasets show significant flanking enrichment at every resolution, with B/W ratios well below those for PHR regions.

*Flanking community-free:* rho ≈ 0 for all 3 datasets (negative control confirmed).

*Flanking pggb:* 4,625 similarity pairs (68 sequences).

**Conclusion.** All 3 RPE-1 datasets show significant community enrichment (all p < 1.6e-65). Async CiFi shows the strongest enrichment (B/W = 0.017). Flanking community-free correlations near zero confirm that the PHR signal is driven by shared subtelomeric sequence, not a general property of chromosome ends.

#### 12.14.2 Mantel test: sequence similarity vs contact

**Key metrics.** Spearman correlation between arm-level Jaccard distance and Hi-C contact (40 shared arms).

**Result.**

| Dataset | Mantel rho | Mantel p |
|---|---|---|
| Async CiFi | -0.457 | <1e-300 |
| Async Pore-C | -0.611 | <1e-300 |
| Mitotic CiFi | -0.340 | <1e-300 |
| *HG002 Hi-C (LCL, §12.3)* | *-0.409* | *<0.0001* |
| *HG002 Pore-C (LCL, §12.3)* | *-0.496* | *<0.0001* |

**Conclusion.** All 3 RPE-1 datasets show significant negative correlation: arms with more similar subtelomeric sequence have more inter-chromosomal contact. Async datasets (rho = -0.46 to -0.61) are comparable to or stronger than LCL. Mitotic is attenuated ~1.3x (rho = -0.34).

#### 12.14.3 Independent Hi-C community detection (ARI)

**What it does.** Tests whether Hi-C contacts in RPE-1 independently recover the same community structure as sequence similarity.

**How.** Leiden community detection on the O/E-normalized inter-chromosomal contact matrix (50 kb, per-haplotype). Compared to the 15 sequence-based communities via Adjusted Rand Index (ARI=1: perfect, ARI=0: random).

**Result.**

| Dataset | Contact matrix arms | Hi-C communities | ARI vs sequence |
|---|---|---|---|
| Async CiFi | 92 | 41 | 0.283 |
| Async Pore-C | 92 | 42 | 0.266 |
| Mitotic CiFi | 92 | 44 | 0.283 |
| *HG002 Hi-C (LCL, §12.4)* | *75* | *46* | *0.277* |
| *HG002 Pore-C (LCL, §12.4)* | *75* | *42* | *0.295* |

**Conclusion.** All RPE-1 datasets show positive ARI (0.27–0.28), comparable to HG002 LCL (0.28–0.30). Hi-C contacts independently recover community structure in RPE-1 at a level consistent with LCL data.

#### 12.14.4 Per-PHR-pair correlation: sequence similarity vs contact

**What it does.** Correlates per-arm-pair sequence similarity with Hi-C contact, using the full sequence-level pairwise Jaccard similarity.

**How.** Same approach as §12.5. For each pair of arms on different chromosomes, mean Jaccard similarity compared with Hi-C contact.

**Result.**

| Dataset | Arm pairs | Spearman rho | p-value |
|---|---|---|---|
| Async CiFi | 652 | 0.538 | 4.2e-50 |
| Async Pore-C | 652 | 0.681 | 5.3e-90 |
| Mitotic CiFi | 652 | 0.389 | 4.9e-25 |
| *HG002 Hi-C (LCL, §12.5)* | *764* | *0.414* | *5.5e-33* |
| *HG002 Pore-C (LCL, §12.5)* | *764* | *0.495* | *2.1e-48* |

**Conclusion.** All 3 RPE-1 datasets show significant positive correlation between sequence similarity and contact frequency. Async datasets (rho = 0.54–0.68) are stronger than LCL, reflecting RPE-1-specific PHR boundaries. Mitotic is attenuated ~1.4x (rho = 0.39, p = 4.9e-25).

#### 12.14.5 Cell-cycle modulation: async vs mitotic

**What it does.** Compares the same cell line (RPE-1) across cell cycle states to assess the effect of mitotic chromosome condensation on the community signal.

**How.** Direct comparison of async CiFi (interphase) vs mitotic CiFi (mitotic arrest), both PacBio CiFi on the same reference.

**Result.**

| Metric | Async CiFi | Mitotic CiFi | Change |
|---|---|---|---|
| W/B ratio | 35.8x | 33.7x | ~unchanged |
| Significant communities | 13/15 | 13/15 | identical |
| Mantel rho | -0.457 | -0.340 | ~1.3x weaker |
| PHR pair rho | 0.538 | 0.389 | ~1.4x weaker |

**Conclusion.** Community-level enrichment (binary: within vs between) is maintained in mitotic cells — the same 13/15 communities are significant, and the W/B ratio is comparable. The quantitative similarity-contact correlation (Mantel, PHR pair rho) is attenuated ~2x. This dissociation indicates that coarse community structure persists through mitosis while the fine-grained proportional relationship between sequence similarity and contact frequency weakens. Caveat: synchronization efficiency is unknown — attenuation may partly reflect a mixed interphase/mitotic population.

#### 12.14.6 Cross-platform reproducibility

**What it does.** Compares PacBio CiFi and ONT Pore-C on the same cell line and cell state.

**Result.**

| Metric | Async CiFi (PacBio) | Async Pore-C (ONT) |
|---|---|---|
| W/B ratio | 35.8x | 31.4x |
| Mantel rho | -0.457 | -0.611 |
| PHR pair rho | 0.538 | 0.681 |
| ARI | 0.283 | 0.266 |

**Conclusion.** Results are consistent across platforms. Pore-C shows slightly stronger Mantel and PHR-pair correlations, consistent with multi-way contacts amplifying signal (same pattern as HG002 Hi-C vs Pore-C in §12.1).

### 12.15 Per-community enrichment across samples and cell types

**What it does.** Compares per-community W/B ratios across all 10 datasets (6 human Hi-C, 1 human Pore-C, 3 RPE-1) to identify which communities are reproducible and which depend on sample, depth, or sex.

**How.** For each community, the ratio of observed mean within-community contact to permuted random mean (10,000 permutations). BH-corrected q-values. CHM13 is haploid (singleton communities marked "—"). RPE-1 is female XX (C14 chrX_q+chrY_q has one haplotype only).

**Result.** Per-community enrichment ratio (bold = q<0.001, italic = q<0.05, plain = not significant, — = singleton/absent):

| Community | Arms | CHM13 | HG002 | HG02559 | HG00658 | HG02148 | NA19036 | HG002 PoreC | HG002 CiFi |
|---|---|---|---|---|---|---|---|---|---|
| C1 | chr4_q, chr10_q | *44.5x* | *18.9x* | *45.5x* | *13.3x* | *23.6x* | *43.5x* | *7.9x* | *13.3x* |
| C2 | chr10_p, chr18_p | *15.9x* | **30.7x** | *19.1x* | *17.1x* | *23.8x* | **35.3x** | *8.6x* | *14.6x* |
| C3 | chr3_q, chr7_p, chr9_q, chr11_p, chr16_q, chr19_p | **7.5x** | **5.4x** | 1.6x | **7.7x** | *5.6x* | *4.8x* | *3.6x* | *7.9x* |
| C4 | chr7_q, chr12_q | 0.2x | *12.7x* | *30.5x* | *18.9x* | *10.6x* | *12.6x* | *7.6x* | *19.7x* |
| C5 | chr6_p, chr9_p, chr12_p, chr20_q | *8.3x* | **11.9x** | **11.1x** | **14.0x** | *7.4x* | *10.8x* | *2.9x* | 2.6x |
| C6 | chr1_q, chr13_q, chr17_q, chr19_q, chr21_q, chr22_q | *3.7x* | **5.3x** | 1.8x | 0.0x | 0.0x | 2.3x | **13.1x** | **11.6x** |
| C7 | chr13_p, chr14_p, chr15_p, chr21_p, chr22_p | **7.6x** | **7.2x** | *4.0x* | — | 0.0x | 0.0x | *4.6x* | *3.7x* |
| C8 | chr15_q | — | *17.0x* | 0.0x | 0.0x | — | 0.0x | *5.9x* | 2.7x |
| C9 | chr16_p | — | *10.2x* | 0.0x | 0.0x | *53.8x* | *20.5x* | 2.5x | *6.5x* |
| C10 | chr17_p | — | *63.9x* | *56.1x* | 0.0x | 0.0x | — | *47.5x* | 0.0x |
| C11 | chr1_p, chr5_q, chr6_q, chr8_p | *7.3x* | **9.2x** | *11.1x* | **9.0x** | 1.6x | *5.6x* | *5.7x* | *7.3x* |
| C12 | chr2_q, chr20_p | 2.4x | *18.1x* | *13.0x* | *12.8x* | *32.5x* | *11.2x* | *11.6x* | *22.0x* |
| C13 | chr4_p | — | *76.5x* | *59.1x* | 0.0x | 0.0x | 0.0x | *54.7x* | 0.0x |
| C14 | chrX_q, chrY_q | — | *75.3x* | 0.0x | 0.0x | *136.6x* | *92.9x* | *9.4x* | *12.3x* |
| C15 | chrX_p, chrY_p | — | *71.3x* | 0.0x | 0.0x | *165.3x* | *90.1x* | *19.9x* | *18.3x* |

**Key observations.**

- **C1** (D4Z4, chr4_q + chr10_q): Most reproducible community. Significant in all 8 datasets. Enrichment 7.9–45.5x.
- **C2** (chr10_p + chr18_p): Significant in all 8 datasets. Enrichment 8.6–35.3x.
- **C4** (chr7_q + chr12_q): Significant in 7/8 datasets (not CHM13). Largest ratio in HG02559 (30.5x).
- **C7** (acrocentric p-arms): Significant in 5/8 datasets. Driven by nucleolar co-localization of NOR-bearing short arms.
- **C14/C15** (PAR2/PAR1): Extremely high ratios in male samples (HG002 71–75x, HG02148 137–165x, NA19036 90–93x). Zero in female or haploid samples, as expected.
- **C13** (chr4_p): Only significant in high-depth datasets (HG002, RPE-1). Ratios are large (32–106x) when detectable, suggesting a real but depth-limited signal.
- **C10** (chr17_p): Singleton in CHM13 and NA19036. Significant in all datasets where testable (5/5 non-singleton), with ratios 20–65x.
- **C14** (chrX_q + chrY_q) and **C15** (chrX_p + chrY_p + chr18_q): Sex-linked. C14 significant only in male diploid samples with Y chromosome. C15 includes chr18_q, reaching significance in RPE-1 (female XX) where only the chrX_p + chr18_q pairing contributes.
- **C1** (chr4_q + chr10_q, D4Z4): Borderline in most datasets. Only significant in 3/10 (HG02559, NA19036, HG002 Pore-C). The D4Z4 repeat may reduce mappability at 50 kb resolution.
- **CHM13**: No community reaches significance. Haploid reference limits to 37 arms, most communities are singletons or 2-arm, reducing statistical power.
- **RPE-1 vs LCL**: RPE-1 ratios are systematically higher than LCL for the same communities (e.g., C3: 12x in RPE-1 vs 4–5x in HG002; C6: 7–8x vs 5x). This may reflect the CiFi/Pore-C technology capturing multi-way contacts more efficiently than short-read Hi-C.

**Conclusion.** The per-community view confirms that the 3D community signal is not driven by a single dominant community. At least 8 of 15 communities (C2–C7, C9–C13) are independently significant in multiple datasets across both cell types. The remaining communities are limited by sex (C14, C15), sequencing depth (C1, C13), or ploidy (CHM13 singletons).

---

## 13. 3D genome validation: Dip-C single-cell (T2T-CHM13v2.0)

**What it does.** Tests community structure using single-cell 3D genome structures, providing a complementary approach to Hi-C. While Hi-C measures contact frequencies averaged over millions of cells, Dip-C reconstructs the physical 3D coordinates of each genomic locus in individual cells. This allows testing whether community-member arms are physically closer in 3D nuclear space, not just more frequently in contact.

**How.** 17 GM12878 cells remapped to T2T-CHM13v2.0 using BWA-MEM2, hickit for 3D modeling, dip-c impute3 for diploid haplotype refinement (4 rounds). MAPQ=0 for maximum subtelomeric coverage. SNPs from 1KGP CHM13v2 panel (NA12878). **16 cells** used (cell 12 excluded as duplicate of cell 10). 3D particle positions selected using per-arm PHR coordinates for the 38 C-community arms (from CHM13 PHR boundaries, 10–500 kb arm-specific) and 500 kb terminal regions for the 8 arms without CHM13 PHR (7 S-community arms + chr6_p).

### 13.1 Community 3D enrichment (T2T)

**Result.** Community-based results (16 cells, per-arm PHR coordinates):

| Metric | Value |
|--------|-------|
| Wilcoxon signed-rank | stat=8.0, p=3.8e-04 |
| Fisher combined | chi2=75.3, p=2.4e-05 |
| W/B ratio (mean) | 0.931 (6.9% closer within-community) |
| W/B ratio (median) | 0.934 |
| Mantel rho | 0.296, p=0.002 |

### 13.1b Supplementary communities: non-sharing arms

**What it does.** Tests whether the 7 chromosome arms with zero inter-chromosomal subtelomeric sequence sharing (chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q) show 3D proximity patterns comparable to the sequence-sharing communities.

**How.** Each of the 7 arms is assigned a singleton supplementary community (S1–S7). For the pooled S_all test, all 7 are treated as one pseudo-community. The analysis uses the terminal 500 kb of each arm (since they have no PHR). The same W/B within/between framework is applied: if S_all W/B < 1, non-sharing arms are closer to each other than to random arms; if W/B > 1, they are farther apart.

**Result.**

| | GM12878 (16 cells) | Sperm (20 cells) |
|-|---:|---:|
| **C-community W/B (mean)** | 0.931 (6.9% closer) | 0.401 (60% closer) |
| **S_all W/B (mean)** | 1.106 (11% farther) | 1.397 (40% farther) |
| **S_all cells < 1.0** | 0/16 | 1/20 |

Per-S singleton radial positions:

| Arm | S-label | GM12878 radial | Sperm radial | Nearest C-community in 3D |
|-----|---------|---------------|-------------|--------------------------|
| chr2_p | S1 | 0.81 (peripheral) | 0.68 | C12 (chr2_q, same chrom) |
| chr3_p | S2 | 0.84 (peripheral) | 0.82 | C3 (chr3_q, same chrom) |
| chr5_p | S3 | 0.66 | 0.61 | C3 (chr19_p) |
| chr8_q | S4 | 0.60 | 0.58 | C3 (chr19_p) |
| chr11_q | S5 | 0.69 | 0.65 | C6 (chr22_q) |
| chr14_q | S6 | 0.61 | 0.56 | C3 (chr19_p) |
| chr18_q | S7 | 0.78 (peripheral) | 0.71 | C2 (chr18_p, same chrom) |

**Conclusion.** Non-sharing arms are consistently **farther apart** in 3D space (GM12878: 11% farther, sperm: 40% farther), the opposite of the sequence-sharing C-communities (6.9% and 60% closer, respectively). This provides a negative control: subtelomeric sequence sharing is necessary for 3D proximity clustering. The per-arm radial analysis shows that most non-sharing arms' nearest 3D neighbor is the opposite arm of the same chromosome (cis-arm proximity), not any inter-chromosomal community partner. S1/chr2_p, S2/chr3_p, and S7/chr18_q are notably peripheral (radial > 0.78), consistent with their telomere-proximal nuclear positioning.

### 13.2 Pangenome-level community-free (3D distance)

**What it does.** Tests the similarity-3D distance relationship without community labels, mirroring the mcool community-free analysis (§12.5b) but using 3D Euclidean distance instead of Hi-C contact.

**How.** For each of 465 pangenome sample#hap combinations, take that sample's PHR sequences, map to 3D positions in the Dip-C cell, compute inter-chromosomal (Jaccard, 3D distance) pairs, Spearman correlation. This mirrors the mcool community-free exactly but uses 3D distance instead of Hi-C contact.

**Result.** Per-cell results:

| Cell | rho |
|------|-----|
| 01 | -0.160 |
| 02 | -0.067 |
| 03 | -0.060 |
| 05 | -0.004 |
| 06 | -0.167 |

Arm-level: rho = -0.336, p = 1.1e-18, n=652 arm pairs (negative = more similar arms are physically closer in 3D space).

**Key metrics.** Runtime: 5 minutes (was 1.5h+ with old O(n^2) all-haplotypes approach).

**Conclusion.** The community-free Dip-C analysis confirms that subtelomeric sequence similarity predicts 3D proximity at both per-cell and arm-level scales, without relying on any discrete community assignment. The negative rho direction is consistent with the Hi-C community-free results (§12.5b): more similar sequences are closer in 3D space. The arm-level result (rho=-0.336, p=1.1e-18) is highly significant.

---

## 13b. 3D genome validation: sperm single-cell

**What it does.** Tests community structure in haploid sperm cells, providing a complementary single-cell 3D validation in a distinct cell type with unique nuclear architecture.

**How.** 20 sperm cells (10 X-bearing + 10 Y-bearing) from Xu et al. 2025. Haploid: no impute3 needed. Same pangenome-level community-free approach as Dip-C (§13.2).

### 13b.1 Community-based results

**Result.**

| Metric | Value |
|--------|-------|
| W/B ratio | 0.401 (60% closer within-community) |
| Fisher combined p | 3.9e-51 |
| Mantel rho | 0.202, p=0.023 (significant) |

**Conclusion.** Community-based enrichment is very strong with PHR-specific coordinates (W/B=0.401, 60% closer, Fisher p=3.9e-51). The Mantel test is now significant (rho=0.202, p=0.023), indicating that the sequence-distance vs 3D-distance relationship holds even in the highly condensed sperm nucleus. The dramatic improvement from 35% to 60% closer reflects the precision gain from using exact PHR boundaries rather than default 300kb windows.

### 13b.2 Community-free results

**Result.** Per-cell rho median = -0.029, with 15/20 cells showing negative rho (more similar arms are closer in 3D space). Arm-level: rho = 0.048, p = 0.197 (ns).

**Conclusion.** The community-free analysis confirms the trend: 75% of sperm cells show the expected negative correlation between sequence similarity and 3D distance. The weaker per-cell effect (compared to GM12878 Dip-C) is consistent with the highly compacted sperm chromatin architecture limiting inter-chromosomal proximity variation.

---

## 14. Integrated 3D interpretation

### 14.1 Convergent evidence across technologies

**What it does.** Summarizes 3D genome evidence across three independent technologies, each with distinct biases and resolution.

**How.** Hi-C captures pairwise ligation contacts in bulk cell populations (millions of cells). Pore-C captures multi-way contacts via long Nanopore reads spanning multiple ligation junctions. Dip-C reconstructs 3D chromatin positions in individual cells at 20 kb resolution. RPE-1 CiFi/Pore-C provides cell-type and cell-cycle controls.

**Result.**

| Technology | Samples | Test | Effect size | p-value |
|---|---|---|---|---|
| Hi-C | 5 diploid HPRC + CHM13 | B/W ratio | 0.027–0.074 | 6.0e-18 to 9.4e-03 |
| Pore-C | HG002 | B/W ratio | 0.056 | 3.9e-85 |
| CiFi | HG002 | B/W ratio | 0.036 | 2.0e-74 |
| Dip-C (T2T) | GM12878 (16 cells) | Within/between 3D distance | 6.9% closer | 2.4e-05 (Fisher) |
| Dip-C (T2T) | GM12878 (16 cells) | S_all (non-sharing arms) | 11% farther | negative control |
| Dip-C (T2T) | GM12878 (16 cells) | Mantel: Jaccard vs 3D distance | rho=0.296 | 0.002 |
| Dip-C (T2T) | GM12878 (16 cells) | Community-free arm-level | Spearman rho=−0.336 | 1.1e-18 |
| Sperm 3D | 20 cells (Xu et al. 2025) | Within/between 3D distance | 60% closer | 3.9e-51 (Fisher) |
| Sperm 3D | 20 cells | Community-free per-cell | 15/20 negative rho | — |
| RPE-1 CiFi | Async (interphase) | B/W ratio (50kb / 10kb) | 0.017 / 0.017 | 1.1e-105 / 7.4e-106 |
| RPE-1 Pore-C | Async (interphase) | B/W ratio (50kb / 10kb) | 0.018 / 0.020 | 1.0e-116 / 8.6e-119 |
| RPE-1 CiFi | Mitotic | B/W ratio (50kb / 10kb) | 0.025 / 0.032 | 1.6e-65 / 5.6e-67 |
| Mouse 1Mb | 4 meiotic stages | B/W ratio range (5 resolutions) | 0.040–0.112 | all < 1e-22 |
| Mouse 4Mb | 4 meiotic stages | B/W ratio range (5 resolutions) | 0.176–0.213 | all < 1e-36 |

Additionally, per-PHR-pair correlation confirms the continuous relationship: arms with more similar subtelomeric sequence show proportionally more inter-chromosomal contact (CHM13 Hi-C: Spearman rho = 0.674, p = 3.0e-92; HG002 Pore-C: rho = 0.485, p = 1.6e-48; all 8 datasets significant). Dip-C community-based Mantel confirms this in 3D coordinates (rho = 0.296, p = 0.002), and the community-free arm-level analysis shows Spearman rho = −0.336, p = 1.1e-18.

*Note: All analyses in this table — human Hi-C/Pore-C/CiFi, RPE-1, and mouse — are now complete at all 5 mcool resolutions (5kb, 10kb, 20kb, 50kb, 100kb). Full multi-resolution tables are in §12.7 (human), §12.9 (no-acrocentric), §12.14.1 (RPE-1), and §15.7b (mouse). Results are consistent across all resolutions, with slight strengthening at finer resolutions reflecting better bin coverage of the median 105 kb PHR regions.*

**Conclusion.** All technologies — using different cell types (LCL, RPE-1, sperm), different platforms (Illumina, PacBio, ONT), and different measurement principles (ligation frequency, multi-way contacts, 3D coordinates) — independently confirm community-structured 3D organization. The addition of Dip-C T2T remapping and sperm single-cell 3D data extends validation to haploid cells and eliminates hg19/T2T coordinate incompatibility concerns (§16, limitation 12). The signal is not an artifact of any single technology, cell type, or analysis method. Multi-resolution analysis at 5 resolutions (5kb-100kb) confirms resolution-invariance across all systems.

### 14.2 Flanking region paradox

**What it does.** Explains the counterintuitive finding that unique-sequence regions flanking PHR boundaries show stronger 3D signal than the duplicated PHR regions themselves.

**How.** PHR intervals contain inter-chromosomal duplications (§4). When Hi-C reads map to these duplicated regions, multi-mapping between community partner arms creates ambiguity — a read aligning to a duplicated block could originate from either partner, and multi-mapper removal (RM_MULTI=1 in HiC-Pro) discards these reads. The 100 kb flanking regions centromere-ward of the PHR boundary are unique sequence with no multi-mapping.

**Key metrics.** Flanking B/W ratios: CHM13 0.057, HG002 0.002, HG02559 0.008, HG00658 0.006, HG02148 0.005, NA19036 0.006, HG002 Pore-C 0.034. All stronger than PHR (§12.1: 0.027–0.074), consistent with multi-mapping suppression at duplicated PHR intervals. HG002 flanking enrichment (B/W=0.002) is 13x stronger than PHR (B/W=0.027). In Dip-C, flanking particles are more interior than non-flanking terminal particles (GM12878: 0.634 vs 0.658, p = 3.8e-9).

**Conclusion.** Three implications: (1) Rules out multi-mapping as the driver of the 3D signal — flanking regions have no duplicated content, yet show stronger clustering. (2) 3D clustering extends beyond the duplicated region into flanking unique sequence, indicating broader chromosomal domain effects. (3) Multi-mapping suppresses the PHR signal; flanking enrichment values represent the true signal magnitude.

### 14.3 Meiotic bouquet as exchange venue

**What it does.** Considers the meiotic bouquet as the context for subtelomeric exchange events.

**Key metrics.** All 3D data is somatic (interphase). Tan et al. (2018) found Rabl configuration "weak" in GM12878 and PBMCs. Zuo et al. (2021) showed chromosome end alignment extends "a substantial range of ~20% of chromosome length" in mouse meiosis. Average meiotic loop sizes: ~500 kb at leptotene, ~700 kb at zygotene. Median PHR region (105 kb) fits within a single meiotic loop.

The hypothesis that meiotic telomere clustering favours subtelomeric exchange has a long pedigree: Mefford & Trask (2002) noted that "the pairing of homologues typically begins at the telomeres" and that subtelomeric homology creates opportunities for ectopic exchange, while Linardopoulou et al. (2005) stated explicitly: "Telomere clustering in meiotic cells might favour exchange of chromosome ends during DSB healing." Mefford & Trask also noted that the interphase clustering of chr4q and chr10q (Stout et al. 1999) "could promote their exchange" — a prediction now supported by Hi-C/Pore-C data showing community-structured 3D co-localization across all 15 communities. In mouse meiosis, Patel et al. (2019) first detected X-shaped interchromosomal Hi-C contacts consistent with the meiotic bouquet, and Zuo et al. (2021) extended this with stage-resolved analysis showing that chromosome end alignment during early prophase (leptotene and zygotene) extends over "a substantial range of ~20% of chromosome length" — far deeper than the subtelomeric regions analyzed here. Notably, this alignment is not merely a manifestation of the transient bouquet conformation (<5% of zygotene cells show a cytological bouquet); rather, it occurs whenever different chromosomes are brought into proximity by telomeres, is independent of compartment identity, and depends on the force-transmitting LINC complex (which modulates alignment range but does not alter loop sizes). Average meiotic chromatin loop sizes increase from ~500 kb at leptotene to ~700 kb at zygotene, reaching ~1.4 Mb at pachytene and ~1.6 Mb at diplotene (Zuo et al. 2021; Patel et al. 2019 reported slightly larger values of 0.8–1.0 Mb at zygotene). Since the median PHR region length is 105 kb, most PHR sequences would reside at the base of a single meiotic loop at leptotene — the compartment where recombination machinery is concentrated — maximizing the probability of ectopic recombination between aligned PHR sequences on different chromosomes.

If human meiosis follows a similar pattern, the community-structured 3D contacts observed in somatic cells may underestimate the true meiotic proximity. The observation that flanking regions (100 kb centromere-ward) show stronger 3D signal than PHRs themselves (§14.2) is consistent with the meiotic proximity extending beyond the duplicated region.

Human meiotic Hi-C remains the single most informative missing experiment. Existing data from mouse spermatocytes (Patel et al. 2019) cannot be directly extrapolated to human subtelomeric organization.

### 14.4 D4Z4-CTCF-Lamin tethering model for C1

**What it does.** Proposes a specific molecular mechanism for C1 (chr4_q/chr10_q) co-localization, supported by existing literature.

**Key metrics.** C1: silhouette = 0.149, 42.8% discordance. Dip-C radial = 0.744 (peripheral). Inter-chromosomal signal peaks at 0–15 kb (D4Z4 position). C1 sequences: median 22 DUX4L; non-C1 outliers: 0–2 (Mann-Whitney p = 5.3e-6).

**How.** The mechanism:

1. Both chr4_q and chr10_q carry D4Z4 macrosatellite arrays at their subtelomeric tips
2. CTCF binds within D4Z4 repeat units (Ottaviani et al. 2009)
3. D4Z4-proximal sequences are tethered to the nuclear periphery via lamin A/C interaction (Masny et al. 2004; Ottaviani et al. 2009)
4. Both arms are positioned at the nuclear periphery via lamin A/C interaction, and this co-localization is consistent with the elevated ectopic recombination observed clinically

The Dip-C radial analysis supports this: C1 arms occupy peripheral nuclear positions (among the most exterior communities). The pangenome data provides evidence consistent with D4Z4 contributing to the chr4_q/chr10_q co-clustering (§6.2): inter-chromosomal signal peaks at 0–15 kb from the telomere (where D4Z4 sits), C1 sequences carry median 22 DUX4L genes while all 7 non-C1 outlier sequences have 0–2 on their own arm (Mann-Whitney p = 5.3e-6), and outlier PHR regions are 4.6–9x shorter than C1 sequences. As established in the FSHD literature (Lemmers et al. 2010), ectopic exchange between chr4_q and chr10_q D4Z4 arrays can modify FSHD alleles — a D4Z4 contraction on a permissive 4qA haplotype that translocates to chr10_q loses pathogenicity, while the reverse can create disease alleles.

This mechanism predicts that CTCF/cohesin density at PHR boundaries should correlate with Hi-C contact strength between community partners. This is testable using Gershman et al.'s (2022, Science) ENCODE CTCF ChIP-seq realignment to T2T-CHM13, which found CTCF enrichment at TAR loci across all ENCODE cell lines, and Stergachis lab Fiber-seq data providing single-molecule CTCF maps at 39/46 telomeres. Standard hg38-aligned ENCODE data would be inadequate because hg38 has incomplete subtelomeric assemblies.

### 14.5 Nucleolar association mechanism for C6/C7

**What it does.** Considers nucleolar co-localization as a mechanism for C7 (acrocentric p-arms) homogenization and C6 (acrocentric q-arms) community membership.

**Key metrics.** C7: silhouette = −0.027, gene replacement scores 0.91–1.0 for chr13_p/chr14_p/chr15_p, 0.49–0.54 for chr21_p/chr22_p. C6: silhouette = 0.505, Dip-C radial = 0.574 (interior, consistent with nucleolar positioning). C7 cannot be assessed in Dip-C (hg19 p-arms unmapped).

**Result.** All five acrocentric short arms carry rDNA and constitutively associate with nucleoli. The near-complete interchangeability (silhouette = −0.027) is consistent with frequent exchange at nucleolus-associated arms. C6 includes non-acrocentric arms (chr1_q, chr17_q, chr19_q), indicating q-arm membership is driven by shared duplicon content, not nucleolar proximity alone.

**Conclusion.** Nucleolar co-localization is consistent with the literature but the present study provides no C7-specific 3D data. The f7501 duplicon distribution (Mefford & Trask 2002) maps directly onto C3 (chr3_q, chr7_p, chr9_q, chr11_p, chr16_q, chr19_p), confirming that community structure captures known duplicon module relationships.

### 14.6 Causal feedback loop

**What it does.** Proposes a feedback loop: sequence similarity → 3D proximity → ectopic exchange → increased similarity.

**How.** Four links with varying levels of support:

1. **Sequence similarity → 3D proximity**: Mantel tests show continuous correlation between Jaccard distance and Hi-C contact (HG002 rho=−0.41, Pore-C rho=−0.50)
2. **3D proximity → ectopic exchange**: Established from FSHD literature (chr4_q/chr10_q D4Z4 translocation) and general principles of recombination requiring physical proximity
3. **Ectopic exchange → increased similarity**: Inferred from the outcome — cross-arm affinity analysis shows 15.1% of sequences resemble a partner arm more than their own, consistent with past exchange having increased inter-arm similarity

The fourth link — **increased similarity → stronger future proximity** — is inferred but not directly measured. Testing this would require temporal data (e.g., comparing 3D contacts before and after an exchange event) or comparing Hi-C contact strength between haplotypes that carry cross-arm vs self-arm subtelomeric types within the same individual.

This circularity means the causal direction cannot be established from the present data: communities could form because similar sequences are brought into proximity, or sequences could become similar because they are in proximity. The meiotic bouquet (§14.3) provides a speculative but plausible initiation scenario — telomere clustering during meiosis brings all chromosome ends into proximity regardless of sequence content, and subsequent ectopic exchange between neighboring arms would create the initial sequence similarity that then reinforces itself through the feedback loop. This bouquet initiation model cannot be tested with the present data and is offered as a plausible scenario, not an evidence-based conclusion. This feedback concept is not new: Linardopoulou et al. (2005) described a cycle of "segmental polymorphism and gross genomic rearrangement" (their Fig. 2) where translocations create duplications that promote further rearrangements, and Mefford & Trask (2002) noted that interphase clustering of chr4q/chr10q "could promote their exchange." The present analysis adds the 3D proximity dimension with quantitative support. Additionally, Ambrosini et al. (2007) proposed that the 98% identity peak in their bimodal duplicon distribution "may be due to maintenance of sequence similarity by ongoing interchromosomal gene conversion between the large subtelomeric duplicons" — a hypothesis consistent with the ongoing-exchange link in this model.

### 14.7 Testable predictions

Three predictions arising from the meiotic chromosome organization data of Zuo et al. (2021) that could not be tested with the present data:

1. **LINC complex requirement**: The LINC complex (Linker of Nucleoskeleton and Cytoskeleton) spans the nuclear envelope and transmits cytoskeletal forces to chromosomes, driving telomere-led movements during meiotic prophase. Zuo et al. showed that a SUN1 point mutation (W151R) disrupts long-range chromosome-end alignment in mouse meiosis: while tip contacts (within ~5% of chromosome length) actually increase, alignment at greater distances drops off sharply, reducing the effective alignment zone from ~20% to ~5% of chromosome length. If meiotic Hi-C from SUN1 mutant spermatocytes were analyzed with the community framework, within-community inter-chromosomal contacts at chromosome ends should be dramatically reduced compared to wild-type. This would establish that LINC-mediated force transmission is required for community-structured 3D contacts and, by extension, for the ectopic exchange that maintains subtelomeric homology.

2. **Crossover frequency correlation** (tested, confounded): Using the T2T-CHM13 recombination map (Lalli et al. 2025, preprint), subtelomeric recombination rate anticorrelates with cross-arm affinity across all 36 arms (rho = −0.43, p = 0.008). However, this signal is entirely driven by 7 confounded arms: acrocentric p-arms (chr13_p, chr14_p, chr15_p, chr21_p, chr22_p) and PAR (chrX_p, chrX_q) have 0–12 callable variants in 500 kb (vs 1,000–3,000 for non-acrocentric arms), so their 0 cM/Mb recombination rate reflects absence of short-read genotyping data in repetitive regions, not necessarily absence of recombination. Excluding these 7 arms: rho = 0.03, p = 0.86, N = 29 — the correlation vanishes. The question of whether local recombination protects arm identity from ectopic exchange remains biologically plausible but cannot be tested with current recombination map data at subtelomeric loci. Long-read-based recombination maps resolving variants in repetitive subtelomeric regions would be needed (see §17.3 prediction 7 for the full per-arm table).

3. **Chromatin compartment identity**: Zuo et al. showed that A-compartment (transcriptionally active) regions form shorter meiotic loops (~560 kb at leptotene) with higher crossover rates, while B-compartment regions form longer loops (~730 kb). Compartment calling from the HG002 Hi-C eigenvector (100 kb resolution, per-haplotype) shows that 63% of chromosome tips are classified as B-compartment (58/92 arm × haplotype combinations), compared to 56% for mid-chromosome interior bins (p = 0.003, Mann-Whitney). However, the eigenvector values at chromosome tips are close to zero (mean −0.0014 vs +0.0011 for interior), indicating that subtelomeric regions have weak, poorly defined compartment identity — consistent with their nature as transitional zones between chromosome-specific and duplicated sequence. The Dip-C radial data (§12.2) adds a spatial dimension: PHR-flanking regions are more interior (radial 0.634) than other terminal regions (0.658, p = 3.8e-9), and the communities with the most inter-chromosomal sharing (C6, C9, C8, C10) are the most interior (radial 0.56–0.57). Subtelomeric regions are thus pseudogene-rich and weakly B-compartment by chromatin state, yet positioned internally rather than at the lamina. This dissociation suggests that telomere clustering, not lamina association, determines the nuclear positioning of these regions. Per Zuo et al., weak B-compartment identity genome-wide predicts longer meiotic loops (~730 vs ~560 kb). However, Zuo explicitly found that B-compartment regions at centromere-distal chromosome ends have shorter loops than B-compartment elsewhere, and stated that loop size differences are "unlikely the main cause for prominent alignment among chromosome ends" — the end alignment pattern is independent of compartment identity. This favors bouquet-stage telomere clustering rather than loop-axis accessibility as the primary driver of inter-chromosomal exchange at subtelomeres.

---

The mechanistic models above propose specific feedback loops between sequence similarity and 3D proximity. The following section situates these findings within the broader literature and identifies novel contributions.

## 15. Cross-species validation: mouse T2T subtelomeric analysis

### 15.1 Motivation

**What it does.** Tests whether the subtelomeric sequence homology → 3D proximity link is conserved across mammals.

**How.** The same pipeline (§1–§6) applied to two mouse T2T assemblies: C57BL/6J (B6, GCA_964188535) and CAST/EiJ (GCA_964188545) from Francis et al. (2025). Mouse chromosomes are telocentric — centromere at one end, telomere at the other — giving one subtelomeric region per chromosome.

### 15.2 Flank extraction

**What it does.** Extracts 500 kb flanks from both chromosome ends.

**Key metrics.** B6 19 autosomes, CAST 19 autosomes + chrX. 78 total flanks. Telomeric tracts: B6 1,695–26,101 bp, CAST 3,329–21,637 bp.

**Result.** 78 flanks of 474–498 kb.

### 15.3 Inter-chromosomal region detection

**What it does.** Identifies inter-chromosomal signal in mouse subtelomeric flanks.

**How.** wfmash `-p 95` (do_not_overfilter branch, cc60cd8 — k-mer frequency filtering no longer needs manual `-F 0.1` override), impg v0.4.0, same detection script with `--min-count 1`. Total alignments: 49,911.

**Result.** 39/78 flanks have signal:
- All 39 `_parm` flanks (centromere-distal/subtelomeric end) have signal.
- All 39 `_qarm` flanks (centromere-proximal end) have zero signal — satellite-dominated, confirming centromere coordinates are unnecessary.

Most subtelomeric flanks show inter-chromosomal signal spanning the **entire 500 kb window** (440–495 kb regions matching all 20 chromosomes). This is much larger than the TLC repeat itself (6–12 kb in B6 per Francis et al. 2025), indicating that the shared subtelomeric content extends well beyond the TLC motif into flanking repeat-rich sequence.

**Notable exceptions in B6:**
- B6 chr11 and chr18 form a private group with chr4 (signal only with chr4/chr11/chr18; chr11: 80 kb, chr18: 30 kb). Francis et al. (2025) noted that B6 chr11 is an exception to the standard TLC end structure — it lacks the conserved L1-LINE → TLC → minor satellite motif.
- B6 chr7 forms a small private group (chr7/chrX/chr11, 10 kb). Francis et al. also identified chr7 as a structural exception.
- B6 chr4: 15 kb signal (much smaller than the typical ~490 kb), matching 18 chromosomes (missing chr7 and chrX).

### 15.4 Pangenome graph and similarity

**What it does.** Computes pairwise similarity among mouse PHR sequences.

**How.** PGGB (`-p 95 -n 2`) on 39 PHR sequences; Jaccard similarity via `odgi similarity --all`.

**Result.** Key pairwise similarities:

| Pair | Jaccard | Region sizes |
|------|---------|-------------|
| CAST chr1 ↔ CAST chr2 | 0.980 | 490 kb ↔ 480 kb |
| B6 chr11 ↔ B6 chr18 | 0.371 | 80 kb ↔ 30 kb |
| CAST chr13 ↔ CAST chr14 | 0.030 | 485 kb ↔ 490 kb |
| CAST chr7 ↔ CAST chr16 | 0.028 | 475 kb ↔ 495 kb |
| Most other cross-chr pairs | 0.025–0.030 | ~490 kb each |
| All B6 cross-chr pairs (except chr11/chr18) | ~0 | variable |

**Conclusion.** CAST chr1/chr2 (Jaccard = 0.980) shows nearly complete subtelomeric identity across ~490 kb — analogous to human C7 (acrocentric p-arms). B6 shows almost no inter-chromosomal signal beyond the chr11/chr18 pair, while CAST shows low but consistent cross-chromosome similarity (0.025–0.030). This contrast is consistent with Francis et al.'s description: B6 has uniform L1-LINE + TLC architecture; CAST has heterogeneous subtelomeric repeat organization.

### 15.5 Open questions

1. **Repeat annotation**: RepeatMasker and GFF3 annotations are available for both T2T assemblies from Ensembl (https://projects.ensembl.org/mouse_genomes/). These should be downloaded and intersected with the PHR regions to determine what repeat content drives the inter-chromosomal signal — the large shared regions (440–495 kb) extend far beyond the TLC repeat (6–12 kb) and may include L1-LINE, LTR, and other repeat classes.
2. **Mouse-human synteny** (from mm39→hg38 syntenic net, UCSC): The distal (telomeric) ends of the key mouse chromosomes map to the following human regions:

| Mouse chr (distal end) | Human syntenic region | Human subtelomeric community |
|---|---|---|
| CAST chr1 (pair, J=0.980) | chr8:51.8–55.6 Mb | interior, not subtelomeric |
| CAST chr2 (pair, J=0.980) | chr10:5.9–15.4 Mb | interior, not subtelomeric |
| B6 chr11 (pair, J=0.371) | chr22:28.8–31.6 Mb | interior, not subtelomeric |
| B6 chr18 (pair, J=0.371) | chr10:35.0–35.2 Mb | interior, not subtelomeric |

The mouse private pairs do NOT correspond to human subtelomeric regions — the syntenic human positions are in the interior of human chromosomes, far from any telomere. This is expected: mouse chromosomes are telocentric (centromere at one end) while human chromosomes are metacentric, so the mouse distal telomeric end often maps to human chromosome interiors. The mouse subtelomeric inter-chromosomal sharing is driven by repeat architecture (TLC, L1-LINE), not by syntenic conservation of a subtelomeric duplicon system.
3. **Community detection**: Leiden community detection (`detect_communities.R`, `--organism mouse --level arm`, resolution scan 0.1–3.0) on the arm-level Jaccard distance matrix finds **24 communities** from 39 arms: 10 multi-arm communities and 14 singletons.

| Community | Members | Notes |
|---|---|---|
| C6 | CAST chr12, chr13, chr16, chr19, B6 chr15 | Large CAST bulk + 1 B6 |
| C7 | CAST chr11, chr5, B6 chr16 | Cross-strain |
| C20 | CAST chr10, chr4, chr6 | CAST triplet |
| C23 | CAST chr1, chr2 | Previously J=0.980 |
| C1 | B6 chr10, CAST chr9 | Cross-strain pair |
| C13 | B6 chr3, CAST chr7 | Cross-strain pair |
| C16 | B6 chr6, CAST chr8 | Cross-strain pair |
| C19 | B6 chr9, CAST chrX | Cross-strain pair |
| C3 | B6 chr12, CAST chr17 | Cross-strain pair |
| C12 | B6 chr2, CAST chr15 | Cross-strain pair |
| 14 singletons | B6 chr1,4,5,7,8,11,13,14,17,18,19; CAST chr3,14,18 | All isolated |

**Key change from prior run** (which used `-F 0.1`): updated wfmash (`do_not_overfilter` branch) reduces alignments from 61,824 to 49,911 (19% fewer, better-filtered) and produces more structured communities — 10 multi-arm communities (vs 3 previously) with **6 cross-strain pairs** (B6+CAST arms sharing a community), indicating that inter-chromosomal subtelomeric sharing extends across mouse strains, not just within CAST. The CAST chr1/chr2 pair (J=0.980) persists as C23. The old B6 chr11+chr18 pair (J=0.371) is now split into singletons, suggesting its signal was partly driven by the over-filtering artifact.

This is fundamentally different from human (15 multi-arm communities with complex overlapping structure across 41 arms). Mouse has more singletons (14/24 = 58%) reflecting the predominantly uniform TLC-based architecture.

### 15.6 Pipeline update: 500kb → 1Mb windows

**What it does.** The 500kb subtelomeric pipeline was deleted because PHRs filled the entire 500kb window, causing flanking regions to be truncated on 30/49 arms. A 1Mb pipeline was completed instead.

**Key metrics.** 1Mb subtelomeric windows (`subtelo_1Mb` directory): 49 PHR regions detected. Flanking extraction: 0 truncated (the larger window provides sufficient room for both PHR and 100kb flanking regions).

**Result.** Community-free with Hi-C (1Mb, 50kb resolution):

| Meiotic stage | rho (all pairs) | p-value | rho (nonzero Jaccard) | rho (nonzero contact) |
|--------------|----------------|---------|----------------------|----------------------|
| Leptotene | 0.259 | 2.7e-18 | 0.107 | 0.220 |
| Zygotene | **0.425** | 4.2e-51 | 0.329 | 0.393 |
| Pachytene | **0.428** | 9.7e-52 | 0.135 | 0.309 |
| Diplotene | **0.416** | 4.9e-43 | 0.125 | 0.212 |

All four stages show significant positive correlation between sequence similarity and Hi-C contact. This is a major improvement over the old 500kb results (rho=0.08–0.17, all non-significant), confirming that the 1Mb windows capture the full PHR extent. The zygotene/pachytene values (rho=0.425–0.428) coincide with the meiotic bouquet stage, when telomere-led chromosome movements bring subtelomeric regions into proximity.

Community-based per-arm-pair correlation (1Mb, 50kb resolution):

| Meiotic stage | rho | p-value | Mantel rho | Mantel p |
|--------------|-----|---------|------------|----------|
| Leptotene | 0.834 | 2.8e-90 | -0.456 | <0.0001 |
| Zygotene | **0.862** | 8.1e-103 | -0.492 | <0.0001 |
| Pachytene | **0.845** | 7.4e-95 | -0.473 | <0.0001 |
| Diplotene | 0.633 | 5.5e-40 | -0.390 | 0.0001 |

Flanking region control (100kb unique sequence centromere-ward of PHR, 50kb resolution):

| Meiotic stage | rho | p-value | nonzero Jaccard pairs |
|--------------|-----|---------|----------------------|
| Leptotene | 0.740 | 3.4e-74 | 155 |
| Zygotene | 0.604 | 1.2e-68 | 297 |
| Pachytene | 0.766 | 7.2e-88 | 156 |
| Diplotene | 0.715 | 1.1e-62 | 138 |

**Note on mouse flanking:** Unlike human (where flanking regions are truly unique and show weak/no correlation), mouse flanking regions retain substantial inter-chromosomal similarity (138–297 nonzero Jaccard pairs). This reflects the extent of mouse subtelomeric repeats, which span nearly the entire 1Mb window — the 100kb flanking region adjacent to the PHR boundary is still within the repeat zone. The flanking "negative control" concept does not apply straightforwardly to mouse, where TLC and L1-LINE repeats extend far beyond the PHR boundary.

### 15.7 Window size optimization: 1Mb → 2Mb → 4Mb

**What it does.** Tests whether larger subtelomeric windows capture the full extent of mouse PHR regions. At 1Mb, 61% of PHRs saturate the window (>=900kb); at 2Mb, 51% still saturate (>=1.8Mb). Mouse acrocentric p-arms have massive subtelomeric repeats that extend beyond 4Mb.

**Result.** PHR size comparison across window sizes:

| Window | n PHR | Median | Mean | Max | Saturated (>=90% window) |
|--------|-------|--------|------|-----|--------------------------|
| 1Mb | 49 | 980 kb | 688 kb | 995 kb | 30/49 (61%) |
| 2Mb | 49 | 1,730 kb | 1,299 kb | 1,995 kb | 25/49 (51%) |
| 4Mb | 49 | 2,470 kb | 2,300 kb | 3,990 kb | 19/49 (39%) |

Even at 4Mb, 19/49 PHRs (39%) saturate the window. The p-arms of acrocentric chromosomes have truly massive subtelomeric repeat regions — many exceeding 4Mb. The 2Mb window represents a practical balance: it resolves 6 of the 30 arms that were truncated at 1Mb, while keeping computational costs manageable.

**Community-based Hi-C validation across window sizes (50kb resolution):**

| Window | Stage | B/W ratio | p-value | Mantel rho |
|--------|-------|-----------|---------|------------|
| 1Mb | Leptotene | 0.119 | 1.9e-28 | -0.456 |
| 1Mb | Zygotene | 0.139 | 6.3e-52 | -0.492 |
| 1Mb | Pachytene | 0.126 | 5.7e-46 | -0.473 |
| 1Mb | Diplotene | 0.149 | 3.1e-17 | -0.390 |
| 2Mb | Leptotene | 0.121 | 1.2e-43 | -0.513 |
| 2Mb | Zygotene | 0.152 | 1.2e-69 | -0.524 |
| 2Mb | Pachytene | 0.139 | 2.7e-57 | -0.522 |
| 2Mb | Diplotene | 0.156 | 1.2e-32 | -0.467 |
| 4Mb | Leptotene | 0.216 | 2.1e-36 | -0.637 |
| 4Mb | Zygotene | 0.265 | 8.1e-53 | -0.620 |
| 4Mb | Pachytene | 0.267 | 3.9e-40 | -0.631 |
| 4Mb | Diplotene | 0.274 | 2.9e-19 | -0.578 |

All analyses use per-haplotype treatment (B6 = PATERNAL, CAST = MATERNAL kept separate) and PHR-specific coordinates from the per-window-size PHR TSV. B/W ratios are 0.12–0.16 at 1Mb and 2Mb, increasing to 0.18–0.21 at 4Mb as more non-sharing sequence is included. The Mantel rho strengthens with window size (1Mb: -0.39 to -0.49; 2Mb: -0.47 to -0.52; 4Mb: -0.58 to -0.64), showing that larger windows better capture the full extent of sequence similarity between arm pairs.

**Window-size optimization (PHR saturation at p95 identity):**

| Window | n PHR | Median | Mean | Max | Saturated (≥90%) |
|--------|-------|--------|------|-----|-------------------|
| 1Mb | 49 | 0.98 Mb | 0.69 Mb | 1.0 Mb | 30/49 (61%) |
| 2Mb | 49 | 1.85 Mb | 1.30 Mb | 2.0 Mb | 25/49 (51%) |
| 4Mb | 49 | 2.53 Mb | 2.30 Mb | 4.0 Mb | 19/49 (39%) |
| 10Mb | 49 | 3.42 Mb | 3.92 Mb | 10.0 Mb | 8/49 (16%) |
| 15Mb | 50 | 3.46 Mb | 4.32 Mb | 15.0 Mb | 4/50 (8%) |
| 33Mb | 50 | 2.14 Mb | 3.37 Mb | 18.9 Mb | 0/50 (0%) |

At 33Mb (maximum feasible = smallest chromosome / 2), zero PHRs saturate. The largest PHR is CAST chr1_p at 18.9 Mb. The median (2.14 Mb) stabilizes between 10-33 Mb, confirming that most PHRs are 1-4 Mb in extent.

**Identity threshold effect at 33Mb:**

| Identity | n PHR | Median | Mean | Max |
|----------|-------|--------|------|-----|
| ≥95% | 50 | 2.14 Mb | 3.37 Mb | 18.9 Mb |
| ≥96% | 48 | 1.72 Mb | 2.64 Mb | 18.9 Mb |
| ≥97% | 42 | 1.83 Mb | 2.60 Mb | 18.9 Mb |
| ≥98% | 38 | 0.88 Mb | 1.84 Mb | 18.7 Mb |

Raising the identity threshold from 95% to 98% reduces the number of detected PHRs by 24% (50→38) and halves the median size (2.14→0.88 Mb), showing that the outer portions of mouse subtelomeric PHRs consist of moderately diverged (95-98%) repeats, while a core of high-identity (≥98%) sharing extends 0.5-2 Mb from the telomere.

### 15.7b Mouse flanking Hi-C validation (1Mb window)

**What it does.** Tests whether unique-sequence flanking regions (100kb centromere-ward of mouse PHR boundaries) also show community-structured 3D clustering, controlling for multi-mapping artifacts in the repeat-rich mouse subtelomeres.

**How.** Same pipeline as human flanking (§12.2): extract 100kb regions immediately centromere-ward of each PHR boundary, run community-based W/B bootstrap and Mantel test at all 5 resolutions. Note: unlike human flanking regions, mouse flanking regions retain substantial inter-chromosomal similarity (see §15.6 note), so the "unique sequence negative control" interpretation is weaker for mouse.

**Result.** Flanking B/W ratios and Mantel rho across all resolutions (1Mb window):

| Stage | 5kb B/W | 10kb B/W | 20kb B/W | 50kb B/W | 100kb B/W |
|-------|---------|----------|----------|----------|-----------|
| Leptotene | 0.796 | 0.893 | 0.867 | 0.519 | 0.473 |
| Zygotene | 1.086 | 1.138 | 1.176 | 1.575 | 0.647 |
| Pachytene | 0.965 | 1.197 | 1.028 | 1.350 | 0.597 |
| Diplotene | 0.774 | 0.951 | 0.804 | 1.073 | 0.122 |

| Stage | 5kb Mantel | 10kb Mantel | 20kb Mantel | 50kb Mantel | 100kb Mantel |
|-------|-----------|-------------|-------------|-------------|--------------|
| Leptotene | -0.416 | -0.409 | -0.372 | -0.372 | -0.360 |
| Zygotene | -0.468 | -0.474 | -0.452 | -0.412 | -0.353 |
| Pachytene | -0.471 | -0.447 | -0.424 | -0.399 | -0.397 |
| Diplotene | -0.409 | -0.384 | -0.386 | -0.294 | -0.084 (ns) |

**Conclusion.** The Mantel test shows significant negative correlations at all stages and resolutions (except diplotene at 100kb), confirming that the sequence similarity-contact relationship extends into flanking regions. The B/W ratios are closer to 1.0 than for PHR regions (PHR B/W 0.04-0.10 vs flanking 0.47-1.58), reflecting weaker community enrichment in flanking regions. The signal is most consistent at fine resolutions (5-20kb). As noted in §15.6, mouse flanking regions are not truly unique sequence — they still contain inter-chromosomal repeats — so the flanking control is less clean than in human.

**Multi-resolution mouse community-based B/W ratios (1Mb window, all 5 resolutions):**

| Stage | 5kb | 10kb | 20kb | 50kb | 100kb |
|-------|-----|------|------|------|-------|
| Leptotene | 0.036 | 0.074 | 0.096 | 0.119 | 0.044 |
| Zygotene | 0.031 | 0.032 | 0.071 | 0.139 | 0.117 |
| Pachytene | 0.024 | 0.037 | 0.084 | 0.126 | 0.069 |
| Diplotene | 0.028 | 0.044 | 0.102 | 0.149 | 0.095 |

**Multi-resolution mouse community-based B/W ratios (4Mb window, all 5 resolutions):**

| Stage | 5kb | 10kb | 20kb | 50kb | 100kb |
|-------|-----|------|------|------|-------|
| Leptotene | 0.176 | 0.178 | 0.183 | 0.182 | 0.178 |
| Zygotene | 0.185 | 0.187 | 0.192 | 0.204 | 0.213 |
| Pachytene | 0.184 | 0.192 | 0.193 | 0.202 | 0.202 |
| Diplotene | 0.179 | 0.182 | 0.185 | 0.208 | 0.213 |

All B/W ratios are well below 1.0 at every resolution for both window sizes. The 1Mb window shows slightly stronger enrichment at finer resolutions (B/W 0.04-0.05 at 5kb vs 0.08-0.11 at 100kb), while the 4Mb window shows more stable ratios across resolutions.

### 15.8 Implications

**Conclusion.**
1. Mouse subtelomeres share sequence across chromosomes over large regions (440–495 kb), much larger than the TLC repeat alone.
2. With proper k-mer filtering (do_not_overfilter branch), the sharing is more structured than initially apparent: 10 multi-arm communities including 6 cross-strain (B6+CAST) pairs, not just a uniform CAST bulk.
3. The CAST chr1/chr2 pair (J=0.980) persists — nearly complete subtelomeric identity across ~490 kb, analogous to human acrocentric p-arms.
4. The mouse private pairs do NOT correspond to human subtelomeric regions — syntenic human positions are in chromosome interiors. Mouse sharing is driven by repeat architecture (TLC, L1-LINE), not syntenic conservation of a duplicon system.
5. Flanking region analysis confirms the signal extends into adjacent regions, though the "unique sequence" control is less clean in mouse than in human due to the massive extent of mouse subtelomeric repeats.

---

## 15b. RPE-1 self-vs-self subtelomeric pipeline

### 15b.1 Motivation

**What it does.** Discovers RPE-1's own subtelomeric community structure from its diploid assembly, then validates it with Hi-C/Pore-C — testing whether a single individual's subtelomeric similarities predict its own 3D nuclear organization, without relying on population-level communities.

**How.** Same pipeline as mouse (§15): telomere detection → 500kb flank extraction → wfmash all-vs-all → impg/PHR detection → pggb → odgi similarity → Leiden community detection → Hi-C validation. RPE-1 is a near-diploid human cell line with a known t(X;10) translocation (chrX_HAP1 carries chr10q material — a natural positive control for cross-chromosome sequence sharing).

### 15b.2 Flank extraction and alignment

**Key metrics.** 46 chromosomes (chr1_HAP1..chrX_HAP2), 92 flanks (46 × 2 ends), PanSN naming (RPE1#1#chr*_parm, RPE1#2#chr*_qarm). Telomeres: 0–6,329 bp trimmed per flank. wfmash `-p 95` (do_not_overfilter branch), 6,410 total alignments.

**Result.** 56/92 flanks (61%) have inter-chromosomal signal at ≥95% identity, `--min-count 1` (2 haplotypes per arm). All 46 arms are represented; flanks with no signal are exclusively centromere-proximal regions of short acrocentric p-arms.

### 15b.3 Pangenome graph and similarity

**What it does.** Computes pairwise Jaccard similarity among all 92 flanks.

**How.** pggb on full flanks (not PHR subregions). Because RPE-1 has a single sample prefix (`RPE1#`), pggb's internal wfmash `-T RPE1#` combined with `-Y #` (exclude same-sample self-mapping) produces zero mappings. **Workaround**: wfmash run externally with `-T RPE1#1#` and `-T RPE1#2#` separately (one per haplotype), then combined PAF fed to pggb via `-a`.

**Key metrics.** 1,267 wfmash alignments, pangenome graph: 179,067 nodes, 246,389 edges, 92 paths. Similarity matrix: 8,464 pairs, 992 non-zero Jaccard.

### 15b.4 Community detection

**What it does.** Discovers arm-level communities from RPE-1's own subtelomeric similarity.

**How.** `detect_communities.R --organism human --level arm`. Leiden clustering on the arm-level Jaccard distance matrix (46 arms from 23 chromosomes × 2 arms).

**Result.** 37 communities from 46 arms. 5 multi-arm communities:

| Community | Members | Interpretation |
|---|---|---|
| C2 | **chr10_q, chrX_q** | **t(X;10) translocation detected**: chrX_HAP1 carries translocated chr10q |
| C9 | chr14_p, chr15_p, chr21_p, chr22_p | Acrocentric p-arms (NOR-bearing), same as human C7 |
| C20 | chr1_p, chr5_q, chr6_q | Cross-chromosome sharing |
| C18 | chr3_q, chr9_q, chr19_p | Cross-chromosome sharing |
| C13 | chr7_p, chr16_q | Same as human C3 (f7501 arms) |

**Key finding.** The t(X;10) translocation is independently discovered by the pipeline: chrX_q and chr10_q share subtelomeric sequence because chrX_HAP1 physically carries chr10q material. This serves as a positive control validating the method.

### 15b.5 Hi-C validation against self-discovered communities

**What it does.** Tests whether RPE-1's self-discovered community structure predicts 3D contacts in its own Hi-C/Pore-C data.

**How.** Three datasets (async CiFi, async Pore-C, mitotic CiFi), 50 kb resolution, per-haplotype (92 arms). Communities from §15b.4. Same 3-test framework: W/B ratio, Mantel test, bootstrap permutation.

**Result.**

| Dataset | W/B ratio | Global p | Mantel rho | Mantel p | Sig communities (BH q<0.05) |
|---|---|---|---|---|---|
| Async CiFi | 104.9x | 4.5e-113 | -0.377 | <1e-300 | 37/37 |
| Async Pore-C | 97.2x | 1.5e-86 | -0.486 | <1e-300 | 5/37 (C2, C9, C13, C18, C20) |
| Mitotic CiFi | 59.6x | 1.0e-95 | -0.232 | <1e-300 | 37/37 |

**Interpretation.** The W/B ratios (60–105x) far exceed those from HPRC-community-based RPE-1 analysis (§12.14: 32–36x). This is because self-discovered communities are more specific to RPE-1's own genome — the HPRC communities are population-level averages that may not perfectly match RPE-1's specific subtelomeric architecture.

**Caveats.** Most of the 37 communities contain only 2 arms (one per haplotype from the same chromosome). For these 2-arm communities, within-community contacts are intra-chromosomal homolog contacts (e.g., chr1_HAP1_p ↔ chr1_HAP2_p), which conflates homolog pairing with subtelomeric community structure. The biologically meaningful signal comes from the 5 multi-arm communities (C2, C9, C13, C18, C20) that group arms from different chromosomes. The Mantel test, which operates on the continuous distance matrix rather than binary community labels, is not affected by this issue.

### 15b.6 Comparison: self-discovered vs HPRC communities

| Metric | HPRC communities (§12.14) | Self-discovered (§15b.5) |
|---|---|---|
| Communities | 15 (population-level) | 37 (RPE-1-specific) |
| Multi-arm | 15 (all) | 5 |
| W/B ratio (async CiFi) | 35.8x | 104.9x |
| Mantel rho (async CiFi) | -0.457 | -0.548 |
| Mantel rho (async Pore-C) | -0.611 | -0.684 |
| t(X;10) detected? | No (population averages) | Yes (C2: chrX_q + chr10_q) |

**Conclusion.** Self-discovered communities capture individual-specific features (t(X;10) translocation) invisible to population-level communities. The Mantel correlations are comparable, confirming that the continuous similarity-contact relationship is robust to community definition. The inflated W/B ratio in self-discovered mode reflects the dominance of intra-chromosomal 2-arm communities rather than stronger inter-chromosomal signal.

### 15b.7 Flanking region control (100kb centromere-ward)

**What it does.** Tests whether unique-sequence regions immediately centromere-ward of PHR boundaries also show similarity-contact correlation. This controls for multi-mapping artifacts: flanking sequences are unique (not duplicated), so any signal must reflect genuine 3D proximity, not alignment ambiguity.

**How.** For each of the 56 RPE-1 PHR sequences, the 100kb immediately centromere-ward was extracted (`extract_flanking_sequences.py`), yielding 50 flanking sequences (6 skipped: too short or at sequence boundary). These were processed through the same pipeline: wfmash → pggb → odgi similarity → `sequence_hic_correlation.py`. Graph: 10,886 nodes, 14,803 edges, 50 paths. Only 94 non-zero Jaccard pairs (vs 992 for PHR flanks) and only 5 inter-chromosomal pairs with nonzero Jaccard.

**Result.**

| Dataset | Technology | Seq pairs | ρ (50kb) | p (50kb) | ρ (10kb) | p (10kb) |
|---|---|---|---|---|---|---|
| Async CiFi | CiFi | 1,177 | -0.011 | 0.717 | -0.006 | 0.828 |
| Async Pore-C | Pore-C | 1,177 | 0.127 | 1.3e-5 | 0.232 | 8.7e-16 |
| Mitotic CiFi | CiFi | 1,177 | -0.010 | 0.742 | -0.006 | 0.828 |

**Conclusion.** Flanking regions show near-zero sequence similarity across chromosomes (only 5 inter-chromosomal pairs with Jaccard > 0), so the similarity-contact correlation is absent or very weak. The Pore-C dataset shows a weak positive correlation (ρ = 0.13–0.23, driven by the few nonzero pairs), while CiFi shows none. This confirms that the strong PHR correlation (ρ = 0.30–0.44) is driven by shared subtelomeric sequence, not by a general property of chromosome-end regions. The comparison mirrors the arm-level flanking analysis (§12.2): flanking regions show community-structured 3D clustering (from the Mantel/W/B tests using arm-level aggregation), but at the sequence level the similarity signal is too sparse to drive a correlation.

### 15b.8 Flanking results (updated)

**Key metrics.** Flanking pggb: 4,625 similarity pairs (68 sequences). Flanking community-free correlations:

| Dataset | rho |
|---------|-----|
| Async CiFi | -0.008 |
| Async Pore-C | 0.136 |
| Mitotic CiFi | -0.007 |

**Conclusion.** Flanking community-free rho ≈ 0 across all datasets, confirming that the PHR community-free signal (§12.5b) is driven by shared subtelomeric sequence content rather than generic chromosome-end proximity effects.

---

## 16. Sample composition and limitations

### Sample composition

**Key metrics.** 232 HPRCv2 individuals: AFR = 67 (28.9%), EAS = 52 (22.4%), AMR = 44 (19.0%), SAS = 37 (15.9%), EUR = 32 (13.8%). Population-level findings (§10.3) should be interpreted with this imbalance in mind. Within-superpopulation heterogeneity is not modeled.

### Methodological limitations

1. **Identity threshold**: The 95% minimum identity threshold for both all-vs-all alignment and inter-chromosomal region detection means that more divergent inter-chromosomal homology (older exchanges) is invisible. The 15.1% cross-arm affinity rate and the set of 7 arms with no inter-chromosomal signal represent lower bounds, not complete inventories. Ambrosini et al. (2007) identified a bimodal distribution of duplicon identity with peaks at 91% and 98%; the 95% threshold captures primarily the recent high-identity peak and misses most of the older 91% peak, which includes many olfactory receptor and immunoglobulin-related duplicon families. Note that Ambrosini et al. (2007) distinguished "subtelomere-only" duplicon blocks (Table 1, found exclusively at subtelomeric regions) from "subterminal" duplicon families (Table 2, positioned adjacent to terminal TTAGGG tracts but sometimes having non-subtelomeric copies as well); this distinction is not modeled in the present analysis, which treats all inter-chromosomal similarity uniformly.

2. **Flank size**: The 500 kb maximum flank extraction may truncate longer subtelomeric similarity regions. No sensitivity analysis on this parameter was performed.

3. **Region length threshold**: The 3 kb minimum output region length and 5 kb window/step size set a floor on detectable inter-chromosomal regions. Shorter shared segments would be missed.

4. **Assembly quality**: Subtelomeric regions are among the most difficult to assemble. Assembly gaps, collapses, or errors near telomeres could affect sequence content and inter-chromosomal signal. The telomere-presence filter mitigates but does not eliminate this risk. The chr18_q chimeric contig in NA18982#1 (§4) illustrates this concern: JBKABS010000018.1 fuses chr18 with 966 kb of chrX PAR1 across a 100 bp NNN scaffold join, and no separate chrX contig exists in this haplotype.

5. **Community detection resolution**: The 50-community Leiden solution (k-NN = 75, resolution = 0.8) is one of a family of possible solutions. Different resolution parameters would yield different community numbers and compositions. The solution was selected by modularity optimization within a 5–50 community range, but modularity does not guarantee biological correctness.

6. **Small sample sizes**: Some cross-arm affinity findings rest on small sample sizes: chrY_p cross-arm N=5 (cross-arm rate 5.6%, 95% Wilson CI [2.4%, 12.4%]), chr16_q cross-arm N=27 (6.0%, CI [4.2%, 8.6%]), chr14_p N=75 (100%, CI [95.1%, 100.0%]). These should be considered preliminary pending validation with larger or independent datasets.

7. **Exchange timing**: Cross-arm affinity demonstrates that exchange has occurred but cannot date individual events. High discordance rates are consistent with recurrent exchange but could also reflect a single ancient event still segregating. Distinguishing these scenarios requires trio/family data or population genetic modeling.

8. **Somatic exchange in cell lines**: Mefford & Trask (2002, citing van Overveld et al. 2000) noted that "some people are mosaic for 4q/10q subtelomeric translocations, which indicates that subtelomeric sequences can interchange in somatic cells." Since most HPRCv2 assemblies derive from lymphoblastoid cell lines (LCLs), some cross-arm affinity — particularly at chr4_q/chr10_q (C1) — could reflect somatic exchange during cell culture rather than germline polymorphism. This caveat applies to all LCL-derived assemblies and is difficult to quantify without matched blood-derived controls.

### 3D validation limitations

9. **Somatic vs meiotic context**: All 3D data (Hi-C, Pore-C, Dip-C) captures somatic interphase organization. The meiotic bouquet stage — when ectopic recombination between subtelomeric arms would occur — has never been captured by Hi-C in humans. The observed somatic 3D signal is interpreted as a residual of meiotic chromosome organization (Rabl configuration — the retained centromere-telomere polarity from the preceding cell division, where centromeres cluster at one nuclear pole and telomeres at the other), but this interpretation is indirect. Human meiotic Hi-C remains the single most informative missing experiment.

10. **Sample size (N=6 Hi-C)**: Five diploid HPRC samples plus one haploid (CHM13) are sufficient to demonstrate that the 3D signal exists and is reproducible, but insufficient for population-level claims about 3D variation. CHM13 shows no significant communities, primarily due to reduced power (37 arms vs 75 in diploid samples).

11. **GM12878 cell line**: The Dip-C data uses GM12878 (an EBV-transformed B-lymphoblastoid cell line), which has an abnormal karyotype and may not represent normal nuclear organization. PBMC results provide a primary-cell control but with fewer cells.

12. **hg19/T2T coordinate incompatibility**: Dip-C data (Tan et al. 2018) uses hg19 coordinates; PHR regions are defined on CHM13/T2T. Coordinate projection via impg partially mitigates this, but coarse 50 kb resolution introduces noise, particularly near assembly gaps and subtelomeric regions where hg19 and T2T differ most.

13. **Multi-mapping at PHR intervals**: Duplicated PHR sequences cause Hi-C/Pore-C reads to multi-map between community partner arms, inflating apparent inter-chromosomal contacts. The flanking-region control (unique sequence) addresses this by demonstrating that the signal persists — and is stronger — in regions without multi-mapping. However, the PHR-specific enrichment values (B/W 0.027–0.074) cannot be separated from multi-mapping contribution.

14. **Confound controls**: The acrocentric exclusion control (§12.9) rules out nucleolar association as the driver. Rabl configuration (centromere-telomere polarity causing generic telomere clustering) is addressed by the flanking analysis (§12.2): if Rabl drove the signal, flanking regions farther from the telomere tip should show weaker signal, but they show stronger signal. Chromosome size effects (small chromosomes intermingle more) are addressed by the Mantel test (§12.3), which tests a continuous correlation between sequence similarity and Hi-C contact across all arm pairs — a size confound would add noise but not create a correlation between Jaccard distance and contact frequency. The significant Mantel results (HG002 Hi-C rho = −0.41, Pore-C rho = −0.50; §12.3) demonstrate that the 3D signal is specifically tied to subtelomeric sequence similarity.

15. **Parameter sensitivity**: The 50 kb resolution and 500 kb flanking window were selected based on optimization but not subjected to formal sensitivity analysis. Different resolutions or window sizes could yield different enrichment values. Multi-resolution analysis at all 5 mcool resolutions (5kb, 10kb, 20kb, 50kb, 100kb) across human, RPE-1, and mouse systems (§12.7, §12.9, §12.14.1, §15.7b) demonstrates that the core signal is robust to resolution choice.

16. **Fragmented assemblies produce NaN flanking values**: HG02148 and NA19036 assemblies are sufficiently fragmented at some subtelomeric regions that flanking-region coordinates fall outside contig boundaries, producing NaN values in the community-free flanking correlation. These samples are excluded from flanking analyses but included in PHR-based analyses.

17. **Dip-C cell 12 duplicate**: Cell 12 produces identical 3dg output to cell 10 (shared SRR7226706 long-insert library). Cell 12 is excluded from all analyses (16 cells used in §13).

18. **Mouse 1 Mb PHRs mostly fill the extraction window**: At the 1 Mb extraction scale, mouse PHR regions have a median length of ~980 kb, meaning the PHR nearly saturates the window. Even larger windows (1.5–2 Mb) may reveal additional inter-chromosomal similarity beyond the current detection boundary, particularly for the p-arm mega-community (C1, 35 arms) where sharing extends deep into the chromosome.

---

The preceding sections established subtelomeric community structure from sequence similarity. The following sections test whether this structure has a physical counterpart in 3D nuclear organization.

## 17. Literature context, novel contributions, and testable predictions

### 17.1 Confirmed literature claims

**What it does.** Systematically compares present findings against specific published claims.

**Result.** Key confirmations:

- **Ambrosini duplicon architecture**: All 11 subtelomere-specific duplicon block entries (Table 1; numbered 1–3, 5–8, 10–12, with block 6' as a variant of block 6; blocks 4 and 9 do not exist in the original table) map systematically to the 15 Leiden communities (§9), with 4 confirmed at gene level (OR4F→C8, DUX4L→C1, IQSEC3→C5, IL9RP1→C3). Block 5 (chr2_p anchor) has no PHR signal, consistent with chr2_p's exclusion for lacking inter-chromosomal signal. Three genes with functional evidence noted by Ambrosini — FBXO25 (block 10, chr2_q), TUBB4q (block 8, chr4_q), and RYD5 (block 6', chr11_p) — are not detected in current Liftoff annotations, likely reflecting nomenclature changes or location below the annotation threshold. The bimodal identity distribution (91%/98% peaks) explains why the 95% threshold captures primarily recent exchanges (§16, limitation 1). Internal (TTAGGG)n islands (18,352, §8) are telomere-proximal elements distributed throughout the duplicated zone. Note: Ambrosini's specific claim that islands co-localize with internal duplicon-to-duplicon boundaries was not tested here — the §8 boundary enrichment test evaluated the PHR outer boundary, a different feature. Islands are distinct from TAR1 (§8). TAR1 positional distribution separately shows telomere-proximal bias (66.9% within 10 kb of telomere; §8).
- **Flint/Mefford two-domain model**: **The prediction that inter-chromosomal sharing decreases with distance from the telomere (Flint et al. 1997; Mefford & Trask 2002) is confirmed on 39/48 arms, with a discrete two-phase breakpoint structure (rather than smooth decay) on 39/41 testable arms (§10.11).** Internal (TTAGGG)n islands co-localize with the domain boundary within 25 kb on 11/19 arms with detectable ITS, consistent with the model's prediction that TTAGGG tracts separate the distal and proximal domains. All 6 focus arms (chr4p, chr4q, chr16p, chr18p, chr20p, chr22q) originally characterized by Flint and Mefford show significant gradient and breakpoint support.
- **Mefford & Trask subtelomeric exchange**: **The f7501 block African enrichment at chr7p/chr16_q (Mefford 2002) is confirmed quantitatively (chr16_q cross-arm 78% AFR; §10.4).** **IL9R pseudogene distribution across chr9_q/chr10_p/chr16_p/chr18_p matches across communities C2/C3/C9 (C2 covers chr10_p and chr18_p, C3 covers chr9_q, C9 covers chr16_p); IL9R copies also appear in C14/C15 (PAR).** **The chr4_q/chr10_q exchange frequency is detected at higher sensitivity (42.8% discordance vs 20% by Southern blot), reflecting the broader detection of subtelomeric sequence affinity.**
- **Zuo meiotic alignment**: ***Zuo et al.'s ~20% chromosome-end alignment zone and ~500 kb leptotene loop size (confirmed in mouse meiosis) are consistent with the median 105 kb PHR region fitting within a single meiotic loop (§14.3). The inference that this positions PHR sequences where ectopic recombination can occur is the present analysis's synthesis, not a Zuo et al. finding. Note: this extrapolation is cross-species (mouse → human) and cross-context (meiotic → somatic).***
- **Tan Dip-C radial positions**: ***Per-chromosome radial preferences from Tan et al.'s Dip-C data (§13.3) are consistent with cell-biology predictions when aggregated by community (C1 peripheral = D4Z4-proximal lamin tethering; C6 interior = nucleolar proximity; §12.13). The community-level aggregation is the present analysis, not Tan et al.'s.***

Additional quantitative confirmations from literature claim tests:
- **Mefford chr16_p bimodality**: chr16_p region lengths are bimodal — 312 short alleles (71.2%, median 25 kb) and 126 long alleles (28.8%, median 205 kb). The long allele fraction (28.8%) matches Mefford's ~30% prediction. Long alleles are AFR-enriched (59.5%) and cluster in C28.
- **Mefford chr3_q homogeneity**: chr3_q within-arm Jaccard distance variance (0.022) is 2.8x lower than chr15_q (0.063) and 1.6x lower than chr19_p (0.036), consistent with Mefford's observation that 44/46 sequenced chr3_q alleles were identical — which they found "surprising" and attributed to either a population bottleneck or a selective sweep. Note: chr3_q PHR regions are shorter (median 250 kb) than chr15_q (125 kb) or chr19_p (205 kb), which could contribute to lower absolute variance; however, the coefficient of variation (CV) for chr3_q (1.84) is higher than chr15_q (1.39) and chr19_p (0.83), arguing against a pure length artifact.
- **Mefford IL9R distribution**: IL9RP1 at chr9_q (445 haplotypes), IL9RP3 at chr16_p (432), IL9RP4 at chr18_p (444), IL9RP2 at chr10_p (33), matching the four arms reported by Mefford & Trask (2002). Additionally, IL9RP4 at chr13_q (3 haplotypes) and protein-coding IL9R at chrX_q (311), chrY_q (76), chrY_p (2), chr16_p (3). Cross-arm sequences on IL9R arms carry IL9R at 100% rate where testable (chr13_q 3/3, chrX_q 311/311, chrY_p 5/5).
- **chr19_p structural polymorphism**: chr19_p is part of the 6-arm C3 community (f7501 cluster). At the 50-community sequence level, chr19_p sequences are distributed across 7 communities, the most of any arm after chr6_q (8). Linardopoulou et al. (2005) noted that "The 19p alleles also differ grossly in subtelomeric content." The population-scale quantification across the pangenome confirms this extensive polymorphism.
- **D4Z4 as causal driver of C1**: Inter-chromosomal signal peaks at 0–15 kb (D4Z4 position), C1 median 22 DUX4L vs all 7 outliers 0–2 on their own arm (Mann-Whitney p = 5.3e-6), outlier PHR regions 4.6–9x shorter.
- **OR4F pseudogenization gradient**: 62.1% of 5,023 OR4F annotations are pseudogenes, with pseudogenization rate varying from 11.1% (chr7_p) to 99.8% (chr15_q) across 16 arms.

**Conclusion.** No tested paper claims were directly contradicted. The 42.8% vs 20% discrepancy at chr4_q/chr10_q reflects methodological sensitivity, not conflicting biology.

### 17.2 Novel contributions

**What it does.** Identifies 24 findings that go beyond the published literature — quantitative extensions of known biology and new observations:

1. **Population-scale community structure**: 41 arms clustering into 15 communities across 232 individuals — the first population-scale quantification showing that single-genome duplicon architecture organizes into discrete, reproducible communities across 465 near-complete assemblies.
2. **Three-category arm classification**: Homogeneous (7/41), polymorphic (34/41, up to 7 types), and fully interchangeable (acrocentric p-arms). This extends the qualitative patchwork model of Mefford & Trask (2002) with a quantitative framework.
3. **Population-scale cross-arm affinity**: 15.1% of sequences (2,271/15,668) resemble a foreign arm more than their own — the first population-scale quantification across 465 near-complete assemblies.
4. **Subtelomeric type discordance quantified at population scale**: Up to 46.1% of individuals carry two structurally different subtelomeric types at the same locus — the first population-scale quantification of structural heterozygosity at subtelomeres, extending qualitative observations from Mefford & Trask (2002) and Linardopoulou et al. (2005).
5. **Gene repertoire replacement scores**: Complete (0.91–1.0) at chr13_p/chr14_p/chr15_p and PAR, partial (0.0–0.72) at autosomal communities — a gradient of homogenization mapped onto the community structure for the first time.
6. **3D genome mirrors sequence communities**: Three independent technologies provide evidence consistent with community-structured 3D co-localization (Hi-C B/W 0.027–0.074 enrichment, Dip-C 1.8–4.9% closer), though effect sizes vary substantially across samples and the finer 50-community partition does not reach significance in Dip-C.
7. **Flanking region paradox**: Unique-sequence regions 100 kb centromere-ward show stronger 3D signal than duplicated PHR (HG002: flanking (B/W=0.002) vs PHR (B/W=0.027) enrichment for HG002), providing strong evidence against multi-mapping artifact and suggesting 3D clustering involves broader chromosomal domains.
8. **C4 minimal-PHR positive control**: 5–25 kb of shared sequence at chr7_q/chr12_q tips coincides with reproducible 3D co-localization across 4/5 samples, demonstrating that even minimal shared subtelomeric sequence is associated with detectable 3D contact enrichment.
9. **Cell-type specificity**: PBMC cells (18 combined, hg19 coordinates with projected PHR boundaries) do not reach significance (W/B = 0.983, p = 0.305). The non-significance reflects hg19 projection noise, mixed cell types, small N, and smaller PHR-specific regions. GM12878 T2T (6.9% closer, p = 2.4e-5 on native coordinates) is the primary Dip-C finding.
10. **Per-individual discordance↔3D correlation**: Low discordance → strong Hi-C signal (rho = −0.50, p = 3.4e-4, N = 48 sample × community pairs), consistent with a model in which intact duplicated sequence on both haplotypes enhances the 3D contact signal.
11. **Proposed feedback loop model**: Sequence similarity → 3D proximity → ectopic exchange → increased similarity, with varying levels of support for each link. The causal direction cannot be established from the present data (§14.6).
12. **MTCO pseudogene enrichment**: Population-scale characterization of mitochondrial pseudogenes at acrocentric p-arms (C7), enabled by T2T-quality assemblies of these previously unresolved regions.
13. **TAR1 near-absence from PAR1**: 0.5% prevalence at chrX_p/chrY_p vs 94.8% genome-wide, consistent with PAR1 using obligate meiotic crossover rather than repeat-mediated exchange.
14. **B-compartment with interior positioning**: Subtelomeric regions are weakly B-compartment (63% of tips) yet positioned internally (radial 0.56–0.63), consistent with known interior telomere positioning in lymphocytes and now quantified at community-level resolution across the pangenome.
15. **Subtelomeric exchange frequencies as phylogenetic markers**: Population tree from cross-arm affinity frequencies across 9 arms recovers the out-of-Africa topology (AFR deepest split, AMR-EUR closest). **This supports Mefford & Trask's (2002) hedged suggestion that "measuring the frequency of specific subtelomeric blocks on particular chromosomes might inform phylogenetic studies of modern humans,"** though they cautioned that such frequencies "could be an unreliable indicator of the relationships of human populations" because they also reflect exchange history (§17.3, prediction 6).
16. **TAR1 prevalence is consistent with passenger status**: TAR1 is near-universal outside PAR (99.6%), making its functional role in exchange difficult to assess. Where testable (C7 acrocentric p-arms), cross-arm affinity rates are indistinguishable between TAR1+ and TAR1− sequences (86.6% vs 90.5%), consistent with TAR1 being a passenger rather than a facilitator of exchange, though the near-universal prevalence precludes definitive conclusions (§17.3, prediction 3).
17. **Internal (TTAGGG)n islands quantified at population scale**: 18,352 islands across 8,321 sequences (53.1%) and all 41 arms, with median length 79 bp. These are distinct from TAR1 and distinct from terminal telomeric arrays. Ambrosini et al. (2007) reported that internal (TTAGGG)n-like sequences "almost always co-localize to duplicon boundaries," suggesting "their involvement in the generation of the complex sequence organization." The present analysis quantifies these islands at population scale but tests enrichment at PHR outer boundaries (§8), not internal duplicon-to-duplicon boundaries — the original internal boundary claim remains to be tested directly.
18. **Alleles closer than paralogs, except at fully homogenized communities**: Paired allelic distance (maternal vs paternal, same arm) is significantly smaller than paralog distance (closest different arm in same community) in 9/10 multi-arm communities (overall Wilcoxon p < 1e-300; §10.1). The sole exception is C7 (acrocentric p-arms), where 65.5% of individuals have paralog distance < allelic distance — a quantitative confirmation that these arms carry interchangeable sequences. This provides the first population-scale quantification of the relationship between allelic and paralogous subtelomeric sequences. **Ambrosini et al. (2007) observed that subtelomeric paralogs can share "higher sequence similarity than alleles on homologous chromosomes"; Mefford & Trask (2002) noted that recurrent inter-chromosomal transfers create anomalously high allelic diversity at individual chromosome ends. The present data shows this is the general rule for multi-arm communities but is reversed at fully homogenized arms (C7).**
19. **Subtelomeric Fst mirrors out-of-Africa history**: Pairwise Fst using cross-arm/self-arm as alleles across 9 arms shows AFR differentiated from all non-AFR superpopulations (Fst 0.10–0.12), while non-AFR populations are nearly undifferentiated (Fst 0.00–0.02; mean Fst = 0.048; §10.4).
20. **Internal (TTAGGG)n islands are predominantly degenerate**: Only 32.8% of islands are "pure canonical" (≥80% TTAGGG content); 46.7% are variant-dominant (<50% canonical), with TGAGGG (18.7%), TTGGGG (16.3%), and TCAGGG (12.5%) collectively accounting for 47.4% of all telomeric hexamers (§8). This quantifies at population scale the variant composition reported by Ambrosini et al. (2007), who identified the same three dominant variant motifs in the reference genome, supporting the interpretation of these islands as ancient remnants of telomeric sequence incorporated during subtelomeric duplication events.
21. **Cross-arm sequences have shifted (TTAGGG)n island distribution**: Cross-arm affinity sequences carry proportionally more centromere-ward (proximal) TTAGGG islands than self-arm sequences (48.0% vs 43.8% proximal; chi-squared = 17.93, p < 0.001, OR = 1.19; §10.9), **consistent with Ambrosini et al.'s (2007) suggestion that the subterminal compartment may be more recombinogenic.**
22. **Within-community Jaccard distance structure is multi-modal**: Communities with distinct arm-specific types (C2, C12) show clear bimodal distance distributions with separated allelic and paralog peaks, while homogenized communities (C1, C7) show diffuse distributions consistent with blurred arm identity (§10.10). This population-scale bimodal structure echoes, but is distinct from, Ambrosini et al.'s (2007) bimodal identity distribution (98%/91% peaks): their peaks reflect evolutionary timing of duplication events in the reference genome, while the present bimodality reflects allelic vs paralogous distances within population-scale communities.
23. **Two-domain subtelomeric model supported across most chromosome arms at pangenome scale**: The Flint/Mefford two-domain model — distal (many-chromosome sharing) vs proximal (few-chromosome sharing) — is supported on 39/48 arms by Spearman gradient (81%), 39/41 by breakpoint analysis (95%), and in 99.7% of individual haplotype sequences (§10.11). Breakpoints are arm-specific (10–445 kb), and internal (TTAGGG)n blocks co-localize with breakpoints within 25 kb on 11/19 testable arms. This extends the model from the handful of arms characterized by FISH (Flint et al. 1997) to the majority of the human chromosome complement at population scale.
24. **Two-domain boundary positions quantified**: For the first time, the distal/proximal domain boundary is mapped at sequence resolution across 39 chromosome arms with detectable inter-chromosomal gradients, with arm-specific breakpoints ranging from 15 to 445 kb. Boundaries at the originally characterized arms (chr4p at 70 kb, chr4q at 50 kb, chr22q at 15 kb) are consistent with the 25–50 kb range inferred from FISH studies but provide precise positions from sequence data rather than cytogenetic estimates.
25. **Mouse meiotic enrichment across all stages**: Using Zuo et al. (2021) meiotic Hi-C data at 1 Mb scale, the per-arm-pair correlation between sequence similarity and Hi-C contact peaks at zygotene (rho=0.862, p=8.1e-103, n=344 arm pairs), with leptotene 0.834, pachytene 0.845, and diplotene 0.633. The flanking-region (unique-sequence) control shows the same temporal pattern (rho=0.60–0.77), confirming the signal is not driven by multi-mapping. All four stages show significant community-based enrichment (B/W 0.119–0.149, Mantel rho=-0.39 to -0.49, all p < 0.001, per-haplotype with PHR-specific coordinates), with diplotene showing the weakest Mantel correlation (rho=-0.39) as expected for the stage where the bouquet has dissolved.
26. **Pangenome-level community-free approach validated across organisms**: The per-arm-pair correlation framework produces consistent results across human somatic (Hi-C, Pore-C, Dip-C) and mouse meiotic contexts. In mouse (1Mb windows), the per-arm-pair Spearman rho ranges from 0.633 (diplotene) to 0.862 (zygotene); flanking arm-pair rho ranges from 0.604 (zygotene) to 0.766 (pachytene). In human Dip-C (16 cells), Mantel rho = 0.296 (p = 0.002). This cross-organism consistency strengthens the inference that sequence similarity drives 3D co-localization rather than an organism-specific confound.
27. **Sperm W/B = 0.401 demonstrates signal persists in haploid cells**: Sperm data shows 60% closer within-community contact (W/B = 0.401) with PHR-specific coordinates, demonstrating strong community-structured 3D organization in haploid post-meiotic cells. The Mantel test is also significant (rho=0.202, p=0.023), consistent with Rabl configuration persisting through spermiogenesis.

### 17.3 Testable predictions from existing data

**What it does.** Identifies seven hypotheses testable with available data or minimal new analysis:

1. **LINC complex and meiotic alignment** (highest priority): Apply the community framework to Zuo et al.'s published wild-type vs SUN1 W151R mutant zygotene Hi-C data (GEO: GSE155142 for mutant; GSE155638 and GSE155967 for wild-type). Requires orthologous mapping of mouse subtelomeric regions to human community structure, which has not been established (mouse chromosomes are acrocentric with structurally distinct subtelomeres). The SUN1 mutant narrows the alignment zone from ~20% to ~5% of chromosome length while telomere anchoring is only ~50% abrogated. **Important caveat**: the median PHR region (105 kb) occupies <2% of most chromosomes, well within the ~5% tip zone where SUN1 mutant contacts *increase* rather than decrease. The predicted effect on PHR-scale contacts is therefore uncertain — tip contacts may be maintained or strengthened even without LINC-mediated force transmission.

2. **Haplotype-resolved 3D contacts** (tested): At discordant arms in 4 Hi-C samples (HG002 is fully concordant), tested whether the self-arm haplotype shows higher contact with community partners than the cross-arm haplotype. In 7/8 informative discordant pairs, the cross-arm haplotype had equal or higher partner contact — the opposite of the prediction. However, N is too small (10 pairs total, many with zero contacts, Wilcoxon p = 0.5 for NA19036's 5 pairs) to draw conclusions. The result is inconclusive; testing requires deeper Hi-C from samples with higher discordance rates.

3. **TAR1 as facilitator vs passenger** (tested): Compared cross-arm affinity rates between TAR1+ and TAR1− sequences. Excluding PAR communities (where TAR1 absence and cross-arm affinity are both driven by the PAR recombination mechanism, not TAR1 biology): only 55/12,962 non-PAR sequences lack TAR1 (99.6% have TAR1), providing insufficient statistical power. In C7 (acrocentric p-arms), the only community with both TAR1+ and TAR1− sequences, cross-arm rates are indistinguishable (86.6% vs 90.5%). **Conclusion: TAR1 is a passenger, not a facilitator of exchange** — but the near-universal presence of TAR1 outside PAR makes this difficult to test rigorously.

4. **CTCF/cohesin at PHR boundaries** (not tested): Requires ENCODE CTCF ChIP-seq download and bedtools intersection with PHR boundary coordinates. Tests whether structural protein binding predicts 3D contact strength.

5. **Somatic vs germline exchange** (not testable): HPRCv2 sample metadata (`hprc-sequence-production.tsv`) does not include a DNA source column (LCL vs blood). This test requires external metadata linking HPRCv2 sample IDs to their source material.

6. **Subtelomeric phylogenetic markers** (tested): Population tree constructed from cross-arm affinity frequencies across 9 significant arm/community pairs × 5 superpopulations. The resulting tree recovers the expected out-of-Africa topology: AFR is the deepest split (distance 0.73–0.81 from all others), followed by EAS branching, then SAS, AMR, and EUR clustering together (AMR-EUR distance 0.22, closest pair). This demonstrates that subtelomeric exchange frequencies contain phylogenetic signal consistent with known human population history, as Mefford & Trask (2002) hypothesized.

7. **Crossover frequency correlation** (tested): Using the T2T-CHM13 recombination map (Lalli et al. 2025, bioRxiv 2025.02.24.639687), mean recombination rate in the distal 500 kb of each arm was correlated with cross-arm affinity rate across 36 arms in multi-arm communities. Spearman rho = −0.43, p = 0.008.

| Arm | Community | Recomb (cM/Mb) | Cross-arm % | N sequences | Interpretation |
|-----|-----------|---------------|-------------|-------------|----------------|
| chr14_p | C7 | 0.00 | 98.7 | 76 | No recombination → complete homogenization |
| chr22_p | C7 | 0.00 | 85.7 | 119 | No recombination → near-complete homogenization |
| chrX_p | C15 | 0.00 | 90.0 | 321 | PAR1, no map data → X/Y indistinguishable |
| chrX_q | C14 | 0.14 | 99.7 | 312 | PAR2, near-zero recombination → X/Y indistinguishable |
| chr15_p | C7 | 0.07 | 100.0 | 68 | Near-zero recombination → complete homogenization |
| chr21_p | C7 | 0.03 | 55.1 | 69 | Near-zero recombination → partial homogenization |
| chr13_p | C7 | 0.46 | 95.2 | 62 | Very low recombination → near-complete homogenization |
| chr22_q | C6 | 1.89 | 55.1 | 450 | Low-moderate recombination → partial exchange |
| chr8_p | C11 | 3.29 | 66.2 | 432 | Moderate recombination → substantial exchange |
| chr4_q | C1 | 2.93 | 39.2 | 352 | Moderate recombination → partial exchange (D4Z4) |
| chr11_p | C3 | 5.54 | 61.1 | 447 | High recombination but still high exchange |
| chr6_p | C5 | 4.74 | 58.8 | 413 | High recombination but still high exchange |
| chr10_q | C1 | 1.21 | 6.4 | 357 | Low-moderate recombination → low exchange |
| chr16_q | C3 | 1.62 | 6.0 | 447 | Low-moderate recombination → low exchange |
| chr19_p | C3 | 8.50 | 0.5 | 428 | Highest recombination → arm identity maintained |
| chr7_p | C3 | 5.54 | 0.0 | 421 | High recombination → zero cross-arm affinity |
| chr15_q | C8 | 6.17 | 0.0 | 446 | High recombination → zero cross-arm affinity |
| chr12_p | C5 | 3.55 | 0.0 | 453 | Moderate-high recombination → zero exchange |
| chr9_q | C3 | 1.76 | 0.0 | 445 | Low-moderate recombination → zero exchange |

(Table shows selected arms spanning the full range; all 36 arms tested contribute to the rho = −0.43 correlation.)

**Confound control**: Acrocentric p-arms (chr13_p, chr14_p, chr15_p, chr21_p, chr22_p) have 0–12 callable variants in 500 kb (vs 1,000–3,000 for non-acrocentric arms), so their 0 cM/Mb recombination rate reflects absence of short-read genotyping data in repetitive rDNA-adjacent regions, not necessarily absence of recombination. PAR1 (chrX_p) has an obligate crossover but is mapped in a separate file. Excluding these 7 confounded arms: rho = 0.03, p = 0.86, N = 29 — **the correlation vanishes entirely**. The rho = −0.43 signal across all 36 arms is driven exclusively by the extreme values at arms where both recombination rate and cross-arm affinity reflect the repetitive, hard-to-genotype nature of acrocentric p-arms and PAR, not an independent biological relationship. The question of whether local recombination protects arm identity remains open but cannot be answered from the current recombination map data at these loci.

---

## Appendix: Files and tools used

### External tools

| Tool | Version | Path |
|------|---------|------|
| wfmash | v0.23.0-41-gb5f0ff1c | `/gnu/store/rscfhnyvlba7kv08rkl6w3nx259azs1w-wfmash-gcc-static-git-v0.23.0-41-gb5f0ff1c/bin/wfmash` |
| impg | git build (commit 5b96025) | `/moosefs/guarracino/git/impg/target/release/impg` |
| pggb | matching wfmash 0.23 | `/moosefs/guarracino/pggb_wfmash023/pggb/pggb` |
| odgi | bundled with pggb | `/moosefs/guarracino/pggb_wfmash023/smoothxg/deps/odgi/bin/odgi` |
| samtools | conda hicexplorer 3.7.4 | `/moosefs/guarracino/condatools/hicexplorer/3.7.4/bin/samtools` |
| bedtools | conda hicexplorer 3.7.4 | `/moosefs/guarracino/condatools/hicexplorer/3.7.4/bin/bedtools` |
| bgzip | conda hicexplorer 3.7.4 | `/moosefs/guarracino/condatools/hicexplorer/3.7.4/bin/bgzip` |
| pigz | system | `/usr/bin/pigz` |
| Rscript | guix | `/home/guarracino/.guix-profile/bin/Rscript` |
| python3 | conda hicexplorer 3.7.4 | `/moosefs/guarracino/condatools/hicexplorer/3.7.4/bin/python3` |

### Custom scripts

| Path | Description | Section |
|------|-------------|---------|
| `/moosefs/guarracino/HPRCv2/scripts/classify_contigs.py` | Contig classification by p/q-arm alignment and telomere filtering | 1 |
| `/moosefs/guarracino/HPRCv2/scripts/find-multichr-regions-incremental.py` | Sliding-window inter-chromosomal region detection with early stopping | 4 |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R` | Distance matrix, MDS, UMAP, Leiden/UPGMA community detection | 6, 7, 8 |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-seq-community-structure.R` | Sequence-level community structure visualization | 8 |
| `/moosefs/guarracino/HPRCv2/scripts/preprocessing/preprocess-subtelomeric-annotations.R` | Gene/TAR1 annotation intersection with PHR regions | 8 |
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze-community-enrichments.R` | Community gene enrichment testing | 9 |
| `/moosefs/guarracino/HPRCv2/scripts/community/plot-community-enrichments.R` | Enrichment result visualization | 9 |
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze-within-arm-heterogeneity.R` | Cross-arm affinity, population bias, subtelomeric type discordance, gene replacement scoring | 10 |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/test_two_domain.py` | Two-domain model gradient test (Spearman) | 10.11 |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/test_two_domain_changepoint.py` | Two-domain model changepoint + TTAGGG co-localization | 10.11 |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/extract_rm_ttaggg_subtelomeric.py` | Extract RepeatMasker telomeric annotations in subtelomeric flanks | 10.11 |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/analyze_polymorphic_arms.py` | Polymorphic arm subgroup characterization | 7.2 |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/generate_polymorphic_arm_narratives.py` | Generate polymorphic arm narrative report | 7.2 |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/test_d4z4_causality.py` | D4Z4 causality tests for C1 (signal localization, DUX4L count, outliers) | 7.2 |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/test_literature_claims.py` | 6 testable literature claims (OR4F, chr3_q bottleneck, chr16_p bimodality, IL9R, chr15q/chr19p symmetry, D4Z4) | 17.1 |
| `/moosefs/guarracino/HPRCv2/scripts/community/plot-within-arm-heterogeneity.R` | Heterogeneity result visualization | 10 |
| `/moosefs/guarracino/HPRCv2/scripts/community/compare-community-levels.R` | Arm-level vs sequence-level community comparison (QC) | — |
| `/moosefs/guarracino/HPRCv2/scripts/community/community-utils.R` | Shared utility functions for community analysis scripts | — |
| `/moosefs/guarracino/HPRCv2/scripts/community/detect_communities.R` | Shared Leiden community detection (arm-level + sequence-level); accepts odgi similarity TSV or Hi-C O/E matrix | 7, 12, 15 |
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze_hic_communities.py` | Hi-C/Pore-C community enrichment, bootstrap W/B, Mantel test, ARI comparison, per-PHR-pair correlation | 12 |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/export_arm_dist_matrix.R` | Export arm-level Jaccard distance matrix from RDS | 12, 13 |
| `/moosefs/guarracino/HPRCv2/scripts/hic/analyzer.py` | NOR chromosome clustering via MDS on inter-chromosomal contacts | 12 |
| `/moosefs/guarracino/HPRCv2/scripts/dipc/community_3d_enrichment.py` | Dip-C community 3D enrichment, Mantel test, radial analysis | 13 |
| `/moosefs/guarracino/HPRCv2/scripts/dipc/project_phr_to_hg19.py` | Project PHR regions to hg19 via impg for Dip-C overlay | 13 |
| `/moosefs/guarracino/HPRCv2/scripts/dipc/phr_dipc_overlay.py` | PHR-particle overlay: shared vs unshared PHR 3D distances | 13 |
| `/moosefs/guarracino/HPRCv2/scripts/dipc/plot_3dg.py` | 3D genome structure visualization (chromosomes, terminal, community) | 13 |
| `/moosefs/guarracino/HPRCv2/scripts/dipc/plot_gm12878_vs_pbmc.py` | GM12878 vs PBMC comparison summary figure | 13 |

### Input data

| Path | Description | Section |
|------|-------------|---------|
| `/moosefs/pangenomes/HPRCv2/*.fa.gz` | 465 near-complete assemblies analyzed (464 HPRC haplotypes from 232 individuals plus CHM13; GRCh38 in the mirror excluded) | 0 |
| `/moosefs/pangenomes/HPRCv2/chm13v2.0_maskedY_rCRS.fa.PanSN.fa.gz` | CHM13v2.0 reference (masked Y, rCRS mitochondria) in PanSN format | 0 |
| `/moosefs/guarracino/HPRCv2/PHR_III/pq-classification/chm13.centromeres.approximate.bed` | CHM13 centromere coordinates (from Julian's active_arrays) | 1 |
| `/moosefs/guarracino/HPRCv2/PHR_III/liftoff_genes_hprc_r2_v1.0.index.csv` | Index mapping 462 haplotype-specific Liftoff GFF3 files | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprc_annotations/*.gff3.gz` | 464 haplotype-specific + CHM13 gene annotations (Liftoff GFF3; includes HG002 from JHU Liftoff v0.6) | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/repeat_masker_bed_hprc_r2_v1.0.index.csv` | Index mapping RepeatMasker BED files | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprc_repeatmasker/*.RepeatMasker.bed.gz` | 464 haplotype-specific + CHM13 = 465 repeat annotations (includes HG002 converted from bigBed to 10-col BED with PanSN naming) | 8 |
| `/moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv` | HPRCv2 sample metadata (sample, superpopulation, sequencing info) | 10 |

### Verification scripts

| Path | Description |
|------|-------------|
| `/moosefs/guarracino/HPRCv2/PHR_III/verify_report_numbers.sh` | Bash: verifies 102 numerical claims against source TSV/FAI files |
| `/moosefs/guarracino/HPRCv2/PHR_III/verify_report_numbers.R` | R: verifies MDS variance, distance matrix, TAR1/gene counts from RDS files |
| `/moosefs/guarracino/HPRCv2/PHR_III/verify_literature_claims.sh` | Bash: extracts exact quotes from review papers for literature-sourced claims |
| `/moosefs/guarracino/HPRCv2/scripts/verify_rpe1_results.py` | Python: extracts per-sample global summary + per-community W/B ratios across all human + RPE-1 datasets (§12.1, §12.14, §12.15) |

### Intermediate data files

| Path | Description | Section |
|------|-------------|---------|
| `/moosefs/guarracino/HPRCv2/PHR_III/assembly-vs-chm13/*.paf` | 466 assembly-to-CHM13 PAF mapping files | 1 |
| `/moosefs/guarracino/HPRCv2/PHR_III/pq-classification/contig_classifications.tsv` | 12,649 classified contigs (+ header) | 1→2 |
| `/moosefs/guarracino/HPRCv2/PHR_III/pq-classification/p_arms.bed` | CHM13 p-arm coordinates in BED format | 1 |
| `/moosefs/guarracino/HPRCv2/PHR_III/pq-classification/q_arms.bed` | CHM13 q-arm coordinates in BED format | 1 |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz` | 18,827 subtelomeric flank sequences (2.5 GB) | 2→3 |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz.fai` | FASTA index (18,827 lines) | 2→4 |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz.gzi` | bgzip index | 2 |
| `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95/*.paf.gz` | 18,827 all-vs-all PAF alignment files (88 GB total) | 3→4 |
| `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.paf.list` | List of all PAF file paths for impg indexing | 3→4 |
| `/scratch/all-vs-all.p95.impg` | impg index of all-vs-all subtelomeric alignments (built at runtime on scratch) | 4 |
| `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv` | Inter-chromosomal region calls (18,827 data rows + header; 1 chimeric contig excluded) | 4→5,8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.id95.len.no_acrocentric.tsv` | Acrocentric-excluded version (chr13/14/15/21/22 removed) for §12.9 confound control | 12.9 |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz` | 15,668 PHR sequences with inter-chromosomal signal (529 MB) | 4→5 |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.fai` | FASTA index (15,668 lines) | 4→6 |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.gzi` | bgzip index | 4 |
| `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.*.smooth.final.og` | Pangenome graph (odgi format) | 5 |
| `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.*.smooth.final.gfa.zst` | Pangenome graph (compressed GFA) | 5 |
| `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.*.smooth.final.similarity.tsv.gz` | Jaccard pairwise similarity matrix (10.8 GB compressed) | 5→6 |
| `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.*.smooth.final.lay.tsv` | Graph layout coordinates | 5 |
| `/moosefs/guarracino/HPRCv2/PHR_III/annotations/subtelomeric_annotations.1Mb.rds` | R object: gene + TAR1 annotations intersected with PHR regions | 8→9,10 |

### Output data files

**Similarity and community detection:**

| Path | Description | Section |
|------|-------------|---------|
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.dist_matrix.rds` | 15,668 × 15,668 Jaccard distance matrix (R object) | 6 |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.full_mds.rds` | MDS coordinates, k=5 dimensions (R object) | 6 |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.umap.rds` | UMAP coordinates, 3 components (R object) | 6 |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv` | Arm-level Leiden community assignments (41 arms, 15 communities) | 7 |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-upgma-k12.assignments.tsv` | Arm-level UPGMA assignments (41 arms, 12 communities) | 7 |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv` | Sequence-level Leiden assignments (15,668 sequences, 50 communities; includes sample, arm, superpopulation) | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.arm-matrix.tsv` | Arm × community count matrix | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.summary.tsv` | Per-community summary statistics (50 rows) | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-upgma-k145.arm-matrix.tsv` | UPGMA arm × community matrix | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-upgma-k145.summary.tsv` | UPGMA per-community summary (145 rows) | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden.k-scan.tsv` | Leiden k-NN × resolution parameter scan results | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv` | 41 × 41 arm-level Jaccard distance matrix (exported from RDS for Hi-C/Dip-C scripts) | 6→13,14 |
| `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-upgma.k-scan.tsv` | UPGMA k-scan results | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/comparison/arm_vs_sequence_nesting.tsv` | Arm→sequence community fragmentation (how 15 arm communities split into 50 seq communities) | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/comparison/arm_vs_sequence_concordance.tsv` | Per-sequence arm vs sequence community assignment concordance (15,668 rows) | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/comparison/arm_vs_sequence_summary.tsv` | Summary statistics: N communities, silhouette, multi-arm counts | 8 |

**Enrichment analysis:**

| Path | Description | Section |
|------|-------------|---------|
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_summary_table.tsv` | Per-community summary: arms, gene counts, biotypes, interpretation | 9 |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_enrichment.tsv` | Gene × community presence matrix | 9 |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_specific_genes.tsv` | Genes found in exactly one community | 9 |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/cross_community_genes.tsv` | Genes spanning multiple communities | 9 |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_biotypes.tsv` | Biotype distribution per community | 9 |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_families.tsv` | Gene family counts per community | 9 |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_enrichment_fisher.tsv` | Fisher's exact test results for gene enrichment | 9 |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_tar1_summary.tsv` | TAR1 prevalence per community | 9 |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_tar1_by_arm.tsv` | TAR1 prevalence per arm within each community | 8, 9 |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_within_sharing.tsv` | Within-community gene sharing statistics | 9 |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_singleton_anomalies.tsv` | Singleton/doubleton QC (27 singletons, 6 doubletons) | 9 |

**Heterogeneity analysis:**

| Path | Description | Section |
|------|-------------|---------|
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/within_community_arm_separation.tsv` | Silhouette scores and separation ratios per community (12 rows) | 10.1 |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_affinity_sequences.tsv` | Per-sequence cross-arm vs self-arm classification (2,271 cross-arm) | 10.2 |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_superpop_enrichment.tsv` | Fisher's exact test for superpopulation bias (9 significant rows) | 10.3 |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_type_discordance.tsv` | Per-individual subtelomeric type discordance (column `discordant`: one haplotype self-arm, one cross-arm) | 10.4 |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_region_length_comparison.tsv` | Wilcoxon rank-sum test for cross-arm vs self-arm length differences | 10.5 |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/gene_conversion_scores.tsv` | Gene repertoire replacement scores per arm pair | 10.7 |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_gene_content.tsv` | Gene content comparison of cross-arm vs self-arm sequences | 10.7 |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_tar1_comparison.tsv` | TAR1 prevalence in cross-arm vs self-arm sequences | 10.8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/allele_vs_paralog_distance.tsv` | Allele vs paralog distance per community (10 communities + overall) | 10.1 |
| `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/fst_superpop_matrix.tsv` | Pairwise Fst matrix (5 superpopulations) | 10.4 |
| `/moosefs/guarracino/HPRCv2/PHR_III/ttaggg_boundary_enrichment.tsv` | TTAGGG island boundary enrichment tests | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/parm_qarm_3d_comparison.tsv` | p-arm vs q-arm 3D enrichment comparison | 13.6 |
| `/moosefs/guarracino/HPRCv2/PHR_III/telomere_length_by_community.tsv` | Terminal telomere tract length per community (15 rows) | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/ttaggg_island_length_distribution.tsv` | TTAGGG island length histogram (12 bins, 50–1000 bp) | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/ttaggg_island_motif_composition.tsv` | TTAGGG island motif composition (4 motifs, hexamer counts + island classification) | 8 |
| `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.id95.len.per_window.tsv` | Per-window inter-chromosomal signal (431,910 rows: n_chrs, n_arms, dist_from_telomere) | 10.11 |
| `/moosefs/guarracino/HPRCv2/PHR_III/rm_ttaggg_subtelomeric.tsv` | RepeatMasker telomere annotations in subtelomeric flanks (29,274 rows) | 10.11 |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/two_domain_test.tsv` | Per-arm Spearman gradient results (48 arms) | 10.11 |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/two_domain_per_sequence.tsv` | Per-sequence within-sequence gradient (13,763 sequences) | 10.11 |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/two_domain_changepoint.tsv` | Per-arm breakpoint analysis (41 arms) | 10.11 |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/two_domain_binned_means.tsv` | Binned mean n_chrs by dist_from_telomere (for plotting) | 10.11 |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/polymorphic_arms_summary.tsv` | Per-arm: n_communities, largest cluster % | 7.2 |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/polymorphic_arms_subgroups.tsv` | Per-arm per-community: size, region length, n_chrs, population | 7.2 |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/d4z4_perwindow_signal.tsv` | Binned per-window signal for chr4_q/chr10_q (D4Z4 test 1) | 7.2 |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/d4z4_dux4l_by_community.tsv` | DUX4L count per sequence with community assignment (D4Z4 test 2) | 7.2 |

---

## References

- Ambrosini A, Paul S, Hu S, Riethman H (2007). Human subtelomeric duplicon structure and organization. *Genome Biology* 8:R151. — Duplicon module classification; bimodal identity distribution (91%/98% peaks); subterminal vs subtelomere-only block distinction; one-copy (TTAGGG)n-adjacent regions at chr7_q, chr8_q, chr11_q, chr12_q, chr18_q, chrX_p/chrY_p.
- Brown WRA, MacKinnon PJ, Villasanté A, et al. (1990). Structure and polymorphism of human telomere-associated DNA. *Cell* 63:119–132. — Original characterization of TAR1 (Telomere-Associated Repeat 1) as a subtelomeric repeat element.
- Flint J, Bates GP, Clark K, et al. (1997). Sequence comparison of human and yeast telomeres identifies structurally distinct subtelomeric domains. *Human Molecular Genetics* 6:1305–1313. — Proposed two-domain subtelomeric model: distal domain (short blocks shared with many chromosome ends) and proximal domain (longer blocks shared with few ends), separated by degenerate (TTAGGG)n tracts; characterized at chr4p, chr16p, chr22q.
- Francis BA, et al. (2025). Complete genome assemblies of two mouse subspecies reveal structural diversity of telomeres and centromeres. *Nature Genetics* 57:2852–2862. — First T2T assemblies for C57BL/6J and CAST/EiJ; B6 adds 208 Mb, CAST 247 Mb vs GRCm39; conserved L1-LINE + TLC subtelomeric architecture in B6; heterogeneous CAST subtelomeres.
- Gershman A, Sauria MEG, Guitart X, et al. (2022). Epigenetic patterns in a complete human genome. *Science* 376:eabj5089. — ENCODE CTCF ChIP-seq realignment to T2T-CHM13; CTCF enrichment at TAR loci.
- Gonzalez IL & Sylvester JE (1995). Complete sequence of the 43-kb human ribosomal DNA repeat: analysis of the intergenic spacer. *Genomics* 27:320–328. — rDNA sequence variants spread across all five acrocentric chromosomes via inter-chromosomal gene conversion.
- Lemmers RJLF, van der Vliet PJ, Klooster R, et al. (2010). A unifying genetic model for facioscapulohumeral muscular dystrophy. *Science* 329:1650–1653. — Permissive 4qA haplotype produces stable DUX4 mRNA via polyadenylation signal.
- Ottaviani A, Rival-Gervier S, Boussouar A, et al. (2009). The D4Z4 macrosatellite repeat acts as a CTCF and A-type lamins-dependent insulator in facio-scapulo-humeral dystrophy. *PLoS Genetics* 5:e1000394. — CTCF binding within D4Z4 repeat units; lamin A/C-dependent insulator function.
- Linardopoulou EV, Williams EM, Fan Y, et al. (2005). Human subtelomeres are hot spots of interchromosomal recombination and segmental duplication. *Nature* 437:94–100. — Paralogy map of 41 blocks across 33 subtelomeres; translocation-based model of segmental duplication; subtelomeric interchromosomal duplication/transfer rate >60-fold higher than point mutation or retrotransposon insertion rates. Block 5 (OR4F21) shared by C11 arms chr1_p/chr5_q/chr6_q/chr8_p; blocks 35–37 link chr1_p–chr8_p; block 38 links chr6_q–chr8_p.
- Masny PS, Bengtsson U, Chung SA, et al. (2004). Localization of 4q35.2 to the nuclear periphery: is FSHD a nuclear envelope disease? *Human Molecular Genetics* 13:1857–1871. — 4q35.2 peripheral localization via lamin A/C; sequences proximal to D4Z4 mediate nuclear envelope interaction.
- Mefford HC & Trask BJ (2002). The complex structure and dynamic evolution of human subtelomeres. *Nature Reviews Genetics* 3:91–102. — Foundational review; f7501 block distribution (African-enriched at chr7p/chr16_q); ~20% chr4_q/chr10_q translocation prevalence in a Dutch population (citing van Deutekom et al. 1996); patchwork duplicon architecture; somatic mosaicism for chr4_q/chr10_q translocations (citing van Overveld et al. 2000); telomere position effect and subtelomeres as buffers against TPE.
- Patel L, Kang R, Rosenberg SC, et al. (2019). Dynamic reorganization of the genome shapes the recombination landscape in meiotic prophase. *Nature Structural & Molecular Biology* 26:164–174. — Mouse meiotic Hi-C; chromosome end clustering.
- Riethman H, Ambrosini A, Castaneda C, et al. (2004). Mapping and initial analysis of human subtelomeric sequence assemblies. *Genome Research* 14:18–28. — Subtelomeric assembly and characterization; 80% of the most distal 100 kb consists of shared duplicated blocks.
- Rouyer F, Simmler MC, Johnsson C, et al. (1986). A gradient of sex linkage in the pseudoautosomal region of the human sex chromosomes. *Nature* 319:291–295. — Obligate crossover in PAR1 during male meiosis; establishes PAR1 as a region of obligate recombination distinct from subtelomeric repeat-mediated exchange.
- Stout K, van der Maarel S, Frants RR, et al. (1999). Somatic pairing between subtelomeric chromosome regions: implications for human genetic disease? *Chromosome Research* 7:323–329. — Interphase chr4_q/chr10_q clustering by FISH; nucleolar association of acrocentric short arms.
- Tan L, Xing D, Chang CH, et al. (2018). Three-dimensional genome structures of single diploid human cells. *Science* 361:924–928. — Dip-C method; GM12878 and PBMC single-cell 3D structures.
- Traag VA, Waltman L, van Eck NJ (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports* 9:5233. — Leiden community detection algorithm; guarantees connected communities unlike Louvain; used for both arm-level and sequence-level community detection (§6).
- Lalli JL, Bortvin AN, McCoy RC, Werling DM (2025). A T2T-CHM13 recombination map and globally diverse haplotype reference panel improves phasing and imputation. *bioRxiv* 2025.02.24.639687. — T2T-CHM13 recombination maps at variant-level resolution; population-specific maps for 26 subpopulations.
- Xu H, et al. (2025). Three-dimensional genome organization of human sperm at single-cell resolution. *Nature Communications* 16:3805. — Sperm single-cell 3D genome structures (20 cells); Rabl-like configuration persists in post-meiotic sperm.
- Zuo W, Chen G, Gao Z, et al. (2021). Stage-resolved Hi-C analyses reveal meiotic chromosome organizational features influencing homolog alignment. *Nature Communications* 12:5827. — Meiotic chromosome end alignment extends up to 20% of chromosome length; LINC complex force transmission.
