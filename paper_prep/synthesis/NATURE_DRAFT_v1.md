---
title: Concerted evolution and unorthodox recombination of human subtelomeres
target_journal: Nature
target_format: Article
abstract_words: 218
main_text_words: 3905
generated: 2026-05-17
inputs:
  outline: paper_prep/synthesis/NATURE_DRAFT_OUTLINE.md (sha: f96d52b294bdf295014751426c803d084cd211cc)
  pptx_outline: paper_prep/synthesis/PPTX_OUTLINE.md (sha: f2f5bba7d00ab0e09dd014873cd3010ef871114a)
  references: paper_prep/synthesis/REFERENCES_v3.bib (sha: 81c6766bbf19af68b771156ada14f17e3fa25535)
---

# Concerted evolution and unorthodox recombination of human subtelomeres

## Abstract

Human subtelomeres are among the most dynamic and structurally complex regions of the genome, yet interchromosomal sequence relationships have resisted systematic analysis owing to assembly incompleteness and per-chromosome alignment frames [@MeffordTrask2002; @Riethman2004; @Linardopoulou2005; @BaileyEichler2006]. Using 465 near-complete haplotype assemblies from 233 individuals of the Human Pangenome Reference Consortium v2 plus CHM13v2.0 (466 haplotype-equivalent units), we build an implicit pangenome graph: the all-vs-all wfmash PAF set, queried by IMPG transitive closure, without chromosomal partitioning [@hprc_hprcv2_2025; @Garrison2024pggb; @pangenome_graphs_impg_GarrisonGuarracino2023]. The graph maps 18,827 telomere-anchored 500 kb flanks to 15,668 pseudohomologous regions (PHRs; median 105 kb) on 41 of 48 chromosome arms, in scale comparable to the canonical PAR2. A neighbour-joining tree of arm-level Jaccard distances recovers PAR1, PAR2, the acrocentric short arms, a 10p/18p clade, a tight {22q, 21q, 19q, 1q, 13q, 17q} q-arm clade and 4q/10q DUX4 with wide copy-number diversity; all six clades map one-to-one to a Leiden partition of 15 arm-level and 50 sequence-level communities. Bulk and single-cell Hi-C, Pore-C, CiFi and Dip-C across six individuals, 20-cell sperm scHi-C and mouse meiotic Hi-C peaking at zygotene (Spearman ρ = 0.715) tie sequence similarity to nuclear-envelope proximity through the meiotic bouquet. A T2T 3-generation pedigree resolves 538 inter-chromosomal patches with 92% within Leiden communities, including 133 gene-conversion-like and 16 crossover-like events. Human subtelomeres are unified by ongoing inter-chromosomal recombination and concerted evolution.

## Main text

The chromosome ends of the human genome were the first regions in which inter-chromosomal sequence exchange was identified [@Brown1990; @Wilkie1991; @Trask1991; @Trask1998]. The pseudoautosomal regions PAR1 (about 2.6 Mb at Xp/Yp) and PAR2 (around 334 kb at Xq/Yq) sustain obligate meiotic crossover between the sex chromosomes [@Rouyer1986; @sexchrompars_acquaviva2020; @sexchrompars_bellott2024]; the acrocentric rDNA-bearing short arms 13p, 14p, 15p, 21p and 22p share large blocks of identity that recombine in pre-meiotic nucleoli [@acrocentric_Altemose2022; @acrocentric_Guarracino2025ape]; the D4Z4 macrosatellite of 4q35 is found in degenerate copies on 10q26 and underlies facioscapulohumeral dystrophy [@dux4_d4z4_fshd_lemmers2010worldwide; @dux4_d4z4_fshd_lemmers2007]. Cytogenetic FISH and BAC-walking efforts in the 1990s and 2000s mapped multi-chromosomal duplicons across roughly half of the 48 chromosome arms [@Mefford2001; @MeffordTrask2002; @Linardopoulou2005; @Riethman2001; @Riethman2004; @Riethman2008; @Flint1997; @Ambrosini2007], but population-scale extension stalled because no reference assembly was complete to the telomeres [@Nurk2022; @Logsdon2021] and because per-chromosome alignment frames hid trans-chromosomal sequence sharing. The Human Pangenome Reference Consortium v2 release delivers near-T2T haplotype-resolved assemblies for 233 individuals from five superpopulations [@hprc_hprcv2_2025; @Liao2023], removing both constraints. We use these data to revisit subtelomeric architecture without chromosomal partitioning, asking three quantitative questions: how extensive is inter-chromosomal sequence sharing across the unbiased population, does that sharing predict three-dimensional nuclear proximity, and can the events that build the population structure be observed directly in human pedigrees? This paper, the subtelomere companion to the HPRC v2 reference assembly publication [@hprc_hprcv2_2025], shows that the answers are yes, yes, and yes.

We treat every haplotype as its own reference. From each of the 233 individuals we used both haplotype assemblies (465 haplotypes total; one HPRC haplotype assembly was excluded for quality) plus CHM13v2.0, then extracted 18,827 telomere-anchored 500 kb flanks across 48 chromosome arms (Extended Data Fig. 1a, 1b) and ran wfmash all-vs-all at 95% minimum identity [@Guarracino2023] without restricting alignment to homologous chromosomes. The resulting PAF set is the implicit pangenome graph: nodes are the 18,827 flanks, edges are wfmash alignments, and the canonical query is the IMPG transitive closure [@pangenome_graphs_impg_GarrisonGuarracino2023; @pangenome_graphs_impg_GuarracinoHeumos2022; @pangenome_graphs_impg_IMPG2023; @pangenome_graphs_impg_Hickey2024]. Because pairwise alignment of 18,827 flanks is C(18,827, 2) = 177 million pairs, full all-to-all alignment is computationally infeasible; wfmash k-mer prefiltering instead evaluates 11.6% of pairs (rounded to 12% throughout). The realised sampling rate sits 230x above the Erdős and Rényi connectivity threshold p* = log(n)/n = 5.2 x 10^-4 for n = 18,827, so transitive closure reaches virtually every subtelomere in the dataset (Methods). The all-vs-all alignment is the substrate; PHR detection is the call. We define a pseudohomologous region as a flank window in which IMPG transitive closure recovers alignment to at least 5 segments on at least 2 different chromosomes at identity ≥ 95%, with total aligned length ≥ 3 kb [@Garrison2018; @Garrison2024pggb]. Applying this filter to the 18,827 flanks yields 15,668 PHRs (83.2% of flanks) on 41 of 48 chromosome arms. The 7 remaining arms (chr5p, chr6q, chr7q, chr12q, chr14q, chr20p, chr20q) carry no detectable inter-chromosomal homology under the same filter; we use these as built-in negative controls throughout the paper. Sequence-level Jaccard similarity between PHRs was computed downstream from a PGGB graph of the 15,668 sequences [@Garrison2024pggb] with odgi similarity [@Guarracino2023], yielding a 15,668 x 15,668 Jaccard matrix that is the input to all subsequent cladistics, community detection and 3D-validation analyses.

Genome-wide stacked identity heatmaps (Fig. 1a) show telomere-anchored high-identity blocks at nearly every chromosome end. Each row is one of the 465 HPRC v2 haplotypes plus CHM13 (466 total), each column is one of 24 chromosomes scanned in 100 kb windows, and each cell is the maximum identity of that window to any matching window on any other chromosome. Inter-chromosomal blocks extend from the very tip inward to between 30 kb and several hundred kilobases. Stacking by chromosome contributor (Fig. 1b) makes the same landscape quantitative: the number of distinct partner chromosomes per flank scales with the depth of sharing, and the called PHR-BED tracks overlay the visible blocks one-to-one. The detected PHRs span tens to hundreds of kilobases (median 105 kb, mean 144 kb, range 5 kb to 500 kb; Extended Data Fig. 1c), with the upper bound determined by the 500 kb flank window rather than by biology. The chr18q chimera control (Extended Data Fig. 1d) removes one false positive (15,669 to 15,668). To put scale in cytogenetic terms, the median PHR (105 kb) is 31% of the length of PAR2 (334 kb), and the mean PHR (144 kb) is 43% of PAR2: the largest known fully homogenised inter-chromosomal region in the human genome [@sexchrompars_bellott2024] is now flanked, on the population scale, by 15,668 examples of structurally comparable architecture [@Bailey2002; @BaileyEichler2006; @RuizHerrera2008; @Vollger2023; @concerted_evolution_nahr_Vollger2023; @Stong2014]. The 41 signal-bearing arms partition into 15 communities at the arm level (Fig. 1c) and into 50 communities at the sequence level (Extended Data Fig. 2a, modularity 0.97; Methods). The 7 silent arms produce zero PHRs under the 95% identity filter and provide the S_all negative control invoked in the 3D analysis below.

A neighbour-joining tree built on the 41 x 41 arm-level Jaccard distance matrix recovers six monophyletic clades that match every known case of inter-chromosomal subtelomere homology in the human genome (paper_prep/figures/nj_tree_arms/; Fig. 1c shows the same matrix as a heatmap with UPGMA dendrogram). PAR1 (Xp, Yp) forms one clade; PAR2 (Xq, Yq) forms a second [@sexchrompars_acquaviva2020]. The acrocentric short arms 13p, 14p, 15p, 21p and 22p form a third clade consistent with their shared rDNA arrays and Robertsonian recombination [@acrocentric_Altemose2022; @acrocentric_Guarracino2025ape]. A 10p/18p clade reproduces the high-identity pair first reported by Linardopoulou and colleagues in 2005 [@Linardopoulou2005]. A tight q-arm clade comprising 22q, 21q, 19q, 1q, 13q and 17q is, to our knowledge, the first cladistic recovery of this composite group at population scale. Finally, 4q and 10q pair through the D4Z4 macrosatellite [@dux4_d4z4_fshd_lemmers2010worldwide; @dux4_d4z4_fshd_lemmers2007; @Cabianca2012], with wide copy-number diversity across the 465 haplotypes (DUX4 family copy counts in Extended Data Fig. 4c, range 0 to 22). The placement of each clade is robust to algorithm choice. Leiden community detection on a kernel-weighted version of the same distance matrix yields 15 arm-level communities (mean silhouette 0.347, k = 15), and the six abstract-named clades map one-to-one to communities C15 (PAR1), C14 (PAR2), C7 (acrocentric p-arms), C2 (10p/18p), C6 ({22q, 21q, 19q, 1q, 13q, 17q}) and C1 (4q/10q DUX4). UPGMA at k = 14 agrees with Leiden on 14 of 15 communities, and a 1,000-replicate perturbation bootstrap (Gaussian noise on off-diagonal distances at sigma = 25% of the IQR) puts the MRCA support of every named clade at 100%. The deeper backbone of the tree, by contrast, has bootstrap support of 32 to 90% across internal edges, mirroring the well-known fact that subtelomere-internal duplicon order is not strictly monophyletic at the chromosome level [@Skaletsky2003; @Rudd2009; @Mefford2001]. The cladistic signal in human subtelomeres is robust to choice of algorithm, robust to perturbation of the input distances, and partitions the genome into a small set of named, interpretable groups.

The 41 signal-bearing arms partition into three architectural categories on the joint axes of cross-arm sequence rate and per-community silhouette (Fig. 1d). Four arms are homogeneous: each is the sole member of its arm-level community (C8, C9, C10, C13) with allele-vs-paralog distance ratios near unity. Twenty-eight arms are polymorphic: they sit in multi-arm communities but retain arm identity, in the sense that within-community sequences cluster by source arm and the allele is closer than the paralog at the median (Wilcoxon paired test, allele-minus-paralog distance significantly negative). The remaining nine arms are fully interchangeable: the five acrocentric p-arms (C7), Xq with Yq (C14, PAR2) and Xp with Yp (C15, PAR1). All three interchangeable communities have negative mean silhouette, paralog-closer-than-allele distance and sample-of-origin labels that do not separate from chromosome-of-origin labels. The cross-arm sequence rate peaks at chrX_q 99.7%, chr21_p 94.0% and chr11_p 74.1% as the highest autosomal value (Extended Data Fig. 2c). These three categories refine the 8/34/7 architectural skeleton suggested by FISH-era inspection [@Linardopoulou2005; @Stong2014; @Eichler2001] into a 4/28/9 split that respects negative-silhouette membership.

The 15 arm-level communities are not internally homogeneous. Across 5,946 paired distances from the nine multi-arm communities, the within-haplotype-pair allelic distance is shorter than the corresponding cross-paralog distance in 8 of 9 communities, with combined p < 10^-300 (Fig. 2a) [@Flint1997; @MeffordTrask2002]. The single exception is the acrocentric p-arms (C7), where 70.5% of pairs have the paralog closer than the allele (p = 2.0 x 10^-7): quantitative confirmation that the rDNA-bearing short arms are homogenised more deeply than within-individual allelic variation [@acrocentric_Altemose2022; @acrocentric_rdna_robertsonian_bandyopadhyay2001]. Within communities C1, C2, C3, C5, C6, C7, C11 and C12, the within-community Jaccard distance distribution is bimodal (Extended Data Fig. 2b): an allelic mode at low distance and a paralog mode at higher distance separate cleanly. The two-domain Flint and Mefford model [@Flint1997; @MeffordTrask2002; @Mefford2001] of a proximal unique segment and a distal duplicon-rich segment generalises to the pangenome. In 39 of 48 arms, the number of distinct chromosome contributors per flank window decreases monotonically with distance from the telomere (Spearman ρ < 0; Fig. 2b top), and in 39 of 41 signal-bearing arms a piecewise linear model with a single breakpoint outperforms a linear model on Akaike's information criterion (Fig. 2b bottom). Sixteen of nineteen arms that contain an internal (TTAGGG)n island also have an ITS within 50 kb of the inferred breakpoint, consistent with telomeric repeat fragments demarcating the proximal-to-distal transition [@subtelstruct_NergadzeITS2007; @subtelstruct_Nergadze2007; @subtelstruct_NergadzeITSReview2007; @Ambrosini2007]. TAR1 prevalence (Extended Data Fig. 3a) is correlated with arm architecture: PAR1 arms are TAR1-free (chrXp 0.3%, chrYp 1.1%), acrocentric p-arms sit at 73 to 79%, and the remaining autosomal arms saturate above 99%. Inter-chromosomal exchange leaves a population-genetic signature. A 2 x 5 Fisher exact test for superpopulation composition of cross-arm vs self-arm sequences is significant after Benjamini and Hochberg correction in 10 of 19 testable arms (Fig. 2c, left). Hudson pairwise F_ST [@subtel_popgen_hudson1992; @subtel_popgen_weir1984; @subtel_popgen_lewontin1972] computed on cross-arm sequences yields F_ST values of 0.10 to 0.15 between AFR and each of AMR, EAS, EUR and SAS, and 0.02 to 0.04 within the non-AFR set (Fig. 2c, right). A UPGMA tree of the F_ST matrix recovers an out-of-Africa topology in which AFR splits first and AMR, EAS, EUR and SAS form a tight non-AFR clade (Fig. 2d) [@Bergstrom2020; @subtel_popgen_anderson2008; @subtel_popgen_bhatia2013; @subtel_popgen_rosenberg2002; @subtel_popgen_levysakin2019; @subtel_popgen_1000g2010].

Sequence-defined communities are physical. We assembled 14 inter-arm tests of three-dimensional proximity across six HPRC v2 individuals and the CHM13 cell line: bulk Hi-C in CHM13 and five HPRC samples, Pore-C and CiFi in HG002, Dip-C in GM12878 [@Tan2018], 20-cell single-sperm scHi-C [@Xu2025] and four mouse meiotic Hi-C stages [@Zuo2021; @Cechova2025]. In every test the within-community contact frequency exceeds the between-community frequency: B/W ratios for Hi-C and W/B ratios for distance-based scHi-C all sit left of unity, with a range of 0.020 to 0.93 across the 14 tests (Fig. 3b forest plot). The strongest single measurement is HG002 Hi-C at 50 kb, B/W = 0.027, p = 4.0 x 10^-66 (Mann-Whitney test on bootstrap-resampled within and between distributions). CHM13 Hi-C yields B/W = 0.071, p = 6.0 x 10^-18; HG002 Pore-C yields B/W = 0.056, p = 3.9 x 10^-85 (Fig. 3a) [@Ulahannan2019; @hic3d_cifi2025]. Mantel correlation between the 41 x 41 arm-level sequence similarity matrix and the corresponding 41 x 41 Hi-C contact matrix is Spearman ρ = 0.66 for both CHM13 and HG002 (10,000 row and column permutations). Per-individual sequence-pair Spearman ρ rises to 0.83 in the lowest-coverage samples where the long-range signal-to-noise is best. The signal is not driven by acrocentric, nucleolar or PAR contacts. Across five mcool resolutions (5, 10, 20, 50, 100 kb) and five exclusion sets (no acrocentric p-arms, no sex chromosomes, no acrocentric p plus sex, no all-acrocentric plus sex, no strongest community), the Mantel correlation strengthens when subset confounds are excluded: HG002 0.66 to 0.80, HG02148 0.15 to 0.21, CHM13 0.66 to 0.85, NA19036 0.27 to 0.49 (Extended Data Fig. 5a, 5b) [@hic3d_dixon2012; @hic3d_imakaev2012; @hic3d_alavattam2019; @hic3d_wolff2018; @hic3d_deshpande2022]. Within-vs-between observed-over-expected enrichment across all trans-arm pairs ranges from 5.9x (HG02559) to 18.4x (HG002 CiFi) and is robust across Hi-C, CiFi and Pore-C platforms (8 of 8 tests; Extended Data Fig. 5c). Across 11 datasets and 15 communities, the within-community contact is enriched over the random-label null in every cell of the 15 x 11 reproducibility heatmap (Extended Data Fig. 5d), with Benjamini and Hochberg q < 0.05 in the majority of cells and q < 0.001 in the strongest.

The 3D signal is not a multi-mapping artefact. We computed B/W ratios separately for the duplicated PHR windows themselves and for the unique-sequence 100 kb regions immediately centromere-ward of the PHR boundaries. The flanking unique-sequence regions, which by construction carry no inter-chromosomal duplication, produce a stronger within-community signal than the PHRs themselves (Fig. 3d). In HG002 Hi-C the PHR B/W of 0.027 falls to flanking B/W of 0.0031, a 9-fold strengthening. Multi-mapping of reads to identical paralogous sequence cannot inflate the within-community signal in regions that contain no paralogous sequence. The 7 silent arms provide the complementary negative control. Pooled as the pseudo-community S_all, they sit systematically farther apart in 3D space: in GM12878 Dip-C, 16 of 16 C-community cells have W/B < 1 while 0 of 16 S_all cells do; in 20-cell sperm scHi-C [@Xu2025], 20 of 20 C-cells have W/B < 1 while only 1 of 20 S_all cells does (Fig. 3c). S_all is 11% farther in GM12878 and 40% farther in sperm than the average C community. The Dip-C radial-position inset (Fig. 3d) shows the same effect on a different axis: flanking unique-sequence particles are more nuclear-interior than non-flanking terminal particles (radial 0.504 vs 0.556, p = 1.6 x 10^-35). A mechanism is available. Meiotic prophase I telomeres are tethered to the nuclear envelope by the MAJIN, TERB2 and TERB1 complex, which drives the bouquet stage of zygotene chromosome organisation [@bouquet_KotaSUN1MAJIN2020; @bouquet_Scherthan2001; @bouquet_Scherthan2003; @bouquet_ShibuyaRPMs2015; @bouquet_ChikashigeTelomere1994; @bouquet_HarperBouquet2004; @bouquet_HornKASH52013; @bouquet_DingSUN12007; @bouquet_MorimotoKASH2012; @bouquet_ZicklerKleckner1999; @ZicklerKleckner1998; @ZicklerKleckner2015]. Telomere clustering during meiosis pre-positions subtelomeric regions in physical proximity at every meiosis, providing structural opportunity for ectopic exchange between non-homologous chromosomes whose terminal sequences happen to share identity. The D4Z4 array on 4q35 binds CTCF and is tethered to the nuclear lamina [@Ottaviani2009; @OttavianiGilson2008; @Masny2004; @Cabianca2012], which fixes the C1 (4q, 10q) co-peripheralisation observed in Dip-C (radial 0.732, Extended Data Fig. 8b).

Indirect inference of inter-chromosomal exchange can be replaced by direct observation. The WashU T2T pedigree comprises four haplotype-resolved T2T-quality genomes from three generations of one family: grandmother PAN010, grandfather PAN011, mother PAN027 and granddaughter PAN028 [@Cechova2025]. For each pair of parent and child haplotypes, we ran `odgi untangle nth-best=1` per flank to call inter-chromosomal patches of recent recombination, then filtered to high-quality patches (minimum patch and alignment score 0.95) and restricted to patches whose query and target arms belong to the same HPRC v2 Leiden community. The filter yields 538 high-quality inter-chromosomal patches. 494 of 538 (92%) sit within a Leiden community: the population-scale partition almost completely predicts where new inter-chromosomal recombination is found in a single family (Fig. 4a). The 494 within-community patches break down into five mutually exclusive patterns. 229 are acros_like: 5 or more inter-chromosomal patches drawn from 3 or more source chromosomes within a single flank, the classical multi-paralog non-allelic homologous recombination signature [@concerted_evolution_nahr_SamonteEichler2002; @concerted_evolution_nahr_Eichler2001; @concerted_evolution_nahr_Hastings2009; @concerted_evolution_nahr_Myers2010; @Sharp2006]. 133 are gene-conversion-like sandwiches in which a query arm flank reads chrN:hX, chrM:hY, chrN:hX, with the central chrM:hY block aligning at 1.000 identity to the source paralog; roughly 90% of these are in C7, the acrocentric p-arm community. 16 are crossover-like: a single reciprocal exchange in which the query haplotype switches its source chromosome at a discrete breakpoint and stays on the new chromosome to the telomere. The largest crossover-like event spans 27.97 kb on the PAN028 maternal chr3q. 115 are sandwich_same_hap (a within-haplotype interchange of patches), and 1 is complex. The inheritance pattern is observed directly: PAN027 inherits her maternal haplotype from PAN010 and her paternal haplotype from PAN011, and PAN028 inherits her maternal haplotype from PAN027 [@Cechova2025; @StankiewiczLupski2002; @StankiewiczLupski2010]. Twelve of the sixteen crossover-like events are in PAN028, confirming that meiotic-resolution inter-chromosomal breakpoints transmit across generations. The CEPH1463 4-generation Platinum Pedigree [@Porubsky2025; @acrocentric_Porubsky2025denovo] provides a stricter test. We required that a parent x chromosome-pair feature be independently called by both hifiasm and verkko assembly pipelines and that the involved arms belong to the same Leiden community. 11 features pass. They include chr10 x chr18 (C2, the Linardopoulou pair) independently in NA12877 paternal and NA12878 maternal, chr12 x chr9 (C5) in NA12889 and NA12890 grandparents, chr6 x chr9 (C5) independently in NA12877 and NA12878, and a chr19 x chr22 feature transmitted via NA12878 (Fig. 4b). Every cross-assembler-validated event in the platinum pedigree sits within an HPRC v2 Leiden community: a second, fully independent family confirms that the partition predicts where new inter-chromosomal exchange is generated.

The methodology generalises to a single diploid genome and to a non-human mammal. We applied the full pipeline to the 46-arm RPE-1 retinal pigment epithelial cell line, the only diploid human cell line with a publicly available T2T assembly [@Francis2025]. Leiden community detection on the RPE-1 distance matrix recovered 37 self-discovered communities, including a Leiden C2 = {chr10_q, chrX_q} community (Fig. 4c). The two arms also carry elevated asynchronous CiFi Hi-C contact at 50 kb: the well-known t(X;10) constitutional translocation present in this karyotypically aneuploid cell line is rediscovered from sequence alone in one individual, and the same arms show enriched 3D contact. The pipeline does not require a population. The mouse genome was processed with the same pipeline on the T2T assemblies of B6 and CAST [@Francis2025], yielding a 2-community arm-level partition and a per-PHR-pair Jaccard matrix. The per-PHR-pair Jaccard between mouse subtelomeres correlates with mouse zygotene Hi-C contact at Spearman ρ = 0.715, p = 4.4 x 10^-55, n = 344 inter-chromosomal pairs (Fig. 4d) [@Zuo2021]. The correlation is present at all four meiotic stages (leptotene, zygotene, pachytene, diplotene; ρ 0.574 to 0.715) and peaks at zygotene, the meiotic bouquet stage [@bouquet_BhattTERBEvolution2020; @Patel2019]. The sequence-to-3D coupling generalises to one of the two species in which the meiotic bouquet has been most extensively characterised.

The five lines of evidence close a four-link causal loop (Extended Data Fig. 8a). Sequence sharing predicts three-dimensional proximity at the meiotic bouquet (Mantel Spearman ρ = 0.66 in CHM13 and HG002 Hi-C, 0.485 in HG002 Pore-C, 0.715 in mouse zygotene Hi-C). Three-dimensional proximity at the bouquet creates structural opportunity for ectopic recombination between non-homologous chromosomes. Ectopic recombination generates new shared sequence between the partner arms: 16 crossover-like and 133 gene-conversion-like patches in the WashU pedigree, 11 cross-assembler-validated parent x chr-pair features in the CEPH1463 pedigree. The new shared sequence enters the population, strengthens sequence-similarity edges and reinforces the original 3D proximity at the next meiosis. Each link is measured, not inferred. We use "concerted evolution" in the loose sense intended by Arnheim [@concerted_evolution_nahr_Arnheim1980], Ohta [@concerted_evolution_nahr_Ohta1984] and Charlesworth and colleagues [@concerted_evolution_nahr_Charlesworth1994]: ongoing inter-chromosomal recombination exchange that homogenises sequence between non-homologous chromosomes. The pedigree results provide the direct empirical anchor that earlier statistical and phylogenetic treatments lacked [@concerted_evolution_nahr_Hillis1991; @concerted_evolution_nahr_Vollger2023; @Vollger2023]. Four limitations bound the inference. First, bulk Hi-C in lymphoblastoid cell lines measures somatic chromatin organisation, not the meiotic bouquet directly; the 3D signal is consistent with envelope tethering and with single-sperm scHi-C evidence [@Xu2025; @hic3d_scnanoHiC2023; @hic3d_scnanoHiC2_2025] but is not a direct measurement of meiotic chromosome conformation in human. Second, we do not yet have genome-wide LAD or Lamin B1 ChIP-seq in matched germline cell types; Dip-C radial position is the current proxy (C1 D4Z4 radial 0.732 peripheral, C14 PAR2 0.840 peripheral, C10 chr17p 0.474 interior) and compartment identity at tips is weakly A-leaning (mean e1 = +0.0073 across HG002 100 kb Hi-C; 63 of 92 arm-by-haplotype windows in A, Extended Data Fig. 8d). Third, Lalli et al. [@Lalli2025] reported a centiMorgan-per-megabase anti-correlation between subtelomeric recombination rate and cross-arm affinity (Spearman ρ = -0.35, n = 46), but the correlation collapses to ρ ≈ 0 (n = 40) once seven low-callability arms are excluded (Extended Data Fig. 8c): short-read recombination maps cannot resolve PHRs and long-read maps are required [@Sasani2019; @Smolka2024]. Fourth, the 12% wfmash sampling rate is justified by Erdős and Rényi connectivity at the chosen 95% identity threshold; a stricter threshold could lose sensitivity to deeper-divergence inter-chromosomal blocks. The principal outlook is that long-read recombination maps in trios, matched germline LAD or HiCAR data, and extension of the cross-assembler pedigree intersection to the full CEPH1463 4-generation T2T set [@acrocentric_Porubsky2025denovo; @Porubsky2025] will close the remaining open links and quantify, for the first time, how much of the human variant landscape at chromosome ends is generated by inter-chromosomal recombination in each generation [@Logsdon2024].

## Methods

**Sample selection and reference frame.** 233 HPRC v2 v1.1 individuals contributed 465 haplotype-resolved assemblies (one HPRC haplotype was excluded for quality) and were complemented by CHM13v2.0 for a reference total of 466 haplotype-equivalent units. Per-superpopulation breakdown: AFR 67, EAS 52, AMR 44, SAS 37, EUR 33 [@hprc_hprcv2_2025; @Nurk2022].

**Telomere-anchored 500 kb flank extraction.** For each haplotype assembly we extracted the terminal 500 kb of every contig classified as a p- or q-arm telomere-bearing contig (contig length minimum 1 Mb), producing 18,827 flanks across 48 arms. Per-arm flank counts are in Extended Data Fig. 1b; pq-classification at `pq-classification/contig_classifications.tsv`.

**wfmash all-vs-all alignment.** `wfmash v0.23.0-41-gb5f0ff1c -p 95 -t 48 --quiet`; each flank serves as target in turn [@Guarracino2023].

**Implicit pangenome graph and IMPG transitive closure.** The all-vs-all PAF set is the implicit pangenome graph; IMPG `query -x` computes reachability [@pangenome_graphs_impg_GarrisonGuarracino2023; @pangenome_graphs_impg_GuarracinoHeumos2022; @pangenome_graphs_impg_IMPG2023; @pangenome_graphs_impg_Hickey2024]. Justification of "no chromosomal partitioning": realised sampling rate 11.6%, rounded to 12%, sits 230x above the Erdős and Rényi connectivity threshold p* = log(n)/n = 5.2 x 10^-4 for n = 18,827.

**PHR detection.** IMPG sliding-window scan with identity ≥ 95%, ≥ 5 alignments per chromosome on ≥ 2 different chromosomes, total aligned length per window ≥ 3 kb; output: 15,668 PHRs (83.2% of flanks), 41 of 48 arms.

**Pangenome graph and Jaccard similarity.** `pggb -p 95 -D /scratch`; `odgi similarity --all -P` over the 15,668 sequences yields the 15,668 x 15,668 Jaccard matrix [@Garrison2024pggb; @Guarracino2023].

**Arm-level distance matrix.** The 41 x 41 arm-level matrix is the mean of pairwise sequence-level Jaccard distances within each arm pair (`hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`).

**Community detection (Leiden).** Arm-level: edge weight w_ij = exp(-d_ij / median(d)); resolution scanned from 0.1 to 3.0 at 0.01 step; selected resolution 1.16 (mean silhouette 0.347, k = 15). Sequence-level: k-NN graph with k in {10, 25, 50, 75, 100, 125} and resolution 0.1 to 3.0; selected k = 75, resolution 0.8 (modularity 0.97, mean silhouette 0.602, k = 50).

**UPGMA dendrogram.** `hclust(..., method = "average")` on the 41 x 41 distance matrix; k = 14 (mean silhouette 0.342); agreement with Leiden 14 of 15.

**Neighbour-joining tree.** `ape::nj()` on the same 41 x 41 distance matrix; rooted at the MRCA of the acrocentric short-arm clade; 1,000-replicate perturbation bootstrap with Gaussian noise at sigma = 25% of the off-diagonal IQR.

**Heterogeneity tests.** Wilcoxon paired test (allele vs paralog distance); per-arm Spearman correlation between number of chromosome contributors and distance from the telomere; piecewise regression with changepoint detection; 2 x 5 Fisher exact test for cross-arm vs self-arm superpopulation composition with Benjamini and Hochberg correction; Hudson pairwise F_ST [@subtel_popgen_hudson1992; @subtel_popgen_weir1984].

**Hi-C, Pore-C and CiFi pipeline.** mcool inputs at 5, 10, 20, 50 and 100 kb; per-haplotype mode preserves maternal and paternal contigs as separate arms; MAPQ filters disabled to retain multi-mappers with one random alignment per read; within-vs-between (W/B) ratio computed by bootstrap on 10,000 permutations; Mann-Whitney global test; Mantel Spearman correlation with 10,000 row-and-column permutations; observed-over-expected inter-chromosomal normalisation [@hic3d_dixon2012; @hic3d_imakaev2012; @Ulahannan2019; @hic3d_cifi2025].

**Exclusion controls.** Five exclusion sets (no acrocentric p-arms, no sex chromosomes, no acrocentric p plus sex, no all-acrocentric plus sex, no strongest community) applied at all five mcool resolutions; results in Extended Data Fig. 5a, 5b.

**Single-cell 3D.** Dip-C in 16 GM12878 cells [@Tan2018] remapped to T2T-CHM13v2.0; sperm scHi-C in 20 cells [@Xu2025; @hic3d_scnanoHiC2023; @hic3d_scnanoHiC2_2025]. Negative control S_all: pooled 7 zero-signal arms (chr5p, chr6q, chr7q, chr12q, chr14q, chr20p, chr20q).

**Mouse pipeline.** B6 and CAST T2T assemblies [@Francis2025]; telomere-anchored flanks; 1, 2 and 4 Mb window scan; 2-community arm-level partition; per-PHR-pair Jaccard vs Zuo et al. 2021 zygotene Hi-C contact at leptotene, zygotene, pachytene and diplotene [@Zuo2021].

**Pedigree odgi-untangle.** `odgi untangle nth-best=1` per flank; high-quality filter: minimum patch score 0.95 and minimum alignment score 0.95; within-Leiden filter applied as credibility constraint; pattern classification (gene_conversion_like, crossover_like, acros_like, sandwich_same_hap, complex) via `scripts/pedigree/analyze-pedigree-recombination.py` [@Cechova2025; @Porubsky2025].

**Cross-assembler intersection (CEPH1463).** Independent calls by hifiasm and verkko; within-community parent x chr-pair features kept; 11 features pass [@acrocentric_Porubsky2025denovo].

**Software versions.** wfmash 0.23.0-41-gb5f0ff1c; impg commit 5b96025; pggb and odgi bundled; samtools, bedtools and bgzip via conda; hicexplorer 3.7.4; R packages `ape`, `vegan` and `cluster`.

**Data and code availability.** GitHub `ekg/phrs`; on-disk roots `/moosefs/guarracino/HPRCv2/PHR_III/` and `/moosefs/erikg/phrs/`; CHM13-coordinate PHR BEDs at repo root (`chm13.phrs.bed`, `chm13.phrs.no_acro.bed`, `CHM13-HG002.sub-telo-phrs.bed`).

**Limitations.** (i) bulk Hi-C is somatic LCL, not germline meiotic; (ii) no genome-wide LAD or Lamin B1 ChIP-seq in matched germline cell types; (iii) the short-read cM/Mb anti-correlation reported by Lalli 2025 collapses once seven low-callability arms are excluded, so the rate-vs-affinity relationship is underpowered; (iv) the 12% wfmash sampling rate is justified by Erdős and Rényi connectivity but rests on the chosen 95% identity threshold; (v) flank size truncation at 500 kb may underestimate PHR length; (vi) Hi-C N is small (six individuals) and confidence intervals are correspondingly wide.

## References

acrocentric_Altemose2022
acrocentric_Guarracino2025ape
acrocentric_Porubsky2025denovo
acrocentric_rdna_robertsonian_bandyopadhyay2001
Ambrosini2007
Bailey2002
BaileyEichler2006
Bergstrom2020
bouquet_BhattTERBEvolution2020
bouquet_ChikashigeTelomere1994
bouquet_DingSUN12007
bouquet_HarperBouquet2004
bouquet_HornKASH52013
bouquet_KotaSUN1MAJIN2020
bouquet_MorimotoKASH2012
bouquet_Scherthan2001
bouquet_Scherthan2003
bouquet_ShibuyaRPMs2015
bouquet_ZicklerKleckner1999
Brown1990
Cabianca2012
Cechova2025
concerted_evolution_nahr_Arnheim1980
concerted_evolution_nahr_Charlesworth1994
concerted_evolution_nahr_Eichler2001
concerted_evolution_nahr_Hastings2009
concerted_evolution_nahr_Hillis1991
concerted_evolution_nahr_Myers2010
concerted_evolution_nahr_Ohta1984
concerted_evolution_nahr_SamonteEichler2002
concerted_evolution_nahr_Vollger2023
dux4_d4z4_fshd_lemmers2007
dux4_d4z4_fshd_lemmers2010worldwide
Eichler2001
Flint1997
Francis2025
Garrison2018
Garrison2024pggb
Guarracino2023
hic3d_alavattam2019
hic3d_cifi2025
hic3d_deshpande2022
hic3d_dixon2012
hic3d_imakaev2012
hic3d_scnanoHiC2023
hic3d_scnanoHiC2_2025
hic3d_wolff2018
hprc_hprcv2_2025
Lalli2025
Liao2023
Linardopoulou2005
Logsdon2021
Logsdon2024
Masny2004
Mefford2001
MeffordTrask2002
Nurk2022
Ottaviani2009
OttavianiGilson2008
pangenome_graphs_impg_GarrisonGuarracino2023
pangenome_graphs_impg_GuarracinoHeumos2022
pangenome_graphs_impg_Hickey2024
pangenome_graphs_impg_IMPG2023
Patel2019
Porubsky2025
Riethman2001
Riethman2004
Riethman2008
Rouyer1986
Rudd2009
RuizHerrera2008
Sasani2019
sexchrompars_acquaviva2020
sexchrompars_bellott2024
Sharp2006
Skaletsky2003
Smolka2024
StankiewiczLupski2002
StankiewiczLupski2010
Stong2014
subtel_popgen_1000g2010
subtel_popgen_anderson2008
subtel_popgen_bhatia2013
subtel_popgen_hudson1992
subtel_popgen_levysakin2019
subtel_popgen_lewontin1972
subtel_popgen_rosenberg2002
subtel_popgen_weir1984
subtelstruct_Nergadze2007
subtelstruct_NergadzeITS2007
subtelstruct_NergadzeITSReview2007
Tan2018
Trask1991
Trask1998
Ulahannan2019
Vollger2023
Wilkie1991
Xu2025
ZicklerKleckner1998
ZicklerKleckner2015
Zuo2021

## Figure list

- Fig. 1: Population-scale subtelomeric communities. (a) Genome-wide stacked identity heatmap of 465 HPRC v2 haplotypes plus CHM13 across 24 chromosomes in 100 kb windows. (b) Genome-wide number-of-chromosome-contributors heatmap with PHR-BED overlay. (c) 41 x 41 arm-level Jaccard distance heatmap, arms ordered by Leiden k = 15 community, with UPGMA k = 14 dendrogram on top. (d) Per-arm architecture (homogeneous, polymorphic, fully interchangeable) and cross-arm sequence rate.
- Fig. 2: Within-community heterogeneity, two-domain model and population history. (a) Per-community Wilcoxon paired allele-vs-paralog distance. (b) Per-arm Spearman gradient and piecewise-regression breakpoints for the two-domain model. (c) Cross-arm superpopulation enrichment and Hudson F_ST. (d) UPGMA out-of-Africa tree from the F_ST matrix.
- Fig. 3: Three-dimensional nuclear organisation mirrors sequence communities. (a) HG002 Pore-C inter-arm contact at 50 kb, B/W = 0.056. (b) Forest plot of 14 inter-arm 3D tests across Hi-C, Pore-C, CiFi, Dip-C, sperm scHi-C and mouse meiotic Hi-C. (c) S_all negative-control comparison in GM12878 Dip-C and 20-cell sperm scHi-C. (d) Flanking paradox: PHR vs flanking 100 kb B/W; Dip-C radial inset.
- Fig. 4: Pedigree-resolved exchanges and cross-species generalisation. (a) WashU 3-generation T2T pedigree: 538 high-quality patches, 494 within Leiden communities. (b) CEPH1463 4-generation pedigree: 11 hifiasm-and-verkko-validated features. (c) RPE-1 self: Leiden C2 = {chr10_q, chrX_q} and elevated CiFi contact. (d) Mouse zygotene Hi-C vs per-PHR-pair Jaccard, Spearman ρ = 0.715, n = 344.
- Extended Data Fig. 1: Pipeline schematic, per-arm flank counts, PHR length distribution, chr18q chimera control.
- Extended Data Fig. 2: Sequence-level 50-community UMAP, within-community Jaccard bimodality, cross-arm affinity chord, arm-vs-sequence confusion matrix.
- Extended Data Fig. 3: TAR1 prevalence per arm, internal (TTAGGG)n island length and motif composition, terminal telomere length by community, per-arm TAR1 distance from telomere.
- Extended Data Fig. 4: PHR-only GO biological process enrichment, copy-weighted vs deduplicated GO, top-15 high-copy gene families (OR4F, IL9R, DUX4, FRG2), OR4F pseudogenisation gradient.
- Extended Data Fig. 5: Multi-resolution W/B robustness, Mantel before-vs-after acrocentric and sex exclusion, observed-over-expected within vs between, 15-community by 11-dataset reproducibility heatmap.
- Extended Data Fig. 8: Causal feedback loop, D4Z4-CTCF-lamin tethering model for C1, Lalli 2025 cM/Mb anti-correlation honest null, compartment diagnostic at tips.
- Extended Data nj_tree_arms: Neighbour-joining tree of arm-level Jaccard distances; six monophyletic clades with 100% MRCA bootstrap support.
