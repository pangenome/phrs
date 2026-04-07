## Integrated 3D interpretation

### Convergent evidence across technologies

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
| Dip-C (T2T) | GM12878 (16 cells) | Spearman: community-free arm-level | rho=0.336 | 1.1e-18 |
| Sperm 3D | 20 cells (Xu et al. 2025) | Within/between 3D distance | 60% closer | 3.9e-51 (Fisher) |
| Sperm 3D | 20 cells | Community-free per-cell | 15/20 positive rho | — |
| RPE-1 CiFi | Async (interphase) | B/W ratio (50kb / 10kb) | 0.024 / 0.033 | 9.1e-102 / 4.2e-79 |
| RPE-1 Pore-C | Async (interphase) | B/W ratio (50kb / 10kb) | 0.031 / 0.048 | 1.3e-95 / 1.7e-95 |
| RPE-1 CiFi | Mitotic | B/W ratio (50kb / 10kb) | 0.008 / 0.030 | 3.3e-68 / 3.0e-56 |
| Mouse 1Mb | 4 meiotic stages | B/W ratio range (5 resolutions) | 0.029–0.122 | all < 1e-22 |
| Mouse 4Mb | 4 meiotic stages | B/W ratio range (5 resolutions) | 0.019–0.057 | all < 1e-36 |

Additionally, per-PHR-pair correlation confirms the continuous relationship: arms with more similar subtelomeric sequence show proportionally more inter-chromosomal contact (CHM13 Hi-C: Spearman rho = 0.674, p = 3.0e-92; HG002 Pore-C: rho = 0.485, p = 1.6e-48; all 8 datasets significant). Dip-C community-based Mantel confirms this in 3D coordinates (rho = 0.296, p = 0.002), and the community-free arm-level analysis shows Spearman rho = 0.336, p = 1.1e-18.

*Note: All analyses in this table — human Hi-C/Pore-C/CiFi, RPE-1, and mouse — are now complete at all 5 mcool resolutions (5kb, 10kb, 20kb, 50kb, 100kb). Full multi-resolution tables are in the resolution sensitivity section (human), the acrocentric exclusion control (no-acrocentric), the RPE-1 validation section.1 (RPE-1), and the mouse flanking Hi-C section (mouse). Results are consistent across all resolutions, with slight strengthening at finer resolutions reflecting better bin coverage of the median 105 kb PHR regions.*

**Conclusion.** All technologies — using different cell types (LCL, RPE-1, sperm), different platforms (Illumina, PacBio, ONT), and different measurement principles (ligation frequency, multi-way contacts, 3D coordinates) — independently confirm community-structured 3D organization. The addition of Dip-C T2T remapping and sperm single-cell 3D data extends validation to haploid cells and eliminates hg19/T2T coordinate incompatibility concerns (the limitations section, limitation 12). The signal is not an artifact of any single technology, cell type, or analysis method. Multi-resolution analysis at 5 resolutions (5kb-100kb) confirms resolution-invariance across all systems.

### Flanking region paradox

**What it does.** Explains the counterintuitive finding that unique-sequence regions flanking PHR boundaries show stronger 3D signal than the duplicated PHR regions themselves.

**How.** PHR intervals contain inter-chromosomal duplications (the inter-chromosomal detection section). When Hi-C reads map to these duplicated regions, multi-mapping between community partner arms creates ambiguity — a read aligning to a duplicated block could originate from either partner, and multi-mapper removal (RM_MULTI=1 in HiC-Pro) discards these reads. The 100 kb flanking regions centromere-ward of the PHR boundary are unique sequence with no multi-mapping.

**Key metrics.** Flanking B/W ratios range from 0.002 (HG002) to 0.057 (CHM13) — all stronger than PHR B/W (0.027–0.074, see the flanking analysis section for per-sample detail), consistent with multi-mapping suppression at duplicated PHR intervals. HG002 flanking enrichment is 13x stronger than PHR (0.002 vs 0.027). In Dip-C, flanking particles are more interior than non-flanking terminal particles (GM12878: 0.503 vs 0.551, p = 7.4e-35).

**Conclusion.** Three implications: (1) Rules out multi-mapping as the driver of the 3D signal — flanking regions have no duplicated content, yet show stronger clustering. (2) 3D clustering extends beyond the duplicated region into flanking unique sequence, indicating broader chromosomal domain effects. (3) Multi-mapping suppresses the PHR signal; flanking enrichment values represent the true signal magnitude.

### Meiotic bouquet as exchange venue

**What it does.** Considers the meiotic bouquet as the context for subtelomeric exchange events.

**Key metrics.** All 3D data is somatic (interphase). Tan et al. (2018) found Rabl configuration "weak" in GM12878 and PBMCs. Zuo et al. (2021) showed chromosome end alignment extends "a substantial range of ~20% of chromosome length" in mouse meiosis. Average meiotic loop sizes: ~500 kb at leptotene, ~700 kb at zygotene. Median PHR region (105 kb) fits within a single meiotic loop.

The hypothesis that meiotic telomere clustering favours subtelomeric exchange has a long pedigree: Mefford & Trask (2002) noted that "the pairing of homologues typically begins at the telomeres" and that subtelomeric homology creates opportunities for ectopic exchange, while Linardopoulou et al. (2005) stated explicitly: "Telomere clustering in meiotic cells might favour exchange of chromosome ends during DSB healing." Mefford & Trask also noted that the interphase clustering of chr4q and chr10q (Stout et al. 1999) "could promote their exchange" — a prediction now supported by Hi-C/Pore-C data showing community-structured 3D co-localization across all 15 communities. In mouse meiosis, Patel et al. (2019) first detected X-shaped interchromosomal Hi-C contacts consistent with the meiotic bouquet, and Zuo et al. (2021) extended this with stage-resolved analysis showing that chromosome end alignment during early prophase (leptotene and zygotene) extends over "a substantial range of ~20% of chromosome length" — far deeper than the subtelomeric regions analyzed here. Notably, this alignment is not merely a manifestation of the transient bouquet conformation (<5% of zygotene cells show a cytological bouquet); rather, it occurs whenever different chromosomes are brought into proximity by telomeres, is independent of compartment identity, and depends on the force-transmitting LINC complex (which modulates alignment range but does not alter loop sizes). Average meiotic chromatin loop sizes increase from ~500 kb at leptotene to ~700 kb at zygotene, reaching ~1.4 Mb at pachytene and ~1.6 Mb at diplotene (Zuo et al. 2021; Patel et al. 2019 reported slightly larger values of 0.8–1.0 Mb at zygotene). Since the median PHR region length is 105 kb, most PHR sequences would reside at the base of a single meiotic loop at leptotene — the compartment where recombination machinery is concentrated — maximizing the probability of ectopic recombination between aligned PHR sequences on different chromosomes.

If human meiosis follows a similar pattern, the community-structured 3D contacts observed in somatic cells may underestimate the true meiotic proximity. The observation that flanking regions (100 kb centromere-ward) show stronger 3D signal than PHRs themselves (the flanking paradox section) is consistent with the meiotic proximity extending beyond the duplicated region.

Human meiotic Hi-C remains the single most informative missing experiment. Existing data from mouse spermatocytes (Patel et al. 2019) cannot be directly extrapolated to human subtelomeric organization.

### D4Z4-CTCF-Lamin tethering model for C1

**What it does.** Proposes a specific molecular mechanism for C1 (chr4_q/chr10_q) co-localization, supported by existing literature.

**Key metrics.** C1: silhouette = 0.147, 43.4% discordance (chr4_q). Dip-C radial = 0.732 (peripheral). Inter-chromosomal signal peaks at 0–15 kb (D4Z4 position). C1 sequences: median 22 DUX4L; non-C1 outliers: 0–2 (Mann-Whitney p = 5.3e-6).

**How.** The mechanism:

1. Both chr4_q and chr10_q carry D4Z4 macrosatellite arrays at their subtelomeric tips
2. CTCF binds within D4Z4 repeat units (Ottaviani et al. 2009)
3. D4Z4-proximal sequences are tethered to the nuclear periphery via lamin A/C interaction (Masny et al. 2004; Ottaviani et al. 2009)
4. Both arms are positioned at the nuclear periphery via lamin A/C interaction, and this co-localization is consistent with the elevated ectopic recombination observed clinically

The Dip-C radial analysis supports this: C1 arms occupy peripheral nuclear positions (among the most exterior communities). The pangenome data provides evidence consistent with D4Z4 contributing to the chr4_q/chr10_q co-clustering (the sequence-level community detection): inter-chromosomal signal peaks at 0–15 kb from the telomere (where D4Z4 sits), C1 sequences carry median 22 DUX4L genes while all 7 non-C1 outlier sequences have 0–2 on their own arm (Mann-Whitney p = 5.3e-6), and outlier PHR regions are 4.6–9x shorter than C1 sequences. As established in the FSHD literature (Lemmers et al. 2010), ectopic exchange between chr4_q and chr10_q D4Z4 arrays can modify FSHD alleles — a D4Z4 contraction on a permissive 4qA haplotype that translocates to chr10_q loses pathogenicity, while the reverse can create disease alleles.

This mechanism predicts that CTCF/cohesin density at PHR boundaries should correlate with Hi-C contact strength between community partners. This is testable using Gershman et al.'s (2022, Science) ENCODE CTCF ChIP-seq realignment to T2T-CHM13, which found CTCF enrichment at TAR loci across all ENCODE cell lines, and Stergachis lab Fiber-seq data providing single-molecule CTCF maps at 39/46 telomeres. Standard hg38-aligned ENCODE data would be inadequate because hg38 has incomplete subtelomeric assemblies.

### Nucleolar association mechanism for C6/C7

**What it does.** Considers nucleolar co-localization as a mechanism for C7 (acrocentric p-arms) homogenization and C6 (acrocentric q-arms) community membership.

**Key metrics.** C7: silhouette = −0.029, gene replacement scores 0.91–1.0 for chr13_p/chr14_p/chr15_p, 0.49–0.54 for chr21_p/chr22_p. C6: silhouette = 0.521 (sequence separation), Dip-C radial = 0.505 (3D interior positioning, consistent with nucleolar co-localization). C7 cannot be assessed in Dip-C (hg19 p-arms unmapped).

**Result.** All five acrocentric short arms carry rDNA and constitutively associate with nucleoli. The near-complete interchangeability (silhouette = −0.029) is consistent with frequent exchange at nucleolus-associated arms. C6 includes non-acrocentric arms (chr1_q, chr17_q, chr19_q), indicating q-arm membership is driven by shared duplicon content, not nucleolar proximity alone.

**Conclusion.** Nucleolar co-localization is consistent with the literature but the present study provides no C7-specific 3D data. The f7501 duplicon distribution (Mefford & Trask 2002) maps directly onto C3 (chr3_q, chr7_p, chr9_q, chr11_p, chr16_q, chr19_p), confirming that community structure captures known duplicon module relationships.

### Causal feedback loop

**What it does.** Proposes a feedback loop: sequence similarity → 3D proximity → ectopic exchange → increased similarity.

**How.** Four links with varying levels of support:

1. **Sequence similarity → 3D proximity**: Mantel tests show continuous positive correlation between Jaccard similarity and Hi-C contact (HG002 rho=0.66, Pore-C rho=0.49)
2. **3D proximity → ectopic exchange**: Established from FSHD literature (chr4_q/chr10_q D4Z4 translocation) and general principles of recombination requiring physical proximity
3. **Ectopic exchange → increased similarity**: Inferred from the outcome — cross-arm affinity analysis shows 15.9% of sequences resemble a partner arm more than their own, consistent with past exchange having increased inter-arm similarity

The fourth link — **increased similarity → stronger future proximity** — is inferred but not directly measured. Testing this would require temporal data (e.g., comparing 3D contacts before and after an exchange event) or comparing Hi-C contact strength between haplotypes that carry cross-arm vs self-arm subtelomeric types within the same individual.

This circularity means the causal direction cannot be established from the present data: communities could form because similar sequences are brought into proximity, or sequences could become similar because they are in proximity. The meiotic bouquet (the meiotic bouquet section) provides a speculative but plausible initiation scenario — telomere clustering during meiosis brings all chromosome ends into proximity regardless of sequence content, and subsequent ectopic exchange between neighboring arms would create the initial sequence similarity that then reinforces itself through the feedback loop. This bouquet initiation model cannot be tested with the present data and is offered as a plausible scenario, not an evidence-based conclusion. This feedback concept is not new: Linardopoulou et al. (2005) described a cycle of "segmental polymorphism and gross genomic rearrangement" (their Fig. 2) where translocations create duplications that promote further rearrangements, and Mefford & Trask (2002) noted that interphase clustering of chr4q/chr10q "could promote their exchange." The present analysis adds the 3D proximity dimension with quantitative support. Additionally, Ambrosini et al. (2007) proposed that the 98% identity peak in their bimodal duplicon distribution "may be due to maintenance of sequence similarity by ongoing interchromosomal gene conversion between the large subtelomeric duplicons" — a hypothesis consistent with the ongoing-exchange link in this model.

### Testable predictions

Three predictions arising from the meiotic chromosome organization data of Zuo et al. (2021) that could not be tested with the present data:

1. **LINC complex requirement**: The LINC complex (Linker of Nucleoskeleton and Cytoskeleton) spans the nuclear envelope and transmits cytoskeletal forces to chromosomes, driving telomere-led movements during meiotic prophase. Zuo et al. showed that a SUN1 point mutation (W151R) disrupts long-range chromosome-end alignment in mouse meiosis: while tip contacts (within ~5% of chromosome length) actually increase, alignment at greater distances drops off sharply, reducing the effective alignment zone from ~20% to ~5% of chromosome length. If meiotic Hi-C from SUN1 mutant spermatocytes were analyzed with the community framework, within-community inter-chromosomal contacts at chromosome ends should be dramatically reduced compared to wild-type. This would establish that LINC-mediated force transmission is required for community-structured 3D contacts and, by extension, for the ectopic exchange that maintains subtelomeric homology.

2. **Crossover frequency correlation** (tested, confounded): Using the T2T-CHM13 recombination map (Lalli et al. 2025, preprint), subtelomeric recombination rate anticorrelates with cross-arm affinity across all 39 shared arms (rho = −0.43, p = 0.006). However, this signal is entirely driven by 7 confounded arms: acrocentric p-arms (chr13_p, chr14_p, chr15_p, chr21_p, chr22_p) and PAR (chrX_p, chrX_q) have 0–12 callable variants in 500 kb (vs 1,000–3,000 for non-acrocentric arms), so their 0 cM/Mb recombination rate reflects absence of short-read genotyping data in repetitive regions, not necessarily absence of recombination. Excluding these 7 arms: rho = 0.00, p = 0.98, N = 32 — the correlation vanishes. The question of whether local recombination protects arm identity from ectopic exchange remains biologically plausible but cannot be tested with current recombination map data at subtelomeric loci. Long-read-based recombination maps resolving variants in repetitive subtelomeric regions would be needed (see the testable predictions section prediction 7 for the full per-arm table).

3. **Chromatin compartment identity**: Zuo et al. showed that A-compartment (transcriptionally active) regions form shorter meiotic loops (~560 kb at leptotene) with higher crossover rates, while B-compartment regions form longer loops (~730 kb). Compartment calling from the HG002 Hi-C eigenvector (100 kb resolution, per-haplotype, GC-oriented) shows that 68% of chromosome tips are classified as A-compartment (63/92 arm × haplotype combinations), indicating a weak gene-rich signature at subtelomeric tips. However, the eigenvector values at tips are close to zero (mean +0.007), indicating that subtelomeric regions have weak, poorly defined compartment identity — consistent with their nature as transitional zones between chromosome-specific and duplicated sequence. The Dip-C radial data (the flanking analysis section) adds a spatial dimension: PHR-flanking regions are more interior than other terminal regions, and the communities with the most inter-chromosomal sharing (C6, C10) are among the most interior (C10 radial=0.474, C6=0.505). Subtelomeric regions are thus pseudogene-rich with weak and poorly defined compartment identity (68% A-compartment by GC-oriented eigenvector, but mean e1 close to zero), yet positioned internally rather than at the lamina. This dissociation suggests that telomere clustering, not lamina association, determines the nuclear positioning of these regions. Per Zuo et al., weak B-compartment identity genome-wide predicts longer meiotic loops (~730 vs ~560 kb). However, Zuo explicitly found that B-compartment regions at centromere-distal chromosome ends have shorter loops than B-compartment elsewhere, and stated that loop size differences are "unlikely the main cause for prominent alignment among chromosome ends" — the end alignment pattern is independent of compartment identity. This favors bouquet-stage telomere clustering rather than loop-axis accessibility as the primary driver of inter-chromosomal exchange at subtelomeres.

---

The mechanistic models above propose specific feedback loops between sequence similarity and 3D proximity. The following section situates these findings within the broader literature and identifies novel contributions.

