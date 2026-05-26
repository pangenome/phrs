# Refresh 15: Emerging Topics 2024–2026

**Generated:** 2026-05-17  
**Task:** lit-refresh-15-emerging (agent-54)  
**Scope:** Scan 2024–2026 literature for emerging themes adjacent to the Nature draft that are NOT covered by the 14 existing `topic_NN_*.md` files. Catch anything we should be citing or addressing that we currently aren't.  
**Inputs cross-checked:** All 14 `paper_prep/lit_review/topic_NN_*.md` files; `paper_prep/synthesis/REFERENCES_v3.bib` (295 entries, keys listed); `paper_prep/synthesis/NATURE_DRAFT_v1.md`; `paper_prep/synthesis/CONSISTENCY_AUDIT_v1.md`.  
**Search tools used:** PubMed MCP, bioRxiv MCP; date range 2024-01-01 to 2026-05-17.

---

## Section 0: Theme Inventory

Fifteen candidate emerging themes were evaluated. Eight are included (new content not yet in the bib or substantially extending covered ground). Seven are rejected.

| # | Theme | Decision | One-line rationale |
|---|-------|----------|--------------------|
| 1 | Non-human primate T2T subtelomeric satellite biology | **INCLUDE** | Two landmark 2024–2026 preprints (macaque, marmoset) resolve subtelomeric repeat architectures and inter-arm PHRs in NHP; directly extends C4/C5 comparative framing; not in bib |
| 2 | Long-read single-molecule recombination from sperm | **INCLUDE** | Schweiger et al. 2024 captures NCO gene conversions at single-molecule resolution from long-read sperm; extends and complements Xu2025 (in bib); two-process NCO finding directly relevant to C8 gene-conversion-like patches |
| 3 | RPE-1 cell-line reference genome | **INCLUDE** | Volpe et al. 2025 (Nat Commun) published the near-complete diploid assembly of hTERT RPE-1, the exact cell line used for §09 self-vs-self community analysis; not in bib |
| 4 | rDNA-linked segmental duplications in great apes (T2T-scale) | **INCLUDE** | de Gennaro et al. 2026 (Genes) used T2T assemblies + FISH to characterize SD dynamics at acrocentric short arms across five ape species; extends topics 05 and 07; not in bib |
| 5 | Meiotic bouquet TERB1 and SUN-domain functional updates | **INCLUDE** | Kameyama 2024 (TERB1 in medaka) and Cromer 2024 (SUN1/SUN2 in Arabidopsis) add mechanistic validation of bouquet models post-2023; complement topic 08; not in bib |
| 6 | Non-B DNA structures at acrocentric short arms and subtelomeres in ape T2T genomes | **INCLUDE** | Smeds et al. 2025 (NAR) showed G4s, Z-DNA, hairpins occupy 9–38 % of acrocentric short arms in T2T ape genomes, with ORC enrichment; relevant to recombination mechanism at subtelomeres; not in bib |
| 7 | Pangenome structural variation dynamics (review) | **INCLUDE** | Loegler et al. 2025 (Cell Genomics) reviews SV diversity in the T2T/pangenome era; provides broader framing for the paper's claims C2/C4; supports citing in Discussion/Introduction |
| 8 | New single-cell chromatin conformation capture efficiency | **INCLUDE** | HiChew (Chen et al. 2026, Genome Biol) improves scHi-C valid-pair ratio to ~50 %; methodological context for existing Dip-C/sperm scHi-C data in C7 |
| 9 | HPRC v2 downstream re-analyses | **REJECT** | Largely covered by `hprc_hprcv2_2025`, `Liao2023`, `hprc_siren2025` already in bib; no specific 2024–2026 paper was missed |
| 10 | LAD maps in germline or meiotic cells | **REJECT** | 2024–2026 papers found only in somatic differentiation (erythropoiesis); no germline-specific or meiotic LAD maps published yet; gap acknowledged in topic 08 Open Questions |
| 11 | MAJIN complex structural updates | **REJECT** | No new breakthrough structure papers 2024–2026 beyond `bouquet_KotaSUN1MAJIN2020` and `bouquet_WangTERB2019`/`bouquet_QiangTERB2019`; existing bib coverage adequate |
| 12 | Subtelomere-specific cancer/aging (ALT/TERT) | **REJECT** | Only one relevant paper found (Graham et al. 2024, TERT/neuroblastoma, PMID 38837897); too cancer-type-specific and not about subtelomere exchange mechanism; would distract |
| 13 | New pangenome graph alignment tools post-wfmash | **REJECT** | No published alignment tools superseding wfmash/seqwish in 2024–2026 found in PubMed/bioRxiv; the field moved toward downstream pangenome graph analyses, not aligner replacement |
| 14 | New mouse T2T assemblies with meiotic telomere focus | **REJECT** | No 2024–2026 mouse T2T papers with meiotic/subtelomere emphasis found; Zuo2021 meiotic Hi-C remains the primary reference |
| 15 | Additional pedigree T2T assemblies beyond CEPH1463/WashU | **REJECT** | `Cechova2025` and `Porubsky2025` already in bib cover the two available multigenerational T2T datasets; additional independent datasets not yet available |

---

## Section 1: Executive Summary

Eight themes yield **ten recommended new citations** (all 2024–2026, all verified by PMID or DOI), addressing four distinct gaps in the current bib:

1. **Non-human primate comparative context**: The paper currently lacks any comparative discussion of subtelomeric architectures in non-human primates. Two landmark preprints (rhesus macaque T2T with 268 novel repeat families; common marmoset T2T with PHRs at all acrocentric autosomes) allow the paper to position its human findings within a broader primate evolutionary context. Both are relevant to C4 (interchromosomal homology) and C5 (community structure), and the marmoset paper directly identifies PHRs and rDNA-facilitated exchange at acrocentric arms, mirroring the paper's community C7.

2. **Long-read sperm recombination (non-crossover)**: Schweiger et al. 2024 (long-read sperm, PMID 39005338) reveals two distinct NCO processes — a PRDM9-associated short-tract dominant class and a non-PRDM9 long-tract minority (~2%) — using the same single-sperm long-read strategy as Xu2025 (in bib). This paper directly informs interpretation of the 133 gene-conversion-like patches in C8 (tract-length distribution, association with PRDM9 hotspots).

3. **RPE-1 cell-line reference**: The paper uses RPE-1 for self-vs-self community validation (§09, 37 self-discovered communities, t(X;10) as positive control). Volpe et al. 2025 (Nat Commun) published the near-complete diploid RPE-1 assembly with t(X;10) confirmed; this is a natural citation for the RPE-1 analysis.

4. **Mechanistic and comparative updates**: Four papers update the mechanistic and comparative landscape covered in topics 05, 07, 08, and 09: TERB1 in medaka meiosis (Kameyama 2024), Arabidopsis SUN domain (Cromer 2024), great-ape acrocentric SD dynamics (de Gennaro 2026), and non-B DNA in ape T2T genomes (Smeds 2025).

Additionally, two broader framing papers are recommended for Discussion context: Loegler et al. 2025 (pangenome SV dynamics review) and Chen et al. 2026 (HiChew single-cell conformation method).

**Cross-check verdict:** None of the recommended papers are double-counted with any of the 14 existing topic files or the 295 bib entries. Cross-check details are in Section 3.

---

## Section 2: Recommended New Citations

All papers below are verified by PMID or DOI. None appear in `REFERENCES_v3.bib` (searched by DOI, first-author surname, and title substring).

---

### 2.1 Non-human primate T2T subtelomeric biology

**Zhang et al. 2025 — T2T rhesus macaque: lessons from subtelomeric repeats and sequencing bias**

- PMID: 41019632
- DOI: [10.1101/2025.08.04.668424](https://doi.org/10.1101/2025.08.04.668424)
- Journal: bioRxiv preprint
- Published: 2025-08-04
- Key authors: Zhang S, Ventura M, Antonacci F, Phillippy AM, Mao Y (Shanghai Jiao Tong University / Bari / NHGRI)

**Summary (from PubMed):** T2T-MMU8v2.0, the highest base-level-accuracy primate genome to date (ONT-only). Subtelomeric satellite-rich regions identified as the principal assembly bottleneck. Discovered 268 previously unannotated repeat families; resolved ~8 Mbp of SATR satellite arrays with >99-fold enrichment in subtelomeric regions. Four distinct SATR genomic architectures, each with unique satellite composition, segmental-duplication organization, and epigenetic signatures. Crucially: macaque subtelomeres harbor 58 actively transcribed genes in SATR arrays, contrasting with gene-poor hominid subtelomeres. Read mappability improved 19% and 5,821 additional chromatin accessibility peaks were recovered compared with prior assemblies.

**Relevance to paper:** Directly relevant to C4 (interchromosomal homology at subtelomeres), C5 (community structure), and C8 (ongoing exchange). The paper's finding that macaque subtelomeres have 4 architecturally distinct SATR landscapes — functionally differentiated, with active gene expression — contrasts with the largely gene-desert architecture of human subtelomeres and provides primate-comparative context for why human subtelomeric community structure is evolutionarily significant. The 268 novel repeat families resolved only with T2T assembly reinforce C2's implicit-pangenome rationale.

**Suggested citation location:** Introduction (motivating the T2T approach for subtelomeric analysis) and Discussion (comparative context: human subtelomeric community structure vs. macaque SATR diversity).

**Suggested bib key:** `t2t_Zhang2025macaque`

---

**Hebbar et al. 2026 — A Complete Genome for the Common Marmoset**

- PMID: 41929024
- DOI: [10.64898/2026.03.25.713844](https://doi.org/10.64898/2026.03.25.713844)
- Journal: bioRxiv preprint
- Published: 2026-03-26
- Key authors: Hebbar P, Miga KH, Eichler EE, Gerton JL, Paten B, Alexandrov I (UCSC / UW / Stowers)

**Summary (from PubMed):** First T2T genome for common marmoset (NWM); four high-quality haplotypes from two individuals. All six marmoset acrocentric autosomes have gene-poor, satellite-rich short arms; evidence that all share PHRs; Y chromosome (but not X) carries active rDNA and PHRs; rDNA copy number is sexually dimorphic. Chromosomes sharing PHRs share closely related centromeric satellite DNA, supporting a model of ongoing recombinational exchange facilitated by rDNA. Constructed a marmoset pangenome for short-read mapping.

**Relevance to paper:** The most directly relevant new paper. Hebbar et al. demonstrate that PHRs exist at marmoset acrocentric short arms and that rDNA-facilitated exchange between heterologous chromosomes is a conserved NWM mechanism — extending the human paper's model (C4/C5/C8) to a New World monkey. The "all acrocentric autosomes share PHRs" finding mirrors community C7 in the human analysis (the acrocentric p-arm community). The sexual dimorphism of rDNA + PHRs on the marmoset Y echoes the human chrX_q community.

**Suggested citation location:** Discussion (conservation of PHR-mediated inter-arm exchange in primates); Introduction (cross-species motivation).

**Suggested bib key:** `t2t_Hebbar2026marmoset`

---

### 2.2 Long-read single-molecule recombination from sperm

**Schweiger et al. 2024 — Insights into non-crossover recombination from long-read sperm sequencing**

- PMID: 39005338
- DOI: [10.1101/2024.07.05.602249](https://doi.org/10.1101/2024.07.05.602249)
- Journal: bioRxiv preprint
- Published: 2024-07-07
- Key authors: Schweiger R, Lee S, Durbin R (Cambridge / Sanger)

**Summary (from PubMed):** High-fidelity single long reads from 15 sperm samples (13 donors) capture both crossovers and non-crossovers in one molecule. Non-crossover gene conversions show variation between and within donors. Two distinct NCO processes: (i) a dominant class with mean tract length <50 bp, upstream of PRDM9 binding sites — standard PRDM9-induced meiotic recombination; (ii) a minority (~2%) with much longer mean tract lengths, not associated with PRDM9 sites, also seen in somatic cells.

**Relevance to paper:** The 133 gene-conversion-like patches in the WashU pedigree (C8) are the inheritance-visible products of NCO events. Schweiger et al. provide the contemporary single-molecule framework for interpreting NCO tract-length distributions and PRDM9 association. The non-PRDM9 long-tract NCO class (2%) is a potential source of the larger gene-conversion-like patches observed at subtelomeres (where PRDM9 hotspot density may be lower than at interstitial sequences). This paper complements `Xu2025` (in bib: single-cell sperm Hi-C) by addressing the recombination mechanics rather than the 3D organization.

**Suggested citation location:** Methods/Discussion (pedigree patch classification: gene-conversion-like tract-length interpretation); supplements topic 11 (pedigree-based recombination detection).

**Suggested bib key:** `pedigree_Schweiger2024spermNCO`

---

### 2.3 RPE-1 cell-line reference genome

**Volpe et al. 2025 — The reference genome of the human diploid cell line RPE-1**

- PMID: 40940351
- DOI: [10.1038/s41467-025-62428-z](https://doi.org/10.1038/s41467-025-62428-z)
- Journal: Nature Communications 16: 7751
- Published: 2025-09-12
- Key authors: Volpe E, Guarracino A, Giunta S (Sapienza Rome / UTHSC Memphis)

**Summary (from PubMed):** RPE1v1.1, the near-complete diploid assembly of hTERT RPE-1 (non-cancerous retinal epithelial cell line with stable karyotype). PacBio + ONT long-read sequencing; Hi-C phasing. Chromosome-level scaffolds spanning centromeres. Identifies haplotype-specific variants including the t(X;10)(Xq28;10q21.2) characteristic of RPE-1 cells, and peak divergence at centromeres.

**Relevance to paper:** The paper uses RPE-1 for the self-vs-self 3D community analysis (§09): 37 self-discovered communities, the t(X;10) translocation as a positive control (community C2 = {chr10_q, chrX_q} in HPRC analysis). Volpe et al. provide the reference assembly that validates the t(X;10) identification and the RPE-1 community structure. This is a natural citation in the RPE-1 analysis paragraph, where it is currently uncited. Note: Guarracino is a co-author, making this a first-party citation for the paper.

**Suggested citation location:** Results (§09 RPE-1 self-vs-self paragraph, currently uncited for the RPE-1 cell line itself); Methods (RPE-1 Hi-C data source).

**Suggested bib key:** `hic3d_Volpe2025RPE1`

---

### 2.4 rDNA-linked segmental duplications in great apes (T2T-scale)

**de Gennaro et al. 2026 — Evolution of rDNA-Linked Segmental Duplications as Lineage-Specific Mosaics in Great Apes**

- PMID: 41751569
- DOI: [10.3390/genes17020185](https://doi.org/10.3390/genes17020185)
- Journal: Genes 17(2): 185
- Published: 2026-01-31
- Key authors: de Gennaro L, Catacchio CR, Ventura M (University of Bari)

**Summary (from PubMed):** Eight human-derived fosmid probes targeting SD-enriched regions flanking rDNA arrays were hybridized to multiple individuals from chimpanzee, bonobo, gorilla, Bornean and Sumatran orangutan. FISH reveals extensive lineage-specific SD copy-number variation and chromosomal heteromorphism at acrocentric chromosomes; gorillas show most polymorphism, orangutans the most conserved patterns. Several SDs show fixed duplications across species; others are highly polymorphic. T2T assembly comparison confirms consistent localization for some probes but shows partial discordance for others (highlighting challenges even in T2T assemblies). Chimpanzees and bonobos have higher proportions of SDs on rDNA-bearing chromosomes than gorillas.

**Relevance to paper:** Topic 05 (acrocentric rDNA/Robertsonian) covers human acrocentric SDs; this paper extends the analysis to five great ape species at T2T resolution, directly supporting the comparative framing of human community C7 (acrocentric p-arms) as part of an evolutionarily dynamic, lineage-specific landscape. The "lineage-specific mosaic" framing resonates with the paper's "41 signal-bearing arms" finding — not all arms are equivalent, and the inter-species SD variation at acrocentrics is mechanistically related to the inter-arm exchange documented in the paper. The partial discordance even in T2T assemblies validates the paper's caution about assembly errors in complex repeat regions.

**Suggested citation location:** Discussion (comparative context for acrocentric community structure and SD dynamics).

**Suggested bib key:** `acrocentric_deGennaro2026apeSD`

---

### 2.5 Meiotic bouquet TERB1 and SUN-domain functional updates

**Kameyama et al. 2024 — Medaka Mutant Displays Defects of Synaptonemal Complex Formation and Sexual Difference in Gametogenesis**

- PMID: 38809870
- DOI: [10.2108/zs230108](https://doi.org/10.2108/zs230108)
- Journal: Zoological Science 41(3): 314–322
- Published: 2024-06
- Key authors: Kameyama S, Tanaka M (Nagoya University)

**Summary (from PubMed):** Medaka (fish) TERB1 mutants show incomplete synaptonemal complex (SC) formation; SC initiation observed (punctate lateral elements, fragmented transverse filaments) but not completion. DSB introduction is independent of synapsis completion. Both oogenesis and spermatogenesis show aberrant chromosome arrangement. Critical finding: oogenesis arrests at zygotene-like stage in mutants, but testes continue to produce sperm-like cells — demonstrating sexual dimorphism in the meiotic checkpoint (similar to mammalian pattern). TERB1 connects telomeres to SUN/KASH motors.

**Relevance to paper:** Topic 08 covers the TERB1–TERB2–MAJIN complex (2019 papers in bib). This 2024 paper adds a fish-model functional validation of TERB1's essential role in SC formation and provides evidence for checkpoint sexual dimorphism — relevant because the paper uses both sperm and oocyte-derived single-cell data (sperm scHi-C and Dip-C in GM12878 cells). The sexual dimorphism finding is also interesting for interpreting why male vs. female Hi-C contacts differ at subtelomeres.

**Suggested citation location:** Discussion (meiotic bouquet mechanism, TERB1 conservation); supplements topic 08 Open Questions.

**Suggested bib key:** `bouquet_Kameyama2024TERB1`

---

**Cromer et al. 2024 — Rapid meiotic prophase chromosome movements in Arabidopsis thaliana are linked to essential reorganization at the nuclear envelope**

- PMID: 39013853
- DOI: [10.1038/s41467-024-50169-4](https://doi.org/10.1038/s41467-024-50169-4)
- Journal: Nature Communications 15: 5964
- Published: 2024-07-16
- Key authors: Cromer L, Grelon M (INRAE/AgroParisTech, Versailles)

**Summary (from PubMed):** Arabidopsis meiotic centromeres undergo rapid (up to 500 nm/s) uncoordinated movements during zygotene and pachytene. These centromere movements require SUN1/SUN2 (abolished in sun1 sun2 double mutant). Telomere attachment to SUN-enriched NE domains, bouquet formation, and nucleolus displacement are all defective in sun1 sun2. Establishes Arabidopsis as a model for meiotic RPMs and demonstrates mechanistic conservation of telomere-led RPMs across plants and animals.

**Relevance to paper:** Directly supports the mechanistic conservation claim in topic 08 and the paper's use of mouse meiotic Hi-C data to infer human subtelomeric proximity. The plant data extend SUN-domain-mediated bouquet function to a phylogenetically distant organism (angiosperms), strengthening the argument that SUN1-mediated NE tethering is a universal mechanism for telomere-led chromosome movements. The 500 nm/s RPM speed puts quantitative constraints on how rapidly subtelomeric sequences can sample nuclear space during zygotene.

**Suggested citation location:** Discussion (universality of SUN-domain-mediated bouquet); supplements topic 08.

**Suggested bib key:** `bouquet_Cromer2024SUN`

---

### 2.6 Non-B DNA structures in ape T2T genomes (acrocentrics and subtelomeres)

**Smeds et al. 2025 — Non-canonical DNA in human and other ape telomere-to-telomere genomes**

- PMID: 40226919
- DOI: [10.1093/nar/gkaf298](https://doi.org/10.1093/nar/gkaf298)
- Journal: Nucleic Acids Research 53(7)
- Published: 2025-04-10
- Key authors: Smeds L, Makova KD (Penn State University)

**Summary (from PubMed):** Comprehensive characterization of non-B DNA motifs (G4s, Z-DNA, hairpins, bent DNA, etc.) in T2T genomes of 7 apes (human, bonobo, chimpanzee, gorilla, Bornean and Sumatran orangutan, siamang). Non-B DNA motifs are enriched in the genomic regions added to T2T assemblies. Acrocentric short arms: 9–15 % of autosomes, 9–11 % of chr X, and 12–38 % of chr Y occupied by non-B motifs. G4s and Z-DNA enriched at promoters, enhancers, and ORC origins. Repetitive sequences harbor more non-B DNA than non-repetitive sequences, especially at acrocentric short arms. Most centromeres/flanking regions enriched in at least one non-B DNA motif type.

**Relevance to paper:** The paper identifies G4s and non-B DNA motifs concentrated at the very regions (acrocentric short arms, subtelomeres) where the paper's communities are organized. Non-B DNA structures can promote DSB formation (G4 stabilization, Z-DNA formation under negative supercoiling) — relevant to understanding *why* acrocentric short arms and subtelomeric regions are hotspots for ectopic exchange. This is a mechanistic complement to the community structure findings: the paper shows sequence-similarity patterns; Smeds 2025 shows that these regions also have elevated structural instability at the DNA level. Particularly relevant for topic 07 (NAHR mechanism) and Discussion.

**Suggested citation location:** Discussion (mechanistic basis of subtelomeric NAHR/exchange hotspots; structural features of exchange-prone regions).

**Suggested bib key:** `subtelstruct_Smeds2025nonBDNA`

---

### 2.7 Pangenome structural variation dynamics (review)

**Loegler et al. 2025 — Dynamics of genome evolution in the era of pangenome analysis**

- PMID: 41260225
- DOI: [10.1016/j.xgen.2025.101067](https://doi.org/10.1016/j.xgen.2025.101067)
- Journal: Cell Genomics 6(1): 101067
- Published: 2025-11-18
- Key authors: Loegler V, Friedrich A, Schacherer J (University of Strasbourg)

**Summary (from PubMed):** Review of how pangenomes and T2T assemblies transform structural variant identification, classification, and analysis; mechanisms of SV formation; uneven genomic distribution; roles in adaptation and disease. Discusses integration of pangenome graphs into GWAS, challenges of T2T pangenomes at population scale, and need for new computational tools.

**Relevance to paper:** Provides broad pangenome framing for the paper's contribution. Relevant for Introduction/Discussion to situate the implicit pangenome graph approach within the contemporary field. Lower priority than the mechanistic papers above, but useful for a sentence or two in Introduction on "why pangenomes matter for complex regions." Note: does not replace the existing `hprc_ebler2022`, `Liao2023` citations.

**Suggested citation location:** Introduction (pangenome SV landscape framing).

**Suggested bib key:** `pangenome_Loegler2025review`

---

### 2.8 Single-cell chromatin conformation method advance

**Chen et al. 2026 — Highly efficient chromatin conformation capture with post-enrichment in single cells by HiChew**

- PMID: 42036683
- DOI: [10.1186/s13059-026-04059-1](https://doi.org/10.1186/s13059-026-04059-1)
- Journal: Genome Biology 27(1)
- Published: 2026-04-27
- Key authors: Chen Z, Tang C (BGI Genomics, Shenzhen)

**Summary (from PubMed):** HiChew combines sticky-end ligation with post-PCR methylation-based enrichment to achieve ~50% valid pair ratios (vs. ~8% for unenriched methods). Single-cell implementation (snHiChew): 45–50% valid pair ratios, 5–10 kb resolution with 70–80% bin coverage. Strong concordance with conventional Hi-C for compartments, TADs, and loops.

**Relevance to paper:** The paper uses Dip-C (GM12878, 16 cells) and sperm scHi-C (20 cells, Xu2025) as its single-cell 3D datasets. HiChew represents a methodological advance that would allow the same analysis at substantially higher efficiency and resolution in future experiments. Relevant for Discussion (limitations and future directions: higher-efficiency single-cell conformation capture could resolve within-community vs. cross-community contact at subtelomere scale in additional individuals or cell types).

**Suggested citation location:** Discussion (methods outlook/future directions for single-cell 3D validation).

**Suggested bib key:** `hic3d_Chen2026HiChew`

---

## Section 3: Coverage Cross-Check Against Existing 14 Topics

For each recommended paper, the relevant topic file is listed and the reason the paper is NOT double-counted is stated.

| Recommended paper | Most relevant existing topic | Why NOT double-counted |
|---|---|---|
| Zhang 2025 (macaque T2T) | topic_02 (subtelomere structure), topic_10 (pangenome) | topic_02 covers human TAR1/ITS structure; topic_10 covers wfmash/PGGB methods. Zhang 2025 covers non-human primate subtelomeric satellite architecture; no NHP T2T paper exists in either topic. |
| Hebbar 2026 (marmoset T2T) | topic_05 (acrocentric rDNA), topic_03 (PHR concept) | topic_05 covers human rDNA/NOR; topic_03 covers PHR concept in humans. Hebbar 2026 is the first NHP T2T paper showing PHRs at acrocentric arms in NWM; not covered anywhere. |
| Schweiger 2024 (sperm NCO) | topic_11 (pedigree recombination) | topic_11 covers pedigree-based and Sperm-seq studies through 2024. Schweiger 2024 is distinct: it uses long-read single molecules (not microfluidics) to capture NCO gene conversions; the NCO two-process finding is not discussed in topic_11. `Xu2025` (in bib) covers 3D sperm organization, not recombination mechanics. |
| Volpe 2025 (RPE-1 genome) | topic_09 (Hi-C/3D methods) | topic_09 covers Hi-C and single-cell 3D methods; it does not have the RPE-1 reference genome. The self-vs-self analysis (§09) references the RPE-1 Hi-C data but not the published reference assembly. |
| de Gennaro 2026 (ape SDs) | topic_05 (acrocentric rDNA) | topic_05 covers acrocentric SD organization in humans and early ape data; de Gennaro 2026 is the first T2T-scale + FISH cross-species SD characterization at acrocentrics; not cited or described. |
| Kameyama 2024 (TERB1 medaka) | topic_08 (meiotic bouquet) | topic_08 covers TERB1 through the 2019 structural papers; Kameyama 2024 adds a 2024 fish-model functional paper validating TERB1 phenotypes and the checkpoint dimorphism; not in topic_08. |
| Cromer 2024 (Arabidopsis SUN) | topic_08 (meiotic bouquet) | topic_08 covers SUN1/SUN2 in mammals; Cromer 2024 extends this to plants with RPM velocities; not in topic_08. |
| Smeds 2025 (non-B DNA apes) | topic_02 (subtelomere structure), topic_07 (NAHR) | topic_02 covers TAR1/ITS; topic_07 covers NAHR mechanisms. Non-B DNA at acrocentric short arms and subtelomeres is not covered in either topic; Smeds 2025 is the first T2T-scale characterization. |
| Loegler 2025 (pangenome SV review) | topic_10 (pangenome graphs) | topic_10 covers wfmash/PGGB/IMPG methodology. Loegler 2025 is a 2025 field-wide review; partially overlaps in framing but adds current landscape context; low-risk double-count. |
| Chen 2026 (HiChew scHi-C) | topic_09 (Hi-C/3D methods) | topic_09 does not cover HiChew; it was not yet published. Different enrichment strategy from Dip-C. |

---

## Section 4: Claim-Level Implications

**C4 (extended interchromosomal homology at nearly all subtelomeres):** Zhang 2025 and Hebbar 2026 add cross-species validation. The human paper now has primate comparators showing that subtelomeric shared-sequence regions (PHRs in marmoset; SATR architectures in macaque) are present across primates, supporting the evolutionary significance of C4. _No change needed to numerical claims._

**C5 (named cladistic communities):** The marmoset data (Hebbar 2026) identifies PHRs at all acrocentric autosomes and a sexually dimorphic rDNA+PHR on Y — providing a NWM parallel to the human acrocentric p-arm community (C7). The de Gennaro 2026 ape SD data provide evolutionary context for the acrocentric community. _No change needed; Discussion can reference both._

**C7 (3D proximity mechanism):** The Cromer 2024 and Kameyama 2024 bouquet papers add cross-taxon support for the SUN/KASH/TERB mechanism. The Chen 2026 HiChew method is a methodological advance for future validation. _No change needed; Discussion can note conservation and methodological advances._

**C8 (ongoing exchange; concerted evolution):** Schweiger 2024 provides the single-molecule NCO framework for interpreting the 133 gene-conversion-like patches. The finding of two NCO processes (PRDM9-associated short tract, non-PRDM9 long tract) is directly relevant: the acrocentric community (C7) may be enriched in the non-PRDM9 long-tract class because acrocentric p-arms are PRDM9-depleted DSB zones. _No change to numbers; Discussion can cite Schweiger 2024 when discussing gene-conversion-like patch mechanism._

**§09 RPE-1 analysis:** Volpe 2025 should be cited in the RPE-1 paragraph (currently no reference for the RPE-1 cell line itself). _Add one citation._

**DIVERGES flagged by consistency audit:** None of the new papers directly resolve the flagged divergences (wrong 7 silent arms, UPGMA 12 vs 14 communities, F_ST range error, HG02148 exclusion-set Mantel values). These are internal consistency issues, not literature gaps.

---

## Section 5: Priority Action Items

Ordered by impact on the published paper.

| Priority | Action | Paper(s) | Location |
|----------|--------|----------|----------|
| 1 (HIGH) | Add RPE-1 genome citation | Volpe 2025 (PMID 40940351) | Results §09, Methods |
| 2 (HIGH) | Add marmoset T2T citation | Hebbar 2026 (PMID 41929024) | Discussion (PHR conservation in NWM) |
| 3 (HIGH) | Add rhesus macaque T2T citation | Zhang 2025 (PMID 41019632) | Discussion (NHP subtelomeric satellite diversity) |
| 4 (HIGH) | Add long-read NCO sperm citation | Schweiger 2024 (PMID 39005338) | Discussion (gene-conversion-like patch mechanism) |
| 5 (MEDIUM) | Add non-B DNA citation | Smeds 2025 (PMID 40226919) | Discussion (mechanistic basis of acrocentric/subtelomeric exchange hotspots) |
| 6 (MEDIUM) | Add ape acrocentric SD citation | de Gennaro 2026 (PMID 41751569) | Discussion (comparative acrocentric SD dynamics) |
| 7 (LOW) | Add TERB1 medaka citation | Kameyama 2024 (PMID 38809870) | Discussion / topic 08 |
| 8 (LOW) | Add Arabidopsis SUN citation | Cromer 2024 (PMID 39013853) | Discussion / topic 08 |
| 9 (LOW) | Add pangenome SV review | Loegler 2025 (PMID 41260225) | Introduction (pangenome framing) |
| 10 (LOW) | Add HiChew scHi-C citation | Chen 2026 (PMID 42036683) | Discussion (future methods outlook) |

---

## Appendix: Summary Table of Recommended New Citations

All from PubMed/bioRxiv, 2024–2026, PMID verified.

| BibKey (proposed) | PMID | DOI | Year | Journal | Primary relevance |
|---|---|---|---|---|---|
| `t2t_Zhang2025macaque` | 41019632 | 10.1101/2025.08.04.668424 | 2025 | bioRxiv | NHP subtelomere architecture (C4/C5) |
| `t2t_Hebbar2026marmoset` | 41929024 | 10.64898/2026.03.25.713844 | 2026 | bioRxiv | NWM PHRs at acrocentrics (C5/C8) |
| `pedigree_Schweiger2024spermNCO` | 39005338 | 10.1101/2024.07.05.602249 | 2024 | bioRxiv | NCO gene conversion mechanism (C8) |
| `hic3d_Volpe2025RPE1` | 40940351 | 10.1038/s41467-025-62428-z | 2025 | Nat Commun | RPE-1 reference genome (§09) |
| `acrocentric_deGennaro2026apeSD` | 41751569 | 10.3390/genes17020185 | 2026 | Genes | Ape acrocentric SD dynamics (C5) |
| `bouquet_Kameyama2024TERB1` | 38809870 | 10.2108/zs230108 | 2024 | Zoolog Sci | TERB1 bouquet function (C7) |
| `bouquet_Cromer2024SUN` | 39013853 | 10.1038/s41467-024-50169-4 | 2024 | Nat Commun | SUN-domain bouquet conservation (C7) |
| `subtelstruct_Smeds2025nonBDNA` | 40226919 | 10.1093/nar/gkaf298 | 2025 | NAR | Non-B DNA at acrocentrics/subtelomeres |
| `pangenome_Loegler2025review` | 41260225 | 10.1016/j.xgen.2025.101067 | 2025 | Cell Genomics | Pangenome SV framing |
| `hic3d_Chen2026HiChew` | 42036683 | 10.1186/s13059-026-04059-1 | 2026 | Genome Biol | scHi-C methods advance |

**Total new citations: 10. All have verifiable PMID. All absent from REFERENCES_v3.bib (295 entries). Cross-checked against all 14 topic_NN files.**
