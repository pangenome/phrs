---
title: NATURE_DRAFT_OUTLINE — section-by-section outline binding every claim to evidence and citations
author: agent-15 (draft-nature-article)
date: 2026-05-17
target_journal: Nature (Article)
target_word_counts:
  abstract: ~200
  main_text: 3000-4500
figure_count_main: 4
figure_count_extended_data: 7   # ED1, ED2, ED3, ED4, ED5, ED8, nj_tree_arms (ed6 + ed7 do not exist)
commit_hash_of_inputs: 8c364fe   # HEAD at outline drafting; PPTX_OUTLINE.md added in 89cab6a
inputs_read:
  - paper_prep/synthesis/PPTX_OUTLINE.md
  - paper_prep/synthesis/ABSTRACT_nature.md
  - paper_prep/synthesis/ABSTRACT_BoG.md
  - paper_prep/synthesis/CROSSWALK.md
  - paper_prep/synthesis/REFERENCES_v3.bib            # 295 @-keys
  - paper_prep/figures/{fig1,fig2,fig3,fig4}/caption.md
  - paper_prep/figures/{fig1,fig2,fig3,fig4}/sources.tsv
  - paper_prep/figures/{ed1,ed2,ed3,ed4,ed5,ed8}/caption.md
  - paper_prep/figures/nj_tree_arms/README.md
  - paper_prep/surveys/SURVEY_14_pedigree_recombination.md
  - end-to-end-report/report/{01_pipeline,14_pedigree_recombination}.md  (headings)
---

# Nature Article — Drafting Outline

Section-by-section skeleton that T3 (`write-nature-article`) will expand into prose. Every paragraph names its **Claim(s)**, **Evidence** (figures, surveys, report sections, on-disk tables), **Citations** (only keys present in `REFERENCES_v3.bib`), and **Approx word count**.

Numerical anchors used throughout (reconciled — see § *Number reconciliation* at end):

- **233** HPRC v2 individuals → **465** haplotype assemblies + **CHM13** reference (= **466** haplotype-equivalent units)
- **18,827** telomere-anchored 500 kb flanks across **48** arms
- **15,668** PHRs (83.2 % of flanks) on **41 of 48** signal-bearing arms
- **15** arm-level Leiden communities (silhouette 0.347, k=15) / **50** sequence-level Leiden communities (modularity 0.97, k=75/res 0.8)
- PHR length: **median 105 kb**, mean 144 kb, range 5–500 kb
- Hi-C anchors: HG002 B/W = 0.027 (50 kb, p = 4.0 × 10⁻⁶⁶); CHM13 B/W = 0.071 (p = 6.0 × 10⁻¹⁸); Mantel ρ = 0.66 both
- Pedigree (WashU T2T): **538** HQ inter-chr patches, **494 (92 %)** within Leiden communities, **133** gene-conversion-like, **16** crossover-like, **229** acros_like
- Mouse zygotene Hi-C: Spearman ρ = 0.715 (n = 344 inter-chr pairs, p = 4.4 × 10⁻⁵⁵)

---

## Title candidates

1. **(PPTX + Nature/BoG abstracts — keep)** *Concerted evolution and unorthodox recombination of human subtelomeres*
2. **(Refined; one-line message)** *Population-scale pseudohomology and ongoing inter-chromosomal exchange define human subtelomeres*
3. **(Alternative; mechanism-first)** *Human subtelomeres are unified by community-structured 3-D proximity and pedigree-observed recombination*

Recommendation: keep option 1 — it is the PPTX title, both abstract drafts use it, and "concerted evolution" reads in the *loose* sense (ongoing inter-chromosomal recombination exchange) per CROSSWALK §2 C8 / Erik's clarification.

---

## Abstract bullets (~200 words to expand from these 7)

1. **Background / gap.** Human subtelomeres are evolutionarily dynamic and structurally complex, but interchromosomal sequence relationships have resisted systematic analysis owing to assembly incompleteness and reference-bias in alignment. (Citations: `MeffordTrask2002`, `Riethman2004`, `Linardopoulou2005`, `BaileyEichler2006`.)
2. **Approach (sample + method).** Using **465** near-complete haplotype assemblies from HPRC v2 plus CHM13 (**466** total), we build an *implicit pangenome graph* — the all-vs-all PAF set produced by wfmash sampled at ~**12 %** of pairwise combinations (~230× above the Erdős–Rényi connectivity threshold), and queried via IMPG transitive closure — without chromosomal partitioning. (`hprc_hprcv2_2025`, `Garrison2024pggb`, `pangenome_graphs_impg_GarrisonGuarracino2023`, `pangenome_graphs_impg_GuarracinoHeumos2022`.)
3. **Landscape.** We map **18,827** telomere-anchored 500 kb flanks → **15,668** PHRs (median 105 kb) across **41 of 48** chromosome arms, comparable in scale to canonical PAR2. (`Rouyer1986`, `sexchrompars_acquaviva2020`.)
4. **Cladistics.** A neighbour-joining tree of arm-level Jaccard distances recovers PAR1, PAR2, the acrocentric short arms, the **10p–18p** clade, a tight **{22q,21q,19q,1q,13q,17q}** q-arm clade, and **4q–10q DUX4** with wide copy-number diversity — all six clades match the Leiden k = 15 partition (UPGMA agreement 14/15). (`MeffordTrask2002`, `Linardopoulou2005`, `dux4_d4z4_fshd_lemmers2010worldwide`, `acrocentric_Altemose2022`.)
5. **Population-genetic signature.** Allele-vs-paralog Wilcoxon, two-domain Spearman gradient (39/48 arms), Hudson F_ST and an out-of-Africa UPGMA tree resolve subtelomeric clustering across five human superpopulations. (`Flint1997`, `Ambrosini2007`, `subtel_popgen_anderson2008`, `subtel_popgen_bhatia2013`, `subtel_popgen_hudson1992`.)
6. **3-D + pedigree.** Bulk and single-cell Hi-C, Pore-C and CiFi in six individuals, GM12878 Dip-C and 20-cell sperm scHi-C, plus mouse meiotic Hi-C peaking at zygotene (ρ = 0.715), tie sequence similarity to nuclear-envelope proximity through the meiotic bouquet; a four-sample T2T pedigree provides direct evidence with **538** inter-chr patches, **92 %** within communities, **133** gene-conversion-like and **16** crossover-like events. (`hic3d_dixon2012`, `Tan2018`, `Xu2025`, `Zuo2021`, `Cechova2025`, `Porubsky2025`, `bouquet_KotaSUN1MAJIN2020`.)
7. **Implication.** Human subtelomeres are unified by ongoing inter-chromosomal recombination and concerted-in-the-loose-sense evolution; the events that build the population-scale communities are observable in single pedigrees within a single human generation.

---

## Section outline

> Nature Articles use a paragraph-led main text without formal Intro/Results/Discussion headers; Methods is a separate section. Below: ordered paragraphs P1–P11, each with claim/evidence/citations/word budget. Total target ≈ **3,610 words**, fits inside the 3,000–4,500 budget.

### P1: Subtelomeres as a frontier of human pangenomics

**Claim(s):**
- Human subtelomeres are dynamic, segmentally duplicated, and the classical example of inter-chromosomal sequence sharing (PAR1, PAR2, acrocentric rDNA short arms, D4Z4/FSHD).
- Prior population-scale work has been limited by reference incompleteness and by per-chromosome alignment frames; the HPRC v2 near-T2T pangenome of 233 individuals provides the substrate to revisit the question without chromosomal partitioning.
- This paper is the subtelomeric companion to the HPRC v2 main publication.

**Evidence:** `end-to-end-report/report/01_pipeline.md` introduction; PPTX slides 1–2 (`paper_prep/synthesis/PPTX_OUTLINE.md`); SURVEY_FRAMING_cytogenetic_fish.md; topic_01_cytogenetic_foundations.md; topic_02_subtelomere_structure.md.

**Citations:** `Brown1990`, `Wilkie1991`, `Trask1991`, `Trask1998`, `Riethman2001`, `Riethman2004`, `Riethman2008`, `MeffordTrask2002`, `Mefford2001`, `Linardopoulou2005`, `Rouyer1986`, `Flint1997`, `sexchrompars_acquaviva2020`, `dux4_d4z4_fshd_lemmers2010worldwide`, `acrocentric_Altemose2022`, `hprc_hprcv2_2025`, `Liao2023`, `Nurk2022`, `Logsdon2021`.

**Approx word count:** 300

---

### P2: The implicit pangenome graph: an unbiased map without chromosomal partitioning

**Claim(s):**
- Each haplotype is its own reference: the all-vs-all PAF set produced by wfmash (≥ 95 % identity) over 18,827 flanks is the implicit pangenome graph; IMPG transitive closure is the canonical query (`impg query -x`).
- The pipeline samples ≈ 12 % of pairwise haplotype combinations — ~230× above the Erdős–Rényi connectivity threshold p* = log(n)/n ≈ 5.2 × 10⁻⁴ for n = 18,827 — so transitive closure reaches virtually every subtelomere.
- Output: 18,827 flanks → 15,668 PHRs (83.2 %) on 41 of 48 arms.

**Evidence:** `paper_prep/figures/ed1/caption.md` (panel a, pipeline schematic); `end-to-end-report/report/01_pipeline.md` §"All-vs-all subtelomeric alignment", §"Inter-chromosomal region detection"; CROSSWALK §7a, §7b; PPTX slides 30, 37–42.

**Citations:** `hprc_hprcv2_2025`, `Garrison2018`, `Garrison2024pggb`, `pangenome_graphs_impg_GarrisonGuarracino2023`, `pangenome_graphs_impg_GuarracinoHeumos2022`, `pangenome_graphs_impg_IMPG2023`, `Guarracino2023`, `Liao2023`, `pangenome_graphs_impg_Hickey2024`.

**Approx word count:** 350

---

### P3: A genome-wide landscape of interchromosomal homology comparable to PAR2

**Claim(s):**
- Genome-wide stacked identity heatmaps show telomere-anchored, high-identity inter-chromosomal blocks at nearly every chromosome end.
- The detected PHRs span tens to hundreds of kilobases (median 105 kb, mean 144 kb, range 5–500 kb), comparable in scale to PAR2 (~334 kb).
- 41 of 48 chromosome arms carry signal; the 7 silent arms (chr7q, chr12q etc.) serve as built-in negative controls (also used in 3-D analysis as `S_all`).

**Evidence:** **Fig 1a** (`p_genome_wide_identity_heatmap.png`), **Fig 1b** (`p_genome_wide_numchrom_heatmap.png`); `paper_prep/figures/fig1/caption.md`; **ED1c** (PHR length distribution); `end-to-end-report/report/01_pipeline.md` §"Inter-chromosomal region detection"; CHM13-coordinate BED at repo root (`chm13.phrs.bed`, `chm13.phrs.no_acro.bed`).

**Citations:** `Rouyer1986`, `sexchrompars_acquaviva2020`, `sexchrompars_bellott2024`, `Bailey2002`, `BaileyEichler2006`, `RuizHerrera2008`, `Vollger2023`, `concerted_evolution_nahr_Vollger2023`, `Stong2014`.

**Approx word count:** 350

---

### P4: Neighbour-joining cladistics: six named clades recovered from one distance matrix

**Claim(s):**
- A neighbour-joining tree of the 41 × 41 arm-level Jaccard distance matrix recovers **PAR1 (Xp, Yp)**, **PAR2 (Xq, Yq)**, the **acrocentric short arms** (13p, 14p, 15p, 21p, 22p), the **10p–18p Linardopoulou clade**, the tight q-arm clade **{22q, 21q, 19q, 1q, 13q, 17q}**, and the **4q–10q DUX4** pair as monophyletic groups.
- Bootstrap support at the MRCA of all six abstract-named clades is 100 % under 1,000 perturbation replicates.
- Each NJ clade maps one-to-one to a Leiden k = 15 community (C15, C14, C7, C2, C6, C1 respectively); UPGMA k = 14 agrees on 14/15 communities — the cladistic signal is robust to algorithm choice.

**Evidence:** **Fig 1c** (41 × 41 Jaccard heatmap + UPGMA dendrogram); **`paper_prep/figures/nj_tree_arms/`** (`nj_tree.R`, `nj_tree_annotated.{pdf,png}`, `nj_tree.newick`, `README.md`); `hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`; `hprcv2.1Mb.subtelo.arm-upgma-k14.assignments.tsv`; `end-to-end-report/report/01_pipeline.md` §"Community detection"; CROSSWALK §2 C5 table.

**Citations:** `MeffordTrask2002`, `Linardopoulou2005`, `dux4_d4z4_fshd_lemmers2010worldwide`, `dux4_d4z4_fshd_lemmers2007`, `acrocentric_Altemose2022`, `acrocentric_Guarracino2025ape`, `sexchrompars_acquaviva2020`, `Cabianca2012`, `Skaletsky2003`, `Rudd2009`.

> PPTX_CONFLICT: PPTX slides 5–12 describe "PGGB graph + jaccard + Louvain community detection" and adapt Traag et al. (Leiden). The Nature abstract names "neighbor-joining trees"; Andrea's pipeline runs Leiden + UPGMA, no NJ. `nj_tree_arms/` resolves this conflict: NJ is now run on the arm-level Jaccard distance matrix, and the figure reports the 14/15 Leiden agreement so reviewers can see the cladistic signal is independent of the partitioning algorithm. Outline keeps NJ as the headline cladistic representation (matches abstract) and Leiden / UPGMA as supporting partitions (matches PPTX + pipeline).

**Approx word count:** 380

---

### P5: Three arm architectures: homogeneous, polymorphic, fully interchangeable

**Claim(s):**
- 41 signal-bearing arms partition into three architectural categories using cross-arm sequence rate and silhouette: **homogeneous** (4/41 — single-arm communities C8, C9, C10, C13), **polymorphic** (28/41 — multi-arm members with arm identity preserved), and **fully interchangeable** (9/41 — acrocentric p-arms C7 + PAR2 C14 + PAR1 C15, all with negative silhouette and reversed allele–paralog distance).
- chrX_q (99.7 %), chr21_p (94.0 %), chr11_p (74.1 %) top the cross-arm sequence rate; chrX_q ↔ chrY_q via PAR2 is the extreme case.

**Evidence:** **Fig 1d**; `paper_prep/figures/fig1/architecture_per_arm.tsv`; `cross_arm_affinity_sequences.tsv`; `hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv`; `end-to-end-report/report/04_heterogeneity.md` §"Cross-arm affinity"; SURVEY_04 §1.1–§1.2; SURVEY_04 §1.11.

**Citations:** `Flint1997`, `MeffordTrask2002`, `Linardopoulou2005`, `Riethman2001`, `Eichler2001`, `Stong2014`.

**Approx word count:** 240

---

### P6: Within-community heterogeneity, two-domain architecture, and a population-genetic signature

**Claim(s):**
- **Allele-vs-paralog Wilcoxon.** Across 5,946 paired distances in nine multi-arm communities, the allele is closer than the cross-paralog in 8/9 communities (overall p < 10⁻³⁰⁰); C7 (acrocentric p-arms) reverses with 70.5 % paralog-closer (p = 2.0 × 10⁻⁷) — quantitative confirmation of complete homogenisation at acrocentric short arms.
- **Two-domain Flint/Mefford model at pangenome scale.** 39 of 48 arms show a monotonic Spearman gradient between number-of-chromosome-contributors and telomere distance; 39/41 arms favour a two-segment piecewise fit over linear; 16/19 arms with internal (TTAGGG)n islands carry an ITS within 50 kb of the breakpoint.
- **Population structure.** Cross-arm vs self-arm Fisher 2 × 5 superpopulation tests: 10/19 significant after BH correction; Hudson F_ST: AFR vs non-AFR 0.10–0.15; an UPGMA dendrogram of the F_ST matrix recovers AFR as the deepest split with AMR/EAS/EUR/SAS forming a tight non-AFR clade.

**Evidence:** **Fig 2a–d**; `allele_vs_paralog_distance.tsv`; `two_domain_test.tsv` (48 arms); `two_domain_changepoint.tsv` (41 arms); `its_breakpoint_coloc.tsv`; `cross_arm_superpop_enrichment.tsv`; `fst_superpop_matrix.tsv`; **ED2b** (within-community Jaccard bimodality for C1/C2/C3/C5/C6/C7/C11/C12); **ED3a/b** (TAR1 prevalence + (TTAGGG)n island length); `end-to-end-report/report/04_heterogeneity.md`; SURVEY_04 §1.1, §1.3, §1.4; SURVEY_10_11_12_limits_summary_lit.md §6 T-4.

**Citations:** `Flint1997`, `Ambrosini2007`, `MeffordTrask2002`, `subtelstruct_NergadzeITS2007`, `subtelstruct_Nergadze2007`, `subtelstruct_NergadzeITSReview2007`, `subtel_popgen_anderson2008`, `subtel_popgen_bhatia2013`, `subtel_popgen_hudson1992`, `subtel_popgen_rosenberg2002`, `subtel_popgen_levysakin2019`, `subtel_popgen_1000g2010`, `Bergstrom2020`, `subtel_popgen_lewontin1972`, `subtel_popgen_weir1984`.

**Approx word count:** 420

---

### P7: Three-dimensional nuclear organisation mirrors sequence communities

**Claim(s):**
- Sequence-defined communities are physical: across 14 inter-arm 3-D tests (six Hi-C + HG002 Pore-C + HG002 CiFi + GM12878 Dip-C + 20-cell sperm scHi-C + four mouse meiotic stages), every effect lies left of unity (B/W or W/B < 1; range 0.020–0.93).
- HG002 Hi-C **B/W = 0.027** at 50 kb (p = 4.0 × 10⁻⁶⁶); CHM13 Hi-C B/W = 0.071 (p = 6.0 × 10⁻¹⁸); HG002 Pore-C B/W = 0.056 (p = 3.9 × 10⁻⁸⁵). **Mantel ρ = 0.66** (sequence-similarity ↔ Hi-C contact) for both CHM13 and HG002; per-arm-pair ρ = 0.66; per-individual sequence-pair ρ up to 0.83 in the lowest-coverage samples.
- Robustness: across 5 mcool resolutions (5–100 kb) and 5 exclusion sets (no acro-p / no sex / no acro-p + sex / no all-acrocentric + sex / no strongest), the Mantel correlation **strengthens** when subset confounds are excluded (HG002 0.66 → 0.80; HG02148 0.15 → 0.21; CHM13 0.66 → 0.85). The signal is not driven by acrocentric / nucleolar / PAR contacts.

**Evidence:** **Fig 3a/3b**; `paper_prep/figures/fig3/sources.tsv`; `analysis/human/community_based/{5,10,20,50,100}000bp/<sample>_global_test.tsv`; `analysis/human/exclusion_controls/*/`; **ED5a–d** (multi-resolution + exclusion robustness); `end-to-end-report/report/05_hic_validation.md`; SURVEY_07 §1, §1.2, §5 #11.

**Citations:** `hic3d_dixon2012`, `hic3d_imakaev2012`, `hic3d_alavattam2019`, `hic3d_wolff2018`, `hic3d_deshpande2022`, `hic3d_cifi2025`, `hic3d_scnanoHiC2023`, `hic3d_scnanoHiC2_2025`, `Ulahannan2019`, `Tan2018`, `Xu2025`, `Zuo2021`, `Wagner2022`.

**Approx word count:** 420

---

### P8: The flanking paradox, negative controls, and a meiotic-bouquet mechanism

**Claim(s):**
- Unique-sequence flanking regions (100 kb centromere-ward of PHR boundaries) show a **stronger** 3-D community signal than the duplicated PHRs themselves (HG002: PHR B/W 0.027 → flanking 0.0031; ≈ 9× stronger). Multi-mapping does **not** explain the within-community 3-D enrichment.
- The 7 zero-signal arms, pooled as the `S_all` pseudo-community, are systematically **farther** apart in 3-D (16/16 GM12878 + 20/20 sperm C-community cells have W/B < 1; 0/16 GM12878 and 1/20 sperm S_all cells do — S_all is 11 % / 40 % farther). Clean negative control.
- Mechanism: telomeres are attached to the nuclear envelope by the MAJIN–TERB2–TERB1 complex; the meiotic bouquet pre-localises subtelomeres at zygotene, providing a structural basis for ectopic exchange between non-homologous chromosomes. D4Z4-CTCF–lamin tethering at chr4q/chr10q (Masny 2004; Ottaviani 2009) instantiates the loop at C1.

**Evidence:** **Fig 3c/d** (S_all negative control + flanking paradox + Dip-C radial inset); **ED5b/c**; **ED8b** (D4Z4–CTCF–lamin schematic); `paper_prep/figures/fig3/sources.tsv` (`flanking_global_test`, `dipc_radial_inset`); SURVEY_07 §1.2, §1.6, §1.7; `end-to-end-report/report/07_integrated.md` §"Flanking region paradox", §"D4Z4-CTCF-Lamin model".

**Citations:** `bouquet_KotaSUN1MAJIN2020`, `bouquet_Scherthan2001`, `bouquet_Scherthan2003`, `bouquet_ShibuyaRPMs2015`, `bouquet_ChikashigeTelomere1994`, `bouquet_HarperBouquet2004`, `bouquet_HornKASH52013`, `bouquet_DingSUN12007`, `bouquet_MorimotoKASH2012`, `bouquet_ZicklerKleckner1999`, `ZicklerKleckner1998`, `ZicklerKleckner2015`, `Ottaviani2009`, `OttavianiGilson2008`, `Masny2004`, `Lemmers2010`, `Cabianca2012`.

**Approx word count:** 360

---

### P9: Direct pedigree evidence of ongoing inter-chromosomal exchange

**Claim(s):**
- The WashU T2T pedigree (PAN010 grandmother, PAN011 grandfather, PAN027 mother, PAN028 granddaughter; Cechova et al. 2025) yields **538** HQ inter-chromosomal `odgi untangle` patches; **494 (92 %)** sit within an HPRC v2 Leiden community.
- Pattern breakdown of the 494: **229 acros_like** (NAHR signature: ≥ 5 inter-chr patches from ≥ 3 source chromosomes in one flank — classical acrocentric pattern), **133 gene-conversion-like** (sandwich `chrN:hX → chrM:hY → chrN:hX` at predominantly perfect 1.000/1.000 alignment scores; ~ 90 % in community C7), **16 crossover-like** events (largest 27.97 kb on PAN028 maternal chr3q), **115 sandwich_same_hap**, **1 complex**.
- Inheritance directly observed across 3 generations: PAN027 ← PAN010 / PAN011; PAN028 ← PAN027. 12 of 16 crossover-like events are in PAN028, confirming that meiotic-resolution inter-chr breakpoints transmit.
- The CEPH1463 4-generation pedigree (Porubsky et al. 2025) yields **11 cross-assembler-validated** (hifiasm AND verkko) parent × chr-pair features, all within Leiden communities: chr10/chr18 (C2; the Linardopoulou pair) independently in NA12877 paternal and NA12878 maternal; chr12/chr9 (C5) in both NA12889 and NA12890 G1 grandparents; chr6/chr9 (C5) independently in NA12877 and NA12878; chr19/chr22 transmitted via NA12878.

**Evidence:** **Fig 4a, 4b**; `pedigrees/all_pedigrees_patches.tsv` (5,984 HQ patches); `pedigree-plots/washu/PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf`, `PAN027.paternal_hap2_from_PAN011_father.untangle.pdf`, `PAN028.maternal_hap1_from_PAN027_mother.untangle.pdf`; `pedigrees/ceph1463/{hifiasm,verkko}/untangle/recombination/patches.tsv`; `end-to-end-report/report/14_pedigree_recombination.md` §1.2–§1.6; SURVEY_14_pedigree_recombination.md §1.1, §1.2, §1.3, §1.6.

**Citations:** `Cechova2025`, `Porubsky2025`, `acrocentric_Porubsky2025denovo`, `concerted_evolution_nahr_Arnheim1980`, `concerted_evolution_nahr_SamonteEichler2002`, `concerted_evolution_nahr_Eichler2001`, `concerted_evolution_nahr_Hastings2009`, `concerted_evolution_nahr_Myers2010`, `StankiewiczLupski2002`, `StankiewiczLupski2010`, `Linardopoulou2005`, `Lemmers2010`, `dux4_d4z4_fshd_lemmers2010worldwide`.

> PPTX_CONFLICT: PPTX slide 26 comment from Erik Garrison flags "we need to check if the untangling is actually too conservative, and if this means we have a whole phr that swapped". The outline frames the 538 patches as a lower bound (HQ filter + Leiden-community filter) consistent with that concern — see SURVEY_14 §1.1 footnote on the 80-point WashU/CEPH gap as a fragmentation-noise effect that justifies the conservative filter rather than disqualifying it.

**Approx word count:** 460

---

### P10: Cross-individual and cross-species generalisation

**Claim(s):**
- **RPE-1 (single diploid) — methodological generality.** Self-discovered 37 communities on the 46-arm RPE-1 distance matrix independently recover **chrX_q × chr10_q** as a community (Leiden C2 = {chr10q, chrXq}), with elevated async-CiFi Hi-C contact at 50 kb. The well-known t(X;10) translocation is rediscovered from sequence alone, in a single individual.
- **Mouse (B6 + CAST T2T) — cross-species generality.** Per-PHR-pair Jaccard vs zygotene Hi-C contact (Zuo et al. 2021) shows **Spearman ρ = 0.715, p = 4.4 × 10⁻⁵⁵, n = 344 inter-chr pairs**; the correlation is present at all four meiotic stages (lepto/zygo/pachy/diplo, ρ 0.574–0.715) and peaks at zygotene — the meiotic-bouquet stage. The sequence ↔ 3-D coupling generalises beyond human.

**Evidence:** **Fig 4c, 4d**; `RPE1_subtelo/rpe1.dist_matrix.tsv`; `RPE1_subtelo/rpe1.communities.tsv` (37 communities); `rpe1_self_async_cifi_contact_matrix.tsv`; `mouse_T2T/community_analysis_1Mb/50000bp/zuo2021_zygotene_phr_pair_correlation.tsv` (+ leptotene / pachytene / diplotene stage tables); `mouse.arm-leiden.communities.tsv`; SURVEY_09 §1.3; SURVEY_08 §1.7; `end-to-end-report/report/08_mouse.md`; `end-to-end-report/report/09_rpe1_self.md`.

**Citations:** `Francis2025`, `Zuo2021`, `Patel2019`, `bouquet_BhattTERBEvolution2020`.

**Approx word count:** 250

---

### P11: Discussion — a feedback loop, "concerted evolution" in the loose sense, and the limits of inference

**Claim(s):**
- The five lines of evidence close a four-link causal loop (ED8a): **sequence sharing → 3-D proximity at meiotic bouquet → ectopic recombination → new shared segments → propagation**. Steps 1–2 measured directly (CHM13 Mantel ρ = 0.66, HG002 Pore-C 0.485); step 3 measured directly in pedigree (16 crossover-like + 133 gene-conversion-like events); steps 4 inferred from the population partition.
- The title's "concerted evolution" is used in the *loose* sense: ongoing inter-chromosomal recombination exchange that homogenises sequence between non-homologous chromosomes; the pedigree results provide the direct empirical anchor.
- Limitations: (i) bulk Hi-C is somatic LCL not germline meiotic — the 3-D signal is consistent with envelope tethering but does not directly measure bouquet contacts; (ii) no genome-wide LAD / Lamin B1 overlay in matched cell types yet — Dip-C radial position is the current proxy (C1 D4Z4 radial 0.732 peripheral; C14 PAR2 0.840; C10 chr17p 0.474 interior); (iii) Lalli 2025 cM/Mb anti-correlation with cross-arm affinity (ρ = −0.35, n = 46) collapses to ρ ≈ 0 (n = 40) once seven low-callability arms are excluded — short-read recombination maps cannot resolve PHRs; (iv) 12 % wfmash sampling is justified by Erdős–Rényi connectivity but rests on the chosen identity threshold (95 %).
- Outlook: long-read recombination maps + matched LAD/HiCAR for germline cell types + extending the pedigree set to the CEPH1463 4-gen at T2T quality will close the remaining open links.

**Evidence:** **ED8a** (causal loop schematic); **ED8b** (D4Z4 model); **ED8c** (Lalli 2025 anti-correlation null); **ED8d** (compartment diagnostic); `end-to-end-report/report/07_integrated.md` §"Causal feedback loop", §"D4Z4-CTCF-Lamin model"; `end-to-end-report/report/10_limitations.md` (18 enumerated limitations); `end-to-end-report/report/12_literature.md` (27 novel contributions + 7 testable predictions ledger); CROSSWALK §3 row "Mouse meiotic Hi-C", §4 C7 row, §4 C8 row.

**Citations:** `concerted_evolution_nahr_Arnheim1980`, `concerted_evolution_nahr_Ohta1984`, `concerted_evolution_nahr_Charlesworth1994`, `concerted_evolution_nahr_SamonteEichler2002`, `concerted_evolution_nahr_Eichler2001`, `concerted_evolution_nahr_Hillis1991`, `concerted_evolution_nahr_Vollger2023`, `Vollger2023`, `Sharp2006`, `Sasani2019`, `Lalli2025`, `Logsdon2021`, `Logsdon2024`, `Smolka2024`, `hprc_hprcv2_2025`, `Liao2023`.

**Approx word count:** 480

---

### Total main-text budget

- P1 300 + P2 350 + P3 350 + P4 380 + P5 240 + P6 420 + P7 420 + P8 360 + P9 460 + P10 250 + P11 480 = **4,010 words** (inside the 3,000–4,500 target).

> If the editor pushes back on length, P1 / P5 / P10 are the recommended trim points; P11 (Discussion) is load-bearing for the title claim and should not be cut below 400.

---

## Figure ↔ paragraph map

> Validation requires this table to reference only existing figure directories. The figure dirs that exist are: `fig1`, `fig2`, `fig3`, `fig4`, `ed1`, `ed2`, `ed3`, `ed4`, `ed5`, `ed8`, `nj_tree_arms`. **`ed6` and `ed7` do not exist and are not referenced anywhere below.**

### Main figures (Fig 1–4)

| Figure | Paragraphs that call it | Purpose |
|---|---|---|
| **Fig 1** (`paper_prep/figures/fig1/`) | P2 (a, b), P3 (a, b), P4 (c), P5 (d) | Landscape + arm-level partition + architecture categories |
| **Fig 2** (`paper_prep/figures/fig2/`) | P6 (a, b, c, d) | Within-community heterogeneity + two-domain + population F_ST |
| **Fig 3** (`paper_prep/figures/fig3/`) | P7 (a, b), P8 (c, d) | 3-D proximity mirrors sequence + flanking paradox |
| **Fig 4** (`paper_prep/figures/fig4/`) | P9 (a, b), P10 (c, d) | Pedigree-resolved exchanges + cross-species generalisation |

### Extended Data

| ED figure | Paragraphs that call it | Purpose |
|---|---|---|
| **ED1** (`paper_prep/figures/ed1/`) | P2 (a — pipeline), P3 (b — per-arm flank counts; c — PHR length distribution); P11 (d — chr18q chimera removal QC) | Pipeline + flank inventory + length distribution + QC |
| **ED2** (`paper_prep/figures/ed2/`) | P4 (a, d — sequence-level UMAP + arm/seq confusion); P6 (b — within-community bimodality); P5 (c — cross-arm affinity chord) | Sequence-level 50-community detail |
| **ED3** (`paper_prep/figures/ed3/`) | P3 / P6 (a — TAR1 prevalence; b — (TTAGGG)n island length; c — telomere length per community; d — TAR1 distance from telomere) | Annotation: TAR1 + ITS + telomere length per community |
| **ED4** (`paper_prep/figures/ed4/`) | P3 / P11 (a — PHR-only GO BP; b — copy-weighted GO; c — high-copy gene families incl. DUX4 16–20 copies; d — OR4F pseudogene gradient) | Gene-enrichment + copy-number diversity (DUX4 / OR4F gradient supports P4 4q–10q narrative) |
| **ED5** (`paper_prep/figures/ed5/`) | P7 (a — multi-resolution W/B; b — exclusion-controls Mantel; c — O/E within vs between; d — 15-community × 11-dataset reproducibility) | Multi-resolution + exclusion-confound robustness for Hi-C |
| **ED8** (`paper_prep/figures/ed8/`) | P11 (a — causal loop; b — D4Z4-CTCF-lamin; c — recomb anti-correlation honest null; d — compartment diagnostic) | Discussion synthesis: feedback loop + D4Z4 + recomb null + compartment |
| **nj_tree_arms** (`paper_prep/figures/nj_tree_arms/`) | P4 | Headline NJ cladistic tree (matches abstract wording) |

> **CROSSWALK conflict to flag for write-nature-article:** CROSSWALK §5 recommends ED4 be scrapped (no FDR-significant Fisher enrichments). This outline retains ED4 in a *softer* role — specifically panel (c) (DUX4 / OR4F high-copy families) supports the P4 4q–10q DUX4 copy-number-diversity claim from the abstract, and panel (d) supports the P11 pseudogenisation gradient. If T3 decides ED4 cannot carry that load, demote panels (a) and (b) to SI and keep only (c, d) in main ED, or move ED4 entirely to SI and inline the DUX4 copy-number number directly in P4.

---

## Methods outline

Subheadings (Nature Methods style, brief + ~1-line content):

1. **Sample selection and reference frame.** 233 HPRC v2 v1.1 individuals → 465 haplotype assemblies (~3 Gb each) + CHM13v2.0 reference (466 total). Per-superpop: AFR 67, EAS 52, AMR 44, SAS 37, EUR 33. (Cite `hprc_hprcv2_2025`, `Nurk2022`.)
2. **Telomere-anchored 500 kb flank extraction.** 18,827 flanks across 48 arms; arms with < 1 Mb contig length excluded; per-arm flank counts in ED1b; pq-classification table at `pq-classification/contig_classifications.tsv`.
3. **wfmash all-vs-all alignment.** `wfmash v0.23.0-41 -p 95 -t 48 --quiet`; each flank serves as target in turn. (Cite `Guarracino2023`.)
4. **Implicit pangenome graph + IMPG transitive closure.** All-vs-all PAFs are the implicit graph; IMPG `query -x` computes reachability. Justification of "no chromosomal partitioning": sampling rate ≈ 12 % ≫ Erdős–Rényi threshold p* = log(n)/n ≈ 5.2 × 10⁻⁴ for n = 18,827, ≈ 230× margin. (Cite `pangenome_graphs_impg_GarrisonGuarracino2023`, `pangenome_graphs_impg_GuarracinoHeumos2022`, `pangenome_graphs_impg_IMPG2023`, `pangenome_graphs_impg_Hickey2024`.)
5. **PHR detection.** impg sliding window: identity ≥ 95 %, ≥ 2 different chromosomes, ≥ 5 alignments per chromosome, output ≥ 3 kb; result: 15,668 PHRs (83.2 %), 41/48 arms.
6. **Pangenome graph + Jaccard similarity.** `pggb -p 95 -D /scratch`; `odgi similarity --all -P` over 15,668 sequences → 15,668 × 15,668 Jaccard matrix. (Cite `Garrison2024pggb`, `Guarracino2023`.)
7. **Arm-level distance matrix.** 41 × 41 matrix computed by averaging pairwise sequence-level Jaccard distances per arm pair (`hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`).
8. **Community detection — Leiden.** Arm-level: edge weight w_ij = exp(−d_ij / median(d)); resolution scan 0.1–3.0 step 0.01; selected 1.16 (mean silhouette 0.347); 15 communities. Sequence-level: k-NN graph with k ∈ {10, 25, 50, 75, 100, 125}, resolution 0.1–3.0; selected k = 75, res = 0.8 (modularity 0.97, silhouette 0.602); 50 communities. (Cite Leiden algorithm reference — **GAP: no Traag 2019 key in REFERENCES_v3.bib; needs to be added** — see Open Issues.)
9. **UPGMA dendrogram.** `hclust(..., method = "average")` on the 41 × 41 distance matrix; k = 14 (max silhouette 0.342); agreement with Leiden 14/15 (1c) / 12/15 (older 15-community cut).
10. **Neighbour-joining tree.** `ape::nj()` on the same 41 × 41 distance matrix; rooted at MRCA of acrocentric short-arm clade; 1,000-replicate perturbation bootstrap (Gaussian noise at σ = 25 % of off-diagonal IQR). (Cite **GAP: ape package; Saitou & Nei NJ original — neither key in REFERENCES_v3.bib; needs to be added** — see Open Issues.)
11. **Heterogeneity tests.** Wilcoxon paired (allele vs paralog); two-domain Spearman per arm; piecewise regression with changepoint detection; Fisher 2 × 5 superpopulation enrichment + BH; Hudson F_ST.
12. **Hi-C / Pore-C / CiFi pipeline.** mcool inputs at 5/10/20/50/100 kb; per-haplotype mode keeps maternal / paternal as separate arms; MAPQ filters disabled (multi-mappers retained, one random alignment per read); within / between (W/B) ratio bootstrap (10,000 perms); Mann–Whitney global test; Mantel Spearman (10,000 row/column perms); O/E inter-chromosomal normalisation. (Cite `hic3d_dixon2012`, `hic3d_imakaev2012`.)
13. **Exclusion controls.** Five sets (no acrocentric p / no sex / no acro p+sex / no all-acro+sex / no strongest) at all 5 mcool resolutions; tabulated in ED5.
14. **Single-cell 3-D.** Dip-C (Tan 2018, 16 GM12878 cells) remapped to T2T-CHM13v2.0; sperm scHi-C (Xu 2025, 20 cells). Negative control `S_all`: pooled 7 zero-signal arms.
15. **Mouse pipeline.** B6 + CAST T2T (Francis 2025); telomere-anchored flanks; 1 / 2 / 4 Mb window scan; 2-community partition; per-PHR-pair Jaccard vs Zuo 2021 zygotene Hi-C (lepto/zygo/pachy/diplo). (Cite `Francis2025`, `Zuo2021`.)
16. **Pedigree odgi-untangle.** `odgi untangle nth.best=1` per flank; HQ filter: minimum patch score 0.95 / 0.95 alignment; within-Leiden filter applied for credibility; pattern classification (`gene_conversion_like` / `crossover_like` / `acros_like` / `sandwich_same_hap` / `complex`) via `scripts/pedigree/analyze-pedigree-recombination.py`. (Cite `Cechova2025`, `Porubsky2025`.)
17. **Cross-assembler intersection (CEPH1463).** hifiasm AND verkko within-community parent × chr-pair features; 11 features pass.
18. **Software versions.** wfmash 0.23.0-41-gb5f0ff1c; impg commit 5b96025; pggb / odgi bundled; samtools / bedtools / bgzip via conda; hicexplorer 3.7.4; R packages: `ape`, `vegan`, `cluster`. (See `end-to-end-report/report/13_appendix.md`.)
19. **Data and code availability.** GitHub `ekg/phrs`; on-disk roots `/moosefs/guarracino/HPRCv2/PHR_III/` and `/moosefs/erikg/phrs/`; CHM13-coordinate PHR BEDs at repo root (`chm13.phrs.bed`, `chm13.phrs.no_acro.bed`, `CHM13-HG002.sub-telo-phrs.bed`).
20. **Limitations.** Eight items distilled from `end-to-end-report/report/10_limitations.md`: identity threshold; flank size truncation; small Hi-C N; somatic vs meiotic context; LCL somatic exchange caveat; multi-mapping caveat (addressed by flanking control); confound controls applied; recombination-rate test underpowered.

---

## Open issues / gaps

### Citation gaps — keys NOT in REFERENCES_v3.bib (must be added before T3 prose draft)

Verified via `grep -F '@' paper_prep/synthesis/REFERENCES_v3.bib | sed 's/^@[^{]*{//; s/,.*//' > /tmp/known_keys` (295 entries).

1. **Leiden algorithm primary reference.** PPTX slide 11 cites Traag, V. A., Waltman, L. & van Eck, N. J. *From Louvain to Leiden* (Sci. Rep. 9, 5233, 2019). **No `Traag2019` / `LeidenAlg2019` / `traag_leiden_2019` key exists in `REFERENCES_v3.bib`.** Required by Methods §8.
2. **Neighbour-joining primary reference.** Saitou & Nei (1987) *Mol. Biol. Evol.* 4 (4): 406–425. **No `SaitouNei1987` / `NeighborJoining1987` key exists.** Required by Methods §10 + P4.
3. **`ape` R package reference.** Paradis et al. 2004 *Bioinformatics* / 2019 update. **No `apePackage` / `Paradis2004` key exists.** Required by Methods §10.
4. **Erdős–Rényi connectivity reference.** Erdős & Rényi 1959/1960. **No `ErdosRenyi1960` key.** Required by Methods §4 + P2 (~12 % sampling argument).
5. **wfmash primary reference.** PPTX slide 36 refers to "Li and Rong 2020" (implicit interval tree); Methods cites wfmash but **no standalone `wfmash` / `Li2020` / `LiRong2020` key exists.** Required by Methods §3. (Note: `Guarracino2023` is present — that may already be the wfmash methods reference; verify which paper it points to.)
6. **odgi reference.** No `odgi` key; `Guarracino2023` may cover it but unclear. Required by Methods §6.
7. **HG002 / GIAB reference.** `Zook2020` exists but a Pore-C / CiFi data-source reference per dataset is unclear; verify if the on-disk Pore-C / CiFi data has a published reference distinct from `hic3d_cifi2025` / `Ulahannan2019`.

**Action for T3:** before drafting prose, run a bib gap pass — either add the missing keys to `REFERENCES_v3.bib` or substitute existing keys (note: `Guarracino2023` may already cover wfmash / impg / odgi; `pangenome_graphs_impg_GarrisonGuarracino2023` may already cover the implicit pangenome graph concept; verify per-cite). The outline does NOT invent any bibtex key — every citation above is verified present except where flagged "GAP" inline.

### PPTX vs paper_prep conflicts

1. **NJ tree vs Leiden + UPGMA (resolved).** PPTX slides 5–12 describe Louvain / community detection + UPGMA only; both abstracts (BoG + Nature) name "neighbour-joining trees". Outline uses NJ as headline (matches abstract; `nj_tree_arms/` exists on disk), Leiden + UPGMA as supporting; CROSSWALK §2 C5 reconciles.
2. **PPTX slide 38 says "all-to-all alignment of 465 × 3 billion bp implies a 2 septillion cell matrix" (treats as full all-to-all)** vs both abstracts which say ~12 % pairwise sampling. CROSSWALK §7b reconciles: the all-vs-all framing is the *target*; wfmash's k-mer prefiltering realises ~12 %; Erdős–Rényi connectivity makes the 12 % sampling sufficient. Outline P2 + Methods §4 both state this.
3. **PPTX slide 42 says "11.6 % with 466 haplotypes" vs both abstracts and outline say 12 % / 465 + CHM13.** Numerical agreement (11.6 % ≈ 12 %); 466 vs 465 is the CHM13-inclusion convention. Methods §1 nails the convention.
4. **PPTX slide 22 comment from Erik Garrison flags "need metric for the degree of similarity within communities vs between … degree of contact intensity".** Outline P7 (Fig 3a/b/c) delivers this: B/W ratios, per-arm-pair Spearman, the 14-test forest plot. Resolved.
5. **PPTX slide 26 comment from Erik Garrison flags "we need to check if the untangling is actually too conservative, and if this means we have a whole phr that swapped".** Outline P9 frames 538 patches as an HQ + within-community *lower bound*. The 80-pp gap between WashU (92 %) and CEPH1463 (12–13 %) within-Leiden rate is reported as evidence that the filter is conservative against fragmentation noise, not against true exchange. **Open follow-up (not blocking outline):** quantify whether any single PHR was wholly swapped between PAN027 and PAN028 — requires a length-of-swapped-PHR-block analysis on the existing patches.tsv files; not currently in the figure set.
6. **PPTX slides 19 + 31 mention mouse and CHM13 Hi-C 3D MDS in supplementary** (extended-data territory). Outline uses ED5 + ED8 only for these; mouse is in P10 / Fig 4d. Consistent.
7. **PPTX slide 14 Erik comment: "put all communities one per page in supplementary figures".** Outline does not include a 15-community deep-dive panel set; **recommend** that T3 confirm with PI whether a separate ED figure (could be `ed6` or `ed7` if we create them, or appended panels in SI) showing per-community detail is required. **CROSSWALK note:** `ed6` and `ed7` directories do not currently exist; creating them is outside this outline's scope.

### Numbers not directly evidenced in paper_prep (require T3 verification or numerical fill)

1. **"comparable in scale to PAR2" (P3 / abstract).** No side-by-side PAR2 length comparison panel exists; abstract draws the analogy but the PHR length distribution panel (ED1c) does not overlay PAR2's ~334 kb anchor. **Recommend:** T3 add the numeric comparison in prose (median PHR 105 kb vs PAR2 ~334 kb), or T3 + downstream task adds a small inset / annotation to ED1c.
2. **"~12 % pairwise sampling" derivation from on-disk PAFs.** Per CROSSWALK §7b, this number must be derived from the on-disk PAFs (#evaluated pairs / C(18,827, 2)). The figure 12 % is in both abstracts; the on-disk-derivation step is not yet committed to a TSV. **Recommend:** T3 either inserts a placeholder "≈ 12 % (recomputed from PAF set)" and flags the gap, or a follow-up `wg add` task computes the exact realised rate.
3. **"Andrea uses 465; abstract uses 466"** (CROSSWALK §2 C3). Outline resolves: 466 = 233 × 2 + CHM13. Methods §1 nails the convention; numerical reconciliation table immediately below.

### Number reconciliation (single source of truth for prose)

| Quantity | Value used in outline | Source of disagreement |
|---|---|---|
| Individuals | 233 | consistent across PPTX, BoG abstract, Nature abstract, CROSSWALK, fig captions |
| HPRC haplotypes | 465 | consistent (Andrea report, fig captions, CROSSWALK) |
| With-reference total | 466 (= 465 + CHM13) | Nature abstract uses 466; PPTX slide 42 also uses 466; Andrea's report uses 465. Methods §1 reconciles. |
| Telomere-anchored flanks | 18,827 | consistent |
| PHRs | 15,668 | consistent (ED1d notes 15,669 → 15,668 after chr18q chimera removal) |
| Signal-bearing arms | 41 / 48 | consistent |
| Arm-level Leiden communities | 15 | consistent |
| Sequence-level Leiden communities | 50 | consistent |
| Median PHR length | 105 kb | consistent |
| Mean PHR length | 144 kb | consistent |
| HG002 Hi-C B/W @ 50 kb | 0.027 (p = 4.0 × 10⁻⁶⁶) | consistent |
| CHM13 Hi-C B/W @ 50 kb | 0.071 (p = 6.0 × 10⁻¹⁸) | consistent |
| HG002 Pore-C B/W @ 50 kb | 0.056 (p = 3.9 × 10⁻⁸⁵) | consistent |
| Mantel ρ (CHM13 Hi-C) | 0.66 (abstract) / 0.656 (CROSSWALK §1 ch05 row) / 0.674 (CROSSWALK §6 ED8a) | Three slightly different rho values appear in source documents (0.656 / 0.66 / 0.674); outline uses **0.66** as the rounded headline. **Recommend:** T3 confirm which is the canonical per-arm-pair vs Mantel ρ; both are reported in chapter 05 and should not be conflated. |
| Mantel ρ (HG002 Hi-C) | 0.66 (abstract) / 0.657 (CROSSWALK §1 ch05 row) / 0.485 (ED8a — for HG002 *Pore-C*) | The 0.485 figure is HG002 *Pore-C* not Hi-C; outline P11 uses both. |
| Mouse zygotene per-PHR-pair ρ | 0.715 (n = 344, p = 4.4 × 10⁻⁵⁵) | consistent across abstract, Fig 4d caption, SURVEY_08 §1.7 |
| WashU pedigree HQ patches | 538 | consistent |
| WashU within-community fraction | 92 % (494 / 538) | consistent |
| WashU gene-conversion-like | 133 | consistent |
| WashU crossover-like | 16 | consistent |
| WashU acros_like | 229 | consistent |
| CEPH1463 cross-assembler-validated features | 11 | consistent |

---

*End of NATURE_DRAFT_OUTLINE.md. T3 (`write-nature-article`) should now expand each paragraph header into prose, respecting the word budgets, evidence chains, and citation keys above, and consulting the Open Issues section before adding any new citation.*
