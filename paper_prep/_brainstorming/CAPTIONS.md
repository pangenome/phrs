# CAPTIONS — Manuscript display items (deduplicated)

**Origin.** Captions distilled from `paper_prep/synthesis/MANUSCRIPT_SKELETON.md` (preliminary captions for 4 main + 8 ED figures) and from per-survey caption stubs (`SURVEY_07 §6 T-1..T-7`, `SURVEY_FIG_inventory.md`). One caption per figure. Numerical values pulled from the named-TSV evidence tables in the cited surveys.

**Convention.**
- *p-values* are uncorrected unless suffixed `(BH-q)`. Multi-test corrections are tracked separately in `STATS_AUDIT.md` (task `validate-statistics-fdr-cis`).
- *n* refers to the unit being counted (haplotypes / pairs / cells / arms) — explicit in each caption.
- *Effect size* = the test statistic listed in the parent survey: B/W ratio (bulk Hi-C / Pore-C; values <1 = within-community contacts more frequent than between), W/B ratio (Dip-C / sperm scHi-C 3D distances; values <1 = within closer), Spearman ρ (Mantel + per-pair correlations), Wilcoxon paired test, Fst, Fisher OR, silhouette s, paralog-vs-allele conversion score.
- *Schematic-only* panels (ED1a, ED8a, ED8b) carry no statistic but cite the underlying source numbers in the caption.

---

## Figure 1 — Population-scale subtelomeric communities (the landscape)

*Inter-chromosomal subtelomeric sharing partitions 41 chromosome arms into 15 communities across 465 HPRCv2 haplotypes from 233 individuals.* **(a)** Genome-wide stacked identity heatmap (n = 465 haplotypes × 24 chromosomes) shows discrete blocks of telomere-anchored shared sequence; 15,668 PHRs called from 18,827 telomere-anchored 500 kb flanks (`SURVEY_01 §1.1–§1.5`). **(b)** Genome-wide n-chromosomes-sharing heatmap with PHR-call BED overlay; 83 % of flanks (n = 15,621/18,827) carry inter-chromosomal sequence sharing. **(c)** 41 × 41 arm-level Jaccard distance heatmap (Leiden, k = 15; modularity Q = 0.97; silhouette s = 0.347 arm-level / 0.602 sequence-level; 12/15 UPGMA-Leiden agreement, `SURVEY_01 §1.8, §1.10`). **(d)** Per-arm architecture-category bar — homogeneous 8/41, polymorphic 34/41, fully interchangeable 7/41 (acrocentric p-arms + PARs). 15.9 % of subtelomeric sequences cross-arm at the arm level; 11.1 % at the sequence level (n = 15,668 sequences; `SURVEY_04 §1.4, §1.11`).

---

## Figure 2 — Within-community heterogeneity, the two-domain model, and population history

*Inter-chromosomal exchange leaves a quantitative population signature.* **(a)** Allele-vs-paralog Wilcoxon paired test across 9 multi-arm communities (n = 5,946 paired distances). Allele closer than paralog in 8/9 communities; overall Wilcoxon p < 1 × 10⁻³⁰⁰; C7 (acrocentric p-arms) reversed — 70.5 % paralog closer (Wilcoxon p = 2.0 × 10⁻⁷); C7 silhouette = −0.029, conversion score = 1.000 toward all four other acrocentric p-arms (`SURVEY_04 §1.1`). **(b)** Per-arm Spearman gradient and piecewise breakpoint for the two-domain model (n = 41 testable arms): 39/48 arms with significant Spearman gradient; 39/41 prefer two-segment over single-segment; gradient detectable in 99.7 % of individual haplotypes; arm-specific breakpoints 15 kb (chr22q) – 445 kb; internal (TTAGGG)n islands within 25 kb of breakpoint on 11/19 testable arms (`SURVEY_04 §1.3`). **(c)** Cross-arm superpopulation enrichment heatmap with Fst overlay (n = 19 cross-arm pairs × 5 superpopulations); 10/19 Fisher BH-significant; AFR vs non-AFR mean Fst = 0.10–0.15; non-AFR/non-AFR Fst ≈ 0 (`SURVEY_04 §2.1`). **(d)** Out-of-Africa population tree from cross-arm-affinity frequencies (n = 9 arms × 5 superpopulations); recovers AFR-deepest topology with AMR–EUR closest pair (`SURVEY_10/11/12 §6 T-4`).

---

## Figure 3 — Three-dimensional nuclear organisation mirrors sequence communities

*Sequence-defined communities are physical across six technologies and two species.* **(a)** HG002 Pore-C inter-chromosomal contact matrix at 50 kb, arms ordered by sequence community (n = 41 arms; B/W = 0.056, p = 3.9 × 10⁻⁸⁵; CHM13 per-pair Spearman ρ = 0.674; `SURVEY_05 §1, §2.1`). **(b)** Convergent-evidence forest plot — 14 independent tests across bulk Hi-C, Pore-C, CiFi, Dip-C, sperm scHi-C, and mouse meiotic Hi-C; all 14 effect sizes on the within-community-closer side of unity; bootstrap 95 % CIs (sperm W/B = 0.401, p = 3.9 × 10⁻⁵¹, n = 20 cells; mouse zygotene per-pair ρ = 0.715, p = 4.4 × 10⁻⁵⁵; `SURVEY_07 §6 T-1`). **(c)** S_all negative-control panel for the 7 non-sharing arms: 0/16 GM12878 Dip-C and 1/20 sperm cells fall below W/B = 1, vs 16/16 + 19/20 for sharing communities; non-sharing arms are 11 % (GM12878) / 40 % (sperm) *farther* in 3D than between-community pairs (`SURVEY_06 §1.3`). **(d)** Flanking-vs-PHR B/W panel across 8 bulk Hi-C/Pore-C samples — flanks 100 kb centromere-ward of PHR are unique sequence yet give a stronger signal (HG002 flanking B/W = 0.002 vs PHR B/W = 0.027; 13×); inset: Dip-C flanking radial 0.503 vs non-flanking terminal 0.551 (p = 7.4 × 10⁻³⁵; `SURVEY_05 §6.1, SURVEY_06 §1.4, SURVEY_07 §6 T-5`).

---

## Figure 4 — Pedigree-resolved exchanges and cross-species generalisation

*Direct observation of the events that build the communities, in human pedigrees and across mammals.* **(a)** WashU 3-generation untangle ribbons (PAN027 maternal hap1 from PAN010, PAN028 maternal hap1 from PAN027) on T2T assemblies (Cechova et al. 2025); off-diagonal inter-chr patches highlighted; 494 / 538 (92 %) HQ inter-chr patches sit in a Leiden community vs 12–13 % at fragmented CEPH1463 (`SURVEY_14 §1.1`). **(b)** CEPH1463 4-generation cross-assembler-validated parent-feature matrix (hifiasm + verkko intersection; n = 11 robust events across 28 samples per Porubsky et al. 2025); chr10/chr18 (community C2 / Linardopoulou pair) detected independently in NA12877 paternal and NA12878 maternal — same canonical exchange in unrelated individuals (`SURVEY_14 §1.6`). **(c)** RPE-1 t(X;10) translocation rediscovery — pipeline run on a single diploid individual partitions chr10_q with chrX_q (community C2; `SURVEY_09 §6 T-1`). **(d)** Mouse cross-species generalisation — zygotene per-PHR-pair similarity vs Hi-C contact (n = 27 arms × 27 arms; ρ = 0.715, p = 4.4 × 10⁻⁵⁵, Mantel ρ = 0.718); same correlation across all 4 meiotic stages (B/W 0.029–0.122, leptotene → diplotene); B6 + CAST T2T from Francis et al. 2025; Hi-C from Zuo et al. 2021 (`SURVEY_08 §1.7, §6 T-1`).

---

## Extended Data Figure 1 — Pipeline and per-arm flank inventory

*From assemblies to communities.* **(a)** Pipeline schematic (no statistic): 465 HPRCv2 haplotype-resolved assemblies → 18,827 telomere-anchored 500 kb flanks → wfmash all-vs-all (asm20, p95, id95, len ≥ 30 kb) → impg projection → Leiden community detection. **(b)** Per-arm flank counts across 48 arms with assembly QC overlay (n = 18,827 flanks; assembly classifications from `contig_classifications.tsv`; `SURVEY_01 §3`). **(c)** PHR length distribution (n = 15,668 PHRs; median 105 kb, mean 144 kb). **(d)** Chr18_q (NA18982#1) chimera evidence — wfmash + minimap2 dotplot, NNN gap, and Flagger annotation flagging the chimera (`SURVEY_01 §1.5, §5 item 6`).

---

## Extended Data Figure 2 — Sequence-level (50-community) detail

*Sequence-level partition and confusion with the arm-level structure.* **(a)** UMAP / force-directed layout of 15,668 PHR sequences coloured by 50-community Leiden partition (silhouette s = 0.602; modularity Q = 0.97; `SURVEY_01 §1.10`). **(b)** Within-community Jaccard distance bimodality across 8 multi-arm communities (C1, C2, C3, C5, C6, C7, C11, C12); separable allelic vs paralogous peaks for C2 and C12 (`SURVEY_04 §1.10`). **(c)** Cross-arm affinity circular plot — 41 arms with edges weighted by absorbed sequences (n = 1,740 cross-arm sequences arm-level; n = 2,484 at the broader cross-arm level; 11.1 % / 15.9 % rates; `SURVEY_01 §6 F5`). **(d)** Confusion matrix of Arm-Leiden vs Sequence-Leiden assignments (15 × 50; ARI = 0.35; NMI = 0.76).

---

## Extended Data Figure 3 — Annotation: TAR1 + internal (TTAGGG)n + telomere length

*Subtelomeric annotation across 465 haplotypes.* **(a)** TAR1 prevalence per arm (PAR1 absence; acrocentric intermediate; autosomal saturation; `SURVEY_02 §6 Fig M1a`). **(b)** Internal (TTAGGG)n island length distribution and canonical-fraction histogram (n = 18,352 islands across 8,321 sequences = 53.1 % of PHR sequences; mode 50–74 bp; 32.2 % canonical, 47.2 % variant-dominant). **(c)** Terminal telomere length by community (n = 15 communities; Kruskal-Wallis H = 100.89, p = 3.2 × 10⁻¹⁵; `SURVEY_02 §6 ED3`). **(d)** Per-arm TAR1 positional distance-from-telomere (n = 41 arms; `SURVEY_02 §6 Fig M1b`).

---

## Extended Data Figure 4 — Gene enrichment, pseudogene gradient, and copy-weighted GO

*Gene content of PHRs.* **(a)** GSEA / GO:BP top terms — snRNP, olfactory, sensory perception (caveat: 1 Mb window; PHR-only re-run flagged as gap; `SURVEY_FIG_inv §3`). **(b)** Copy-weighted vs deduplicated comparison; olfactory family fold-enrichment = 598× when copies counted (`SURVEY_DATA §4`). **(c)** High-copy gene families bar — DUX4 ×18 (community C1), BAGE2, MTCO, RPL23A, SEPTIN14P22, OR4F (n = 6 families; `SURVEY_DATA §2`). **(d)** OR4F pseudogenisation gradient across 9 arms — 62.1 % overall pseudogene rate; 11.1 % at chr7_p → 99.8 % at chr15_q (`SURVEY_10/11/12 C12`).

---

## Extended Data Figure 5 — Multi-resolution + confound-exclusion robustness for Hi-C

*Robustness across resolutions and confound exclusions.* **(a)** B/W ratio across 5 mcool resolutions (5/10/20/50/100 kb) for 8 bulk Hi-C / Pore-C datasets (n = 8 samples × 5 resolutions = 40 tests; per-test p-values from `*_global_test.tsv`; all p < 0.01; `SURVEY_05 §6.1`). **(b)** Acrocentric / sex / strong-community exclusion grid — Mantel ρ before/after for 5 exclusion conditions × 8 samples; ρ remains within ±0.05 of full-data values (no exclusion *weakens* the signal; `SURVEY_05 §6.2 SI 5.5`). **(c)** O/E-normalised within vs between contact (n = 8 samples; 8.6× – 34.4× separation; `SURVEY_05 §1.10`). **(d)** Per-community reproducibility heatmap (15 communities × 10 datasets; * = BH-q < 0.05; `SURVEY_05 §6.2 SI 5.8`).

---

## Extended Data Figure 6 — Single-cell 3D (Dip-C + sperm) full panel + RPE-1 self

*Single-cell 3D validation.* **(a)** GM12878 Dip-C Mantel scatter (n = 16 cells) + radial-by-community panel; Dip-C method from Tan et al. 2018 (`SURVEY_06 §2`). **(b)** Sperm scHi-C Mantel scatter + per-cell W/B + radial (n = 20 cells, 10 X-bearing + 10 Y-bearing; Xu et al. 2025; W/B = 0.401, p = 3.9 × 10⁻⁵¹). **(c)** Per-arm radial positions across 46 arms — S1/S2/S7 peripheral; cis-arm proximity for non-sharing arms (`SURVEY_06 §1.4`). **(d)** RPE-1 self-discovered cell-cycle modulation: async vs mitotic across 3 datasets and all 5 resolutions; mitotic 3× stronger global W/B but 1.4× weaker per-arm-pair ρ (`SURVEY_09 §1.4`).

---

## Extended Data Figure 7 — Cross-species mouse meiotic generalisation

*Mouse meiosis recovers the same architecture.* **(a)** Mouse 27-arm Leiden 2-community map vs human 41-arm 15-community structure (architectural contrast; `SURVEY_08 §6 T-5`). **(b)** Mouse window-size scaling across 6 windows (1/2/4/10/15/33 Mb) — PHR size + saturation fraction (`SURVEY_08 §6 T-3`). **(c)** Mouse B/W & Mantel ρ across 4 meiotic stages × 5 resolutions × 3 windows (n = 60 cells/stage tests; B/W = 0.029–0.122; ρ = 0.61–0.72 across leptotene → diplotene; Hi-C from Zuo et al. 2021). **(d)** Mouse private-pair → human syntenic-net (mm39 → hg38; all in human interiors; `SURVEY_08 §1.6`).

---

## Extended Data Figure 8 — Discussion synthesis: feedback loop, recombination null, mechanistic models

*Discussion synthesis.* **(a)** Causal feedback loop schematic (4 links coloured by support level — direct, established literature, inferred, untested). **(b)** D4Z4-CTCF-lamin tethering schematic for community C1 — D4Z4 macrosatellite → CTCF (Ottaviani et al. 2009) → lamin A/C (Masny et al. 2004) — with 0–15 kb sharing-peak inset and DUX4L copy number (median 22, 0–2 outliers; p = 5.3 × 10⁻⁶ for outlier deviation; FSHD context Lemmers et al. 2010). **(c)** Recombination map vs cross-arm affinity per arm (Lalli et al. 2025 T2T-CHM13 cM/Mb; n = 39 arms full ρ = −0.43, p = 0.006; n = 32 after excluding 7 acrocentric/PAR arms with 0–12 callable variants ρ = 0.00, p = 0.98 — the explicit honesty figure for the Discussion). **(d)** Compartment-identity-at-tips diagnostic (n = 92 arm × haplotype HG002 tips; e1 distribution 68 % A; mean e1 = +0.007; tips A-leaning by GC but interior-positioned by Dip-C radial — telomere clustering, not lamina association, drives end alignment).

---

## Caption conventions and gaps

- **Effect-size convention** is included in each caption (B/W or W/B for 3D distances, ρ for correlations, Wilcoxon p for paired distance comparisons, Fst for population structure). For pure schematic panels (ED1a, ED8a, ED8b) no statistic is associated; the underlying TSV path is cited instead.
- **Headline numbers** (n = 233 individuals, n = 465 haplotypes, n = 18,827 flanks, n = 15,668 PHRs, n = 41 arms / 15 communities / 50 sequence-communities) are quoted in the Results main text and in `MANUSCRIPT_SKELETON.md` headline-numbers section, not repeated in every caption.
- **Multi-test correction** state per p-value is the responsibility of `validate-statistics-fdr-cis` (task 10) — captions here flag uncorrected p-values where the underlying survey did so.
- **Schematics** (ED1a pipeline; ED8a feedback loop; ED8b D4Z4 tethering) intentionally carry no test statistic — the validation task `validate-captions-references` flags these as schematic-only and not as missing statistics.
