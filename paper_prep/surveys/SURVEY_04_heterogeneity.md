# Survey 04 — Within-arm heterogeneity

Source: `end-to-end-report/report/04_heterogeneity.md`
Survey scope: extract findings, figures, CSVs, methods, gaps, suggested figures, and talk takeaways for the within-arm heterogeneity section.

---

## 1. Key findings (with metrics)

The section contains **eleven** sub-analyses on the same underlying question: are arm-of-origin signatures preserved within each multi-arm community, or does inter-chromosomal exchange erase them? The following are the load-bearing findings for the manuscript and the talk.

### 1.1 Allele vs paralog distance — arm identity is preserved in 8 of 9 multi-arm communities
- **Test:** Wilcoxon signed-rank, paired (allelic Jaccard distance vs paralog Jaccard distance) within each individual.
- **Overall:** 5,946 paired observations; median allelic 0.069 vs median paralog 0.387; only 20.4 % of pairs have a paralog closer than the allele; combined Wilcoxon p < 1e-300.
- **Eight communities significantly favour the allele** (C1–C6, C11, C12), Wilcoxon p between 6.0e-5 (C1, D4Z4) and 3.1e-159 (C6).
- **The lone reversal is C7 (acrocentric p-arms):** N = 156 pairs, median allelic 0.481 vs median paralog 0.353, **70.5 %** of individuals have a paralog closer than the allele, p = 2.0e-7. This is the population-scale signature of complete inter-chromosomal homogenization of acrocentric p-arms.
- **C1 (D4Z4) is the autosomal weak spot:** 44.2 % of pairs have a paralog closer; consistent with ongoing chr4_q ↔ chr10_q exchange.

### 1.2 Silhouette / arm separation — two-tier outcome across communities
- **Excellent separation (sil > 0.8):** C2 (chr10_p + chr18_p, sil = 0.888) and C12 (chr2_q + chr20_p, sil = 0.820) — arms share duplicon content but retain distinct sequence content.
- **Negative-silhouette outcomes — arms indistinguishable or anti-clustered:**
  - **C14 (PAR2: chrX_q + chrY_q), sil = −0.163** — X and Y sequences are more similar to each other than to same-chromosome sequences.
  - **C7 (acrocentric p-arms: chr13_p, chr14_p, chr15_p, chr21_p, chr22_p), sil = −0.029** — five-arm near-interchangeability.
  - **C15 (PAR1: chrX_p + chrY_p), sil = −0.025** — pseudoautosomal as expected.
- **C1 (D4Z4: chr4_q + chr10_q), sil = 0.147** — poor separation; supports the C1 reversal in §1.1.
- Single-arm communities (C8/C9/C10/C13) excluded from silhouette computation.

### 1.3 Two-domain subtelomeric model — Flint/Mefford supported at all three test levels
The most quantitatively load-bearing result in the section.

- **Test 1 (Spearman gradient on per-window n_chrs vs distance from telomere):** **39 / 48 arms** show a significant negative correlation; in the focus arm set: chr4p ρ = −0.85, chr4q ρ = −0.68, chr16p ρ = −0.32, chr18p ρ = −0.31, chr20p ρ = −0.47, chr22q ρ = −0.95. Per-sequence: **99.7 %** of haplotypes with ≥ 5 windows have negative within-sequence ρ (median ρ = −0.79; 13,728 of 13,763).
- **Test 2 (Piecewise vs single linear):** **39 / 41 testable arms** prefer the two-segment model. Focus-arm breakpoints / R² gain:
  - chr4p — 70 kb, R² 0.13 → 0.81 (+0.68)
  - chr4q — 50 kb, R² 0.14 → 0.92 (+0.77)
  - chr16p — 295 kb, R² 0.44 → 0.64 (+0.20)
  - chr18p — 120 kb, R² 0.05 → 0.30 (+0.24)
  - chr20p — 165 kb, R² 0.24 → 0.47 (+0.23)
  - chr22q — 15 kb, R² 0.81 → 1.00 (+0.19)
  - Two non-significant arms: chr12q (p = 0.48, already linear) and chr19p (p = 0.10).
- **Test 3 (Internal (TTAGGG)n co-localization):** Of the 19 arms with internal TTAGGG blocks > 5 kb from telomere, **11 / 19 arms** have an ITS block within 25 kb of the breakpoint, **16 / 19** within 50 kb. Best-aligned: chr8p (0.1 kb), chr20p (0.6 kb), chr13q (1.5 kb), chr7p (3.5 kb). Overall Spearman of breakpoint vs median ITS position is marginal (ρ = 0.42, p = 0.08) — colocalization is per-arm-specific, not universal.
- **Synthesis:** the pangenome-scale data extends Flint et al. 1997 / Mefford & Trask 2002 from a handful of arms to the whole human chromosome complement; breakpoints are arm-specific (10–445 kb), not at a fixed distance from the telomere.

### 1.4 Population structure of cross-arm exchange (Fisher + Hudson Fst)
- **15.9 % of all sequences (2,484 / 15,668) carry cross-arm affinity** at the arm level.
- **10 / 19 arm-community pairs show significant superpopulation bias** (Fisher exact, BH-corrected p_adj < 0.05). Strongest: chr4_q-C1, chr16_q-C3, chr6_p-C5, chrX_p-C15 (all p_adj = 4.7e-4).
- **Hudson Fst across the 10 strongest pairs averaged:** mean Fst = 0.044; AFR vs non-AFR Fst = 0.10–0.15; non-AFR vs non-AFR Fst between −0.05 and 0.01. Pattern parallels out-of-Africa demography; mechanism (population-specific exchange vs drift / ILS) cannot be distinguished without an outgroup.

### 1.5 Subtelomeric type discordance (heterozygosity for haplotype class)
- **9 / 35 arm-community pairs show > 20 % discordance.** Top: chr22_q-C6 47.5 %, chr22_p-C7 46.3 %, chr9_q-C3 45.2 %, chr4_q-C1 43.4 %, chr6_p-C5 41.7 %, chr11_p-C3 37.5 %.
- **chr4_q discordance 43.4 % is roughly twice the ~20 % from a Dutch cohort (Mefford & Trask 2002 citing van Deutekom 1996)** — but this work captures the broader PHR, so the comparison is not apples-to-apples.

### 1.6 Region-length dichotomy between cross-arm and self-arm
- **14 / 18 arm-community pairs differ significantly** (Wilcoxon rank-sum, BH-corrected p_adj < 0.05).
- **Largest effect:** chrX_p (C15) cross-arm 500 kb vs self 20 kb (25× longer; nearly all chrX_p clusters with chrY_p, the rare "self" haplotypes are truncated).
- **Cross-shorter examples:** chr22_q (C6) 20 vs 25 kb; chr6_p (C5) 100 vs 153 kb.
- **Cross-longer examples:** chr16_q (C3) 210 vs 155 kb; chr9_q (C3) 305 vs 225 kb.

### 1.7 Gene-repertoire replacement (cross-arm "conversion score")
- **Score 1.000 (complete gene-repertoire overlap with affinity arm)** on the acrocentric chr13_p (vs chr14_p / chr15_p / chr21_p / chr22_p) and at PAR1 (chrX_p ↔ chrY_p, 10 genes) / PAR2 (chrX_q ↔ chrY_q, 22 genes).
- **f7501 cluster:** chr11_p → chr19_p score 0.698 (75 affinity-specific genes); chr11_p → chr3_q 0.677.
- **D4Z4 (C1):** balanced bidirectional exchange (chr4_q → chr10_q 0.500; chr10_q → chr4_q 0.494).

### 1.8 TAR1 in cross- vs self-arm
- TAR1 prevalence > 85 % in both classes for most communities. Only 3 / 19 pairs differ significantly: chr16_q-C3 (100 % cross vs 71.1 % self, p_adj = 2.4e-10), chr13_p-C7 (88.1 % vs 0 %, p_adj = 1.6e-6), chr14_p-C7 (72.1 % vs 94.9 %, p_adj = 0.010). PAR1 is essentially TAR1-free (< 2 %).
- A causal role for TAR1 in exchange cannot be inferred from prevalence alone.

### 1.9 Internal TTAGGG island position is exchange-status-invariant
- Cross-arm proximal fraction = 45.9 %; self-arm = 45.9 %; χ² = 0.00, p = 0.99, OR = 1.00. Island distribution reflects duplicon architecture, not exchange history.

### 1.10 Within-community Jaccard distance is bimodal in well-separated communities
- C2 peaks ≈ 0.02 / 0.37; C12 peaks ≈ 0.04 / 0.73; bimodality is qualitatively consistent with Ambrosini et al. 2007. C1 / C7 are diffuse — direct visual confirmation of homogenization.

### 1.11 Cross-arm affinity per arm (top of the cross-arm-rate ranking)
- chrX_q (C14) 99.7 % toward chrY_q; chrX_p (C15) 93.3 % toward chrY_p; chr21_p (C7) 94.0 %; chr11_p (C3) 74.1 % — the highest autosomal rate, toward chr3_q / chr19_p / chr9_q (f7501).

---

## 2. Existing figures referenced

The report file does not embed image links. The following figures already exist on disk and align directly with the analyses above. PDF versions are listed; matching `.png` files exist alongside each `.pdf` in the same directory.

### 2.1 Within-arm heterogeneity panel suite
Path: `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/`

| File | Used for |
|---|---|
| `within_arm_heterogeneity_01_separation_overview.pdf` | Silhouette / arm separation across communities (§1.2) |
| `within_arm_heterogeneity_02_affinity_heatmap.pdf` | Per-arm cross-arm affinity (§1.11, §1.4) |
| `within_arm_heterogeneity_03_mds_zoom.pdf` | Within-community MDS (qualitative companion to §1.10) |
| `within_arm_heterogeneity_04_superpop_composition.pdf` | Superpop composition of cross-arm vs self-arm (§1.4) |
| `within_arm_heterogeneity_05_heterozygosity_heatmap.pdf` | Type discordance heatmap (§1.5) |
| `within_arm_heterogeneity_06_region_length.pdf` | Cross- vs self-arm length distributions (§1.6) |
| `within_arm_heterogeneity_07_gene_content.pdf` | Gene-content overlap (§1.7) |
| `within_arm_heterogeneity_08_conversion_evidence.pdf` | Conversion evidence figure (§1.7) |
| `within_arm_heterogeneity_09_conversion_score.pdf` | Cross-arm gene-conversion score (§1.7) |
| `within_arm_heterogeneity_10_affinity_network.pdf` | Cross-arm affinity network (§1.4, §1.11) |
| `within_arm_heterogeneity_11_tar1_prevalence.pdf` | TAR1 prevalence by class (§1.8) |
| `within_arm_heterogeneity_12_tar1_cross_vs_self.pdf` | Cross- vs self-arm TAR1 contrast (§1.8) |

### 2.2 Sequence-level subset (used in §1.2, §1.4, §1.6, §1.8)
Path: `/moosefs/guarracino/HPRCv2/PHR_III/sequence_level/heterogeneity/`
- `within_arm_heterogeneity_01_separation_overview.pdf`
- `within_arm_heterogeneity_04_superpop_composition.pdf`
- `within_arm_heterogeneity_06_region_length.pdf`
- `within_arm_heterogeneity_11_tar1_prevalence.pdf`

### 2.3 Companion community / MDS figures (context for §1.1–§1.2)
Path: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/`
- `hprcv2.1Mb.subtelo.heatmap.arm.clustered.pdf` / `.heatmap.arm.ordered.pdf`
- `hprcv2.1Mb.subtelo.heatmap.chrom.clustered.pdf` / `.heatmap.chrom.ordered.pdf`
- `hprcv2.1Mb.subtelo.mds.arm-leiden-k15.communities.pdf` / `.superpop-hulls.pdf`
- `hprcv2.1Mb.subtelo.mds.arm-upgma-k14.communities.pdf` / `.superpop-hulls.pdf`
- `hprcv2.1Mb.subtelo.umap.color-by-arm.pdf` / `.color-by-chromosome.pdf` / `.color-by-superpopulation.pdf`
- `hprcv2.1Mb.subtelo.mds.color-by-arm.pdf` / `.color-by-chromosome.pdf` / `.color-by-superpopulation.pdf`

### 2.4 No image yet for the two-domain results
The `/moosefs/guarracino/HPRCv2/PHR_III/plots/` directory contains the two-domain TSVs (`two_domain_test.tsv`, `two_domain_changepoint.tsv`, `two_domain_per_sequence.tsv`, `two_domain_binned_means.tsv`, `its_breakpoint_coloc.tsv`) but **no rendered figure**. See §6 for the suggested figure.

---

## 3. Existing CSVs / TSVs

All paths verified by directory listing.

### 3.1 Heterogeneity outputs (`/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/`)
- `allele_vs_paralog_distance.tsv` — per-community Wilcoxon results (§1.1)
- `cross_arm_affinity_sequences.tsv` — per-sequence cross-arm classification (§1.4, §1.11)
- `cross_arm_type_discordance.tsv` — per-individual concordant/discordant calls (§1.5)
- `cross_arm_region_length_comparison.tsv` — Wilcoxon length comparison (§1.6)
- `cross_arm_gene_content.tsv` — gene-content overlap matrix (§1.7)
- `gene_conversion_scores.tsv` — affinity-arm replacement scores (§1.7)
- `cross_arm_tar1_comparison.tsv` — TAR1 Fisher / Wilcoxon results (§1.8)
- `cross_arm_superpop_enrichment.tsv` — per-arm Fisher exact + BH (§1.4)
- `fst_superpop_matrix.tsv` — Hudson pairwise Fst (§1.4)
- `island_exchange_status.tsv` — proximal/distal × cross/self χ² inputs (§1.9)
- `within_community_arm_separation.tsv` — silhouette / separation ratio (§1.2)

### 3.2 Two-domain analysis (`/moosefs/guarracino/HPRCv2/PHR_III/plots/`)
- `two_domain_test.tsv` — per-arm Spearman (48 arms; §1.3 Test 1)
- `two_domain_per_sequence.tsv` — per-sequence ρ (99.7 % statistic; §1.3 Test 1)
- `two_domain_changepoint.tsv` — per-arm piecewise breakpoint + R² (41 arms; §1.3 Test 2)
- `two_domain_binned_means.tsv` — binned mean n_chrs vs distance (input for changepoint plots)
- `its_breakpoint_coloc.tsv` — closest internal-TTAGGG block to each arm's breakpoint (§1.3 Test 3)

### 3.3 Sample metadata
- `/moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv` — superpopulation labels for the 232 individuals (§1.4)

---

## 4. Methods (extracted)

| Sub-analysis | Statistical method | Multiple-testing | Notes |
|---|---|---|---|
| Allele vs paralog distance (§1.1) | Wilcoxon signed-rank, paired allele-vs-paralog Jaccard distance per individual per community | None reported across communities (each community tested independently) | Jaccard distance from k-mer set comparison; "paralog" = closest sequence from a different arm in the same community |
| Arm separation (§1.2) | Silhouette score on within-community Jaccard distances; "separation ratio" = inter-arm / intra-arm distance | n/a | Only multi-arm communities are testable (single-arm communities excluded) |
| Cross-arm affinity (§1.11) | Per-sequence rule: nearest neighbour in same community is "foreign-arm" | n/a (descriptive) | 2,484 / 15,668 sequences cross-arm |
| Population structure (§1.4) | 2 × 5 contingency Fisher exact per arm-community pair (cross/self × AFR/AMR/EAS/EUR/SAS) | BH-correction across 19 pairs | All 232 individuals have superpop annotation |
| Pairwise Fst (§1.4) | Hudson estimator on binary site (cross-arm = 1, self-arm = 0) per pair, averaged across the 10 strongest pairs | n/a | Mean Fst = 0.044 |
| Type discordance (§1.5) | Per-individual classification of the two haplotypes; rate = discordant / (concordant + discordant) | n/a | 35 arm-community pairs reported; threshold > 20 % flagged |
| Region length (§1.6) | Wilcoxon rank-sum on PHR length, cross vs self per arm-community | BH across 18 pairs | 14 / 18 significant |
| Gene replacement (§1.7) | Fraction of cross-arm genes matching affinity arm vs own arm; reported as "score" | n/a | Score ∈ [0, 1]; affinity-specific gene count reported alongside |
| TAR1 by class (§1.8) | Fisher exact for prevalence; Wilcoxon for length | BH across 19 pairs | 3 / 19 pairs significant for prevalence |
| Island position by class (§1.9) | χ² on proximal vs distal × cross vs self | n/a (single test) | χ² = 0.00, p = 0.99 |
| Within-community Jaccard distribution (§1.10) | Pairwise Jaccard sampling (≤ 4,950 pairs/community); visual bimodality assessment | n/a | Qualitative; cross-references Ambrosini 2007 |
| Two-domain — gradient (§1.3 Test 1) | Spearman ρ between per-window n_chrs and distance from telomere; per-arm + per-sequence | None across 48 arms (each tested independently) | 5 kb windows, 5 kb step, 500 kb max; 18,827 haplotypes |
| Two-domain — breakpoint (§1.3 Test 2) | Piecewise (two-segment) linear regression vs single linear, F-test for fit improvement | n/a | 41 testable arms (≥ 5 unique distance bins) |
| Two-domain — ITS colocalization (§1.3 Test 3) | Distance from arm's breakpoint to nearest internal (TTAGGG)n block > 5 kb from telomere; Spearman across 19 arms | n/a (descriptive correlation) | 467 per-sample BED RepeatMasker files; 29,274 annotations |

### 4.1 Scripts
Path: `/moosefs/guarracino/HPRCv2/scripts/`
- `community/analyze-within-arm-heterogeneity.R` — affinity, population bias, type discordance, gene replacement
- `community/plot-within-arm-heterogeneity.R` — plotting (the 12 panel suite)
- `community/allele_vs_paralog_distance.R` — §1.1
- `community/compute_fst_superpop.py` — Hudson Fst (§1.4)
- `similarity/analyze_polymorphic_arms.py` — polymorphic-arm characterization
- `similarity/test_d4z4_causality.py` — D4Z4-specific tests for C1
- `similarity/test_two_domain.py` — Spearman gradient (§1.3 Test 1)
- `similarity/test_two_domain_changepoint.py` — piecewise + colocalization (§1.3 Tests 2–3)
- `similarity/extract_rm_ttaggg_subtelomeric.py` — RepeatMasker (TTAGGG)n / TAR1 / (ACCCTA)n extraction
- `similarity/compute_its_breakpoint_coloc.py` — ITS-to-breakpoint distances (§1.3 Test 3)

---

## 5. Gaps

1. **No rendered figure for the two-domain analysis.** The TSVs exist (`two_domain_test.tsv`, `two_domain_changepoint.tsv`, `two_domain_per_sequence.tsv`, `two_domain_binned_means.tsv`, `its_breakpoint_coloc.tsv`) but no PDF/PNG. This is the most quantitatively powerful claim in the section and currently has zero visual support.
2. **Mechanism caveat is acknowledged but never resolved.** The section opens with "cannot distinguish among gene conversion, unequal crossover, reciprocal exchange, or other processes". The pedigree analysis (file `14_pedigree_recombination.md`) presumably addresses this; this section should at minimum cross-reference whatever direct mechanism evidence the paper carries.
3. **No outgroup for the Fst story.** The conclusion in §1.4 explicitly flags that population-specific exchange cannot be separated from drift / ILS without a non-human outgroup. A chimpanzee or gorilla comparison is not in scope here.
4. **C7 (acrocentric) "complete homogenization" is supported by silhouette = −0.029 and 70.5 % paralog-closer — but no per-individual time-since-exchange estimate.** The data establish that homogenization has occurred; they do not establish how recently.
5. **chr4_q discordance (43.4 %) exceeds the published Dutch ~20 %** but the comparison is to a different assay (Southern / PFGE for full D4Z4 array translocations vs broader PHR sequence affinity). A proper apples-to-apples re-analysis using the original D4Z4 boundary would clarify whether the discrepancy is methodological or biological.
6. **The cross-arm affinity classifier itself is not benchmarked.** A sequence is cross-arm if its nearest community neighbour is from a foreign arm. Robustness to neighbour-count k and to alternative distance metrics (other than Jaccard) is not reported.
7. **Single-arm communities (C8/C9/C10/C13) are excluded from silhouette and most cross-arm tests by design** — but this leaves four communities un-described in the heterogeneity story even though the genes they contain (e.g. chr15_q, chr16_p, chr17_p, chr4_p) are biologically interesting.
8. **No interaction analysis between the population-structure signal (§1.4) and the gene-content signal (§1.7).** Are AFR-enriched cross-arm haplotypes also gene-content-shifted, or only sequence-shifted?
9. **The two-domain breakpoint range (10–445 kb) is reported but not biologically grouped.** Are short-breakpoint arms enriched for any structural feature (TAR1 distance, segdup density, GC) vs long-breakpoint arms?
10. **Jaccard distance is the only metric used.** Reproducing the headline claim with an alignment-based or compositional metric would harden the conclusion (Ambrosini-style identity scores, or k = 21–31 minimizer Jaccard sweep).

---

## 6. Suggested figures (produced-vs-todo)

Each entry below is "title — caption — production status".

### Produced — ready for the paper as-is

**F1. Allele vs paralog distance, all multi-arm communities.**
Caption: Per-individual paired Wilcoxon distance (allelic vs paralog) across nine multi-arm communities; allele closer in 8 / 9 communities (overall p < 1e-300, 5,946 pairs). C7 (acrocentric p-arms) is the lone reversal (70.5 % paralog closer, p = 2.0e-7).
Status: produced. Likely backed by `within_arm_heterogeneity_03_mds_zoom.pdf` + table from `allele_vs_paralog_distance.tsv`. A clean two-panel figure (boxplot + reversal table) may need a small re-render from the TSV.

**F2. Silhouette / arm separation across communities.**
Caption: Per-community silhouette score and separation ratio. C2 / C12 (sil > 0.8) retain arm identity; C7, C14, C15 (sil ≤ 0) are interchangeable; C1 (sil = 0.147) is the autosomal weak spot (D4Z4).
Status: produced — `within_arm_heterogeneity_01_separation_overview.pdf`.

**F3. Cross-arm affinity heatmap.**
Caption: Per-arm cross-arm affinity rate; chrX_q (99.7 %), chrX_p (93.3 %), acrocentric p-arms (43–94 %), chr11_p (74.1 %) lead the autosomal ranking.
Status: produced — `within_arm_heterogeneity_02_affinity_heatmap.pdf` and `within_arm_heterogeneity_10_affinity_network.pdf`.

**F4. Superpopulation composition of cross-arm haplotypes.**
Caption: Cross-arm vs self-arm super-population composition for the 10 / 19 significant arm-community pairs (Fisher BH p_adj < 0.05); AFR strongly differentiated (mean Fst 0.10–0.15 vs non-AFR).
Status: produced — `within_arm_heterogeneity_04_superpop_composition.pdf` + `fst_superpop_matrix.tsv`.

**F5. Subtelomeric type discordance heterozygosity heatmap.**
Caption: Per-arm fraction of individuals heterozygous for haplotype class. chr22_q-C6 (47.5 %) and chr4_q-C1 (43.4 %) lead.
Status: produced — `within_arm_heterogeneity_05_heterozygosity_heatmap.pdf`.

**F6. PHR region-length distribution by exchange class.**
Caption: Cross- vs self-arm PHR length per arm-community. chrX_p cross 25× longer than self; chr22_q cross shorter (20 vs 25 kb); chr16_q cross longer (210 vs 155 kb). 14 / 18 pairs significant.
Status: produced — `within_arm_heterogeneity_06_region_length.pdf`.

**F7. Gene-repertoire conversion-score plot.**
Caption: Cross-arm gene-content match to affinity arm. chr13_p, PAR1, PAR2 saturate at 1.000; D4Z4 balanced ≈ 0.5; f7501 ≈ 0.7.
Status: produced — `within_arm_heterogeneity_07_gene_content.pdf` + `_09_conversion_score.pdf`.

**F8. TAR1 in cross- vs self-arm.**
Caption: TAR1 prevalence and length difference by class. Three pairs differ significantly; PAR1 essentially TAR1-free.
Status: produced — `within_arm_heterogeneity_11_tar1_prevalence.pdf` + `_12_tar1_cross_vs_self.pdf`.

### Todo — to be produced

**F9 (NEW). Two-domain model: per-arm gradient + breakpoint composite.**
Caption: (a) Per-window mean n_chrs vs distance from telomere for the six focus arms (chr4p, chr4q, chr16p, chr18p, chr20p, chr22q) with two-segment fits overlaid (R² gain noted); (b) per-arm Spearman ρ across all 48 arms (39 negative, 99.7 % of individual haplotypes negative); (c) breakpoint position vs nearest internal (TTAGGG)n block, highlighting the 11 / 19 arms within 25 kb (chr8p 0.1 kb, chr20p 0.6 kb, chr13q 1.5 kb, chr7p 3.5 kb).
Status: **TODO.** Inputs ready: `two_domain_binned_means.tsv`, `two_domain_changepoint.tsv`, `two_domain_per_sequence.tsv`, `its_breakpoint_coloc.tsv`. No rendered figure currently in `/moosefs/guarracino/HPRCv2/PHR_III/plots/`.

**F10 (NEW). Within-community Jaccard distance bimodality, per community.**
Caption: Pairwise Jaccard distance density plots for C1 / C2 / C3 / C5 / C6 / C7 / C11 / C12. Bimodal in C2 (peaks 0.02 / 0.37) and C12 (0.04 / 0.73); diffuse in C1 / C7 — direct visual signature of homogenization.
Status: **TODO.** Plotting script (`plot-within-arm-heterogeneity.R`) likely already supports this; not present in the panel suite.

**F11 (NEW). C7 acrocentric homogenization summary.**
Caption: Three-panel figure for the talk's punchline: (a) C7 allele-vs-paralog reversal (70.5 % paralog closer); (b) C7 silhouette = −0.029 vs the rest of the communities; (c) chr13_p conversion score = 1.000 toward all four other acrocentric p-arms.
Status: TODO (composite of existing data; no current single-figure version).

**F12 (NEW). Cross-arm vs self-arm length-bimodality on chrX_p (PAR1) and chr22_q (C6).**
Caption: Two side-by-side density plots showing the 25× length difference at chrX_p (cross 500 kb vs self 20 kb) and the 5 kb length gap at chr22_q (cross 20 kb vs self 25 kb) — the two extremes of the §1.6 distribution.
Status: TODO.

---

## 7. Talk slide takeaways (15-min talk)

Suggested allocation: 3 slides for this section.

### Slide A — "Arm identity is preserved (almost) everywhere"
- Headline: **5,946 paired observations, allele closer than paralog in 8 of 9 multi-arm communities, p < 1e-300.**
- Visual: F1 (allele-vs-paralog) + F2 (silhouette).
- Punchline: arm-of-origin signature survives despite cross-arm exchange. **One exception: acrocentric p-arms (C7) — paralog closer in 70.5 % of pairs, silhouette = −0.029. Acrocentric homogenization is complete at population scale.**
- Secondary: PAR1 / PAR2 (C14, C15) behave as expected (silhouette ≤ 0). D4Z4 (C1) is the autosomal weak spot (44.2 % paralog closer, silhouette 0.147).

### Slide B — "Two-domain subtelomere extends genome-wide"
- Headline: **Flint–Mefford two-domain model holds across 39 / 48 arms; 99.7 % of individual haplotypes show the gradient.**
- Visual: F9 (TODO) — per-arm gradient panel + breakpoint composite, internal (TTAGGG)n co-localization.
- Numbers: chr4p ρ = −0.85, chr22q ρ = −0.95; piecewise > linear on 39 / 41 arms; ITS within 25 kb of breakpoint on 11 / 19 testable arms (chr8p 0.1 kb, chr20p 0.6 kb).
- Punchline: the model originally inferred from a handful of arms (Flint 1997 / Mefford 2002) is the rule, not the exception, across all 48 chromosome ends — but the breakpoint position is arm-specific (10–445 kb).

### Slide C — "Cross-arm exchange leaves a population signature"
- Headline: **15.9 % of subtelomeric sequences are cross-arm; AFR vs non-AFR Fst = 0.10–0.15 vs ~0 within non-AFR.**
- Visual: F4 (superpop composition) + F5 (heterozygosity).
- Numbers: 10 / 19 arm-community pairs show significant superpop bias (Fisher BH p_adj < 0.05). chr22_q discordance 47.5 %, chr4_q 43.4 % — subtelomeric exchange is a segregating polymorphism, not a rare event.
- Mechanism caveat (one line): the analyses detect outcomes, not mechanisms — gene conversion vs unequal crossover vs reciprocal exchange are not distinguished without pedigree data (cross-reference 14_pedigree_recombination).
