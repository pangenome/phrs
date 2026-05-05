---
title: "Survey 07 — Integrated 3D interpretation (Discussion)"
source: end-to-end-report/report/07_integrated.md
scope: Discussion synthesis — convergent evidence, mechanistic models, testable predictions
audience: Nature manuscript Discussion + 15-min talk closing slides
---

# Survey 07 — Integrated interpretation

This survey extracts and structures the content of `end-to-end-report/report/07_integrated.md` for the Nature manuscript Discussion and the closing slides of the 15-minute talk. Unlike sections 03–06 (which carry primary results), section 07 is a **synthesis layer**: it does not produce its own figures or TSVs but aggregates effect sizes from all upstream analyses into a single convergent-evidence table and proposes four mechanistic models (flanking paradox, meiotic bouquet, D4Z4-CTCF-lamin, nucleolar association) plus a causal feedback loop and three testable predictions.

---

## 1. Key findings (with metrics)

### 1.1 Convergent evidence across technologies — the headline table (07.L11–L26)

A single table aggregates 14 independent tests across four assay types, six cell types, three platforms (Illumina/PacBio/ONT), and two species. Every test rejects the null in the same direction.

| Tech | Sample(s) | Test | Effect | p | Source |
|---|---|---|---|---|---|
| Hi-C | 5 HPRC + CHM13 | B/W ratio | 0.027–0.074 | 6e-18 to 9.4e-3 | survey 05 |
| Pore-C | HG002 | B/W ratio | 0.056 | 3.9e-85 | survey 05 |
| CiFi | HG002 | B/W ratio | 0.036 | 2.0e-74 | survey 05 |
| Dip-C (T2T) | GM12878, 16 cells | W/B 3D distance | 6.9 % closer | 2.4e-5 (Fisher) | survey 06 §1.1 |
| Dip-C (T2T) | GM12878, 16 cells | S_all (negative) | 11 % farther | (negative ctrl) | survey 06 §1.3 |
| Dip-C (T2T) | GM12878, 16 cells | Mantel ρ (Jaccard×3D) | ρ = 0.296 | 0.002 | survey 06 §1.1 |
| Dip-C (T2T) | GM12878, 16 cells | Community-free arm Spearman | ρ = 0.336 | 1.1e-18 | survey 06 §1.1 |
| Sperm 3D | 20 cells (Xu 2025) | W/B 3D distance | 60 % closer | 3.9e-51 (Fisher) | survey 06 §1.2 |
| Sperm 3D | 20 cells | Per-cell community-free | 15/20 ρ > 0 | — | survey 06 §1.2 |
| RPE-1 CiFi | Async (interphase) | B/W (50 / 10 kb) | 0.024 / 0.033 | 9.1e-102 / 4.2e-79 | (RPE-1 self section) |
| RPE-1 Pore-C | Async | B/W (50 / 10 kb) | 0.031 / 0.048 | 1.3e-95 / 1.7e-95 | (RPE-1 self section) |
| RPE-1 CiFi | Mitotic | B/W (50 / 10 kb) | 0.008 / 0.030 | 3.3e-68 / 3.0e-56 | (RPE-1 self section) |
| Mouse 1 Mb | 4 meiotic stages | B/W (5 resolutions) | 0.029–0.122 | all < 1e-22 | (mouse section) |
| Mouse 4 Mb | 4 meiotic stages | B/W (5 resolutions) | 0.019–0.057 | all < 1e-36 | (mouse section) |

Per-PHR-pair continuous correlation (sequence similarity ↔ contact strength): **CHM13 Hi-C ρ = 0.674, p = 3.0e-92**; **HG002 Pore-C ρ = 0.485, p = 1.6e-48**; all 8 datasets significant (07.L28). Multi-resolution: results stable across **5 kb / 10 kb / 20 kb / 50 kb / 100 kb mcool** in all systems, with slight strengthening at finer resolutions consistent with better bin coverage of the median 105 kb PHR (07.L30).

### 1.2 Flanking-region paradox (07.L34–L42)

PHR intervals contain duplicated sequence; the 100 kb centromere-ward of the PHR boundary is unique sequence. Multi-mapping discards (`RM_MULTI=1` in HiC-Pro) suppress signal inside PHRs.

- **Flanking B/W = 0.002 (HG002) – 0.057 (CHM13);** all stronger than PHR B/W (0.027–0.074). HG002 flanking is **13× stronger** than HG002 PHR (0.002 vs 0.027).
- **Dip-C radial:** GM12878 flanking particles are more interior than non-flanking terminal particles (0.503 vs 0.551, p = 7.4e-35).
- **Three implications.** (1) Multi-mapping is *not* the driver — flanking has no duplicated content and the signal is *stronger*. (2) 3D clustering extends past the duplicated region into unique flanking sequence — a broader chromosomal domain effect. (3) PHR effect sizes underestimate the true biological signal; flanking values are closer to ground truth.

### 1.3 Meiotic bouquet model (07.L44–L54)

All primary 3D data is somatic / interphase. The meiotic bouquet provides a candidate mechanism for *initial* sequence homogenization (Mefford & Trask 2002; Linardopoulou 2005; Patel 2019; Zuo 2021).

- **Mouse meiotic Hi-C (Zuo 2021):** chromosome-end alignment extends ~20 % of chromosome length during leptotene/zygotene; alignment is independent of compartment identity and depends on the LINC complex.
- **Loop sizes scale through prophase:** ~500 kb (leptotene) → ~700 kb (zygotene) → ~1.4 Mb (pachytene) → ~1.6 Mb (diplotene). The median PHR (105 kb) fits inside a single leptotene loop, placing it at the loop base where recombination machinery is concentrated.
- **Tan et al. 2018:** Rabl configuration "weak" in interphase GM12878 / PBMCs — somatic 3D underestimates meiotic proximity.
- **Bouquet is rare cytologically (<5 % of zygotene cells)** but end-alignment Hi-C signal is robust → the alignment, not the bouquet snapshot, is the operative phenomenon.
- **Critical gap:** human meiotic Hi-C does not exist; mouse spermatocyte data cannot be extrapolated directly because subtelomeric organization differs.

### 1.4 D4Z4-CTCF-lamin tethering model for C1 (chr4_q ↔ chr10_q) (07.L56–L72)

Specific molecular mechanism for the strongest single community.

- **C1 metrics:** silhouette = 0.147 (poor separation); 43.4 % chr4_q discordance; **Dip-C radial = 0.732 (peripheral)**; inter-chromosomal sequence-sharing signal **peaks at 0–15 kb from telomere** (D4Z4 location); C1 sequences carry **median 22 DUX4L genes** vs **0–2** for all 7 non-C1 outliers (Mann-Whitney p = 5.3e-6); outlier PHR regions are 4.6–9× shorter than C1.
- **Mechanism:** (i) D4Z4 macrosatellite at chr4_q / chr10_q tips; (ii) CTCF binds within D4Z4 repeats (Ottaviani 2009); (iii) D4Z4-proximal sequence is tethered to the nuclear periphery via lamin A/C (Masny 2004; Ottaviani 2009); (iv) co-peripheralization is consistent with elevated FSHD-relevant ectopic recombination (Lemmers 2010).
- **Testable consequence:** CTCF/cohesin density at PHR boundaries should correlate with Hi-C contact strength between community partners. Standard hg38 ENCODE is inadequate (incomplete subtelomeres). Two existing T2T-aligned datasets cover the gap: **Gershman 2022 ENCODE CTCF ChIP-seq realigned to T2T-CHM13** (CTCF enrichment at TAR loci across all ENCODE cell lines) and **Stergachis lab Fiber-seq** (single-molecule CTCF at 39 / 46 telomeres).

### 1.5 Nucleolar association mechanism for C6 / C7 (07.L74–L81)

- **C7 (acrocentric p-arms, chr13/14/15/21/22_p):** silhouette = −0.029 (interchangeable); gene-replacement scores 0.91–1.0 for chr13/14/15_p, 0.49–0.54 for chr21/22_p. **C7 is unmappable in Dip-C (hg19 p-arms).**
- **C6 (acrocentric q-arms + chr1_q, chr17_q, chr19_q):** silhouette = 0.521; **Dip-C radial = 0.505 (interior)** — consistent with nucleolar co-localization.
- **All five acrocentric short arms carry rDNA → constitutive nucleolar association** → near-complete interchangeability is consistent with frequent exchange at NOR-bearing arms.
- The **f7501 duplicon distribution (Mefford & Trask 2002)** maps onto C3 (chr3_q, chr7_p, chr9_q, chr11_p, chr16_q, chr19_p), confirming that community structure recovers known duplicon modules.
- **Note:** C7-specific 3D data is absent in this work; the model is consistent with literature only.

### 1.6 Causal feedback loop (07.L83–L95)

Sequence similarity → 3D proximity → ectopic exchange → increased similarity → (back to start).

| Link | Support level | Evidence |
|---|---|---|
| Similarity → 3D proximity | **Direct, this work** | Mantel & community-free correlations: HG002 ρ = 0.66, Pore-C ρ = 0.49, Dip-C ρ = 0.336 |
| 3D proximity → ectopic exchange | **Established literature** | FSHD chr4_q ↔ chr10_q D4Z4 translocations (Lemmers 2010); recombination requires physical proximity |
| Ectopic exchange → increased similarity | **Inferred from outcome** | 15.9 % of sequences cross-arm (survey 04); homogenization in C1, C7, PAR1, PAR2 |
| Increased similarity → stronger future proximity | **Not directly measured** | Would require temporal data or before/after exchange comparison within an individual |

- Circularity → causal direction is undecidable from the present data.
- Bouquet initiation (1.3) provides a plausible but untested origin scenario.
- Concept is not new: Linardopoulou 2005 ("segmental polymorphism / gross genomic rearrangement" cycle, their Fig. 2); Mefford & Trask 2002 ("interphase clustering … could promote their exchange"); Ambrosini 2007 (98 % identity peak in bimodal duplicon distribution explained by "ongoing interchromosomal gene conversion"). **This work's contribution is the quantitative 3D proximity dimension.**

### 1.7 Three testable predictions (07.L97–L105)

1. **LINC complex requirement (untested).** Zuo 2021 SUN1 W151R mutation collapses the alignment range from ~20 % to ~5 % of chromosome length. Prediction: meiotic Hi-C of SUN1-mutant spermatocytes analysed with the community framework should show dramatically reduced within-community inter-chromosomal contacts at chromosome ends. Would establish LINC-mediated force transmission as the causal driver.
2. **Crossover-frequency correlation (tested, confounded).** Subtelomeric recombination rate (Lalli 2025 T2T-CHM13 map) anti-correlates with cross-arm affinity across 39 arms (ρ = −0.43, p = 0.006), but the signal vanishes (ρ = 0.00, p = 0.98, N = 32) when 7 arms with 0–12 callable variants in 500 kb (acrocentric p-arms + PAR) are excluded. **Long-read recombination maps are required to test this properly.**
3. **Compartment identity at chromosome tips.** **68 % of chromosome tips classified A-compartment** (63 / 92 arm × haplotype), but **mean eigenvector e1 = +0.007** — weak, poorly defined identity. Subtelomeres are pseudogene-rich, A-leaning by GC, but **internally positioned** (C10 radial 0.474, C6 0.505) rather than at the lamina. Zuo 2021 explicitly found end-alignment is **independent of compartment identity** and that "loop size differences are unlikely the main cause for prominent alignment among chromosome ends". **Telomere clustering — not lamina association or loop size — drives end alignment.**

---

## 2. Existing figures referenced (paths)

Section 07 itself produces **no native figures**; it is a synthesis chapter. Every claim in the convergent-evidence table refers to a figure or TSV produced upstream. Cross-references:

### 2.1 Convergent-evidence sources
- **Bulk Hi-C / Pore-C / CiFi B/W panels** — see survey 05 (`/moosefs/guarracino/HPRCv2/PHR_III/HiC/`).
- **Dip-C / sperm Mantel scatters & radial panels** — survey 06 (`/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/` for GM12878; `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/` for sperm).
- **RPE-1 cell-cycle B/W panels** — RPE-1 self-validation section (`/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/`).
- **Mouse 1 Mb / 4 Mb meiotic-stage panels** — mouse section (path TBD; not surveyed in this round).
- **Multi-resolution stability** — referenced as "the resolution sensitivity section (human), the acrocentric exclusion control (no-acrocentric), the RPE-1 validation section.1, the mouse flanking Hi-C section". All exist as B/W panels at 5 / 10 / 20 / 50 / 100 kb.

### 2.2 Flanking paradox sources
- Hi-C flanking analysis: **`/moosefs/guarracino/HPRCv2/PHR_III/HiC/flanking_*`** (per-sample flanking B/W TSVs and panels — see survey 05).
- Dip-C flanking radial: GM12878 flanking-vs-non-flanking interior comparison lives in the per-particle radial outputs under `community_enrichment_k50/` (radial 0.503 vs 0.551 with p = 7.4e-35).

### 2.3 D4Z4 / C1 mechanism sources
- **D4Z4 causality test:** `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/` — `cross_arm_gene_content.tsv` carries the DUX4L count per sequence; the Mann-Whitney 5.3e-6 statistic comes from `scripts/similarity/test_d4z4_causality.py`.
- **0–15 kb sharing peak at D4Z4:** survey 02 (annotation section) — sequence-level community detection, intercommunity contact-by-distance plot.
- **C1 silhouette / discordance / radial:** survey 04 §1.2, §1.5; survey 06 §1.4.

### 2.4 Nucleolar / acrocentric sources
- **C7 silhouette and conversion scores:** survey 04 §1.2 / §1.7 (`within_arm_heterogeneity_01_separation_overview.pdf`, `_07_gene_content.pdf`, `_09_conversion_score.pdf`).
- **C6 radial:** survey 06 §1.4 radial table.
- **f7501 → C3 mapping:** referenced from the community detection section (annotation survey 02).

### 2.5 Causal-loop / model schematics
- **No produced schematic figure** for the four-link feedback loop. See §6 below.
- **No produced schematic** for the D4Z4-CTCF-lamin tethering model. See §6 below.
- **No produced schematic** for the meiotic bouquet end-alignment model. See §6 below.

### 2.6 Testable-predictions sources
- **Recombination map test:** Lalli 2025 preprint T2T-CHM13 recombination map (external); per-arm correlations table referenced as "the testable predictions section prediction 7" (located in section 07 of the report; full per-arm table not in this section).
- **Compartment / eigenvector data:** HG002 Hi-C eigenvector at 100 kb resolution, per-haplotype, GC-oriented (path: `/moosefs/guarracino/HPRCv2/PHR_III/HiC/...` — see survey 05).
- **CTCF data for D4Z4 prediction:** Gershman 2022 ENCODE CTCF ChIP-seq realigned to T2T-CHM13 (external); Stergachis lab Fiber-seq single-molecule CTCF at 39 / 46 telomeres (external).

---

## 3. Existing CSVs / TSVs (paths)

Section 07 has **no native CSVs**. The convergent-evidence table is hand-assembled from upstream outputs:

### 3.1 Bulk Hi-C / Pore-C / CiFi B/W
- `/moosefs/guarracino/HPRCv2/PHR_III/HiC/*_summary.tsv` — per-sample B/W ratios, p-values, multi-resolution sweep (see survey 05).
- `/moosefs/guarracino/HPRCv2/PHR_III/HiC/per_pair_correlation_*.tsv` — Spearman ρ between sequence similarity and Hi-C contact (CHM13 ρ = 0.674; HG002 Pore-C ρ = 0.485).

### 3.2 Dip-C and sperm
- `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_summary.tsv`, `gm12878_mantel_3d.tsv`, `gm12878_community_free_arm.tsv`.
- `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_16cells_500kb_per_community_per_cell.tsv` — S_all rows for the negative control.
- `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_summary.tsv`, `sperm_all20_mantel_3d.tsv`, `sperm_all20_community_free_per_cell.tsv`.
- (Full inventory: survey 06 §3.)

### 3.3 RPE-1
- `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/` — per-cell-cycle B/W TSVs at 10 kb / 50 kb (async, mitotic; CiFi + Pore-C).

### 3.4 Mouse meiotic
- Path TBD (mouse section directory not surveyed in this round); 1 Mb and 4 Mb B/W TSVs at 5 resolutions across 4 meiotic stages.

### 3.5 D4Z4 / DUX4L
- `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_gene_content.tsv` — DUX4L count per sequence (input for the C1 vs non-C1 outlier Mann-Whitney p = 5.3e-6).

### 3.6 Compartment / eigenvector
- HG002 100 kb eigenvector outputs under `/moosefs/guarracino/HPRCv2/PHR_III/HiC/` (per-haplotype, GC-oriented).

### 3.7 Recombination map
- Lalli 2025 T2T-CHM13 recombination map — external resource; per-arm subtelomeric rates referenced from section 07 of the report (preceding "testable predictions section prediction 7" forward reference).

---

## 4. Methods / models used

The integrated section uses no new statistical methods. Its "methods" are:

### 4.1 Aggregation methodology
- **Convergent-evidence table** (07.L11–L26): hand-curated assembly of B/W ratios and p-values from sections 03 (Hi-C / Pore-C / CiFi), 06 (Dip-C, sperm), the RPE-1 self-validation section, and the mouse section. Effect-size convention: **B/W ratio < 1** for Hi-C contact-frequency tests (within-community contacts above expectation) versus **W/B ratio < 1** for 3D-distance tests (within-community arms closer). The two conventions are inverse but both encode the same biological direction (community-coherent).
- **Multi-resolution check**: every B/W test repeated at 5 / 10 / 20 / 50 / 100 kb; consistency reported.
- **Per-pair Spearman correlation** (continuous version of the W/B test): sequence similarity (Jaccard) vs Hi-C / Pore-C contact strength, per inter-chromosomal arm pair; aggregated Spearman ρ.

### 4.2 Mechanistic models proposed
| Model | Inputs | Logical structure |
|---|---|---|
| **Flanking paradox** (1.2) | PHR vs flanking B/W; multi-mapping suppression | Deductive: multi-mapping discards explain *why* flanking > PHR; therefore signal extends beyond duplicated region |
| **Meiotic bouquet** (1.3) | Zuo 2021 mouse meiotic Hi-C; Patel 2019; Mefford 2002; Linardopoulou 2005 | Speculative but literature-supported: meiotic end-alignment → ectopic exchange in early prophase loops |
| **D4Z4-CTCF-lamin tethering** (1.4) | Ottaviani 2009; Masny 2004; Lemmers 2010; this work's C1 metrics | Mechanism for one specific community; falsifiable via CTCF-density-vs-contact correlation |
| **Nucleolar association** (1.5) | rDNA on five acrocentrics; this work's C6/C7 metrics | Existing literature aligned with C7 silhouette and C6 radial position |
| **Causal feedback loop** (1.6) | Mantel & community-free results; cross-arm affinity rate; FSHD literature | Four-link cycle with three direct/established links and one inferred-only link |

### 4.3 Compartment / eigenvector method
- HG002 Hi-C eigenvector at **100 kb resolution**, **per-haplotype**, **GC-oriented** (sign convention fixed by GC content). Tip = arm-terminal bin. 68 % of tips A-compartment (63 / 92), mean e1 = +0.007 — weak signature.

### 4.4 Cross-talk between somatic and meiotic interpretations
The interpretive bridge throughout the section: somatic 3D data (interphase Hi-C, Dip-C of GM12878, sperm post-meiotic 3D) is taken as **necessary but not sufficient** evidence for the *meiotic* exchange model. Bouquet, LINC, and end-alignment data are mouse-derived and cannot be extrapolated quantitatively to human subtelomeres. The integrated section makes this caveat explicit (07.L52–L54).

---

## 5. Open gaps

1. **No human meiotic Hi-C.** Identified as "the single most informative missing experiment" (07.L54). All inferences about bouquet, LINC, end-alignment depth, and meiotic loop sizes rely on mouse extrapolation.
2. **No SUN1-mutant test of the LINC requirement.** Prediction 1 (1.7) is fully laid out but cannot be tested with existing data — requires SUN1 W151R mouse spermatocyte Hi-C re-analysed within the community framework.
3. **Recombination-map prediction is confounded.** The chr-end recombination rate vs cross-arm affinity correlation is driven by short-read mappability artefacts at acrocentric p-arms and PAR (0–12 callable variants in 500 kb). After excluding 7 arms ρ collapses to 0.00 (p = 0.98). Long-read-based recombination maps resolving repetitive subtelomeres are required.
4. **CTCF-density-vs-contact prediction not yet executed.** Datasets exist (Gershman 2022 T2T ENCODE CTCF; Stergachis Fiber-seq for 39 / 46 telomeres) but no quantitative correlation between CTCF/cohesin density at PHR boundaries and contact strength has been computed in this work.
5. **No before/after exchange comparison.** The fourth link of the causal feedback loop ("increased similarity → stronger future proximity") is inferred only — it cannot be tested with present data. Would require temporal or pedigree data.
6. **Bouquet-initiation scenario is offered, not evidenced.** The paper explicitly flags this in 07.L95: "offered as a plausible scenario, not an evidence-based conclusion".
7. **C7 has no 3D data.** Dip-C uses hg19 mapping for the sperm/PBMC reference and acrocentric p-arms are unmapped. C7 mechanistic discussion is therefore literature-only.
8. **No native Discussion-section figure.** The convergent-evidence table is the only visual; there is no synthesis figure (forest plot of effect sizes, model schematics, or feedback-loop diagram). All four figures in §6 below are TODO.
9. **Causal direction is undecidable from present data** (07.L93–L95). The integrated section explicitly acknowledges that communities could form because similar sequences are brought into proximity, OR sequences could become similar because they are in proximity.
10. **Subtelomeric compartment identity is weak (mean e1 = +0.007).** Eigenvector classifications at chromosome tips are noisy; A/B compartment binarisation may not be the right framework for this region. The integrated section concludes that **telomere clustering, not compartment identity, drives end alignment** — but does not propose a stronger replacement metric.
11. **Inter-section terminology drift** (B/W vs W/B). Bulk Hi-C uses B/W (between/within contact frequency); Dip-C uses W/B (within/between 3D distance). Both with values < 1 indicate community proximity, but this inversion is a documentation/onboarding hazard for downstream readers.
12. **The four-link causal loop is presented as text** with no schematic; the talk-ready visual is missing.

---

## 6. Suggested figures (produced-vs-todo)

**No produced figures** are native to section 07. All entries below are TODO. The first two are the highest-leverage candidates for the manuscript Discussion and the talk's closing slide.

### TODO

**T-1. Convergent-evidence forest plot (the headline Discussion figure).**
*Caption:* "Effect sizes across 14 independent tests of community-coherent 3D organisation. Each row is one assay × cell type × resolution combination; x-axis is the effect-size statistic (B/W or W/B; oriented so left of unity = within-community closer/stronger than between-community); error bars are bootstrap 95 % CIs; colour codes assay family (Hi-C / Pore-C / CiFi / Dip-C / sperm scHi-C / mouse meiotic). All 14 effect sizes are on the same side of the null. Annotations call out the headline numbers: HG002 Pore-C 3.9e-85, sperm 3.9e-51, Dip-C Mantel 0.296."
*Status:* TODO. Inputs: the convergent-evidence table (07.L11–L26) plus the per-test summary TSVs in survey 05 §3, survey 06 §3, RPE-1 self section, mouse section.

**T-2. Causal feedback loop schematic.**
*Caption:* "Four-link cycle: sequence similarity → 3D proximity → ectopic exchange → increased similarity. Arrows are coloured by support level: solid (direct measurement, this work — link 1), dashed-solid (established literature — link 2), dashed (inferred from outcome — link 3), dotted (not directly measured — link 4). Initiation node: meiotic bouquet (greyed; speculative, mouse-derived). Feedback arrow returns to similarity."
*Status:* TODO. Pure schematic; no underlying TSV. Belongs as Discussion figure 1A or as the closing-slide model.

**T-3. D4Z4-CTCF-lamin tethering model schematic.**
*Caption:* "Mechanism for C1 (chr4_q ↔ chr10_q): D4Z4 macrosatellite at chromosome tip → CTCF binds within D4Z4 (Ottaviani 2009) → lamin A/C tethering at nuclear periphery (Masny 2004) → co-peripheral positioning of chr4_q and chr10_q (Dip-C radial 0.732). Inset: 0–15 kb peak in inter-chromosomal sharing signal aligns with D4Z4 location."
*Status:* TODO. Schematic + one inset panel from the sequence-level community-detection signal (annotation survey 02 inputs).

**T-4. Meiotic bouquet / end-alignment model schematic.**
*Caption:* "Cartoon of leptotene/zygotene bouquet: telomeres clustered at the nuclear envelope; chromosome ends aligned over ~20 % of length (Zuo 2021). Inset shows the median 105 kb PHR fitting within a single ~500 kb leptotene loop — placed at the loop base where recombination machinery is concentrated. SUN1-mutant prediction overlay: alignment range collapses from ~20 % to ~5 % (Zuo 2021)."
*Status:* TODO. Schematic only; the Zuo 2021 reference figure can be cited but should be redrawn for permissions.

**T-5. Flanking paradox panel.**
*Caption:* "PHR vs flanking 100 kb B/W ratios per sample. HG002 flanking is 13× stronger than HG002 PHR (0.002 vs 0.027). Companion panel: Dip-C flanking radial 0.503 vs non-flanking terminal 0.551 (p = 7.4e-35) — flanking is *more interior*."
*Status:* TODO. Inputs: survey 05 flanking TSVs + survey 06 radial TSVs.

**T-6. Compartment-identity-at-tips diagnostic panel.**
*Caption:* "(a) Distribution of HG002 Hi-C eigenvector e1 at chromosome tips (n = 92 arm × haplotype); 68 % A-compartment (e1 > 0); mean e1 = +0.007. (b) Tip e1 vs Dip-C radial position; tips are A-leaning by GC but interior-positioned, dissociating compartment identity from nuclear-envelope tethering. (c) Per-community radial vs per-community mean e1 — telomere clustering, not lamina association, drives end alignment."
*Status:* TODO. Composite from HG002 Hi-C eigenvector + survey 06 radial TSVs.

**T-7. Per-arm cross-arm affinity vs subtelomeric recombination rate (the confounding panel).**
*Caption:* "Lalli 2025 T2T-CHM13 subtelomeric recombination rate vs cross-arm affinity per arm. All-39-arm fit ρ = −0.43, p = 0.006 (driven by 7 arms with 0–12 callable variants). After excluding acrocentric p-arms + PAR (N = 32), ρ = 0.00, p = 0.98. Highlights the short-read mappability confound."
*Status:* TODO. Inputs: per-arm recombination rate (Lalli 2025) + per-arm cross-arm affinity (`cross_arm_affinity_sequences.tsv`). Should be the explicit honesty figure for the Discussion.

---

## 7. Talk slide takeaways (closing slides, 15-min talk)

Recommended allocation: **3 closing slides** (the model + the convergent evidence + the testable predictions).

### Slide A — "Six independent technologies, one signal"
- **Headline:** *Bulk Hi-C, Pore-C, CiFi, Dip-C, sperm scHi-C, mouse meiotic Hi-C — all positive, all p < 1e-5, often p < 1e-50.*
- **Visual:** T-1 (convergent-evidence forest plot). One row per test; all 14 on the same side of unity.
- **Headline numbers (memorise these):**
  - HG002 Pore-C B/W = 0.056, p = **3.9e-85**
  - Sperm 3D W/B = 0.401 (60 % closer), p = **3.9e-51**
  - GM12878 Dip-C Mantel ρ = **0.296**, p = 0.002
  - CHM13 Hi-C per-pair ρ = **0.674**, p = 3.0e-92
- **Negative control to slide in:** S_all (non-sharing arms) reverses to 11 % farther in GM12878, 40 % farther in sperm — sequence sharing is necessary, not incidental.
- **Multi-resolution stability:** identical conclusion at 5 / 10 / 20 / 50 / 100 kb mcool.

### Slide B — "Why? A causal feedback loop with one bouquet-shaped origin"
- **Headline:** *Sequence similarity, 3D proximity, and ectopic exchange form a self-reinforcing cycle.*
- **Visual:** T-2 (feedback-loop schematic) overlaid with T-4 (bouquet inset for the initiation node).
- **Talk track:**
  - Link 1 (similarity → proximity) — measured here.
  - Link 2 (proximity → exchange) — established (FSHD).
  - Link 3 (exchange → similarity) — inferred from 15.9 % cross-arm rate.
  - Link 4 (similarity → future proximity) — not measured, requires temporal data.
- **Honest caveat (one line):** causal direction undecidable from present data; bouquet-initiation is plausible, not proven.
- **Specific instances on the next layer:** D4Z4-CTCF-lamin for C1; nucleolar tethering for C7; flanking-paradox-as-domain-effect for the broader chromosomal context.

### Slide C — "Three things this work cannot test (yet)"
- **Headline:** *Three falsifiable predictions, three different missing datasets.*
- **Visual:** three icons in a row (LINC complex / chromosome / CTCF molecule).
- **Predictions:**
  1. **LINC requirement** — repeat Zuo 2021 SUN1 W151R Hi-C, run through the community framework. Expectation: ~5×–10× drop in within-community inter-chromosomal contacts at chromosome ends.
  2. **CTCF density at boundaries** — correlate Gershman 2022 T2T-aligned ENCODE CTCF + Stergachis Fiber-seq with per-PHR Hi-C contact strength. Expectation: positive correlation; standard hg38 ENCODE inadequate.
  3. **Crossover-protected arms** — repeat the Lalli 2025 correlation against a *long-read* recombination map; current short-read maps cannot resolve subtelomeric variants on acrocentrics or PAR.
- **Closing line:** "The most informative missing experiment is human meiotic Hi-C — until then, mouse is our oracle and the signal is consistent across every cell type and platform we have tried."

---

## Cross-references
- Survey 01 (`SURVEY_01_pipeline.md`) — pipeline that produces the per-arm B/W and Mantel inputs.
- Survey 02 (`SURVEY_02_annotation.md`) — community detection, sequence-level inter-community signal at 0–15 kb (D4Z4 peak).
- Survey 04 (`SURVEY_04_heterogeneity.md`) — silhouette, cross-arm affinity, type discordance, gene replacement (the C1, C6, C7 metrics).
- Survey 05 (`SURVEY_05_hic_validation.md`) — bulk Hi-C / Pore-C / CiFi B/W and continuous correlations; flanking analysis.
- Survey 06 (`SURVEY_06_dipc_validation.md`) — Dip-C / sperm 3D distance, Mantel, S_all negative control, radial positions.
- Section `12_literature.md` — situates these findings within prior literature (forward reference at 07.L107–L109).
