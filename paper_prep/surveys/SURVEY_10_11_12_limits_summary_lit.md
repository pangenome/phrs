---
title: "Survey 10/11/12 — Limitations, summary, literature & novelty"
sources:
  - end-to-end-report/report/10_limitations.md
  - end-to-end-report/report/11_summary.md
  - end-to-end-report/report/12_literature.md
scope: Caveats, the 10–12 anchoring biological findings, and the 24+7 novel/predictive claim ledger
audience: Nature manuscript and 15-min talk
---

# Survey 10/11/12 — Limitations, key findings, literature

This survey extracts and structures the content of three short closing
sections of `end-to-end-report/report/`:

- `10_limitations.md` — sample composition + 18 numbered methodological /
  3D-validation caveats.
- `11_summary.md` — the 12 anchoring biological findings (numbered 1–12)
  with metrics that the paper and talk are built on.
- `12_literature.md` — confirmed prior-literature claims, 27 novel-to-this-paper
  contributions, and 7 explicit testable predictions (some already tested in
  the report).

These sections do not introduce new figures or new TSVs of their own — they
are synthesis prose that points back to artefacts produced upstream
(sections 01–09 and 14). Where a figure or CSV is named below it is the
underlying source for the claim, not a new artefact created in 10/11/12.

---

## 1. Key findings / claims with metrics

### 1.1 From `10_limitations.md` — caveats and lower-bound qualifiers

The limitations section makes 18 numbered points that constrain how every
finding in 11 and 12 should be read. They split into three groups.

**Sample composition (preface, no number):**

| Superpop | N | % |
|---|---|---|
| AFR | 67 | 28.8 |
| EAS | 52 | 22.3 |
| AMR | 44 | 18.9 |
| SAS | 37 | 15.9 |
| EUR | 33 | 14.2 |
| **Total** | **233** | **100.0** |

Within-superpopulation heterogeneity is not modelled. Population-scale
findings (the cross-arm-affinity / Fst section) inherit this imbalance.

**Methodological limitations (1–8):**

1. **Identity threshold = 95 %.** Lower bound — Ambrosini et al. (2007)
   reported a bimodal duplicon-identity distribution at 91 % and 98 %; the
   95 % cut captures the recent peak and misses much of the older peak
   (many olfactory and immunoglobulin duplicons). The 15.9 % cross-arm
   rate and 7 no-signal arms are floors, not censuses.
2. **Flank size = 500 kb.** May truncate longer subtelomeric similarity;
   no sensitivity analysis performed.
3. **Region-length floor = 3 kb output, 5 kb window/step.** Shorter
   shared segments are missed.
4. **Assembly quality.** Subtelomeres are among the hardest regions to
   assemble; the chr18_q chimera in `NA18982#1` (JBKABS010000018.1 fuses
   chr18 with 966 kb of chrX PAR1 across a 100 bp NNN scaffold join)
   illustrates how assembly errors can distort inter-chromosomal signal.
5. **Community-detection resolution.** The 50-community Leiden solution
   (k-NN = 75, resolution = 0.8) is one of a family; modularity does not
   guarantee biological correctness.
6. **Small N.** chrY_p has only 10 cross-arm sequences out of 92 in the
   population enrichment test (p_adj = 0.028); chr15_p discordance is
   based on 6/22 individuals (27.3 %). These should be treated as
   preliminary.
7. **Exchange timing.** Cross-arm affinity demonstrates that exchange
   has occurred but cannot date individual events. High discordance is
   consistent with recurrent exchange or a single ancient event still
   segregating; trio/family or population-genetics modelling is
   required to distinguish these.
8. **Somatic exchange in cell lines.** Most HPRCv2 assemblies derive
   from LCLs; some cross-arm affinity (especially chr4_q/chr10_q in C1)
   could reflect somatic exchange during culture rather than germline
   polymorphism (van Overveld et al. 2000 via Mefford & Trask 2002).

**3D-validation limitations (9–18):**

9. **Somatic vs meiotic context.** All Hi-C / Pore-C / Dip-C captures
    somatic interphase. The meiotic bouquet — when ectopic recombination
    actually occurs — has never been imaged by Hi-C in humans. The
    somatic 3D signal is interpreted as a residual of meiotic
    organisation (Rabl configuration), but the inference is indirect.
    Human meiotic Hi-C is the single most informative missing experiment.
10. **N = 6 Hi-C** (5 diploid HPRC + 1 haploid CHM13). Sufficient for
    reproducibility, insufficient for population-scale 3D claims. CHM13
    shows no significant communities, primarily because of reduced
    power (37 arms vs 75 in diploid samples).
11. **GM12878.** EBV-transformed B-LCL with abnormal karyotype; PBMC
    provides primary-cell control but with fewer cells.
12. **hg19 / T2T coordinate incompatibility** in Tan et al. 2018 Dip-C;
    impg projection at 50 kb resolution introduces noise near assembly
    gaps and subtelomeres.
13. **Multi-mapping at PHR intervals.** All four 3D technologies (Hi-C,
    Pore-C, CiFi, Dip-C/sperm) disable MAPQ filters; each multimapped
    read keeps exactly one randomly chosen alignment (no `-k`/`--all`),
    so noise is symmetric. The flanking-region control rules out
    multi-mapping as the source of the signal, but the absolute PHR
    enrichment values (B/W 0.027–0.074) cannot be cleanly separated
    from a multi-mapping contribution.
14. **Confound controls (positive note).** Five exclusion sets
    (acrocentric p-arms, sex chromosomes, both, all acrocentric p+q +
    sex, strongest C1+C7+C14/C15) all *increase* Mantel ρ when
    excluded — HG002 0.657 → 0.790, HG02148 0.152 → 0.720 —
    ruling out nucleolar association, PAR sharing, and D4Z4 as drivers.
    Rabl is addressed by the flanking control (signal *strengthens*
    away from the tip, not weakens).
15. **Parameter sensitivity (positive note).** Multi-resolution mcool
    analysis at 5/10/20/50/100 kb across human, RPE-1 and mouse shows
    the core signal is robust to resolution choice.
16. **Fragmented assemblies.** HG02148 and NA19036 produce NaN flanking
    values at some subtelomeres → excluded from flanking-only analyses
    but kept for PHR-based analyses.
17. **Dip-C cell 12 = duplicate of cell 10** (shared SRR7226706
    library); excluded → 16 GM12878 cells used.
18. **Mouse 1 Mb PHRs ≈ window size.** Median mouse PHR ≈ 980 kb at
    1 Mb extraction → window saturation. Larger windows (1.5–2 Mb) may
    expose deeper sharing, especially in C1 (16 arms).

### 1.2 From `11_summary.md` — twelve anchoring findings

The summary numbers 12 findings (the source numbers 1–12; the original
brief mentions "10" because items 5/11/12 are stated more tersely).
Each is one paragraph in the source; metrics are reproduced verbatim
where the paper-prep team needs them as numbers to quote.

| # | Finding | Anchor metric(s) |
|---|---|---|
| 1 | Subtelomeric regions form discrete inter-chromosomal communities. | 41 arms → 15 arm-level / 50 sequence-level communities; 233 individuals × 465 haplotypes × 15,668 PHR sequences. |
| 2 | Three categories of subtelomeric architecture. | (a) Homogeneous: 8/41 arms with 0 % cross-arm. (b) Polymorphic: 34/41 spanning ≥ 2 seq-level communities. (c) Fully interchangeable: 7/41 with 100 % cross-arm. |
| 3 | Recurrent inter-chromosomal exchange. | Arm-level cross-arm 15.9 % (2 484 / 15 668); seq-level 11.1 % (1 740 / 15 668). Discordance up to 47.5 %. |
| 4 | Extensive gene-repertoire overlap at acrocentric and sex chromosome subtelomeres. | chr13_p replacement score 1.000 toward all other acrocentric p-arms; chr14_p (N = 229) 83.0 % cross-arm with other acrocentric p-arms; PAR1/PAR2 score 1.000. |
| 5 | Population-specific exchange histories. | 10 arm/community pairs FDR-significant (Fisher's exact, BH); chr16_q cross-arm 70 % AFR (60/86); chrX_p self-arm 82 % AFR (18/22). |
| 6 | Subtelomeric gene content dominated by pseudogenes. | Protein-coding 4–9 % per community; PAR1 exception 32.1 %. RPL23AP45 spans 10 communities / 21 arms; SEPTIN14P22 spans 9 / 22; DUX4L → C1, MTCO → C7, SHOX → C15. |
| 7 | TAR1 as subtelomeric marker. | Present in 94.6 % of sequences, all 41 arms; near-absent from PAR1. |
| 8 | 3D genome organisation mirrors sequence communities. | Human Hi-C B/W 0.027–0.074 (8 datasets); RPE-1 HPRC B/W 0.005–0.052 (3 datasets, PHR coords); Dip-C T2T 6.9 % closer (16 cells, Fisher p = 2.4 × 10⁻⁵, Mantel ρ = 0.30); sperm W/B = 0.401 (60 % closer, 20 cells); mouse Hi-C per-arm-pair ρ = 0.57–0.72, Mantel ρ = 0.58–0.72 (4 stages, 1 Mb, per-haplotype PHR coords); 2 Mb 0.61–0.70; 4 Mb 0.65–0.73. Negative control: 7 non-sharing arms 11 % farther (GM12878), 40 % farther (sperm). Community-free arm-level Dip-C ρ = 0.34, human Hi-C ρ = 0.66–0.83. |
| 9 | Flanking unique-sequence regions show equal or stronger 3D signal. | HG002 100 kb flanking B/W = 0.002 vs PHR B/W = 0.027 — multi-mapping ruled out. |
| 10 | C4 (chr7_q / chr12_q) as minimal-PHR positive control. | Significant in 4/5 diploid Hi-C samples despite only 5–25 kb PHR and zero gene annotations. |
| 11 | Community-specific 3D predictions confirmed. | C7 ↔ nucleolar association; C1 peripheral (D4Z4-proximal lamin tethering, Masny et al. 2004); C14/C15 strongest in male samples; singletons C8/C10/C13 enriched in highest-depth samples (homolog–homolog contacts). |
| 12 | Two-domain subtelomeric model supported genome-wide. | Flint/Mefford gradient confirmed on 39/48 arms; breakpoint structure on 39/41 testable arms; 99.7 % of individual haplotype sequences. TTAGGG co-localisation within 25 kb on 11/19 arms. Breakpoint range 10–445 kb. |

These twelve are the load-bearing claims of the paper. The talk should
quote items 1, 3, 4, 8, 9, 12 verbatim (see § 7).

### 1.3 From `12_literature.md` — claim ledger

`12_literature.md` is a three-part ledger over the same biology:

- **Confirmed literature claims** — six headline confirmations plus six
  quantitative-only sub-confirmations.
- **Novel contributions** — 27 numbered claims. (The source heading
  reads "24 findings" but the numbered list runs 1–27; treat the
  numbered list as authoritative.)
- **Testable predictions** — 7 numbered hypotheses (3 already tested in
  the report, 1 not testable, 3 not yet tested).

§ 4 of this survey itemises all three lists. Headline metrics with
direct quotability for the manuscript:

- Sequence/3D feedback loop is the proposed unifying mechanism (novel
  contribution #11).
- Out-of-Africa topology recovered from cross-arm-affinity frequencies
  (novel contribution #15; testable prediction #6, *tested*) supports
  Mefford & Trask's (2002) hedged phylogenetic suggestion.
- Recombination map vs cross-arm affinity correlation **vanishes
  (ρ = 0.00, p = 0.98, N = 32) when 7 confounded acrocentric/PAR arms
  are excluded** (testable prediction #7, *tested* — null result;
  important caveat for the manuscript discussion).

---

## 2. Existing figures (paths)

`10/11/12` themselves do not introduce new figures. Every figure cited
by these sections is produced upstream. Below is the minimal cross-link
table back to the upstream surveys (paths reproduced from those
surveys; verify in the upstream survey if you need newer copies).

### 2.1 Figures cited from §10 caveats

| Caveat | Upstream artefact | Survey |
|---|---|---|
| 4. chr18_q chimera in NA18982#1 | Inter-chromosomal detection panel | SURVEY_01_pipeline.md |
| 5. 50-community Leiden solution | Community structure plots (`plot-seq-community-structure.R` outputs) | SURVEY_01_pipeline.md / 02_annotation.md |
| 6. Population enrichment chrY_p / chr15_p | Population structure heatmap | SURVEY_04_heterogeneity.md |
| 9. Somatic vs meiotic | Mouse meiotic Hi-C panels (Zuo et al.) | SURVEY_05_hic_validation.md / 08_mouse |
| 10. Hi-C N = 6 panel set | Per-sample W/B + Mantel scatter | SURVEY_05_hic_validation.md |
| 11. GM12878 / PBMC | Dip-C scatter + radial | SURVEY_06_dipc_validation.md |
| 13. Multi-mapping at PHRs | MAPQ-filter ablation panel | SURVEY_05_hic_validation.md |
| 14. Exclusion-set Mantel ρ | `output_q0_XX/community_enrichment_k50/exclusion_no_*` PDFs | SURVEY_05_hic_validation.md / 06 |
| 15. Multi-resolution mcool | Resolution-sensitivity multipanel | SURVEY_05_hic_validation.md |

### 2.2 Figures cited from §11 anchoring findings

Each numbered finding maps to one or more upstream PDFs. The single
most economical "key findings" panel for the paper does not yet exist
(see § 6, **T-1**).

| Finding | Existing artefact (representative) | Survey |
|---|---|---|
| 1. 15 / 50 communities | `plot-seq-community-structure.R` PDF | SURVEY_01 / 02 |
| 2. Three architecture categories | Polymorphic-arm narratives (`generate_polymorphic_arm_narratives.py`) | SURVEY_04 |
| 3. Cross-arm 15.9 % / 11.1 % | Pangenome-scale cross-arm bar/heatmap | SURVEY_01 |
| 4. Acrocentric and PAR overlap | Replacement-score heatmap | SURVEY_03_gene_enrichment.md |
| 5. Population-specific exchange | Fisher exact heatmap | SURVEY_04 |
| 6. Pseudogene dominance | Gene-class composition stacked bars | SURVEY_03 |
| 7. TAR1 distribution | TAR1 prevalence per community | SURVEY_02 |
| 8. 3D mirrors communities | `gm12878_mantel_scatter.pdf`, sperm Mantel, RPE-1 Mantel, mouse meiotic Hi-C | SURVEY_05 / 06 |
| 9. Flanking paradox | Flanking-vs-PHR W/B panel | SURVEY_05 |
| 10. C4 minimal-PHR control | C4 Hi-C zoom | SURVEY_05 |
| 11. Community-specific 3D predictions | Radial-position figure(s); cell-biology overlays | SURVEY_06 |
| 12. Two-domain model | `test_two_domain.py` / `test_two_domain_changepoint.py` outputs (gradient + changepoint plots) | SURVEY_01 (or dedicated) |

### 2.3 Figures cited from §12 literature

| Claim | Artefact | Survey |
|---|---|---|
| Mefford chr16_p bimodality | chr16_p length-distribution histogram (long-allele tail) | SURVEY_04 |
| Mefford chr3_q homogeneity | Within-arm Jaccard variance per arm | SURVEY_04 |
| IL9R distribution | Liftoff-projected gene heatmap | SURVEY_03 |
| D4Z4 / DUX4L causality | `test_d4z4_causality.py` outputs | SURVEY_05 (D4Z4 panel) |
| OR4F pseudogenisation gradient | Per-arm pseudogene fraction | SURVEY_03 |
| Out-of-Africa from cross-arm freqs | Population-tree dendrogram + Fst heatmap (testable prediction #6 result) | SURVEY_04 |
| Recombination ρ = −0.43 → 0.00 | Cross-arm % vs cM/Mb scatter, with confound annotation | SURVEY_04 / new (see § 6, T-3) |
| Two-domain model | gradient plot, changepoint plot per arm | SURVEY_01 / 14_pedigree_recombination |

---

## 3. Existing CSVs / TSVs (paths)

`10/11/12` cite no new TSVs. The numbers quoted in those sections are
recomputed from upstream tables. The minimal pointer set:

### 3.1 §10 (limitations)

- Sample composition counts (233 / 5 superpops): `hprc-sequence-production.tsv` (referenced in §12 caveat 5 but resolved in §10).
- Cross-arm 15.9 % / arm-level cross-arm rates: per-arm cross-arm TSV from `find-multichr-regions-incremental.py` outputs (SURVEY_01).
- Discordance 47.5 % cap: per-arm discordance TSV (SURVEY_04).
- Mantel ρ exclusion-set values (HG002 0.657 → 0.790; HG02148 0.152 → 0.720): `output_q0_XX/community_enrichment_k50/exclusion_no_*/<sample>_mantel_3d.tsv` and the matched flanking files (SURVEY_05 / 06).

### 3.2 §11 (summary)

- Same upstream TSVs as the corresponding row in § 2.2. The summary
  section does not introduce new tables.
- One TSV that should be assembled (does not yet exist as a single
  file): the **paper-key-findings table** — twelve rows, one per
  finding, with columns for finding label / anchor metric / source
  TSV / source figure. See § 6, **T-1**.

### 3.3 §12 (literature)

- `/moosefs/guarracino/HPRCv2/scripts/similarity/test_literature_claims.py`
  — runs the six testable literature claims listed under "Confirmed
  literature claims" (cited explicitly in §12, end of the section).
- Companion scripts in the same directory used for §12 results:
  - `test_two_domain.py`, `test_two_domain_changepoint.py` — Flint/Mefford two-domain confirmation (novel contribution #23, #24).
  - `test_d4z4_causality.py` — D4Z4-driver confirmation.
  - `compute_its_breakpoint_coloc.py` — internal (TTAGGG)n co-localisation with breakpoints (Ambrosini et al. 2007; novel contribution #17 and confirmation row).
  - `analyze_polymorphic_arms.py`, `generate_polymorphic_arm_narratives.py` — arm-classification engine for novel contributions #2, #4.
- Recombination-map test (testable prediction #7): T2T-CHM13 Lalli et al. 2025 recombination map intersected with arm-level cross-arm-affinity TSV (SURVEY_04). The intermediate per-arm table is reproduced inline in §12 lines 75–96.
- Cross-arm-affinity Fst / population tree (testable prediction #6 + novel contribution #19): per-superpopulation cross-arm count TSV from SURVEY_04.
- IL9R distribution counts (confirmation row): Liftoff annotation TSV from SURVEY_03.

---

## 4. §12 ledger — explicit lists

This section satisfies requirement (4) of the task: an explicit list of
(a) confirmed prior literature, (b) novel-to-this-paper contributions,
(c) testable predictions. All three lists are reproduced from
`12_literature.md` so that the manuscript and talk can quote without
re-reading the source.

### 4.1 Confirmed prior literature

Six headline confirmations:

C1. **Ambrosini duplicon architecture (Ambrosini et al. 2007).** All 11
    subtelomere-specific duplicon block entries (Table 1, blocks 1–3,
    5–8, 10–12, plus block 6′; blocks 4 and 9 do not exist in the
    original) map systematically to the 15 Leiden communities. Four
    confirmed at gene level: OR4F → C8, DUX4L → C1, IQSEC3 → C5,
    IL9RP1 → C3. Block 5 (chr2_p anchor) has no PHR signal, consistent
    with chr2_p exclusion. Three Ambrosini-noted genes (FBXO25, TUBB4q,
    RYD5) are not detected in current Liftoff annotations — likely
    nomenclature drift. The 91 %/98 % bimodal-identity peak structure
    explains why the 95 % threshold under-represents older exchanges
    (limitation 1). **Caveat reproduced in source:** Ambrosini's
    specific claim about TTAGGG island co-localisation with internal
    duplicon-to-duplicon boundaries is **not** tested — the present
    test evaluates the PHR outer boundary, a different feature.

C2. **Flint/Mefford two-domain model (Flint et al. 1997; Mefford &
    Trask 2002).** Inter-chromosomal sharing decreases with distance
    from the telomere on 39/48 arms; discrete two-phase breakpoint
    structure on 39/41 testable arms; internal (TTAGGG)n islands
    co-localise with the boundary within 25 kb on 11/19 arms with
    detectable ITS. All 6 originally FISH-characterised arms (chr4p,
    chr4q, chr16p, chr18p, chr20p, chr22q) confirm.

C3. **Mefford & Trask subtelomeric exchange.** (a) f7501 block African
    enrichment at chr7p/chr16_q — chr16_q cross-arm 70 % AFR (60/86).
    (b) IL9R pseudogene distribution chr9_q/chr10_p/chr16_p/chr18_p
    matches across communities C2/C3/C9 (and PAR copies in C14/C15).
    (c) chr4_q/chr10_q exchange detected at higher sensitivity
    (43.4 % discordance vs 20 % by Southern blot) — methodological
    sensitivity, not a contradiction.

C4. **Zuo meiotic alignment (Zuo et al. 2021).** Median 105 kb PHR fits
    within Zuo's ~500 kb leptotene loop and the ~20 % chromosome-end
    alignment zone, positioning PHR sequences where ectopic
    recombination can occur. **Caveat:** the inference from Zuo's
    measurements to PHR positioning is the *present* analysis's
    synthesis, and the extrapolation is cross-species (mouse → human)
    and cross-context (meiotic → somatic).

C5. **Tan Dip-C radial positions (Tan et al. 2018).** Per-chromosome
    radial preferences are consistent with cell-biology predictions
    when aggregated by community: C1 peripheral = D4Z4-proximal lamin
    tethering; C6 interior = nucleolar proximity. **Caveat:** the
    community-level aggregation is the present analysis, not Tan et
    al.'s.

C6. **Linardopoulou et al. (2005).** "Human subtelomeres are hot spots
    of interchromosomal recombination and segmental duplication" — the
    paper provides the population-scale quantification (15.9 %
    cross-arm, 11.1 % seq-level cross-arm, discordance up to 47.5 %).

Six quantitative-only sub-confirmations (also explicitly numbered in
the source, lines 16–21):

C7. **Mefford chr16_p bimodality** — 71.2 % short alleles (median 25 kb), 28.8 % long (median 205 kb); long-allele AFR enrichment 59.5 %, cluster in C28; matches Mefford's "~30 %" prediction.
C8. **Mefford chr3_q homogeneity** — within-arm Jaccard variance 0.022, 2.8× lower than chr15_q. CV 1.84 (chr3_q) > 1.39 (chr15_q) > 0.83 (chr19_p), arguing against a pure length artefact.
C9. **Mefford IL9R distribution** — IL9RP1 chr9_q (445 hap), IL9RP3 chr16_p (432), IL9RP4 chr18_p (444), IL9RP2 chr10_p (33); plus chr13_q (3), chrX_q (311), chrY_q (76), chrY_p (2), chr16_p (3). Cross-arm sequences carry IL9R at 100 % rate where testable.
C10. **chr19_p structural polymorphism (Linardopoulou)** — chr19_p in 7 seq-level communities, second only to chr6_q.
C11. **D4Z4 as causal driver of C1** — inter-chromosomal signal peaks 0–15 kb (D4Z4 position), C1 median 22 DUX4L vs all 7 outliers 0–2 on their own arm (Mann-Whitney p = 5.3 × 10⁻⁶); outlier PHRs 4.6–9× shorter.
C12. **OR4F pseudogenisation gradient** — 62.1 % of 5 023 OR4F annotations are pseudogenes; pseudogenisation rate 11.1 % (chr7_p) → 99.8 % (chr15_q) across 16 arms.

**Conclusion (verbatim from §12):** "No tested paper claims were
directly contradicted. The 43.4 % vs 20 % discrepancy at
chr4_q/chr10_q reflects methodological sensitivity, not conflicting
biology."

### 4.2 Novel contributions (27)

Numbering matches §12 lines 29–55 (the source heading reads "24
findings" but the enumerated list runs 1–27).

N1.  **Population-scale community structure** — first quantification of 41 arms → 15 communities across 233 individuals / 465 haplotypes.
N2.  **Three-category arm classification** — homogeneous / polymorphic / fully interchangeable; quantitative extension of Mefford & Trask's qualitative patchwork model.
N3.  **Population-scale cross-arm affinity** — 15.9 % (2 484 / 15 668), first across 465 haplotypes.
N4.  **Subtelomeric type discordance at population scale** — up to 47.5 % structural heterozygosity per locus.
N5.  **Gene-repertoire replacement scores** — complete (0.91–1.0) at chr13_p / chr14_p / chr15_p / PAR; partial (0.0–0.72) at autosomal communities — first homogenisation gradient mapped onto community structure.
N6.  **3D genome mirrors sequence communities** — three independent technologies; effect sizes vary across samples; finer 50-community partition does not reach significance in Dip-C.
N7.  **Flanking-region paradox** — unique-sequence flanks 100 kb centromere-ward show stronger 3D signal than duplicated PHR (HG002: B/W 0.002 vs 0.027). Strong evidence against multi-mapping artefact.
N8.  **C4 minimal-PHR positive control** — 5–25 kb at chr7_q / chr12_q tips → 3D co-localisation in 4/5 samples.
N9.  **Cell-type specificity** — PBMC W/B = 0.983, p = 0.305 (n.s.; mixed cells, hg19 projection, smaller PHRs); GM12878 T2T 6.9 % closer is the primary Dip-C finding.
N10. **Per-individual discordance ↔ 3D correlation** — low discordance → strong Hi-C signal (ρ = −0.50, p = 3.4 × 10⁻⁴, N = 48 sample × community pairs).
N11. **Proposed feedback loop model** — sequence similarity → 3D proximity → ectopic exchange → increased similarity. Causal direction not established.
N12. **MTCO pseudogene enrichment** at acrocentric p-arms (C7); enabled by T2T-quality assemblies.
N13. **TAR1 near-absence from PAR1** (0.5 % at chrX_p / chrY_p vs 94.8 % genome-wide).
N14. **Weak A-compartment with interior positioning** — 68 % A-compartment but mean e1 = +0.007; radial 0.50–0.63.
N15. **Subtelomeric exchange frequencies as phylogenetic markers** — out-of-Africa topology recovered across 9 arms. Supports Mefford & Trask's hedged 2002 suggestion.
N16. **TAR1 prevalence consistent with passenger status** — C7 cross-arm rates 86.6 % (TAR1+) vs 90.5 % (TAR1−); near-universal prevalence outside PAR (99.6 %) precludes a definitive functional test.
N17. **Internal (TTAGGG)n islands at population scale** — 18 352 islands / 8 321 sequences (53.1 %) / all 41 arms; median 79 bp. Distinct from TAR1 and from terminal arrays.
N18. **Alleles closer than paralogs, except at C7** — Wilcoxon p < 1 × 10⁻³⁰⁰ in 9/10 multi-arm communities; C7 reversed in 65.5 % of individuals (paralog distance < allelic distance).
N19. **Subtelomeric Fst mirrors out-of-Africa history** — AFR Fst 0.10–0.12 vs all non-AFR; non-AFR 0.00–0.02; mean Fst = 0.048 across 9 arms.
N20. **Internal (TTAGGG)n islands predominantly degenerate** — only 32.2 % "pure canonical"; TGAGGG 19.0 %, TTGGGG 16.0 %, TCAGGG 12.7 %. Quantifies Ambrosini variant motifs at population scale.
N21. **Cross-arm and self-arm sequences share identical (TTAGGG)n island distribution** — 45.9 % vs 45.9 % proximal; χ² = 0.00, p = 0.99, OR = 1.00.
N22. **Within-community Jaccard distance is multi-modal** — C2/C12 bimodal (allelic vs paralog peaks); C1/C7 diffuse (homogenised).
N23. **Two-domain subtelomeric model supported across most arms at pangenome scale** — gradient on 39/48 arms (81 %), breakpoint on 39/41 (95 %), 99.7 % of haplotype sequences. Internal (TTAGGG)n co-localises with breakpoint within 25 kb on 11/19 testable arms.
N24. **Two-domain boundary positions quantified** — 15–445 kb arm-specific; chr4p 70 kb, chr4q 50 kb, chr22q 15 kb; consistent with FISH-inferred 25–50 kb range but at sequence resolution.
N25. **Mouse meiotic enrichment across all stages** — per-arm-pair ρ peaks at zygotene (ρ = 0.715, p = 4.4 × 10⁻⁵⁵, n = 344); leptotene 0.680, pachytene 0.677, diplotene 0.574. Flanking-region control matches temporal pattern (ρ = 0.60–0.77).
N26. **Pangenome community-free framework validated cross-organism** — mouse 1 Mb per-arm-pair ρ 0.574–0.715; flanking 0.604–0.766. Human Dip-C 16 cells Mantel ρ = 0.296, p = 0.002.
N27. **Sperm signal in haploid post-meiotic cells** — W/B = 0.401 (60 % closer); Mantel ρ = 0.202, p = 0.023. Consistent with Rabl persistence through spermiogenesis.

### 4.3 Testable predictions

T1. **LINC complex and meiotic alignment** (highest priority, *not yet
    tested in this paper*). Apply community framework to Zuo et al.
    SUN1 W151R mutant Hi-C (GEO: GSE155142 mutant; GSE155638 +
    GSE155967 WT). **Caveat:** median PHR ≈ 105 kb occupies < 2 % of
    most chromosomes, well within the ~5 % tip zone where SUN1 mutant
    contacts *increase* — predicted effect on PHR-scale contacts is
    therefore uncertain.

T2. **Haplotype-resolved 3D contacts** (*tested in this paper, result
    inconclusive*). At discordant arms in 4 Hi-C samples (HG002 fully
    concordant), 7/8 informative pairs show the cross-arm haplotype with
    equal or higher partner contact — opposite of prediction.
    N too small (10 pairs total, NA19036 Wilcoxon p = 0.5).

T3. **TAR1 facilitator vs passenger** (*tested, conclusion: passenger*).
    Outside PAR, only 55 / 12 962 sequences lack TAR1 (99.6 % present).
    In C7 (only community with both TAR1+ and TAR1−), cross-arm rates
    are indistinguishable (86.6 % vs 90.5 %).

T4. **CTCF/cohesin at PHR boundaries** (*not yet tested*). Requires
    ENCODE CTCF ChIP-seq + bedtools intersection with PHR boundary
    coordinates.

T5. **Somatic vs germline exchange** (*not testable from public data*).
    `hprc-sequence-production.tsv` does not include a DNA-source
    column (LCL vs blood); requires external metadata linking sample
    IDs to source material.

T6. **Subtelomeric phylogenetic markers** (*tested, supports Mefford &
    Trask 2002 hypothesis*). Tree across 9 arm/community × 5 superpops
    recovers out-of-Africa topology: AFR deepest split (0.73–0.81 from
    all others), then EAS, then SAS/AMR/EUR clustering (AMR–EUR 0.22).

T7. **Crossover frequency correlation** (*tested, important null
    result*). Lalli et al. 2025 T2T-CHM13 recombination map vs cross-arm
    affinity across 39 shared arms: ρ = −0.43, p = 0.006. **But:**
    excluding 7 confounded acrocentric/PAR arms (0–12 callable variants
    per 500 kb due to short-read genotyping limitations in repetitive
    rDNA-adjacent regions), **ρ = 0.00, p = 0.98, N = 32 — the
    correlation vanishes entirely**. The signal across all 39 arms is
    driven exclusively by the extreme values at confounded arms; the
    biological question of whether local recombination protects arm
    identity remains open.

---

## 5. Gaps for the manuscript

The closing trio of sections is mostly synthesis; the gaps below are
what 10/11/12 *don't* yet do, and what the manuscript will need.

1. **No single "key findings" master panel.** §11 enumerates 12
   findings but the paper has no one-figure visual summary that places
   all 12 anchors next to their numbers. See § 6, **T-1**.
2. **Novel contributions table missing.** §12 lists 27 novel claims
   but as prose. Reviewers and the talk will both want a Table 1-style
   ledger (claim / metric / supporting figure / supporting TSV / prior
   literature it extends). See § 6, **T-2**.
3. **Recombination null result needs its own panel.** Testable
   prediction T7's all-39-arm ρ = −0.43 and N = 32 ρ = 0.00 are the
   manuscript's most quotable null result and the easiest to
   misinterpret in a single sentence. See § 6, **T-3**.
4. **Out-of-Africa tree from cross-arm affinity** is a striking
   independent recovery of population history but currently lives only
   inside §12 as a paragraph (testable prediction T6 + novel
   contribution N15/N19). It deserves its own Fst heatmap + dendrogram
   panel in the paper (and a single talk slide). See § 6, **T-4**.
5. **Quantitative confound disclosure for §10 limitations** is uneven:
   limitations 1 (95 % threshold), 6 (small N), 8 (LCL somatic
   exchange) and 13 (multi-mapping) cap interpretations of specific
   findings — but the paper does not yet provide a one-line "subject
   to limitation X" annotation per anchoring finding. See § 6, **T-5**.
6. **Two-domain model output is split across scripts.** §12
   confirmation C2 / N23 / N24 are produced by `test_two_domain.py`
   and `test_two_domain_changepoint.py` but no consolidated
   per-arm-table / panel exists. See § 6, **T-6**.
7. **Sperm + GM12878 + Hi-C consistency claim** is asserted in
   §11 finding 8 but no single panel shows W/B / Mantel ρ across all
   four 3D technologies side-by-side. SURVEY_06 § 6 already proposes
   a similar T-6; this should be promoted to a paper main-text figure
   (see § 6, **T-7**).
8. **Talk does not yet have a "what we did *not* claim" slide.** §10's
   18 caveats include several that anticipate reviewer questions
   (chimeric NA18982 contig, GM12878 abnormal karyotype, hg19/T2T
   incompatibility for Tan et al. 2018, LCL somatic exchange, lack of
   human meiotic Hi-C). A single slide enumerating these pre-empts
   the obvious objection round. See § 7, slide 6.
9. **Liftoff missing-gene caveat** (Ambrosini's FBXO25, TUBB4q, RYD5
   not detected in Liftoff annotations) is mentioned only in §12 line
   9. The manuscript should record this as a Methods caveat with the
   proposed remediation (re-annotate with current Ensembl).
10. **No glossary linking community labels (C1…C50) to canonical
    cytogenetic / locus shorthand** (e.g., C1 ≈ D4Z4 / 4q35 family,
    C7 ≈ acrocentric p-arm rDNA-adjacent, C14/C15 ≈ PAR1/PAR2). The
    talk especially needs this for non-specialist audiences. See
    § 6, **T-8**.

---

## 6. Suggested figures (with captions)

Two of the eight items below are explicitly required by the task
brief: **T-1** (master "key findings" summary panel) and **T-2**
(novel contributions table). The other six are derived from the gaps
in § 5.

**T-1. Master "key findings" summary panel (paper main text + talk).**
*Caption:* "Twelve anchoring biological findings of the paper.
Numbered 1–12, each labelled with its single most quotable metric:
(1) 41 arms → 15 communities, 233 individuals; (2) 8/41 homogeneous,
34/41 polymorphic, 7/41 fully interchangeable; (3) 15.9 % cross-arm
arm-level (11.1 % at sequence level); (4) chr13_p replacement score
1.000, chr14_p 83.0 % cross-arm; (5) 10 arm/community pairs FDR-
significant, chr16_q 70 % AFR; (6) 4–9 % protein-coding (32.1 % at
PAR1); (7) TAR1 in 94.6 % of sequences, all 41 arms, near-absent
from PAR1; (8) Hi-C / Pore-C / CiFi / Dip-C / sperm / mouse all
positive (Mantel ρ 0.20–0.83); (9) HG002 flanking B/W 0.002 vs PHR
B/W 0.027; (10) C4 chr7_q/chr12_q significant in 4/5 samples with
5–25 kb PHR; (11) C7 nucleolar, C1 peripheral, C14/C15 male-strongest;
(12) gradient on 39/48 arms, breakpoint on 39/41 testable arms,
99.7 % of haplotypes." Build by composing thumbnails from upstream
SURVEY_01..SURVEY_09 figures listed in § 2.2; new artefact = the
panel composition + one TSV (`paper_key_findings.tsv`, 12 rows).

**T-2. Novel contributions ledger table (paper Table 1 + talk back-up
slide).** *Caption:* "Twenty-seven novel claims, each annotated with
(a) anchor metric, (b) supporting figure, (c) supporting TSV, (d)
prior literature the claim extends. Includes a column distinguishing
quantitative extension (e.g. N3, N5, N17) from observation new to
this paper (e.g. N9, N11, N15)." Source: § 4.2 above. New artefact =
`paper_prep/tables/novel_contributions.tsv` (27 rows × 5 columns).

**T-3. Recombination map vs cross-arm affinity (T7 null result).**
*Caption:* "Mean recombination rate (cM/Mb, Lalli et al. 2025
T2T-CHM13 map) vs cross-arm affinity across 39 shared arms.
All-arms Spearman ρ = −0.43 (p = 0.006, N = 39); excluding 7
acrocentric p-arms + PAR1 (chr13_p, chr14_p, chr15_p, chr21_p,
chr22_p, chrX_p, chrX_q) where 0–12 callable variants per 500 kb
reflect short-read genotyping limitations rather than recombination
biology, ρ = 0.00, p = 0.98, N = 32. The all-arm correlation is
driven exclusively by extreme values at confounded loci." Source
TSV: §12 inline table lines 75–96. New artefact = scatter PDF.

**T-4. Out-of-Africa from subtelomeric exchange (T6 + N15 + N19).**
*Caption:* "Pairwise Fst (cross-arm/self-arm as alleles, 9 arms ×
5 superpopulations) heatmap and the resulting population-tree
dendrogram. AFR is the deepest split (Fst 0.10–0.12 from all
non-AFR superpopulations), followed by EAS, with SAS/AMR/EUR
clustering (AMR–EUR Fst 0.22, closest pair). Mean Fst = 0.048."
Source TSVs: per-superpopulation cross-arm count TSV from
SURVEY_04. New artefact = Fst-heatmap + dendrogram PDF.

**T-5. Limitations × findings cross-reference (paper supplementary
table).** *Caption:* "Each row is one of the 12 anchoring findings
(rows of T-1) annotated with the §10 limitation numbers that bound
its interpretation. Highlights the four cross-cutting caveats —
LCL somatic exchange (#8), multi-mapping (#13), 95 % identity
threshold (#1), small N at chrY_p / chr15_p (#6)." New artefact =
`paper_prep/tables/limitations_x_findings.tsv` (12 rows × 18
columns or a long format).

**T-6. Two-domain model panel (C2 + N23 + N24).** *Caption:* "Per-
arm gradient (Spearman ρ between sequence sharing and distance from
telomere) and breakpoint position (kb) for all 41 arms. 39/48 arms
significant by gradient; 39/41 testable arms significant by
breakpoint; 99.7 % of haplotype sequences support the model.
Breakpoint range 10–445 kb. Internal (TTAGGG)n islands co-localise
with breakpoints within 25 kb on 11/19 testable arms. The 6
originally FISH-characterised arms (chr4p, chr4q, chr16p, chr18p,
chr20p, chr22q) are highlighted." Source: outputs of
`test_two_domain.py` + `test_two_domain_changepoint.py` +
`compute_its_breakpoint_coloc.py`. New artefact = composite per-arm
panel.

**T-7. Cross-technology W/B / Mantel ρ summary (promote SURVEY_06
T-6 to main text).** *Caption:* "Within/between 3D distance ratio
(W/B) and Mantel ρ across four single-cell and bulk technologies —
GM12878 Dip-C (W/B 0.93, ρ 0.30, n = 16), sperm scHi-C (W/B 0.40,
ρ 0.20, n = 20), bulk Hi-C (B/W 0.027–0.074, n = 8 datasets), mouse
meiotic Hi-C (Mantel ρ 0.58–0.72, 4 stages). Negative-control 7
non-sharing arms overlaid (11 % farther GM12878, 40 % farther
sperm)." Source TSVs already exist in SURVEY_06 § 3. New artefact =
single composite figure for paper main text.

**T-8. Community-to-locus glossary (talk slide; paper supplementary
figure).** *Caption:* "Each of the 15 arm-level communities
labelled with its dominant cytogenetic / locus identity: C1 D4Z4 /
4q35 family (chr4q, chr10q, plus 14 additional arms in mouse-like
extended sharing), C3 f7501 / IL9R-related, C5 IQSEC3, C6 C-tail
(chr22q etc.), C7 acrocentric p-arms (rDNA-adjacent, MTCO
pseudogenes, PAR-absent), C8 OR4F (olfactory family), C14 PAR2
(chrX_q / chrY_q), C15 PAR1 (chrX_p / chrY_p, SHOX). Also: 50-
community sequence-level partition collapses to 15 at arm level."
New artefact = annotated network/community plot.

---

## 7. Talk slide takeaways (15-min talk)

The talk should treat 10/11/12 as the **closing 3 slides** of the
biology section — the talk's "what we found / what's new / what we
didn't claim" beats. Suggested order:

1. **One-slide headline (slide ~9 of ~15).** "Across 233 individuals
   and 465 haplotypes we resolve 41 chromosome arms into 15 inter-
   chromosomal subtelomeric communities, with 15.9 % of sequences
   resembling a foreign arm more than their own — and 3D nuclear
   organisation mirrors this structure across four technologies, two
   cell types, and two organisms." (Use § 4 of T-1 + finding #8.)

2. **The three architecture categories slide.** Homogeneous (8/41) /
   polymorphic (34/41) / fully interchangeable (7/41) — quantitative
   extension of the qualitative Mefford & Trask patchwork model.
   (Finding #2; novel contribution #2.)

3. **Cross-technology 3D summary slide.** One panel: W/B and Mantel ρ
   across GM12878 Dip-C / sperm scHi-C / bulk Hi-C / mouse meiotic
   Hi-C. Negative control on the same panel. (Finding #8 / #9 / #10;
   T-7.)

4. **Out-of-Africa-from-subtelomeres slide.** Fst heatmap +
   dendrogram. "Subtelomeric exchange frequencies recover human
   population history — a 24-year-old hedged hypothesis (Mefford &
   Trask 2002), now testable from one pangenome." (Novel contribution
   #15 / #19; testable prediction #6 — *tested*; T-4.)

5. **What we *didn't* find slide (recombination null).** "Cross-arm
   affinity correlates with cM/Mb at ρ = −0.43 across 39 arms — but
   the correlation collapses to ρ = 0.00 once 7 confounded
   acrocentric/PAR arms are removed. Local recombination ≠ arm-
   identity protector at this resolution." (Testable prediction #7;
   T-3.)

6. **What we *can't* claim slide (top 5 caveats).**
   (a) No human meiotic Hi-C exists; somatic 3D is interpreted as a
       Rabl residual, *indirectly* (limitation 9).
   (b) 95 % identity threshold misses the 91 % older-exchange peak
       (limitation 1).
   (c) GM12878 has an abnormal karyotype; PBMC primary-cell control
       is hg19-projected and underpowered (limitations 11–12).
   (d) Most HPRCv2 assemblies are LCL-derived; chr4_q/chr10_q (C1)
       cross-arm signal could include somatic exchange in culture
       (limitation 8).
   (e) The recombination ρ = 0.00 null is real and important, not a
       failure to support — the data simply do not test the
       hypothesis at confounded loci (testable prediction #7).

7. **The "what comes next" slide (3 highest-priority asks).**
   (i) Human meiotic Hi-C (LINC mutants if possible).
   (ii) Re-annotate with current Ensembl/GENCODE to recover the
        Ambrosini-noted FBXO25 / TUBB4q / RYD5 absences.
   (iii) Source-stratified validation (LCL vs blood) of cross-arm
        affinity to bound the LCL somatic-exchange contribution.

Recall numbers, prioritised for memorisation in the talk:

| Number | Anchor |
|---|---|
| 233 / 465 | individuals / haplotypes |
| 41 → 15 | arms → arm-level communities |
| 15.9 % | cross-arm at arm level (11.1 % seq-level) |
| 47.5 % | maximum discordance |
| 6.9 % closer (Dip-C) / 60 % closer (sperm) | W/B-derived 3D proximity |
| ρ = 0.00, N = 32 | the recombination null after confound removal |
| 39/48 / 39/41 / 99.7 % | the three two-domain model coverages |
