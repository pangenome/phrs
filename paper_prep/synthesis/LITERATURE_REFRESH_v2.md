---
title: Round-2 literature refresh — executive synthesis
date: 2026-05-17
agent: r2-integrator-merge (agent-136)
inputs:
  - paper_prep/lit_review/REFRESH_R2_foundational.md (topics 01, 02, 03; agent-125)
  - paper_prep/lit_review/REFRESH_R2_recombination.md (topics 05, 06, 07, 11; agent-124)
  - paper_prep/lit_review/REFRESH_R2_3d-bouquet.md (topics 08, 09; agent-122)
  - paper_prep/lit_review/REFRESH_R2_methods.md (topics 10, 12; agent-121)
  - paper_prep/lit_review/REFRESH_R2_frontier.md (topics 13, 14, 15; agent-123)
references_base: paper_prep/synthesis/REFERENCES_v4.bib (364 entries)
references_output: paper_prep/synthesis/REFERENCES_v5.bib (372 entries; +8 new, deduped by DOI)
upgrade_plan: paper_prep/synthesis/CITATION_UPGRADE_PLAN_v2.md (round-2 ADDs only)
---

# Round-2 literature refresh — executive synthesis

## 0. Executive summary

The five R2 refresh sweeps cover every topic 01–15 partitioned into five non-overlapping bundles. Combined query volume: **~90 PubMed queries, 1,500+ bioRxiv preprints screened by title/abstract preview, OpenAlex sweeps for 4 topics.** Coverage window: 2024-06 to 2026-05 with explicit 60-day preprint frontier (2026-03-17 to 2026-05-17).

**Headline:** 8 genuinely new bibkeys identified (2 STRONG, 6 MEDIUM). One proposed key (`rhie2026rob` from R2_methods) is a duplicate of the existing `acrocentric_Cechova2026` entry (identical DOI `10.64898/2026.03.08.710242`) and is therefore not added to REFERENCES_v5.bib. The two STRONG papers are Porsborg 2025 (long-read primate testis recombination; calibrates the 133 gene-conversion-like patch interpretation in P9) and Ataei 2025 (LINE1 mechanism for acrocentric distal-junction nucleolar anchoring; first molecular mechanism for the C7 co-localisation model).

**No new contradictions** were introduced or resolved by R2 across any topic. The 60-day preprint frontier is essentially silent for our subtelomere niche.

The R2 sweep also produced explicit "no new strong findings" verdicts for 7 of 14 sub-bundles (R2_foundational topics 01, 03; R2_recombination topics 06, 07; R2_3d-bouquet topic 09 primary data; R2_methods topics 10, 12; R2_frontier topics 13, 14). Each null is documented with the queries that establish it (Section 5 of each refresh).

---

## 1. Per-topic synopsis (5 blocks)

### Block A — REFRESH_R2_foundational (topics 01, 02, 03; agent-125)

- **Topic 01 (Cytogenetic Foundations).** No new STRONG papers. R2_AUDIT_PLAN marked TRUST; verified by deeper sweep (yang2025chr2fusion / poszewiecka2023phasedancer already cover the 2024-2026 window). Section covers 1916–2005 historical material where no 2024+ substitution is possible.
- **Topic 02 (Subtelomere Structure — TAR1/ITS/TERRA).** DEEPEN per R2_AUDIT_PLAN P-1 angle. One STRONG find: `santagostino2025terra` (RNA 32(1):97–112; PMID 41193243) — T2T-CHM13v2.0-based map of TERRA promoters at **39 of 46** human subtelomeres, identifies **106 intrachromosomal TERRA-like promoters** at ITS loci, and **205 ITS loci transcribed in ≥1 cell line**. Upgrades the "roughly half" of Gershman 2022 / "more than half" of Rodrigues 2024 to a precise 39/46 count, and provides the first systematic T2T-based ITS transcription census. Author group is Nergadze/Giulotto — same lab as the foundational `subtelstruct_Nergadze2007` entries.
- **Topic 03 (PHR Concept).** No new STRONG papers. R2_AUDIT_PLAN marked TRUST; confirmed.
- **Two MEDIUM dropped** (Hsieh 2025 TERRA-QUANT/aging; Schmidt 2024 Telo-seq) — complementary to existing R1 STRONG cites but no new biological claim for the draft.

### Block B — REFRESH_R2_recombination (topics 05, 06, 07, 11; agent-124)

- **Topic 05 (Acrocentric / rDNA / Robertsonian).** Two NEW: `ataei2025line1dj` STRONG and `hao2024snul` MEDIUM.
  - `ataei2025line1dj` (Genes & Dev 39(3-4):280–298; PMID 39797762): primate-specific full-length LINE1 conserved at all five acrocentric distal junctions; CRISPR removal disrupts nucleolar positioning, transcriptional output, self-renewal. First molecular mechanism assigned to DJ nucleolar anchoring — directly upgrades the C7 co-localisation explanation in P3/EDFig.4.
  - `hao2024snul` (eLife 13; PMID 38240312): SNUL ncRNAs form monoallelically-expressed sub-nucleolar territories on individual NOR-containing chromosomes — addressable per-chromosome sub-domains *within* the shared nucleolus.
- **Topic 06 (DUX4/D4Z4/FSHD).** No new STRONG. R1 was comprehensive (7 STRONG papers).
- **Topic 07 (Concerted Evolution / NAHR).** **Critical null result:** no follow-up to Schweiger 2024 (P-3 audit angle). 8 targeted PubMed queries + OpenAlex returned zero new long-tract NCO mechanism papers in 2024-2026. Schweiger 2024 remains the sole primary cite for the two-process NCO model.
- **Topic 11 (Pedigree-Based Recombination).** One NEW MEDIUM: `sasani2026kfam` (bioRxiv 2026-03-06; PMID 41959501) — uses the **same CEPH1463/K1463 four-generation pedigree** as the PHR study; HiFi long-read sequencing of 8 million TR loci across 20 children; 1,270 expansions/contractions; 43 hyper-mutable loci. Porubsky and Eichler are co-authors. Contextualises PHR exchange within broader pedigree mutagenesis.

### Block C — REFRESH_R2_3d-bouquet (topics 08, 09; agent-122)

- **Topic 08 (Meiotic Bouquet / Envelope).** No new STRONG. The 2023-2025 mechanistic anchors (KASH5, SpdyA/SUN1, SUN1 NOA variant, CTCF loops, centromere-trigger) all already in R1. One MEDIUM: `bouquet_ChenCEP164Cilia2025` (bioRxiv 2025-12-04; PMID 41415467) — CEP164 KO eliminates zygotene cilia + sperm flagella; **meiotic chromosome pairing and DSB repair proceed normally**. Functional test of the bouquet–zygotene-cilium connection. Does not change a draft claim (the draft does not assert cilium necessity), so MEDIUM not STRONG.
- **Topic 09 (Hi-C 3D Methods).** **Key gap remains unfilled:** no human stage-resolved meiotic Hi-C paper has appeared in 2024-2026. Cheng 2024 NIDDK Micro-C preprint not yet published. Marín-Gual 2025 (already R1) and He 2023 (already R1) remain the most recent stage-resolved 3D meiotic data.
- 600 bioRxiv preprints (6 categories × 100) reviewed in the 60-day window: zero relevant.

### Block D — REFRESH_R2_methods (topics 10, 12; agent-121)

- **Topic 10 (Pangenome Graphs / IMPG).** No new STRONG. Two persistent R1 gaps confirmed in R2: (1) IMPG still has no peer-reviewed methods paper — the `pangenome_graphs_impg_*` keys remain `@misc`; (2) the 12% wfmash sampling rate / Erdős-Rényi argument is uncited in any external 2024-2026 publication. Helmy 2025 mouse pangenome (PMID 41330379, Garrison+Guarracino as co-authors) is WEAK — applies existing methods, no PHR-relevant claim. Bao & Weigel 2025 plant pangenome review — WEAK.
- **Topic 12 (HPRC Population Pangenomes).** No new STRONG. HPRC v2 main paper still not deposited as preprint or indexed by PubMed as of 2026-05-17 (14 months after data release) — `hprc_hprcv2_2025` must remain `@misc`. Schloissnig 2025 / Logsdon 2025 / Jeong 2025 already in v4.
- One MEDIUM dropped: Rhie et al. 2026 (bioRxiv 2026-03-10) "Biobank-scale genotyping of Robertsonian translocations." This is a **duplicate of the existing `acrocentric_Cechova2026` entry** (DOI `10.64898/2026.03.08.710242`). NOT added.

### Block E — REFRESH_R2_frontier (topics 13, 14, 15; agent-123)

- **Topic 13 (Subtelomere Pop-gen / F_ST).** No new STRONG. 12+ PubMed queries + 3 OpenAlex queries; absence of papers on FST at subtelomeric haplotypes confirms the manuscript's novelty.
- **Topic 14 (Olfactory / OR4F).** No new STRONG. Graham 2025 (altitude OR reduction, PMID 40562037) is WEAK. The REDO flag for REFRESH_14 (`@chi2025primate` hedge for OR4F paragraph) is a **citation-placement** issue addressed by `R2_AUDIT_PLAN.md` §5 edit 3, not a literature gap.
- **Topic 15 (Emerging).** Three NEW: one STRONG (`pedigree_Porsborg2025primaterecom`) and two MEDIUM (`subtelstruct_Lee2026SEPTIN14`, `pedigree_Joseph2024PRDM9indep`).
  - `pedigree_Porsborg2025primaterecom` (Nat Commun 16:10337; PMID 41285744): long-read testis sequencing from 6 primate species + 3 human sperm samples; 2,881 crossovers, 2,314 simple gene conversions, 555 complex events; **NCO mean tracts 22–95 bp, CO-associated mean tracts 318–688 bp**. Calibrates the 133 gene-conversion-like patch interpretation in v2 P9.
  - `subtelstruct_Lee2026SEPTIN14` (Mobile DNA 17(1); PMID 41699652): SEPTIN14P pseudogene dispersal via subtelomeric SD blocks across great apes; mechanistic basis for the SEPTIN14P22 hub-span named in P11.
  - `pedigree_Joseph2024PRDM9indep` (PNAS 121:e2401973121; PMID 38809707): gBGC-substitution estimator quantifies PRDM9-independent hotspot activity in 52 placental mammals; humans show a deficit at these sites — contextualises the rare non-PRDM9 NCO class in human sperm.

---

## 2. Cross-cutting themes

1. **C7 (acrocentric p-arms) is the most-strengthened community.** Two of the eight new papers (`ataei2025line1dj`, `hao2024snul`) provide molecular and sub-nucleolar mechanism for the rDNA-linked co-localisation pattern that drives community C7's negative silhouette and paralog-closer-than-allele signal. Together with `acrocentric_rdna_robertsonian_floutsakou2013` / `mcstay2016` / `acrocentric_Altemose2022`, the C7 mechanism story is now complete from sequence identity (>98% DJ) → physical anchor (LINE1 + SNUL ncRNA) → recombination output.

2. **P9 (gene-conversion-like patches) now has a tract-length calibrator.** `pedigree_Porsborg2025primaterecom` provides primate-scale tract-length distributions: NCO tracts 22–95 bp, CO-associated tracts 318–688 bp. This calibrates the mechanistic interpretation: 133 gene-conversion-like patches are likely a mixture of NCO and CO-associated events depending on tract length, not all NCO. Combined with `pedigree_Schweiger2024spermNCO` (already in v4) and `pedigree_Joseph2024PRDM9indep` (Topic 15 contextualiser), P9's recombination model is now triangulated by three orthogonal cites: tract-length distribution, sperm long-read NCO mechanism, and the PRDM9-independent hotspot landscape.

3. **TERRA / ITS axis is the most-strengthened structural story.** `santagostino2025terra` provides the precise T2T count (39 of 46) for TERRA-promoter-bearing subtelomeres and the first systematic ITS transcription census (205 transcribed ITS loci). Combined with the existing `rodrigues2024terra`, `gershman2022telomeres` and the Nergadze 2007 trilogy, the TAR1 prevalence sentence in v2 L40 can now be supplemented with quantitative TERRA-context anchors.

4. **The 60-day preprint frontier (2026-03-17 to 2026-05-17) is essentially silent for our niche.** Across all five sweeps, 1,500+ bioRxiv preprints were screened; the only directly relevant 60-day preprint is `sasani2026kfam`. This is informative — the manuscript occupies a genuinely sparse research niche, and the v2 citation backbone is not at risk of being outdated within typical reviewer-turn-around windows.

5. **Two persistent methods-citation gaps remain unfilled.** (a) IMPG has no peer-reviewed methods paper as of 2026-05; the four `pangenome_graphs_impg_*` keys remain `@misc`. (b) The 12% wfmash sampling rate / Erdős-Rényi connectivity argument is uncited externally. R2 found no candidate to close either gap. Both are listed in v2 Limitations (clause iv).

6. **No new contradictions.** Every R2 sweep explicitly confirmed §3 of its parent file: zero new 2025-11 to 2026-05 contradictions across topics 01–15. The Salsi 2026 D4Z4 hedge (Topic 06), Chi 2025 OR4F hedge (Topic 14), and Lalli 2025 cM/Mb collapse (Topic 11) — all R1-flagged contradictions — remain the only ones in play; all are either landed in v2 already or scheduled for r2-fix-apply.

---

## 3. Top-10 must-add citations (round-2 only)

Each entry below is **uniquely R2** (absent from REFERENCES_v4.bib and from any R1 REFRESH file) and lists its source REFRESH_R2 file. Numbers 1–8 are genuinely new; the list closes at 8 because there are no additional defensible MUST-ADD candidates from the R2 sweep — see "Why only 8" below.

| # | Bibkey | Score | Topic | Target paragraph | Source REFRESH_R2 file |
|---|--------|-------|-------|------------------|------------------------|
| 1 | `pedigree_Porsborg2025primaterecom` | **STRONG** | 15 / P-3 | v2 P9 (L46) — supplement `@pedigree_Schweiger2024spermNCO` with tract-length calibration | `REFRESH_R2_frontier.md` §1 STRONG-1 |
| 2 | `ataei2025line1dj` | **STRONG** | 05 | v2 P3 / EDFig.4 — DJ anchoring mechanism sentence alongside `@acrocentric_rdna_robertsonian_floutsakou2013` | `REFRESH_R2_recombination.md` §1 |
| 3 | `santagostino2025terra` | STRONG | 02 | v2 L40 (TAR1 prevalence) and Extended Data Fig. 3 caption — supplement `@Gershman2022` / `@rodrigues2024terra` | `REFRESH_R2_foundational.md` §1.1 |
| 4 | `bouquet_ChenCEP164Cilia2025` | MEDIUM | 08 | EDFig.8 caption (causal-loop) or open-questions Methods note — cilium dispensability for pairing/DSB | `REFRESH_R2_3d-bouquet.md` §1 |
| 5 | `hao2024snul` | MEDIUM | 05 | Methods M14 (3D / Dip-C) or EDFig.4 caption — sub-nucleolar territory mechanism | `REFRESH_R2_recombination.md` §1.2 |
| 6 | `sasani2026kfam` | MEDIUM | 11 | v2 P9 (Limitations / pedigree context) — K1463 TR mutagenesis on same pedigree | `REFRESH_R2_recombination.md` §2 |
| 7 | `subtelstruct_Lee2026SEPTIN14` | MEDIUM | 15 | v2 P11 (gene-enrichment, SEPTIN14P22 hub-span sentence) | `REFRESH_R2_frontier.md` §1 MEDIUM-1 |
| 8 | `pedigree_Joseph2024PRDM9indep` | MEDIUM | 15 / P-3 | v2 P9 supporting cite OR Methods (recombination model) | `REFRESH_R2_frontier.md` §1 MEDIUM-2 |

**Why only 8.** The R2 sweep produced 8 genuinely new bibkeys. The single additional proposal (`rhie2026rob` from REFRESH_R2_methods) was found to be a duplicate of the already-present `acrocentric_Cechova2026` entry (same DOI). All other 2024-2026 candidate papers surfaced by the R2 queries were already in REFERENCES_v4.bib, already in an R1 REFRESH file, or judged WEAK / off-topic (~22 papers dropped with rationale across the five files). The top-10 cap is reserved for genuinely round-2 MUST-ADD findings and is not padded with WEAK/MEDIUM dropped candidates.

---

## 4. Top-5 contradictions (round-2)

**No new contradictions were found in any of the five R2 sweeps for any of topics 01–15.** Each R2 file's §3 ("Contradiction Follow-Up") affirmatively confirms this null. The top-5 list of pre-existing contradictions (carried from R1) and their R2 status:

| # | R1 contradiction | Round-2 status | Where in draft |
|---|-----------------|----------------|----------------|
| 1 | Salsi 2026 — degenerate D4Z4 on ≥10 additional chromosomes (vs canonical 4q/10q-only framing) | No new 2025-11 to 2026-05 paper challenges the Salsi 2026 direction. Hedge already landed at v2 L30. | v2 L30 `@Salsi2026fshd` |
| 2 | Chi 2025 — OR pseudogenisation as sensory-reallocation (vs Gilad 2004 trichromacy trade-off) | No new 2024-2026 paper deepens the Gilad/Chi divergence. Hedge required by R2_AUDIT_PLAN.md §5 edit 3 but **not yet landed** in v2 (REDO scheduled for r2-fix-apply). | v2 L50 — pending |
| 3 | Lalli 2025 — cM/Mb anti-correlation ρ=−0.43 collapses to ρ≈0 after low-callability arm exclusion | No new 2024-2026 paper offers an alternative cM/Mb relationship at PHRs. Limitations clause already cites the honest null. | v2 L52 limitations (iii) |
| 4 | Porubsky 2025 — no whole-genome crossover-SV correlation in CEPH1463 | No new 2024-2026 paper revisits whole-genome crossover-SV correlations. Already addressed in v2 P9 closing sentence. | v2 P9 close |
| 5 | Rodrigues 2024 ("more than half" of subtelomeres) vs Gershman 2022 ("roughly half") on TERRA promoter count | `santagostino2025terra` resolves to a precise count: 39 of 46 (84.8%). This is *greater precision*, not a contradiction — but it shifts the numeric ground for the v2 TAR1/TERRA sentence. | v2 L40 — supplementable |

**No 6th contradiction exists.** R2 reaffirms: the v2 hedge surface is the right one.

---

## 5. Round-2 bibtex additions (full DOIs and PMIDs)

All 8 bibkeys below are now in `paper_prep/synthesis/REFERENCES_v5.bib` (entry count: 372). Each PMID was verified in the source REFRESH_R2 file via PubMed MCP; each DOI was confirmed absent from REFERENCES_v4.bib.

| Bibkey | DOI | PMID | Source REFRESH_R2 |
|--------|-----|------|-------------------|
| `santagostino2025terra` | 10.1261/rna.080790.125 | 41193243 | `REFRESH_R2_foundational.md` |
| `ataei2025line1dj` | 10.1101/gad.351979.124 | 39797762 | `REFRESH_R2_recombination.md` |
| `hao2024snul` | 10.7554/eLife.80684 | 38240312 | `REFRESH_R2_recombination.md` |
| `sasani2026kfam` | 10.64898/2026.03.06.710071 | 41959501 | `REFRESH_R2_recombination.md` |
| `bouquet_ChenCEP164Cilia2025` | 10.64898/2025.12.04.692363 | 41415467 | `REFRESH_R2_3d-bouquet.md` |
| `pedigree_Porsborg2025primaterecom` | 10.1038/s41467-025-65248-3 | 41285744 | `REFRESH_R2_frontier.md` |
| `subtelstruct_Lee2026SEPTIN14` | 10.1186/s13100-026-00394-z | 41699652 | `REFRESH_R2_frontier.md` |
| `pedigree_Joseph2024PRDM9indep` | 10.1073/pnas.2401973121 | 38809707 | `REFRESH_R2_frontier.md` |

**Duplicate suppressed:** `rhie2026rob` (proposed by REFRESH_R2_methods §4) is identical-DOI (10.64898/2026.03.08.710242) to existing v4 entry `acrocentric_Cechova2026`. NOT added.

---

## 6. Hand-off to r2-fix-apply

The companion plan `paper_prep/synthesis/CITATION_UPGRADE_PLAN_v2.md` lists every per-paragraph ADD action for v3 drafting. r2-fix-apply should apply that plan plus the R2_AUDIT_PLAN.md §5 line edits (1-10) to produce `NATURE_DRAFT_v3.md`. The round-1 CITATION_UPGRADE_PLAN.md has already been applied during Pass B (per audit §3); v2 plan is round-2-only and does NOT duplicate any v1 ADD action.

*End of LITERATURE_REFRESH_v2.md. Generated 2026-05-17 by r2-integrator-merge (agent-136).*
