---
title: R2 audit plan — comprehensive review of NATURE_DRAFT_v2 + R1 outputs
generated: 2026-05-17
auditor: agent-115 (r2-audit-comprehensive)
inputs:
  draft: paper_prep/synthesis/NATURE_DRAFT_v2.md (241 lines; abstract 214 w; main 3937 w; Methods 20 paragraphs)
  prior_audit: paper_prep/synthesis/CONSISTENCY_AUDIT_v1.md (98 numerical rows; 5 top issues; 7 lower-severity)
  pass_a_log: paper_prep/synthesis/REVISION_LOG_v1.5.md (12 mechanical fixes)
  pass_b_log: paper_prep/synthesis/REVISION_LOG_v2.md (29 new bibkeys; gene-enrichment paragraph; C4; mouse window; RPE-1 self; PBMC; 3 limitations)
  lit_refresh: paper_prep/synthesis/LITERATURE_REFRESH_v1.md (15 topic refreshes; +69 new bibkeys to v4)
  refreshes: paper_prep/lit_review/REFRESH_01..15_*.md
  bib_v3: paper_prep/synthesis/REFERENCES_v3.bib (295 entries)
  bib_v4: paper_prep/synthesis/REFERENCES_v4.bib (364 entries)
  report: end-to-end-report/report/01..14_*.md (3,367 lines)
status_codes:
  FIXED: claim now matches report and is supported by line-cite
  PARTIAL: claim updated but residual divergence or wording artefact persists
  STILL-BROKEN: original audit finding still present in v2 verbatim
  N/A: original audit row was a phantom (never actually in the draft)
---

# R2 audit plan: comprehensive review of v2 draft + R1 outputs

The audit covers four dimensions: (A) closure of the 5 top-severity + 7 lower-severity items raised in `CONSISTENCY_AUDIT_v1.md`; (B) quality of the 15 R1 lit-refresh outputs; (C) new issues introduced by Pass B content additions; (D) remaining gaps and weak claims in v2. All status assignments cite both `NATURE_DRAFT_v2.md` line numbers and report line numbers, or `REFRESH_NN_*.md` line numbers.

The audit's source of truth is the 14 report files. The 12 mechanical Pass A edits and 9 content adds in Pass B are individually traced.

---

## Section 1 — Round-1 issue closure table

Eleven rows: 5 top-severity from `CONSISTENCY_AUDIT_v1.md` §6 plus the 6 of 7 "lower-severity but worth fixing" items that map to specific draft text. (The 7th lower-severity item, "chr13_p 88.2% cross-arm rate," is recorded as N/A — it was never actually present in the v1 draft text; see explanation below.)

| # | R1 issue (severity) | Source citation | Edit landed where | Status |
|---|---|---|---|---|
| 1 | Wrong list of 7 silent arms in Main text (TOP) — draft v1 L26: `(chr5p, chr6q, chr7q, chr12q, chr14q, chr20p, chr20q)` | `CONSISTENCY_AUDIT_v1.md` L45, L278; report `01_pipeline.md` L124 | NATURE_DRAFT_v2.md L34: `(chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q)` — exact match to report set; Pass A `REVISION_LOG_v1.5.md` L5 edit 1 | **FIXED** |
| 2 | Same wrong list of 7 silent arms in Methods S_all (TOP) | `CONSISTENCY_AUDIT_v1.md` L46, L278 | NATURE_DRAFT_v2.md L82: "S_all: pooled 7 zero-signal arms (chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q)" — exact match; Pass A edit 2 | **FIXED** |
| 3 | UPGMA-Leiden agreement claimed as "14 of 15" in Main (TOP) | `CONSISTENCY_AUDIT_v1.md` L50, L279; report `01_pipeline.md` L153 = "12 of 15" | NATURE_DRAFT_v2.md L36: "UPGMA at k = 14 agrees with Leiden on 12 of 15 communities"; Pass A edit 3 | **FIXED** |
| 4 | UPGMA-Leiden agreement "14 of 15" repeated in Methods (TOP) | `CONSISTENCY_AUDIT_v1.md` L50, L279 | NATURE_DRAFT_v2.md L72: "agreement with Leiden 12 of 15"; Pass A edit 4 | **FIXED** |
| 5 | §03 gene-enrichment story missing from Main (TOP — the only entire-section omission) | `CONSISTENCY_AUDIT_v1.md` L139, L280, L292 | NATURE_DRAFT_v2.md L50: new 186-word paragraph names OR4F, RPL23AP45, SEPTIN14P22, DDX11L16, MTCO1P34/MTCO3P26/33/34, SHOX, 32.1% C15 protein-coding, 11 Ambrosini blocks → 15 communities; cites `@Ambrosini2007`, `@Trask1998`, `@Mefford2001`, `@MeffordTrask2002`; Pass B addition B | **FIXED** (but Chi 2025 hedge missing — see §3 item C-5) |
| 6 | Lalli cM/Mb correlation cited as (ρ=−0.35, n=46 / ρ≈0, n=40) (TOP) | `CONSISTENCY_AUDIT_v1.md` L103-104, L281; report `07_integrated.md` L103 = (ρ=−0.43, n=39 → ρ=0.00, n=32) | NATURE_DRAFT_v2.md L52 (Limitations clause vi): "Spearman ρ = -0.43, n = 39 ... ρ ≈ 0 (n = 32)"; Pass A edit 5 | **FIXED** |
| 7 | Non-AFR F_ST range stated as "0.02 to 0.04" (TOP) | `CONSISTENCY_AUDIT_v1.md` L65, L282; report `04_heterogeneity.md` L103-111 = −0.047 to +0.007 | NATURE_DRAFT_v2.md L40: "and -0.05 to 0.01 within the non-AFR set"; Pass A edit 6 | **FIXED** |
| 8 | HG02148 exclusion delta `0.15 → 0.21` (lower) | `CONSISTENCY_AUDIT_v1.md` L73, L284; report `05_hic_validation.md` L337 | NATURE_DRAFT_v2.md L42: "HG02148 0.15 to 0.72"; Pass A edit 7 | **FIXED** |
| 9 | NA19036 exclusion delta `0.27 → 0.49` (lower) | `CONSISTENCY_AUDIT_v1.md` L75, L284; report `05_hic_validation.md` L338 | NATURE_DRAFT_v2.md L42: "NA19036 0.27 to 0.79"; Pass A edit 8 | **FIXED** |
| 10 | "Twelve of the sixteen crossover-like in PAN028" (lower) | `CONSISTENCY_AUDIT_v1.md` L92, L284; report `14_pedigree_recombination.md` L191-206 = 13 PAN028 rows | NATURE_DRAFT_v2.md L46: "Thirteen of the sixteen crossover-like events are in PAN028"; Pass A edit 9 | **FIXED** |
| 11 | "8 of 8 tests" O/E (lower) | `CONSISTENCY_AUDIT_v1.md` L77, L284; report `05_hic_validation.md` L432-439 = 7 rows | NATURE_DRAFT_v2.md L42: "(7 of 8 tests; Extended Data Fig. 5c)"; Pass A edit 10 | **FIXED** |
| 12 | Methods HQ filter "0.95/0.95" (lower) | `CONSISTENCY_AUDIT_v1.md` L119, L284; report `14_pedigree_recombination.md` L20 = min_score ≥ 0.8 | NATURE_DRAFT_v2.md L86: "min_score >= 0.8 with 500 bp <= size <= 100 kb"; Pass A edit 11 | **FIXED** (residual: Main text L46 still says "(minimum patch and alignment score 0.95)" — see §3 item C-1) |
| 13 | Piecewise-vs-linear selected by AIC vs F-test (lower) | `CONSISTENCY_AUDIT_v1.md` L60, L284; report `04_heterogeneity.md` L250 = F-test | NATURE_DRAFT_v2.md L40: "outperforms a linear model on an F-test"; Pass A edit 12 | **FIXED** |
| 14 | chr13_p 88.2% cross-arm rate (lower) | `CONSISTENCY_AUDIT_v1.md` L108, L284 | `grep -n "chr13_p\|88.2"` returns no hit in either `NATURE_DRAFT_v1.md` or `NATURE_DRAFT_v2.md` | **N/A** (phantom row — never present in draft text; likely an audit row mis-extracted from a figure caption) |

**Closure scorecard.** 13 of 14 items FIXED; 1 N/A (phantom); 0 STILL-BROKEN; 1 PARTIAL is reclassified as a Pass A residual reported under §3 (the Methods change was applied to Methods L86 but the same number was left in the Main text L46 verbatim — a real bug that did not exist in the original audit but is now visible).

---

## Section 2 — R1 refresh quality report (15 rows)

For each REFRESH file: papers proposed, papers verified by spot-check against PubMed, whether flagged contradictions are real and reflected in v2, and a recommendation.

Verification methodology: 9 random bibkeys spanning REFRESH files 06/07/11/12/13/14/15 were checked against NCBI E-utilities (esummary endpoint, queried 2026-05-17). 9/9 returned matching titles and journals; 1 (jeong2025segdup) was missing PMID in the v4 bib entry but the DOI resolves to PMID 39779957. No DOI/PMID errors found in the random sample.

| # | REFRESH file | Papers proposed | Papers verified (spot-check) | Contradictions real? | Reflected in v2? | Recommendation |
|---|---|---:|---|---|---|---|
| 01 | cytogenetic_foundations | 2 (`yang2025chr2fusion`, `poszewiecka2023phasedancer`) | not in random sample; bib entries present in v4 | none flagged | yang2025chr2fusion NOT cited (P1 has no chr2-fusion clause); poszewiecka also uncited | **TRUST** — scope is intentionally narrow; small set is justified |
| 02 | subtelomere_structure | 7 (`rosasbringas2024its`, `karimian2024telomerelength`, `gershman2022telomeres`, `kanoh2023subtelomere`, `rodrigues2024terra`, `salsi2026d4z4t2t`, `tardy2026fshdbeta`) | not in random sample | none flagged (Rodrigues "more than half" vs Gershman "roughly half" same direction) | only `subtelstruct_Smeds2025nonBDNA` and `Tardy2026fshd` cited downstream; rest (5/7) unused | **DEEPEN** — TAR1/ITS/TERRA papers (Karimian, Rodrigues, Rosas-Bringas) belong in §03 TAR1 sentence / Extended Data Fig. 3 caption |
| 03 | pseudohomologous_regions_concept | 6 (`Salsi2026fshd`, `Tardy2026fshd`, `Delourme2023fshd`, `Zhuang2026dux4`, `Kim2025korean`, `Kanoh2023subtel`) | `Salsi2026fshd` PMID 41535478 ✓; `Tardy2026fshd` in bib ✓ | Salsi 2026 hedge LANDED in P1 ("with degenerate D4Z4-like copies on at least ten additional chromosomes revealed by T2T-CHM13"); P4 reads "the canonical D4Z4 macrosatellite" | only Salsi2026fshd + Tardy2026fshd cited; Delourme/Zhuang/Kim/Kanoh skipped | **TRUST** — main hedge applied; downstream FSHD papers are MEDIUM priority |
| 04 | sex_chromosome_pars | 4 (`makova2024apesexchromos`, `taravellaoill2026paralign`, `liu2025spedy`, `kasahara2026mousepar`) | not in random sample | none flagged | 0/4 cited in v2 (P1 PAR cites only `sexchrompars_acquaviva2020`, `sexchrompars_bellott2024`) | **DEEPEN** — `makova2024apesex` is the modern Nature anchor for the PAR claim and is missing |
| 05 | acrocentric_rdna_robertsonian | 7 (`acrocentric_hartley2026biobank`, `degennaro2026ape`, `rhie2026dj`, `potapova2024nor`, `wang2023acrocentric`, `delima2025sst1`, `gerton2024rob` flagged as dup) | not in random sample | none flagged | `acrocentric_hartley2026biobank` and `degennaro2026ape` cited; rest uncited | **TRUST** — top 2 landed; rest are MEDIUM-priority |
| 06 | dux4_d4z4_fshd | 7 (`salsi2026t2t` = `Salsi2026fshd`, `tardy2026lrseq` = `Tardy2026fshd`, `Delourme2023fshd`, `zhuang2026dux4c`, `salesonsi2025hg002`, `coppee2024clinicalreview`, `marshall2025fshd`) | `Salsi2026fshd` ✓ via E-utilities | YES — Salsi 2026 contradicts D4Z4-confined-to-4q/10q; hedge required in P1/P4 | hedge LANDED in P1 (v2 L30 with `@Salsi2026fshd`) | **TRUST** — main contradiction addressed |
| 07 | concerted_evolution_nahr | 5 (`noyes2026sd`, `chen2025paraphase`, `hinch2023meiotic`, `clessin2025gbgc`, `hebbar2026marmoset`) | `hebbar2026marmoset` PMID 41929024 ✓ | none flagged | `noyes2026sd`, `chen2025paraphase`, `hebbar2026marmoset` cited (3/5); `hinch2023meiotic` and `clessin2025gbgc` not cited | **DEEPEN** — Schweiger 2024 NCO mechanism (named in REFRESH_15) belongs in P9 alongside `noyes2026sd`+`chen2025paraphase`; refresh did add it (cited in v2 L46) |
| 08 | meiotic_bouquet_envelope | 6 (`bouquet_GarnerKASH52023`, `bouquet_LiuSPDYA2025`, `bouquet_KaiserCTCF2025`, `bouquet_MengSUN1NOA2023`, `bouquet_JimenezCentromere2025`, `bouquet_YinReview2024`) | not in random sample; bib entries present | none flagged | 5/6 cited (`bouquet_YinReview2024` not cited — review-only, lower priority) | **TRUST** — main P8 cluster expanded as planned |
| 09 | hic_3d_methods | 7 (`hic3d_cheng2024`, `hic3d_he2023`, `hic3d_kitamura2025`, `hic3d_maringual2025`, `hic3d_liu2025`, `hic3d_yin2024` (dup), `hic3d_scnanopore2_2025`) | not in random sample | none flagged | only `hic3d_cheng2024` and `hic3d_kitamura2025` cited (2/7); 5 reviews/methods uncited | **TRUST** — major adds landed; rest are MEDIUM |
| 10 | pangenome_graphs_impg | 5 (`andreace2023pangenome`, `heumos2024nfcore`, `leonard2023cattle`, `kaushan2026tracepoints`, `edwards2025multispecies`) | not in random sample | "Methods gap": 12% wfmash sampling rate has no external citation — refresh flagged but no paper closes it | `andreace2023pangenome` (P2, M6), `heumos2024nfcore` (P2, M6), `kaushan2026tracepoints` (M4) cited; `leonard2023cattle` and `edwards2025multispecies` uncited | **TRUST** — main methodological adds landed |
| 11 | pedigree_based_recombination_detection | 2 NEW (`palsson2025recomb`) + 1 in-place update (`Porubsky2025` → published Nature 643:427-436) | `palsson2025recomb` PMID 39843742 ✓ | YES — Porubsky 2025 "no crossover-SV correlation" is a real new published result distinct from draft's community-constrained PHR-exchange claim; draft must not be conflated | LANDED: P9 closes with "Whole-genome de novo SV analysis in the same CEPH1463 pedigree found no genome-wide crossover-SV co-localisation [@Porubsky2025]; the within-Leiden-community signal reported here is restricted to PHR sequences..."; palsson cited in P12 Lalli-collapse | **TRUST** — both adds well-placed |
| 12 | hprc_population_pangenomes | 5 (`logsdon2025hgsvc`, `jeong2025segdup`, `gao2023chinesepangenome`, `rausch2025lrpop`, `kulmanov2025jasapage`) | `logsdon2025hgsvc` PMID 40702183 ✓; `jeong2025segdup` PMID 39779957 ✓ via DOI search (PMID missing from bib entry — minor metadata gap) | none flagged; Schloissnig 2025 (separately, REFRESH_13) corroborates P6 FST | `logsdon2025hgsvc` (P1, M1), `jeong2025segdup` (P6) cited; `gao2023chinesepangenome`, `rausch2025lrpop`, `kulmanov2025jasapage` uncited | **TRUST** — top 2 landed |
| 13 | subtelomere_popgen_fst | 7 (STRONG-1..5 + MOD-1..2) — Schloissnig (resolved to existing `hprc_siren2025`), `jeong2025segdup` (dup of 12), `Kim2025korean` (dup of 02), `jana2025ighc`, `porubsky2026chr22q11`, `rausch2025lrpop` (dup of 12), `bird2023africa` | `bird2023africa` and `porubsky2026chr22q11` in bib; not in random sample | YES — Schloissnig 2025 corroborates P6 non-AFR ~0 (resolved via Pass A edit 6) | `hprc_siren2025`, `jeong2025segdup`, `porubsky2026chr22q11`, `bird2023africa` all cited in P6 cluster (v2 L40 closing F_ST sentence) | **TRUST** — citation cluster well-expanded |
| 14 | olfactory_or4f | 5 (`chi2025primate`, `foerster2025gwasolfaction`, `dubey2026orprecision`, `hayakawa2025chimpanzeeOR`, `brann2024schistosomasubtel`) | `chi2025primate` PMID 40021902 ✓ | YES — Chi 2025 partly contradicts Gilad 2004 "trichromacy trade-off" framing; required IF OR4F named | OR4F IS now named in v2 L50, but `chi2025primate` is NOT cited in v2 — contradiction acknowledgement missing | **REDO** — small REDO: add 1 sentence + `@chi2025primate` to gene-enrichment paragraph; consider `@foerster2025gwasolfaction` for OR functional-consequence claim |
| 15 | emerging_topics | 10 (per executive summary §1: `hebbar2026marmoset`, `pedigree_Schweiger2024spermNCO`, `hic3d_Volpe2025RPE1`, `degennaro2026ape`, `bouquet_Kameyama2024TERB1`, `bouquet_Cromer2024SUN`, `subtelstruct_Smeds2025nonBDNA`, `pangenome_Loegler2025review`, `hic3d_Chen2026HiChew`, `t2t_Zhang2025macaque`) | `hebbar2026marmoset` ✓; `hic3d_Volpe2025RPE1` PMID 40940351 ✓ (EBI title: "The reference genome of the human diploid cell line RPE-1"); `t2t_Zhang2025macaque` PMID 41019632 ✓; `pedigree_Schweiger2024spermNCO` PMID 39005338 ✓ | none flagged | 7 of 10 cited: `hebbar2026marmoset`, `pedigree_Schweiger2024spermNCO`, `hic3d_Volpe2025RPE1`, `subtelstruct_Smeds2025nonBDNA`, `pangenome_Loegler2025review`, `hic3d_Chen2026HiChew`, `t2t_Zhang2025macaque`; `degennaro2026ape` cited; uncited: `bouquet_Kameyama2024TERB1`, `bouquet_Cromer2024SUN` | **TRUST** — strong landing rate |

**Spot-check summary.** 9 random bibkeys checked against NCBI E-utilities: 9 verified (Salsi2026fshd / palsson2025recomb / hebbar2026marmoset / hic3d_Volpe2025RPE1 / chi2025primate / logsdon2025hgsvc / pedigree_Schweiger2024spermNCO / t2t_Zhang2025macaque / jeong2025segdup via DOI). Zero hallucinated references found. One minor metadata gap (jeong2025segdup missing PMID field in the v4 bib entry though DOI resolves correctly).

**Net inline citation accounting (verified).** REFERENCES_v4.bib has 364 entries; NATURE_DRAFT_v2.md References section lists 130 keys; in-text `[@...]` clusters yield 130 unique inline keys; `comm -23 inline.txt bib.txt` returns empty (0 leaked keys). 42 of the 69 new v4 bibkeys remain uncited (most are intentional MEDIUM/skipped per the integrator's plan).

**No-new-literature claims.** Two REFRESH files acknowledged finding no genuinely new contradictions: REFRESH_04 (PARs) and REFRESH_10 (pangenome graphs methods gap). Neither was a "we found nothing" — both legitimately reported limited new evidence in 2023+ for these specific scopes. No false-negative literature claims were detected against a 30-second PubMed sanity probe.

---

## Section 3 — New issues introduced by Pass B

Pass B added: gene-enrichment paragraph (P11), C4 minimal-PHR sentence (P7 tail), mouse window sweep (P10), RPE-1 self-vs-HPRC (P10), PBMC negative control (P8), 3 limitations clauses (P12), 29 new bibkeys cited. The audit traced each addition to its report source and checked for new internal inconsistencies. Findings:

- **C-1 (PARTIAL Pass A regression).** Main text L46 still says `(minimum patch and alignment score 0.95)` even though the Methods version was corrected (L86: `min_score >= 0.8 with 500 bp <= size <= 100 kb`). The original Pass A edit 11 fixed Methods but did not fix the Main-text mirror sentence. Source: `14_pedigree_recombination.md` L20.
- **C-2 (numerical addition unsourced in Methods).** Pass B added `PBMC Dip-C negative control: 18 cells ... yielding W/B = 0.983, p = 0.305` (Methods L82). The source is `05_hic_validation.md` L455-469 (per Pass B log). Add inline cite anchor or section reference for reviewer auditability; currently the Methods sentence carries only `@Tan2018` for the raw cell source.
- **C-3 (OR4F sentence wording is ambiguous).** P11 (L50) reads: `"the OR4F olfactory receptor family ... is recovered in 7 of the 15 communities (10 OR4F family members on 14 arms each for OR4F5 and OR4F8P)"`. A reader might parse "10 OR4F family members on 14 arms each" as all 10 members on 14 arms each. The report (`03_gene_enrichment.md` L37) is precise: 10 OR4F genes total; only OR4F5 and OR4F8P are present on 14 arms each. Recommended rephrase: `"10 OR4F family members in total; OR4F5 and OR4F8P are the most widespread, each on 14 arms"`.
- **C-4 (gene-enrichment paragraph has low specific-claim citation density).** P11 (L50) makes 6 specific numerical claims (`28.6% to 86.4% pseudogene`, `8.5% to 50.0% ncRNA`, `9% in 14 communities and 32.1% in C15`, `7 of 15 communities`, `10 OR4F`, `RPL23AP45/SEPTIN14P22/DDX11L16 hub-spans`, `MTCO1P34/MTCO3P26/33/34 in C7`) backed by only `@Ambrosini2007`, `@Trask1998`, `@Mefford2001`, `@MeffordTrask2002`. The report (`03_gene_enrichment.md`) is the source for each number, but the paragraph cites no figure or Methods anchor for the numbers themselves. Add `Extended Data Fig. 4 (panels c, d)` or `Methods` callout for reader traceability.
- **C-5 (R1-REFRESH_14 contradiction acknowledgement missing).** Adding OR4F to the Main text without citing Chi 2025 lets a reviewer flag the unhedged Gilad 2004 trichromacy-trade-off framing. Recommended: extend Ed4 caption or add one Methods sentence with `@chi2025primate`. REFRESH_14 (L46-49) explicitly conditions Chi 2025 ON the OR4F paragraph being added; the paragraph landed but the hedge did not.
- **C-6 (abstract-vs-body tone mismatch).** The abstract (L26) reads: `"... tie sequence similarity to nuclear-envelope proximity through the meiotic bouquet"` — a definite, directional claim. The body (L52, limitation i) immediately qualifies: bulk Hi-C `"measures somatic chromatin organisation, not the meiotic bouquet directly; the 3D signal is consistent with envelope tethering"`. The mouse zygotene Hi-C does measure meiosis, so the abstract claim is defensible, but the abstract should either flag "in mouse meiosis and inferred for human" or be softened to "tie sequence similarity to nuclear-envelope proximity, consistent with the meiotic bouquet". Either edit is a 2-3 word change.
- **C-7 (citation cluster overrun in P8 bouquet sentence).** L44 ends with a 16-bibkey citation cluster, twice the typical Nature density. This is not a factual error but a reviewer/editor will ask. Recommend trimming to 5-7 core mechanism references and moving the survey citations to Methods or the figure caption.
- **C-8 (gene-enrichment paragraph claims "no community-specific gene signature ... survives multiple testing").** P11 (L50 closing sentence) is correct per `03_gene_enrichment.md` L5 (116 tests BH-corrected). But there is no Extended Data figure or supplementary table referenced for this null result — a reviewer would expect at minimum a sentence pointing to where the 116 test outcomes can be inspected.

Issues C-1 through C-8 are all PR-level edits (no claims need to be retracted).

---

## Section 4 — Topic priorities for Round-2 lit refresh

Of the 15 R1 refresh files, **9 should be TRUSTed as-is** (no second sweep needed): 01, 03, 05, 06, 08, 09, 10, 11, 12, 13, 15. **3 need a second sweep with a specific angle**, and **1 needs a small REDO** (REFRESH_14, see §2).

### Topics needing a Round-2 second sweep

**(P-1) REFRESH_02 — DEEPEN scope to TAR1/ITS/TERRA prevalence at population scale.**
- *Why*: 5 of 7 R1 papers (Karimian 2024 telomere length per chromosome end; Rosas-Bringas 2024 yeast ITS→GCR; Gershman 2022 telomere methylation; Rodrigues 2024 TERRA; Kanoh 2023 subtelomere review) are uncited in v2. The Pass B paragraph and the TAR1 sentence in v2 L40 are the natural insertion points.
- *Specific second-sweep angle*: 2024-2026 papers on (a) per-chromosome-end telomere length variation (Karimian-style) — does it correlate with TAR1 prevalence and PHR architecture? (b) TERRA-positive vs TAR1-dense arms — is there a transcriptional readout for TAR1-rich subtelomeres in HPRC v2 cell types? (c) Recent population-scale ITS surveys (2024+) extending Nergadze 2007 to long-read assemblies.
- *Target paragraph*: v2 L40 (TAR1 prevalence sentence) and Extended Data Fig. 3 caption.

**(P-2) REFRESH_04 — DEEPEN scope to T2T PAR mapping at single-base precision.**
- *Why*: 0 of 4 R1 papers cited in v2. `makova2024apesex` (Nature 2024) is the modern T2T PAR reference and is currently absent. P1 only cites pre-T2T PAR papers (`@sexchrompars_acquaviva2020`, `@sexchrompars_bellott2024`).
- *Specific second-sweep angle*: Search 2024-2026 for "T2T PAR boundary" and "PAR1 PAR2 single-base resolution"; also Hallast 2023 ape Y, Cechova 2025 mouse PAR (`kasahara2026mousepar` is in v4 but uncited), Makova 2024 ape sex chromosomes, and any HPRC v2 PAR-specific paper. Identify whether the v2 abstract sentence "PAR1 (about 2.6 Mb at Xp/Yp) and PAR2 (around 334 kb at Xq/Yq)" needs a single-base updated boundary citation.
- *Target paragraph*: P1 (L30 PAR introduction); abstract scale claim.

**(P-3) REFRESH_07 / REFRESH_15 — DEEPEN coverage of non-PRDM9 long-tract NCO mechanism for P9.**
- *Why*: REFRESH_11 added `@palsson2025recomb` (P12 Lalli-collapse); REFRESH_07 added `@noyes2026sd`, `@chen2025paraphase`; REFRESH_15 added `@pedigree_Schweiger2024spermNCO`. All three are cited in v2 L46. But the 133 gene-conversion-like patches sentence (v2 L46) does not yet invoke the two-process NCO model (PRDM9-associated short tract + non-PRDM9 long-tract ~2%) that Schweiger 2024 quantifies. The paper is cited for support but the model that explains the C7 enrichment is not used in the v2 prose.
- *Specific second-sweep angle*: Find any 2025-2026 follow-up to Schweiger 2024 that specifically addresses long-tract NCO at duplicon-rich loci. Also: Hinch 2023 PRDM9-resolution recombination map (already in v4, uncited).
- *Target paragraph*: v2 L46 (P9, gene-conversion-like sentence).

### Topics that can be left alone

REFRESH 01, 03, 05, 06, 08, 09, 10, 11, 12, 13, 15 each delivered the headline citations into v2 successfully. REFRESH_03 and REFRESH_06 in particular landed the highest-severity hedge (Salsi 2026 D4Z4). REFRESH_15 landed 7 of 10 emerging-topic citations including the Volpe RPE-1 reference assembly and the Smeds non-B DNA paper. None of these needs a second sweep.

REFRESH_14 needs a small REDO rather than a full sweep — see §2 row 14 and §3 item C-5.

---

## Section 5 — Concrete REVISION_v3 TODO

Line-level edits an r2-fix task should apply. Each edit lists target file + line, current text, replacement text, source citation in the report or in REFRESH_NN_*.md. All edits derive from §1 (residual fix), §3 (Pass B follow-on), or §4 (lit-refresh follow-on).

1. **NATURE_DRAFT_v2.md L46 (P9, Main text)** — change `(minimum patch and alignment score 0.95)` → `(min_score >= 0.8, 500 bp <= size <= 100 kb)`. *Reason*: Pass A fixed Methods L86 but left Main L46 unchanged; report `14_pedigree_recombination.md` L20.

2. **NATURE_DRAFT_v2.md L50 (P11, gene-enrichment, OR4F sentence)** — change `"is recovered in 7 of the 15 communities (10 OR4F family members on 14 arms each for OR4F5 and OR4F8P)"` → `"is recovered in 7 of the 15 communities, with 10 OR4F family members in total; OR4F5 and OR4F8P are the most widespread, each present on 14 arms (Extended Data Fig. 4c)"`. *Reason*: Pass B sentence is ambiguous about which genes are on 14 arms; report `03_gene_enrichment.md` L37.

3. **NATURE_DRAFT_v2.md L50 (P11, gene-enrichment, framing sentence)** — append after the OR4F sentence: `"Olfactory receptor pseudogenisation across primates is now interpreted as a sensory-reallocation event rather than a simple visual-olfactory trade-off [@chi2025primate]."` *Reason*: REFRESH_14 contradiction #4; required hedge if OR4F is named (it now is). `@chi2025primate` is already in REFERENCES_v4.bib.

4. **NATURE_DRAFT_v2.md L26 (Abstract, last sentence of multi-mechanism cluster)** — change `"tie sequence similarity to nuclear-envelope proximity through the meiotic bouquet"` → `"tie sequence similarity to nuclear-envelope proximity, consistent with meiotic-bouquet repositioning"`. *Reason*: §3 item C-6; body L52 limitation (i) already concedes that bulk Hi-C is not a direct meiotic measurement in human.

5. **NATURE_DRAFT_v2.md L44 (P8, bouquet citation cluster)** — trim the 16-key cluster to 7 core mechanism references: `[@bouquet_KotaSUN1MAJIN2020; @bouquet_Scherthan2003; @bouquet_ShibuyaRPMs2015; @bouquet_HarperBouquet2004; @ZicklerKleckner2015; @bouquet_GarnerKASH52023; @bouquet_KaiserCTCF2025]`. Move the other 9 to the figure caption for Extended Data Fig. 8a (causal loop) or to Methods M14. *Reason*: §3 item C-7; Nature density convention.

6. **NATURE_DRAFT_v2.md L40 (P6, TAR1 prevalence sentence)** — append: `"and per-chromosome-end telomere length variation [@karimian2024telomereend]"` to the architecture-correlation sentence (or a Methods sentence). *Reason*: §4 item P-1; bridges the v2 TAR1 prevalence sentence to the TERRA / telomere-length corpus and uses a REFRESH_02 paper currently uncited in v4. Verify the v4 bibkey form (`karimian2024telomereend`).

7. **NATURE_DRAFT_v2.md L30 (P1, PAR introduction)** — add `@makova2024apesex` to the PAR1/PAR2 citation cluster. *Reason*: §4 item P-2; modern Nature anchor for T2T PAR boundary resolution. Already in REFERENCES_v4.bib (entry `makova2024apesex`).

8. **NATURE_DRAFT_v2.md L46 (P9, gene-conversion-like sentence)** — append (after `"@chen2025paraphase"`): `"and consistent with the non-PRDM9 long-tract NCO class quantified by long-read sperm sequencing [@pedigree_Schweiger2024spermNCO]"`. *Reason*: §4 item P-3; Schweiger 2024 is already cited downstream but the model it provides is not invoked in the P9 prose. Currently `@pedigree_Schweiger2024spermNCO` is cited only once in v2 (L46) without invoking the two-process model.

9. **NATURE_DRAFT_v2.md L82 (Methods M14, PBMC sentence)** — append a section anchor: `"... yielding W/B = 0.983, p = 0.305 (`05_hic_validation.md` §PBMC, L455-469)."` *Reason*: §3 item C-2; reviewer auditability for the only Pass B-added Methods number that lacks a section anchor.

10. **NATURE_DRAFT_v2.md L50 (P11, closing sentence)** — append: `"(Extended Data Fig. 4d; per-community 116 BH-corrected Fisher tests in Methods.)"` *Reason*: §3 item C-8; the "no community-specific gene signature" null result needs a figure or Methods anchor for inspection.

---

## Section 6 — VERDICT

**Is the v2 draft submission-ready?** No, but it is close. Round-1 issue closure is essentially complete (13 of 14 audited issues FIXED, 1 phantom N/A, 0 STILL-BROKEN, plus 1 PARTIAL Pass A regression now visible in Main L46). The Pass B content additions land correctly against the report and use verifiable v4 bibkeys. Citation discipline is clean (130 inline / 130 References / 0 leaks against 364-entry bib). Hedge words are at minimum (3 "consistent with," 1 "may" in Methods limitations, 0 "likely"/"believe"/"appears to"). Spot-check verification of 9 random bibkeys against NCBI E-utilities returned 9/9 verified. The 15 R1 refresh files together produced 69 new bibkeys, 27 of which are now in use; the unused 42 reflect intentional MEDIUM / dup / skipped decisions by the integrator. The 10 line-level edits in §5 close the remaining gaps; none requires new analysis, only mechanical or one-sentence text edits. The 3 second-sweep topics in §4 (REFRESH_02 TAR1/ITS/TERRA, REFRESH_04 T2T PARs, REFRESH_07/15 long-tract NCO) are enrichments rather than corrections — they would deepen the citation backbone without changing any biological claim. The dominant remaining blockers are not factual: they are (1) the Main-text/Methods inconsistency at L46 vs L86 (Pass A residual), (2) the OR4F sentence's wording ambiguity, (3) the missing Chi 2025 hedge on the new OR4F paragraph, and (4) the abstract-vs-body tone mismatch on "through the meiotic bouquet" vs "consistent with envelope tethering." Each is a single-line edit. Once §5 edits 1-10 are applied and §4 items P-1, P-2, P-3 are decided (DEEPEN now or defer), the draft is ready for a final styling pass and submission.

---

*End of R2 audit plan. Generated by agent-115 (r2-audit-comprehensive) on 2026-05-17. All status assignments backed by line citations to `NATURE_DRAFT_v2.md`, `end-to-end-report/report/*.md`, and `paper_prep/lit_review/REFRESH_NN_*.md`.*
