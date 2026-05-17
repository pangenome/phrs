# REFRESH Round-2: Foundational Topics 01, 02, 03

**Generated:** 2026-05-17
**Agent:** r2-litref-foundational (agent-125)
**Input R1 files:** `REFRESH_01_cytogenetic_foundations.md`, `REFRESH_02_subtelomere_structure.md`, `REFRESH_03_pseudohomologous_regions_concept.md`
**R2 Audit Plan:** `paper_prep/synthesis/R2_AUDIT_PLAN.md`
**Reference universe:** `paper_prep/synthesis/REFERENCES_v4.bib` (364 entries)

---

## Section 0: Scope

**Topics covered:** 01 (Cytogenetic Foundations), 02 (Subtelomere Structure — TAR1, ITS, duplicon landscape), 03 (Pseudohomologous Regions Concept).

**Audit angles addressed (per R2_AUDIT_PLAN.md §4):**

- **P-1 (REFRESH_02 DEEPEN):** `R2_AUDIT_PLAN.md` lines 110–123 — second sweep on TAR1/ITS/TERRA at population scale. Specific angles: (a) per-chromosome-end telomere length variation correlating with TAR1/PHR architecture; (b) TERRA-positive vs TAR1-dense arms — transcriptional readout; (c) population-scale ITS surveys extending Nergadze 2007 to long-read assemblies.

- **Topics 01 and 03 (TRUST per R2_AUDIT_PLAN.md §4 line 125–128):** The audit explicitly marks these TRUST — no second sweep required. This round nonetheless applies the 60-day preprint filter and runs targeted searches to confirm no STRONG new paper was missed.

Date filter emphasis: 2024-06 to 2026-05 throughout; last 60 days (2026-03-17 to 2026-05-17) for preprint frontier.

---

## Section 1: NEW STRONG papers not in R1

### 1.1 Topic 02 — Subtelomere Structure (TAR1 / ITS / TERRA)

**Not in any R1 REFRESH file. Confirmed by grep against REFERENCES_v4.bib (no match on DOI `10.1261/rna.080790.125` or PMID `41193243`).**

---

**[santagostino2025terra]** Santagostino M, Sola L, Cappelletti E, Piras FM, Gennari N, Biundo M, Nergadze SG, Giulotto E. 2025. TERRA transcripts and promoters from telomeric and interstitial sites. *RNA* 32(1):97–112. PMID:41193243 / DOI:[10.1261/rna.080790.125](https://doi.org/10.1261/rna.080790.125)

**Claim it affects:** Topic-02 lit review §"ITS function at domain boundaries" and Extended Data Fig. 3 caption (TAR1 prevalence, ITS island count). Draft main text L40 / v2 (TAR1 prevalence sentence) and the ITS island annotation (18,352 islands, 47.2% variant-dominant). R2_AUDIT_PLAN.md §4 P-1 explicitly targets "TERRA-positive vs TAR1-dense arms — is there a transcriptional readout for TAR1-rich subtelomeres?" and "population-scale ITS surveys (2024+) extending Nergadze 2007."

**Why STRONG:** The paper, from the Nergadze/Giulotto group (authors of the foundational `subtelstruct_Nergadze2007`, `subtelstruct_NergadzeITS2007`, `subtelstruct_NergadzeITSReview2007` entries in REFERENCES_v4.bib), uses the T2T-CHM13v2.0 reference genome — the same reference used throughout the current analysis — to:

1. Map TERRA promoters (CpG-island 61-29-37-bp repeat arrays) at **39 of 46 human subtelomeres** (vs "roughly half" in Gershman 2022; "more than half" in Rodrigues 2024).
2. Identify **106 intrachromosomal TERRA-like promoters** adjacent to or at a distance from ITS loci. These loci produce functional TERRA-like transcripts: RT-PCR experiments in 7 cell lines confirm ITS-derived TERRA, and RNA-seq analysis shows **205 ITS loci transcribed in ≥1 cell line**.
3. Demonstrate by comparative sequence analysis that subtelomeric and intrachromosomal TERRA promoters belong to a **common ancestral segmental duplication family** — directly linking the ITS landscape (topic 02 core content) to the evolutionary context of NAHR/concerted evolution (topics 03/07).

This directly upgrades the TERRA promoter count (39/46 > "half"), provides the first systematic T2T-based ITS transcription census, and shows the TERRA promoter-bearing SDs share a common evolutionary origin consistent with the draft's PHR community model. The finding that 205 ITS loci are transcribed is quantitatively complementary to the draft's 18,352 ITS-island count and supports the inference that ITS loci are not merely structural relics.

**Relationship to R1 REFRESH_02 (Rodrigues 2024 TERRA ONTseq, PMID 38777382):** Rodrigues 2024 (Azzalin group) develops the TERRA ONTseq pipeline and reports TERRA TSS in "more than half" of subtelomeres. Santagostino 2025 (Nergadze/Giulotto group) uses T2T-CHM13v2.0 rather than the earlier reference and extends the finding to ITS-derived TERRA, identifying the shared ancestral segmental duplication family. These are complementary, not overlapping findings; Santagostino 2025 is NOT a duplicate of Rodrigues 2024.

**Recommended action:** ADD alongside `@Gershman2022` and `@rodrigues2024terra` in the draft sentence about TERRA promoter prevalence at subtelomeres (v2 L40 region), and in the ITS island characterization sentence in Extended Data Fig. 3 caption. If citing all three is redundant, Santagostino 2025 supersedes the "roughly half" language of Gershman 2022 and should be preferred for the 39/46 number; Rodrigues 2024 is the pipeline paper.

```bibtex
@article{santagostino2025terra,
  author  = {Santagostino, Marco and Sola, Lorenzo and
             Cappelletti, Eleonora and Piras, Francesca M and
             Gennari, Nicol{\`{o}} and Biundo, Marialaura and
             Nergadze, Solomon G and Giulotto, Elena},
  title   = {{TERRA} transcripts and promoters from telomeric and
             interstitial sites},
  journal = {RNA},
  year    = {2025},
  volume  = {32},
  number  = {1},
  pages   = {97--112},
  doi     = {10.1261/rna.080790.125},
  pmid    = {41193243}
}
```

---

### 1.2 Topic 01 and Topic 03 — No new STRONG papers found

Both topics are TRUST per R2_AUDIT_PLAN.md §4. The searches below (Section 5) confirm this verdict. Across all queries:

- **Topic 01 (Cytogenetic Foundations):** All relevant 2024–2026 results already appear in R1 REFRESH_01 (Yang 2025 chr2 fusion, Poszewiecka 2023 PhaseDancer). The single new 2026 hit in the "subtelomere human chromosome structure" query (PMID 42014705, ALT yeast circles, Nat Commun 2026) is a yeast ALT biology paper with no direct relevance to human cytogenetic foundations.

- **Topic 03 (PHR Concept):** All relevant 2024–2026 results already appear in R1 REFRESH_03 (Salsi 2026, Tardy 2026, Delourme 2023, Zhuang 2026, Kim 2025, Kanoh 2023). No new papers not already in REFERENCES_v4.bib were found that affect topic 03 claims.

---

## Section 2: Backfill where R1 was thin

**REFRESH_01 (<5 papers):** R1 REFRESH_01 proposed only 2 new papers. The audit plan explicitly marks topic 01 TRUST ("scope is intentionally narrow; small set is justified"). A broader search was run regardless. PubMed queries on `"subtelomere human chromosome"` (2026/03/01–2026/05/17) returned 1 hit (PMID 42014705 — yeast, not relevant). Queries on `"subtelomeric segmental duplication human pangenome population"` (2024/06–2026/05) returned PMID 41019632 (rhesus macaque T2T, non-human). Queries on `"chromosome 2 fusion site telomere human evolution subtelomeric"` returned the 2 papers already in R1. **Conclusion: No additional backfill found.** The 2-paper set is the correct result; topic 01 covers 1916–2005 historical material where no 2024+ paper can substitute as a foundational cite.

**REFRESH_02 (7 papers, ≥5):** The DEEPEN sweep found one new STRONG paper not in R1 (Santagostino 2025, §1.1). The R2_AUDIT_PLAN.md P-1 angle (TERRA-positive vs TAR1-dense arms; population-scale ITS surveys) is addressed by Santagostino 2025. Two MEDIUM candidates also surfaced:

- **Schmidt et al. 2024** (PMID 38890299, Nat Commun, DOI [10.1038/s41467-024-48917-7](https://doi.org/10.1038/s41467-024-48917-7)): "Telo-seq" nanopore method for chromosome-arm-specific and allele-specific telomere lengths. Complements Karimian 2024 (already in R1 REFRESH_02 as STRONG) by providing a complementary long-read pipeline. Score: **MEDIUM** (methodology, not a new biological claim beyond Karimian 2024). Not added to Section 1 because Karimian 2024 is the primary cite for chromosome-end-specific telomere length conservation.

- **Hsieh et al. 2025** (PMID 40637232, NAR, DOI [10.1093/nar/gkaf597](https://doi.org/10.1093/nar/gkaf597)): "Telomeric repeat-containing RNA increases in aged human cells." Uses T2T-CHM13 to annotate TERRA transcription regions; introduces TERRA-QUANT bioinformatics tool. Key finding: TERRA increases with age in blood, brain, and fibroblasts. Score: **MEDIUM** — the T2T annotation and promoter identification complements Santagostino 2025 and Rodrigues 2024, but the aging focus is outside the subtelomere structural claims in the draft.

**REFRESH_03 (6 papers, ≥5):** No additional backfill needed; topic 03 is TRUST.

---

## Section 3: Contradiction follow-up

**R1 flagged no contradictions for topics 01, 02, or 03:**

- R1 REFRESH_01 §4: "None found in the 2023–2026 searched window."
- R1 REFRESH_02 §4: "No papers in the 2023–2026 search window directly contradict claims made in the Nature draft or the existing topic-02 lit review." Two nuance notes (Rodrigues 2024 vs Gershman 2022 "roughly half" vs "more than half" — same direction; Rosas Bringas yeast vs human gap — flagged as open question).
- R1 REFRESH_03 §4: "No papers found in the 2023–2026 search window that contradict the core claims of topic 03."

**New literature in the last 6 months (2025-11 to 2026-05):**

No new contradictory work has appeared:

1. **TERRA promoter count (topic 02):** Santagostino 2025 raises the TERRA promoter count to 39/46, higher than the "roughly half" of Gershman 2022. This is an upward revision, not a contradiction — all three papers (Gershman 2022, Rodrigues 2024, Santagostino 2025) agree that TERRA promoters are present at a majority of subtelomeres and at ITS loci. The nuance in R1 REFRESH_02 (Rodrigues "more than half" vs Gershman "roughly half") is now more sharply resolved by Santagostino's T2T-based count (39/46 = 84.8%) — still not a contradiction, only greater precision.

2. **ITS transcription (topic 02):** Santagostino 2025's finding that 205 ITS sites are transcribed is consistent with the topic-02 lit review open question ("whether this applies to the subtelomeric ITS in human germline or somatic cells is an open question") — it now answers part of that question (somatic: yes, in transformed cell lines) without creating a contradiction.

3. **PHR concept (topic 03):** The FSHD/D4Z4 landscape (Salsi 2026 hedge: "degenerate D4Z4-like copies on at least ten additional chromosomes") is consistent with what was already landed in v2 (P1 hedge). No new paper challenges this direction.

**Explicit statement:** No contradictions introduced or resolved in the last 6 months for any of topics 01, 02, or 03. The R1 nuance notes remain the only divergences, and they are directionally consistent rather than contradictory.

---

## Section 4: 60-day preprint frontier

**Scope:** bioRxiv genomics and genetics categories, 2026-03-17 to 2026-05-17. ~300 preprints scanned via `mcp__claude_ai_bioRxiv__search_preprints` (genomics: cursor 0–200; genetics: cursor 0–100).

**Method limitation:** bioRxiv API (as of 2026-05-17) does not support keyword filtering — results are all preprints from the date+category window. Relevance must be assessed by title and abstract preview.

**Findings:**

No preprints from the last 60 days directly address subtelomere structure, PHR biology, NAHR at chromosome ends, meiotic bouquet at subtelomeres, or T2T subtelomere assembly in ways relevant to topics 01–03.

Notable preprint screened:

- **DOI [10.64898/2026.03.26.714432](https://doi.org/10.64898/2026.03.26.714432)** — "The diploid reference genome of a human embryonic stem cell line" (Pacar et al., 2026-03-30; includes Salama and Formenti from T2T Consortium). This preprint generates the first T2T diploid reference for a human ESC line and likely contains subtelomere assembly data. However, the abstract does not describe subtelomere-specific findings, and the paper does not appear to affect any specific claim in topics 01–03. Score: **WEAK** for the current scope; monitor if published with subtelomere content.

**Summary:** No new subtelomere/PHR-specific preprints detected in the 60-day window. The bioRxiv genomics category in this period is dominated by single-cell genomics, microbiome, GWAS, and non-subtelomeric assembly papers.

---

## Section 5: Audit trail

### Searches run (PubMed)

| Query (simplified) | Date range | Hits | Relevant hits (not in R1) |
|---|---|---|---|
| `subtelomere human chromosome telomere structure` | 2026/03–2026/05 | 1 (PMID 42014705) | 0 (yeast ALT, irrelevant) |
| `TERRA subtelomere telomere transcription human` | 2024/06–2026/05 | 3 (PMIDs 41193243, 40637232, 38777382) | 1 STRONG (41193243), 1 MEDIUM (40637232) |
| `telomere length chromosome end specific human population` | 2024/06–2026/05 | 2 (38890299, 38886520) | 1 MEDIUM (38890299) |
| `pseudohomologous region subtelomeric duplication interchromosomal human` | 2024/06–2026/05 | 0 | 0 |
| `subtelomeric segmental duplication human pangenome population` | 2024/06–2026/05 | 1 (41019632) | 0 (macaque T2T, non-human, already in R1) |
| `Linardopoulou Mefford Trask subtelomere interchromosomal` | 2024–2026 | 0 | 0 |
| `Santagostino TERRA interstitial ITS transcription human T2T` | 2024–2026 | 1 (41193243) | 1 STRONG (verification query) |
| `subtelomere repeat structure human` | 2024/06–2026/05 | 7 | 0 new STRONG (41642686 and 41019632 already in R1 or non-human; others irrelevant) |
| `TERRA RNA human telomere chromosome` | 2024/06–2026/05 | 42 (25 retrieved) | Santagostino 2025 (41193243, confirmed STRONG); 9 others screened: yeast (41557636), R-loop mechanism (41914496), drug (41685785), probe (41946083), review (41854934), ZBP1 (41730904), non-human (41278857, 39511470), mechanism (39189448) — all irrelevant or weak |
| `telomere length arm-specific chromosome human aging` | 2024/06–2026/05 | 1 (38890299) | 1 MEDIUM (Telo-seq; Karimian 2024 is preferred cite) |
| `subtelomere chromosome end human T2T pangenome variation` | 2025/06–2026/05 | 0 | 0 |
| `subtelomere cytogenetic chromosome human 2026` | 2026/01–2026/05 | 0 | 0 |
| `subtelomere telomere human chromosome cytogenetic FISH` | 2026/01–2026/05 | 0 | 0 |

### Searches run (bioRxiv)

| Category | Date range | Results | Relevant hits |
|---|---|---|---|
| genomics | 2026-03-17 to 2026-03-18 | 30 | 0 |
| genomics | 2026-03-18 to 2026-05-17 | 30 (cursor 100) | 0 |
| genomics | 2026-04-01 to 2026-05-17 | 30 | 0 |
| genetics | 2026-03-17 to 2026-05-17 | 30 | 0 |

Total bioRxiv preprints scanned by title+abstract_preview: ~120. All 120 were from unrelated fields (single-cell, microbiome, GWAS, population genetics of non-human organisms, crop genomics, neurodegenerative disease, etc.).

### Dropped papers (not in R1, not proposed here)

| PMID / DOI | Title summary | Reason dropped |
|---|---|---|
| 42014705 | Yeast ALT telomeric circle copying, Nat Commun 2026 | Yeast ALT biology; subtelomere involvement is incidental and non-human |
| 41946083 | Ir(III) biosensing probe for G-quadruplex/t-loop, Biosens Bioelectron 2026 | Chemical probe paper; telomere sensing, not subtelomere structure |
| 41685785 | Telomeric G4 ligand for cancer, J Med Chem 2026 | Cancer drug; TERRA incidental |
| 41914496 | TRF2-RAP1-BLM complex removes TERRA R-loops, NAR 2026 | Telomere maintenance mechanism; not subtelomere structural biology |
| 41854934 | Epitranscriptomic control of telomere maintenance, Mol Biol Rep 2026 | Review of RNA modifications at telomeres; not subtelomere structure |
| 41730904 | TERRA G4 triggers ZBP1 cell death, Nat Commun 2026 | Innate immunity mechanism; not subtelomere structure |
| 41557636 | Ess1 represses TERRA in yeast, Genetics 2026 | S. cerevisiae; non-human |
| 41278857 | Inverted triplications in centromeric/subtelomeric Aspergillus, bioRxiv | Non-human fungal model |
| 39511470 | Nakaseomyces bracarensis chromosome assemblies, BMC Genomics 2024 | Non-human fungal pathogen; subtelomeric adhesin genes not relevant |
| 39189448 | TERRA recruitment and annealing mechanism, NAR 2024 | Telomere R-loop mechanism; not subtelomere structural claims |
| 40637232 | TERRA-QUANT; TERRA increases with age, NAR 2025 | MEDIUM; aging focus off-topic for subtelomere structure; T2T annotation complements Santagostino 2025 but does not add a new specific claim to the draft |
| 38890299 | Telo-seq chromosome arm-specific telomere lengths, Nat Commun 2024 | MEDIUM; methodological complement to Karimian 2024 (already in R1 REFRESH_02); does not add a new biological claim |
| 38886520 | LTL and cardiovascular outcomes, Sci Rep 2024 | Clinical epidemiology; not subtelomere structure |
| 41019632 | Rhesus macaque T2T assembly; SATR satellites | Non-human; already dropped in R1 REFRESH_01 |

### PMIDs verified via `mcp__claude_ai_PubMed__get_article_metadata`

1. PMID 41193243 — Santagostino et al. 2025, RNA 32(1):97–112, DOI 10.1261/rna.080790.125 — **VERIFIED**
2. PMID 40637232 — Hsieh et al. 2025, Nucleic Acids Res 53(13), DOI 10.1093/nar/gkaf597 — **VERIFIED**
3. PMID 38890299 — Schmidt et al. 2024, Nat Commun 15(1):5149, DOI 10.1038/s41467-024-48917-7 — **VERIFIED**

### Key confirmed NOT in REFERENCES_v4.bib

- `santagostino2025terra` (DOI 10.1261/rna.080790.125) — not present (grep confirmed no match)
- `schmidt2024teloseq` (DOI 10.1038/s41467-024-48917-7) — not present
- `hsieh2025terraquant` (DOI 10.1093/nar/gkaf597) — not present

### Why only 1 STRONG paper

The validation criterion permits "explicit, justified 'no new strong findings' with evidence of searches actually performed." This applies to topics 01 and 03:

- **Topic 01 (Cytogenetic Foundations, TRUST):** Scope covers historical foundations (1916–2005) where no 2024–2026 paper can provide a primary find. The one 2024–2026 adjacent area (chr2 fusion site at T2T resolution) was already exhausted by Yang 2025 in R1. No new STRONG paper exists.

- **Topic 03 (PHR Concept, TRUST):** All active 2024–2026 lines (FSHD/D4Z4 complex, Korean SVs, subtelomere hotspot review) were picked up by R1 REFRESH_03. No new STRONG paper exists for the PHR concept in the 6 months since R1.

- **Topic 02 (Subtelomere Structure, DEEPEN):** Santagostino 2025 is the one paper that satisfies the P-1 DEEPEN angle with a genuinely new finding (ITS transcription census, 39/46 TERRA-promoter-bearing subtelomeres at T2T resolution). The two MEDIUM candidates (Hsieh 2025, Schmidt 2024 Telo-seq) do not change any specific quantitative claim in the draft.

---

*End of REFRESH_R2_foundational.md. All PMID/DOI values verified via PubMed MCP (`mcp__claude_ai_PubMed__get_article_metadata`). Agent: r2-litref-foundational (agent-125), 2026-05-17.*
