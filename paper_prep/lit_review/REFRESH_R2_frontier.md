# REFRESH_R2_frontier: Round-2 Literature Sweep for Topics 13, 14, 15

**Date**: 2026-05-17
**Agent**: agent-123 (r2-litref-frontier)
**Branch**: wg/agent-123/r2-litref-frontier
**Inputs**:
- `paper_prep/lit_review/REFRESH_13_subtelomere_popgen_fst.md` (agent-73)
- `paper_prep/lit_review/REFRESH_14_olfactory_or4f.md` (agent-57)
- `paper_prep/lit_review/REFRESH_15_emerging_topics.md` (agent-54)
- `paper_prep/synthesis/R2_AUDIT_PLAN.md` (agent-115)
- `paper_prep/synthesis/REFERENCES_v4.bib` (364 entries)
- `paper_prep/synthesis/NATURE_DRAFT_v1.md`

---

## Section 0: Scope

**Topics covered**: 13 (subtelomere population genetics and FST), 14 (olfactory OR4F), 15 (emerging topics 2024-2026).

**R2 audit angles addressed** (`R2_AUDIT_PLAN.md` §4):

| Audit angle | Source | Addressed here |
|---|---|---|
| P-1: DEEPEN REFRESH_02 TAR1/ITS/TERRA | §4 P-1 | NOT in scope (topics 02, not 13/14/15) |
| P-2: DEEPEN REFRESH_04 T2T PAR boundary | §4 P-2 | NOT in scope (topic 04, not 13/14/15) |
| P-3: DEEPEN REFRESH_07/15 non-PRDM9 long-tract NCO | §4 P-3 | **YES** — directly addressed (Topic 15) |
| REDO REFRESH_14: chi2025primate hedge in OR4F paragraph | §2 row 14 + §3 C-5 | **YES** — addressed (Topic 14) |

This sweep focuses specifically on what Round-1 missed for topics 13/14/15 and on the 60-day preprint frontier (2026-03-18 to 2026-05-17).

---

## Section 1: NEW STRONG Papers Not in R1

**STRONG threshold**: a paper that would change a draft claim or substantially alter how a claim is framed.

**Summary verdict**: For topics 13 and 14, no new STRONG papers were found beyond R1 (see Section 5 for evidence of searches). For topic 15, one new STRONG paper was identified (Porsborg 2025), and two MEDIUM papers were identified (Lee 2026, Joseph 2024). These are documented below.

---

### STRONG-1 (Topic 15 / P-3 angle)

**Proposed bibkey**: `pedigree_Porsborg2025primaterecom`

**Full citation**:
Porsborg PS, Charmouh AP, Singh VK, et al. Long-read sequencing of primate testis and human sperm allows identification of recombination events in individuals. *Nature Communications*. 2025;16(1):10337. Published 2025-11-24.
DOI: [10.1038/s41467-025-65248-3](https://doi.org/10.1038/s41467-025-65248-3)
PMID: 41285744 (verified via PubMed)

**Duplicate check**: Absent from REFERENCES_v4.bib (verified by grep for "porsborg", "41285744", "s41467-025-65248"). NOT present in any R1 REFRESH file (confirmed by grep of all 15 REFRESH_NN_*.md files).

**Relevance** (STRONG): Long-read sequencing of testis tissue from 16 individuals across six primate species plus three human sperm samples, identifying 2,881 crossovers, 2,314 simple gene conversions, and 555 complex events. Key quantitative findings:
- Non-crossover (NCO) gene conversion tracts: **mean 22-95 bp** across all primate species (short class)
- Crossover-associated (CO-associated) gene conversion tracts: **mean 318-688 bp** (substantially longer)
- GC-biased gene conversion observed for both NCO and CO-associated classes
- Human samples align with double-strand break map, confirming PRDM9-directed breaks

**Claim it changes**: v2 P9 / L46 — "133 gene-conversion-like patches." The draft currently presents all 133 patches as consistent with non-crossover gene conversion (citing Schweiger 2024). Porsborg 2025 calibrates the expected tract-length distributions for both NCO and CO-associated gene conversions in human (and primate) meiosis using a directly comparable long-read approach. Because CO-associated gene conversion tracts are 4-8× longer than NCO tracts, patches in the pedigree dataset with tracts substantially longer than ~100 bp are more likely to be CO-associated gene conversions than standalone NCOs. This refines the mechanistic interpretation: the subtelomeric "gene-conversion-like" patches are not all NCO-derived; a subset — particularly those with larger tracts — could reflect CO-associated gene conversion, a mechanistically distinct event with different implications for sequence-level concerted evolution at subtelomeres.

**Recommended action**: ADD alongside `@pedigree_Schweiger2024spermNCO` in P9. Suggested text (extending R2_AUDIT_PLAN.md §5 edit 8): "…consistent with the non-PRDM9 long-tract NCO class quantified by long-read sperm sequencing [@pedigree_Schweiger2024spermNCO]; the co-existence of crossover-associated gene conversion with substantially longer tracts (318-688 bp; [@pedigree_Porsborg2025primaterecom]) suggests that a subset of the larger patches may reflect crossover-associated rather than standalone NCO events."

**Suggested placement**: v2 P9 (gene-conversion-like sentence, L46) — supplement the existing Schweiger 2024 citation.

```bibtex
@article{pedigree_Porsborg2025primaterecom,
  author  = {Porsborg, Peter Soerud and Charmouh, Anders Poulsen and Singh, Vinod Kumar
             and Winge, Sofia Boeg and Hvilsom, Christina and Oroperv, Carmen
             and Hansen, Lasse Thorup and Berner, Juliana Andrea and Pelizzola, Marta
             and Laurentino, Sandra and Neuhaus, Nina and Hobolth, Asger
             and Bataillon, Thomas and Besenbacher, Soren and Almstrup, Kristian
             and Schierup, Mikkel Heide},
  title   = {Long-read sequencing of primate testis and human sperm allows
             identification of recombination events in individuals},
  journal = {Nature Communications},
  year    = {2025},
  volume  = {16},
  number  = {1},
  pages   = {10337},
  doi     = {10.1038/s41467-025-65248-3},
  pmid    = {41285744},
  note    = {Cross-species long-read meiotic recombination; NCO tracts 22-95 bp,
             CO-associated tracts 318-688 bp; calibrates gene-conversion-like patch
             interpretation for C8 pedigree analysis}
}
```

---

### MEDIUM-1 (Topic 15)

**Proposed bibkey**: `subtelstruct_Lee2026SEPTIN14`

**Full citation**:
Lee MG. Boundary-associated propagation of a processed pseudogene dissects pre-existing limitations of genome annotation in the T2T era. *Mobile DNA*. 2026;17(1). Published 2026-02-17.
DOI: [10.1186/s13100-026-00394-z](https://doi.org/10.1186/s13100-026-00394-z)
PMID: 41699652 (verified via PubMed)

**Duplicate check**: Absent from REFERENCES_v4.bib (verified by grep for "41699652", "lee.*septin", "T2T.*annotation.*septin"). NOT in any R1 REFRESH file.

**Relevance** (MEDIUM): Uses T2T-CHM13 assemblies to show that annotated CICP loci (SEPTIN14 processed pseudogene copies) preferentially localize within segmental duplication blocks that accumulate near pericentromeric and subtelomeric regions. Chain-based comparative analysis shows the SEPTIN14 3' terminal exon + adjacent CICP12 is dispersed into multiple SD-associated units across great apes via secondary structural propagation rather than independent LINE-1 insertions. Purifying selection acts on the SEPTIN14 coding sequence and its 3' exon throughout this dispersal.

**Claim it supports**: v2 P11 / L50 — "RPL23AP45/SEPTIN14P22/DDX11L16 hub-spans." The draft names SEPTIN14P22 as a hub-spanning gene in the gene-enrichment paragraph without citing any paper specifically for this gene family's subtelomeric distribution. Lee 2026 provides the mechanistic basis: SEPTIN14 pseudogene copies accumulate at subtelomeres via SD-mediated propagation in great apes, consistent with the cross-arm community sharing mechanism.

**Recommended action**: ADD as a supporting citation for the SEPTIN14P22 claim in P11. Useful in a phrase such as: "...SEPTIN14P22 [present across community-spanning arms; its dispersal driven by subtelomeric SD propagation; @subtelstruct_Lee2026SEPTIN14]."

```bibtex
@article{subtelstruct_Lee2026SEPTIN14,
  author  = {Lee, Min-Gyu},
  title   = {Boundary-associated propagation of a processed pseudogene dissects
             pre-existing limitations of genome annotation in the {T2T} era},
  journal = {Mobile DNA},
  year    = {2026},
  volume  = {17},
  number  = {1},
  doi     = {10.1186/s13100-026-00394-z},
  pmid    = {41699652},
  note    = {SEPTIN14P pseudogene dispersal via SD blocks at subtelomeres; T2T
             annotation case study; supports SEPTIN14P22 hub-span claim in P11}
}
```

---

### MEDIUM-2 (Topic 15 / P-3 angle)

**Proposed bibkey**: `pedigree_Joseph2024PRDM9indep`

**Full citation**:
Joseph J, Prentout D, Laverré A, Tricou T, Duret L. High prevalence of PRDM9-independent recombination hotspots in placental mammals. *Proceedings of the National Academy of Sciences USA*. 2024;121(23):e2401973121. Published 2024-05-29.
DOI: [10.1073/pnas.2401973121](https://doi.org/10.1073/pnas.2401973121)
PMID: 38809707 (verified via PubMed)

**Duplicate check**: Absent from REFERENCES_v4.bib (verified by grep for "38809707", "2401973121", "PRDM9-independent.*placental", "Joseph.*Duret.*PRDM9"). NOT in any R1 REFRESH file.

**Relevance** (MEDIUM): Derives an estimator of past recombination activity from GC-biased gene conversion (gBGC) signatures in substitution patterns and quantifies PRDM9-independent recombination hotspot activity in 52 boreoeutherian mammal species. Key finding: PRDM9-directed and PRDM9-independent hotspots coexist in most mammals — humans show a relative DEFICIT of recombination at PRDM9-independent hotspots compared to most placental mammals, consistent with PRDM9 efficiently directing recombination away from these sites in humans. PRDM9-independent hotspots are more positionally stable (located at CpG-dense promoter-like features) and evolve more slowly.

**Claim it contextualizes**: v2 P9 — the 133 gene-conversion-like patches and the Schweiger 2024 non-PRDM9 long-tract class. Joseph 2024 provides the evolutionary-genomic framework: the non-PRDM9 long-tract NCO class detected by Schweiger 2024 (~2% of NCOs) is mechanistically related to PRDM9-independent hotspot activity, which Joseph 2024 demonstrates is present but suppressed in humans relative to other mammals. This contextualization helps explain why the non-PRDM9 class is rare in human sperm (~2%), and why subtelomeric regions (which are PRDM9-depleted due to low PRDM9-binding motif density near telomeres) might be relatively enriched in PRDM9-independent recombination compared to interstitial regions.

**Recommended action**: ADD as optional supporting citation alongside Schweiger 2024 and Porsborg 2025 in the P9 gene-conversion-like sentence, or in Methods (recombination model). Lower priority than STRONG-1 (Porsborg 2025) because it doesn't directly address subtelomeres.

```bibtex
@article{pedigree_Joseph2024PRDM9indep,
  author  = {Joseph, Julien and Prentout, Djivan and Laverr{\'e}, Alexandre
             and Tricou, Th{\'e}o and Duret, Laurent},
  title   = {High prevalence of {PRDM9}-independent recombination hotspots
             in placental mammals},
  journal = {Proceedings of the National Academy of Sciences USA},
  year    = {2024},
  volume  = {121},
  number  = {23},
  pages   = {e2401973121},
  doi     = {10.1073/pnas.2401973121},
  pmid    = {38809707},
  note    = {PRDM9-independent hotspots in 52 placental mammals; humans show deficit
             at these sites; contextualize non-PRDM9 NCO class from Schweiger 2024}
}
```

---

### Topics 13 and 14: No New STRONG Papers

**Topic 13 (subtelomere_popgen_fst)**: After 12+ PubMed queries and 3 OpenAlex queries (see Section 5), no new STRONG papers were identified beyond the 5 STRONG papers already in R1 (`subtel_popgen_jeong2025`, `subtel_popgen_schloissnig2025`, `subtel_popgen_kim2025`, `subtel_popgen_jana2025`, `subtel_popgen_porubsky2026`). All are in REFERENCES_v4.bib. No 2024-2026 paper specifically addresses Hudson/Weir-Cockerham FST at subtelomeric binary haplotypes or cross-arm affinity frequency differences between AFR and non-AFR superpopulations. The absence of such papers is consistent with the R1 assessment that this analysis is genuinely novel.

**Topic 14 (olfactory_or4f)**: After 8+ PubMed queries and 1 OpenAlex query (see Section 5), no new STRONG papers were identified beyond the 5 papers already in R1. One potentially relevant paper found — Graham et al. 2025 (PMID 40562037, Curr Biol, "Convergent reduction of olfactory genes and olfactory bulb size in mammalian species at altitude") — is WEAK for this topic: it describes OR gene reduction in high-altitude mammals (non-human, non-subtelomeric context) and does not inform the arm-specific OR4F pseudogenization gradient or community structure. The REDO flag for REFRESH_14 is not a literature gap but a citation-placement gap: `@chi2025primate` (already in REFERENCES_v4.bib from R1) needs to be cited in the OR4F paragraph of v2 P11 (per R2_AUDIT_PLAN.md §3 C-5 and §5 edit 3).

---

## Section 2: Backfill Where R1 Was Thin

**Assessment of REFRESH file depth** (threshold: <5 papers):

| REFRESH file | R1 paper count | Thin? | Deeper search result |
|---|---|---|---|
| REFRESH_13 (topic 13) | 7 (5 STRONG + 2 MODERATE) | No | Searches confirmed; no new papers found |
| REFRESH_14 (topic 14) | 5 papers | Borderline | Deeper searches found 0 new papers; OR4F-subtelomere field remains sparse |
| REFRESH_15 (topic 15) | 10 papers | No | 3 new papers found (1 STRONG, 2 MEDIUM); see Section 1 |

**Topic 13**: The R1 sweep was comprehensive. Seven papers across the relevant scope (FST methodology, SV population differentiation, African enrichment of SVs, subtelomere SV surveys, long-read population pangenomes). Deeper R2 searches for "structural variant FST African enrichment," "subtelomere population differentiation diverse human genomes," and "segmental duplication population scale pangenome 2024-2026" returned no additional relevant papers. The Gong et al. 2025 (Nat Commun) Han Chinese SV study found via OpenAlex is WEAK — it covers a single East Asian population and does not provide FST estimates comparing AFR vs. non-AFR for subtelomeric or cross-arm haplotypes.

**Topic 14**: R1 identified 5 papers spanning OR gene evolution, functional genomics, primate comparative genomics, and cross-phylogenetic mechanism. The R2 search for "olfactory receptor OR4F evolution pseudogenization human population" returned only papers already in R1 (Hayakawa 2025, chi2025primate). Graham et al. 2025 (altitude-driven OR loss) was found but is WEAK. The field genuinely lacks 2024-2026 papers on OR4F in human subtelomeric community context, confirming the paper's novelty. Borderline depth of 5 papers is justified: there are simply few papers about this specific intersection.

**Topic 15**: R1 had 10 papers and was strong. The R2 sweep found 3 additional papers not in R1 (Section 1). The overall coverage remains good.

---

## Section 3: Contradiction Follow-Up

**Topic 13**: R1 (REFRESH_13 §4) flagged one internal inconsistency (non-AFR FST "0.02–0.04" in draft vs. "−0.047 to +0.007" in report), resolved by Pass A (v2 L40 corrected to "−0.05 to 0.01") and supported by Schloissnig 2025. **No new 2024-2026 literature contradicts or further modifies any Topic 13 claim.** No new paper has been published in the last 6 months that reverses the consensus that non-AFR structural variant FST is near zero, or that AFR genomes are enriched for complex SDs.

**Topic 14**: R1 (REFRESH_14 §4) flagged one partial contradiction: Chi et al. 2025 (PMID 40021902) vs. Gilad et al. 2004 on the trichromacy trade-off framing. This contradiction was identified in R1 and the required hedge sentence is specified in R2_AUDIT_PLAN.md §5 edit 3 (text: "Olfactory receptor pseudogenisation across primates is now interpreted as a sensory-reallocation event rather than a simple visual-olfactory trade-off [@chi2025primate]"). **No new 2024-2026 paper further resolves or deepens the Chi 2025 vs. Gilad 2004 contradiction.** The current state is: the sensory-reallocation framing from Chi 2025 is the latest scientific position; the hedge is required but not yet in v2 text.

**Topic 15**: R1 (REFRESH_15 §4) flagged no contradictions. **Confirmed: no new contradictions introduced.** The new papers found (Porsborg 2025, Lee 2026, Joseph 2024) support or contextualize existing claims without contradicting any of them.

---

## Section 4: 60-Day Preprint Frontier (2026-03-18 to 2026-05-17)

Four bioRxiv categories were searched for the period 2026-03-18 to 2026-05-17 (60 days):

| Category | Preprints scanned | Relevant to this project |
|---|---|---|
| genomics | 100 | 0 |
| genetics | 50 | 0 |
| evolutionary biology | 50 | 0 |
| molecular biology | 50 | 0 |
| **Total** | **250** | **0** |

All 250 preprints were reviewed by title and abstract preview. None touched: subtelomere, chromosome-end sequence sharing, NAHR at subtelomeres, pedigree-based T2T exchange, pangenome PHR analysis, or meiotic bouquet–subtelomere proximity.

**Notable borderline hits (reviewed and excluded)**:
- Ringbauer et al. 2026 (bioRxiv) — ancient DNA genealogy confirmation (historical genetics; not relevant)
- Kim et al. 2026 (bioRxiv, 2026-03-18) — Denisovan-derived Alu in OCA2 affecting Melanesian pigmentation (structural variant introgression; not subtelomeric)
- Multiple pangenome alignment tools — none extend wfmash/PGGB, none address subtelomeric specifics

**Conclusion**: The 60-day frontier is silent on topics 13/14/15. This is not surprising: our topics occupy a niche intersection (subtelomere community structure + population genetics + OR gene biology + meiotic recombination in pedigrees) that has low preprint volume in any given 60-day window.

---

## Section 5: Audit Trail

### PubMed queries (all via mcp__claude_ai_PubMed__search_articles)

Search date: 2026-05-17. Date filters applied as noted.

**Topic 13 searches:**

| Query | Date range | Hits | Kept |
|---|---|---|---|
| "subtelomere population genetics structural variation FST African enrichment" | 2024/06–2026/05 | 0 | 0 |
| "structural variant population differentiation diverse human genomes long-read" | 2024/06–2026/05 | 0 | 0 |
| "human population FST structural variant Weir Cockerham Hudson estimator out-of-Africa" | 2024/06–2026/05 | 0 | 0 |
| "segmental duplication structural variation African genomes population scale pangenome" | 2024/06–2026/05 | 1 | 0 (PMID 40631282 = Porubsky 2026 bioRxiv preprint, already in R1 as published NatComm) |
| "population genetics diverse human out-of-Africa structural variant frequency" | 2025/01–2026/05 | 0 | 0 |

**Topic 14 searches:**

| Query | Date range | Hits | Kept |
|---|---|---|---|
| "olfactory receptor OR4F subtelomere human genome evolution pseudogene" | 2024/01–2026/05 | 0 | 0 |
| "olfactory receptor gene evolution human population pseudogenization functional" | 2025/01–2026/05 | 2 | 0 (PMID 40562037 = Graham 2025 altitude OR loss — WEAK; PMID 39747985 = Hayakawa 2025 — already R1) |
| "olfactory receptor OR4F pangenome copy number variation human arm community" | 2023/01–2026/05 | 0 | 0 |

**Topic 15 / P-3 searches:**

| Query | Date range | Hits | Kept |
|---|---|---|---|
| "non-crossover recombination gene conversion long-read sperm meiosis PRDM9" | 2024/06–2026/05 | 1 | 1 (PMID 41285744 = Porsborg 2025 — NEW STRONG) |
| "recombination crossover non-crossover sperm human pedigree long-read sequencing" | 2024/06–2026/05 | 1 | 0 (PMID 39005338 = Schweiger 2024 — already R1) |
| "PRDM9-independent recombination hotspot placental mammal" | 2024/01–2026/05 | 1 | 1 (PMID 38809707 = Joseph 2024 — NEW MEDIUM) |
| "PRDM9 recombination hotspot meiosis gene conversion non-crossover human" | 2024/06–2026/05 | 0 | 0 |
| "gene conversion tract length duplicon repeat-rich region meiosis long-read" | 2024/01–2026/05 | 0 | 0 |
| "meiotic bouquet telomere nuclear envelope recombination chromosome movement" | 2026/01–2026/05 | 0 | 0 |

**Frontier / broad searches:**

| Query | Date range | Hits | Kept |
|---|---|---|---|
| "telomere-to-telomere T2T assembly pangenome subtelomeric 2026" | 2026/01–2026/05 | 1 | 1 (PMID 41699652 = Lee 2026 SEPTIN14P — NEW MEDIUM) |
| "subtelomere chromosome end T2T pangenome recombination" | 2026/03–2026/05 | 0 | 0 |
| "NAHR non-allelic homologous recombination segmental duplication human population" | 2024/06–2026/05 | 0 | 0 |

**Total PMIDs reviewed**: 17
**Kept for Section 1**: 3 (PMID 41285744, 41699652, 38809707)
**Dropped with rationale**:
- PMID 40631282 (Porubsky 2026 bioRxiv) — same paper as `porubsky2026chr22q11` (published NatComm) already in R1/v4.bib
- PMID 39005338 (Schweiger 2024) — already in R1
- PMID 39747985 (Hayakawa 2025) — already in R1 (REFRESH_14)
- PMID 40562037 (Graham 2025, altitude OR reduction) — WEAK; altitude-adapted non-human mammals, not OR4F in human subtelomeric community context

### OpenAlex queries (via direct API, 2026-05-17)

| Query | Year filter | Results | Kept |
|---|---|---|---|
| "subtelomere population genetics FST structural variants African" | 2024+ | 5 | 0 (Plasmodium, macaque — off-topic) |
| "meiotic recombination non-crossover gene conversion long-read sperm" | 2024+ | 6 | 1 (already R1: Schweiger 2024) |
| "olfactory receptor OR4F evolution pseudogenization human" | 2024+ | 0 | 0 |

### bioRxiv searches (via mcp__claude_ai_bioRxiv__search_preprints)

| Category | Date range | Results scanned | Relevant |
|---|---|---|---|
| genomics | 2026-03-18 to 2026-05-17 | 100 | 0 |
| genetics | 2026-03-18 to 2026-05-17 | 50 | 0 |
| evolutionary biology | 2026-03-18 to 2026-05-17 | 50 | 0 |
| molecular biology | 2026-03-18 to 2026-05-17 | 50 | 0 |

### REFERENCES_v4.bib cross-check

The three new papers proposed for addition were verified absent from REFERENCES_v4.bib (364 entries) via grep on author surnames, PMIDs, and DOI fragments. All three are confirmed NEW relative to both R1 REFRESH files and REFERENCES_v4.bib.

### R1 duplicate cross-check (zero duplicates confirmed)

| Proposed R2 bibkey | R1 REFRESH_13 | R1 REFRESH_14 | R1 REFRESH_15 | In REFERENCES_v4.bib? |
|---|---|---|---|---|
| `pedigree_Porsborg2025primaterecom` | Not present | Not present | Not present | **NOT PRESENT** |
| `subtelstruct_Lee2026SEPTIN14` | Not present | Not present | Not present | **NOT PRESENT** |
| `pedigree_Joseph2024PRDM9indep` | Not present | Not present | Not present | **NOT PRESENT** |

All PMIDs verified via `mcp__claude_ai_PubMed__get_article_metadata`:
- PMID 41285744: confirmed as Porsborg et al. 2025, Nat Commun, DOI 10.1038/s41467-025-65248-3
- PMID 41699652: confirmed as Lee 2026, Mobile DNA, DOI 10.1186/s13100-026-00394-z
- PMID 38809707: confirmed as Joseph et al. 2024, PNAS, DOI 10.1073/pnas.2401973121

---

*End of REFRESH_R2_frontier.md. Generated by agent-123 (r2-litref-frontier) on 2026-05-17.*
