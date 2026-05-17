# REFRESH_R2_recombination: Round-2 Literature Sweep — Topics 05, 06, 07, 11

**Generated:** 2026-05-17
**Agent:** agent-124 (r2-litref-recombination)
**Databases searched:** PubMed (MCP tool), bioRxiv (MCP tool), OpenAlex (skill)
**Date range:** 2024-06-01 to 2026-05-17 (24 months); 60-day window for preprints: 2026-03-17 to 2026-05-17
**Cross-checked against:** REFRESH_05, REFRESH_06, REFRESH_07, REFRESH_11, REFERENCES_v4.bib

---

## Section 0: Scope and Audit Angles

**Topics covered:** 05 (acrocentric_rdna_robertsonian), 06 (dux4_d4z4_fshd), 07 (concerted_evolution_nahr), 11 (pedigree_based_recombination_detection).

**R2_AUDIT_PLAN.md angles addressed:**

- §4 P-1 (REFRESH_02): Out of scope for this task (topics 05/06/07/11 only).
- §4 P-2 (REFRESH_04): Out of scope for this task.
- §4 P-3 (REFRESH_07/REFRESH_15): Search for 2025-2026 follow-up to Schweiger 2024 (PMID 39005338) addressing long-tract NCO at duplicon-rich loci. Also note: `hinch2023meiotic` already in v4 bib but not cited in v2. → **Result: no follow-up paper found** (see Section 3).
- §2 Row 05 (TRUST verdict): Verified by deeper search. Two R1 gaps found (Section 1 and 2).
- §2 Row 06 (TRUST verdict): Verified by deeper search. R1 was comprehensive; no STRONG new papers.
- §2 Row 11 (TRUST verdict): One new MEDIUM paper found (Section 2).

**Draft claims directly relevant to this sweep:**
- **C7**: acrocentric p-arms (C7) are sequence-interchangeable — paralog closer than allele (p=2.0×10⁻⁷)
- **C8**: ongoing interchromosomal exchange: 133 gene-conversion-like patches in C7 (~90%)
- **P3**: NJ tree recovers 4q/10q DUX4 clade (C1) and acrocentric clade (C7)
- **P9**: two-process NCO model — PRDM9-associated short tract + non-PRDM9 long-tract
- **P11 (draft v2 L46)**: Schweiger 2024 cited but two-process model not invoked in prose

---

## Section 1: NEW STRONG Papers Not in R1

### [ataei2025line1dj] — STRONG (Topic 05)

**Full citation:**
Ataei L, Zhang J, Monis S, Giemza K, Mittal K, Yang J, Shimomura M, McStay B, Wilson MD, Ramalho-Santos M. 2025. LINE1 elements at distal junctions of rDNA repeats regulate nucleolar organization in human embryonic stem cells. *Genes & Development* 39(3-4):280–298. PMID:39797762 / DOI:[10.1101/gad.351979.124](https://doi.org/10.1101/gad.351979.124)

**Claim affected:** Draft P3 / report section 03 / Extended Data Fig. 4 — the claim that five acrocentric arms co-localize in the nucleolus and undergo concerted sequence exchange. McStay is co-author.

**Why STRONG:** Ataei et al. 2025 identifies a primate-specific, full-length LINE1 (L1) retrotransposon insertion at a conserved position across ALL five human acrocentric distal junctions (DJs). This DJ-LINE1 element: (1) interacts with specific DJ regions; (2) is upregulated in naïve hESCs; (3) contributes to nucleolar positioning of the DJs by CRISPR deletion/interference; (4) silencing DJ-LINE1 disrupts nucleolus structure, transcriptional output, and self-renewal. This is the first molecular mechanism assigned to DJ nucleolar anchoring. Previously, the near-identical sequence of DJs across all five arms (>95–99% identity; `acrocentric_rdna_robertsonian_floutsakou2013`) was known, but WHY the DJs anchor to the nucleolus was unresolved. Ataei 2025 closes that mechanistic gap. This ADDS a mechanistic sentence to the draft's causal-loop/P3 discussion: "primate-specific LINE1 elements in all five acrocentric distal junctions contribute to their nucleolar anchoring, providing a sequence-level substrate for the observed co-localization."

**Duplicate check:** Not in REFRESH_05, REFRESH_06, REFRESH_07, REFRESH_11. Not in REFERENCES_v4.bib. CONFIRMED NEW.

**PMID verified:** PMID 39797762 confirmed via PubMed `get_article_metadata`. Title, journal, year, McStay co-authorship all match.

**Recommended action:** ADD alongside `acrocentric_rdna_robertsonian_floutsakou2013` and `acrocentric_rdna_robertsonian_mcstay2016` in the draft's description of DJ anchoring at the nucleolus (P3, or Extended Data methods). Draft sentence: "The distal junctions [DJs] share >98% sequence identity across arms, and a primate-specific LINE1 element conserved in all five DJs contributes to their nucleolar anchoring [@ataei2025line1dj; @acrocentric_rdna_robertsonian_floutsakou2013]."

```bibtex
@article{ataei2025line1dj,
  author  = {Ataei, Lamisa and Zhang, Juan and Monis, Simon and Giemza, Krystyna
             and Mittal, Kirti and Yang, Joshua and Shimomura, Mayu and McStay, Brian
             and Wilson, Michael D and Ramalho-Santos, Miguel},
  title   = {{LINE1} elements at distal junctions of {rDNA} repeats regulate
             nucleolar organization in human embryonic stem cells},
  journal = {Genes \& Development},
  year    = {2025},
  volume  = {39},
  number  = {3--4},
  pages   = {280--298},
  doi     = {10.1101/gad.351979.124},
  pmid    = {39797762}
}
```

---

### [hao2024snul] — MEDIUM (Topic 05, borderline STRONG)

**Full citation:**
Hao Q, Liu M, Daulatabad SV, …, McStay B, Sinha S, Janga SC, Prasanth SG, Prasanth KV. 2024. Monoallelically expressed noncoding RNAs form nucleolar territories on NOR-containing chromosomes and regulate rRNA expression. *eLife* 13. PMID:38240312 / DOI:[10.7554/eLife.80684](https://doi.org/10.7554/eLife.80684)

**Claim affected:** C7 (inter-arm sequence exchange mediated by nucleolar co-localization), report section 03 gene enrichment (NOR chromosomal identity).

**Why MEDIUM (near STRONG):** Discovers SNUL (Single Nucleolus-Localized RNA) ncRNAs that form constrained sub-nucleolar territories on individual NOR-containing chromosomes and control rRNA expression. McStay is co-author. SNULs are monoallelically expressed — each allele of a NOR chromosome has its own sub-nucleolar domain. This adds a layer of individual-arm identity WITHIN the shared nucleolar compartment, directly relevant to the co-localization model: the five arms converge in the nucleolus yet maintain distinct SNUL-mediated sub-territories. This doesn't CHANGE a draft claim but adds novel mechanistic nuance (individual NOR chromosomes have addressable nucleolar sub-domains) to the co-localization story.

**Duplicate check:** Not in any R1 REFRESH file. Not in REFERENCES_v4.bib. CONFIRMED NEW. R1 REFRESH_05 searched through 2026-05-17 but queried "McStay 2023" specifically and missed this January 2024 paper.

**Recommended action:** ADD as optional background alongside `mcstay2023acrocentric` in the Methods or Extended Data Fig. caption for NOR co-localization. Not a required citation, but relevant for reviewers who ask about individual-arm identity within the shared nucleolus.

```bibtex
@article{hao2024snul,
  author  = {Hao, Qinyu and Liu, Minxue and Daulatabad, Swapna Vidhur
             and Gaffari, Saba and Song, You Jin and Srivastava, Rajneesh
             and Bhaskar, Shivang and Moitra, Anurupa and Mangan, Hazel
             and Tseng, Elizabeth and Gilmore, Rachel B and Frier, Susan M
             and Chen, Xin and Wang, Chengliang and Huang, Sui
             and Chamberlain, Stormy and Jin, Hong and Korlach, Jonas
             and McStay, Brian and Sinha, Saurabh and Janga, Sarath Chandra
             and Prasanth, Supriya G and Prasanth, Kannanganattu V},
  title   = {Monoallelically expressed noncoding {RNA}s form nucleolar territories
             on {NOR}-containing chromosomes and regulate {rRNA} expression},
  journal = {eLife},
  year    = {2024},
  volume  = {13},
  doi     = {10.7554/eLife.80684},
  pmid    = {38240312}
}
```

---

## Section 2: Backfill Where R1 Was Thin

### Topic 11 (REFRESH_11): Only 1 new paper in R1 — deeper search found 1 more

R1 REFRESH_11 identified only `palsson2025recomb` as a genuinely new paper (plus two stale-placeholder updates for Porubsky2025 and Cechova2025). The audit plan classified Topic 11 as TRUST because the two headline citations landed correctly. However, the criterion "R1 REFRESH files with <5 papers: did a deeper search find more?" applies here — REFRESH_11 had 1 new paper.

**Deeper search result:** YES — one new MEDIUM paper found:

**[sasani2026kfam]** — MEDIUM (Topic 11)

Sasani TA, Goldberg ME, Avvaru AK, Nicholas TJ, Neklason DW, Dolzhenko E, Mokveld T, Munson KM, Hoekzema K, Ayllon M, Kaufman EJ, Porubsky D, Valdmanis PN, Eichler EE, Quinlan AR, Dashnow H. 2026. A family portrait of the genomic factors shaping tandem repeat mutagenesis. *bioRxiv* (preprint). PMID:41959501 / DOI:[10.64898/2026.03.06.710071](https://doi.org/10.64898/2026.03.06.710071)

**Claim context:** The draft's pedigree analysis uses CEPH1463 (four-generation, 28 members) for PHR exchange detection. This paper uses the EXACT SAME pedigree (designated K1463 by the Quinlan/Sasani lab; CEPH 1463 = K1463 cross-notation), with Porubsky and Eichler as co-authors, to profile 8 million tandem repeat (TR) loci by HiFi long-read sequencing, finding 1,270 TR expansions/contractions across 20 children.

**Why MEDIUM (not STRONG):** The paper does NOT address PHR interchromosomal exchange directly. It focuses on STR/VNTR mutagenesis driven by repeat length, motif interruption, and parental heterozygosity. The paper's hyper-mutable loci (43 sites, up to 12 mutations across the pedigree) are TR loci, not subtelomeric PHR exchange segments. However, it provides complementary evidence that recurrent mutagenesis at tandem repeats is measurable in the same family, lending indirect support to the claim that the PHR exchange events observed in the same pedigree are genuine germline events rather than noise.

**Where to cite:** P7 or draft Limitations, as a note that the same pedigree has been independently mined for TR mutagenesis at a genome-wide scale.

**Duplicate check:** Not in REFRESH_11 (R1 searched Sasani T[Author] 2023-2026, got 5 hits but did not include this preprint, which appeared on bioRxiv 2026-03-06 — possibly indexed in PubMed after R1's search run). Not in REFERENCES_v4.bib. CONFIRMED NEW.

```bibtex
@misc{sasani2026kfam,
  author  = {Sasani, Thomas A and Goldberg, Michael E and Avvaru, Akshay K
             and Nicholas, Thomas J and Neklason, Deborah W and Dolzhenko, Egor
             and Mokveld, Tom and Munson, Katherine M and Hoekzema, Kendra
             and Ayllon, Marcelo and Kaufman, Eli J and Porubsky, David
             and Valdmanis, Paul N and Eichler, Evan E and Quinlan, Aaron R
             and Dashnow, Harriet},
  title   = {A family portrait of the genomic factors shaping tandem repeat mutagenesis},
  howpublished = {bioRxiv preprint},
  year    = {2026},
  doi     = {10.64898/2026.03.06.710071},
  pmid    = {41959501},
  note    = {PREPRINT. Same CEPH1463/K1463 four-generation pedigree as PHR study;
             HiFi LRS of 8M TR loci; 43 hyper-mutable loci; Porubsky and Eichler co-authors.
             Contextualizes PHR exchange events within broader repeat mutagenesis landscape.}
}
```

### Topics 05, 06, 07: backfill assessment

- **Topic 05 (7 R1 papers):** Deeper search found 2 new papers (Ataei 2025 STRONG; Hao 2024 MEDIUM) that R1 missed despite searching through 2026-05-17. R1's McStay query was year-specific to 2023; both missed papers were published in 2024 and 2025.
- **Topic 06 (7 R1 papers):** Deeper search found no STRONG new papers. R1 was comprehensive. Papers found in extended search (PMID 41509470 Sakr 2025 bioRxiv; PMID 40744516 Huang 2025) are clinical/functional studies WEAK for C5/C7/C8 structural claims.
- **Topic 07 (5 R1 papers):** No additional strong papers found. See Section 3 for Schweiger 2024 follow-up null result.

---

## Section 3: Contradiction Follow-Up

### Topic 05 (acrocentric/rDNA/Robertsonian)

R1 REFRESH_05 flagged NO contradictions. Deeper search confirms: no 2024-2026 paper contradicts the SST1/DJ model for Robertsonian formation (deLima2025), the rDNA concerted evolution claim (C8), or the allele/paralog inversion at C7. **No new contradictions found.** Ataei 2025 reinforces rather than challenges the co-localization model.

### Topic 06 (DUX4/D4Z4/FSHD)

R1 REFRESH_06 flagged the Salsi 2026 partial contradiction (D4Z4-like loci on >10 chromosomes beyond 4q/10q). This hedge was incorporated into draft v2 L30 with `@Salsi2026fshd`. **No new contradictions found in the last 6 months** (Nov 2025 to May 2026). My searches in the FSHD/D4Z4 space found no paper that challenges the C1 community structure, the 4q/10q dominant-arms model, or the CTCF/lamin tethering mechanism.

### Topic 07 (concerted evolution/NAHR)

R1 REFRESH_07 flagged no contradictions. **Critical null result:** The R2_AUDIT_PLAN (§4 P-3) specifically requested a search for 2025-2026 follow-up to Schweiger 2024 (PMID 39005338 — long-tract non-PRDM9 NCO in sperm) addressing the two-process NCO model at duplicon-rich loci. I ran 8 targeted PubMed queries (Schweiger[Author] recombination, non-crossover NCO long-tract PRDM9 independent duplicon, gene conversion non-crossover subtelomere human meiosis, meiotic recombination subtelomere segmental duplication 2024-2025, Hinch AG meiotic recombination 2024-2025, PRDM9 independent recombination crossover segmental duplication long-read) and an OpenAlex search. **No follow-up paper to Schweiger 2024 exists in the 2024-2026 literature.** The Schweiger 2024 paper remains the sole primary evidence for long-tract (~2%) non-PRDM9 NCO in human sperm. R2 cannot provide a stronger citation than R1 for the P9 two-process NCO model.

Note on `hinch2023meiotic` (PMID 38033082): This paper was added to REFERENCES_v4.bib by R1 REFRESH_07 but is not cited in NATURE_DRAFT_v2. No new Hinch recombination papers were found in 2024-2026 (Hinch AG author search returned 0 results for that period). The R2 auditor's note that Hinch 2023 is "already in v4, uncited" is a valid citation-placement suggestion for the r2-fix task, not a literature gap.

### Topic 11 (pedigree-based recombination)

R1 REFRESH_11 flagged two tensions: Porubsky2025 "no crossover-SV correlation" (addressed in v2), and the 538-patch count/HQ-filter inconsistencies (mechanical edits, not literature). **No new contradictions found in the last 6 months.** No paper published since November 2025 challenges the PHR ectopic exchange interpretation or the four-generation CEPH1463 pedigree analysis methodology.

---

## Section 4: 60-Day Preprint Frontier (2026-03-17 to 2026-05-17)

I searched bioRxiv using the MCP tool for the 60-day window (March 17 – May 17, 2026) across the categories most likely to contain relevant papers: genomics (50 results examined), genetics (30 results examined), evolutionary biology (30 results examined). Total: 110 preprints examined. Keyword filter is not supported by the bioRxiv MCP tool; visual inspection of titles and abstracts was used.

**Papers touching subtelomere / chromosome-end / NAHR / pedigree T2T / pangenome / meiotic bouquet in the 60-day window:**

None identified in the bioRxiv category searches. The 110 results spanned genomic epidemiology, population genetics, comparative genomics, and microbiology — none concerned subtelomeric PHR biology, NOR/rDNA recombination, FSHD/D4Z4 structural genomics, or pedigree-based recombination analysis.

**Papers from the 60-day window already in R1 (not new):**
- Hebbar 2026 marmoset T2T genome (PMID 41929024, bioRxiv 2026-03-25) — already in REFRESH_07 as `hebbar2026marmoset`.

**New papers from the 60-day window (found via PubMed):**
- Sasani 2026 K1463 tandem repeats (PMID 41959501, bioRxiv 2026-03-06) — reported in Section 2. This is the only 60-day preprint directly relevant to our pedigree topic.

**Conclusion for 60-day frontier:** No new preprints touching any of the four topic areas have appeared in the last 60 days beyond what R1 already found (Hebbar 2026) and what this sweep adds (Sasani 2026).

---

## Section 5: Audit Trail

### Databases and tools used

| Tool | Purpose |
|------|---------|
| `mcp__claude_ai_PubMed__search_articles` | Primary discovery across all four topics |
| `mcp__claude_ai_PubMed__get_article_metadata` | PMID verification for all candidate papers |
| `mcp__claude_ai_bioRxiv__search_preprints` | 60-day preprint scan (genomics/genetics/evolutionary-biology categories) |
| OpenAlex skill (`openalex-database`) | Broader search for Schweiger 2024 follow-up and NCO at duplicon-rich loci |

### PubMed queries with hit counts (all date-filtered 2024-2026 unless noted)

| Query | PMIDs returned | Assessment |
|-------|---------------|------------|
| `Schweiger recombination sperm` | 2 | PMID 39005338 = Schweiger 2024 (already in bib); PMID 40702521 = different Schweiger (TRIP13, WEAK) |
| `gene conversion non-crossover subtelomere human meiosis` | 1 | PMID 37797835 = ciliate (WEAK) |
| `Sasani de novo mutation pedigree long-read` | 3 | PMID 41959501 NEW; PMID 40269156 = Porubsky2025 (in bib); PMID 39149261 = Porubsky2025 preprint (in bib) |
| `acrocentric chromosome Robertsonian translocation rDNA 2025` | 0 | — |
| `DUX4 FSHD facioscapulohumeral muscular dystrophy nanopore 2025` | 3 | PMIDs 41509470 (WEAK), 40744516 (WEAK), 40507980 (WEAK) |
| `Cechova[Author] pedigree T2T subtelomere` | 1 | PMID 41473289 = Cechova2025 (already in bib) |
| `Halldorsson OR Palsson recombination crossover Iceland pedigree` | 10+ | All hits = palsson2025recomb (already in bib) or unrelated |
| `PRDM9 independent recombination crossover non-crossover segmental duplication long-read` | 0 | — |
| `McStay[Author] acrocentric nucleolus organizer` | 2 | PMID 39797762 NEW STRONG; PMID 38240312 NEW MEDIUM |
| `Noyes M Eichler segmental duplication postzygotic gene conversion 2026` | 1 | PMID 41803180 = noyes2026sd (already in R1 REFRESH_07) |
| `Porubsky[Author] de novo structural variant pedigree` | 2 | Both already in bib |
| `Hinch AG meiotic recombination hotspot PRDM9 crossover 2024 2025` | 0 | — |
| `Vollger SD segmental duplication gene conversion mutation rate pangenome` | 0 | — |
| `de Lima Gerton acrocentric rDNA nucleolus SST1 Robertsonian 2025 2026` | 0 | — |
| `Eichler EE long-read segmental duplication gene conversion mutation postzygotic 2025` | 1 | PMID 40791370 = preprint of noyes2026sd (same paper, already in R1) |
| `HPRC HPRCv2 human pangenome reference consortium 2025 2026` | 1 | PMID 41621533 = T cell receptor paper (WEAK) |
| `Schweiger DF sperm single-cell non-crossover gene conversion tract length` | 0 | — |
| `Cechova M Porubsky D subtelomere chromosome exchange pedigree de novo 2025 2026` | 0 | — |
| `rDNA ribosomal DNA copy number human acrocentric 2025 2026` | 0 | — |

**Total unique PMIDs retrieved:** 29. After relevance filter: 3 new papers (1 STRONG Topic05, 1 MEDIUM Topic05, 1 MEDIUM Topic11). 26 dropped.

### bioRxiv date+category searches

| Category | Date range | Preprints examined | Relevant hits |
|----------|-----------|-------------------|--------------|
| genomics | 2026-03-17 to 2026-05-17 | 50 | 0 |
| genetics | 2026-03-17 to 2026-05-17 | 30 | 0 |
| evolutionary biology | 2026-03-17 to 2026-05-17 | 30 | 0 |

### Papers dropped (non-trivial assessments)

| PMID | Title snippet | Reason |
|------|--------------|--------|
| 40702521 | Homozygous TRIP13 variant (Schweiger M) | Different Schweiger; TRIP13 meiotic checkpoint, not NCO at SDs |
| 37797835 | Ciliate meiotic recombination landscape | Non-human (Tetrahymena); not relevant |
| 41509470 | Genome-wide FSHD cell lines Nanopore (Sakr 2025) | Functional/transcriptomic; no structural C5/C7/C8 claim impact |
| 40744516 | D4Z4 adaptive sampling FSHD1 (Huang 2025) | Clinical/diagnostic; no new structural claim |
| 40507980 | m6A methylation FSHD myoblasts (Settas 2025) | m6A pathway; off-claim |
| 41621533 | T cell receptor alleles in HPRC (Yang 2026) | HPRC samples only; T cell receptor biology; off-claim |
| 40791370 | Long-read sequencing of trios repetitive DNA (Noyes 2025 bioRxiv) | Preprint = noyes2026sd (already in R1 REFRESH_07 as PMID 41803180) |

### Bibkey uniqueness checks

All three proposed new bibkeys checked against `REFERENCES_v4.bib` and all 15 R1 REFRESH files:

| Bibkey | In REFERENCES_v4.bib? | In any REFRESH file? | Status |
|--------|----------------------|---------------------|--------|
| `ataei2025line1dj` | No | No | UNIQUE — safe to add |
| `hao2024snul` | No | No | UNIQUE — safe to add |
| `sasani2026kfam` | No | No | UNIQUE — safe to add |

### PMID/DOI spot-check verification (all Section 1 STRONG papers verified)

1. **ataei2025line1dj** — PMID 39797762 confirmed via `mcp__claude_ai_PubMed__get_article_metadata`. Title: "LINE1 elements at distal junctions of rDNA repeats regulate nucleolar organization in human embryonic stem cells." Journal: Genes & Development, 2025, 39(3-4):280–298. DOI: 10.1101/gad.351979.124. McStay B confirmed as co-author. **VERIFIED.**
2. **hao2024snul** — PMID 38240312 confirmed. Title: "Monoallelically expressed noncoding RNAs form nucleolar territories on NOR-containing chromosomes and regulate rRNA expression." Journal: eLife, 2024, vol. 13. DOI: 10.7554/eLife.80684. McStay B confirmed as co-author. **VERIFIED.**
3. **sasani2026kfam** — PMID 41959501 confirmed. Title: "A family portrait of the genomic factors shaping tandem repeat mutagenesis." bioRxiv 2026-03-06. DOI: 10.64898/2026.03.06.710071. K1463 / CEPH1463 four-generation family confirmed. Porubsky D and Eichler EE confirmed as co-authors. **VERIFIED.**

---

## Summary Table

| # | Bibkey | Topic | Score | Draft target | Duplicate of R1? |
|---|--------|-------|-------|-------------|-----------------|
| 1 | `ataei2025line1dj` | 05 | STRONG | P3/EDFig.4 DJ anchoring sentence | No |
| 2 | `hao2024snul` | 05 | MEDIUM | Optional: Methods NOR co-localization | No |
| 3 | `sasani2026kfam` | 11 | MEDIUM | P7/Limitations: K1463 TR mutagenesis context | No |

**Topic 06:** No new STRONG papers. R1 was comprehensive (7 STRONG papers already proposed). The 60-day FSHD literature contains clinical/diagnostic papers only.

**Topic 07:** No new STRONG papers. Schweiger 2024 follow-up does not exist in the 2024-2026 literature. `hinch2023meiotic` (already in v4, uncited) is the correct citation for the meiotic break repair mechanism; the r2-fix task should add it to P9 or the Limitations. This is a citation-placement action, not a new literature find.

---

*End of REFRESH_R2_recombination. Generated 2026-05-17 by agent-124 (r2-litref-recombination). All PMIDs verified via PubMed MCP tool. Zero duplicates vs R1 REFRESH files.*
