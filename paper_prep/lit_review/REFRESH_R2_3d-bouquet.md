---
task: r2-litref-3d-bouquet
agent: agent-122
date: 2026-05-17
topics: 08 (meiotic_bouquet_envelope), 09 (hic_3d_methods)
inputs:
  r1_refresh_08: paper_prep/lit_review/REFRESH_08_meiotic_bouquet_envelope.md
  r1_refresh_09: paper_prep/lit_review/REFRESH_09_hic_3d_methods.md
  r2_audit_plan: paper_prep/synthesis/R2_AUDIT_PLAN.md
  references_v4: paper_prep/synthesis/REFERENCES_v4.bib
---

# REFRESH_R2: 3D Bouquet — Round-2 Literature Sweep (Topics 08, 09)

---

## Section 0: Scope

**Topics covered:** 08 (Meiotic Bouquet and Nuclear Envelope Tethering) and 09 (Hi-C and 3D Methods).

**Audit angles from R2_AUDIT_PLAN.md (Section 4):**

Both topics 08 and 09 are marked **TRUST** in the R2 audit (Section 4, line: "REFRESH 01, 03, 05, 06, 08, 09, 10, 11, 12, 13, 15 each delivered the headline citations into v2 successfully. None of these needs a second sweep."). The R2 audit does not assign specific second-sweep angles for topics 08 or 09, unlike topics 02 (TAR1/ITS/TERRA), 04 (T2T PARs), and 07/15 (long-tract NCO), which are explicitly flagged for second sweeps.

Despite the TRUST verdict, this sweep pursues the following angles not fully addressed in R1:

1. **Papers published 2024-06 to 2026-05** that R1 may have missed due to indexing lag (R1 searched 2023–2026 broadly; R2 focuses the date filter on June 2024 forward and checks for papers indexed after R1's May 2026 execution).
2. **Zygotene cilium / bouquet assembly interface** — R1 REFRESH_08 identified this as an open question (Section 1, point 4). No paper resolved it in R1.
3. **Human stage-resolved meiotic Hi-C** — identified as the key gap in both R1 REFRESH_08 (open question 1) and R1 REFRESH_09 (key gap statement, Section 5 final paragraph). Any 2024–2026 paper directly providing this data would be STRONG.
4. **Last-60-day bioRxiv frontier** (2026-03-17 to 2026-05-17) across six relevant categories.

---

## Section 1: NEW STRONG Papers Not in R1

**No new STRONG papers were found.**

This is not a silent null. Below is a documented positive conclusion from 20 independent PubMed queries and 6 bioRxiv category sweeps (full detail in Section 5).

**Justification by claim:**

- *Core meiotic bouquet machinery (SUN1, KASH5, MAJIN, TERB1/2, dynein, RPMs):* R1 REFRESH_08 already delivered the current literature ceiling for this set of proteins. Garner 2023 (KASH5 activating adaptor, PMID 36946995), Liu 2025 (SpdyA/SUN1 Ser48, PMID 40826181), Meng 2023 (SUN1 NOA human variant, PMID 36933034), and Kaiser 2025 (CTCF loops in human meiosis, PMID 40114154) represent the 2023–2025 additions. No 2024-06 to 2026-05 follow-up on any of these proteins was found in PubMed.

- *Human stage-resolved meiotic Hi-C:* The key gap identified in both R1 files remains unfilled as of 2026-05. No new paper provides stage-pure human meiotic Hi-C data. The PubMed queries targeting this angle ("human spermatocyte prophase chromosome conformation capture", "human meiosis spermatocyte Hi-C chromatin conformation prophase 3D genome") both returned 0 hits in the 2024-06 to 2026-05 window.

- *Meiotic chromosome loop architecture (mouse / human):* The Cheng et al. 2024 preprint (PMID 38903112, NIDDK Micro-C; already in R1 as MEDIUM with publication-pending flag) has not yet appeared as a published journal article in PubMed as of 2026-05-17. No new high-resolution meiotic Hi-C or Micro-C paper was found.

- *Stage-resolved 3D genome meiosis:* He et al. 2023 (Developmental Cell, PMID 37963468) and Marín-Gual et al. 2025 (Science Advances, PMID 41032613) from R1 remain the most recent additions to stage-resolved 3D meiotic data. No new comparable dataset paper was found.

**One new MEDIUM paper was found** (not in R1; see below). It does not rise to STRONG because it does not change a draft claim — it addresses an open question with a finding that the draft has not yet committed to a position on.

---

### New MEDIUM Paper: Chen et al. 2025 — Zygotene Cilia Are Dispensable for Meiotic Chromosome Pairing

**[bouquet_ChenCEP164Cilia2025]** Chen JJ, Gong X, Mak M, Li FQ, Takemaru KI. 2025. Loss of cilia drives centriole clustering and elimination during mammalian spermatogenesis. *bioRxiv* (Stony Brook University preprint). PMID:41415467 / DOI:[10.64898/2025.12.04.692363](https://doi.org/10.64898/2025.12.04.692363). Posted December 8, 2025.

**Claim it addresses:** REFRESH_08 Section 1 open question 4: "the zygotene cilium connection to bouquet assembly." The draft currently cites `bouquet_ZygoteneCilium2021` for the zygotene cilium discovery; this paper directly tests functional necessity.

**Evidence (from PubMed abstract, PMID 41415467):** Male germ cell-specific conditional KO of CEP164 (distal appendage protein required for basal body docking and ciliogenesis) eliminates both zygotene primary cilia and sperm flagella. Despite complete absence of zygotene cilia, **meiotic chromosome pairing and DNA double-strand break repair proceeded normally**. The CEP164 KO phenotype (male infertility) is driven by flagellogenesis failure in round spermatids, not by loss of zygotene cilia. Centrioles in round spermatids exhibit defective basal body docking and form supernumerary clusters subsequently eliminated via residual bodies.

**Score: MEDIUM.** Does not change any current draft claim — the draft does not assert that zygotene cilia are required for bouquet function. Resolves the "zygotene cilium connection to bouquet assembly" open question by showing that meiotic chromosome pairing and DSB repair are cilium-independent in mammals. However: (a) the paper does not directly measure telomere clustering at the NE or bouquet stage timing; (b) it is a preprint, not yet peer-reviewed; (c) the Mytlis 2023 paper (already in R1 REFRESH_08 notes) established the bouquet-MTOC/zygotene cilium association but also did not show functional necessity.

**Duplicate check:** Confirmed NOT present in REFRESH_08_meiotic_bouquet_envelope.md, REFRESH_09_hic_3d_methods.md, or any R1 REFRESH file (grep on PMID, DOI, "CEP164," and "zygotene cilia" returned no hits across all REFRESH files).

**Recommended action:** ADD as background note in REFRESH_08 Section 3 update or directly in the open-questions appendix of the draft. Cite alongside `bouquet_ZygoteneCilium2021` as evidence that the cilium is dispensable for the core meiotic chromosome events (pairing, DSB repair). Do not use to claim cilium dispensability for telomere clustering at the NE — that was not measured.

```bibtex
@article{bouquet_ChenCEP164Cilia2025,
  author    = {Chen, Jun Jie and Gong, Xiangyu and Mak, Michael and Li, Feng-Qian and Takemaru, Ken-Ichi},
  title     = {Loss of cilia drives centriole clustering and elimination during mammalian spermatogenesis},
  year      = {2025},
  journal   = {bioRxiv},
  doi       = {10.64898/2025.12.04.692363},
  pmid      = {41415467},
  note      = {Preprint (Stony Brook). CEP164 KO eliminates zygotene cilia and sperm flagella; meiotic chromosome pairing and DSB repair proceed normally — zygotene cilia are dispensable for these processes in mammals. Male infertility driven by flagellogenesis defect, not by absence of zygotene cilia. Keywords: CEP164, cilia, flagella, meiotic bouquet, spermatogenesis, zygotene.}
}
```

---

## Section 2: Backfill Where R1 Was Thin

**R1 REFRESH_08** proposed 6 new papers; R1 REFRESH_09 proposed 7 new papers. Neither falls below the 5-paper threshold for a "thin R1" requiring backfill. Both refreshes are classified as TRUST in the R2 audit.

**Did a deeper search find more papers?**

For Topic 08 (meiotic bouquet): Searches returned 0–2 relevant PMIDs per query in the 2024-06 to 2026-05 window (detailed in Section 5). No papers were found beyond those already in R1. The search confirms that the meiotic bouquet machinery literature published after June 2024 is genuinely sparse — the primary mechanistic papers (KASH5 2023, SpdyA 2025, CTCF 2025, centromere-trigger 2025) were published in 2023–2025 and captured in R1. No cluster of missed papers was identified.

For Topic 09 (Hi-C 3D methods): The key gap remains human stage-resolved meiotic Hi-C (confirmed absent from the 2024-06 to 2026-05 window). The Cheng 2024 NIH Micro-C preprint has not been published. The review landscape added Kitamura 2025 (Trends in Genetics, PMID 41407613) and Thadani 2026 (Cold Spring Harbor Perspectives, PMID 41419316) — both already in R1 REFRESH_09. No new primary data paper was found.

**Conclusion:** No backfill needed. R1 coverage is consistent with the actual publication rate in these topics.

---

## Section 3: Contradiction Follow-Up

**R1 flagged contradictions:** Neither REFRESH_08 nor REFRESH_09 flagged any contradictions in the 2023-2026 window. REFRESH_08 Section 4 ("Contradictions") is empty; REFRESH_09 Section 4 ("Contradictions") is explicit: "No papers found in the searched window (2023–2026) that directly contradict the draft's claims in topic 09."

**New contradictory literature in the last 6 months (2025-11 to 2026-05):**

No new papers were found that contradict the following draft claims:

1. **MAJIN–TERB2–TERB1 bridge telomeres to the inner nuclear membrane:** Yin 2024 (Andrology), Liu 2025 (EMBO J SpdyA/SUN1 Ser48), and Garner 2023 (KASH5 activating adaptor) all remain the most current literature; no 2025-2026 papers challenge or revise the TTM complex mechanism.

2. **SUN1–KASH5 spans the NE; dynein drives RPMs:** No new structural or functional work in 2025-2026 was found to contradict Garner 2023. The Liu 2025 SpdyA/CDK2 paper extends the SUN1-SPDYA axis to pachytene NE remodeling and is additive, not contradictory.

3. **Zygotene peaks the Mantel ρ = 0.715 mouse Hi-C result (Zuo 2021):** No new mouse meiotic Hi-C dataset appeared in 2025-2026 to challenge the stage-specific result. He et al. 2023 (Developmental Cell) confirmed the leptotene-to-zygotene chromatin reorganization is independent; the Cheng 2024 preprint is consistent with the Zuo 2021 loop-size picture.

4. **Human bulk Hi-C (GM12878, sperm) ties sequence similarity to nuclear-envelope proximity:** The Marín-Gual 2025 RAD21L paper (Science Advances) is consistent (disrupting bouquet cohesin increases non-specific telomeric contacts rather than eliminating the community-structured signal). No paper challenges this.

The new MEDIUM paper (Chen 2025 CEP164 KO) shows cilia are dispensable for chromosome pairing and DSB repair. This does not contradict any draft claim since the draft does not assert that zygotene cilia are mechanistically required for bouquet function.

**Conclusion: No contradictions, old or new. Section 3 is affirmatively silent.**

---

## Section 4: 60-Day Preprint Frontier

**Search window:** 2026-03-17 to 2026-05-17 (60 days before 2026-05-17).
**Categories searched on bioRxiv:** genetics, cell biology, molecular biology, genomics, developmental biology, evolutionary biology (100 results each = 600 preprints reviewed by title/abstract).
**Note:** bioRxiv API does not support keyword search; results are filtered by date + category only. Titles and abstract previews were scanned manually for: subtelomere, chromosome-end, NAHR, pedigree T2T, pangenome, meiotic bouquet, Hi-C, Pore-C, CiFi, zygotene, KASH5, SUN1, MAJIN, TERB.

**Results:**

No preprints in the last 60 days directly address:
- Meiotic bouquet biology (MAJIN, TERB, SUN1, KASH5, RPMs, zygotene, bouquet)
- Hi-C or 3D genome architecture in meiotic spermatocytes
- CiFi or Pore-C applied to subtelomeric/repetitive regions
- Human meiotic chromosome organization

**One preprint of marginal relevance found:**

**[Lorber et al. 2026, bioRxiv]** "Chromatin tethering to the nuclear envelope enhances its accessibility to RNAPII and promotes chromatin asymmetric organization." Lorber D, Azuri I, Kumar A, Rotkopf R, Safran S, Volk T. DOI: [10.64898/2026.03.06.710131](https://doi.org/10.64898/2026.03.06.710131). Posted 2026-03-07. Category: molecular biology. Weizmann Institute.

**Assessment: WEAK / DROP.** Uses high-resolution live imaging in *Drosophila* larval muscle nuclei to show that LINC complex–mediated chromatin tethering to the NE counteracts chromatin self-attraction and maintains 3D organization. Relevant mechanistic principle (LINC complex mediates chromatin-NE proximity) but in somatic Drosophila cells, not mammalian meiosis. Does not provide new meiotic or human data. Does not affect any draft claim or add to the meiotic bouquet mechanism. Dropped.

**Conclusion:** No preprints from the last 60 days affect topics 08 or 09.

---

## Section 5: Audit Trail

### PubMed Queries (all 2024-06 to 2026-05 unless noted)

| Query | Hits (raw) | Relevant new (not in R1) |
|---|---:|---:|
| `meiotic bouquet telomere nuclear envelope zygotene` | 2 | 0 (both plant/fish drops from R1) |
| `SUN1 KASH5 meiosis spermatocyte` | 0 | 0 |
| `human meiosis spermatocyte Hi-C chromatin conformation prophase 3D genome` | 0 | 0 |
| `single sperm Hi-C Dip-C 3D genome structure conformation chromatin` | 1 | 0 (PMID 42036683 = hic3d_Chen2026HiChew, R1 REFRESH_15 dup) |
| `TERB1 TERB2 MAJIN meiosis telomere` (2024-01 to 2026-05) | 1 | 0 (PMID 38511802 = bouquet_YinReview2024, R1 dup) |
| `subtelomere chromatin architecture 3D nuclear organization Hi-C` | 0 | 0 |
| `meiotic telomere bouquet chromosome pairing homolog prophase` | 5 | 0 (all plant/fish R1 drops + hic3d_xie2025 R1 dup) |
| `LINC complex nuclear envelope chromosome movement meiosis 2025` | 1 | 0 (PMID 40715639 = plant drop) |
| `meiotic chromatin loop synaptonemal complex TAD zygotene prophase I` | 0 | 0 |
| `spermatocyte meiosis Hi-C chromatin organization homolog pairing` | 0 | 0 |
| `meiosis oogenesis 3D genome chromosome conformation prophase female` | 0 | 0 |
| `germline 3D genome subtelomeric chromosome end nuclear lamina` | 0 | 0 |
| `meiosis mouse chromosome movement rapid prophase I cilia zygotene` | 0 | 0 |
| `meiotic bouquet zygotene chromosome` (2025-03 to 2026-05) | **1** | **1** (PMID 41415467 = Chen 2025 CEP164/cilia, **NEW MEDIUM**) |
| `human spermatocyte prophase chromosome conformation capture` | 0 | 0 |
| `SUN1 meiosis spermatocyte telomere azoospermia male infertility` | 0 | 0 |
| `Pore-C long-read chromatin conformation meiosis germline 3D` | 0 | 0 |
| `meiosis spermatogenesis 3D genome organization nuclear architecture 2025 2026` | 1 | 0 (PMID 41407613 = hic3d_kitamura2025, R1 dup) |
| `meiosis chromosomes hazards hubs telomere Cooper 2026` | 1 | 0 (PMID 41419316 = hic3d_thadani2026, R1 dup) |
| `Cheng Pratto Brick Camerini-Otero mouse meiosis chromatin Micro-C Hi-C` | 0 | 0 (preprint not yet published) |
| `sperm 3D genome chromatin conformation single cell Hi-C` | 1 | 0 (PMID 42036683, R1 dup) |
| `telomere clustering homolog pairing meiosis nuclear organization` | 1 | 0 (PMID 39283979 = plant drop) |

**Total unique PMIDs reviewed:** 16 non-duplicate records.
**New papers not in R1:** 1 (PMID 41415467, MEDIUM).

### bioRxiv Category Sweeps (2026-03-17 to 2026-05-17, 100 results each)

| Category | Results scanned | Relevant finds |
|---|---:|---|
| genetics | 100 | 0 |
| cell biology | 100 | 0 |
| molecular biology | 100 | Lorber 2026 (LINC/NE in somatic Drosophila — WEAK, dropped) |
| genomics | 100 | 0 |
| developmental biology | 100 | 0 |
| evolutionary biology | 100 | 0 |

**Total bioRxiv preprints reviewed:** 600 titles/abstract previews.
**New preprints: 0 relevant** (1 tangential/dropped).

### Duplicate Cross-Check for the One New Paper

PMID 41415467 (Chen 2025, CEP164/cilia) was verified NOT present in:
- REFRESH_08_meiotic_bouquet_envelope.md (grep: "41415467", "692363", "CEP164", "cilia.*zygot", "zygot.*cilia" — no hits)
- REFRESH_09_hic_3d_methods.md (same grep — no hits)
- REFERENCES_v4.bib (same grep — no hits)

### Rationale for Absence of STRONG Papers

The meiotic bouquet (topic 08) literature in 2024-06 to 2026-05 is thin because the major 2023-2025 mechanistic papers were:
- **KASH5 activating adaptor** (Garner 2023) — already in R1
- **SpdyA/SUN1 Ser48 phosphorylation** (Liu 2025) — already in R1
- **SUN1 human NOA variant** (Meng 2023) — already in R1
- **CTCF-loop dynamics human meiosis** (Kaiser 2025) — already in R1
- **Centromere-trigger bouquet** (Jiménez-Martín 2025) — already in R1

The Hi-C 3D methods (topic 09) literature in 2024-06 to 2026-05 added only reviews (Kitamura 2025, Thadani 2026) and one methods advance (HiChew 2026 = hic3d_Chen2026HiChew from R1 REFRESH_15). The key primary-data gap (human stage-resolved meiotic Hi-C) remains unfilled. This absence is informative: the draft's cross-species inference from mouse zygotene to human prophase-I remains the only available path as of 2026-05.

---

*Generated by agent-122 (r2-litref-3d-bouquet) on 2026-05-17.*
*All PMIDs verified via PubMed MCP retrieval. The one new paper (PMID 41415467) has its metadata confirmed through mcp__claude_ai_PubMed__get_article_metadata. DOIs are provided as required for proper attribution.*
