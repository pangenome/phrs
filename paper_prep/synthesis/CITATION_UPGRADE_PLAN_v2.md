---
title: Per-paragraph citation upgrade plan (round-2)
draft_target: paper_prep/synthesis/NATURE_DRAFT_v1.md (line range 18-84; line numbers as in v1)
bib_source: paper_prep/synthesis/REFERENCES_v5.bib (372 entries; +8 over v4)
round1_plan_already_applied: paper_prep/synthesis/CITATION_UPGRADE_PLAN.md (37 ADDs; landed during Pass B → NATURE_DRAFT_v2.md)
generated: 2026-05-17
agent: r2-integrator-merge (agent-136)
note: This plan adds only the 8 round-2 bibkeys produced by the R2 sweep. Every action below is round-2-only; round-1 ADDs are already applied and are not duplicated here. Target paragraph references use NATURE_DRAFT_v1.md line numbers and the same P1..P11 / M1..M17 anchors as the v1 plan.
---

# Per-paragraph citation upgrade plan (round-2)

This plan is the round-2 companion to `paper_prep/synthesis/CITATION_UPGRADE_PLAN.md` (round-1; already applied). It contains **only** the ADD / REPLACE / REMOVE actions arising from round-2 findings — i.e. the 8 new bibkeys integrated into `REFERENCES_v5.bib`. No action below appears in the round-1 plan; the two plans are fully disjoint in their action set.

Action verbs (same vocabulary as v1):
- **ADD**: insert a new citation; existing block remains intact.
- **REPLACE x → y**: swap a bib key (none in this round; see Summary).
- **REMOVE**: drop without replacement (none in this round; see Summary).

---

## Summary table — round-2 totals

| Paragraph | ADDs | REPLACEs | REMOVEs | Net | Source REFRESH_R2 |
|-----------|-----:|---------:|--------:|----:|-------------------|
| Abstract | 0 | 0 | 0 | 0 | — |
| P1 (history, founding refs) | 0 | 0 | 0 | 0 | — |
| P2 (methods substrate) | 0 | 0 | 0 | 0 | — |
| P3 (heatmaps, scale) | 0 | 0 | 0 | 0 | — |
| P4 (NJ tree) | 0 | 0 | 0 | 0 | — |
| P5 (3-architecture split) | 0 | 0 | 0 | 0 | — |
| P6 (heterogeneity, FST, TAR1) | 1 | 0 | 0 | +1 | foundational |
| P7 (3D Hi-C panel) | 0 | 0 | 0 | 0 | — |
| P8 (flanking paradox, bouquet) | 1 | 0 | 0 | +1 | recombination |
| P9 (pedigree exchange) | 3 | 0 | 0 | +3 | frontier (×2), recombination (×1) |
| P10 (RPE-1, mouse) | 0 | 0 | 0 | 0 | — |
| P11 (causal loop, limits) | 1 | 0 | 0 | +1 | frontier |
| EDFig.4 caption (gene enrichment) | 1 | 0 | 0 | +1 | frontier |
| EDFig.8 caption (causal loop / bouquet) | 1 | 0 | 0 | +1 | 3d-bouquet |
| Methods M14 (single-cell 3D) | 1 | 0 | 0 | +1 | recombination |
| **TOTAL** | **9** | **0** | **0** | **+9 ADD actions; 8 unique bibkeys** |

**Why 9 actions on 8 bibkeys.** `pedigree_Porsborg2025primaterecom` appears in two locations (P9 main-text sentence and a Methods/recombination-model anchor) because the R2 audit P-3 angle is split across the prose and the Methods recombination model. All other bibkeys are cited exactly once.

**Zero REPLACEs and zero REMOVEs.** Every R2 paper is additive: each fills a gap (mechanism, calibration, sub-nucleolar detail, dispersal mechanism, hotspot context) rather than superseding a v4 cite. The R1 plan also had zero true REPLACEs, so the pattern continues.

**No overlap with round-1 plan.** Cross-checked: every bibkey in this plan is unique to round-2 (not in `CITATION_UPGRADE_PLAN.md`) and was added to REFERENCES_v5 by this integrator task (not to v4 by R1).

---

## P6: Within-community heterogeneity, F_ST, TAR1 prevalence (NATURE_DRAFT_v1.md line 32)

**Already-cited TAR1/TERRA cluster (post-v1 + Pass B context):** `@Ambrosini2007`, `@subtelstruct_NergadzeITS2007`, `@subtelstruct_Nergadze2007`, `@subtelstruct_NergadzeITSReview2007`.

**Actions:**
- **ADD `@santagostino2025terra`** on the TAR1 prevalence + TERRA-promoter clause (v1 L32 — TAR1 prevalence sentence, or in Extended Data Fig. 3 caption). Suggested phrasing: "TAR1 prevalence is correlated with arm architecture (PAR1 0.3–1.1%, acrocentric p-arms 73–79%, other autosomal arms >99%), consistent with the T2T-based TERRA promoter map of subtelomeres (39 of 46 subtelomeres carry a TERRA promoter; [@santagostino2025terra])." Source: `REFRESH_R2_foundational.md` §1.1 (PMID 41193243; DOI 10.1261/rna.080790.125).

**Rationale.** This is the R2_AUDIT_PLAN P-1 angle landed at the natural insertion point. Quantitatively complementary to the draft's 18,352 ITS-island count (Santagostino reports 205 ITS sites transcribed in ≥1 cell line; the draft does not need to add this number but it supports the inference that ITS loci are not merely structural relics).

---

## P8: Flanking paradox + S_all negative control + bouquet mechanism (NATURE_DRAFT_v1.md line 36)

**Already-cited bouquet cluster (post-v1 + Pass B context):** `@bouquet_KotaSUN1MAJIN2020`, `@bouquet_Scherthan2001`, `@bouquet_Scherthan2003`, `@bouquet_ShibuyaRPMs2015`, `@bouquet_ChikashigeTelomere1994`, `@bouquet_HarperBouquet2004`, `@bouquet_HornKASH52013`, `@bouquet_DingSUN12007`, `@bouquet_MorimotoKASH2012`, `@bouquet_ZicklerKleckner1999`, `@ZicklerKleckner1998`, `@ZicklerKleckner2015`, plus Pass B additions (`@bouquet_GarnerKASH52023`, `@bouquet_LiuSPDYA2025`, `@bouquet_MengSUN1NOA2023`, `@bouquet_JimenezCentromere2025`, `@bouquet_KaiserCTCF2025`).

**Actions:**
- **ADD `@ataei2025line1dj`** on the C7 acrocentric / DJ-anchoring clause near v1 L36 (or move to EDFig.4 caption — see EDFig.4 entry below). Suggested phrasing for main text: "The five acrocentric distal junctions share >98% sequence identity across arms; a primate-specific LINE1 element conserved at all five DJs contributes to their nucleolar anchoring [@ataei2025line1dj]." Source: `REFRESH_R2_recombination.md` §1 (PMID 39797762; DOI 10.1101/gad.351979.124).

**Rationale.** This is the first molecular mechanism for the C7 co-localisation pattern observed in Dip-C and sperm scHi-C. The current P8 paragraph asserts mechanism via meiotic-bouquet tethering for autosomal arms; ataei2025line1dj closes the equivalent mechanism gap for the acrocentric p-arms specifically. Placement choice (P8 vs EDFig.4) is at integrator discretion; both are valid.

---

## P9: Pedigree-resolved exchanges (NATURE_DRAFT_v1.md line 38)

**Already-cited cluster (post-v1 + Pass B context):** `@Cechova2025`, `@concerted_evolution_nahr_SamonteEichler2002`, `@concerted_evolution_nahr_Eichler2001`, `@concerted_evolution_nahr_Hastings2009`, `@concerted_evolution_nahr_Myers2010`, `@Sharp2006`, `@StankiewiczLupski2002`, `@StankiewiczLupski2010`, `@Porubsky2025`, `@acrocentric_Porubsky2025denovo`, plus Pass B additions (`@pedigree_Schweiger2024spermNCO`, `@noyes2026sd`, `@chen2025paraphase`, `@Tardy2026fshd`).

**Actions:**
- **ADD `@pedigree_Porsborg2025primaterecom`** on the 133 gene-conversion-like clause (v1 L38; the sandwich-pattern sentence). Suggested phrasing (extending R2_AUDIT_PLAN.md §5 edit 8): "…consistent with the non-PRDM9 long-tract NCO class quantified by long-read sperm sequencing [@pedigree_Schweiger2024spermNCO]; the co-existence of crossover-associated gene conversion with substantially longer tracts (318–688 bp; [@pedigree_Porsborg2025primaterecom]) suggests that a subset of the larger patches may reflect crossover-associated rather than standalone NCO events." Source: `REFRESH_R2_frontier.md` §1 STRONG-1 (PMID 41285744; DOI 10.1038/s41467-025-65248-3). **HIGH PRIORITY** — R2 STRONG.
- **ADD `@pedigree_Joseph2024PRDM9indep`** as supporting cite alongside the Schweiger 2024 anchor (either in P9 prose or in Methods recombination-model commentary). Suggested phrasing: "PRDM9-independent recombination hotspots are present but suppressed in humans relative to other placental mammals [@pedigree_Joseph2024PRDM9indep], consistent with the rarity (~2%) of the long-tract non-PRDM9 NCO class detected in human sperm." Source: `REFRESH_R2_frontier.md` §1 MEDIUM-2 (PMID 38809707; DOI 10.1073/pnas.2401973121).
- **ADD `@sasani2026kfam`** in P9 Limitations/pedigree-context clause, or in v1 L42 limitation cluster (P11), to acknowledge that the **same CEPH1463/K1463 pedigree** has been independently mined for tandem-repeat mutagenesis at genome-wide scale by the Quinlan/Eichler/Porubsky group. Suggested phrasing for a P9 Limitations addendum: "The same CEPH1463 four-generation pedigree has been independently analysed for tandem-repeat expansion/contraction at 8 million STR/VNTR loci, identifying 1,270 expansions/contractions and 43 hyper-mutable sites [@sasani2026kfam]; the PHR-specific exchange events reported here occupy a complementary class of pedigree variation." Source: `REFRESH_R2_recombination.md` §2 (PMID 41959501; DOI 10.64898/2026.03.06.710071).

**Rationale.** P9 is the most-strengthened paragraph in the round-2 refresh. The three additions provide (a) tract-length calibration of the 133 gene-conversion-like patches (Porsborg), (b) evolutionary-genomic context for the rarity of the non-PRDM9 class in humans (Joseph), and (c) pedigree-level cross-validation that the same family is a productive substrate for de novo mutation discovery beyond PHRs (Sasani). Together they triangulate the recombination model in v2 P9 without changing the underlying claim count.

---

## P11: Causal loop, limitations, outlook (NATURE_DRAFT_v1.md line 42)

**Already-cited cluster (post-v1 + Pass B context):** 15 concerted-evolution / methodology bibkeys plus Pass B additions (`@palsson2025recomb`, `@hic3d_Chen2026HiChew`, `@pangenome_Loegler2025review`, plus the SEPTIN14P22/OR4F gene-enrichment cluster added in P11).

**Actions:**
- **ADD `@subtelstruct_Lee2026SEPTIN14`** on the SEPTIN14P22 hub-span sentence in P11 gene-enrichment paragraph (v1 L42, the SEPTIN14P22 named-gene clause). Suggested phrasing: "...the SEPTIN14P pseudogene family, whose subtelomeric dispersal across great apes proceeds via SD-mediated propagation rather than independent retrotransposition [@subtelstruct_Lee2026SEPTIN14]." Source: `REFRESH_R2_frontier.md` §1 MEDIUM-1 (PMID 41699652; DOI 10.1186/s13100-026-00394-z).

**Rationale.** P11's gene-enrichment paragraph names SEPTIN14P22 as a hub-spanning pseudogene but currently lacks a mechanism citation for why that pseudogene family is enriched at subtelomeres. Lee 2026 supplies the mechanistic basis.

---

## Extended Data Fig. 4 caption (gene enrichment)

**Action:**
- **ADD `@hao2024snul`** as an optional supporting cite in the EDFig.4 caption (NOR co-localisation / individual-arm-identity panel, if present). Suggested phrasing for the caption: "Sub-nucleolar territories are individually addressable per NOR-containing chromosome via monoallelically-expressed SNUL ncRNAs [@hao2024snul]." Source: `REFRESH_R2_recombination.md` §1.2 (PMID 38240312; DOI 10.7554/eLife.80684).

**Rationale.** Optional MEDIUM that adds nuance to the C7 co-localisation panel — each acrocentric arm has its own SNUL-defined sub-territory within the shared nucleolus. Skippable if EDFig.4 word budget is tight.

---

## Extended Data Fig. 8 caption (causal feedback loop / bouquet)

**Action:**
- **ADD `@bouquet_ChenCEP164Cilia2025`** in the EDFig.8a caption (causal loop) as a note in the bouquet / cilium open-question summary. Suggested phrasing: "Loss of zygotene cilia (CEP164 KO) eliminates ciliogenesis without disrupting meiotic chromosome pairing or DSB repair [@bouquet_ChenCEP164Cilia2025], indicating that the cilium is dispensable for the core meiotic chromosome events." Source: `REFRESH_R2_3d-bouquet.md` §1 (PMID 41415467; DOI 10.64898/2025.12.04.692363).

**Rationale.** Resolves the zygotene-cilium connection open question without changing any draft claim. Belongs in EDFig.8 caption (or open-questions appendix) rather than main text, since the draft does not currently assert cilium necessity. Caveat: telomere clustering at the NE was not directly measured by Chen 2025, so the placement should not be used to claim cilium dispensability for the bouquet specifically.

---

## Methods M14 — Single-cell 3D / pedigree context

**Already-cited (M14 line 72-76):** `@Tan2018`, `@Xu2025`, `@hic3d_scnanoHiC2023`, `@hic3d_scnanoHiC2_2025`, plus Pass B addition `@hic3d_kitamura2025`.

**Action:**
- **ADD `@pedigree_Porsborg2025primaterecom`** in Methods M16 (`pedigree odgi-untangle`, line 76) as a secondary anchor for the gene-conversion-like patch classifier. Suggested addendum: "Patch tract-length distributions are interpreted in the framework of long-read primate meiotic recombination, in which NCO tracts average 22–95 bp and CO-associated tracts average 318–688 bp [@pedigree_Porsborg2025primaterecom]."

**Rationale.** This is the second Porsborg 2025 anchor — Methods-only — to make the tract-length interpretation explicit in the methodological substrate, not just the prose.

---

## Verification checklist (every recommended ADD bibkey is present in REFERENCES_v5.bib)

```
ADD bibkey                            v5 presence  Source REFRESH_R2 file
santagostino2025terra                 OK           REFRESH_R2_foundational.md §1.1
ataei2025line1dj                      OK           REFRESH_R2_recombination.md §1
hao2024snul                           OK           REFRESH_R2_recombination.md §1.2
sasani2026kfam                        OK           REFRESH_R2_recombination.md §2
bouquet_ChenCEP164Cilia2025           OK           REFRESH_R2_3d-bouquet.md §1
pedigree_Porsborg2025primaterecom     OK           REFRESH_R2_frontier.md §1 STRONG-1
subtelstruct_Lee2026SEPTIN14          OK           REFRESH_R2_frontier.md §1 MEDIUM-1
pedigree_Joseph2024PRDM9indep         OK           REFRESH_R2_frontier.md §1 MEDIUM-2
```

All 8 bibkeys exist in `paper_prep/synthesis/REFERENCES_v5.bib` (entry count 372 = 364 v4 + 8). pybtex parse exit code 0.

---

## Notes for r2-fix-apply (downstream)

1. **No overlap with round-1 plan.** Every action in this file targets a bibkey introduced by REFERENCES_v5 (not present in v4); the round-1 plan (now landed in NATURE_DRAFT_v2.md) is fully orthogonal. Confirm by intersecting the action set: 0 bibkeys appear in both plans.
2. **Apply order.** Apply R2_AUDIT_PLAN.md §5 line edits (1–10) FIRST, then layer this plan's ADDs on top — the audit edits adjust prose wording at the same insertion points used here. Applying in the reverse order will cause merge conflicts on P9 L38 and EDFig.4 caption.
3. **Word-count delta.** Estimated +90 words across main text (P6, P8, P9, P11 combined) and +50 words across captions/Methods. Total: ~140 words. Below the 500-word budget headroom for a Nature Article.
4. **Bib housekeeping deferred.** R1 plan §3 flagged the `hprc_siren2025` naming bug (carries Schloissnig DOI) and the `subtelstruct_Sholes2022` wrong-authors issue. Neither was within the R2 sweep scope. Leave these to a follow-up bib-cleanup task.

---

*End of CITATION_UPGRADE_PLAN_v2.md. Generated 2026-05-17 by r2-integrator-merge (agent-136).*
