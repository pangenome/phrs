---
title: NATURE_DRAFT v3 revision log
generated: 2026-05-17
agent: r2-fix-apply (agent-140)
inputs:
  base: paper_prep/synthesis/NATURE_DRAFT_v2.md (v2; Pass A + Pass B applied)
  audit: paper_prep/synthesis/R2_AUDIT_PLAN.md (§5 REVISION_v3 TODO; 10 line edits)
  upgrade: paper_prep/synthesis/CITATION_UPGRADE_PLAN_v2.md (9 ADD actions on 8 unique bibkeys)
  contradictions: paper_prep/synthesis/LITERATURE_REFRESH_v2.md (§4; 5 contradictions, R2 status)
  references: paper_prep/synthesis/REFERENCES_v5.bib (372 entries; sha 8f777fec749c096977de84e52f19b34affc9e304)
output:
  draft: paper_prep/synthesis/NATURE_DRAFT_v3.md (244 lines; abstract 214 w; main 4158 w)
  delta_words: +221 main text (3937 -> 4158)
  delta_refs: +1 net (130 -> 131; +11 ADD bibkeys; -10 unused bouquet keys after P8 cluster trim)
---

# NATURE_DRAFT v3 revision log

Each entry lists action, source, location in v3, and status. Order: R2_AUDIT_PLAN §5 line edits 1-10 first, then CITATION_UPGRADE_PLAN_v2 ADD actions, then contradiction resolutions, then hedge audit. The plan author's instruction "apply audit edits FIRST, then layer upgrade plan ADDs" was followed exactly.

---

## R2_AUDIT_PLAN.md §5 line edits

| # | Edit | Source | v3 location | Status |
|---|------|--------|-------------|--------|
| 1 | "(minimum patch and alignment score 0.95)" -> "(min_score >= 0.8, 500 bp <= size <= 100 kb)" in Main text P9 (Pass A had fixed Methods L86 but missed Main mirror) | R2_AUDIT_PLAN §5 #1; report `14_pedigree_recombination.md` L20 | v3 L46 P9 (WashU pedigree filter) | DONE |
| 2 | OR4F sentence rephrase: "(10 OR4F family members on 14 arms each for OR4F5 and OR4F8P)" -> "with 10 OR4F family members in total; OR4F5 and OR4F8P are the most widespread, each present on 14 arms (Extended Data Fig. 4c)" | R2_AUDIT_PLAN §5 #2; report `03_gene_enrichment.md` L37 | v3 L50 P11 gene-enrichment | DONE |
| 3 | Append Chi 2025 hedge after OR4F sentence: "Olfactory receptor pseudogenisation across primates is now interpreted as a sensory-reallocation event rather than a simple visual-olfactory trade-off [@chi2025primate]." | R2_AUDIT_PLAN §5 #3; REFRESH_14 contradiction #4 | v3 L50 P11 gene-enrichment | DONE |
| 4 | Abstract softening: "tie sequence similarity to nuclear-envelope proximity through the meiotic bouquet" -> "tie sequence similarity to nuclear-envelope proximity, consistent with meiotic-bouquet repositioning" | R2_AUDIT_PLAN §5 #4; §3 item C-6 (abstract-vs-body tone mismatch) | v3 L26 Abstract last sentence | DONE |
| 5 | Trim P8 bouquet citation cluster from 16 -> 7 keys: keep KotaSUN1MAJIN2020, Scherthan2003, ShibuyaRPMs2015, HarperBouquet2004, ZicklerKleckner2015, GarnerKASH52023, KaiserCTCF2025. Remove (now uncited): Scherthan2001, ChikashigeTelomere1994, HornKASH52013, DingSUN12007, MorimotoKASH2012, ZicklerKleckner1999, ZicklerKleckner1998, LiuSPDYA2025, MengSUN1NOA2023, JimenezCentromere2025 | R2_AUDIT_PLAN §5 #5; §3 item C-7 (Nature citation density) | v3 L44 P8 bouquet mechanism cluster | DONE |
| 6 | Append TAR1/TERRA/telomere length anchor to P6 TAR1 prevalence sentence: "consistent with the T2T-based TERRA promoter map in which 39 of 46 subtelomeres carry a TERRA promoter [@santagostino2025terra] and with per-chromosome-end telomere length variation [@karimian2024telomereend]." | R2_AUDIT_PLAN §5 #6; §4 item P-1 (REFRESH_02 DEEPEN) | v3 L41 P6 TAR1 prevalence | DONE (bibkey resolved to `karimian2024telomereend`; `@santagostino2025terra` added during this edit covers CITATION_UPGRADE_PLAN_v2 P6 ADD) |
| 7 | Add `@makova2024apesex` to P1 PAR1/PAR2 citation cluster | R2_AUDIT_PLAN §5 #7; §4 item P-2 (REFRESH_04 DEEPEN, T2T PAR anchor) | v3 L30 P1 (PAR introduction) | DONE |
| 8 | Append two-process NCO model to P9 gene-conversion-like sentence: "and consistent with the non-PRDM9 long-tract NCO class quantified by long-read sperm sequencing [@pedigree_Schweiger2024spermNCO]" | R2_AUDIT_PLAN §5 #8; §4 item P-3 (REFRESH_07/15 long-tract NCO) | v3 L47 P9 gene-conversion-like | DONE (extended further via CITATION_UPGRADE_PLAN_v2 P9 ADDs: Porsborg 2025 tract length + Joseph 2024 PRDM9-indep context) |
| 9 | Append section anchor to Methods M14 PBMC sentence: "(`05_hic_validation.md` §PBMC, L455-469)" | R2_AUDIT_PLAN §5 #9; §3 item C-2 (PBMC number unsourced) | v3 L82 Methods M14 (single-cell 3D) | DONE |
| 10 | Append Methods anchor to P11 closing null-result sentence: "(Extended Data Fig. 4d; per-community 116 BH-corrected Fisher tests in Methods.)" | R2_AUDIT_PLAN §5 #10; §3 item C-8 (no community-specific gene signature) | v3 L50 P11 gene-enrichment closing | DONE |

**Audit edits scorecard.** 10 of 10 DONE; 0 SKIPPED.

---

## CITATION_UPGRADE_PLAN_v2.md ADD actions

| # | ADD | Target | v3 location | Status |
|---|-----|--------|-------------|--------|
| U-1 | `@santagostino2025terra` on TAR1 prevalence + TERRA promoter clause | P6 (v1 L32 -> v3 L41) | v3 L41 P6 (TAR1 prevalence sentence) | DONE (combined with audit edit #6 in single edit) |
| U-2 | `@ataei2025line1dj` on C7 acrocentric/DJ-anchoring (P8 or EDFig.4) | P8 (v1 L36) | v3 L44 P8 (post-bouquet mechanism sentence: "For the C7 acrocentric p-arm community a complementary somatic anchor is available...") | DONE (placed in P8 main text alongside the bouquet mechanism; v3 prose anchors the C7 mechanism story to LINE1 DJ) |
| U-3 | `@pedigree_Porsborg2025primaterecom` on P9 gene-conversion-like clause (tract-length calibration) | P9 (v1 L38 -> v3 L47) | v3 L47 P9 ("the co-existence of crossover-associated gene conversion at substantially longer tracts (318 to 688 bp; [@pedigree_Porsborg2025primaterecom])") | DONE |
| U-4 | `@pedigree_Joseph2024PRDM9indep` as supporting cite alongside Schweiger 2024 anchor | P9 prose | v3 L47 P9 ("while PRDM9-independent recombination hotspots remain suppressed in humans relative to other placental mammals [@pedigree_Joseph2024PRDM9indep]") | DONE |
| U-5 | `@sasani2026kfam` on P9 Limitations / pedigree-context clause | P9 Limitations / P11 outlook | v3 L47 P9 close ("The same CEPH1463 four-generation pedigree has been independently analysed for tandem-repeat expansion and contraction at roughly 8 million STR and VNTR loci...") | DONE |
| U-6 | `@subtelstruct_Lee2026SEPTIN14` on SEPTIN14P22 hub-span sentence | P11 (v1 L42 -> v3 L50) | v3 L50 P11 gene-enrichment ("with SEPTIN14P pseudogene dispersal across great apes proceeding via SD-mediated propagation rather than independent retrotransposition [@subtelstruct_Lee2026SEPTIN14]") | DONE |
| U-7 | `@hao2024snul` on Extended Data Fig. 4 caption | EDFig.4 caption | v3 L238 Figure list ("Sub-nucleolar territories are individually addressable per NOR-containing chromosome via monoallelically-expressed SNUL ncRNAs [@hao2024snul].") | DONE |
| U-8 | `@bouquet_ChenCEP164Cilia2025` on Extended Data Fig. 8 caption | EDFig.8 caption | v3 L240 Figure list ("Loss of zygotene cilia (CEP164 KO) eliminates ciliogenesis without disrupting meiotic chromosome pairing or DSB repair [@bouquet_ChenCEP164Cilia2025]...") | DONE |
| U-9 | `@pedigree_Porsborg2025primaterecom` second anchor in Methods M16 | Methods M16 (pedigree odgi-untangle) | v3 L86 Methods M16 ("Patch tract-length distributions are interpreted in the framework of long-read primate meiotic recombination, in which NCO tracts average 22 to 95 bp and CO-associated tracts average 318 to 688 bp [@pedigree_Porsborg2025primaterecom].") | DONE |

**Upgrade plan scorecard.** 9 of 9 ADD actions DONE; 8 unique bibkeys (Porsborg appears at both P9 and Methods M16 anchors as planned). 0 REPLACEs, 0 REMOVEs.

---

## Contradictions (LITERATURE_REFRESH_v2.md §4)

| # | Contradiction | R2 status | v3 disposition |
|---|--------------|-----------|----------------|
| 1 | Salsi 2026 (D4Z4 on >=10 chromosomes vs canonical 4q/10q) | Hedge already landed at v2 L30 with `@Salsi2026fshd` | No additional edit; carried forward unchanged into v3 L30 |
| 2 | Chi 2025 OR pseudogenisation as sensory-reallocation (vs Gilad 2004) | REDO scheduled | RESOLVED in v3 via R2_AUDIT_PLAN §5 edit #3 (added Chi 2025 hedge sentence in P11) |
| 3 | Lalli 2025 cM/Mb anti-correlation collapse | Honest null already cited in v2 P11 limitations (vi) | No additional edit; carried into v3 L53 |
| 4 | Porubsky 2025 no whole-genome crossover-SV correlation | Already addressed in v2 P9 close | No additional edit; carried into v3 L47 |
| 5 | Rodrigues 2024 "more than half" vs Gershman 2022 "roughly half" TERRA promoter count | Resolved by Santagostino 2025 to precise 39/46 | RESOLVED in v3 via R2_AUDIT_PLAN §5 edit #6 / CITATION_UPGRADE_PLAN_v2 U-1 (Santagostino added on TAR1 prevalence sentence; precise count cited) |

**Contradictions scorecard.** 2 of 5 newly resolved in v3; 3 carried forward (already addressed in v2). No active unresolved contradictions remain.

---

## Hedge audit

R2_AUDIT_PLAN §6 VERDICT recorded the v2 hedge state as: "3 'consistent with', 1 'may' in Methods limitations, 0 'likely'/'believe'/'appears to'". v3 adds 3 net "consistent with" occurrences (TERRA promoter map; non-PRDM9 NCO model; abstract softening per audit edit #4) and 1 "may" ("may reflect crossover-associated rather than standalone NCO events" in P9; this is a real tract-length inference based on Porsborg 2025 data, not a softener). All instances are scientifically appropriate hedges where the data permits but does not directly demonstrate the mechanism. No "likely", "we believe", or "appears to" occurrences introduced.

| Hedge | v2 count | v3 count | Notes |
|-------|---------:|---------:|-------|
| consistent with | 3 | 6 | 3 added: P6 TERRA (audit/upgrade-driven), P9 NCO model (audit-driven), abstract softening (audit-driven) |
| may | 1 | 2 | 1 added: P9 "may reflect crossover-associated" (data-supported by Porsborg tract-length distribution; not over-claim) |
| likely | 0 | 0 | not introduced |
| we believe | 0 | 0 | not introduced |
| appears to | 0 | 0 | not introduced |

No hedge-tightening edits were required by R2 audit (its verdict noted hedges already at minimum). The 3 new "consistent with" occurrences are mandated by the audit's own edits (the audit explicitly wrote "consistent with meiotic-bouquet repositioning" as the replacement for the abstract). Accepting the audit's own wording.

---

## YAML frontmatter changes

| Field | v2 | v3 |
|-------|----|----|
| `version` | v2 (Pass B content additions + REFERENCES_v4 integration) | v3 (R2 audit fixes + REFERENCES_v5 round-2 citation upgrades) |
| `abstract_words` | 214 | 214 (unchanged) |
| `main_text_words` | 3937 | 4158 (+221) |
| `pass_history` | 3 entries | 4 entries (added v3) |
| `references` | `paper_prep/synthesis/REFERENCES_v4.bib` | `paper_prep/synthesis/REFERENCES_v5.bib (sha: 8f777fec749c096977de84e52f19b34affc9e304)` |
| `citation_plan` | `paper_prep/synthesis/CITATION_UPGRADE_PLAN.md` | `paper_prep/synthesis/CITATION_UPGRADE_PLAN_v2.md` |
| `lit_refresh` | `paper_prep/synthesis/LITERATURE_REFRESH_v1.md` | `paper_prep/synthesis/LITERATURE_REFRESH_v2.md` |
| `audit` | `paper_prep/synthesis/CONSISTENCY_AUDIT_v1.md` | `paper_prep/synthesis/R2_AUDIT_PLAN.md` |
| `generated` | 2026-05-17 | 2026-05-17 (same day) |

---

## References section reconciliation

Inline citations and References section verified identical after edits:
- v2 had 130 inline + 130 References.
- v3 has 131 inline + 131 References (delta +1 net = +11 added - 10 unused after P8 cluster trim).
- Added (11): `ataei2025line1dj`, `bouquet_ChenCEP164Cilia2025`, `chi2025primate`, `hao2024snul`, `karimian2024telomereend`, `makova2024apesex`, `pedigree_Joseph2024PRDM9indep`, `pedigree_Porsborg2025primaterecom`, `santagostino2025terra`, `sasani2026kfam`, `subtelstruct_Lee2026SEPTIN14`.
- Removed from References (10; no longer cited inline after R2_AUDIT §5 #5 P8 cluster trim): `bouquet_ChikashigeTelomere1994`, `bouquet_DingSUN12007`, `bouquet_HornKASH52013`, `bouquet_JimenezCentromere2025`, `bouquet_LiuSPDYA2025`, `bouquet_MengSUN1NOA2023`, `bouquet_MorimotoKASH2012`, `bouquet_Scherthan2001`, `bouquet_ZicklerKleckner1999`, `ZicklerKleckner1998`.

Validation (run during this task):
- `comm -23 inline.txt bib_v5.txt` returned empty (0 leaked keys against REFERENCES_v5.bib).
- `diff inline_sorted.txt references_section_sorted.txt` returned empty (perfect 1-to-1 mapping inline <-> References).

---

## Validation summary (task validation checklist)

- [x] `NATURE_DRAFT_v3.md` exists, will be committed.
- [x] Every REVISION_v3 TODO has a corresponding REVISION_LOG_v3 entry marked DONE (10 of 10 DONE; 0 SKIPPED).
- [x] Every CITATION_UPGRADE_PLAN_v2.md action has a log entry (9 of 9 DONE).
- [x] grep verifies: zero citation keys in v3 draft that don't exist in REFERENCES_v5.bib (`comm -23` returned 0 lines).
- [x] Main text word count in [3000, 5000]: 4158 words.
- [x] Abstract retains its claims, has zero citations (verified with grep).
- [x] Pass A and Pass B fixes still present (not regressed): diff against v2 shows only the documented v3 edits, all v2 prose retained or extended with new ADDs.
- [x] Zero em-dashes.
- [x] Zero `---` outside YAML (verified: only L1 and L21).

*End of REVISION_LOG_v3.md. Generated 2026-05-17 by r2-fix-apply (agent-140).*
