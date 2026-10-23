# response/ — reviewer-response material

## Why this directory exists, and a correction about what is in it

The real peer review is **`PEER_REVIEW.md` at the repository root**: the Nature editor's letter
plus **five referee reports** (Referee #1 evolutionary genomics/recombination; #2; #3 genome
organization; #4; #5 human genomics), dated as received, with 6–8 week-scale revision
expectations. **That is the document that governs everything here.**

For part of this work the response material was drafted against the wrong document:
`paper_prep/synthesis/PEER_REVIEW_v1.md`, an internal **simulated** critique written by
agent-147 on 2026-05-17 against `NATURE_DRAFT_v3.md`. That file has now been renamed
`paper_prep/synthesis/SIMULATED_PEER_REVIEW_v1.md` and carries a SIMULATED warning at the top.

Consequence: the cross-references in the drafts below (labels like `M1`–`M12`, `m1`–`m18`,
`Q1`–`Q3`) refer to the **simulated** document's numbering, **not** to the real referees'
numbering. The technical content is still valid — the analyses were real and the numbers are
ours — but every "we answered point M4" claim must be re-mapped onto the actual referee
reports before anything is sent.

## Contents

| file | what it is |
|---|---|
| `01_evidence_log_pedigree.txt` | evidence log for the two pedigree objections raised outside the formal review (Slack thread with E. Garrison & M. Čechová, 2026-10-15): flagger QC, recurring-candidate argument, all identity/coverage numbers, revisions 1–14 of the reasoning |
| `02_letter_draft_v1.md` | assembled response letter (overview + 12 points + files changed). **Drafted against the simulated review.** |
| `02b_letter_draft_pedigree_points.md` | earlier, longer prose version of the pedigree points |
| `03_change_log.md` | item-by-item log of what was changed and why, with the numbers behind each change. This one is review-independent and is the most reliable record. |
| `04_answers_mapq0_and_rna.md` | answers to the multi-mapping (MAPQ0) and RNA-proximity questions |
| `05_answers_arm_asymmetry_and_architecture.md` | arm-type asymmetry (P/P, Q/Q, P-Q), module orientation, module order, module maps |
| `06_answers_stats_and_bib.md` | statistical controls (permutation null, Wilson/jackknife CIs, character bootstrap, arm-level Mantel), threshold sensitivity, bibliography hygiene |

Material that stayed elsewhere because it is canonical for the project rather than
response-specific: `PEDIGREE_NAMING.md` (naming and per-chromosome copy rules),
`scripts/ci/audit_response_numbers.py` (checks every quoted number against its source),
`submission/supporting_material/Fig5_margin_table.tsv`.

## Status of the re-mapping

`REAL_REVIEW_MAP.md` in this directory is the tracker keyed to the **actual** referee reports.
It starts as a skeleton extracted from `PEER_REVIEW.md`; each referee point gets a status
(`answered` / `partially` / `not addressed` / `disputed`) and a pointer to the artefact that
answers it. Until that file is filled in, no claim of "we addressed X" should be treated as
verified.

## Analysts' note on what is genuinely settled regardless of numbering

These are review-independent and already in the manuscript:

- pedigree label correction (per-chromosome maternal copy; child-side calls −73%) and the
  withdrawal of the child-side exchange claims;
- one retained ectopic exchange (chr9q terminus from chr3q) with heterozygous-donor evidence,
  threshold-invariance across FastGA frequency 10–8192, and HiFi flagger cleanliness;
- chr1p/chr5q reclassified as ancestral homology;
- mouse analysis replaced by an arm-level Mantel (ρ 0.41/0.58/0.63/0.33; pachytene peak);
- permutation null for the within-community fraction (0.68/0.84 nulls vs 0.872/0.945 observed);
- character-level bootstrap with the consequence that TIGHT_q and 4q/10q are **not** clades;
- matched F_ST baseline (ratios 0.72–1.44, all equivalent) plus block heterogeneity;
- module architecture (98.6% same telomere-relative orientation; repertoire 0.83 vs 0.57;
  order not conserved) as Extended Data Fig. 3.
