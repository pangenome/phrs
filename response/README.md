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

## Locus-level report (read-level + multi-alignment)

- **`notes/9q3q_locus_full_report.md`** — the complete record for the one retained ectopic
  exchange: pipeline and metric definitions (dice/jaccard/cosine, `estimated.identity`, gap-compressed
  `dv`, length-penalised cols 10/11), the seven candidate explanations with verdicts, ONT read
  support at every involved end (father chr9 hap1/hap2, chr3 hap1/hap2; daughter paternal/maternal
  chr9) including the junction-spanning reads that refute assembly over-duplication, the
  per-window multi-alignment dice profile showing the homolog→paralog switch at the tract
  boundary (and its absence at 5q/1p), flagger status with corrected `Dup`/`Col` semantics, a
  corrections log of claims retracted during this analysis, and open items.
  Companion: `notes/read_support.json`, `scripts/pedigree/ont_read_support_9q3q.py`.

## Heterolog / coalescence / telomere analyses (2026-11)

Six notes added after the first response round. Read them in this order:

| note | what it establishes |
|---|---|
| `notes/hprc_heterolog_divergence.md` | cohort-scale allelic vs paralogous divergence across 465 HPRCv2 assemblies (allelic 0.15%, paralog 0.89%, ratio 5.93x; intra vs inter chromosomal indistinguishable at 1.07x); **read the "READ THIS FIRST" banner** - early sections are coarse-bin and superseded by the FINAL NUMBERS section |
| `notes/community_tier_table.md` | the three regimes: allelic-like (C7 acrocentric, PAR1, PAR2), near-allelic cohesive (C1, C11, C3), weak, and residual (C4) - source material for Extended Data Table `ED_community_tier_table` |
| `notes/exchange_rate_estimates.md` | exchange rate and transfer size implied by the divergences (1.35e-6 per locus pair per generation; tract ~50 kb; ~17% of the obligate rate for the average block) |
| `notes/community_coalescence_times.md` | per-community coalescence times, distributions, skew (5-12 Myr upper bound, 2-5 Myr at mu = 3e-8) |
| `notes/telomere_distance_profile.md` | ratio vs distance from the telomere-subtelomere junction: step structure with breakpoints near 10 and 20 kb |
| `notes/telomere_gene_vs_intergenic.md` | gene vs intergenic split; the terminal ~6 kb is **unmeasurable** because the allelic denominator fails there |
| `notes/lit_Tan2026_methods_extract.md` | the companion paper (Tan et al., bioRxiv 10.64898/2026.09.03.748401): verified full-text extract; they do **not** date by molecular clock, and both places where an early draft let their paper soften ours are retracted in-section |

Figures: `notes/heterolog_scope.pdf`, `notes/telomere_profile.pdf`, `notes/coalescence_distribution.pdf`,
`notes/telomere_end_curves.pdf`. Code: `scripts/ci/hprc_heterolog_divergence.py`,
`telomere_distance_profile.py`, `telomere_gene_profile.py`, `telomere_end_ratio_curves.py`,
`make_heterolog_figure.py`, `make_telomere_profile_figure.py`, `make_end_curve_figure.py`,
`scripts/cladistics/community_coalescence.py`.

**Outstanding (see `ROADMAP_analysis_integration.md`):** Methods and Results text for all of the above
(patch files, not yet applied); the tip control on long/high-identity records; the primate extension.
