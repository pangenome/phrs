---
title: REVISION_LOG_v6 — deferred-analyses integration + bib hygiene
date: 2026-05-18
agent: draft-v6-integrate (agent-239)
prior_draft: paper_prep/synthesis/NATURE_DRAFT_v5.md
new_draft: paper_prep/synthesis/NATURE_DRAFT_v6.md
inputs:
  - paper_prep/synthesis/ANALYSIS_D_M4.md
  - paper_prep/synthesis/ANALYSIS_D_M5.md
  - paper_prep/synthesis/ANALYSIS_D_M6.md
  - paper_prep/synthesis/ANALYSIS_D_M9.md   (read from main: commit 7d8e378)
  - paper_prep/synthesis/ANALYSIS_D_M12.md  (read from branch: commit 72865c1)
  - paper_prep/synthesis/ANALYSIS_D_PEERQ1.md
  - paper_prep/synthesis/ANALYSIS_D_PEERQ3.md
  - paper_prep/synthesis/ANALYSIS_F34.md
  - paper_prep/synthesis/BIB_MERGE_v6_LOG.md
  - paper_prep/synthesis/REFERENCES_v6.bib
  - paper_prep/synthesis/RENDERED_REFERENCES_v6.md
word_budget:
  abstract: 200 / 200 (hard cap)
  main: 3299 / 3300 (hard cap)
  methods: 1591 (v5 1391; growth 200, within allowed 100-200)
---

# Revision log for NATURE_DRAFT_v6

## A. Analysis integrations

| Analysis | Status | Where applied | What changed |
|---|---|---|---|
| **D-M4** Monte Carlo null for pedigree 92% | APPLIED | P11 (pedigree paragraph) + Methods §Pedigree | Inserted null mean 77.0% (95% CI 75.4-78.8%; B = 10,000; one-sided p < 1e-4) and CEPH1463 null mean 30.7%, p ≈ 1e-4. Removed v5 sentence "The 92% lacks a published null baseline; a Monte Carlo permutation comparison ... is deferred". Methods §Pedigree gains the permutation-null procedure sentence. |
| **D-M5** Arm-level mouse Mantel per stage | APPLIED | P8 (mouse Hi-C paragraph) + Methods §Mouse pipeline + Fig. 4d caption | Replaced "Spearman ρ = 0.715 (n = 344 ... arm-level Mantel pending)" with "Mantel ρ = 0.718 (10,000 row+column permutations, p < 1 × 10⁻⁴; 95% CI 0.47-0.86; n_arms = 27)". Added per-stage trajectory (lepto 0.687, pachy 0.683, diplo 0.577). Methods sentence "Mantel ... is pending" removed. |
| **D-M6** Matched F_ST control | APPLIED | P6 (heterogeneity paragraph) + Methods §F_ST + Fig. 2 caption | Replaced "matched control deferred" with the matched 1000G/HGDP comparison + block-jackknife 95% CIs per pair (D-M12 §d numbers used: AFR-AMR 0.108 [0.026, 0.191]; AFR-EAS 0.155 [0.076, 0.234]; AFR-EUR 0.128 [-0.020, 0.275]; AFR-SAS 0.112 [-0.001, 0.224]). Verdict per pair (equivalent/elevated/depressed) embedded in Methods. |
| **D-M9** Character-level NJ bootstrap | APPLIED | Abstract + P4 (NJ paragraph) + Methods §Neighbour-joining tree + ED nj_tree_arms caption | Abstract: "robust under distance-matrix perturbation" → "with PAR1, PAR2, 10p/18p and 4q/10q recovered at >= 99% character-level bootstrap support" (compressed via grouping rewrite). P4 + Methods + ED caption updated with per-grouping support (PAR1/PAR2/10p18p 100%, 4q/10q 99.4%, acrocentric p 52% NJ / 87% UPGMA; tight q-arm 0.1% under chromosome-granular surrogate, deferred to follow-up). |
| **D-M12** Bootstrap CIs for headline correlations | APPLIED | P8 (mouse), P9 (human Hi-C Mantel), P11 (pedigree), P6 (F_ST) + Methods | Five CIs added in-text: Mantel CHM13 0.43-0.81, HG002 0.44-0.80; mouse zygotene 0.47-0.86; pedigree Wilson 89.2-93.9%; F_ST AFR-vs-non-AFR block-jackknife 0.026-0.234. §b (5 exclusion x 5 resolution Mantel CI grid) deferred — Open Reviewer Concerns. |
| **D-PeerQ1** Hi-C MAPQ-strict | APPLIED (partial) | Methods §Hi-C | Sentence added: strict-MAPQ re-binning is expected to reproduce flanking B/W within Poisson noise while PHR-internal contact collapses to noise floor. The numerical MAPQ-strict B/W table was not produced (`.allValidPairs` blocked on moosefs mount); the flanking 1.25-13.5x B/W ratio across 7/7 datasets is the falsification documented in-text. |
| **D-PeerQ3** P-arm/Q-arm asymmetry | APPLIED | P13 (gene-content paragraph) + Supp Table S-PQ reference | Replaced "a systematic P/Q orientation audit per community is deferred" with the D-PeerQ3 result: 17/75 within-community pairs are P-Q (vs 51% expected); C6, C7 significantly intra-orientation (p_BH < 0.005). Supplementary Table S-PQ referenced. |
| **F34** Per-meiosis crossover rate | APPLIED | P11 (pedigree paragraph) + Methods §Pedigree | Appended after the "Thirteen of sixteen crossover-like events are in PAN028" sentence: 0.76 / Mb / meiosis (95% CI 0.43-1.14; 16 events / 3 transmissions / 7 Mb of PHR per transmission), ~76x the genome-wide rate [@Halldorsson2019]. Methods §Pedigree gains denominator definition + Poisson-resampling CI procedure. |

## B. Narrative-match items (F-track) applied this round

| Item | Status | Where applied | What changed |
|---|---|---|---|
| **F20** DUX4 cancer / oncofetal-programme clause | APPLIED | P1 | Inserted after "degenerate D4Z4-like copies on at least ten additional chromosomes revealed by T2T-CHM13": "and is reactivated as an oncofetal programme in diverse solid cancers [@dux4_cancer_Chew2019]". New bibkey added in v6 bib pass. |
| **F21** Hi-C rare-contact justification | APPLIED | Methods §Hi-C | Justification of observed-over-expected normalisation for the rare inter-chromosomal regime (~2-5% of read pairs) added alongside the D-PeerQ1 strict-MAPQ sentence. |
| **F30** End P14 on directionality | DEFERRED (author judgement) | — | The directionality (chicken-and-egg) sentence sits mid-paragraph in v5 P14. Re-ordering to end on it would require restructuring around the Nature limitations-section convention; deferred. Not load-bearing for v6 scientific content. |
| **F32** "Simulate the full graph" closing clause | DEFERRED | — | Low-value polish (existing P2 closing sentence already conveys the substance). |

## C. New bibkeys cited in v6 vs v5

| Bibkey | Where cited in v6 | Source bib patch |
|---|---|---|
| `cech2004chromosomeend` | P7 (bouquet model), next to `bouquet_KotaSUN1MAJIN2020` and `ZicklerKleckner2015` | `BIB_PATCH_d_cech2004.bib` |
| `dux4_cancer_Chew2019` | P1 (oncofetal programme clause) | `BIB_PATCH_f20_dux4_cancer.bib` |
| `Halldorsson2019` | P11 (per-meiosis crossover-rate comparison) | already in REFERENCES_v5.bib (not a new bib entry, but newly cited in v6) |

`smith1976crossover` patch was skipped (duplicate paper of existing `Smith1976`); `Smith1976` remains in the bib but is not cited in v6 main text. The v6 NJ/Methods text uses 76 unique bibkeys, up from 73 in v5.

## D. Hard-gate verification

- Abstract: 200 / 200.
- Main text: 3299 / 3300.
- Methods: 1591 (v5 baseline 1391, allowed growth 100-200; growth = 200 ✓).
- Forbidden phrases (grep, zero hits): "arm-level Mantel test is pending"; "Monte Carlo permutation comparison ... is deferred"; "matched control deferred"; "arm-level Mantel pending".
- Em-dashes: 0 in body.
- Standalone `---` outside YAML frontmatter: 0.
- ED Fig. 6 numbering preserved (3 occurrences in figure list + main text).
- Leaked bibkeys vs REFERENCES_v6.bib: 0.
- In-text CIs added: 5 (Mantel HG002, Mantel CHM13, mouse zygotene Mantel, pedigree 92% Wilson, F_ST AFR-vs-non-AFR block-jackknife). Requirement: >= 4 ✓.

## E. Outputs

- `paper_prep/synthesis/NATURE_DRAFT_v6.md` (committed).
- `paper_prep/synthesis/REVISION_LOG_v6.md` (this file).
- `paper_prep/synthesis/RENDERED_REFERENCES_v6.md` (renumbered for 76 cited entries).
- `paper_prep/synthesis/OPEN_REVIEWER_CONCERNS.md` (appended with v6 state).

## F. Still-open items carried forward to v7+

See `OPEN_REVIEWER_CONCERNS.md` §v6-residual.

- D-M9 q-grouping at character level: requires per-arm-pair PGGB closure data (chromosome-granular surrogate in current bootstrap collapses the q-grouping to 0.1%).
- D-M12 §b: Mantel ρ trajectory under 5 exclusions x 5 resolutions; bootstrap CI grid (blocked on /moosefs mount).
- D-PeerQ1 numerical MAPQ-strict B/W table: blocked on `.allValidPairs` upstream files; flanking control falsification documented in-text.
- D-PeerQ2 RNA-soup speculation (cover letter only).
- D-PeerQ2b pan-segmental-duplication Mantel: separate paper.
- F30 P14 ending restructuring (author judgement).
- F32 "simulate the full graph" verbatim phrasing (low-value).
