# Reviewer-response drafts: statistical controls and bibliography (D-M4, D-M5, D-M9, D-M12, D-Bib, D-PeerQ2b)

Response-letter material. Not manuscript prose.

## D-M4 — the within-Leiden statistic now has a null

`scripts/pedigree/monte_carlo_null_d_m4.py` run with 10,000 permutations that preserve the
per-arm patch-count marginals (so the arm-size-distribution objection the reviewer raised is
built into the null):

| set | observed | null mean | null 95% | p (enrichment) |
|---|---|---|---|---|
| published, 538 patches | 469/538 = 0.872 | 0.677 | [0.664, 0.693] | < 1e-4 |
| **corrected copies, 310 patches** | **293/310 = 0.945** | 0.840 | [0.832, 0.852] | < 1e-4 |

The reviewer's arithmetic worry ("if 70% of coverage sits in three communities the random null
could be 50–70%") is confirmed as the right magnitude for the null — 68% for the published set
and 84% for the corrected set — and the observation still sits at the 100th percentile of it.
Wilson 95% CIs: published 0.918 [0.892, 0.939]; corrected 0.945 [0.914, 0.965].
Outputs: `paper_prep/synthesis/ANALYSIS_D_M4.json`, `paper_prep/figures/fig4/null_distribution_d_m4.{pdf,png}`.

## D-M5 — arm-level Mantel replaces the non-independent per-pair Spearman

`scripts/hic/mantel_mouse_d_m5.py`: arm-level sequence similarity (1 − Jaccard distance from
`mouse_T2T/similarity_1Mb/mouse.dist_matrix.tsv`) against arm-level mean normalised Hi-C
contact, aggregated from the per-PHR-pair tables; row-and-column permutation, 10,000
replicates; all four meiotic stages, two window sizes:

| stage | 20 kb | 50 kb |
|---|---|---|
| leptotene | ρ = +0.41, p < 1e-4 | +0.29, p < 1e-4 |
| **zygotene** | **ρ = +0.58, p < 1e-4** | +0.34, p < 1e-4 |
| **pachytene** | **ρ = +0.63, p < 1e-4** | +0.31, p < 1e-4 |
| diplotene | ρ = +0.33, p < 1e-4 | +0.30, p < 1e-4 |

Permutation null is centred on 0 with a 95% range of ±0.11, so every stage is far outside it.
The correlation is *attenuated* relative to the reported per-pair value (0.715 → 0.58 at
zygotene, 0.63 at pachytene) exactly as expected once the non-independence the reviewer
objected to is respected. We report the arm-level Mantel values and drop the per-pair p-value.
Outputs: `scripts/hic/mantel_mouse_d_m5.json`.

## D-M9 — a real character-level bootstrap, and what it does to the tree claims

`scripts/cladistics/char_bootstrap_d_m9.R`, B = 1,000 replicates resampling the 15,089
signal-bearing PHRs with replacement, recomputing the arm-level Jaccard matrix per replicate,
rebuilding NJ and UPGMA trees, and scoring exact clade recovery:

| clade | members | monophyletic in full tree | NJ bootstrap | UPGMA bootstrap | previous "sensitivity" |
|---|---|---|---|---|---|
| PAR1 | Xp, Yp | yes | **100.0%** | **100.0%** | 100% |
| PAR2 | Xq, Yq | yes | **100.0%** | 99.9% | 100% |
| P10_18 | 10p, 18p | yes | **100.0%** | **100.0%** | 100% |
| ACRO_p | 13p, 14p, 15p, 21p, 22p | yes | 54.4% | **87.8%** | 100% |
| TIGHT_q | 22q, 21q, 19q, 1q, 13q, 17q | **no** | 0% | 0% | 100% |
| DUX4 | 4q, 10q | **no** | 0% | 0% | 100% |

Two things follow, and they are different in kind:

1. For the clades that are monophyletic in the reference tree, the character bootstrap is
   strong for PAR1, PAR2 and the 10p/18p pair (100%), moderate for the acrocentric p-arm group
   (54% NJ, 88% UPGMA). This is the honest replacement for the "100% bootstrap" language.
2. TIGHT_q and DUX4 are **not monophyletic in the reference NJ/UPGMA tree at all**, so their
   clade support is 0 by construction. The earlier 100% "sensitivity support" measured
   recovery of a *partition*, not monophyly of a *clade*; the two numbers were never
   comparable. In the manuscript these groups must be described as community assignments
   from the Leiden partition, and the DUX4/4q–10q relationship must rest on its direct
   sequence evidence (and the community statistics), not on tree-support language.

Outputs: `scripts/cladistics/d_m9/` (support summary, per-replicate matrices, meta).

## D-M12 — confidence intervals

`scripts/ci/results_d_m12.json` now carries, for every headline statistic:

| statistic | point estimate | 95% CI | method |
|---|---|---|---|
| pedigree within-Leiden, published | 0.918 (494/538) | [0.892, 0.939] | Wilson score |
| pedigree within-Leiden, corrected copies | 0.945 (293/310) | [0.914, 0.965] | Wilson score |
| F_ST, per superpopulation pair | 0.102–0.152 | block jackknife over arms (10 blocks) | jackknife |
| Mantel ρ, human Hi-C | 0.66 | bootstrap/permutation distribution | `mantel_bootstrap_ci.tsv` |
| mouse zygotene ρ (arm-level Mantel) | +0.58 | permutation 95% outside ±0.11 → p < 1e-4 | D-M5 above |

No point estimate in the pedigree or 3D sections is now reported without an interval.

## D-M6 — matched F_ST baseline (now analysed)

Subtelomeric Hudson F_ST against matched genome-wide continental values for the same five
superpopulation labels and the same estimator (Bhatia 2013 Table 1; Patterson 2012; block
jackknife over ten arm/community blocks):

| pair | subtelomeric F_ST (95% jk) | matched genome-wide | ratio | verdict |
|---|---|---|---|---|
| AFR-AMR | +0.102 [+0.020, +0.184] | 0.071 | 1.44 | equivalent |
| AFR-EAS | +0.153 [+0.063, +0.242] | 0.144 | 1.06 | equivalent |
| AFR-EUR | +0.108 [-0.026, +0.242] | 0.150 | 0.72 | equivalent |
| AFR-SAS | +0.103 [-0.003, +0.210] | 0.110 | 0.94 | equivalent |

No subtelomere-specific signature: differentiation is not elevated above genome-wide. Per-block
values span -0.08 to +0.35 (s.d. 0.12) and do not correlate with within-community sequence
sharing (r = +0.07, ten blocks), so the spread is heterogeneity rather than a gradient.

Reproduce offline: `python3 scripts/popgen/matched_fst_from_blocks.py` writes
`scripts/popgen/matched_fst_d_m6.tsv` from the repository's per-block F_ST table plus the
published baselines (the upstream `matched_fst_d_m6.py` needs a moosefs table and is the version
that also computes per-arm contributions).

LIMITATION WE STATE: local-ancestry inference is not validated in duplicated subtelomeric
sequence, so we do not build any argument on ancestry-deconvolved genotypes; the comparison is
at arm resolution with continental labels, which is what these regions support.

## D-Bib — bibliography hygiene

1. **Vollger2023 is not a duplicate.** `Vollger2023` is *Segmental duplications and their
   variation in a complete human genome*, Science 376:eabj6965 (2023), and
   `concerted_evolution_nahr_Vollger2023` is *Increased mutation and gene conversion within
   human segmental duplications*, Nature 617:325–334 (2023). Different papers, different DOIs;
   the reviewer's suspicion is resolved rather than a deletion made.
2. **Smith 1976 is already present** in `REFERENCES_v3.bib` (`Smith1976`, *Evolution of
   repeated DNA sequences by unequal crossover*, Science 191:528–535).
3. **Cech 2004 was missing from the synthesis bibliography** — added as `Cech2004`
   (*Beginning to understand the end of the chromosome*, Cell 116:273–279,
   doi:10.1016/s0092-8674(04)00038-8). It already exists in the submission bibliography as
   `cech2004chromosomeend`, so the manuscript can cite it without a build change.

## D-PeerQ2b — segmental-duplication generality (roadmap statement, not an analysis)

The question is whether the sequence-to-3D correlation is specific to subtelomeric PHRs or
general to segmental duplications. The honest answer for this revision: not yet measured, and
the comparison is on the roadmap. What it requires: the same arm-level Mantel machinery
applied to all non-subtelomeric HPRCv2 segmental duplications (~50–100 Mb) using the same
Hi-C pipeline, then comparing that ρ with the subtelomeric value. We give the reviewer the
number that would settle it (subtelomeric arm-level Mantel ρ = 0.66 human, +0.58/−0.63 mouse)
and state that the all-SD counterpart is pending.
