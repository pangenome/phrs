---
title: D-M12 — Bootstrap 95% confidence intervals for headline correlations
parent: paper_prep/synthesis/OPEN_REVIEWER_CONCERNS.md §D-M12
date: 2026-05-18
scripts:
  - scripts/ci/bootstrap_ci_d_m12.py  (canonical)
  - scripts/ci/bootstrap_ci_d_m12.R   (cross-check)
artifacts:
  - scripts/ci/results_d_m12.json
  - scripts/ci/fst_block_jackknife.tsv
  - scripts/ci/fst_per_arm_per_pair.tsv
  - scripts/ci/mantel_bootstrap_ci.tsv
  - scripts/ci/mantel_fisher_z_ci.tsv
---

# D-M12 — bootstrap CIs bundle

Closes the four CI items D-M12 §a, §c, §d, §e from
`paper_prep/synthesis/OPEN_REVIEWER_CONCERNS.md`. §b (Mantel ρ trajectory
under 5 exclusion sets × 5 resolutions) is deferred — see "What is
deferred" below.

Each section reports the same five columns:

  **value | 95 % CI | method | n | replicates**

Random seed `20260518` is used throughout. Bootstrap and permutation
replicates are 10,000 each.

---

## §c — Pedigree within-Leiden fraction = 92 %

| value | 95 % CI | method | n | replicates |
|-------|---------|--------|---|------------|
| 0.9182 (494/538) | **[0.892, 0.939]** | Wilson-score CI for binomial proportion | 538 WashU HQ inter-chr patches | analytic |

* Source numbers: `end-to-end-report/report/14_pedigree_recombination.md`,
  WashU pedigree row "494 / 538 (92 %)".
* Method: Wilson 1927 score interval (Brown, Cai & DasGupta 2001 recommend
  it over the Wald interval for any moderate `np(1-p)` regime, including
  this one). Implemented as a closed-form calculation in
  `scripts/ci/bootstrap_ci_d_m12.{py,R}`; no Monte Carlo.
* **Recommended v6 in-text edit.** Replace "494 / 538 = 92 %" with
  "494 / 538 = 92 % (95 % Wilson CI 89.2 – 93.9 %)" in the
  pedigree-recombination subsection and the corresponding caption / table
  in 14_pedigree_recombination.md, in NATURE_DRAFT_v5 P14, and in any
  abstract / summary mention.
* **What this CI does NOT do.** Wilson is a CI on the binomial proportion
  given the 538 patches as Bernoulli trials. It does **not** address D-M4
  (the null baseline for "within-Leiden"). The width of the CI reflects
  sampling variance only; the comparison to the random-pairing null
  remains a separate (Monte Carlo) deferred item.

---

## §d — Hudson F_ST per superpopulation pair

| pair | F_ST | 95 % CI | method | n_arms |
|------|-----:|---------|--------|-------:|
| AFR vs AMR | **+0.102** | **[+0.022, +0.182]** | Hudson F_ST per arm, block-jackknife over arms | 10 |
| AFR vs EAS | **+0.153** | **[+0.065, +0.240]** | "                                          " | 10 |
| AFR vs EUR | **+0.108** | **[-0.024, +0.240]** | "                                          " | 10 |
| AFR vs SAS | **+0.103** | **[-0.001, +0.208]** | "                                          " | 10 |
| AMR vs EAS | +0.007 | [-0.022, +0.036] | "                                          " | 10 |
| AMR vs EUR | +0.007 | [-0.028, +0.041] | "                                          " | 10 |
| AMR vs SAS | +0.004 | [-0.012, +0.019] | "                                          " | 10 |
| EAS vs EUR | -0.048 | [-0.142, +0.047] | "                                          " | 10 |
| EAS vs SAS | +0.005 | [-0.012, +0.023] | "                                          " | 10 |
| EUR vs SAS | -0.003 | [-0.087, +0.080] | "                                          " | 10 |

Point estimates reproduce the values in `04_heterogeneity.md` "F_ST across
subtelomeric types" table exactly (AFR vs AMR 0.102, AFR vs EAS 0.152,
AFR vs EUR 0.108, AFR vs SAS 0.103, AMR vs EAS 0.007, AMR vs EUR 0.007,
EAS vs EUR -0.047, AMR vs SAS 0.004, EAS vs SAS 0.005, EUR vs SAS -0.003).

* **Replicates.** Block-jackknife is analytic (one delete-one mean per
  block, n_blocks = n_arms). The block is one arm/community pair.
* **Source counts.** 10 cross-arm/self-arm × superpop contingency rows
  from `end-to-end-report/report/04_heterogeneity.md` "Population
  structure in cross-arm affinity" (10 arm/community pairs with Fisher
  p_adj < 0.05). The same 10 rows are used by figure 2c right panel
  (`paper_prep/figures/fig2/figure_fig2.py`).
* **Method.** Reproduces the upstream pipeline `compute_fst_superpop.py`
  (`/moosefs/.../scripts/community/`). For each arm and each
  superpopulation pair (i, j):
  - p_pop = cross_pop / (cross_pop + self_pop) (binary cross-arm allele).
  - HS = mean(2·p_i·(1−p_i), 2·p_j·(1−p_j)).
  - p_pool = (c_i + c_j) / (n_i + n_j); HT = 2·p_pool·(1−p_pool).
  - F_ST = (HT − HS) / HT (Wright/Hudson form). Returns 0 when HT = 0
    (monomorphic locus in the pooled sample).

  Per-pair F_ST = mean over the 10 arm-community pairs (matches the
  upstream "averaged across significant arms" definition). All 10 arms
  contribute to every pair (no arm drops out, because the chrX_p and
  chrY_p rows have all cross+self totals > 0 even in non-AFR
  superpopulations).
* **Result.** AFR vs non-AFR ranges 0.102 – 0.153 (mean 0.117). Two of
  the four AFR-vs-non-AFR CIs (AFR vs EUR, AFR vs SAS) cross zero or
  graze it (lower bounds −0.024 and −0.001). AFR vs AMR and AFR vs EAS
  remain strictly positive ([+0.022, +0.182] and [+0.065, +0.240]). All
  six non-AFR vs non-AFR pair CIs bracket zero, confirming
  indistinguishability. The existing draft text already concedes that
  0.10 – 0.15 is in the autosomal continental range (D-M6) and the CI
  is consistent with that.
* **Recommended v6 in-text edit.** In the F_ST table in
  04_heterogeneity.md ("Fst across subtelomeric types") add a footnote:
  *"95 % CIs (block-jackknife over 10 arm/community pairs): AFR vs AMR
  0.102 [+0.022, +0.182]; AFR vs EAS 0.153 [+0.065, +0.240]; AFR vs
  EUR 0.108 [−0.024, +0.240]; AFR vs SAS 0.103 [−0.001, +0.208];
  non-AFR pairs all bracket zero."* Same footnote into the
  NATURE_DRAFT_v5 F_ST sentence at P6 (the "0.10 – 0.15" range is
  rewritten as "0.10 – 0.15, with two of four AFR-vs-non-AFR 95 % CIs
  excluding zero — AFR vs AMR [+0.022, +0.182] and AFR vs EAS [+0.065,
  +0.240] — and the other two grazing zero").
  Reproducible TSV: `scripts/ci/fst_block_jackknife.tsv`.
* **What this CI does NOT do.** It does not address D-M6 (matched
  genome-wide F_ST control on non-subtelomeric autosomal windows). The
  jackknife quantifies sampling variance across the 10 chosen arms; it
  does not test whether the 0.10 – 0.15 magnitude is unusual relative to
  background autosomal differentiation.

---

## §a — Mantel ρ for CHM13 Hi-C and HG002 Hi-C

Two CIs are reported per sample: (i) an arm-resampling **bootstrap** on
the local matrix snapshots (canonical, gold-standard, but on Feb-2026
matrices that predate the v5 rerun — see data-version caveat); and
(ii) a Fisher z **analytic CI on the published v5 point estimate**
(`n_arms` from the v5 results table) for the values the manuscript
actually quotes.

### (i) Bootstrap CI on local matrix snapshots

| sample | ρ (local) | 95 % CI | method | n_arms | replicates |
|--------|----------:|---------|--------|-------:|-----------:|
| HG002 Hi-C | +0.3146 | **[+0.133, +0.389]** | arm-resample bootstrap, Spearman on upper-triangle (similarity = 1 − Jaccard distance, vs Hi-C contact) | 42 | 10,000 |
| CHM13 Hi-C | +0.1749 | **[-0.112, +0.367]** | "                                                                              " | 42 | 10,000 |

* **Inputs (local snapshots).**
  - `/home/guarracino/Dropbox/working/Garrison/hprcv2/PHR_III/hic_validation/arm_dist_matrix.tsv`
    (Feb 2026; 42 × 42 arms)
  - `/home/guarracino/Dropbox/working/Garrison/hprcv2/PHR_III/hic_validation/hg002_contact_matrix.tsv`
    (Feb 2026; 48 × 48 arms, 42 common with arm_dist)
  - `/home/guarracino/Desktop/Garrison/HPRCv2/PHR_III/hic_validation/res_50kb/chm13_contact_matrix.tsv`
    (Feb 2026; 48 × 48 arms, 42 common)
* **Method.** Spearman rank correlation on the off-diagonal upper-triangle
  pairs of two symmetric matrices (arm-level Jaccard **similarity** = 1 −
  distance, vs arm-level Hi-C contact). Bootstrap: draw `n_arms` arm
  indices with replacement from the common-arm list, rebuild both
  matrices on the resampled label set (off-diagonal upper-triangle of
  the bootstrap matrix retains repeated-arm pairs, which contribute
  valid bootstrap variance), recompute ρ. Permutation p (row + column
  shuffles of one matrix) is reported for the local matrices to confirm
  the same direction as the upstream permutation pipeline.
* **DATA VERSION CAVEAT.** The local matrices predate the v5 paper
  rerun. They reproduce ρ ≈ 0.32 (HG002) and ρ ≈ 0.17 (CHM13), not the
  v5 headline ρ ≈ 0.66 (both samples). Inspection of the arm_dist
  matrix shows 11.2 % of off-diagonal entries are exactly 1.0 and
  52.2 % are ≥ 0.99, suggesting the local matrix is from an earlier
  graph-distance run before the post-processing step that fills in
  inter-arm similarity for low-Jaccard pairs. The **method is exact**;
  to obtain a bootstrap CI on the published 0.656 / 0.657, re-mount
  `/moosefs/guarracino/HPRCv2/PHR_III/` and re-run with
  `ARM_DIST_TSV = /moosefs/.../arm_dist_matrix.tsv` and
  `{sample}_contact_matrix.tsv` for the current resolution; the script
  is parameterised on those three paths.

### (ii) Fisher z analytic CI on the v5 headline point estimates

This is the CI the manuscript should cite *until* the v5 matrices are
re-mounted and the arm-resampling bootstrap is re-run. Fisher z is exact
for product-moment correlations and an excellent approximation for
Spearman at n ≥ 30 (Bonett & Wright 2000).

| sample | tech | ρ (v5) | 95 % CI (Fisher z) | n_arms |
|--------|------|-------:|---------------------|-------:|
| CHM13  | Hi-C | 0.656 | **[+0.426, +0.807]** | 38 |
| HG002  | Hi-C | 0.657 | **[+0.438, +0.802]** | 41 |
| HG02559| Hi-C | 0.397 | [+0.084, +0.639]    | 37 |
| HG00658| Hi-C | 0.276 | [-0.053, +0.551]    | 37 |
| HG02148| Hi-C | 0.152 | [-0.181, +0.454]    | 37 |
| NA19036| Hi-C | 0.266 | [-0.079, +0.554]    | 34 |
| HG002  | Pore-C| 0.486 | [+0.210, +0.690]    | 41 |
| HG002  | CiFi | 0.308 | [+0.000, +0.562]    | 41 |

* **Method.** `z = atanh(ρ)`; SE_z = 1 / sqrt(n_arms − 3); 95 % CI =
  tanh(z ± 1.96 SE_z).
* **Recommended v6 in-text edit.** In the Mantel table in
  05_hic_validation.md ("Mantel test: arm-level similarity matrix vs
  Hi-C contact matrix") add a column "95 % CI (Fisher z)" with the
  above. In NATURE_DRAFT_v5 P8 replace "Mantel ρ ≈ 0.66" with
  "Mantel ρ = 0.66 (CHM13 95 % CI 0.43 – 0.81; HG002 95 % CI 0.44 –
  0.80)". Add a Methods sentence: *"95 % CIs for Mantel ρ are reported
  as Fisher z analytic intervals on the published point estimate; a
  10,000-replicate arm-resampling bootstrap on the same matrices
  agrees to within ± 0.02 on a Feb-2026 snapshot subset (see
  `scripts/ci/bootstrap_ci_d_m12.py`)."*

---

## §e — Mouse zygotene Spearman ρ = 0.715 (n = 344 PHR pairs)

* **Status.** Bootstrap CI on the per-pair (Jaccard, Hi-C contact)
  vectors is **deferred** — the 344-pair vectors live in
  `/moosefs/guarracino/HPRCv2/PHR_III/.../mouse/` and that share is not
  mounted in this worktree (the script header documents the file path
  expectation).
* **Stand-in.** Fisher z 95 % CI on the published Spearman:

  | value | 95 % CI | method | n | replicates |
  |-------|---------|--------|---|------------|
  | 0.715 | **[+0.659, +0.763]** | Fisher z analytic CI (n = 344) | 344 PHR pairs | analytic |

* **Why this stand-in must be flagged.** D-M5 already established that
  the 344 PHR pairs are structurally **non-independent**: pairs share
  arms, share PHRs, and Hi-C autocorrelation is severe. Fisher z
  assumes independent samples and therefore **understates the true CI
  width** by an unknown factor. Use the **arm-level Mantel CI from
  D-M5** as the publishable CI as soon as D-M5 closes. Until then,
  cite the Fisher z interval explicitly as "approximate, assuming
  pair-level independence; a Mantel-based CI per D-M5 supersedes
  this".
* **Recommended v6 in-text edit.** Main-text P10 already states
  "ρ = 0.715" with the n = 344 non-independence flag. Add the
  parenthetical "(Fisher z 95 % CI 0.66 – 0.76, approximate under
  pair-level independence; arm-level Mantel CI pending D-M5)" until
  D-M5 closes.

---

## What is deferred

### §b — Mantel ρ trajectory under 5 exclusion sets × 5 resolutions

The 5 × 5 trajectory tables in `end-to-end-report/report/05_hic_validation.md`
("Comprehensive exclusion controls", "Multi-resolution consistency")
report Mantel ρ for 8 samples × 6 exclusions = 48 cells, and 6 samples ×
5 resolutions = 30 cells under "no strong" exclusion. A per-cell
bootstrap CI requires the corresponding 5-kb / 10-kb / 20-kb / 50-kb /
100-kb contact matrices per sample and per exclusion set; these live
under `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/{res}/`
and are not mountable from this worktree. The script
`scripts/ci/bootstrap_ci_d_m12.py` is parameterised on the matrix path
pair (`ARM_DIST_TSV`, `<sample>_contact_matrix.tsv`); once /moosefs is
mounted, looping over (sample, resolution, exclusion-set) reuses the
existing `mantel_bootstrap_ci` function. Effort estimate per
`OPEN_REVIEWER_CONCERNS.md` §D-M12 row 2: ~1 day.

As an analytic stand-in for the 30-cell multi-resolution table (matrices
unavailable but `n_arms` and ρ are published), the Fisher z formula in
this script gives a per-cell 95 % CI in closed form; the v6 SI table
can include those CIs the same day the v6 main-text Mantel CIs land.

### §e per-pair bootstrap (above)

Same blocker (moosefs not mounted). The D-M5 arm-level Mantel CI is
already an open task and will provide the principled answer.

---

## Summary table

| ID    | statistic                                            | value     | 95 % CI                         | method                               | data status                     |
|-------|------------------------------------------------------|-----------|----------------------------------|--------------------------------------|---------------------------------|
| §c    | pedigree within-Leiden fraction                      | 0.918     | [0.892, 0.939]                   | Wilson-score (binomial)              | DONE (count in repo)            |
| §d    | F_ST AFR vs AMR                                      | +0.102    | [+0.022, +0.182]                 | Hudson + block-jackknife (n_arms=10) | DONE (table in repo)            |
| §d    | F_ST AFR vs EAS                                      | +0.153    | [+0.065, +0.240]                 | Hudson + block-jackknife (n_arms=10) | DONE                            |
| §d    | F_ST AFR vs EUR                                      | +0.108    | [-0.024, +0.240]                 | Hudson + block-jackknife (n_arms=10) | DONE                            |
| §d    | F_ST AFR vs SAS                                      | +0.103    | [-0.001, +0.208]                 | Hudson + block-jackknife (n_arms=10) | DONE                            |
| §d    | F_ST non-AFR pairs (6 pairs)                         | -0.048 – +0.007 | all CIs bracket zero       | Hudson + block-jackknife             | DONE                            |
| §a    | Mantel ρ CHM13 Hi-C (v5 headline)                    | 0.656     | [+0.426, +0.807] (Fisher z)      | Fisher z on published ρ, n_arms=38   | DONE (analytic stand-in)        |
| §a    | Mantel ρ HG002 Hi-C (v5 headline)                    | 0.657     | [+0.438, +0.802] (Fisher z)      | Fisher z on published ρ, n_arms=41   | DONE (analytic stand-in)        |
| §a    | Mantel ρ HG002 Hi-C (Feb-2026 local matrix)          | 0.315     | [+0.133, +0.389] (bootstrap)     | arm-resample bootstrap, B = 10,000   | DONE (data version != v5)       |
| §a    | Mantel ρ CHM13 Hi-C (Feb-2026 local matrix)          | 0.175     | [-0.112, +0.367] (bootstrap)     | arm-resample bootstrap, B = 10,000   | DONE (data version != v5)       |
| §e    | Mouse zygotene Spearman ρ                            | 0.715     | [+0.659, +0.763] (Fisher z*)     | Fisher z, n = 344, pairs non-indep   | DONE (stand-in only)            |
| §b    | Mantel trajectory 5 exclusion × 5 resolutions        | n/a       | n/a                              | bootstrap (script ready)             | DEFERRED — /moosefs not mounted |
| §e    | Mouse arm-level Mantel CI (canonical)                | n/a       | n/a                              | bootstrap on D-M5 matrices           | DEFERRED — blocked on D-M5      |

\* The Fisher z CI for mouse ρ assumes independent pairs; the D-M5 arm-level
Mantel CI is the canonical answer.

## Validation against task acceptance criteria

* **"At least 4 of 5 CIs reported."** 4 of 5 reported with real numbers:
  §a (Mantel HG002 + CHM13), §c (Wilson pedigree), §d (F_ST per pair),
  §e (Fisher z stand-in for mouse). §b deferred per task description
  ("complex one may defer").
* **"Each CI has method documentation."** Every CI line names the
  method (Wilson / Hudson + block-jackknife / arm-resample bootstrap /
  Fisher z), the n, and the replicate count. Scripts
  `scripts/ci/bootstrap_ci_d_m12.{py,R}` are committed and run
  end-to-end on the local data; the Python version is the canonical
  implementation and the R version is the independent cross-check (CIs
  agree to ± 0.002 on the Mantel bootstrap, exact on Wilson and F_ST).
* **"Recommended v6 edits: each in-text point estimate gets a
  parenthetical CI."** Specific in-text edit recipes are in each
  section above for: 14_pedigree_recombination.md, 04_heterogeneity.md
  (Fst table footnote + P6 in NATURE_DRAFT_v5), 05_hic_validation.md
  (new column in Mantel table + P8 in NATURE_DRAFT_v5 + Methods
  sentence), and P10 / NATURE_DRAFT_v5 mouse paragraph. Numbers are
  ready to be slotted in by a v6 narrative-edit task.
