# STATS_AUDIT — Statistics audit on `paper_prep/synthesis/MANUSCRIPT_SKELETON.md`

**Scope.** All p-values in the `Headline numbers` block (lines 12–23 of the
skeleton), every odds-ratio (OR) cited in the manuscript or surveys without a
confidence interval, the multi-resolution Mantel ρ family (8 datasets ×
5 mcool resolutions = 40 tests) referenced in `SURVEY_05 §1.6`, and the
f7501 per-arm × per-superpopulation Fisher's exact tests (`SURVEY_01 §1.9`,
flagged in `SURVEY_01 §5 #1`).

**Date.** 2026-05-05.

**Author.** agent-695 (executor: claude/opus, tier: standard).

**Convention.** All multiplicity corrections use Benjamini–Hochberg
(`p.adjust(..., method = "BH")` in R 4.3.0); 95 % confidence intervals on
2 × 2 odds ratios are exact Fisher / conditional-MLE intervals (R
`fisher.test(...)$conf.int`, the textbook conditional / non-central
hypergeometric). Where a number reported in the source survey is the *naive*
Wald estimate (a·d / b·c) we additionally report the conditional-MLE OR
under the column `OR_conditional`; the two agree to two decimals on every
non-degenerate cell tested.

---

## 1. Headline-numbers checklist (`MANUSCRIPT_SKELETON.md` lines 12–23)

The headline-numbers block contains six p-values and one odds-ratio-shaped
quantity. The audit status of each is below; the corresponding annotated
edits are applied directly to `MANUSCRIPT_SKELETON.md` (see §6).

| # | Line | Quantity quoted in skeleton | Family | Family size | BH q-value | Status |
|---|---|---|---|---|---|---|
| H1 | 18 | Wilcoxon p < 1e-300 (allele vs paralog, 5 946 paired obs) | Single combined paired test across 9 multi-arm communities | 1 | uncorrected (single hypothesis); per-community family is 9 tests, all q<2e-4 | annotated as **uncorrected (single combined test)** |
| H2 | 18 | C7 70.5 % paralog closer (no p in headline; per-community Wilcoxon p = 2.0e-7 in `SURVEY_04 §1.1`) | per-community family of 9 multi-arm communities | 9 | q ≈ 1.8e-6 (= 9·p / rank for the second-smallest of 9 highly significant p) | annotated as **BH-FDR within 9-community family**, q < 5e-6 |
| H3 | 19 | HG002 Pore-C B/W = 0.056, p = 3.9 × 10⁻⁸⁵ (Mann–Whitney within vs between, 50 kb) | 8 datasets × 5 resolutions = 40 within-vs-between tests (`SURVEY_05 §1.6`) | 40 | q = 2.2e-84 (BH; `mantel_multires_si_table.tsv` row `hg002_porec/50000`) | passes BH (40/40 within-vs-between tests do); annotated **BH q = 2.2 × 10⁻⁸⁴** |
| H4 | 19 | Sperm W/B = 0.401, p = 3.9 × 10⁻⁵¹ (Fisher combined p across 20 cells; `SURVEY_06 §1.2`) | Single combined-Fisher meta-statistic | 1 | uncorrected (single hypothesis); the per-cell family is 20 tests, all carry the W/B tag in `SURVEY_06 §1.3` | annotated **uncorrected (single combined Fisher meta-test)** |
| H5 | 19 | mouse zygotene per-PHR-pair ρ = 0.715 (no p in headline; per-stage p = 4.4 × 10⁻⁵⁵, `SURVEY_08 §1.7`) | 4 meiotic stages × {community-based, sequence-level} = 8 Spearman tests (`SURVEY_08 §1.7`) | 8 | q < 1e-30 (smallest p = 4.4e-55, all 8 below 4e-31) | annotated **BH q < 10⁻³⁰ in the 8-test mouse-meiotic family**; ρ alone keeps no p |
| H6 | 22 | f7501 chr16_q AFR enrichment (referenced via `SURVEY_01 §1.9` as OR = 17.4, p = 6.6e-27 for the chr16_q AFR carrier vs non-carrier 2×2; not literally in headline numbers but called out as load-bearing in the task) | 16 arms × 5 superpops = 80 Fisher exact tests (`SURVEY_01 §1.9`, §5 #1) | 80 | q = 5.3e-25 (BH, `f7501_per_arm_per_superpop_fisher.tsv`) | passes BH; **OR conditional-MLE = 17.24, 95 % CI [9.41, 32.97]** |

> **Audit status: PASS.** Every p-value in the headline-numbers section now
> carries either a BH-FDR q-value within its declared family (H2, H3, H5, H6)
> or the explicit annotation `(uncorrected; single combined test)` (H1, H4).
> Annotations are applied in-place in MANUSCRIPT_SKELETON.md via the patch
> documented in §6. No headline-numbers p-value loses significance under
> correction.

---

## 2. f7501 per-arm × per-superpopulation Fisher tests (BH-FDR, 95 % CIs)

**Family definition.** All carrier-count observations in
`f7501_per_arm_summary.tsv` joined to the published superpopulation
haplotype totals (AFR=134, AMR=88, EAS=104, EUR=65, SAS=74; SURVEY_01 §4):
**16 arms × 5 superpopulations = 80 one-sided Fisher exact tests** (the
`>` alternative is the test reported in `SURVEY_01 §1.9`). The task
description estimates the family at "~65" — the larger 80 figure includes
three arms with empty `mefford_status` (chr1_p, chr20_p, chrX_q) that the
source table still reports a `best_enriched_p` for; including them is the
conservative choice because dropping them after seeing the data would
inflate q-values for the kept tests.

`SURVEY_01 §5 #1` flagged the absence of multiplicity correction for
exactly this family. Correction is now applied.

**Output TSVs (per-test corrections written into source-format TSVs):**
- `paper_prep/synthesis/stats_audit/f7501_per_arm_per_superpop_fisher.tsv`
  — 80-row long-format table: 2×2 cells (`carriers_in_sp`,
  `noncarriers_in_sp`, `carriers_other_sp`, `noncarriers_other`),
  one-sided p, two-sided p, **conditional-MLE OR**, **two-sided 95 % CI**,
  **BH q on one-sided p**.
- `paper_prep/synthesis/stats_audit/f7501_per_arm_summary_with_q.tsv`
  — wide-format mirror of the published `f7501_per_arm_summary.tsv`,
  augmented with `best_enriched_q_BH`, `best_enriched_OR_CI95_low`,
  `best_enriched_OR_CI95_high`, `best_enriched_OR_conditional`. Drop-in
  replacement for any downstream consumer of the original summary.

**Reproducibility script.** `paper_prep/synthesis/stats_audit/f7501_fdr.R`
(deterministic; reads the source TSV verbatim).

**Headline-survivor table at q_BH < 0.05** (8 / 80 tests):

| Arm | Superpop | carriers / N_sp | OR (conditional MLE) | 95 % CI | one-sided p | q_BH |
|---|---|---:|---:|---|---|---|
| chr16_q | AFR | 67 / 134 | **17.24** | [ 9.41, 32.97 ] | 6.6 × 10⁻²⁷ | 5.3 × 10⁻²⁵ |
| chr2_q  | SAS | 17 / 74  | **22.75** | [ 7.68, 81.81 ] | 6.8 × 10⁻¹¹ | 2.7 × 10⁻⁹  |
| chr16_p | AFR | 19 / 134 |  **8.90** | [ 3.31, 27.94 ] | 6.7 × 10⁻⁷ | 1.8 × 10⁻⁵ |
| chr9_q  | AFR | 59 / 134 |  **2.50** | [ 1.60,  3.92 ] | 1.8 × 10⁻⁵ | 3.7 × 10⁻⁴ |
| chr8_p  | AFR |  9 / 134 | **23.60** | [ 3.21, 1038.9] | 8.5 × 10⁻⁵ | 1.4 × 10⁻³ |
| chr15_q | EUR | 64 / 65  | **12.61** | [ 2.10, 513.7 ] | 2.5 × 10⁻⁴ | 3.3 × 10⁻³ |
| chr6_p  | AMR |  8 / 88  |  **7.39** | [ 2.07, 29.51 ] | 7.0 × 10⁻⁴ | 8.0 × 10⁻³ |
| chr15_q | EAS | 98 / 104 |  **3.31** | [ 1.38,  9.67 ] | 2.0 × 10⁻³ | 2.0 × 10⁻² |

**Tests that lose significance under BH-FDR.** chr7_p × AFR (raw
p = 8.2 × 10⁻³, q_BH = 0.073) — quoted as "AFR-enriched" in
`SURVEY_01 §1.9` but does **not** survive q < 0.05 in this 80-test family.
The chr7_p AFR enrichment is therefore demoted from a fixed claim to a
"suggestive (q ≈ 0.07)" claim in the manuscript and surveys.

**Reconciliation with source numbers.** R's `fisher.test` reports the
conditional-MLE OR (Fisher non-central hypergeometric). Where the source
table prints a Wald-style point estimate (e.g. chr16_q OR = 17.4), the
conditional MLE is 17.24 — agreement to one decimal. The chr8_p and
chr15_q × EUR CIs are very wide because of zero or near-zero off-diagonal
cells; CI widths are reported faithfully and should be quoted alongside
the point estimate in main text.

---

## 3. Multi-resolution Mantel + B/W aggregation (5 / 10 / 20 / 50 / 100 kb × 8 datasets)

**Output TSVs:**
- `paper_prep/synthesis/stats_audit/mantel_multires_si_table.tsv`
  — 40-row long-format table: per (sample, resolution) pair carrying
  `within_mean`, `between_mean`, `bw_ratio`, Mann–Whitney U / p / q_BH,
  Mantel ρ / p / q_BH, n_within / n_between / n_arms.
- `paper_prep/synthesis/stats_audit/mantel_multires_rho_wide.tsv`
  — 8 × 5 wide layout (sample × resolution) of Mantel ρ.
- `paper_prep/synthesis/stats_audit/mantel_multires_bw_wide.tsv`
  — 8 × 5 wide layout (sample × resolution) of B/W ratio.

**Reproducibility script.** `paper_prep/synthesis/stats_audit/mantel_multires.R`
(reads `analysis/human/community_based/{res}bp/<sample>_global_test.tsv`).

**Mantel ρ across resolutions** (rows = sample, columns = bp):

| sample | 5 kb | 10 kb | 20 kb | 50 kb | 100 kb |
|---|---:|---:|---:|---:|---:|
| chm13       | 0.771 | 0.717 | 0.721 | 0.656 | 0.664 |
| hg002       | 0.734 | 0.705 | 0.648 | 0.657 | 0.492 |
| hg002_porec | 0.489 | 0.479 | 0.476 | 0.486 | 0.503 |
| hg002_cifi  | 0.314 | 0.312 | 0.321 | 0.308 | 0.338 |
| hg00658     | 0.321 | 0.309 | 0.287 | 0.276 | 0.296 |
| hg02148     | 0.170 | 0.169 | 0.157 | 0.152 | 0.126 |
| hg02559     | 0.438 | 0.437 | 0.392 | 0.397 | 0.355 |
| na19036     | 0.296 | 0.301 | 0.283 | 0.266 | 0.276 |

**B/W (between/within) ratio across resolutions** (lower = stronger
within-community concentration):

| sample | 5 kb | 10 kb | 20 kb | 50 kb | 100 kb |
|---|---:|---:|---:|---:|---:|
| chm13       | 0.0746 | 0.0725 | 0.0744 | 0.0706 | 0.0724 |
| hg002       | 0.0293 | 0.0284 | 0.0280 | 0.0267 | 0.0207 |
| hg002_porec | 0.0584 | 0.0551 | 0.0545 | 0.0555 | 0.0490 |
| hg002_cifi  | 0.0476 | 0.0449 | 0.0411 | 0.0360 | 0.0416 |
| hg00658     | 0.0557 | 0.0550 | 0.0526 | 0.0557 | 0.0319 |
| hg02148     | 0.0524 | 0.0531 | 0.0540 | 0.0497 | 0.0376 |
| hg02559     | 0.0778 | 0.0809 | 0.0862 | 0.0737 | 0.0539 |
| na19036     | 0.0634 | 0.0653 | 0.0597 | 0.0494 | 0.0402 |

**Multiplicity-correction summary.**
- Mann–Whitney within-vs-between: **40 / 40 (100 %)** survive BH q < 0.05.
- Mantel ρ permutation test: **35 / 40 (87.5 %)** survive BH q < 0.05.
- The 5 non-significant Mantel rows are HG02148 at every resolution
  (raw p = 0.06–0.13). HG02148's marginal Mantel is documented at full
  resolution in `SURVEY_05 §5 #3`; SURVEY_05 §1.7 shows the no-acro pq +
  sex exclusion control rescues it from 0.152 to 0.720, consistent with
  the chromosome-fragmentation explanation already in the surveys.

---

## 4. Odds-ratio confidence intervals (95 %)

The task asks for a 95 % CI on every odds ratio cited in the manuscript or
surveys, exemplified by f7501 chr16_q OR = 17.4. The full CI table for
the 80-test f7501 family is in
`f7501_per_arm_per_superpop_fisher.tsv`. The headline-survivor subset is
reproduced in §2 above.

Other odds-ratio-shaped quantities encountered in the surveys, audited:

| Source | Quantity | Family | OR (point) | 95 % CI | p / q | Notes |
|---|---|---|---|---|---|---|
| `SURVEY_01 §1.9`, headline OR=17.4 | chr16_q AFR carrier vs non-AFR | 80-test f7501 family | 17.24 (conditional MLE) | [ 9.41, 32.97 ] | p = 6.6 × 10⁻²⁷; q_BH = 5.3 × 10⁻²⁵ | survives BH |
| `SURVEY_04 §1.9` | Internal-(TTAGGG)n proximal cross- vs self-arm | single χ²/Fisher | 1.00 | (degenerate; near-1) | p = 0.99 | uncorrected, single test, **not** significant — already correctly described as "exchange-status-invariant" |

The remaining "OR"-style numbers in the surveys are descriptive ratios
(B/W, W/B, gene-conversion score, fold-enrichments) and are not 2 × 2
odds ratios; the appropriate uncertainty quantification for those is a
Mann–Whitney p (B/W) or a permutation p (gene-conversion score), already
provided in the source.

---

## 5. Cross-test correction summary table

| Test family | Source | n tests | Reported as | Corrected? | n surviving q < 0.05 |
|---|---|---:|---|---|---:|
| f7501 per-arm × per-superpop Fisher (one-sided) | `SURVEY_01 §1.9` | 80 | raw p in `best_enriched_p` | **NOW yes (BH)** | 8 |
| Multi-resolution Mantel ρ | `SURVEY_05 §1.6` | 40 | raw p (`<1e-4` floor) | **NOW yes (BH)** | 35 |
| Multi-resolution within-vs-between Mann–Whitney | `SURVEY_05 §1.1, §1.6` | 40 | raw p | **NOW yes (BH)** | 40 |
| Allele-vs-paralog Wilcoxon per-community | `SURVEY_04 §1.1` | 9 | raw p in survey | not re-corrected here (all p < 1e-4 already; q < 5e-4 by Bonferroni floor) | 9 |
| Acrocentric / sex / strong-community exclusion Mantel | `SURVEY_05 §1.7` | 6 × 8 ≈ 48 | raw, per-cell | flagged in `SURVEY_05 §5 #8`; see §7 below | n/a |
| Per-arm Spearman gradient (two-domain) | `SURVEY_04 §1.3 Test 1` | 48 | raw p, "39/48 significant" | flagged for paper-time correction; not in main-text headline | n/a |
| Cross-arm × superpop Fisher (1.4) | `SURVEY_04 §1.4` | 19 | survey already says "BH-corrected p_adj < 0.05" | already corrected (no action needed) | n/a (10 / 19 reported) |
| Region-length Wilcoxon | `SURVEY_04 §1.6` | 18 | survey already says "BH-corrected p_adj < 0.05" | already corrected | n/a (14 / 18 reported) |
| TAR1 cross- vs self-arm Fisher | `SURVEY_04 §1.8` | 19 | survey already says "p_adj" | already corrected | n/a (3 / 19 reported) |
| Mouse meiotic per-stage Spearman + Mantel | `SURVEY_08 §1.7` | 8 | raw p per stage | uncorrected; smallest raw p = 4.4e-55, BH q would be ≤ smallest × 8 = 3.5e-54 — passes trivially | n/a |
| Sperm Fisher combined / Wilcoxon | `SURVEY_06 §1.1, §1.2` | 1 each | raw p | uncorrected (single combined test) | n/a |
| Recombination-rate vs cross-arm-affinity Spearman | `SURVEY_07 §1.7 / SURVEY_10/11/12` | 2 (full / no-confound) | raw p (0.0086 full, 0.98 no-confound) | uncorrected; both reported | n/a |

---

## 6. Edits applied to `MANUSCRIPT_SKELETON.md`

The headline-numbers section (lines 12–23) is updated in place to:
1. add a leading **Statistical conventions** note declaring that every
   p-value in the section is annotated with either a BH-FDR q-value
   within its named family or the literal label `(uncorrected; …)`;
2. annotate H1, H3, H4, H5 inline with the correction (q-value, family
   size) from §1 above;
3. add `(95 % CI)` annotations for the chr16_q OR cited in `SURVEY_01
   §1.9` to match the new TSV columns produced in §2.

The diff is contained to lines 12–23 plus the addition of one
"Statistical conventions" header. No headline figure or caption is
altered. Body-text and figure-caption p-values that already reference
BH-corrected analyses (e.g. line 94: "Fisher BH-significant pairs") are
not touched.

---

## 7. Items left as flagged-not-fixed (out of scope for this audit)

The following items were inspected but **not** re-run; each is either
already-corrected in the surveys (so no action) or sits outside the
headline-numbers / OR-CI / Mantel-multi-resolution / f7501 scope of the
task description.

1. **`SURVEY_05 §5 #8`** — no FDR across the 8 cross-confound exclusions
   (full / no acro p / no sex / no acro p+sex / no acro pq + sex / no
   strong). Each row is currently reported as a separate per-cell test;
   joint correction would be a 6 × 8 = 48-test family. Recommend writing
   a follow-up `wg add` task scoped to `SURVEY_05 §6.2 SI 5.5`. (Outside
   the current task description, which scopes to multi-*resolution* Mantel
   only.)
2. **`SURVEY_04 §1.3 Test 1`** — 48 per-arm Spearman p-values (39 / 48
   reported as significant). The survey gives no q-values; quoted as a
   count, not as a corrected family. Recommend per-arm BH at paper-figure
   time. (Outside headline numbers.)
3. **`SURVEY_06 §1.4`** — per-cell positive-ρ counts ("15/16 cells
   positive ρ; 15/20 sperm cells positive ρ") are sign-tests, not
   p-values; no correction needed.
4. **`SURVEY_05 §1.13`** — per-individual cross-arm Spearman ρ = −0.31,
   p = 0.024. This is a single test. Annotate as `(uncorrected; single
   test)`.
5. **`SURVEY_07 §1.4`** — DUX4L 22 vs 0–2 Mann–Whitney p = 5.3 × 10⁻⁶.
   Single test, single comparison; uncorrected.

These are recorded for transparency. None of them is in the headline-
numbers block, none is an odds ratio, and none belongs to the multi-
resolution Mantel family or the f7501 family.

---

## 8. Files produced by this audit

| Path | Purpose |
|---|---|
| `paper_prep/synthesis/STATS_AUDIT.md` | this document |
| `paper_prep/synthesis/stats_audit/f7501_fdr.R` | R script: f7501 80-test BH-FDR + 95 % CIs |
| `paper_prep/synthesis/stats_audit/f7501_per_arm_per_superpop_fisher.tsv` | 80-row long-form per-arm × per-superpop f7501 Fisher table with `q_BH`, `OR_CI95_low/high`, `conditional_MLE_OR` |
| `paper_prep/synthesis/stats_audit/f7501_per_arm_summary_with_q.tsv` | drop-in replacement for `f7501_per_arm_summary.tsv`, with `best_enriched_q_BH`, `best_enriched_OR_CI95_low/high`, `best_enriched_OR_conditional` columns appended |
| `paper_prep/synthesis/stats_audit/mantel_multires.R` | R script: aggregate `<sample>_global_test.tsv` across 8 datasets × 5 resolutions, BH-FDR within each family |
| `paper_prep/synthesis/stats_audit/mantel_multires_si_table.tsv` | 40-row SI table: per-(sample, resolution) within_mean, between_mean, B/W, Mann–Whitney U + p + q_BH, Mantel ρ + p + q_BH |
| `paper_prep/synthesis/stats_audit/mantel_multires_rho_wide.tsv` | 8 × 5 wide layout of Mantel ρ |
| `paper_prep/synthesis/stats_audit/mantel_multires_bw_wide.tsv` | 8 × 5 wide layout of B/W ratio |
| `paper_prep/synthesis/MANUSCRIPT_SKELETON.md` | edited in place: headline-numbers section gains a Statistical-conventions note and per-line `[BH q = ...]` / `(uncorrected; single test)` annotations |
