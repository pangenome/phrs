# Statistical Validation Report: Copy-Number Weighted phyper() Modifications

**Date:** 2026-04-01 (updated by verification task)
**Task:** validate-weighted-phyper (verified by .verify-validate-weighted-phyper)
**Objective:** Validate statistical properties of copy-number weighted hypergeometric test modifications

## Executive Summary

This report presents comprehensive statistical validation of copy-number weighted modifications to R's `phyper()` hypergeometric test function for Over-Representation Analysis (ORA). The validation covers null distributions, Type I error rates, power analysis, ROC analysis, FDR correction behavior, positive controls, copy-number stratified analyses, and different background models.

### Key Findings

- **Mathematical Correctness**: 100% equivalence with instance expansion approach (200/200 test cases, all parameters match exactly)
- **Anti-Conservative Behavior**: Under gene-level sampling, weighted phyper() has inflated Type I error (~0.22 at alpha=0.05) because the hypergeometric model assumes instance-level independence while gene selection brings all copies as a cluster
- **FDR Compromised**: BH correction does not control FDR for weighted p-values (actual FDR ~0.66 vs target 0.05) due to anti-conservative null p-values
- **Standard phyper() Well-Calibrated**: Unweighted test controls Type I error at all alpha levels and across all background models

## Validation Methods and Results

### 1. Null Distribution Validation (KS Test for Uniformity)

Under gene-level null (random gene sampling, 1000 simulations):

| Method | KS p-value | Mean p-value | Type I rate (α=0.05) | Uniform? |
|--------|-----------|--------------|---------------------|----------|
| Weighted phyper | ~0 | 0.558 | 0.217 | NO |
| Standard phyper | ~0 | 0.599 | 0.024 | NO (conservative) |

Both methods fail the KS uniformity test due to discreteness of the hypergeometric distribution. However:
- Standard phyper is conservative (Type I < α), which is expected and safe
- Weighted phyper is **anti-conservative** (Type I >> α), which is problematic

### 2. Type I Error Rate Control

Type I error rates at multiple alpha levels (2000 null simulations):

| Alpha | Weighted | Standard |
|-------|----------|----------|
| 0.01 | 0.142 | 0.004 |
| 0.05 | 0.225 | 0.018 |
| 0.10 | 0.263 | 0.071 |
| 0.20 | 0.321 | 0.181 |

**At alpha=0.05**: Weighted 0.225 (95% CI: 0.207–0.244), Standard 0.018 (95% CI: 0.013–0.025).

Weighted phyper() does **not** control Type I error at any conventional alpha level.

#### Copy Number Magnitude Effect

| Scenario | Type I (weighted) | Type I (standard) | Mean CN |
|----------|------------------|-------------------|---------|
| All CN=1 | 0.019 | 0.016 | 1.0 |
| All CN=5 | 0.153 | 0.022 | 5.0 |
| CN 1–3 | 0.104 | 0.019 | 2.0 |
| CN 1–8 | 0.198 | 0.021 | 4.7 |
| CN 1–20 | 0.287 | 0.017 | 10.5 |

When all copy numbers = 1, weighted test is properly calibrated (reduces to standard test). Inflation increases monotonically with copy number magnitude.

### 3. Power Analysis vs Standard Approaches

Power comparison (500 simulations per enrichment level, alpha = 0.05):

| Enrichment Factor | Power (weighted) | Power (standard) | Ratio |
|-------------------|-----------------|-----------------|-------|
| 1.0 (null) | 0.210 | 0.024 | 8.75 |
| 1.5 | 0.480 | 0.170 | 2.82 |
| 2.0 | 0.728 | 0.414 | 1.76 |
| 3.0 | 0.954 | 0.860 | 1.11 |
| 5.0 | 1.000 | 0.998 | 1.00 |

The apparent power advantage of the weighted test at low enrichment levels is an artifact of its anti-conservative behavior. At high enrichment (≥3x), both methods approach 100% power.

### 4. ROC Analysis for Power Comparison

ROC analysis with 500 null + 500 enriched (3x) samples:

| Metric | Weighted | Standard |
|--------|----------|----------|
| AUC | 0.970 | 0.978 |

True positive rate at controlled false positive rates:

| FPR | TPR (weighted) | TPR (standard) |
|-----|---------------|----------------|
| 0.01 | 0.626 | 0.740 |
| 0.05 | 0.818 | 0.850 |
| 0.10 | 0.920 | 0.930 |
| 0.20 | 0.970 | 0.970 |

Standard phyper() achieves **higher** TPR at every controlled FPR level, demonstrating better calibrated sensitivity. The weighted test's ROC is distorted by its anti-conservative p-value distribution.

### 5. FDR Correction Validation

Testing 100 pathways (10 truly enriched), BH correction at FDR < 0.05:

| Method | Mean discoveries | True positives | False positives | Actual FDR |
|--------|-----------------|----------------|-----------------|------------|
| Weighted | 22.4 | 7.4 | 15.0 | **0.661** |
| Standard | 0.3 | 0.3 | 0.1 | 0.026 |

**Critical finding**: BH correction fails to control FDR for weighted phyper() because the procedure assumes valid (uniform or conservative) null p-values. The inflated null p-values lead to an actual FDR of 66% vs the target 5%.

### 6. Positive Controls (Known Enrichment)

Both methods detect strong enrichment:

| Scenario | p (weighted) | p (standard) | Both significant? |
|----------|-------------|-------------|-------------------|
| Strong enrichment, small pathway | 1.28e-37 | 2.41e-08 | Yes |
| Moderate enrichment, small pathway | 1.99e-17 | 2.29e-03 | Yes |
| Strong enrichment, medium pathway | 1.40e-50 | 1.60e-08 | Yes |
| Moderate enrichment, medium pathway | 6.21e-07 | 6.86e-03 | Yes |
| Strong enrichment, large pathway | 1.38e-24 | 2.77e-05 | Yes |
| Moderate enrichment, large pathway | 0.096 | 0.156 | No |

Both methods successfully detect strongly enriched pathways. Weighted p-values are consistently much smaller due to effective sample size inflation.

### 7. Copy-Number Stratified Analyses

Type I error by copy-number stratum (bimodal background: low CN 1–3, high CN 10–30):

| Query stratum | Type I (weighted) | Type I (standard) | Mean CN |
|---------------|-------------------|-------------------|---------|
| Low CN | 0.188 | 0.016 | 2.0 |
| High CN | 0.320 | 0.014 | 19.8 |
| Mixed | 0.284 | 0.024 | 10.9 |

Anti-conservative behavior is substantially worse for high-CN strata, where the effective sample size inflation is greatest.

### 8. Different Background Models

| Background Model | Mean CN | Var(CN) | Type I (weighted) | Type I (standard) |
|-----------------|---------|---------|-------------------|-------------------|
| Diploid (all CN=2) | 2.0 | 0.0 | 0.070 | 0.022 |
| Uniform CN 1–5 | 2.9 | 2.0 | 0.172 | 0.028 |
| Geometric CN | 3.5 | 9.4 | 0.164 | 0.030 |
| Bimodal CN | 8.4 | 52.6 | 0.308 | 0.018 |
| Heavy-tail CN | 5.5 | 24.0 | 0.210 | 0.012 |

Anti-conservative behavior correlates with copy number magnitude. Even the diploid model (all CN=2) shows mild inflation (0.070 vs 0.022). Standard phyper is well-calibrated across all background models.

### 9. Mathematical Equivalence Confirmation

200 random test cases with varying copy numbers:
- P-value equivalence: **100%**
- Parameter q match: **100%**
- Parameter m match: **100%**
- Parameter n match: **100%**
- Parameter k match: **100%**

The weighted parameter calculation is mathematically identical to full instance expansion.

## Root Cause Analysis

The anti-conservative behavior has a clear mathematical explanation:

1. **Model assumption**: The hypergeometric test (whether via parameter weighting or instance expansion) models draws as independent samples from the population of gene *instances*
2. **Actual sampling**: In ORA, genes are selected as units, bringing all their copies together. This creates **clustering** — instances of the same gene are perfectly correlated
3. **Effect**: Clustering inflates the effective sample size relative to what the hypergeometric model expects, reducing relative variance and producing systematically small p-values
4. **Scaling**: Even with uniform copy numbers, scaling all parameters by constant c produces a hypergeometric distribution with smaller relative variance than the original, so p-values shrink

This is analogous to the "design effect" in survey statistics, where cluster sampling produces larger variance than simple random sampling, making tests that assume independence anti-conservative.

## Recommendations

### When Weighted phyper() Is Appropriate

The weighted hypergeometric test is statistically valid when the **instance-level independence assumption** holds — i.e., when individual gene copies can independently be part of the query set. This may apply when:

- Individual genomic copies at different locations independently overlap with regions of interest (e.g., PHR regions)
- Selection is truly at the instance level, not the gene level
- The analysis specifically concerns gene dosage effects

### When to Prefer Standard phyper()

Use unweighted hypergeometric when:
- Genes are selected as units (e.g., differential expression analysis)
- Copy number information is primarily used for annotation, not selection
- Type I error control is critical
- FDR guarantees are needed

### Mitigation Strategies

For cases where copy-number weighting is desired but gene-level selection applies:

1. **Permutation-based p-values**: Permute gene labels (preserving copy number structure) and compute empirical p-values. This correctly accounts for clustering.
2. **Report both methods**: Always provide standard and weighted results for comparison.
3. **Conservative correction**: Apply additional Bonferroni-type correction to compensate for anti-conservative behavior.
4. **Calibration**: Use simulation under the actual null model to establish critical thresholds.

## Technical Implementation

### Core Function

```r
weighted_hypergeometric_test <- function(query_df, pathway_genes, background_df) {
  valid_query_genes <- intersect(query_df$gene, background_df$gene)
  query_filtered <- query_df[query_df$gene %in% valid_query_genes, ]
  query_with_bg_copies <- merge(query_filtered[, "gene", drop = FALSE],
                               background_df, by = "gene", all.x = TRUE)

  query_in_pathway <- query_with_bg_copies[query_with_bg_copies$gene %in% pathway_genes, ]
  pathway_in_background <- background_df[background_df$gene %in% pathway_genes, ]

  k_weighted <- sum(query_with_bg_copies$copy_number)
  q_weighted <- sum(query_in_pathway$copy_number)
  m_weighted <- sum(pathway_in_background$copy_number)
  n_weighted <- sum(background_df$copy_number) - m_weighted

  pvalue <- phyper(q_weighted-1, m_weighted, n_weighted, k_weighted, lower.tail=FALSE)
  fold_enrichment <- (q_weighted/k_weighted) / (m_weighted/(m_weighted+n_weighted))

  list(pvalue = pvalue, fold_enrichment = fold_enrichment,
       overlap_instances = q_weighted, query_instances = k_weighted,
       pathway_instances = m_weighted, background_instances = m_weighted + n_weighted)
}
```

### Statistical Framework

**Transformation**: Gene counts → Gene instance counts
- k_weighted = Σ(copy_number) for query genes
- q_weighted = Σ(copy_number) for query ∩ pathway genes
- m_weighted = Σ(copy_number) for pathway genes in background
- n_weighted = total background instances − m_weighted

**Test**: `phyper(q_weighted−1, m_weighted, n_weighted, k_weighted, lower.tail=FALSE)`

## Artifacts

| File | Description |
|------|-------------|
| `comprehensive_statistical_validation.R` | Full validation suite covering all 9 test categories |
| `debug_weighted_phyper.R` | Core implementation and equivalence testing |
| `theoretical_validation.R` | Mathematical property validation |
| `corrected_statistical_validation.R` | Earlier corrected null distribution tests |
| `phr_specific_validation.R` | PHR-specific scenario testing |
| `statistical_best_practices_weighted_ora.md` | Usage guidelines |
| `comprehensive_validation_results.RData` | Saved R validation results |

## Conclusion

The copy-number weighted phyper() is **mathematically correct** — it is exactly equivalent to instance expansion. However, it has **significant statistical limitations** under gene-level sampling:

- Type I error is inflated 4–6× above the nominal level
- FDR correction is ineffective (actual FDR ~66% at target 5%)
- The standard (unweighted) hypergeometric outperforms the weighted version at every controlled FPR threshold

**Recommendation**: Use weighted phyper() only when instance-level independence is justified. For gene-level analyses, prefer standard phyper() or permutation-based approaches. Always report both weighted and standard results.

---

*Statistical validation completed 2026-04-01. Updated with comprehensive verification results by .verify-validate-weighted-phyper.*
