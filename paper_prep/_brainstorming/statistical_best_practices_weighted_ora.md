# Statistical Best Practices for Copy-Number Weighted ORA

## Critical Statistical Caveat

**The weighted hypergeometric test is anti-conservative under gene-level sampling.** When genes are selected as units (the typical ORA scenario), Type I error is inflated ~4–6× above the nominal level. This occurs because the hypergeometric model assumes instance-level independence, but gene selection brings all copies as a cluster.

See the full validation report (`weighted_phyper_statistical_validation_report.md`) for details.

## When to Use Each Method

### Standard (Unweighted) Hypergeometric — Preferred Default
- Gene-level selection (differentially expressed genes, genes in a region, etc.)
- When Type I error control and FDR guarantees are important
- When copy numbers are used for annotation, not as part of the test

### Weighted Hypergeometric — Use With Caution
- Instance-level selection (individual genomic copies independently selected)
- Exploratory analysis where sensitivity matters more than specificity
- **Always pair with standard results for comparison**

### Permutation-Based — Gold Standard
- When copy number effects matter AND gene-level selection applies
- Permute gene labels, preserving copy number structure
- Computationally expensive but correctly accounts for clustering

## Implementation

```r
weighted_hypergeometric_test <- function(query_df, pathway_genes, background_df) {
  # Ensure query genes are subset of background
  valid_query_genes <- intersect(query_df$gene, background_df$gene)
  query_filtered <- query_df[query_df$gene %in% valid_query_genes, ]

  # Use background copy numbers for consistency
  query_with_bg_copies <- merge(query_filtered[, "gene", drop = FALSE],
                               background_df, by = "gene", all.x = TRUE)

  # Calculate weighted parameters
  query_in_pathway <- query_with_bg_copies[query_with_bg_copies$gene %in% pathway_genes, ]
  pathway_in_background <- background_df[background_df$gene %in% pathway_genes, ]

  k_weighted <- sum(query_with_bg_copies$copy_number)
  q_weighted <- sum(query_in_pathway$copy_number)
  m_weighted <- sum(pathway_in_background$copy_number)
  n_weighted <- sum(background_df$copy_number) - m_weighted

  # Validate constraints
  stopifnot(q_weighted <= k_weighted)
  stopifnot(q_weighted <= m_weighted)
  stopifnot(k_weighted <= (m_weighted + n_weighted))

  # Hypergeometric test
  pvalue <- phyper(q_weighted-1, m_weighted, n_weighted, k_weighted, lower.tail = FALSE)

  list(
    pvalue = pvalue,
    overlap_instances = q_weighted,
    query_instances = k_weighted,
    pathway_instances = m_weighted,
    background_instances = m_weighted + n_weighted,
    fold_enrichment = (q_weighted/k_weighted) / (m_weighted/(m_weighted+n_weighted))
  )
}
```

## Multiple Testing Correction

**For standard phyper**: Use BH (FDR) correction as usual.
```r
standard_pvalues_fdr <- p.adjust(standard_pvalues, method = "fdr")
```

**For weighted phyper**: BH correction is **unreliable** because the null p-values are not uniform. Options:
1. Permutation-based FDR (preferred)
2. Very conservative Bonferroni correction
3. Calibration via simulation under the actual null model

## Validation Summary

| Test | Standard phyper | Weighted phyper |
|------|----------------|-----------------|
| Type I error control | 0.018–0.024 | 0.15–0.31 |
| FDR control (BH at 5%) | 0.026 actual | 0.661 actual |
| ROC AUC (3x enrichment) | 0.978 | 0.970 |
| Positive control detection | Yes | Yes |
| Mathematical correctness | N/A | 100% equiv. |

## Reporting Guidelines

**Required**:
1. State whether gene-level or instance-level selection applies
2. Report both standard and weighted results when copy numbers vary
3. Document the copy number source and methodology
4. Apply and document the multiple testing correction method
5. Acknowledge the anti-conservative behavior of weighted tests

**Recommended**:
1. Include sensitivity analysis with different copy number thresholds
2. Compare standard vs weighted results explicitly
3. Use permutation-based p-values for formal inference
4. Discuss biological relevance of copy number effects for the pathways tested

## Quality Control Checklist

- [ ] Selection mechanism identified (gene-level vs instance-level)
- [ ] Copy number data consistent between query and background
- [ ] Extreme copy number outliers handled appropriately
- [ ] Background gene set represents appropriate universe
- [ ] Multiple testing correction appropriate for the p-value type
- [ ] Both standard and weighted results reported when different
- [ ] Anti-conservative behavior acknowledged if using weighted test
- [ ] Results interpreted in biological context

---

*Updated 2026-04-01 based on comprehensive validation findings (see weighted_phyper_statistical_validation_report.md)*
