# R phyper() Modifications for Copy-Number-Weighted ORA: Technical Research

## Executive Summary

This research investigates how to modify R's `phyper()` hypergeometric test function to handle copy-number-weighted gene instances for Over-Representation Analysis (ORA). The key insight is that copy-number weighting transforms the standard hypergeometric model from unique gene counts to gene instance counts, requiring careful parameter recalculation while preserving statistical validity.

**Core Finding:** Copy-number weighted ORA can be implemented using standard `phyper()` by transforming input parameters from gene counts to instance counts, maintaining mathematical equivalence with instance expansion approaches while improving computational efficiency.

## Background: Standard vs Weighted Hypergeometric Testing

### Standard ORA with phyper()

R's `phyper(q, m, n, k, lower.tail=FALSE)` implements the hypergeometric distribution:

```r
# Standard hypergeometric test parameters:
# q = observed overlap - 1 (for lower.tail=FALSE)
# m = genes in pathway (white balls in urn)
# n = genes not in pathway (black balls in urn) 
# k = genes in query set (balls drawn)

pvalue <- phyper(overlap-1, pathway_size, background_size-pathway_size, query_size, lower.tail=FALSE)
```

**Mathematical Model:**
```
P(X ≥ k) = Σ [C(m,i) × C(n,k-i)] / C(m+n,k)  for i = k to min(m,k)

Where C(n,r) = n! / (r! × (n-r)!)
```

### Copy-Number Weighted ORA Modification

Copy-number weighting transforms the model from "genes" to "gene instances":

```r
# Weighted hypergeometric parameters:
# q_weighted = observed instance overlap - 1
# m_weighted = pathway instances in background  
# n_weighted = non-pathway instances in background
# k_weighted = query instances (total copy number)

pvalue <- phyper(weighted_overlap-1, weighted_pathway_size, 
                 weighted_background_size-weighted_pathway_size, 
                 weighted_query_size, lower.tail=FALSE)
```

**Key Transformation:** Each parameter scales from gene counts to instance counts based on copy numbers.

## Parameter Mapping Methodology

### 1. Query Size (k → k_weighted)

**Standard:** Count unique genes in query set
**Weighted:** Sum of copy numbers for all genes in query

```r
# Standard approach
k_standard <- length(query_genes)

# Weighted approach  
k_weighted <- sum(query_copy_numbers)

# Example: PHR dataset
k_standard <- 35      # unique genes
k_weighted <- 1189    # total instances
```

### 2. Pathway Size in Background (m → m_weighted)

**Standard:** Count unique pathway genes in background
**Weighted:** Sum copy numbers for pathway genes in background

```r
# Standard approach
pathway_genes_in_bg <- intersect(pathway_genes, background_genes)
m_standard <- length(pathway_genes_in_bg)

# Weighted approach
pathway_in_bg_df <- background_df[background_df$gene %in% pathway_genes, ]
m_weighted <- sum(pathway_in_bg_df$copy_number)

# Example: Olfactory receptors
m_standard <- 400     # unique OR genes in genome
m_weighted <- 800     # OR gene instances (avg 2x copies)
```

### 3. Non-Pathway Background Size (n → n_weighted)

**Standard:** Background genes not in pathway
**Weighted:** Background instances not in pathway

```r
# Standard approach
n_standard <- length(background_genes) - m_standard

# Weighted approach  
n_weighted <- sum(background_df$copy_number) - m_weighted

# Example: Human genome
n_standard <- 20000 - 400 = 19600    # non-OR genes
n_weighted <- 2000000 - 800 = 1999200 # non-OR instances
```

### 4. Observed Overlap (q → q_weighted)

**Standard:** Count unique genes in query ∩ pathway
**Weighted:** Sum copy numbers for genes in query ∩ pathway

```r
# Standard approach
overlap_genes <- intersect(query_genes, pathway_genes)
q_standard <- length(overlap_genes)

# Weighted approach
query_pathway_genes <- query_df[query_df$gene %in% pathway_genes, ]
q_weighted <- sum(query_pathway_genes$copy_number)

# Example: PHR olfactory genes
q_standard <- 4      # OR4F17, OR4F29, OR4F3, OR4F5
q_weighted <- 56     # 4 genes × 14 copies each
```

## Mathematical Equivalence Analysis

### Instance Expansion vs Parameter Weighting

**Hypothesis:** Weighted parameter approach is mathematically equivalent to instance expansion.

**Instance Expansion Method:**
```r
# Expand all datasets by copy number
query_expanded <- rep(query_df$gene, query_df$copy_number)
background_expanded <- rep(background_df$gene, background_df$copy_number)

# Standard phyper on expanded sets
q_exp <- sum(query_expanded %in% pathway_genes)
m_exp <- sum(background_expanded %in% pathway_genes) 
n_exp <- length(background_expanded) - m_exp
k_exp <- length(query_expanded)

pval_expansion <- phyper(q_exp-1, m_exp, n_exp, k_exp, lower.tail=FALSE)
```

**Parameter Weighting Method:**
```r
# Calculate weighted parameters directly
q_weight <- sum(query_df[query_df$gene %in% pathway_genes, "copy_number"])
m_weight <- sum(background_df[background_df$gene %in% pathway_genes, "copy_number"])
n_weight <- sum(background_df$copy_number) - m_weight
k_weight <- sum(query_df$copy_number)

pval_weighted <- phyper(q_weight-1, m_weight, n_weight, k_weight, lower.tail=FALSE)
```

**Mathematical Proof:** Both methods generate identical hypergeometric parameters:
- Both count the same total instances
- Both preserve the same overlap structure
- Both maintain the same background composition

**Verification:** `pval_expansion == pval_weighted` (within numerical precision)

## Computational Efficiency Implications

### Memory Usage Comparison

| Method | Memory Complexity | PHR Example |
|--------|------------------|-------------|
| Instance expansion | O(total_instances) | ~1.2K elements |
| Parameter weighting | O(unique_genes) | ~35 elements |
| **Reduction** | **97% less memory** | **34x smaller** |

### Runtime Performance

**Instance Expansion:**
```r
# Time complexity: O(copy_expansion + hypergeometric_calculation)
system.time({
  query_expanded <- rep(genes, copies)        # O(total_copies)
  result <- phyper(...)                       # O(1)
})
```

**Parameter Weighting:**
```r  
# Time complexity: O(aggregation + hypergeometric_calculation)
system.time({
  weighted_params <- calculate_weights(...)   # O(unique_genes)
  result <- phyper(...)                       # O(1)
})
```

**Scaling Analysis:**
- Small datasets: Minimal difference
- Large datasets (>50K instances): Parameter weighting significantly faster
- Memory-constrained environments: Parameter weighting essential

## Statistical Validation Framework

### Null Distribution Properties

**Requirement:** Under null hypothesis (no true enrichment), p-values should be uniformly distributed.

**Test Implementation:**
```r
validate_null_distribution <- function(n_simulations = 1000) {
  null_pvalues <- replicate(n_simulations, {
    # Generate random query with no true enrichment
    null_query <- sample(background_genes, size = query_size)
    
    # Calculate weighted parameters
    weighted_params <- calculate_weighted_params(null_query, pathway, background)
    
    # Run weighted phyper
    phyper(weighted_params$q-1, weighted_params$m, 
           weighted_params$n, weighted_params$k, lower.tail=FALSE)
  })
  
  # Test uniformity
  ks_result <- ks.test(null_pvalues, punif)
  return(ks_result$p.value > 0.05)  # Should NOT reject uniformity
}
```

### Type I Error Control

**Requirement:** False positive rate should equal nominal α level.

**Test Implementation:**
```r
test_type_i_error <- function(alpha = 0.05, n_tests = 10000) {
  false_positives <- sum(replicate(n_tests, {
    null_result <- run_weighted_ora(random_query, pathway, background)
    return(null_result$pvalue < alpha)
  }))
  
  observed_rate <- false_positives / n_tests
  expected_rate <- alpha
  
  # Should be approximately equal
  return(abs(observed_rate - expected_rate) < 0.01)
}
```

## Implementation Recommendations

### 1. Recommended R Implementation

```r
weighted_hypergeometric_test <- function(query_df, pathway_genes, background_df) {
  # Input validation
  stopifnot(all(c("gene", "copy_number") %in% names(query_df)))
  stopifnot(all(c("gene", "copy_number") %in% names(background_df)))
  
  # Calculate weighted parameters
  query_in_pathway <- query_df[query_df$gene %in% pathway_genes, ]
  pathway_in_background <- background_df[background_df$gene %in% pathway_genes, ]
  
  k_weighted <- sum(query_df$copy_number)                    # query instances
  q_weighted <- sum(query_in_pathway$copy_number)            # overlap instances  
  m_weighted <- sum(pathway_in_background$copy_number)       # pathway instances
  n_weighted <- sum(background_df$copy_number) - m_weighted  # non-pathway instances
  
  # Parameter validation
  stopifnot(q_weighted <= k_weighted)
  stopifnot(q_weighted <= m_weighted)
  stopifnot(k_weighted <= (m_weighted + n_weighted))
  
  # Hypergeometric test
  pvalue <- phyper(q_weighted-1, m_weighted, n_weighted, k_weighted, 
                   lower.tail = FALSE)
  
  # Return comprehensive results
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

## Best Practices and Recommendations

### 1. Statistical Guidelines

1. **Always validate null distributions** for your specific background model
2. **Use appropriate multiple testing correction** (FDR recommended)
3. **Report both standard and weighted results** when they differ substantially
4. **Include confidence intervals** for fold-enrichment estimates
5. **Validate with positive controls** known to be enriched

### 2. Implementation Guidelines

1. **Prefer parameter weighting over instance expansion** for computational efficiency
2. **Implement comprehensive input validation** for data quality
3. **Use vectorized operations** for multiple pathway testing
4. **Handle edge cases gracefully** with informative warnings
5. **Provide both gene-level and instance-level result interpretations**

## Conclusion

Copy-number weighted ORA using modified R `phyper()` parameters provides a statistically sound and computationally efficient approach to account for gene dosage effects in pathway enrichment analysis. The key insights are:

1. **Mathematical Equivalence:** Parameter weighting is mathematically equivalent to instance expansion but computationally superior
2. **Statistical Validity:** Standard hypergeometric theory applies after parameter transformation
3. **Practical Implementation:** Standard R functions can be used with careful parameter calculation
4. **Performance Benefits:** Significant memory and runtime improvements over naive approaches

The approach is particularly valuable for datasets with extreme copy number variation (like PHRs) where gene dosage effects may substantially influence biological pathways. This research provides the foundation for robust, efficient implementation of copy-number weighted ORA in R while maintaining statistical rigor and computational practicality.
