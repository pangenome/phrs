# R phyper() Parameter Modifications for Copy-Number Weighted ORA: Technical Analysis

## Research Question

**How can R's `phyper(q, m, n, k)` parameters be modified to handle copy-number weights for gene instance counting rather than unique gene counting?**

## Core Parameter Transformation

### Standard phyper() Usage
```r
phyper(q, m, n, k, lower.tail = FALSE)
```

Where:
- `q` = observed overlap count - 1 (for P(X ≥ q+1))
- `m` = items of interest in population (pathway genes)
- `n` = items not of interest in population (non-pathway genes)
- `k` = sample size (query genes)

### Weighted phyper() Modification
```r
phyper(q_weighted, m_weighted, n_weighted, k_weighted, lower.tail = FALSE)
```

Where each parameter transforms from gene counts to gene instance counts based on copy numbers.

## Detailed Parameter Mapping

### 1. Query Sample Size: k → k_weighted

**Standard Gene Counting:**
```r
k_standard <- length(unique(query_genes))
# Counts each gene symbol once, regardless of copy number
```

**Copy-Number Weighted Counting:**
```r
k_weighted <- sum(query_df$copy_number)
# Counts each gene copy as separate instance
# Example: Gene with 14 copies contributes 14 to k_weighted
```

**Mathematical Relationship:**
```
k_weighted = Σ(copy_number_i) for all i in query
k_standard = |{gene_i : gene_i in query}|

Ratio: k_weighted/k_standard = average copy number in query
```

### 2. Pathway Population Size: m → m_weighted

**Standard Gene Counting:**
```r
pathway_genes_in_background <- intersect(pathway_genes, background_genes)
m_standard <- length(pathway_genes_in_background)
# Counts unique pathway genes present in background
```

**Copy-Number Weighted Counting:**
```r
pathway_bg_df <- background_df[background_df$gene %in% pathway_genes, ]
m_weighted <- sum(pathway_bg_df$copy_number)
# Sums copy numbers for all pathway genes in background
```

**Mathematical Relationship:**
```
m_weighted = Σ(copy_number_j) for all j in (pathway ∩ background)
m_standard = |{gene_j : gene_j in pathway AND gene_j in background}|
```

### 3. Non-Pathway Population Size: n → n_weighted

**Standard Gene Counting:**
```r
n_standard <- length(background_genes) - m_standard
# Total background genes minus pathway genes
```

**Copy-Number Weighted Counting:**
```r
n_weighted <- sum(background_df$copy_number) - m_weighted
# Total background instances minus pathway instances
```

**Mathematical Relationship:**
```
n_weighted = Σ(copy_number_all) - m_weighted
n_standard = |background| - m_standard

Total population: m_weighted + n_weighted = total instances in background
```

### 4. Observed Overlap: q → q_weighted

**Standard Gene Counting:**
```r
overlap_genes <- intersect(query_genes, pathway_genes)
q_standard <- length(overlap_genes)
# Counts unique genes in both query and pathway
```

**Copy-Number Weighted Counting:**
```r
query_pathway_df <- query_df[query_df$gene %in% pathway_genes, ]
q_weighted <- sum(query_pathway_df$copy_number)
# Sums copy numbers for genes in both query and pathway
```

**Mathematical Relationship:**
```
q_weighted = Σ(copy_number_k) for all k in (query ∩ pathway)
q_standard = |{gene_k : gene_k in query AND gene_k in pathway}|
```

## Parameter Validation and Constraints

### Hypergeometric Distribution Constraints

The modified parameters must satisfy standard hypergeometric constraints:

```r
# Constraint validation for weighted parameters
validate_weighted_params <- function(q_weighted, m_weighted, n_weighted, k_weighted) {
  # 1. Non-negative integers
  all_params <- c(q_weighted, m_weighted, n_weighted, k_weighted)
  stopifnot(all(all_params >= 0))
  stopifnot(all(all_params == floor(all_params)))  # integers
  
  # 2. Logical constraints
  stopifnot(q_weighted <= k_weighted)              # overlap ≤ sample size
  stopifnot(q_weighted <= m_weighted)              # overlap ≤ population of interest
  stopifnot(k_weighted <= m_weighted + n_weighted) # sample ≤ total population
  
  # 3. Feasibility constraint
  min_possible_overlap <- max(0, k_weighted - n_weighted)
  max_possible_overlap <- min(k_weighted, m_weighted)
  stopifnot(q_weighted >= min_possible_overlap)
  stopifnot(q_weighted <= max_possible_overlap)
  
  return(TRUE)
}
```

### Copy Number Model Assumptions

**Critical Assumption:** Copy numbers in background model should reflect the same biological context as query copy numbers.

```r
# Ensure consistent copy number models
validate_copy_model_consistency <- function(query_df, background_df) {
  # 1. Same genes should have similar copy number distributions
  common_genes <- intersect(query_df$gene, background_df$gene)
  
  if (length(common_genes) > 5) {  # Sufficient overlap for comparison
    query_copies <- query_df[query_df$gene %in% common_genes, ]
    bg_copies <- background_df[background_df$gene %in% common_genes, ]
    
    # Merge and compare
    merged <- merge(query_copies, bg_copies, by = "gene", suffixes = c("_query", "_bg"))
    correlation <- cor(merged$copy_number_query, merged$copy_number_bg)
    
    if (correlation < 0.5) {
      warning("Low correlation between query and background copy numbers")
    }
  }
  
  # 2. Check for extreme outliers that might skew results
  max_query_copies <- max(query_df$copy_number)
  max_bg_copies <- max(background_df$copy_number)
  
  if (max_query_copies > 10 * max_bg_copies || max_bg_copies > 10 * max_query_copies) {
    warning("Extreme copy number differences between query and background")
  }
}
```

## Concrete Implementation Example

### Using PHR Dataset Parameters

Based on actual PHR data from dependency analysis:

```r
# PHR dataset characteristics (from gene_copy_summary.csv)
phr_genes <- 35       # unique protein-coding genes
phr_instances <- 1189 # total gene instances (copy-weighted)

# Example: Olfactory receptor pathway testing
olfactory_genes <- c("OR4F17", "OR4F29", "OR4F3", "OR4F5")

# Standard approach parameters
k_standard <- 35                    # PHR unique genes
q_standard <- 4                     # OR genes in PHR (unique)
m_standard <- 400                   # OR genes in genome (unique)
n_standard <- 20000 - 400          # Non-OR genes in genome

# Calculate standard p-value
pval_standard <- phyper(q_standard-1, m_standard, n_standard, k_standard, 
                       lower.tail = FALSE)

# Weighted approach parameters (based on copy numbers)
k_weighted <- 1189                  # PHR total instances
q_weighted <- 4 * 14               # OR genes in PHR (4 genes × 14 copies each)
m_weighted <- 400 * 2              # OR instances in genome (est. 2x copies avg)
n_weighted <- 2000000 - 800        # Non-OR instances in genome

# Calculate weighted p-value
pval_weighted <- phyper(q_weighted-1, m_weighted, n_weighted, k_weighted, 
                       lower.tail = FALSE)

# Compare results
print(paste("Standard p-value:", format(pval_standard, scientific = TRUE)))
print(paste("Weighted p-value:", format(pval_weighted, scientific = TRUE)))
print(paste("Fold change:", pval_weighted / pval_standard))
```

## Mathematical Equivalence Verification

### Proof of Equivalence with Instance Expansion

**Theorem:** Parameter weighting method produces identical results to instance expansion method.

**Instance Expansion Method:**
```r
# Create expanded vectors
query_expanded <- rep(query_df$gene, query_df$copy_number)
background_expanded <- rep(background_df$gene, background_df$copy_number)

# Apply standard hypergeometric test to expanded data
q_exp <- sum(query_expanded %in% pathway_genes)
m_exp <- sum(background_expanded %in% pathway_genes)
n_exp <- length(background_expanded) - m_exp
k_exp <- length(query_expanded)

pval_expansion <- phyper(q_exp-1, m_exp, n_exp, k_exp, lower.tail = FALSE)
```

**Parameter Weighting Method:**
```r
# Calculate weighted parameters directly
q_weight <- sum(query_df[query_df$gene %in% pathway_genes, "copy_number"])
m_weight <- sum(background_df[background_df$gene %in% pathway_genes, "copy_number"])
n_weight <- sum(background_df$copy_number) - m_weight
k_weight <- sum(query_df$copy_number)

pval_weighted <- phyper(q_weight-1, m_weight, n_weight, k_weight, lower.tail = FALSE)
```

**Proof:**
1. `k_weight == k_exp` (both count total query instances)
2. `q_weight == q_exp` (both count overlap instances)  
3. `m_weight == m_exp` (both count pathway instances in background)
4. `n_weight == n_exp` (both count non-pathway instances)

Therefore: `pval_weighted == pval_expansion` (within numerical precision)

### Empirical Verification

```r
# Test mathematical equivalence empirically
test_equivalence <- function(query_df, pathway_genes, background_df, tolerance = 1e-12) {
  # Method 1: Instance expansion
  query_expanded <- rep(query_df$gene, query_df$copy_number)
  background_expanded <- rep(background_df$gene, background_df$copy_number)
  
  q_exp <- sum(query_expanded %in% pathway_genes)
  m_exp <- sum(background_expanded %in% pathway_genes)
  n_exp <- length(background_expanded) - m_exp
  k_exp <- length(query_expanded)
  
  pval_expansion <- phyper(q_exp-1, m_exp, n_exp, k_exp, lower.tail = FALSE)
  
  # Method 2: Parameter weighting
  q_weight <- sum(query_df[query_df$gene %in% pathway_genes, "copy_number"])
  m_weight <- sum(background_df[background_df$gene %in% pathway_genes, "copy_number"])
  n_weight <- sum(background_df$copy_number) - m_weight
  k_weight <- sum(query_df$copy_number)
  
  pval_weighted <- phyper(q_weight-1, m_weight, n_weight, k_weight, lower.tail = FALSE)
  
  # Verify parameter equivalence
  params_match <- all(c(
    q_weight == q_exp,
    m_weight == m_exp,
    n_weight == n_exp,
    k_weight == k_exp
  ))
  
  # Verify p-value equivalence
  pval_diff <- abs(pval_weighted - pval_expansion)
  pvals_match <- pval_diff < tolerance
  
  return(list(
    parameters_match = params_match,
    pvalues_match = pvals_match,
    pvalue_difference = pval_diff,
    expansion_pval = pval_expansion,
    weighted_pval = pval_weighted
  ))
}
```

## Edge Cases and Numerical Considerations

### 1. Zero Copy Numbers

**Issue:** Genes with 0 copies in background (annotation errors)

```r
# Handle zero-copy genes
filter_zero_copies <- function(df) {
  zero_copy_genes <- df[df$copy_number == 0, "gene"]
  if (length(zero_copy_genes) > 0) {
    warning(paste("Removing", length(zero_copy_genes), "genes with zero copies"))
    df <- df[df$copy_number > 0, ]
  }
  return(df)
}
```

### 2. Large Copy Numbers

**Issue:** Extremely large copy numbers (>500) may cause numerical overflow

```r
# Cap extreme copy numbers
cap_extreme_copies <- function(df, max_copies = 500) {
  extreme_genes <- df[df$copy_number > max_copies, ]
  if (nrow(extreme_genes) > 0) {
    warning(paste("Capping", nrow(extreme_genes), "genes at", max_copies, "copies"))
    df$copy_number[df$copy_number > max_copies] <- max_copies
  }
  return(df)
}
```

### 3. Small Sample Sizes

**Issue:** Very small query sets may lack statistical power

```r
# Check minimum sample size requirements
validate_sample_size <- function(k_weighted, min_instances = 10) {
  if (k_weighted < min_instances) {
    warning(paste("Small sample size:", k_weighted, "instances. Results may be unreliable."))
  }
  
  # For hypergeometric tests, we also need sufficient background
  # Rule of thumb: background should be at least 10x larger than sample
  return(k_weighted >= min_instances)
}
```

## Performance Optimization

### Memory Efficiency

**Advantage of Parameter Weighting:**

```r
# Memory usage comparison
memory_usage_comparison <- function(query_df, background_df) {
  # Instance expansion memory requirement
  expansion_memory <- (sum(query_df$copy_number) + sum(background_df$copy_number)) * 8  # bytes per string
  
  # Parameter weighting memory requirement  
  weighting_memory <- (nrow(query_df) + nrow(background_df)) * 16  # bytes per row (gene + copy_number)
  
  reduction_factor <- expansion_memory / weighting_memory
  
  return(list(
    expansion_mb = expansion_memory / (1024^2),
    weighting_mb = weighting_memory / (1024^2),
    reduction_factor = reduction_factor
  ))
}
```

### Computational Speed

**Benchmark Results (estimated for PHR-scale data):**

| Method | Memory (MB) | Runtime (ms) | Scalability |
|--------|-------------|--------------|-------------|
| Instance expansion | ~50 | ~100 | O(total_copies) |
| Parameter weighting | ~1 | ~5 | O(unique_genes) |
| **Improvement** | **50x less** | **20x faster** | **Much better** |

## Conclusion

R's `phyper()` function can be effectively modified for copy-number weighted ORA by:

1. **Parameter Transformation:** Convert all parameters from gene counts to gene instance counts
2. **Mathematical Equivalence:** Results are identical to instance expansion but computationally superior  
3. **Statistical Validity:** Standard hypergeometric theory applies after parameter transformation
4. **Practical Benefits:** Significant improvements in memory usage and computational speed

**Key Implementation Points:**
- Carefully validate parameter constraints
- Handle edge cases (zero copies, extreme values)
- Ensure copy number model consistency between query and background
- Provide comprehensive error checking and warnings

This approach enables efficient, statistically rigorous copy-number weighted pathway enrichment analysis using standard R statistical functions.
