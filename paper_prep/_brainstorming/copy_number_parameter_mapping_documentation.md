# Copy-Number Parameter Mapping Documentation

## Overview

This document provides comprehensive documentation for the copy-number-weighted hypergeometric parameter mapping implementation in `copy_number_phyper_mapping.R`. This implementation addresses the failed `map-copy-number` task by providing concrete parameter transformation from copy-number weighted gene counts to R phyper(q,m,n,k) parameters.

## Mathematical Foundation

### Standard Hypergeometric Model

The standard hypergeometric distribution models sampling without replacement:

```
P(X = k | N, K, n) = C(K,k) × C(N-K,n-k) / C(N,n)
```

**Standard phyper() parameters:**
- `q`: observed overlap count - 1 (for P(X ≥ q+1))
- `m`: items of interest in population (pathway genes)  
- `n`: items not of interest in population (non-pathway genes)
- `k`: sample size (query genes)

### Copy-Number-Weighted Model

Copy-number weighting transforms the model from gene counts to gene instance counts:

```
P(X = k_w | N_w, K_w, n_w) = C(K_w,k_w) × C(N_w-K_w,n_w-k_w) / C(N_w,n_w)
```

**Weighted phyper() parameters:**
- `q_weighted`: observed instance overlap - 1
- `m_weighted`: pathway instances in background
- `n_weighted`: non-pathway instances in background  
- `k_weighted`: query instances (total copy number)

## Parameter Transformation Logic

### 1. Query Size: k_standard → k_weighted

**Transformation:**
```
k_standard = |{unique genes in query}|
k_weighted = Σ(copy_number_i) for all i in query
```

**Example:**
```r
# Query with 3 genes
query <- data.frame(
  gene_name = c("GENE1", "GENE2", "GENE3"),
  copy_number = c(5, 14, 3)
)
k_standard = 3      # unique genes
k_weighted = 22     # 5 + 14 + 3 instances
```

### 2. Overlap: q_standard → q_weighted

**Transformation:**
```
q_standard = |{genes in query ∩ pathway}|
q_weighted = Σ(copy_number_j) for all j in (query ∩ pathway)
```

**Example:**
```r
# Pathway genes: OR4F17, OR4F29, OR4F3
# Query overlap: OR4F17 (14 copies), OR4F29 (14 copies)
q_standard = 2      # unique overlap genes
q_weighted = 28     # 14 + 14 instances
```

### 3. Pathway Size: m_standard → m_weighted

**Transformation:**
```
m_standard = |{pathway genes in background}|
m_weighted = Σ(copy_number_k) for all k in (pathway ∩ background)
```

**Example:**
```r
# Pathway has 400 olfactory genes in genome
# Each has ~2 copies on average
m_standard = 400    # unique pathway genes
m_weighted = 800    # pathway instances
```

### 4. Non-Pathway Background: n_standard → n_weighted

**Transformation:**
```
n_standard = |{non-pathway genes in background}|
n_weighted = Σ(all background copy numbers) - m_weighted
```

**Example:**
```r
# Background: 20,000 genes, 2M total instances
# Pathway: 800 instances
n_standard = 19600  # non-pathway genes
n_weighted = 1999200 # non-pathway instances
```

## Function Reference

### calculate_weighted_phyper_params()

**Purpose:** Calculate copy-number-weighted hypergeometric parameters

**Parameters:**
- `query_df`: Data frame with `gene_name` and `copy_number` columns
- `pathway_genes`: Character vector of pathway gene names
- `background_df`: Data frame with `gene_name` and `copy_number` columns
- `validate_params`: Logical, whether to validate constraints

**Returns:**
- `k_weighted`: Query instances
- `q_weighted`: Overlap instances
- `m_weighted`: Pathway instances in background
- `n_weighted`: Non-pathway instances in background
- `parameters_standard`: Standard parameters for comparison
- `fold_enrichment_weighted`: Copy-weighted enrichment ratio
- `validation`: Parameter validation results
- `metadata`: Additional analysis information

**Example Usage:**
```r
# Load required function
source("copy_number_phyper_mapping.R")

# Prepare data
query_df <- data.frame(
  gene_name = c("OR4F17", "OR4F29", "GENE1"),
  copy_number = c(14, 14, 5)
)

pathway_genes <- c("OR4F17", "OR4F29", "OR4F3", "OR4F5")

background_df <- data.frame(
  gene_name = paste0("GENE", 1:1000),
  copy_number = rep(2, 1000)
)

# Calculate parameters
params <- calculate_weighted_phyper_params(
  query_df, pathway_genes, background_df
)

# View results
print(params$k_weighted)     # 33 (14+14+5)
print(params$q_weighted)     # 28 (14+14)
print(params$m_weighted)     # Pathway instances in background
print(params$n_weighted)     # Non-pathway instances
```

### run_weighted_hypergeometric_test()

**Purpose:** Complete copy-number-weighted hypergeometric test

**Parameters:**
- `query_df`: Query data with gene names and copy numbers
- `pathway_genes`: Pathway gene names
- `background_df`: Background data with gene names and copy numbers
- `alpha`: Significance level (default 0.05)

**Returns:**
- `pvalue`: Hypergeometric test p-value
- `significant`: Logical, whether result is significant
- `observed_overlap_weighted`: Weighted overlap count
- `expected_overlap_weighted`: Expected overlap under null
- `fold_enrichment`: Enrichment ratio
- `method`: Test method identifier
- `parameters`: Full parameter calculation results
- `warnings`: Any validation warnings

**Example Usage:**
```r
# Run complete test
result <- run_weighted_hypergeometric_test(
  query_df, pathway_genes, background_df
)

print(paste("P-value:", format(result$pvalue, scientific = TRUE)))
print(paste("Significant:", result$significant))
print(paste("Fold enrichment:", round(result$fold_enrichment, 2)))
```

### verify_equivalence_with_expansion()

**Purpose:** Validate mathematical equivalence with instance expansion method

**Parameters:**
- `query_df`: Query data
- `pathway_genes`: Pathway genes
- `background_df`: Background data  
- `tolerance`: Numerical tolerance for comparison

**Returns:**
- `parameters_equivalent`: Whether parameters match exactly
- `pvalues_equivalent`: Whether p-values match within tolerance
- `pvalue_difference`: Absolute difference between p-values
- `weighted_pvalue`: P-value from parameter weighting
- `expansion_pvalue`: P-value from instance expansion
- `memory_reduction_factor`: Memory savings factor

**Example Usage:**
```r
# Verify equivalence
equiv <- verify_equivalence_with_expansion(
  query_df, pathway_genes, background_df
)

print(paste("Parameters equivalent:", equiv$parameters_equivalent))
print(paste("P-values equivalent:", equiv$pvalues_equivalent))
print(paste("Memory reduction:", equiv$memory_reduction_factor, "x"))
```

## Validation Examples

### Example 1: PHR Olfactory Receptor Analysis

```r
# Load PHR data (if available)
if (file.exists("gene_copy_summary.csv")) {
  copy_data <- read.csv("gene_copy_summary.csv")
  phr_query <- copy_data[copy_data$gene_biotype == "protein_coding", 
                        c("gene_name", "total_copies")]
  names(phr_query)[2] <- "copy_number"
  
  # Test olfactory receptor pathway
  or_pathway <- c("OR4F17", "OR4F29", "OR4F3", "OR4F5")
  
  # Create mock background
  background <- data.frame(
    gene_name = paste0("GENE", 1:20000),
    copy_number = sample(1:5, 20000, replace = TRUE, 
                        prob = c(0.7, 0.15, 0.1, 0.03, 0.02))
  )
  
  # Add OR genes to background
  or_bg <- data.frame(
    gene_name = or_pathway,
    copy_number = c(2, 2, 3, 1)
  )
  background <- rbind(background, or_bg)
  
  # Run test
  result <- run_weighted_hypergeometric_test(
    phr_query, or_pathway, background
  )
  
  print(result)
}
```

### Example 2: Parameter Validation

```r
# Test with edge cases
test_edge_cases <- function() {
  
  # Case 1: No overlap
  query_no_overlap <- data.frame(
    gene_name = c("GENE1", "GENE2"),
    copy_number = c(5, 3)
  )
  pathway_no_overlap <- c("GENE3", "GENE4")
  background <- data.frame(
    gene_name = paste0("GENE", 1:100),
    copy_number = rep(2, 100)
  )
  
  params_no_overlap <- calculate_weighted_phyper_params(
    query_no_overlap, pathway_no_overlap, background
  )
  print("No overlap case:")
  print(params_no_overlap$q_weighted)  # Should be 0
  
  # Case 2: Complete overlap
  query_complete <- data.frame(
    gene_name = c("GENE1", "GENE2"),
    copy_number = c(5, 3)
  )
  pathway_complete <- c("GENE1", "GENE2")
  
  params_complete <- calculate_weighted_phyper_params(
    query_complete, pathway_complete, background
  )
  print("Complete overlap case:")
  print(params_complete$q_weighted)  # Should equal k_weighted
  
  # Case 3: Very high copy numbers
  query_high_copy <- data.frame(
    gene_name = c("REPEAT1"),
    copy_number = c(100)
  )
  pathway_high_copy <- c("REPEAT1")
  
  params_high_copy <- calculate_weighted_phyper_params(
    query_high_copy, pathway_high_copy, background
  )
  print("High copy case:")
  print(params_high_copy$k_weighted)  # Should be 100
}

test_edge_cases()
```

### Example 3: Mathematical Equivalence Verification

```r
# Comprehensive equivalence test
test_equivalence <- function() {
  set.seed(42)  # Reproducible results
  
  # Create test data
  query <- data.frame(
    gene_name = paste0("Q_GENE", 1:10),
    copy_number = sample(1:20, 10, replace = TRUE)
  )
  
  pathway <- paste0("Q_GENE", c(1, 3, 5, 7, 15, 17))  # Some overlap
  
  background <- data.frame(
    gene_name = paste0("GENE", 1:1000),
    copy_number = sample(1:5, 1000, replace = TRUE)
  )
  background <- rbind(background, query)  # Include query in background
  
  # Test equivalence
  equiv <- verify_equivalence_with_expansion(query, pathway, background)
  
  print("Equivalence Test Results:")
  print(paste("Parameters match:", equiv$parameters_equivalent))
  print(paste("P-values match:", equiv$pvalues_equivalent))
  print(paste("P-value difference:", equiv$pvalue_difference))
  print(paste("Weighted p-value:", equiv$weighted_pvalue))
  print(paste("Expansion p-value:", equiv$expansion_pvalue))
  
  # Performance comparison
  library(microbenchmark)
  
  print("\nPerformance Comparison:")
  timing <- microbenchmark(
    weighted = run_weighted_hypergeometric_test(query, pathway, background),
    expansion = {
      q_exp <- rep(query$gene_name, query$copy_number)
      b_exp <- rep(background$gene_name, background$copy_number)
      phyper(sum(q_exp %in% pathway) - 1, 
             sum(b_exp %in% pathway),
             length(b_exp) - sum(b_exp %in% pathway),
             length(q_exp), lower.tail = FALSE)
    },
    times = 50
  )
  print(timing)
}

# Run if microbenchmark is available
if (requireNamespace("microbenchmark", quietly = TRUE)) {
  test_equivalence()
}
```

## Integration with Existing Analysis

### Using with clusterProfiler

```r
# Convert copy-weighted results to clusterProfiler format
library(clusterProfiler)

# Expand gene list by copy number for compatibility
expand_for_clusterprofiler <- function(query_df) {
  expanded_genes <- rep(query_df$gene_name, query_df$copy_number)
  return(unique(expanded_genes))  # Remove duplicates for standard ORA
}

# Use with enrichGO
expanded_genes <- expand_for_clusterprofiler(phr_query)
ego_result <- enrichGO(
  gene = expanded_genes,
  universe = background$gene_name,
  OrgDb = org.Hs.eg.db::org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH"
)
```

### Comparison with Standard ORA

```r
# Compare standard vs weighted approaches
compare_approaches <- function(query_df, pathway_genes, background_df) {
  
  # Standard ORA
  standard_result <- list(
    overlap = length(intersect(query_df$gene_name, pathway_genes)),
    query_size = nrow(query_df),
    pathway_size = length(intersect(pathway_genes, background_df$gene_name))
  )
  standard_pval <- phyper(
    standard_result$overlap - 1,
    standard_result$pathway_size,
    nrow(background_df) - standard_result$pathway_size,
    standard_result$query_size,
    lower.tail = FALSE
  )
  
  # Weighted ORA
  weighted_result <- run_weighted_hypergeometric_test(
    query_df, pathway_genes, background_df
  )
  
  # Comparison
  comparison <- data.frame(
    method = c("Standard", "Copy-Weighted"),
    pvalue = c(standard_pval, weighted_result$pvalue),
    overlap = c(standard_result$overlap, weighted_result$observed_overlap_weighted),
    query_size = c(standard_result$query_size, weighted_result$query_size_weighted),
    fold_enrichment = c(
      standard_result$overlap / standard_result$query_size / 
      (standard_result$pathway_size / nrow(background_df)),
      weighted_result$fold_enrichment
    )
  )
  
  return(comparison)
}
```

## Best Practices

### 1. Data Preparation

- Ensure `gene_name` columns use consistent gene identifiers
- Verify `copy_number` values are positive integers
- Remove genes with zero copies (annotation errors)
- Use representative background that matches query context

### 2. Parameter Validation

- Always run validation checks
- Address warnings about small sample sizes
- Verify pathway genes exist in background
- Check for extreme copy numbers (>500) that might cause overflow

### 3. Statistical Interpretation

- Report both standard and weighted results when they differ
- Use appropriate multiple testing correction
- Consider biological relevance alongside statistical significance
- Document copy number model assumptions

### 4. Performance Optimization

- Use parameter weighting instead of instance expansion
- Cache background calculations for multiple pathway tests
- Consider normal approximation for very large parameters
- Monitor memory usage with large copy numbers

## Troubleshooting

### Common Issues

1. **Parameter validation failure**
   - Check for negative or zero copy numbers
   - Ensure pathway genes exist in background
   - Verify query-pathway overlap exists

2. **Numerical precision issues**
   - Use tolerance settings for equivalence tests
   - Consider normal approximation for huge parameters
   - Check for integer overflow with extreme copy numbers

3. **Performance problems**
   - Avoid instance expansion for large datasets
   - Use vectorized operations where possible
   - Cache repeated calculations

### Error Messages

- `"query_df must have columns: gene_name, copy_number"`: Check column names
- `"Parameter validation failed"`: Review constraint violations
- `"No valid genes in query_df after filtering zero copies"`: Check copy number data

## References

- `r_phyper_modifications_research.md`: Mathematical foundation
- `phyper_parameter_modification_analysis.md`: Detailed parameter analysis  
- `mathematical_framework_copy_weighted_ora.md`: Statistical framework
- `copy_number_phyper_mapping.R`: Implementation code