# Copy-Number Weighted Hypergeometric Testing: Documentation and Usage Guide

## Overview

This implementation provides robust, production-ready functions for performing copy-number weighted over-representation analysis (ORA) using modified hypergeometric parameters. The approach transforms standard hypergeometric testing from unique gene counts to gene instance counts, accounting for copy number variation while maintaining statistical validity and computational efficiency.

## Key Features

- **Mathematical Equivalence**: Produces identical results to instance expansion but with superior computational efficiency
- **Comprehensive Validation**: Extensive input validation and constraint checking
- **Error Handling**: Graceful handling of edge cases (zero copies, extreme values, missing data)
- **Performance Optimization**: Significant memory and speed improvements over naive approaches
- **Statistical Rigor**: Maintains proper statistical properties and supports standard multiple testing correction

## Installation and Setup

```r
# Source the implementation
source("copy_weighted_hypergeometric.R")

# For comprehensive testing
source("test_copy_weighted_hypergeometric.R")

# For performance benchmarking
source("benchmark_copy_weighted_hypergeometric.R")
```

## Basic Usage

### Simple Example

```r
# Create query dataset (genes with copy numbers)
query_df <- data.frame(
  gene = c("GENE1", "GENE2", "GENE3", "GENE4"),
  copy_number = c(1, 5, 2, 8)
)

# Define pathway of interest
pathway_genes <- c("GENE1", "GENE2", "GENE5", "GENE6")

# Create background dataset
background_df <- data.frame(
  gene = paste0("GENE", 1:100),
  copy_number = rpois(100, lambda = 3) + 1  # Random copy numbers
)

# Perform copy-weighted hypergeometric test
result <- weighted_hypergeometric_test(
  query_df = query_df,
  pathway_genes = pathway_genes,
  background_df = background_df
)

# View results
print(paste("P-value:", format(result$pvalue, scientific = TRUE)))
print(paste("Fold enrichment:", round(result$fold_enrichment, 2)))
print(paste("Overlap instances:", result$overlap_instances))
```

### PHR-like Dataset Example

```r
# Simulate PHR-like data (high copy numbers)
set.seed(42)

# PHR genes: 35 genes with high copy numbers
phr_genes <- paste0("PHR", 1:35)
phr_copies <- rpois(35, lambda = 30) + 10  # High copy numbers

query_phr <- data.frame(
  gene = phr_genes,
  copy_number = phr_copies
)

# Human genome background with typical copy numbers
genome_genes <- paste0("GENE", 1:20000)
genome_copies <- rpois(20000, lambda = 2) + 1

background_genome <- data.frame(
  gene = genome_genes,
  copy_number = genome_copies
)

# Olfactory receptor pathway
or_pathway <- paste0("GENE", 1:400)  # First 400 genes are OR genes

# Test for OR enrichment in PHRs
or_result <- weighted_hypergeometric_test(
  query_df = query_phr,
  pathway_genes = or_pathway,
  background_df = background_genome
)

# Compare with standard approach
comparison <- compare_weighted_vs_standard(
  query_phr, or_pathway, background_genome
)

print("=== COPY-WEIGHTED vs STANDARD COMPARISON ===")
print(paste("Standard p-value:", format(comparison$standard$pvalue, scientific = TRUE)))
print(paste("Weighted p-value:", format(comparison$weighted$pvalue, scientific = TRUE)))
print(paste("Fold change in significance:", round(comparison$comparison$pvalue_ratio, 3)))
```

## Function Reference

### Main Functions

#### `weighted_hypergeometric_test()`

Primary function for copy-weighted hypergeometric testing.

**Parameters:**
- `query_df`: Data frame with columns 'gene' and 'copy_number' for query set
- `pathway_genes`: Vector of gene symbols defining the pathway of interest
- `background_df`: Data frame with columns 'gene' and 'copy_number' for background
- `validate_inputs`: Logical, whether to perform comprehensive input validation (default: TRUE)
- `handle_zeros`: Logical, whether to automatically remove zero-copy genes (default: TRUE)
- `max_copies`: Integer, maximum allowed copy number to cap extremes (default: 500)
- `min_instances`: Integer, minimum total instances required for reliable results (default: 10)

**Returns:**
- `pvalue`: Hypergeometric test p-value
- `overlap_instances`: Copy-weighted overlap count
- `query_instances`: Total query instances
- `pathway_instances`: Total pathway instances in background
- `background_instances`: Total background instances
- `fold_enrichment`: Fold enrichment ratio
- `parameters`: Named vector of (q, m, n, k) hypergeometric parameters
- `diagnostics`: Validation and diagnostic information

#### `compare_weighted_vs_standard()`

Utility function to compare copy-weighted and standard results.

**Parameters:**
- `query_df`: Query data frame
- `pathway_genes`: Pathway gene vector
- `background_df`: Background data frame

**Returns:**
Comprehensive comparison including both p-values, parameter scaling factors, and significance differences.

### Data Preparation

#### Required Data Format

All input datasets must be data frames with exactly these columns:

```r
# Correct format
correct_df <- data.frame(
  gene = c("GENE1", "GENE2", "GENE3"),        # Character: gene symbols
  copy_number = c(1, 5, 2)                    # Numeric: copy counts (positive integers)
)

# Common mistakes to avoid
wrong_df <- data.frame(
  symbol = c("GENE1", "GENE2"),              # Wrong: column should be "gene"
  copies = c(1, 5)                           # Wrong: column should be "copy_number"
)
```

#### Data Quality Checks

The implementation automatically checks for:

- **Missing columns**: Must have 'gene' and 'copy_number'
- **Data types**: copy_number must be numeric
- **Invalid values**: No negative, zero, infinite, or non-integer copy numbers
- **Duplicated genes**: Each gene should appear only once per dataset
- **Empty datasets**: Must have at least 2 genes in query
- **Pathway coverage**: At least some pathway genes must be in background

## Statistical Interpretation

### Understanding Results

```r
result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)

# P-value interpretation
if (result$pvalue < 0.05) {
  cat("Significant enrichment detected\n")
} else {
  cat("No significant enrichment\n")
}

# Effect size interpretation  
if (result$fold_enrichment > 2) {
  cat("Strong enrichment (>2x expected)\n")
} else if (result$fold_enrichment > 1.5) {
  cat("Moderate enrichment (1.5-2x expected)\n")
} else {
  cat("Weak or no enrichment\n")
}

# Sample size adequacy
if (result$query_instances < 10) {
  warning("Small sample size - results may be unreliable")
}
```

### Multiple Testing Correction

```r
# Test multiple pathways
pathways <- list(
  "Olfactory" = or_genes,
  "Immune" = immune_genes,
  "Metabolic" = metabolic_genes
)

# Run tests
pvalues <- sapply(pathways, function(pathway) {
  result <- weighted_hypergeometric_test(query_df, pathway, background_df)
  return(result$pvalue)
})

# Apply FDR correction
adjusted_pvals <- p.adjust(pvalues, method = "fdr")

# Results table
results_table <- data.frame(
  Pathway = names(pvalues),
  Raw_Pvalue = pvalues,
  FDR_Adjusted = adjusted_pvals,
  Significant = adjusted_pvals < 0.05
)

print(results_table)
```

## Performance Optimization

### When to Use Copy-Weighted Testing

**Strong Benefits:**
- Datasets with high copy number variation (CV > 1.0)
- Large datasets (>1000 genes, >10,000 total instances)
- PHR-like data with extreme copy numbers
- Memory-constrained environments

**Moderate Benefits:**
- Medium-sized datasets with moderate copy variation
- Production pipelines requiring robust error handling

**Minimal Benefits:**
- Small datasets with uniform copy numbers
- Single-use analyses where development time matters more than performance

### Performance Tuning

```r
# For maximum performance (production pipelines)
fast_result <- weighted_hypergeometric_test(
  query_df, pathway_genes, background_df,
  validate_inputs = FALSE,    # Skip validation if data is pre-validated
  handle_zeros = FALSE,       # Skip zero filtering if unnecessary  
  max_copies = NULL          # Skip copy capping if not needed
)

# For maximum safety (exploratory analysis)
safe_result <- weighted_hypergeometric_test(
  query_df, pathway_genes, background_df,
  validate_inputs = TRUE,     # Full validation
  handle_zeros = TRUE,        # Auto-remove problematic genes
  max_copies = 100,          # Cap extreme values
  min_instances = 20         # Require larger sample sizes
)
```

### Benchmarking Your Data

```r
# Run performance benchmark on your specific datasets
source("benchmark_copy_weighted_hypergeometric.R")

# Create test data matching your scale
test_data <- list(
  query_df = your_query_data,
  pathway_genes = your_pathway_genes,
  background_df = your_background_data
)

# Quick timing comparison
timing <- microbenchmark(
  weighted = weighted_hypergeometric_test(
    test_data$query_df, test_data$pathway_genes, test_data$background_df
  ),
  times = 10
)

print(summary(timing))
```

## Troubleshooting

### Common Error Messages

**"Input validation failed: query_df missing columns"**
- Check that your data frame has columns named exactly "gene" and "copy_number"
- Ensure column names have no extra spaces or special characters

**"Hypergeometric constraint violation: overlap cannot exceed query instances"**
- This suggests data inconsistency between query and background
- Check that copy numbers are consistent between datasets
- Verify gene names match exactly (case-sensitive)

**"No pathway genes found in background dataset"**
- Gene names don't match between pathway and background
- Check for naming convention differences (e.g., "GENE1" vs "gene1")
- Verify pathway genes actually exist in your background model

**Warning: "Low correlation between query and background copy numbers"**
- Query and background may represent different biological contexts
- Consider if the background model is appropriate for your query
- May indicate batch effects or systematic differences

### Debugging Tips

```r
# Inspect intermediate calculations
debug_result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)

# Check parameter calculations
print("=== PARAMETER BREAKDOWN ===")
print(paste("Query instances:", debug_result$query_instances))
print(paste("Overlap instances:", debug_result$overlap_instances))  
print(paste("Pathway instances:", debug_result$pathway_instances))
print(paste("Background instances:", debug_result$background_instances))

# Check diagnostics
print("=== DIAGNOSTICS ===")
print(debug_result$diagnostics)

# Validate gene name matching
query_genes_in_bg <- sum(query_df$gene %in% background_df$gene)
pathway_genes_in_bg <- sum(pathway_genes %in% background_df$gene)
pathway_genes_in_query <- sum(pathway_genes %in% query_df$gene)

print(paste("Query genes in background:", query_genes_in_bg, "/", nrow(query_df)))
print(paste("Pathway genes in background:", pathway_genes_in_bg, "/", length(pathway_genes)))
print(paste("Pathway genes in query:", pathway_genes_in_query, "/", length(pathway_genes)))
```

## Validation and Quality Assurance

### Running the Test Suite

```r
# Run comprehensive tests
source("test_copy_weighted_hypergeometric.R")
test_results <- run_all_tests()

if (test_results$all_passed) {
  cat("✓ All tests passed - implementation validated\n")
} else {
  cat("✗ Some tests failed - check implementation\n")
}
```

### Statistical Validation

```r
# Test null distribution properties (simplified)
set.seed(123)
n_simulations <- 100

null_pvals <- replicate(n_simulations, {
  # Generate random query with no true enrichment
  random_genes <- sample(background_df$gene, size = 20)
  random_query <- background_df[background_df$gene %in% random_genes, ]
  
  result <- weighted_hypergeometric_test(random_query, pathway_genes, background_df)
  return(result$pvalue)
})

# Check uniform distribution under null hypothesis
ks_test <- ks.test(null_pvals, "punif")
if (ks_test$p.value > 0.05) {
  cat("✓ Null distribution appears uniform (good)\n")
} else {
  cat("! Null distribution may be biased\n")
}

# Check Type I error rate
alpha <- 0.05
false_positive_rate <- mean(null_pvals < alpha)
cat(paste("Observed Type I error rate:", round(false_positive_rate, 3),
          "vs expected", alpha, "\n"))
```

## Advanced Usage

### Custom Background Models

```r
# Create tissue-specific background
tissue_background <- create_tissue_background(
  tissue_name = "brain",
  expressed_genes = brain_expressed_genes,
  copy_number_source = "genomic_coordinates"
)

# Use with copy-weighted testing
brain_result <- weighted_hypergeometric_test(
  query_df = brain_phr_genes,
  pathway_genes = neurotransmitter_pathway,
  background_df = tissue_background
)
```

### Batch Processing

```r
# Process multiple datasets
datasets <- list(
  "Dataset1" = list(query = query1, background = bg1),
  "Dataset2" = list(query = query2, background = bg2)
)

pathway_results <- lapply(names(datasets), function(dataset_name) {
  data <- datasets[[dataset_name]]
  
  results <- lapply(pathway_list, function(pathway) {
    weighted_hypergeometric_test(data$query, pathway, data$background)
  })
  
  return(results)
})

names(pathway_results) <- names(datasets)
```

### Integration with Existing Workflows

```r
# Convert from common formats
convert_from_gene_list <- function(gene_list, copy_numbers) {
  data.frame(
    gene = gene_list,
    copy_number = copy_numbers
  )
}

# Export results to standard format
export_to_csv <- function(result, filename) {
  result_df <- data.frame(
    pvalue = result$pvalue,
    fold_enrichment = result$fold_enrichment,
    overlap_instances = result$overlap_instances,
    query_instances = result$query_instances,
    significant = result$pvalue < 0.05
  )
  
  write.csv(result_df, filename, row.names = FALSE)
}
```

## Best Practices

1. **Always validate inputs** unless you're certain about data quality
2. **Use appropriate background models** that match your biological context
3. **Apply multiple testing correction** when testing multiple pathways
4. **Check copy number model consistency** between query and background
5. **Report both standard and weighted results** when they differ substantially
6. **Include diagnostic information** in results for reproducibility
7. **Benchmark on your specific data** to understand performance characteristics

## Citation and References

This implementation is based on comprehensive research demonstrating mathematical equivalence between parameter weighting and instance expansion approaches for copy-number weighted hypergeometric testing. Key insights include:

- Parameter transformation maintains statistical validity
- Computational efficiency gains are substantial for high copy number data
- Standard hypergeometric theory applies after parameter transformation

For technical details, see the research documents:
- `r_phyper_modifications_research.md`
- `phyper_parameter_modification_analysis.md`
- `mathematical_framework_copy_weighted_ora.md`

## Version History

- **v1.0 (2026-04-01)**: Initial production release
  - Complete implementation with validation
  - Comprehensive test suite
  - Performance benchmarking
  - Full documentation

---

*This documentation provides comprehensive guidance for using copy-number weighted hypergeometric testing in R. For additional questions or issues, refer to the test suite and benchmark scripts for examples.*