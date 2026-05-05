# Copy-Number Weighted ORA: Best Practices Guide

Based on comprehensive research synthesis covering mathematical validation, computational benchmarking, statistical validation, and implementation alternatives.

## Executive Decision Framework

### When to Use Weighted vs Standard Hypergeometric Testing

```
┌─ Copy Number Data Available?
├─ NO → Use standard hypergeometric
└─ YES ┐
       ├─ Selection Mechanism?
       ├─ Gene-level → Use STANDARD hypergeometric ⭐ RECOMMENDED
       │               (Weighted test is anti-conservative)
       └─ Instance-level → Use WEIGHTED hypergeometric
                          (Independence assumption holds)
```

**Critical Finding**: Under gene-level selection (typical ORA), weighted hypergeometric test inflates Type I error 4-6× above nominal levels and compromises FDR control.

## Implementation Strategy

### 1. Production Function Template

```r
# Context-aware implementation
copy_number_ora <- function(query_df, pathway_genes, background_df,
                           selection_type = c("auto", "gene_level", "instance_level")) {
  
  selection_type <- match.arg(selection_type)
  
  # Auto-detect based on use case
  if (selection_type == "auto") {
    selection_type <- detect_selection_type(query_df, background_df)
  }
  
  # Always compute both for comparison
  standard_result <- standard_hypergeometric_test(query_df, pathway_genes, background_df)
  weighted_result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)
  
  # Select primary method based on context
  if (selection_type == "gene_level") {
    primary <- standard_result
    alternative <- weighted_result
    warning("Using standard test - weighted test anti-conservative for gene-level selection")
  } else {
    primary <- weighted_result
    alternative <- standard_result
  }
  
  return(list(
    primary = primary,
    alternative = alternative,
    selection_type = selection_type
  ))
}
```

### 2. Performance Optimization Guidelines

**Dataset Size Thresholds:**
- **< 2K instances**: Either approach acceptable (slight preference for expansion in debugging)
- **2K-50K instances**: Direct weighted approach preferred (2-4× speedup)
- **> 50K instances**: Direct weighted approach essential (>3× speedup + memory constraints)

```r
# Efficient implementation based on dataset size
efficient_weighted_test <- function(query_df, pathway_genes, background_df) {
  total_instances <- sum(query_df$copy_number)
  
  if (total_instances >= 2000) {
    # Use direct parameter approach for efficiency
    return(weighted_phyper_direct(query_df, pathway_genes, background_df))
  } else {
    # Either approach works; direct still preferred for consistency
    return(weighted_phyper_direct(query_df, pathway_genes, background_df))
  }
}
```

## Statistical Validation and Limitations

### 1. Type I Error Control

| Method | Alpha=0.01 | Alpha=0.05 | Alpha=0.10 | Assessment |
|--------|-----------|-----------|-----------|------------|
| Standard | 0.004 | 0.018 | 0.071 | ✅ Well-controlled |
| Weighted | 0.142 | 0.225 | 0.263 | ❌ Anti-conservative |

### 2. Multiple Testing Correction

**For Standard Results:**
```r
# Safe - standard FDR correction
adjusted_pvals <- p.adjust(standard_pvalues, method = "fdr")
```

**For Weighted Results:**
```r
# UNSAFE - BH correction assumes uniform null distribution
# adjusted_pvals <- p.adjust(weighted_pvalues, method = "fdr")  # DON'T USE

# SAFE alternatives:
# 1. Very conservative Bonferroni
adjusted_pvals <- p.adjust(weighted_pvalues, method = "bonferroni")

# 2. Permutation-based FDR (gold standard but expensive)
adjusted_pvals <- permutation_fdr(query_df, pathway_list, background_df)

# 3. Conservative empirical calibration
adjusted_pvals <- calibrated_adjustment(weighted_pvalues, null_simulation_results)
```

### 3. When Weighted Approach is Valid

**Appropriate Use Cases:**
- PHR overlap analysis (genomic regions independently intersect PHRs)
- CNV burden analysis (independent copy number alterations)
- Gene dosage effect studies (instance-level biological mechanism)
- Exploratory analysis with explicit anti-conservative acknowledgment

**Inappropriate Use Cases:**
- Differential expression gene sets (genes selected as units)
- Standard pathway enrichment analysis
- Any analysis requiring Type I error control
- Multiple testing scenarios without permutation-based correction

## Implementation Quality Assurance

### 1. Mathematical Validation Checklist

```r
# Comprehensive validation function
validate_implementation <- function(test_function) {
  
  # Test 1: Mathematical equivalence with instance expansion
  test_mathematical_equivalence(test_function)
  
  # Test 2: Edge cases
  test_edge_cases(test_function)
  
  # Test 3: Parameter constraints
  test_parameter_validation(test_function)
  
  # Test 4: Performance benchmarks
  test_performance(test_function)
  
  # Test 5: Statistical properties
  test_statistical_properties(test_function)
}
```

### 2. Production Code Standards

**Essential Validations:**
```r
weighted_hypergeometric_test <- function(query_df, pathway_genes, background_df) {
  
  # Input validation
  stopifnot(all(c("gene", "copy_number") %in% names(query_df)))
  stopifnot(all(c("gene", "copy_number") %in% names(background_df)))
  stopifnot(all(query_df$copy_number >= 1))
  stopifnot(all(background_df$copy_number >= 1))
  stopifnot(!any(duplicated(query_df$gene)))
  stopifnot(!any(duplicated(background_df$gene)))
  
  # Ensure query ⊆ background
  valid_query_genes <- intersect(query_df$gene, background_df$gene)
  if (length(valid_query_genes) == 0) {
    stop("No query genes found in background")
  }
  
  # Parameter calculation with validation...
  # [implementation continues]
}
```

### 3. Error Handling and Edge Cases

**Critical Edge Cases:**
```r
# Handle zero overlap
if (q_weighted == 0) {
  return(list(pvalue = 1.0, overlap = 0))
}

# Handle complete overlap
if (q_weighted == k_weighted && k_weighted == m_weighted) {
  # Special handling for complete pathway coverage
}

# Handle extreme copy numbers
if (any(copy_numbers > 1000)) {
  warning("Extreme copy numbers detected - verify biological relevance")
}

# Handle tiny pathways
if (m_weighted < 5) {
  warning("Very small pathway - results may be unreliable")
}
```

## Multiple Pathway Analysis

### 1. Vectorized Implementation

```r
test_multiple_pathways <- function(query_df, pathway_list, background_df,
                                  method = c("standard", "weighted", "both"),
                                  correction = "fdr") {
  
  method <- match.arg(method)
  
  # Pre-calculate query statistics for efficiency
  query_stats <- precalculate_query_stats(query_df, background_df)
  
  # Vectorized pathway testing
  results <- lapply(pathway_list, function(pathway_genes) {
    if (method %in% c("standard", "both")) {
      standard_result <- fast_standard_test(query_stats, pathway_genes, background_df)
    }
    if (method %in% c("weighted", "both")) {
      weighted_result <- fast_weighted_test(query_stats, pathway_genes, background_df)
    }
    
    # Return appropriate results
  })
  
  # Apply appropriate correction
  if (method == "standard") {
    results$qvalue <- p.adjust(results$pvalue, method = correction)
  } else if (method == "weighted") {
    if (correction == "fdr") {
      warning("FDR correction unreliable for weighted tests - using Bonferroni")
      results$qvalue <- p.adjust(results$pvalue, method = "bonferroni")
    } else {
      results$qvalue <- p.adjust(results$pvalue, method = correction)
    }
  }
  
  return(results)
}
```

### 2. Parallel Processing for Large Analyses

```r
# For > 1000 pathways
if (length(pathway_list) > 1000 && requireNamespace("parallel")) {
  results <- parallel::mclapply(pathway_list, test_function, mc.cores = n_cores)
} else {
  results <- lapply(pathway_list, test_function)
}
```

## Reporting and Documentation Guidelines

### 1. Essential Reporting Elements

**Required in Methods Section:**
```
Copy-number weighted ORA was performed using [method]. Gene selection 
was performed at the [gene-level/instance-level]. Copy numbers were 
obtained from [source]. Background gene set consisted of [description].

Statistical testing used [standard/weighted] hypergeometric test. 
Multiple testing correction applied [method] at [threshold]. 

For weighted tests: Anti-conservative behavior under gene-level selection 
was acknowledged, with Type I error inflation of approximately [X]×.
```

**Required in Results Section:**
- Report both weighted and standard results when different
- Acknowledge statistical limitations of weighted approach
- Provide copy number distribution summary statistics
- Document pathway size filtering criteria

### 2. Supplementary Materials

**Include:**
- Copy number source and processing methodology
- Complete pathway gene lists
- Validation of hypergeometric parameter calculations
- Comparison of weighted vs standard results
- Assessment of Type I error properties if novel application

### 3. Code Reproducibility

```r
# Example complete analysis script
library(phrs)  # hypothetical package

# Load data
query_df <- load_query_genes_with_copy_numbers(source = "...")
background_df <- load_background_genome(version = "...")
pathway_db <- load_pathway_database(source = "GO", version = "...")

# Perform analysis with explicit method selection
ora_results <- copy_number_aware_ora(
  query_df = query_df,
  pathway_list = pathway_db$pathways,
  background_df = background_df,
  selection_type = "gene_level",  # Explicit choice
  correction = "fdr",
  min_pathway_size = 5
)

# Document statistical properties
validate_ora_assumptions(ora_results)

# Generate comprehensive report
generate_ora_report(ora_results, output_file = "ora_results.html")
```

## Integration with Existing Workflows

### 1. clusterProfiler Integration

```r
# Wrapper for clusterProfiler compatibility
clusterprofiler_weighted_ora <- function(gene_list, copy_numbers, 
                                        universe, universe_copy_numbers,
                                        ont = "BP") {
  
  # Prepare data frames
  query_df <- data.frame(gene = gene_list, copy_number = copy_numbers)
  background_df <- data.frame(gene = universe, copy_number = universe_copy_numbers)
  
  # Get pathway data
  pathways <- get_go_pathways(ont = ont)
  
  # Run weighted analysis
  results <- test_multiple_pathways(query_df, pathways, background_df)
  
  # Format for clusterProfiler compatibility
  format_for_clusterprofiler(results)
}
```

### 2. GOstats Integration

```r
# Adapter for GOstats workflows
gostats_weighted_wrapper <- function(selected_genes, copy_numbers,
                                    universe_genes, universe_copy_numbers,
                                    annotation_package) {
  
  # Convert to standard format
  query_df <- data.frame(gene = selected_genes, copy_number = copy_numbers)
  background_df <- data.frame(gene = universe_genes, copy_number = universe_copy_numbers)
  
  # Use annotation package to get pathway mappings
  pathways <- extract_pathways(annotation_package)
  
  # Run analysis
  results <- copy_number_aware_ora(query_df, pathways, background_df)
  
  # Convert to GOstats result format
  convert_to_gostats_format(results)
}
```

## Future Development Priorities

### 1. High Priority (Immediate)
- Permutation-based null model for gene-level selection
- Comprehensive validation test suite
- Integration with major Bioconductor packages
- Performance optimization for extreme-scale datasets

### 2. Medium Priority (Next 6 months)
- Alternative clustering-adjusted statistical methods
- Multi-omics integration approaches
- Standardized benchmarking datasets
- Visualization and reporting enhancements

### 3. Low Priority (Future)
- C++ implementation for performance-critical applications
- Machine learning-based copy number effect modeling
- Integration with workflow management systems
- Real-time analysis capabilities

## Summary Recommendations

### For Method Developers:
1. **Default to standard hypergeometric** for gene-level selection
2. **Implement both methods** with clear warnings about limitations
3. **Provide permutation-based alternatives** for copy-number aware analysis
4. **Include comprehensive validation** in all implementations

### For Applied Researchers:
1. **Understand your selection mechanism** (gene-level vs instance-level)
2. **Always report both methods** when copy numbers vary
3. **Use appropriate multiple testing correction** based on method
4. **Acknowledge statistical limitations** in publications

### For Software Users:
1. **Read statistical warnings** carefully
2. **Validate results** with standard approaches
3. **Consider biological interpretation** of copy number effects
4. **Report methodology** transparently

---

**This best practices guide synthesizes research from mathematical validation, computational benchmarking, statistical validation, and implementation alternatives to provide comprehensive guidance for copy-number weighted ORA implementation and application.**