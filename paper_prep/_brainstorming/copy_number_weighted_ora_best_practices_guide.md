# Copy-Number-Weighted ORA Best Practices Guide

**Task:** synthesize-parameter-mapping  
**Date:** 2026-04-01  
**Version:** 1.0  

## Overview

This guide provides best practices for implementing and using copy-number-weighted over-representation analysis (ORA) based on the comprehensive validation results from mathematical verification, performance benchmarking, and integration testing.

## Quick Start Checklist

**Pre-Analysis Requirements:**
- [ ] Copy number data for all genes (query + background)
- [ ] Pathway definitions with consistent gene identifiers
- [ ] Minimum 15,000 background genes for statistical stability
- [ ] Copy number consistency validation completed

**Quality Control Standards:**
- [ ] Zero copy number genes filtered out
- [ ] Gene identifier consistency verified across datasets
- [ ] Background represents appropriate null distribution
- [ ] Copy number source documented and validated

## Data Preparation Best Practices

### 1. Copy Number Data Quality

**Recommended Sources (in order of preference):**
1. **Experimental validation:** qPCR, digital PCR, or targeted sequencing
2. **High-resolution genomic data:** Long-read sequencing assemblies
3. **Computational predictions:** RefSeq or Ensembl annotations (with validation)

**Quality Control Steps:**
```r
# Example quality control workflow
validate_copy_data <- function(copy_data) {
  # Remove invalid entries
  clean_data <- copy_data[copy_data$copy_number > 0, ]
  
  # Check for extreme values
  extreme_copies <- clean_data[clean_data$copy_number > 50, ]
  if (nrow(extreme_copies) > 0) {
    warning("Found genes with >50 copies - verify these are correct")
  }
  
  # Validate completeness
  missing_fraction <- sum(is.na(clean_data$copy_number)) / nrow(clean_data)
  if (missing_fraction > 0.05) {
    stop("More than 5% missing copy number data")
  }
  
  return(clean_data)
}
```

### 2. Background Construction

**Genome-wide Background:**
- Include all protein-coding genes with known copy numbers
- Minimum 15,000 genes for statistical stability
- Use same copy number source as query genes
- Filter out genes not expressed in relevant tissues (optional)

**Tissue-specific Background:**
- Include genes expressed in relevant cell types/tissues
- Maintain breadth to avoid statistical bias
- Document filtering criteria for reproducibility

### 3. Copy Number Consistency Validation

**Critical Issue:** Query genes may have higher copy numbers than the same genes in background data.

**Detection:**
```r
detect_copy_inconsistencies <- function(query_df, background_df) {
  overlaps <- merge(query_df, background_df, by = "gene_name", suffixes = c("_query", "_bg"))
  violations <- overlaps[overlaps$copy_number_query > overlaps$copy_number_bg, ]
  
  if (nrow(violations) > 0) {
    message(paste("Found", nrow(violations), "copy number inconsistencies"))
    print(violations[, c("gene_name", "copy_number_query", "copy_number_bg")])
    return(violations)
  }
  
  message("No copy number inconsistencies found")
  return(NULL)
}
```

**Recommended Solutions:**

**Option A: Background Adjustment (Recommended)**
```r
adjust_background_copies <- function(query_df, background_df) {
  # Increase background copy numbers to match query
  for (gene in query_df$gene_name) {
    query_copies <- query_df[query_df$gene_name == gene, "copy_number"]
    bg_idx <- background_df$gene_name == gene
    
    if (sum(bg_idx) > 0) {
      background_df[bg_idx, "copy_number"] <- pmax(
        background_df[bg_idx, "copy_number"], 
        query_copies
      )
    }
  }
  
  return(background_df)
}
```

**Option B: Query Capping (Conservative)**
```r
cap_query_copies <- function(query_df, background_df) {
  # Reduce query copy numbers to match background
  for (gene in query_df$gene_name) {
    bg_copies <- background_df[background_df$gene_name == gene, "copy_number"]
    
    if (length(bg_copies) > 0) {
      query_df[query_df$gene_name == gene, "copy_number"] <- pmin(
        query_df[query_df$gene_name == gene, "copy_number"],
        bg_copies[1]
      )
    }
  }
  
  return(query_df)
}
```

## Implementation Guidelines

### 1. Recommended Production Workflow

**Complete Analysis Pipeline:**
```r
run_production_copy_weighted_ora <- function(query_genes, pathway_db, background_data, 
                                           consistency_method = "adjust_background") {
  
  # Step 1: Data validation and cleaning
  query_df <- validate_copy_data(query_genes)
  background_df <- validate_copy_data(background_data)
  
  # Step 2: Copy consistency handling
  inconsistencies <- detect_copy_inconsistencies(query_df, background_df)
  
  if (!is.null(inconsistencies)) {
    if (consistency_method == "adjust_background") {
      background_df <- adjust_background_copies(query_df, background_df)
      message("Applied background adjustment for copy consistency")
    } else {
      query_df <- cap_query_copies(query_df, background_df)
      message("Applied query capping for copy consistency")
    }
  }
  
  # Step 3: Pathway enrichment analysis
  results <- list()
  for (pathway_name in names(pathway_db)) {
    pathway_genes <- pathway_db[[pathway_name]]
    
    # Calculate weighted parameters
    params <- calculate_weighted_phyper_params(
      query_df, pathway_genes, background_df, 
      validate_params = TRUE
    )
    
    # Run hypergeometric test
    pvalue <- phyper(params$q - 1, params$m, params$n, params$k, lower.tail = FALSE)
    
    # Calculate enrichment metrics
    fold_enrichment <- (params$q / params$k) / (params$m / (params$m + params$n))
    
    results[[pathway_name]] <- list(
      pvalue = pvalue,
      fold_enrichment = fold_enrichment,
      query_instances = params$q,
      total_query_instances = params$k,
      pathway_background_instances = params$m,
      parameters = params
    )
  }
  
  # Step 4: Multiple testing correction
  pvalues <- sapply(results, function(x) x$pvalue)
  fdr_values <- p.adjust(pvalues, method = "fdr")
  
  for (i in seq_along(results)) {
    results[[i]]$fdr <- fdr_values[i]
  }
  
  return(results)
}
```

### 2. Parameter Calculation Best Practices

**Validation Enabled:**
```r
calculate_weighted_phyper_params <- function(query_df, pathway_genes, background_df, 
                                          validate_params = TRUE) {
  
  # Calculate weighted parameters
  k_weighted <- sum(query_df$copy_number)  # Total query instances
  
  # Find overlapping genes and their copy numbers
  overlap_genes <- query_df[query_df$gene_name %in% pathway_genes, ]
  q_weighted <- sum(overlap_genes$copy_number)  # Overlap instances
  
  # Background pathway instances
  pathway_background <- background_df[background_df$gene_name %in% pathway_genes, ]
  m_weighted <- sum(pathway_background$copy_number)
  
  # Non-pathway background instances
  total_background <- sum(background_df$copy_number)
  n_weighted <- total_background - m_weighted
  
  params <- list(
    q = q_weighted,
    m = m_weighted,
    n = n_weighted,
    k = k_weighted
  )
  
  # Validation
  if (validate_params) {
    validation <- validate_hypergeometric_params(params)
    if (!validation$passed) {
      stop(paste("Parameter validation failed:", 
                paste(validation$errors, collapse = "; ")))
    }
    params$validation <- validation
  }
  
  return(params)
}

validate_hypergeometric_params <- function(params) {
  errors <- character(0)
  
  # Check constraints
  if (params$q > params$k) {
    errors <- c(errors, "q > k: overlap instances exceed total query instances")
  }
  
  if (params$q > params$m) {
    errors <- c(errors, "q > m: overlap instances exceed pathway instances")
  }
  
  if (params$k > params$m + params$n) {
    errors <- c(errors, "k > m+n: query instances exceed total background")
  }
  
  # Check for negative values
  if (any(sapply(params[1:4], function(x) x < 0))) {
    errors <- c(errors, "Negative parameter values detected")
  }
  
  list(passed = length(errors) == 0, errors = errors)
}
```

### 3. Error Handling Best Practices

**Comprehensive Error Checking:**
```r
robust_copy_weighted_test <- function(query_df, pathway_genes, background_df) {
  
  # Input validation
  if (nrow(query_df) == 0) {
    stop("Empty query dataset provided")
  }
  
  if (length(pathway_genes) == 0) {
    warning("Empty pathway provided")
    return(list(pvalue = 1.0, warning = "Empty pathway"))
  }
  
  if (nrow(background_df) < 1000) {
    warning("Small background size may affect statistical power")
  }
  
  # Check for required columns
  required_cols <- c("gene_name", "copy_number")
  if (!all(required_cols %in% colnames(query_df))) {
    stop(paste("Missing required columns:", 
              paste(setdiff(required_cols, colnames(query_df)), collapse = ", ")))
  }
  
  # Calculate parameters with validation
  tryCatch({
    params <- calculate_weighted_phyper_params(query_df, pathway_genes, background_df)
    
    # Handle edge cases
    if (params$q == 0) {
      return(list(pvalue = 1.0, note = "No overlap between query and pathway"))
    }
    
    if (params$m == 0) {
      return(list(pvalue = 1.0, note = "Pathway not present in background"))
    }
    
    # Run test
    pvalue <- phyper(params$q - 1, params$m, params$n, params$k, lower.tail = FALSE)
    
    return(list(pvalue = pvalue, parameters = params))
    
  }, error = function(e) {
    stop(paste("Analysis failed:", e$message))
  })
}
```

## Performance Optimization Guidelines

### 1. Memory Management

**For Large Datasets (>20K genes):**
```r
# Process pathways in batches to manage memory
batch_pathway_analysis <- function(query_df, pathway_list, background_df, batch_size = 100) {
  
  pathway_names <- names(pathway_list)
  n_pathways <- length(pathway_names)
  
  results <- list()
  
  for (i in seq(1, n_pathways, batch_size)) {
    end_idx <- min(i + batch_size - 1, n_pathways)
    batch_pathways <- pathway_list[i:end_idx]
    
    # Process batch
    batch_results <- lapply(batch_pathways, function(pathway_genes) {
      robust_copy_weighted_test(query_df, pathway_genes, background_df)
    })
    
    results <- c(results, batch_results)
    
    # Memory cleanup
    gc()
    
    message(paste("Processed", end_idx, "of", n_pathways, "pathways"))
  }
  
  return(results)
}
```

### 2. Performance Monitoring

**Key Metrics to Track:**
```r
performance_monitor <- function(analysis_function, ...) {
  
  start_time <- Sys.time()
  start_memory <- as.numeric(object.size(ls(envir = .GlobalEnv)))
  
  # Run analysis
  results <- analysis_function(...)
  
  end_time <- Sys.time()
  end_memory <- as.numeric(object.size(ls(envir = .GlobalEnv)))
  
  # Performance summary
  performance <- list(
    duration = as.numeric(end_time - start_time, units = "secs"),
    memory_change = end_memory - start_memory,
    timestamp = end_time
  )
  
  message(sprintf("Analysis completed in %.2f seconds", performance$duration))
  message(sprintf("Memory change: %.2f MB", performance$memory_change / 1e6))
  
  attr(results, "performance") <- performance
  return(results)
}
```

## Integration Workflows

### 1. g:Profiler Integration

**Complete g:Profiler Workflow:**
```r
gprofiler_copy_weighted_analysis <- function(gprofiler_json_file, copy_data_file) {
  
  # Load g:Profiler request
  gprofiler_request <- fromJSON(gprofiler_json_file)
  
  # Load copy number data
  copy_data <- read.csv(copy_data_file)
  
  # Create weighted query
  query_df <- data.frame(
    gene_name = gprofiler_request$query,
    copy_number = copy_data[match(gprofiler_request$query, copy_data$gene_name), "total_copies"]
  )
  
  # Remove genes without copy data
  query_df <- query_df[!is.na(query_df$copy_number), ]
  
  # Create background from copy data (all protein-coding genes)
  background_df <- copy_data[copy_data$gene_biotype == "protein_coding", 
                            c("gene_name", "total_copies")]
  colnames(background_df)[2] <- "copy_number"
  
  # Load pathway databases (example for GO)
  go_pathways <- load_go_pathways(gprofiler_request$sources)
  
  # Run copy-weighted analysis
  results <- run_production_copy_weighted_ora(
    query_df, go_pathways, background_df
  )
  
  # Format results for g:Profiler compatibility
  formatted_results <- format_for_gprofiler(results)
  
  return(formatted_results)
}
```

### 2. Python Workflow Bridge

**R-Python Parameter Mapping:**
```r
export_for_python <- function(query_df, pathway_genes, background_df, output_file) {
  
  # Calculate weighted parameters
  params <- calculate_weighted_phyper_params(query_df, pathway_genes, background_df)
  
  # Create Python-compatible format
  python_params <- list(
    M = params$m + params$n,  # Total population for scipy
    n = params$m,             # Success states  
    N = params$k,             # Sample size
    q = params$q              # Overlap (for sf calculation: q-1)
  )
  
  # Save for Python
  writeLines(jsonlite::toJSON(python_params, auto_unbox = TRUE), output_file)
  
  message(paste("Parameters exported for Python analysis:", output_file))
  return(python_params)
}
```

**Corresponding Python Code:**
```python
import json
from scipy.stats import hypergeom

def import_from_r(params_file):
    with open(params_file, 'r') as f:
        params = json.load(f)
    
    # Calculate p-value using scipy
    pvalue = hypergeom.sf(params['q'] - 1, params['M'], params['n'], params['N'])
    
    return {
        'pvalue': pvalue,
        'parameters': params
    }
```

## Quality Assurance and Validation

### 1. Automated Testing

**Unit Test Framework:**
```r
test_copy_weighted_implementation <- function() {
  
  # Test 1: Mathematical equivalence
  test_mathematical_equivalence()
  
  # Test 2: Parameter constraints
  test_parameter_constraints()
  
  # Test 3: Edge cases
  test_edge_cases()
  
  # Test 4: Performance benchmarks
  test_performance_benchmarks()
  
  message("All validation tests completed")
}

test_mathematical_equivalence <- function() {
  # Create test data
  query_df <- data.frame(
    gene_name = paste0("GENE", 1:100),
    copy_number = sample(1:5, 100, replace = TRUE)
  )
  
  pathway_genes <- paste0("GENE", sample(1:100, 20))
  
  background_df <- data.frame(
    gene_name = paste0("GENE", 1:1000),
    copy_number = sample(1:3, 1000, replace = TRUE)
  )
  
  # Parameter weighting result
  weighted_result <- robust_copy_weighted_test(query_df, pathway_genes, background_df)
  
  # Instance expansion result (for validation)
  expansion_result <- instance_expansion_test(query_df, pathway_genes, background_df)
  
  # Compare results
  diff <- abs(weighted_result$pvalue - expansion_result$pvalue)
  
  if (diff < 1e-12) {
    message("✅ Mathematical equivalence test PASSED")
  } else {
    stop("❌ Mathematical equivalence test FAILED: p-value difference = ", diff)
  }
}
```

### 2. Production Monitoring

**Continuous Validation:**
```r
monitor_production_quality <- function(results, reference_method = NULL) {
  
  # Check for statistical anomalies
  pvalues <- sapply(results, function(x) x$pvalue)
  
  # P-value distribution checks
  if (sum(pvalues == 0) / length(pvalues) > 0.1) {
    warning("High proportion of p-values = 0 detected")
  }
  
  if (sum(pvalues == 1) / length(pvalues) > 0.5) {
    warning("High proportion of p-values = 1 detected")
  }
  
  # Enrichment ratio checks
  fold_enrichments <- sapply(results, function(x) x$fold_enrichment)
  
  if (any(fold_enrichments < 0)) {
    warning("Negative fold enrichment detected - check implementation")
  }
  
  # Cross-validation with reference method
  if (!is.null(reference_method)) {
    reference_results <- reference_method(...)
    
    pvalue_correlation <- cor(pvalues, reference_results$pvalues, use = "complete.obs")
    
    if (pvalue_correlation < 0.95) {
      warning(paste("Low correlation with reference method:", pvalue_correlation))
    } else {
      message(paste("✅ High correlation with reference method:", pvalue_correlation))
    }
  }
}
```

## Troubleshooting Guide

### Common Issues and Solutions

**1. Constraint Violation Errors**
```
Error: q > m (overlap instances exceed pathway instances)
```
**Solution:** Apply copy number consistency validation before analysis:
```r
background_df <- adjust_background_copies(query_df, background_df)
```

**2. Memory Errors with Large Datasets**
```
Error: cannot allocate vector of size X GB
```
**Solution:** Use batched processing:
```r
results <- batch_pathway_analysis(query_df, pathway_list, background_df, batch_size = 50)
```

**3. Zero P-values**
```
Warning: Many p-values = 0 detected
```
**Solution:** Check for extreme enrichment and consider log-space calculations:
```r
log_pvalue <- phyper(params$q - 1, params$m, params$n, params$k, 
                    lower.tail = FALSE, log.p = TRUE)
```

**4. Performance Issues**
```
Analysis taking too long for large datasets
```
**Solution:** Profile and optimize:
```r
# Use vectorized operations
# Reduce pathway database size
# Consider parallel processing
```

### Diagnostic Functions

**Comprehensive Diagnostics:**
```r
diagnose_copy_weighted_analysis <- function(query_df, pathway_db, background_df) {
  
  cat("=== Copy-Weighted ORA Diagnostics ===\n")
  
  # Dataset statistics
  cat("Query genes:", nrow(query_df), "\n")
  cat("Background genes:", nrow(background_df), "\n")
  cat("Pathways:", length(pathway_db), "\n")
  
  # Copy number statistics
  cat("Query copy range:", min(query_df$copy_number), "-", max(query_df$copy_number), "\n")
  cat("Background copy range:", min(background_df$copy_number), "-", max(background_df$copy_number), "\n")
  
  # Check consistency
  inconsistencies <- detect_copy_inconsistencies(query_df, background_df)
  cat("Copy inconsistencies:", ifelse(is.null(inconsistencies), "None", nrow(inconsistencies)), "\n")
  
  # Memory estimate
  estimated_memory <- estimate_memory_usage(query_df, pathway_db, background_df)
  cat("Estimated memory usage:", round(estimated_memory / 1e6, 2), "MB\n")
  
  # Pathway overlap statistics
  overlaps <- sapply(pathway_db, function(pathway) {
    sum(query_df$gene_name %in% pathway)
  })
  
  cat("Pathway overlaps - Mean:", round(mean(overlaps), 1), 
      "Range:", min(overlaps), "-", max(overlaps), "\n")
  
  cat("=== Diagnostics Complete ===\n")
}
```

## Conclusion

This best practices guide provides comprehensive guidance for implementing copy-number-weighted ORA in production environments. Following these guidelines ensures:

- ✅ Mathematical correctness through proper validation
- ✅ Computational efficiency through optimized implementation  
- ✅ Integration compatibility with existing workflows
- ✅ Robust error handling and quality assurance

**Key Success Factors:**
1. Proper copy number consistency handling
2. Comprehensive input validation
3. Performance monitoring and optimization
4. Continuous quality assurance

The copy-number-weighted approach enables genome-scale enrichment analysis while maintaining statistical rigor and computational feasibility.

---

**Generated by:** synthesize-parameter-mapping task  
**Based on:** Comprehensive validation from mathematical-verification-of, performance-benchmarking-parameter, and integration-testing-with tasks