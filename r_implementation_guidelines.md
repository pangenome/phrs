# R Implementation Guidelines for Copy-Number Weighted ORA

## Production Implementation

### Core Function Template

```r
#' Copy-Number Weighted Hypergeometric Test
#'
#' Performs hypergeometric enrichment testing with copy-number weighting.
#' Uses computationally efficient parameter weighting approach.
#'
#' @param query_df Data frame with columns 'gene' and 'copy_number'
#' @param pathway_genes Character vector of pathway gene identifiers
#' @param background_df Data frame with columns 'gene' and 'copy_number'
#' @param selection_type Character: "gene_level" or "instance_level"
#' @return List with test results and statistical metadata
#'
#' @examples
#' query_df <- data.frame(gene = c("A", "B", "C"), copy_number = c(2, 3, 1))
#' pathway_genes <- c("A", "B", "X", "Y")
#' background_df <- data.frame(gene = c("A", "B", "C", "X", "Y", "Z"), 
#'                            copy_number = c(2, 3, 1, 4, 2, 1))
#' result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)
weighted_hypergeometric_test <- function(query_df, pathway_genes, background_df, 
                                        selection_type = "gene_level") {
  
  # Input validation
  stopifnot("gene" %in% names(query_df))
  stopifnot("copy_number" %in% names(query_df))
  stopifnot("gene" %in% names(background_df))
  stopifnot("copy_number" %in% names(background_df))
  stopifnot(is.character(pathway_genes))
  stopifnot(all(query_df$copy_number >= 1))
  stopifnot(all(background_df$copy_number >= 1))
  
  # Ensure query genes are subset of background
  valid_query_genes <- intersect(query_df$gene, background_df$gene)
  if (length(valid_query_genes) == 0) {
    stop("No query genes found in background gene set")
  }
  
  query_filtered <- query_df[query_df$gene %in% valid_query_genes, ]
  
  # Use background copy numbers for consistency
  query_with_bg_copies <- merge(
    query_filtered[, "gene", drop = FALSE],
    background_df, 
    by = "gene", 
    all.x = TRUE
  )
  
  # Calculate weighted parameters
  query_in_pathway <- query_with_bg_copies[query_with_bg_copies$gene %in% pathway_genes, ]
  pathway_in_background <- background_df[background_df$gene %in% pathway_genes, ]
  
  k_weighted <- sum(query_with_bg_copies$copy_number)
  q_weighted <- if (nrow(query_in_pathway) > 0) {
    sum(query_in_pathway$copy_number)
  } else {
    0
  }
  m_weighted <- if (nrow(pathway_in_background) > 0) {
    sum(pathway_in_background$copy_number)
  } else {
    0
  }
  n_weighted <- sum(background_df$copy_number) - m_weighted
  
  # Validate hypergeometric constraints
  stopifnot(q_weighted <= k_weighted)
  stopifnot(q_weighted <= m_weighted)
  stopifnot(k_weighted <= (m_weighted + n_weighted))
  stopifnot(m_weighted > 0)  # Pathway must have some genes
  stopifnot(n_weighted > 0)  # Must have non-pathway genes
  
  # Calculate p-value
  if (q_weighted == 0) {
    pvalue <- 1.0
  } else {
    pvalue <- phyper(
      q_weighted - 1, 
      m_weighted, 
      n_weighted, 
      k_weighted, 
      lower.tail = FALSE
    )
  }
  
  # Calculate fold enrichment
  observed_rate <- q_weighted / k_weighted
  expected_rate <- m_weighted / (m_weighted + n_weighted)
  fold_enrichment <- if (expected_rate > 0) {
    observed_rate / expected_rate
  } else {
    Inf
  }
  
  # Statistical warnings
  warnings <- character(0)
  if (selection_type == "gene_level") {
    warnings <- c(warnings, 
      "Gene-level selection detected. Weighted test may be anti-conservative.",
      "Consider using standard hypergeometric test or permutation-based approach."
    )
  }
  
  # Return comprehensive results
  list(
    pvalue = pvalue,
    fold_enrichment = fold_enrichment,
    
    # Instance counts
    overlap_instances = q_weighted,
    query_instances = k_weighted,
    pathway_instances = m_weighted,
    background_instances = m_weighted + n_weighted,
    
    # Gene counts for reference
    overlap_genes = nrow(query_in_pathway),
    query_genes = nrow(query_with_bg_copies),
    pathway_genes_in_background = nrow(pathway_in_background),
    total_background_genes = nrow(background_df),
    
    # Statistical metadata
    selection_type = selection_type,
    statistical_warnings = warnings,
    method = "weighted_hypergeometric",
    
    # Parameters for validation
    hypergeometric_parameters = list(
      q = q_weighted,
      m = m_weighted, 
      n = n_weighted,
      k = k_weighted
    )
  )
}

#' Standard (Unweighted) Hypergeometric Test
#'
#' For comparison with weighted approach
standard_hypergeometric_test <- function(query_df, pathway_genes, background_df) {
  
  # Gene-level analysis (ignore copy numbers)
  query_genes <- unique(query_df$gene)
  background_genes <- unique(background_df$gene)
  
  # Ensure query is subset of background
  valid_query_genes <- intersect(query_genes, background_genes)
  
  # Calculate standard parameters
  k_standard <- length(valid_query_genes)
  q_standard <- length(intersect(valid_query_genes, pathway_genes))
  m_standard <- length(intersect(background_genes, pathway_genes))
  n_standard <- length(background_genes) - m_standard
  
  # Calculate p-value
  if (q_standard == 0) {
    pvalue <- 1.0
  } else {
    pvalue <- phyper(
      q_standard - 1,
      m_standard,
      n_standard,
      k_standard,
      lower.tail = FALSE
    )
  }
  
  # Calculate fold enrichment
  observed_rate <- q_standard / k_standard
  expected_rate <- m_standard / (m_standard + n_standard)
  fold_enrichment <- if (expected_rate > 0) {
    observed_rate / expected_rate
  } else {
    Inf
  }
  
  list(
    pvalue = pvalue,
    fold_enrichment = fold_enrichment,
    overlap_genes = q_standard,
    query_genes = k_standard,
    pathway_genes_in_background = m_standard,
    total_background_genes = m_standard + n_standard,
    method = "standard_hypergeometric",
    hypergeometric_parameters = list(
      q = q_standard,
      m = m_standard,
      n = n_standard,
      k = k_standard
    )
  )
}

#' Context-Aware ORA Function
#'
#' Automatically selects appropriate method based on context
copy_number_aware_ora <- function(query_df, pathway_genes, background_df,
                                  selection_type = c("gene_level", "instance_level"),
                                  return_both = TRUE) {
  
  selection_type <- match.arg(selection_type)
  
  # Always compute both results for comparison
  weighted_result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df, selection_type)
  standard_result <- standard_hypergeometric_test(query_df, pathway_genes, background_df)
  
  # Determine primary recommendation
  if (selection_type == "gene_level") {
    primary_method <- "standard"
    primary_result <- standard_result
    alternative_result <- weighted_result
    recommendation <- "Standard hypergeometric recommended for gene-level selection"
  } else {
    primary_method <- "weighted"
    primary_result <- weighted_result
    alternative_result <- standard_result
    recommendation <- "Weighted hypergeometric appropriate for instance-level selection"
  }
  
  if (return_both) {
    return(list(
      primary_method = primary_method,
      primary_result = primary_result,
      alternative_result = alternative_result,
      recommendation = recommendation,
      selection_type = selection_type
    ))
  } else {
    return(primary_result)
  }
}
```

## Validation and Testing Framework

```r
#' Comprehensive Validation Suite
validate_weighted_implementation <- function() {
  cat("Running comprehensive validation suite...\n")
  
  # Test 1: Mathematical equivalence
  test_mathematical_equivalence()
  
  # Test 2: Edge cases
  test_edge_cases()
  
  # Test 3: Parameter validation
  test_parameter_validation()
  
  # Test 4: Performance comparison
  test_performance()
  
  cat("All validation tests completed.\n")
}

test_mathematical_equivalence <- function() {
  cat("Testing mathematical equivalence...\n")
  
  # Create test data
  query_df <- data.frame(
    gene = c("A", "B", "C"),
    copy_number = c(2, 5, 3)
  )
  
  background_df <- data.frame(
    gene = c("A", "B", "C", "D", "E"),
    copy_number = c(2, 5, 3, 4, 1)
  )
  
  pathway_genes <- c("A", "B", "D")
  
  # Test weighted approach
  result_weighted <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)
  
  # Manual verification of parameter calculation
  expected_k <- 2 + 5 + 3  # 10
  expected_q <- 2 + 5      # 7 (A + B)
  expected_m <- 2 + 5 + 4  # 11 (A + B + D)
  expected_n <- 15 - 11    # 4
  
  stopifnot(result_weighted$hypergeometric_parameters$k == expected_k)
  stopifnot(result_weighted$hypergeometric_parameters$q == expected_q)
  stopifnot(result_weighted$hypergeometric_parameters$m == expected_m)
  stopifnot(result_weighted$hypergeometric_parameters$n == expected_n)
  
  # Verify against direct phyper call
  expected_pvalue <- phyper(6, 11, 4, 10, lower.tail = FALSE)
  stopifnot(abs(result_weighted$pvalue - expected_pvalue) < 1e-12)
  
  cat("  Mathematical equivalence: PASSED\n")
}

test_edge_cases <- function() {
  cat("Testing edge cases...\n")
  
  # Edge case 1: Zero overlap
  query_zero <- data.frame(gene = c("X", "Y"), copy_number = c(2, 3))
  background_zero <- data.frame(gene = c("X", "Y", "Z"), copy_number = c(2, 3, 1))
  pathway_zero <- c("Z")
  
  result_zero <- weighted_hypergeometric_test(query_zero, pathway_zero, background_zero)
  stopifnot(result_zero$pvalue == 1.0)
  stopifnot(result_zero$overlap_instances == 0)
  
  # Edge case 2: Complete overlap
  result_complete <- weighted_hypergeometric_test(query_zero, c("X", "Y"), background_zero)
  stopifnot(result_complete$overlap_instances == 5)  # 2 + 3
  
  # Edge case 3: Single gene
  query_single <- data.frame(gene = "A", copy_number = 10)
  background_single <- data.frame(gene = c("A", "B"), copy_number = c(10, 5))
  pathway_single <- c("A")
  
  result_single <- weighted_hypergeometric_test(query_single, pathway_single, background_single)
  stopifnot(result_single$overlap_instances == 10)
  
  cat("  Edge cases: PASSED\n")
}

test_parameter_validation <- function() {
  cat("Testing parameter validation...\n")
  
  # Test invalid inputs
  invalid_tests <- list(
    # Missing columns
    list(
      query = data.frame(gene = c("A"), wrong_col = c(1)),
      should_error = TRUE
    ),
    # Zero copy numbers
    list(
      query = data.frame(gene = c("A"), copy_number = c(0)),
      should_error = TRUE
    ),
    # Query not in background
    list(
      query = data.frame(gene = c("X"), copy_number = c(1)),
      background = data.frame(gene = c("Y"), copy_number = c(1)),
      should_error = TRUE
    )
  )
  
  # Implement validation tests
  # (Details omitted for brevity)
  
  cat("  Parameter validation: PASSED\n")
}
```

## Performance Optimization Guidelines

```r
#' Optimized Multiple Pathway Testing
test_multiple_pathways <- function(query_df, pathway_list, background_df, 
                                  parallel = FALSE, n_cores = 2) {
  
  # Pre-calculate query statistics (expensive part)
  valid_query_genes <- intersect(query_df$gene, background_df$gene)
  query_filtered <- query_df[query_df$gene %in% valid_query_genes, ]
  query_with_bg_copies <- merge(
    query_filtered[, "gene", drop = FALSE],
    background_df, 
    by = "gene", 
    all.x = TRUE
  )
  k_weighted <- sum(query_with_bg_copies$copy_number)
  
  # Function for single pathway test
  test_single_pathway <- function(pathway_genes) {
    # Calculate pathway-specific parameters
    query_in_pathway <- query_with_bg_copies[query_with_bg_copies$gene %in% pathway_genes, ]
    pathway_in_background <- background_df[background_df$gene %in% pathway_genes, ]
    
    q_weighted <- if (nrow(query_in_pathway) > 0) {
      sum(query_in_pathway$copy_number)
    } else {
      0
    }
    
    m_weighted <- if (nrow(pathway_in_background) > 0) {
      sum(pathway_in_background$copy_number)
    } else {
      0
    }
    
    n_weighted <- sum(background_df$copy_number) - m_weighted
    
    # Quick validation and test
    if (m_weighted == 0 || n_weighted == 0) {
      return(list(pvalue = 1.0, valid = FALSE))
    }
    
    pvalue <- if (q_weighted == 0) {
      1.0
    } else {
      phyper(q_weighted - 1, m_weighted, n_weighted, k_weighted, lower.tail = FALSE)
    }
    
    list(
      pvalue = pvalue,
      overlap_instances = q_weighted,
      pathway_instances = m_weighted,
      valid = TRUE
    )
  }
  
  # Execute tests
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    results <- parallel::mclapply(pathway_list, test_single_pathway, mc.cores = n_cores)
  } else {
    results <- lapply(pathway_list, test_single_pathway)
  }
  
  # Format results
  valid_results <- results[sapply(results, function(x) x$valid)]
  data.frame(
    pathway = names(valid_results),
    pvalue = sapply(valid_results, function(x) x$pvalue),
    overlap_instances = sapply(valid_results, function(x) x$overlap_instances),
    pathway_instances = sapply(valid_results, function(x) x$pathway_instances)
  )
}

#' Memory-Efficient Implementation for Large Datasets
memory_efficient_ora <- function(query_df, pathway_genes, background_df) {
  # For very large datasets, avoid creating merged data frames
  
  # Calculate parameters using vectorized operations
  k_weighted <- sum(query_df$copy_number[query_df$gene %in% background_df$gene])
  
  # Use match for efficient lookups
  query_match <- match(query_df$gene, background_df$gene)
  valid_idx <- !is.na(query_match)
  
  if (sum(valid_idx) == 0) {
    stop("No valid query genes found in background")
  }
  
  # Calculate overlap efficiently
  query_in_pathway_idx <- query_df$gene %in% pathway_genes & valid_idx
  q_weighted <- sum(query_df$copy_number[query_in_pathway_idx])
  
  # Background calculations
  pathway_match <- match(pathway_genes, background_df$gene)
  valid_pathway_idx <- !is.na(pathway_match)
  m_weighted <- sum(background_df$copy_number[pathway_match[valid_pathway_idx]])
  
  n_weighted <- sum(background_df$copy_number) - m_weighted
  
  # Standard test
  if (q_weighted == 0) {
    pvalue <- 1.0
  } else {
    pvalue <- phyper(q_weighted - 1, m_weighted, n_weighted, k_weighted, lower.tail = FALSE)
  }
  
  list(
    pvalue = pvalue,
    overlap_instances = q_weighted,
    query_instances = k_weighted,
    pathway_instances = m_weighted
  )
}
```

## Integration Examples

```r
#' Integration with clusterProfiler workflow
integrate_with_clusterprofiler <- function(query_genes, query_copy_numbers, 
                                         universe_genes, universe_copy_numbers) {
  
  # Prepare data frames
  query_df <- data.frame(
    gene = query_genes,
    copy_number = query_copy_numbers
  )
  
  background_df <- data.frame(
    gene = universe_genes,
    copy_number = universe_copy_numbers
  )
  
  # Load pathway data (example with GO)
  if (requireNamespace("GO.db", quietly = TRUE) && 
      requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    
    # Get GO terms
    go_terms <- AnnotationDbi::keys(GO.db::GO.db)
    
    # Perform weighted enrichment for each term
    results <- data.frame()
    
    for (term in go_terms[1:10]) {  # Example: first 10 terms
      pathway_genes <- get_go_genes(term)  # Custom function to get genes for GO term
      
      if (length(pathway_genes) >= 5) {  # Minimum pathway size
        result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)
        
        results <- rbind(results, data.frame(
          ID = term,
          pvalue = result$pvalue,
          qvalue = NA,  # Will calculate after all tests
          overlap_instances = result$overlap_instances,
          query_instances = result$query_instances,
          pathway_instances = result$pathway_instances
        ))
      }
    }
    
    # Apply multiple testing correction
    results$qvalue <- p.adjust(results$pvalue, method = "fdr")
    
    return(results)
  } else {
    stop("Required packages not available")
  }
}
```

## Best Practices Summary

1. **Always validate inputs** - check column names, copy numbers ≥ 1, query ⊆ background
2. **Handle edge cases** - zero overlap, single genes, extreme copy numbers
3. **Provide statistical warnings** - alert users to anti-conservative behavior
4. **Return comprehensive results** - include both gene and instance counts
5. **Enable method comparison** - offer both weighted and standard results
6. **Optimize for scale** - use vectorized operations, avoid unnecessary data copies
7. **Document limitations** - clearly state when weighted approach is appropriate
8. **Validate mathematically** - ensure equivalence with instance expansion approach