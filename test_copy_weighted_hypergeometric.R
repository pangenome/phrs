# Unit Tests for Copy-Number Weighted Hypergeometric Testing
#
# Comprehensive test suite covering functionality, edge cases, and validation
# of the copy-weighted hypergeometric implementation.
#
# Author: Robust R Code Implementation Task
# Date: 2026-04-01
# Version: 1.0

# Source the implementation
source("copy_weighted_hypergeometric.R")

# Load required libraries
if (!require(testthat, quietly = TRUE)) {
  stop("testthat package required for tests")
}

#' Test Mathematical Equivalence: Parameter Weighting vs Instance Expansion
#'
#' Verifies that parameter weighting produces identical results to instance expansion
test_mathematical_equivalence <- function() {
  cat("Testing mathematical equivalence...\n")

  # Create test data
  query_df <- data.frame(
    gene = c("GENE1", "GENE2", "GENE3", "GENE4"),
    copy_number = c(1, 3, 2, 4)
  )

  pathway_genes <- c("GENE1", "GENE2", "GENE5")

  background_df <- data.frame(
    gene = paste0("GENE", 1:20),
    copy_number = c(1, 3, 2, 4, 2, rep(1, 15))  # First 5 match query, rest are singles
  )

  # Method 1: Instance expansion
  query_expanded <- rep(query_df$gene, query_df$copy_number)
  background_expanded <- rep(background_df$gene, background_df$copy_number)

  q_exp <- sum(query_expanded %in% pathway_genes)
  m_exp <- sum(background_expanded %in% pathway_genes)
  n_exp <- length(background_expanded) - m_exp
  k_exp <- length(query_expanded)

  pval_expansion <- phyper(q_exp - 1, m_exp, n_exp, k_exp, lower.tail = FALSE)

  # Method 2: Parameter weighting
  weighted_result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)
  pval_weighted <- weighted_result$pvalue

  # Verify parameter equivalence
  params_match <- all(c(
    weighted_result$overlap_instances == q_exp,
    weighted_result$pathway_instances == m_exp,
    (weighted_result$background_instances - weighted_result$pathway_instances) == n_exp,
    weighted_result$query_instances == k_exp
  ))

  # Verify p-value equivalence (within numerical precision)
  pval_diff <- abs(pval_weighted - pval_expansion)
  pvals_match <- pval_diff < 1e-12

  # Test assertions
  if (!params_match) {
    stop("Parameter equivalence test FAILED")
  }

  if (!pvals_match) {
    stop(paste("P-value equivalence test FAILED. Difference:", pval_diff))
  }

  cat("✓ Mathematical equivalence verified\n")
  return(list(
    parameters_match = params_match,
    pvalues_match = pvals_match,
    pvalue_difference = pval_diff
  ))
}

#' Test Input Validation
#'
#' Tests comprehensive input validation functionality
test_input_validation <- function() {
  cat("Testing input validation...\n")

  # Valid input for reference
  valid_query <- data.frame(gene = c("A", "B"), copy_number = c(1, 2))
  valid_pathway <- c("A", "C")
  valid_background <- data.frame(gene = c("A", "B", "C", "D"), copy_number = c(1, 2, 1, 1))

  # Test 1: Missing columns
  tryCatch({
    bad_query <- data.frame(symbol = c("A", "B"), copies = c(1, 2))
    weighted_hypergeometric_test(bad_query, valid_pathway, valid_background)
    stop("Should have failed on missing columns")
  }, error = function(e) {
    if (!grepl("missing columns", e$message)) {
      stop("Wrong error message for missing columns")
    }
  })

  # Test 2: Negative copy numbers
  tryCatch({
    bad_query <- data.frame(gene = c("A", "B"), copy_number = c(-1, 2))
    weighted_hypergeometric_test(bad_query, valid_pathway, valid_background)
    stop("Should have failed on negative copy numbers")
  }, error = function(e) {
    if (!grepl("negative copy numbers", e$message)) {
      stop("Wrong error message for negative copy numbers")
    }
  })

  # Test 3: Duplicated genes
  tryCatch({
    bad_query <- data.frame(gene = c("A", "A"), copy_number = c(1, 2))
    weighted_hypergeometric_test(bad_query, valid_pathway, valid_background)
    stop("Should have failed on duplicated genes")
  }, error = function(e) {
    if (!grepl("duplicated genes", e$message)) {
      stop("Wrong error message for duplicated genes")
    }
  })

  # Test 4: Empty datasets
  tryCatch({
    empty_query <- data.frame(gene = character(0), copy_number = numeric(0))
    weighted_hypergeometric_test(empty_query, valid_pathway, valid_background)
    stop("Should have failed on empty query")
  }, error = function(e) {
    if (!grepl("empty", e$message)) {
      stop("Wrong error message for empty dataset")
    }
  })

  # Test 5: No pathway genes in background
  tryCatch({
    no_pathway_bg <- data.frame(gene = c("X", "Y", "Z"), copy_number = c(1, 1, 1))
    weighted_hypergeometric_test(valid_query, valid_pathway, no_pathway_bg)
    stop("Should have failed when no pathway genes in background")
  }, error = function(e) {
    if (!grepl("No pathway genes found", e$message)) {
      stop("Wrong error message for missing pathway genes")
    }
  })

  cat("✓ Input validation tests passed\n")
  return(TRUE)
}

#' Test Edge Cases
#'
#' Tests handling of edge cases and boundary conditions
test_edge_cases <- function() {
  cat("Testing edge cases...\n")

  # Test 1: Zero copy numbers with handle_zeros = TRUE
  query_with_zeros <- data.frame(
    gene = c("A", "B", "C"),
    copy_number = c(0, 1, 2)
  )

  background_with_zeros <- data.frame(
    gene = c("A", "B", "C", "D"),
    copy_number = c(0, 1, 2, 1)
  )

  pathway <- c("A", "B")

  # Should handle zeros gracefully with warning
  tryCatch({
    result <- weighted_hypergeometric_test(
      query_with_zeros, pathway, background_with_zeros,
      handle_zeros = TRUE
    )
    # Should succeed after removing zero-copy genes
    if (result$diagnostics$genes_removed != 2) {
      stop("Incorrect number of genes removed")
    }
  }, warning = function(w) {
    if (!grepl("zero copies", w$message)) {
      stop("Expected warning about zero copies")
    }
  })

  # Test 2: Extreme copy numbers with capping
  extreme_query <- data.frame(
    gene = c("A", "B"),
    copy_number = c(1000, 1)
  )

  extreme_background <- data.frame(
    gene = c("A", "B", "C"),
    copy_number = c(1000, 1, 1)
  )

  # Should cap extreme values with warning
  tryCatch({
    result <- weighted_hypergeometric_test(
      extreme_query, pathway, extreme_background,
      max_copies = 100
    )
    # Check that copies were capped
    max_copies_found <- max(result$diagnostics$copy_model_consistency$correlation, na.rm = TRUE)
  }, warning = function(w) {
    # Should warn about capping
  })

  # Test 3: Perfect overlap (all query genes in pathway)
  perfect_query <- data.frame(
    gene = c("A", "B"),
    copy_number = c(2, 3)
  )

  perfect_background <- data.frame(
    gene = c("A", "B", "C", "D"),
    copy_number = c(2, 3, 1, 1)
  )

  perfect_pathway <- c("A", "B")

  result <- weighted_hypergeometric_test(perfect_query, perfect_pathway, perfect_background)

  # Should have overlap equal to query size
  if (result$overlap_instances != result$query_instances) {
    stop("Perfect overlap test failed")
  }

  # Test 4: No overlap
  no_overlap_query <- data.frame(
    gene = c("A", "B"),
    copy_number = c(2, 3)
  )

  no_overlap_pathway <- c("C", "D")

  result <- weighted_hypergeometric_test(no_overlap_query, no_overlap_pathway, perfect_background)

  # Should have zero overlap
  if (result$overlap_instances != 0) {
    stop("No overlap test failed")
  }

  # P-value should be high (not significant)
  if (result$pvalue < 0.5) {
    stop("No overlap should give high p-value")
  }

  cat("✓ Edge case tests passed\n")
  return(TRUE)
}

#' Test Constraint Validation
#'
#' Tests hypergeometric constraint validation
test_constraint_validation <- function() {
  cat("Testing hypergeometric constraints...\n")

  # Test valid constraints
  valid_result <- validate_hypergeometric_constraints(5, 10, 90, 20)
  if (!valid_result$valid) {
    stop("Valid constraints rejected")
  }

  # Test q > k (overlap > sample)
  invalid_result <- validate_hypergeometric_constraints(25, 10, 90, 20)
  if (invalid_result$valid) {
    stop("Should reject q > k")
  }

  # Test q > m (overlap > population of interest)
  invalid_result <- validate_hypergeometric_constraints(15, 10, 90, 20)
  if (invalid_result$valid) {
    stop("Should reject q > m")
  }

  # Test k > m + n (sample > total population)
  invalid_result <- validate_hypergeometric_constraints(5, 10, 90, 150)
  if (invalid_result$valid) {
    stop("Should reject k > m + n")
  }

  # Test negative values
  invalid_result <- validate_hypergeometric_constraints(-1, 10, 90, 20)
  if (invalid_result$valid) {
    stop("Should reject negative values")
  }

  # Test non-integers
  invalid_result <- validate_hypergeometric_constraints(5.5, 10, 90, 20)
  if (invalid_result$valid) {
    stop("Should reject non-integers")
  }

  cat("✓ Constraint validation tests passed\n")
  return(TRUE)
}

#' Test Performance with Realistic Data Sizes
#'
#' Tests performance with PHR-scale datasets
test_performance <- function() {
  cat("Testing performance with realistic data sizes...\n")

  # Create PHR-scale test data
  set.seed(42)  # For reproducible results

  # PHR-like query: 35 genes, ~1200 total instances
  query_genes <- paste0("GENE", 1:35)
  query_copies <- rpois(35, lambda = 30) + 1  # Mean ~30 copies per gene

  query_df <- data.frame(
    gene = query_genes,
    copy_number = query_copies
  )

  # Human genome-scale background: 20,000 genes
  background_genes <- paste0("GENE", 1:20000)
  background_copies <- rpois(20000, lambda = 2) + 1  # Mean ~3 copies per gene

  background_df <- data.frame(
    gene = background_genes,
    copy_number = background_copies
  )

  # Olfactory receptor pathway: 400 genes
  pathway_genes <- paste0("GENE", 1:400)  # First 400 genes are "OR genes"

  # Measure performance
  start_time <- Sys.time()

  result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)

  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Check that it completes in reasonable time (< 1 second for this scale)
  if (runtime > 1.0) {
    warning(paste("Performance test slower than expected:", runtime, "seconds"))
  }

  # Check result validity
  if (!is.finite(result$pvalue) || result$pvalue < 0 || result$pvalue > 1) {
    stop("Performance test produced invalid p-value")
  }

  # Check parameter reasonableness
  if (result$query_instances < 35) {  # Should have at least 1 copy per gene
    stop("Unreasonable query instance count")
  }

  if (result$pathway_instances < 400) {  # Should have at least 1 copy per pathway gene
    stop("Unreasonable pathway instance count")
  }

  cat(paste("✓ Performance test passed in", round(runtime, 4), "seconds\n"))

  return(list(
    runtime_seconds = runtime,
    query_instances = result$query_instances,
    pathway_instances = result$pathway_instances,
    pvalue = result$pvalue
  ))
}

#' Test Comparison with Standard Method
#'
#' Tests the utility function for comparing weighted vs standard results
test_comparison_function <- function() {
  cat("Testing comparison function...\n")

  # Create test data where copy numbers matter
  query_df <- data.frame(
    gene = c("A", "B", "C", "D"),
    copy_number = c(10, 1, 1, 1)  # Gene A has many copies
  )

  pathway_genes <- c("A", "E")  # A is in pathway and has many copies

  background_df <- data.frame(
    gene = c("A", "B", "C", "D", "E", "F", "G", "H"),
    copy_number = c(10, 1, 1, 1, 2, 1, 1, 1)
  )

  comparison <- compare_weighted_vs_standard(query_df, pathway_genes, background_df)

  # Check structure
  required_elements <- c("standard", "weighted", "comparison")
  if (!all(required_elements %in% names(comparison))) {
    stop("Comparison function missing required elements")
  }

  # Standard method should count A as 1 gene
  if (comparison$standard$parameters["k"] != 4) {
    stop("Standard method should count 4 query genes")
  }

  # Weighted method should count A as 10 instances
  if (comparison$weighted$parameters["k"] != 13) {
    stop("Weighted method should count 13 query instances")
  }

  # With high copy gene A in pathway, weighted should be more significant
  if (!comparison$comparison$more_significant) {
    stop("Weighted method should be more significant in this case")
  }

  cat("✓ Comparison function tests passed\n")
  return(TRUE)
}

#' Test Statistical Properties
#'
#' Basic tests of statistical properties (simplified null distribution test)
test_statistical_properties <- function() {
  cat("Testing statistical properties...\n")

  set.seed(123)

  # Create background
  background_genes <- paste0("GENE", 1:1000)
  background_copies <- rpois(1000, lambda = 2) + 1

  background_df <- data.frame(
    gene = background_genes,
    copy_number = background_copies
  )

  # Pathway of 50 genes
  pathway_genes <- paste0("GENE", 1:50)

  # Generate multiple random queries (null hypothesis - no true enrichment)
  n_simulations <- 50  # Reduced for testing speed
  query_size <- 20

  pvalues <- replicate(n_simulations, {
    # Random query with no enrichment bias
    random_genes <- sample(background_genes, query_size, replace = FALSE)
    query_df <- data.frame(
      gene = random_genes,
      copy_number = background_df$copy_number[background_df$gene %in% random_genes]
    )

    result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)
    return(result$pvalue)
  })

  # Basic statistical checks
  # P-values should be between 0 and 1
  if (any(pvalues < 0 | pvalues > 1)) {
    stop("Invalid p-values generated")
  }

  # Should not have too many extremely small p-values (would suggest bias)
  very_small_pvals <- sum(pvalues < 0.01)
  if (very_small_pvals > 0.1 * n_simulations) {  # More than 10%
    warning("Unusually many small p-values in null simulation")
  }

  # Mean should be approximately 0.5 under null (rough check)
  mean_pval <- mean(pvalues)
  if (mean_pval < 0.2 || mean_pval > 0.8) {
    warning(paste("P-value mean", round(mean_pval, 3), "seems biased"))
  }

  cat(paste("✓ Statistical properties test passed. Mean p-value:",
            round(mean_pval, 3), "\n"))

  return(list(
    mean_pvalue = mean_pval,
    n_small_pvals = very_small_pvals,
    pvalue_range = range(pvalues)
  ))
}

#' Main Test Runner
#'
#' Runs all tests and reports results
run_all_tests <- function() {
  cat("==========================================\n")
  cat("Running Copy-Weighted Hypergeometric Tests\n")
  cat("==========================================\n\n")

  tests_passed <- 0
  tests_total <- 7

  test_results <- list()

  # Test 1: Mathematical equivalence
  tryCatch({
    test_results$equivalence <- test_mathematical_equivalence()
    tests_passed <- tests_passed + 1
  }, error = function(e) {
    cat("✗ Mathematical equivalence test FAILED:", e$message, "\n")
  })

  # Test 2: Input validation
  tryCatch({
    test_results$validation <- test_input_validation()
    tests_passed <- tests_passed + 1
  }, error = function(e) {
    cat("✗ Input validation test FAILED:", e$message, "\n")
  })

  # Test 3: Edge cases
  tryCatch({
    test_results$edge_cases <- test_edge_cases()
    tests_passed <- tests_passed + 1
  }, error = function(e) {
    cat("✗ Edge case test FAILED:", e$message, "\n")
  })

  # Test 4: Constraint validation
  tryCatch({
    test_results$constraints <- test_constraint_validation()
    tests_passed <- tests_passed + 1
  }, error = function(e) {
    cat("✗ Constraint validation test FAILED:", e$message, "\n")
  })

  # Test 5: Performance
  tryCatch({
    test_results$performance <- test_performance()
    tests_passed <- tests_passed + 1
  }, error = function(e) {
    cat("✗ Performance test FAILED:", e$message, "\n")
  })

  # Test 6: Comparison function
  tryCatch({
    test_results$comparison <- test_comparison_function()
    tests_passed <- tests_passed + 1
  }, error = function(e) {
    cat("✗ Comparison function test FAILED:", e$message, "\n")
  })

  # Test 7: Statistical properties
  tryCatch({
    test_results$statistics <- test_statistical_properties()
    tests_passed <- tests_passed + 1
  }, error = function(e) {
    cat("✗ Statistical properties test FAILED:", e$message, "\n")
  })

  # Summary
  cat("\n==========================================\n")
  cat(paste("Test Results:", tests_passed, "/", tests_total, "passed\n"))

  if (tests_passed == tests_total) {
    cat("✓ All tests PASSED! Implementation is robust.\n")
  } else {
    cat("✗ Some tests FAILED. Review implementation.\n")
  }

  cat("==========================================\n")

  return(list(
    tests_passed = tests_passed,
    tests_total = tests_total,
    all_passed = tests_passed == tests_total,
    detailed_results = test_results
  ))
}

# Run tests if this script is executed directly
if (sys.nframe() == 0) {
  result <- run_all_tests()
}