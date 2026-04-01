#!/usr/bin/env Rscript
#
# Empirical Verification of Copy-Number Weighted phyper() Mathematical Equivalence
#
# This script provides comprehensive empirical verification that parameter weighting
# and instance expansion methods for copy-number weighted hypergeometric testing
# yield mathematically equivalent results.
#
# Author: Generated for PHR copy-weighted ORA analysis
# Date: 2026-04-01

library(stats)

# Comprehensive verification function
verify_weighted_phyper_equivalence <- function(query_df, pathway_genes, background_df,
                                               tolerance = 1e-14) {

  # INPUT VALIDATION
  stopifnot(all(c("gene", "copy_number") %in% names(query_df)))
  stopifnot(all(c("gene", "copy_number") %in% names(background_df)))
  stopifnot(is.character(pathway_genes))

  # METHOD 1: INSTANCE EXPANSION
  cat("Method 1: Instance Expansion\n")
  start_time_exp <- Sys.time()

  # Create expanded vectors
  query_expanded <- rep(query_df$gene, query_df$copy_number)
  background_expanded <- rep(background_df$gene, background_df$copy_number)

  # Calculate expanded parameters
  k_exp <- length(query_expanded)
  q_exp <- sum(query_expanded %in% pathway_genes)
  m_exp <- sum(background_expanded %in% pathway_genes)
  n_exp <- length(background_expanded) - m_exp

  # Calculate p-value
  if (q_exp == 0) {
    pval_expansion <- 1.0
  } else {
    pval_expansion <- phyper(q_exp-1, m_exp, n_exp, k_exp, lower.tail = FALSE)
  }

  end_time_exp <- Sys.time()

  # METHOD 2: PARAMETER WEIGHTING
  cat("Method 2: Parameter Weighting\n")
  start_time_weight <- Sys.time()

  # Calculate weighted parameters
  k_weight <- sum(query_df$copy_number)

  query_pathway <- query_df[query_df$gene %in% pathway_genes, ]
  q_weight <- if(nrow(query_pathway) > 0) sum(query_pathway$copy_number) else 0

  background_pathway <- background_df[background_df$gene %in% pathway_genes, ]
  m_weight <- if(nrow(background_pathway) > 0) sum(background_pathway$copy_number) else 0

  n_weight <- sum(background_df$copy_number) - m_weight

  # Calculate p-value
  if (q_weight == 0) {
    pval_weighted <- 1.0
  } else {
    pval_weighted <- phyper(q_weight-1, m_weight, n_weight, k_weight, lower.tail = FALSE)
  }

  end_time_weight <- Sys.time()

  # VERIFICATION TESTS

  # Test 1: Parameter Equivalence
  params_identical <- all(c(
    k_exp == k_weight,
    q_exp == q_weight,
    m_exp == m_weight,
    n_exp == n_weight
  ))

  # Test 2: P-value Equivalence
  pval_diff <- abs(pval_expansion - pval_weighted)
  pvals_equivalent <- pval_diff < tolerance

  # Test 3: Performance Comparison
  time_expansion <- as.numeric(end_time_exp - start_time_exp)
  time_weighted <- as.numeric(end_time_weight - start_time_weight)
  speedup_factor <- time_expansion / time_weighted

  # COMPREHENSIVE RESULTS
  return(list(
    # Verification results
    mathematical_equivalence = params_identical,
    statistical_equivalence = pvals_equivalent,
    pvalue_difference = pval_diff,

    # Method 1 results
    expansion_parameters = list(k=k_exp, q=q_exp, m=m_exp, n=n_exp),
    expansion_pvalue = pval_expansion,
    expansion_time = time_expansion,

    # Method 2 results
    weighted_parameters = list(k=k_weight, q=q_weight, m=m_weight, n=n_weight),
    weighted_pvalue = pval_weighted,
    weighted_time = time_weighted,

    # Performance metrics
    speedup_factor = speedup_factor,
    memory_reduction_estimate = k_exp / nrow(query_df),

    # Statistical metrics
    fold_enrichment = ifelse(m_weight > 0,
                           (q_weight/k_weight) / (m_weight/(m_weight+n_weight)),
                           NA),
    expected_overlap = k_weight * m_weight / (m_weight + n_weight)
  ))
}

# Test Case 1: Simple Controlled Example
test_simple_case <- function() {
  cat("=== Test Case 1: Simple Controlled Example ===\n")

  # Create minimal test data
  query_df <- data.frame(
    gene = c("GENE1", "GENE2", "GENE3"),
    copy_number = c(2, 3, 1)
  )

  background_df <- data.frame(
    gene = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
    copy_number = c(2, 3, 1, 4, 2)
  )

  pathway_genes <- c("GENE1", "GENE2", "GENE4")

  # Manual calculation for verification
  # Expected: k=6, q=5, m=9, n=3

  result <- verify_weighted_phyper_equivalence(query_df, pathway_genes, background_df)

  cat("Parameters match:", result$mathematical_equivalence, "\n")
  cat("P-values match:", result$statistical_equivalence, "\n")
  cat("P-value difference:", result$pvalue_difference, "\n")
  cat("Speedup factor:", round(result$speedup_factor, 2), "x\n")

  return(result)
}

# Test Case 2: PHR-Realistic Data Scale
test_phr_scale <- function() {
  cat("=== Test Case 2: PHR-Scale Realistic Data ===\n")

  # Simulate PHR-like data: 35 genes, ~1200 total copies
  set.seed(42)  # Reproducible

  query_genes <- paste0("PHR_", 1:35)
  query_copies <- sample(10:50, 35, replace=TRUE)  # High copy variation like PHRs
  query_copies[1:4] <- c(14, 14, 14, 14)  # OR genes with specific copies

  query_df <- data.frame(
    gene = query_genes,
    copy_number = query_copies
  )

  # Simulate genome background: 20,000 genes
  background_genes <- c(query_genes, paste0("BACKGROUND_", 1:19965))
  background_copies <- c(query_copies, sample(1:5, 19965, replace=TRUE))

  background_df <- data.frame(
    gene = background_genes,
    copy_number = background_copies
  )

  # Olfactory receptor pathway (includes first 4 PHR genes)
  or_pathway <- c(query_genes[1:4], paste0("OR_", 1:396))  # 400 OR genes total

  result <- verify_weighted_phyper_equivalence(query_df, or_pathway, background_df)

  cat("Dataset size: Query =", nrow(query_df), "genes,", sum(query_df$copy_number), "instances\n")
  cat("              Background =", nrow(background_df), "genes,", sum(background_df$copy_number), "instances\n")
  cat("              Pathway =", length(or_pathway), "genes\n")
  cat("Parameters match:", result$mathematical_equivalence, "\n")
  cat("P-values match:", result$statistical_equivalence, "\n")
  cat("P-value difference:", result$pvalue_difference, "\n")
  cat("Speedup factor:", round(result$speedup_factor, 2), "x\n")
  cat("P-value:", format(result$weighted_pvalue, scientific=TRUE), "\n")

  return(result)
}

# Test Case 3: Edge Cases
test_edge_cases <- function() {
  cat("=== Test Case 3: Edge Cases ===\n")

  edge_results <- list()

  # Edge Case 3a: No overlap
  cat("\nEdge Case 3a: Zero Overlap\n")
  query_df <- data.frame(gene = c("A", "B"), copy_number = c(5, 3))
  background_df <- data.frame(gene = c("A", "B", "C", "D"), copy_number = c(5, 3, 2, 1))
  pathway_genes <- c("C", "D")  # No overlap with query

  result_3a <- verify_weighted_phyper_equivalence(query_df, pathway_genes, background_df)
  edge_results$zero_overlap <- result_3a
  cat("  P-value:", result_3a$weighted_pvalue, "(should be 1.0)\n")
  cat("  Equivalence:", result_3a$mathematical_equivalence, "\n")

  # Edge Case 3b: Complete overlap
  cat("\nEdge Case 3b: Complete Overlap\n")
  pathway_genes <- c("A", "B")  # Complete overlap

  result_3b <- verify_weighted_phyper_equivalence(query_df, pathway_genes, background_df)
  edge_results$complete_overlap <- result_3b
  cat("  P-value:", format(result_3b$weighted_pvalue, scientific=TRUE), "\n")
  cat("  Equivalence:", result_3b$mathematical_equivalence, "\n")

  # Edge Case 3c: Extreme copy numbers
  cat("\nEdge Case 3c: Extreme Copy Numbers\n")
  query_extreme <- data.frame(gene = c("HIGH1", "HIGH2"), copy_number = c(500, 300))
  background_extreme <- data.frame(
    gene = c("HIGH1", "HIGH2", "LOW1", "LOW2"),
    copy_number = c(500, 300, 1, 1)
  )
  pathway_extreme <- c("HIGH1", "LOW1")

  result_3c <- verify_weighted_phyper_equivalence(query_extreme, pathway_extreme, background_extreme)
  edge_results$extreme_copies <- result_3c
  cat("  P-value:", format(result_3c$weighted_pvalue, scientific=TRUE), "\n")
  cat("  Equivalence:", result_3c$mathematical_equivalence, "\n")
  cat("  Speedup factor:", round(result_3c$speedup_factor, 2), "x\n")

  # Edge Case 3d: Single gene query
  cat("\nEdge Case 3d: Single Gene Query\n")
  single_query <- data.frame(gene = "SINGLE", copy_number = 10)
  single_background <- data.frame(
    gene = c("SINGLE", paste0("BG_", 1:100)),
    copy_number = c(10, rep(2, 100))
  )
  single_pathway <- c("SINGLE", paste0("BG_", 1:20))

  result_3d <- verify_weighted_phyper_equivalence(single_query, single_pathway, single_background)
  edge_results$single_gene <- result_3d
  cat("  P-value:", format(result_3d$weighted_pvalue, scientific=TRUE), "\n")
  cat("  Equivalence:", result_3d$mathematical_equivalence, "\n")

  return(edge_results)
}

# Test numerical precision across different scales
test_numerical_precision <- function() {
  cat("=== Numerical Precision Analysis ===\n")

  precision_results <- data.frame()

  # Test across different scales
  scales <- list(
    small = list(genes=5, max_copies=3),
    medium = list(genes=50, max_copies=20),
    large = list(genes=500, max_copies=100),
    extreme = list(genes=1000, max_copies=1000)
  )

  for (scale_name in names(scales)) {
    scale <- scales[[scale_name]]
    cat("\nTesting", scale_name, "scale:", scale$genes, "genes, max", scale$max_copies, "copies\n")

    # Generate test data
    set.seed(123)
    n_genes <- scale$genes
    query_genes <- paste0("G", 1:n_genes)
    query_copies <- sample(1:scale$max_copies, n_genes, replace=TRUE)

    query_df <- data.frame(gene = query_genes, copy_number = query_copies)

    # Background: 10x larger
    bg_genes <- c(query_genes, paste0("BG", 1:(n_genes*10)))
    bg_copies <- c(query_copies, sample(1:scale$max_copies, n_genes*10, replace=TRUE))

    background_df <- data.frame(gene = bg_genes, copy_number = bg_copies)

    # Pathway: ~30% of genes
    pathway_size <- max(1, round(n_genes * 0.3))
    pathway_genes <- sample(query_genes, pathway_size)

    # Run verification
    result <- verify_weighted_phyper_equivalence(query_df, pathway_genes, background_df,
                                                tolerance = .Machine$double.eps * 100)

    # Record results
    precision_results <- rbind(precision_results, data.frame(
      scale = scale_name,
      n_genes = n_genes,
      max_copies = scale$max_copies,
      total_instances = sum(query_df$copy_number),
      params_match = result$mathematical_equivalence,
      pvals_match = result$statistical_equivalence,
      pval_diff = result$pvalue_difference,
      speedup = result$speedup_factor
    ))

    cat("  Parameters match:", result$mathematical_equivalence, "\n")
    cat("  P-values match:", result$statistical_equivalence, "\n")
    cat("  Difference:", format(result$pvalue_difference, scientific=TRUE), "\n")
    cat("  Speedup:", round(result$speedup_factor, 2), "x\n")
  }

  return(precision_results)
}

# Analyze behavior near machine precision limits
test_machine_precision <- function() {
  cat("=== Machine Precision Limits ===\n")

  cat("Machine precision (.Machine$double.eps):", .Machine$double.eps, "\n")
  cat("Largest finite double (.Machine$double.xmax):", .Machine$double.xmax, "\n")
  cat("Smallest positive double (.Machine$double.xmin):", .Machine$double.xmin, "\n")

  # Test with parameters that might cause precision issues
  extreme_cases <- list(

    # Very small p-values (highly significant)
    high_significance = list(
      query = data.frame(gene = paste0("SIG", 1:20), copy_number = rep(50, 20)),
      background = data.frame(gene = paste0("ALL", 1:1000), copy_number = rep(2, 1000)),
      pathway = paste0("SIG", 1:20)  # All query genes in pathway
    ),

    # Very large parameters
    large_numbers = list(
      query = data.frame(gene = paste0("BIG", 1:10), copy_number = rep(1000, 10)),
      background = data.frame(gene = paste0("HUGE", 1:100), copy_number = rep(1000, 100)),
      pathway = paste0("BIG", 1:5)
    )
  )

  precision_limits <- list()

  for (case_name in names(extreme_cases)) {
    cat("\nTesting", case_name, "case:\n")
    case_data <- extreme_cases[[case_name]]

    tryCatch({
      result <- verify_weighted_phyper_equivalence(
        case_data$query, case_data$pathway, case_data$background,
        tolerance = .Machine$double.eps * 10
      )

      precision_limits[[case_name]] <- result

      cat("  Success: Equivalence =", result$mathematical_equivalence, "\n")
      cat("  P-value =", format(result$weighted_pvalue, scientific=TRUE), "\n")
      cat("  Difference =", format(result$pvalue_difference, scientific=TRUE), "\n")

    }, error = function(e) {
      cat("  Error encountered:", e$message, "\n")
      precision_limits[[case_name]] <<- list(error = e$message)
    })
  }

  return(precision_limits)
}

# Main verification runner
run_comprehensive_verification <- function() {
  cat("######################################################\n")
  cat("# COMPREHENSIVE WEIGHTED PHYPER() VERIFICATION SUITE #\n")
  cat("######################################################\n\n")

  # Initialize results container
  verification_results <- list()

  # Run all test suites
  verification_results$simple <- test_simple_case()
  verification_results$phr_scale <- test_phr_scale()
  verification_results$edge_cases <- test_edge_cases()
  verification_results$precision <- test_numerical_precision()
  verification_results$machine_limits <- test_machine_precision()

  # Summary analysis
  cat("\n######################################################\n")
  cat("# VERIFICATION SUMMARY\n")
  cat("######################################################\n")

  # Count successes
  all_param_matches <- c(
    verification_results$simple$mathematical_equivalence,
    verification_results$phr_scale$mathematical_equivalence,
    sapply(verification_results$edge_cases, function(x) x$mathematical_equivalence),
    verification_results$precision$params_match
  )

  all_pval_matches <- c(
    verification_results$simple$statistical_equivalence,
    verification_results$phr_scale$statistical_equivalence,
    sapply(verification_results$edge_cases, function(x) x$statistical_equivalence),
    verification_results$precision$pvals_match
  )

  param_success_rate <- mean(all_param_matches, na.rm=TRUE) * 100
  pval_success_rate <- mean(all_pval_matches, na.rm=TRUE) * 100

  cat("Parameter Equivalence Success Rate:", round(param_success_rate, 1), "%\n")
  cat("P-value Equivalence Success Rate:", round(pval_success_rate, 1), "%\n")

  # Performance summary
  speedups <- c(
    verification_results$simple$speedup_factor,
    verification_results$phr_scale$speedup_factor,
    sapply(verification_results$edge_cases, function(x) x$speedup_factor)
  )

  cat("Average Speedup Factor:", round(mean(speedups, na.rm=TRUE), 2), "x\n")
  cat("Max Speedup Factor:", round(max(speedups, na.rm=TRUE), 2), "x\n")

  # Final verdict
  if (param_success_rate == 100 && pval_success_rate == 100) {
    cat("\n✅ VERIFICATION PASSED: Mathematical equivalence confirmed\n")
  } else {
    cat("\n❌ VERIFICATION FAILED: Equivalence issues detected\n")
  }

  return(verification_results)
}

# Run the verification if script is executed directly
if (!interactive()) {
  cat("Starting comprehensive verification of weighted phyper() equivalence...\n\n")

  # Run comprehensive verification
  results <- run_comprehensive_verification()

  # Save results to file for later analysis
  saveRDS(results, "weighted_phyper_verification_results.rds")
  cat("\nResults saved to: weighted_phyper_verification_results.rds\n")
}