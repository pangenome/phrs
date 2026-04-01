# Type I Error Rate Validation Functions for Copy-Number Weighted phyper()
# Comprehensive validation of false positive rates and error control
#
# Author: Generated for phrs project
# Task: type-i-error
#
# This file implements comprehensive Type I error control validation for
# copy-number weighted phyper() parameters as specified in task requirements:
# - Test false positive rates at multiple α levels (0.01, 0.05, 0.1)
# - Validate error rate control across different copy number distributions
# - Test with varying background and pathway sizes
# - Ensure error rates stay within 1% of nominal levels

library(tidyverse)
source("debug_weighted_phyper.R")

# ============================================================
# CORE TYPE I ERROR VALIDATION FUNCTIONS
# ============================================================

#' Validate Type I error rates at multiple alpha levels
#'
#' Tests false positive rates for both weighted and standard hypergeometric tests
#' under true null hypothesis (no enrichment). Validates that error rates stay
#' within acceptable bounds of nominal alpha levels.
#'
#' @param n_sims Number of simulations to run (default: 2000)
#' @param alphas Vector of alpha levels to test (default: c(0.01, 0.05, 0.1))
#' @param n_genes Number of background genes (default: 1000)
#' @param n_pathway Number of pathway genes (default: 100)
#' @param n_query Number of query genes (default: 50)
#' @param copy_distribution Function to generate copy numbers (default: sample(1:8, n, replace=TRUE))
#' @param tolerance Acceptable deviation from nominal alpha (default: 0.01 = 1%)
#' @param verbose Whether to print detailed results (default: TRUE)
#' @return List with validation results and pass/fail status
validate_type_i_error_rates <- function(n_sims = 2000,
                                        alphas = c(0.01, 0.05, 0.1),
                                        n_genes = 1000,
                                        n_pathway = 100,
                                        n_query = 50,
                                        copy_distribution = function(n) sample(1:8, n, replace = TRUE),
                                        tolerance = 0.01,
                                        verbose = TRUE) {

  if (verbose) {
    cat("============================================\n")
    cat("TYPE I ERROR RATE VALIDATION\n")
    cat("============================================\n")
    cat(sprintf("Simulations: %d\n", n_sims))
    cat(sprintf("Background genes: %d, Pathway genes: %d, Query genes: %d\n",
                n_genes, n_pathway, n_query))
    cat(sprintf("Alpha levels: %s\n", paste(alphas, collapse = ", ")))
    cat(sprintf("Tolerance: ±%.1f%%\n\n", tolerance * 100))
  }

  set.seed(42)  # Reproducible results

  # Generate background dataset
  background_df <- data.frame(
    gene = paste0("GENE", 1:n_genes),
    copy_number = copy_distribution(n_genes),
    stringsAsFactors = FALSE
  )

  # Random pathway (null hypothesis: no true enrichment)
  pathway_genes <- sample(background_df$gene, n_pathway)

  # Run simulations under null hypothesis
  if (verbose) cat("Running null hypothesis simulations...\n")

  # Weighted phyper p-values
  weighted_pvals <- replicate(n_sims, {
    query_genes <- sample(background_df$gene, n_query)
    query_df <- background_df[background_df$gene %in% query_genes, ]
    weighted_hypergeometric_test_fixed(query_df, pathway_genes, background_df)$pvalue
  })

  # Standard phyper p-values
  standard_pvals <- replicate(n_sims, {
    query_genes <- sample(background_df$gene, n_query)
    overlap <- sum(query_genes %in% pathway_genes)
    phyper(overlap - 1, n_pathway, n_genes - n_pathway, n_query, lower.tail = FALSE)
  })

  # Calculate error rates for each alpha level
  results <- data.frame(
    alpha = alphas,
    weighted_rate = sapply(alphas, function(a) mean(weighted_pvals < a)),
    standard_rate = sapply(alphas, function(a) mean(standard_pvals < a)),
    stringsAsFactors = FALSE
  )

  # Calculate confidence intervals
  results$weighted_ci_lower <- sapply(1:nrow(results), function(i) {
    binom.test(sum(weighted_pvals < results$alpha[i]), n_sims)$conf.int[1]
  })
  results$weighted_ci_upper <- sapply(1:nrow(results), function(i) {
    binom.test(sum(weighted_pvals < results$alpha[i]), n_sims)$conf.int[2]
  })
  results$standard_ci_lower <- sapply(1:nrow(results), function(i) {
    binom.test(sum(standard_pvals < results$alpha[i]), n_sims)$conf.int[1]
  })
  results$standard_ci_upper <- sapply(1:nrow(results), function(i) {
    binom.test(sum(standard_pvals < results$alpha[i]), n_sims)$conf.int[2]
  })

  # Check if rates are within tolerance
  results$weighted_controlled <- abs(results$weighted_rate - results$alpha) <= tolerance
  results$standard_controlled <- abs(results$standard_rate - results$alpha) <= tolerance

  # Print results
  if (verbose) {
    cat(sprintf("\n%-8s %-12s %-12s %-15s %-15s\n",
                "Alpha", "Weighted", "Standard", "W_Controlled", "S_Controlled"))
    cat(sprintf("%-8s %-12s %-12s %-15s %-15s\n",
                "------", "--------", "--------", "-----------", "-----------"))

    for (i in 1:nrow(results)) {
      cat(sprintf("%-8.3f %-12.4f %-12.4f %-15s %-15s\n",
                  results$alpha[i],
                  results$weighted_rate[i],
                  results$standard_rate[i],
                  results$weighted_controlled[i],
                  results$standard_controlled[i]))
    }

    cat("\nWith 95% confidence intervals:\n")
    for (i in 1:nrow(results)) {
      cat(sprintf("α=%.3f: Weighted: %.4f (%.3f-%.3f), Standard: %.4f (%.3f-%.3f)\n",
                  results$alpha[i],
                  results$weighted_rate[i], results$weighted_ci_lower[i], results$weighted_ci_upper[i],
                  results$standard_rate[i], results$standard_ci_lower[i], results$standard_ci_upper[i]))
    }

    # Overall assessment
    all_weighted_controlled <- all(results$weighted_controlled)
    all_standard_controlled <- all(results$standard_controlled)

    cat(sprintf("\nOVERALL ASSESSMENT:\n"))
    cat(sprintf("Weighted phyper Type I error control: %s\n",
                ifelse(all_weighted_controlled, "PASS", "FAIL")))
    cat(sprintf("Standard phyper Type I error control: %s\n",
                ifelse(all_standard_controlled, "PASS", "FAIL")))
  }

  return(list(
    results = results,
    weighted_pvals = weighted_pvals,
    standard_pvals = standard_pvals,
    parameters = list(
      n_sims = n_sims,
      n_genes = n_genes,
      n_pathway = n_pathway,
      n_query = n_query,
      tolerance = tolerance
    ),
    pass_weighted = all(results$weighted_controlled),
    pass_standard = all(results$standard_controlled)
  ))
}

#' Validate Type I error across different copy number distributions
#'
#' Tests how different copy number distributions affect false positive rates.
#' Critical for understanding robustness of the weighted approach.
#'
#' @param copy_scenarios Named list of copy number generation functions
#' @param n_sims Number of simulations per scenario (default: 1000)
#' @param alpha Alpha level to test (default: 0.05)
#' @param n_genes Number of background genes (default: 500)
#' @param n_pathway Number of pathway genes (default: 50)
#' @param n_query Number of query genes (default: 30)
#' @param tolerance Acceptable deviation from nominal alpha (default: 0.01)
#' @param verbose Whether to print results (default: TRUE)
#' @return List with scenario results and overall assessment
validate_type_i_across_copy_distributions <- function(copy_scenarios = NULL,
                                                      n_sims = 1000,
                                                      alpha = 0.05,
                                                      n_genes = 500,
                                                      n_pathway = 50,
                                                      n_query = 30,
                                                      tolerance = 0.01,
                                                      verbose = TRUE) {

  # Default copy number scenarios if not provided
  if (is.null(copy_scenarios)) {
    copy_scenarios <- list(
      "All_CN_1" = function(n) rep(1, n),
      "All_CN_5" = function(n) rep(5, n),
      "Uniform_1_3" = function(n) sample(1:3, n, replace = TRUE),
      "Uniform_1_8" = function(n) sample(1:8, n, replace = TRUE),
      "Uniform_1_20" = function(n) sample(1:20, n, replace = TRUE),
      "High_Variance" = function(n) sample(c(1, 1, 1, 1, 1, 10, 20, 50), n, replace = TRUE),
      "Bimodal" = function(n) sample(c(rep(1, 5), rep(8, 3)), n, replace = TRUE)
    )
  }

  if (verbose) {
    cat("============================================\n")
    cat("TYPE I ERROR ACROSS COPY DISTRIBUTIONS\n")
    cat("============================================\n")
    cat(sprintf("Testing %d copy number scenarios\n", length(copy_scenarios)))
    cat(sprintf("Alpha level: %.3f (tolerance: ±%.1f%%)\n\n", alpha, tolerance * 100))
  }

  results <- data.frame(
    scenario = character(),
    weighted_rate = numeric(),
    standard_rate = numeric(),
    mean_cn = numeric(),
    var_cn = numeric(),
    max_cn = numeric(),
    weighted_controlled = logical(),
    standard_controlled = logical(),
    stringsAsFactors = FALSE
  )

  for (scenario_name in names(copy_scenarios)) {
    if (verbose) cat(sprintf("Testing scenario: %s...\n", scenario_name))

    set.seed(123 + which(names(copy_scenarios) == scenario_name))  # Reproducible per scenario

    # Generate background with this copy distribution
    copy_nums <- copy_scenarios[[scenario_name]](n_genes)
    background_df <- data.frame(
      gene = paste0("GENE", 1:n_genes),
      copy_number = copy_nums,
      stringsAsFactors = FALSE
    )

    pathway_genes <- sample(background_df$gene, n_pathway)

    # Run simulations
    weighted_pvals <- replicate(n_sims, {
      query_genes <- sample(background_df$gene, n_query)
      query_df <- background_df[background_df$gene %in% query_genes, ]
      weighted_hypergeometric_test_fixed(query_df, pathway_genes, background_df)$pvalue
    })

    standard_pvals <- replicate(n_sims, {
      query_genes <- sample(background_df$gene, n_query)
      overlap <- sum(query_genes %in% pathway_genes)
      phyper(overlap - 1, n_pathway, n_genes - n_pathway, n_query, lower.tail = FALSE)
    })

    # Calculate error rates
    weighted_rate <- mean(weighted_pvals < alpha)
    standard_rate <- mean(standard_pvals < alpha)

    # Check control
    weighted_controlled <- abs(weighted_rate - alpha) <= tolerance
    standard_controlled <- abs(standard_rate - alpha) <= tolerance

    # Store results
    results <- rbind(results, data.frame(
      scenario = scenario_name,
      weighted_rate = weighted_rate,
      standard_rate = standard_rate,
      mean_cn = mean(copy_nums),
      var_cn = var(copy_nums),
      max_cn = max(copy_nums),
      weighted_controlled = weighted_controlled,
      standard_controlled = standard_controlled,
      stringsAsFactors = FALSE
    ))
  }

  if (verbose) {
    cat("\nRESULTS BY COPY NUMBER SCENARIO:\n")
    cat(sprintf("%-15s %-10s %-10s %-8s %-8s %-8s %-11s %-11s\n",
                "Scenario", "Weighted", "Standard", "Mean_CN", "Var_CN", "Max_CN", "W_Control", "S_Control"))
    cat(sprintf("%-15s %-10s %-10s %-8s %-8s %-8s %-11s %-11s\n",
                "--------", "--------", "--------", "-------", "-------", "-------", "---------", "---------"))

    for (i in 1:nrow(results)) {
      cat(sprintf("%-15s %-10.4f %-10.4f %-8.1f %-8.1f %-8.0f %-11s %-11s\n",
                  results$scenario[i],
                  results$weighted_rate[i],
                  results$standard_rate[i],
                  results$mean_cn[i],
                  results$var_cn[i],
                  results$max_cn[i],
                  results$weighted_controlled[i],
                  results$standard_controlled[i]))
    }

    # Summary
    pass_weighted <- all(results$weighted_controlled)
    pass_standard <- all(results$standard_controlled)

    cat(sprintf("\nSUMMARY:\n"))
    cat(sprintf("Weighted method controls Type I error across all scenarios: %s\n",
                ifelse(pass_weighted, "PASS", "FAIL")))
    cat(sprintf("Standard method controls Type I error across all scenarios: %s\n",
                ifelse(pass_standard, "PASS", "FAIL")))

    # Identify problematic scenarios
    if (!pass_weighted) {
      problematic <- results[!results$weighted_controlled, "scenario"]
      cat(sprintf("Weighted method failed in: %s\n", paste(problematic, collapse = ", ")))
    }
    if (!pass_standard) {
      problematic <- results[!results$standard_controlled, "scenario"]
      cat(sprintf("Standard method failed in: %s\n", paste(problematic, collapse = ", ")))
    }
  }

  return(list(
    results = results,
    pass_weighted = all(results$weighted_controlled),
    pass_standard = all(results$standard_controlled),
    parameters = list(n_sims = n_sims, alpha = alpha, tolerance = tolerance)
  ))
}

#' Validate Type I error with varying background and pathway sizes
#'
#' Tests robustness across different dataset sizes - critical for ensuring
#' the method works across different experimental scales.
#'
#' @param size_scenarios List with background_size, pathway_size, query_size vectors
#' @param n_sims Number of simulations per scenario (default: 800)
#' @param alpha Alpha level to test (default: 0.05)
#' @param tolerance Acceptable deviation from alpha (default: 0.01)
#' @param verbose Whether to print results (default: TRUE)
#' @return List with size scenario results
validate_type_i_across_sizes <- function(size_scenarios = NULL,
                                         n_sims = 800,
                                         alpha = 0.05,
                                         tolerance = 0.01,
                                         verbose = TRUE) {

  # Default size scenarios
  if (is.null(size_scenarios)) {
    size_scenarios <- list(
      small = list(background = 200, pathway = 20, query = 15),
      medium = list(background = 1000, pathway = 100, query = 50),
      large = list(background = 5000, pathway = 500, query = 200),
      wide_pathway = list(background = 1000, pathway = 300, query = 50),
      narrow_pathway = list(background = 1000, pathway = 25, query = 50),
      large_query = list(background = 1000, pathway = 100, query = 150),
      small_query = list(background = 1000, pathway = 100, query = 20)
    )
  }

  if (verbose) {
    cat("============================================\n")
    cat("TYPE I ERROR ACROSS DATASET SIZES\n")
    cat("============================================\n")
    cat(sprintf("Testing %d size scenarios\n", length(size_scenarios)))
    cat(sprintf("Alpha: %.3f, Tolerance: ±%.1f%%\n\n", alpha, tolerance * 100))
  }

  results <- data.frame(
    scenario = character(),
    background_size = numeric(),
    pathway_size = numeric(),
    query_size = numeric(),
    weighted_rate = numeric(),
    standard_rate = numeric(),
    weighted_controlled = logical(),
    standard_controlled = logical(),
    stringsAsFactors = FALSE
  )

  for (scenario_name in names(size_scenarios)) {
    params <- size_scenarios[[scenario_name]]
    n_bg <- params$background
    n_pw <- params$pathway
    n_q <- params$query

    if (verbose) {
      cat(sprintf("Scenario %s: %d background, %d pathway, %d query...\n",
                  scenario_name, n_bg, n_pw, n_q))
    }

    set.seed(456 + which(names(size_scenarios) == scenario_name))

    # Generate background
    background_df <- data.frame(
      gene = paste0("GENE", 1:n_bg),
      copy_number = sample(1:8, n_bg, replace = TRUE),
      stringsAsFactors = FALSE
    )

    pathway_genes <- sample(background_df$gene, n_pw)

    # Run simulations
    weighted_pvals <- replicate(n_sims, {
      query_genes <- sample(background_df$gene, n_q)
      query_df <- background_df[background_df$gene %in% query_genes, ]
      weighted_hypergeometric_test_fixed(query_df, pathway_genes, background_df)$pvalue
    })

    standard_pvals <- replicate(n_sims, {
      query_genes <- sample(background_df$gene, n_q)
      overlap <- sum(query_genes %in% pathway_genes)
      phyper(overlap - 1, n_pw, n_bg - n_pw, n_q, lower.tail = FALSE)
    })

    weighted_rate <- mean(weighted_pvals < alpha)
    standard_rate <- mean(standard_pvals < alpha)

    results <- rbind(results, data.frame(
      scenario = scenario_name,
      background_size = n_bg,
      pathway_size = n_pw,
      query_size = n_q,
      weighted_rate = weighted_rate,
      standard_rate = standard_rate,
      weighted_controlled = abs(weighted_rate - alpha) <= tolerance,
      standard_controlled = abs(standard_rate - alpha) <= tolerance,
      stringsAsFactors = FALSE
    ))
  }

  if (verbose) {
    cat("\nRESULTS BY SIZE SCENARIO:\n")
    cat(sprintf("%-12s %-6s %-7s %-6s %-10s %-10s %-11s %-11s\n",
                "Scenario", "Bg", "Pathway", "Query", "Weighted", "Standard", "W_Control", "S_Control"))
    cat(sprintf("%-12s %-6s %-7s %-6s %-10s %-10s %-11s %-11s\n",
                "--------", "--", "-------", "-----", "--------", "--------", "---------", "---------"))

    for (i in 1:nrow(results)) {
      cat(sprintf("%-12s %-6.0f %-7.0f %-6.0f %-10.4f %-10.4f %-11s %-11s\n",
                  results$scenario[i],
                  results$background_size[i],
                  results$pathway_size[i],
                  results$query_size[i],
                  results$weighted_rate[i],
                  results$standard_rate[i],
                  results$weighted_controlled[i],
                  results$standard_controlled[i]))
    }

    pass_weighted <- all(results$weighted_controlled)
    pass_standard <- all(results$standard_controlled)

    cat(sprintf("\nSUMMARY:\n"))
    cat(sprintf("Weighted method Type I control across sizes: %s\n",
                ifelse(pass_weighted, "PASS", "FAIL")))
    cat(sprintf("Standard method Type I control across sizes: %s\n",
                ifelse(pass_standard, "PASS", "FAIL")))
  }

  return(list(
    results = results,
    pass_weighted = all(results$weighted_controlled),
    pass_standard = all(results$standard_controlled),
    parameters = list(n_sims = n_sims, alpha = alpha, tolerance = tolerance)
  ))
}

# ============================================================
# ERROR RATE MEASUREMENT AND REPORTING FUNCTIONS
# ============================================================

#' Comprehensive Type I error measurement across all scenarios
#'
#' Runs all Type I error validation tests and generates comprehensive report
#'
#' @param quick_run If TRUE, uses reduced simulation counts (default: FALSE)
#' @param output_file File to save detailed results (default: NULL)
#' @return List with all validation results
comprehensive_type_i_validation <- function(quick_run = FALSE, output_file = NULL) {

  n_sims <- ifelse(quick_run, 500, 2000)

  cat("####################################################\n")
  cat("COMPREHENSIVE TYPE I ERROR VALIDATION\n")
  cat("####################################################\n")
  cat(sprintf("Mode: %s\n", ifelse(quick_run, "Quick", "Full")))
  cat(sprintf("Simulations per test: %d\n\n", n_sims))

  # Test 1: Multiple alpha levels
  cat("TEST 1: Multiple Alpha Levels\n")
  cat("------------------------------\n")
  result1 <- validate_type_i_error_rates(n_sims = n_sims)

  cat("\n\n")

  # Test 2: Copy number distributions
  cat("TEST 2: Copy Number Distributions\n")
  cat("----------------------------------\n")
  result2 <- validate_type_i_across_copy_distributions(n_sims = ifelse(quick_run, 400, 1000))

  cat("\n\n")

  # Test 3: Dataset sizes
  cat("TEST 3: Dataset Sizes\n")
  cat("---------------------\n")
  result3 <- validate_type_i_across_sizes(n_sims = ifelse(quick_run, 300, 800))

  cat("\n\n")

  # Overall summary
  cat("####################################################\n")
  cat("OVERALL TYPE I ERROR VALIDATION SUMMARY\n")
  cat("####################################################\n")

  all_tests_pass_weighted <- result1$pass_weighted && result2$pass_weighted && result3$pass_weighted
  all_tests_pass_standard <- result1$pass_standard && result2$pass_standard && result3$pass_standard

  cat(sprintf("✓ Multiple alpha levels - Weighted: %s, Standard: %s\n",
              ifelse(result1$pass_weighted, "PASS", "FAIL"),
              ifelse(result1$pass_standard, "PASS", "FAIL")))
  cat(sprintf("✓ Copy distributions    - Weighted: %s, Standard: %s\n",
              ifelse(result2$pass_weighted, "PASS", "FAIL"),
              ifelse(result2$pass_standard, "PASS", "FAIL")))
  cat(sprintf("✓ Dataset sizes         - Weighted: %s, Standard: %s\n",
              ifelse(result3$pass_weighted, "PASS", "FAIL"),
              ifelse(result3$pass_standard, "PASS", "FAIL")))

  cat("\nFINAL VERDICT:\n")
  cat(sprintf("Weighted phyper() Type I error control: %s\n",
              ifelse(all_tests_pass_weighted, "✓ PASS", "✗ FAIL")))
  cat(sprintf("Standard phyper() Type I error control: %s\n",
              ifelse(all_tests_pass_standard, "✓ PASS", "✗ FAIL")))

  # Compile results
  comprehensive_results <- list(
    alpha_levels = result1,
    copy_distributions = result2,
    dataset_sizes = result3,
    overall_pass_weighted = all_tests_pass_weighted,
    overall_pass_standard = all_tests_pass_standard,
    timestamp = Sys.time(),
    session_info = sessionInfo()
  )

  # Save to file if requested
  if (!is.null(output_file)) {
    cat(sprintf("\nSaving detailed results to: %s\n", output_file))
    saveRDS(comprehensive_results, output_file)
  }

  return(comprehensive_results)
}

#' Generate Type I error validation report
#'
#' Creates a formatted report of Type I error validation results
#'
#' @param results Results from comprehensive_type_i_validation()
#' @param output_file Markdown file to save report (default: "type_i_error_report.md")
generate_type_i_error_report <- function(results, output_file = "type_i_error_report.md") {

  report_lines <- c(
    "# Type I Error Rate Validation Report",
    "",
    sprintf("Generated: %s", as.character(results$timestamp)),
    "",
    "## Executive Summary",
    "",
    sprintf("**Overall Weighted phyper() Type I Error Control: %s**",
            ifelse(results$overall_pass_weighted, "✓ PASS", "✗ FAIL")),
    sprintf("**Overall Standard phyper() Type I Error Control: %s**",
            ifelse(results$overall_pass_standard, "✓ PASS", "✗ FAIL")),
    "",
    "## Test Results",
    "",
    "### 1. Multiple Alpha Levels",
    ""
  )

  # Alpha levels table
  alpha_results <- results$alpha_levels$results
  report_lines <- c(report_lines,
    "| Alpha | Weighted Rate | Standard Rate | W_Controlled | S_Controlled |",
    "|-------|---------------|---------------|--------------|--------------|"
  )

  for (i in 1:nrow(alpha_results)) {
    report_lines <- c(report_lines,
      sprintf("| %.3f | %.4f | %.4f | %s | %s |",
              alpha_results$alpha[i],
              alpha_results$weighted_rate[i],
              alpha_results$standard_rate[i],
              alpha_results$weighted_controlled[i],
              alpha_results$standard_controlled[i]))
  }

  # Copy distributions
  copy_results <- results$copy_distributions$results
  report_lines <- c(report_lines,
    "",
    "### 2. Copy Number Distributions",
    "",
    "| Scenario | Weighted Rate | Standard Rate | Mean CN | W_Controlled | S_Controlled |",
    "|----------|---------------|---------------|---------|--------------|--------------|"
  )

  for (i in 1:nrow(copy_results)) {
    report_lines <- c(report_lines,
      sprintf("| %s | %.4f | %.4f | %.1f | %s | %s |",
              copy_results$scenario[i],
              copy_results$weighted_rate[i],
              copy_results$standard_rate[i],
              copy_results$mean_cn[i],
              copy_results$weighted_controlled[i],
              copy_results$standard_controlled[i]))
  }

  # Size scenarios
  size_results <- results$dataset_sizes$results
  report_lines <- c(report_lines,
    "",
    "### 3. Dataset Sizes",
    "",
    "| Scenario | Background | Pathway | Query | Weighted Rate | Standard Rate | W_Controlled | S_Controlled |",
    "|----------|------------|---------|-------|---------------|---------------|--------------|--------------|"
  )

  for (i in 1:nrow(size_results)) {
    report_lines <- c(report_lines,
      sprintf("| %s | %.0f | %.0f | %.0f | %.4f | %.4f | %s | %s |",
              size_results$scenario[i],
              size_results$background_size[i],
              size_results$pathway_size[i],
              size_results$query_size[i],
              size_results$weighted_rate[i],
              size_results$standard_rate[i],
              size_results$weighted_controlled[i],
              size_results$standard_controlled[i]))
  }

  # Add conclusions
  report_lines <- c(report_lines,
    "",
    "## Conclusions",
    "",
    sprintf("- Type I error rates were tested at tolerance level of ±1%%"),
    sprintf("- Weighted method overall pass rate: %s", ifelse(results$overall_pass_weighted, "100%", "Failed")),
    sprintf("- Standard method overall pass rate: %s", ifelse(results$overall_pass_standard, "100%", "Failed"))
  )

  # Write report
  writeLines(report_lines, output_file)
  cat(sprintf("Type I error validation report saved to: %s\n", output_file))

  return(output_file)
}

# ============================================================
# QUICK VALIDATION FUNCTION
# ============================================================

#' Quick Type I error validation check
#'
#' Runs a streamlined validation for fast checks during development
#'
#' @param alpha Alpha level to test (default: 0.05)
#' @param n_sims Number of simulations (default: 500)
#' @return Boolean indicating if validation passes
quick_type_i_check <- function(alpha = 0.05, n_sims = 500) {
  cat("Quick Type I Error Validation Check\n")
  cat("===================================\n")

  result <- validate_type_i_error_rates(
    n_sims = n_sims,
    alphas = alpha,
    n_genes = 500,
    n_pathway = 50,
    n_query = 25,
    tolerance = 0.01,
    verbose = FALSE
  )

  pass_result <- result$pass_weighted

  cat(sprintf("Result: %s\n", ifelse(pass_result, "✓ PASS", "✗ FAIL")))
  cat(sprintf("Observed rate: %.4f (expected: %.3f)\n",
              result$results$weighted_rate, alpha))

  return(pass_result)
}