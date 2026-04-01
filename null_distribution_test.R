# Null Distribution Validation Functions for Copy-Number Weighted phyper()
#
# This module implements comprehensive validation that p-values from copy-number
# weighted phyper() follow Uniform(0,1) under the null hypothesis.
#
# Author: AI Assistant
# Date: 2026-04-01
# Task: null-distribution-validation

library(tidyverse)

# Try to load nortest, but continue without it if unavailable
nortest_available <- requireNamespace("nortest", quietly = TRUE)
if (!nortest_available) {
  warning("nortest package not available - Anderson-Darling test will be skipped")
}

# Source the weighted phyper implementation
source("debug_weighted_phyper.R")

# ==============================================================================
# CORE NULL DISTRIBUTION SIMULATION FUNCTIONS
# ==============================================================================

#' Generate background datasets with different copy number distributions
#'
#' @param n_genes Number of genes
#' @param distribution Type of copy number distribution: "uniform", "skewed", "realistic"
#' @param seed Random seed for reproducibility
#' @return Data frame with gene names and copy numbers
generate_background_with_copy_distribution <- function(n_genes, distribution = "uniform", seed = 42) {
  set.seed(seed)

  gene_names <- paste0("GENE", sprintf("%05d", 1:n_genes))

  copy_numbers <- switch(distribution,
    "uniform" = sample(1:8, n_genes, replace = TRUE),
    "skewed" = sample(1:10, n_genes, replace = TRUE, prob = c(0.4, 0.25, 0.15, 0.1, 0.05, 0.03, 0.015, 0.005, 0.003, 0.002)),
    "realistic" = {
      # Based on real copy number variation patterns
      # Most genes have 2 copies, some have CNVs
      base_copies <- rep(2, n_genes)
      # 10% have copy gains (3-8 copies)
      gain_indices <- sample(n_genes, round(0.1 * n_genes))
      base_copies[gain_indices] <- sample(3:8, length(gain_indices), replace = TRUE, prob = c(0.5, 0.25, 0.125, 0.0625, 0.0625, 0.0))
      # 5% have copy losses (1 copy)
      loss_indices <- sample(setdiff(1:n_genes, gain_indices), round(0.05 * n_genes))
      base_copies[loss_indices] <- 1
      base_copies
    },
    stop("Unknown distribution type: ", distribution)
  )

  data.frame(
    gene = gene_names,
    copy_number = copy_numbers,
    stringsAsFactors = FALSE
  )
}

#' Simulate null hypothesis p-values with various copy number scenarios
#'
#' @param n_simulations Number of simulation replicates
#' @param background_size Number of genes in background
#' @param pathway_size Number of genes in pathway
#' @param query_size Number of genes in query set
#' @param copy_distribution Copy number distribution type
#' @param seed Random seed for reproducibility
#' @return Vector of p-values under null hypothesis
simulate_null_pvalues <- function(n_simulations = 1000,
                                  background_size = 500,
                                  pathway_size = 50,
                                  query_size = 30,
                                  copy_distribution = "uniform",
                                  seed = 42) {

  set.seed(seed)

  # Generate stable background dataset
  background_df <- generate_background_with_copy_distribution(
    n_genes = background_size,
    distribution = copy_distribution,
    seed = seed
  )

  # Define pathway genes (fixed for all simulations)
  pathway_genes <- sample(background_df$gene, pathway_size)

  # Simulate null queries by random sampling
  pvalues <- replicate(n_simulations, {
    query_genes <- sample(background_df$gene, query_size)
    query_df <- background_df[background_df$gene %in% query_genes, ]

    result <- weighted_hypergeometric_test_fixed(query_df, pathway_genes, background_df)
    result$pvalue
  })

  return(pvalues)
}

# ==============================================================================
# P-VALUE UNIFORMITY TESTS
# ==============================================================================

#' Test p-value uniformity using multiple statistical tests
#'
#' @param pvalues Vector of p-values to test
#' @param alpha Significance level for tests
#' @return List with test results and uniformity assessment
test_pvalue_uniformity <- function(pvalues, alpha = 0.05) {

  # Remove any NAs or invalid p-values
  valid_pvals <- pvalues[!is.na(pvalues) & pvalues >= 0 & pvalues <= 1]
  n_valid <- length(valid_pvals)
  n_removed <- length(pvalues) - n_valid

  if (n_valid < 10) {
    stop("Insufficient valid p-values for uniformity testing")
  }

  # Kolmogorov-Smirnov test
  ks_result <- ks.test(valid_pvals, punif)

  # Anderson-Darling test (more sensitive to tail deviations)
  ad_result <- if (nortest_available) {
    tryCatch({
      nortest::ad.test(valid_pvals)
    }, error = function(e) {
      list(statistic = NA, p.value = NA, method = "Anderson-Darling failed")
    })
  } else {
    list(statistic = NA, p.value = NA, method = "Anderson-Darling not available")
  }

  # Type I error rate (should be approximately alpha)
  type_i_rate <- mean(valid_pvals < alpha)
  type_i_expected <- alpha
  type_i_se <- sqrt(alpha * (1 - alpha) / n_valid)
  type_i_z <- (type_i_rate - type_i_expected) / type_i_se
  type_i_pvalue <- 2 * (1 - pnorm(abs(type_i_z)))

  # Summary statistics
  mean_pval <- mean(valid_pvals)
  median_pval <- median(valid_pvals)

  # Overall uniformity assessment
  uniformity_ks <- ks_result$p.value > 0.05
  uniformity_ad <- if (!is.na(ad_result$p.value)) ad_result$p.value > 0.05 else FALSE
  uniformity_type_i <- type_i_pvalue > 0.05

  overall_uniform <- uniformity_ks && (uniformity_ad || is.na(ad_result$p.value)) && uniformity_type_i

  return(list(
    n_pvalues = n_valid,
    n_removed = n_removed,
    ks_test = list(
      statistic = ks_result$statistic,
      p_value = ks_result$p.value,
      uniform = uniformity_ks
    ),
    ad_test = list(
      statistic = ad_result$statistic,
      p_value = ad_result$p.value,
      uniform = uniformity_ad
    ),
    type_i_error = list(
      observed = type_i_rate,
      expected = type_i_expected,
      z_score = type_i_z,
      p_value = type_i_pvalue,
      correct = uniformity_type_i
    ),
    summary_stats = list(
      mean = mean_pval,
      median = median_pval,
      min = min(valid_pvals),
      max = max(valid_pvals)
    ),
    overall_uniform = overall_uniform
  ))
}

# ==============================================================================
# COMPREHENSIVE VALIDATION SCENARIOS
# ==============================================================================

#' Run comprehensive null distribution validation across multiple scenarios
#'
#' @param scenarios List of scenario configurations
#' @param n_simulations Number of simulations per scenario
#' @return List of validation results for each scenario
validate_null_distribution_comprehensive <- function(scenarios = NULL, n_simulations = 1000) {

  # Default scenarios if not provided
  if (is.null(scenarios)) {
    scenarios <- list(
      # Different background sizes
      list(name = "small_bg", background_size = 200, pathway_size = 20, query_size = 15, copy_distribution = "uniform"),
      list(name = "medium_bg", background_size = 500, pathway_size = 50, query_size = 30, copy_distribution = "uniform"),
      list(name = "large_bg", background_size = 1000, pathway_size = 100, query_size = 50, copy_distribution = "uniform"),

      # Different pathway sizes
      list(name = "small_pathway", background_size = 500, pathway_size = 25, query_size = 30, copy_distribution = "uniform"),
      list(name = "large_pathway", background_size = 500, pathway_size = 100, query_size = 30, copy_distribution = "uniform"),

      # Different copy number distributions
      list(name = "uniform_copies", background_size = 500, pathway_size = 50, query_size = 30, copy_distribution = "uniform"),
      list(name = "skewed_copies", background_size = 500, pathway_size = 50, query_size = 30, copy_distribution = "skewed"),
      list(name = "realistic_copies", background_size = 500, pathway_size = 50, query_size = 30, copy_distribution = "realistic")
    )
  }

  cat("Running comprehensive null distribution validation\n")
  cat("==============================================\n\n")

  results <- list()

  for (i in seq_along(scenarios)) {
    scenario <- scenarios[[i]]
    cat("Scenario", i, ":", scenario$name, "\n")
    cat("  Background:", scenario$background_size, "genes\n")
    cat("  Pathway:", scenario$pathway_size, "genes\n")
    cat("  Query:", scenario$query_size, "genes\n")
    cat("  Copy distribution:", scenario$copy_distribution, "\n")

    # Simulate p-values under null
    pvalues <- simulate_null_pvalues(
      n_simulations = n_simulations,
      background_size = scenario$background_size,
      pathway_size = scenario$pathway_size,
      query_size = scenario$query_size,
      copy_distribution = scenario$copy_distribution,
      seed = 42 + i  # Different seed per scenario
    )

    # Test uniformity
    uniformity_result <- test_pvalue_uniformity(pvalues)

    # Store results
    results[[scenario$name]] <- list(
      scenario = scenario,
      pvalues = pvalues,
      uniformity = uniformity_result
    )

    # Print summary
    cat("  Results:\n")
    cat("    KS p-value:", format(uniformity_result$ks_test$p_value, digits = 4),
        " (uniform:", uniformity_result$ks_test$uniform, ")\n")
    cat("    AD p-value:", format(uniformity_result$ad_test$p_value, digits = 4),
        " (uniform:", uniformity_result$ad_test$uniform, ")\n")
    cat("    Type I rate:", round(uniformity_result$type_i_error$observed, 3),
        " (correct:", uniformity_result$type_i_error$correct, ")\n")
    cat("    Overall uniform:", uniformity_result$overall_uniform, "\n\n")
  }

  # Overall summary
  uniform_scenarios <- sum(sapply(results, function(x) x$uniformity$overall_uniform))
  total_scenarios <- length(results)

  cat("OVERALL VALIDATION SUMMARY\n")
  cat("===========================\n")
  cat("Scenarios passing uniformity tests:", uniform_scenarios, "/", total_scenarios, "\n")
  cat("Overall validation success:", uniform_scenarios == total_scenarios, "\n")

  return(results)
}

# ==============================================================================
# DIAGNOSTIC AND VISUALIZATION FUNCTIONS
# ==============================================================================

#' Generate diagnostic plots for p-value distribution
#'
#' @param pvalues Vector of p-values
#' @param title Plot title
#' @return ggplot object with diagnostic plots
plot_pvalue_diagnostics <- function(pvalues, title = "P-value Distribution Diagnostics") {

  # Create data frame for plotting
  plot_data <- data.frame(pvalue = pvalues)

  # Histogram
  p1 <- ggplot(plot_data, aes(x = pvalue)) +
    geom_histogram(bins = 20, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_hline(yintercept = length(pvalues) / 20, linetype = "dashed", color = "red") +
    labs(x = "P-value", y = "Frequency", title = paste(title, "- Histogram")) +
    theme_minimal()

  # Q-Q plot against uniform distribution
  p2 <- ggplot(plot_data, aes(sample = pvalue)) +
    stat_qq(distribution = qunif) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(x = "Theoretical Uniform Quantiles", y = "Observed P-value Quantiles",
         title = paste(title, "- Q-Q Plot vs Uniform")) +
    theme_minimal()

  # Cumulative distribution
  p3 <- ggplot(plot_data, aes(x = pvalue)) +
    stat_ecdf() +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(x = "P-value", y = "Cumulative Probability",
         title = paste(title, "- Empirical CDF vs Uniform")) +
    theme_minimal()

  # Combine plots (requires gridExtra or similar)
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
  } else {
    # Return first plot if gridExtra not available
    p1
  }
}

#' Generate validation report
#'
#' @param validation_results Results from validate_null_distribution_comprehensive
#' @return Character vector with formatted report
generate_validation_report <- function(validation_results) {

  report <- c()
  report <- c(report, "COPY-NUMBER WEIGHTED PHYPER() NULL DISTRIBUTION VALIDATION REPORT")
  report <- c(report, paste(rep("=", 70), collapse = ""))
  report <- c(report, paste("Generated:", Sys.time()))
  report <- c(report, "")

  # Summary statistics
  n_scenarios <- length(validation_results)
  n_uniform <- sum(sapply(validation_results, function(x) x$uniformity$overall_uniform))

  report <- c(report, "EXECUTIVE SUMMARY")
  report <- c(report, paste(rep("-", 20), collapse = ""))
  report <- c(report, paste("Total scenarios tested:", n_scenarios))
  report <- c(report, paste("Scenarios passing uniformity tests:", n_uniform))
  report <- c(report, paste("Overall validation success:", n_uniform == n_scenarios))
  report <- c(report, "")

  # Detailed results for each scenario
  for (scenario_name in names(validation_results)) {
    result <- validation_results[[scenario_name]]
    scenario <- result$scenario
    uniformity <- result$uniformity

    report <- c(report, paste("SCENARIO:", toupper(scenario_name)))
    report <- c(report, paste(rep("-", 30), collapse = ""))
    report <- c(report, paste("Background genes:", scenario$background_size))
    report <- c(report, paste("Pathway genes:", scenario$pathway_size))
    report <- c(report, paste("Query genes:", scenario$query_size))
    report <- c(report, paste("Copy distribution:", scenario$copy_distribution))
    report <- c(report, paste("P-values tested:", uniformity$n_pvalues))
    report <- c(report, "")

    report <- c(report, "Statistical Tests:")
    report <- c(report, paste("  Kolmogorov-Smirnov test p-value:", format(uniformity$ks_test$p_value, digits = 4)))
    report <- c(report, paste("  Anderson-Darling test p-value:", format(uniformity$ad_test$p_value, digits = 4)))
    report <- c(report, paste("  Type I error rate:", format(uniformity$type_i_error$observed, digits = 3),
                             " (expected: 0.050)"))
    report <- c(report, "")

    report <- c(report, "Uniformity Assessment:")
    report <- c(report, paste("  KS test passed:", uniformity$ks_test$uniform))
    report <- c(report, paste("  AD test passed:", uniformity$ad_test$uniform))
    report <- c(report, paste("  Type I rate correct:", uniformity$type_i_error$correct))
    report <- c(report, paste("  Overall uniform:", uniformity$overall_uniform))
    report <- c(report, "")
  }

  # Conclusions
  report <- c(report, "CONCLUSIONS")
  report <- c(report, paste(rep("-", 15), collapse = ""))

  if (n_uniform == n_scenarios) {
    report <- c(report, "✓ VALIDATION SUCCESSFUL: Copy-number weighted phyper() produces")
    report <- c(report, "  uniform p-values under the null hypothesis across all tested scenarios.")
    report <- c(report, "✓ The implementation correctly maintains Type I error control.")
    report <- c(report, "✓ Statistical behavior is consistent across different background sizes,")
    report <- c(report, "  pathway sizes, and copy number distributions.")
  } else {
    failed_scenarios <- names(validation_results)[!sapply(validation_results, function(x) x$uniformity$overall_uniform)]
    report <- c(report, "✗ VALIDATION FAILED: Non-uniform p-values detected in scenarios:")
    for (failed in failed_scenarios) {
      report <- c(report, paste("  -", failed))
    }
    report <- c(report, "✗ Implementation requires further debugging and correction.")
  }

  return(report)
}

# ==============================================================================
# MAIN VALIDATION EXECUTION
# ==============================================================================

#' Main function to run complete null distribution validation
#'
#' @param n_simulations Number of simulations per scenario
#' @param output_file Optional file to save results
#' @return Validation results and printed report
run_complete_null_validation <- function(n_simulations = 1000, output_file = NULL) {

  cat("STARTING COMPREHENSIVE NULL DISTRIBUTION VALIDATION\n")
  cat("===================================================\n\n")

  # Run validation across all scenarios
  validation_results <- validate_null_distribution_comprehensive(
    scenarios = NULL,  # Use default scenarios
    n_simulations = n_simulations
  )

  # Generate and print report
  report_lines <- generate_validation_report(validation_results)

  cat("\n")
  for (line in report_lines) {
    cat(line, "\n")
  }

  # Save to file if requested
  if (!is.null(output_file)) {
    writeLines(report_lines, output_file)
    cat("\nReport saved to:", output_file, "\n")
  }

  return(validation_results)
}

# ==============================================================================
# SCRIPT EXECUTION (when not sourced interactively)
# ==============================================================================

if (!interactive()) {
  cat("Executing null distribution validation...\n\n")

  # Run with default parameters
  results <- run_complete_null_validation(
    n_simulations = 500,  # Reduced for faster execution
    output_file = "null_distribution_validation_report.txt"
  )

  cat("\nNull distribution validation complete.\n")
}