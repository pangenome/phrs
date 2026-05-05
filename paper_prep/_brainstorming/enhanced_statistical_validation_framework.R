# Enhanced Statistical Validation Framework for Copy-Number Weighted Methods
#
# This framework addresses the critical statistical issues identified in initial validation:
# 1. Independence assumption violations in gene-level sampling
# 2. Anti-conservative p-values due to cluster sampling effects
# 3. Type I error inflation with copy number weighting
#
# Author: AI Assistant
# Task: statistical-validation-framework
# Date: 2026-04-01

library(tidyverse)
library(parallel)

# ==============================================================================
# CORRECTED IMPLEMENTATIONS
# ==============================================================================

#' Permutation-based copy-number weighted test
#'
#' Uses gene-level permutation to generate appropriate null distribution
#' that accounts for clustering effects of copy numbers.
#'
#' @param query_df Data frame with gene and copy_number columns
#' @param pathway_genes Vector of pathway gene names
#' @param background_df Data frame with all background genes and copy numbers
#' @param n_permutations Number of permutations for null distribution
#' @return List with p-value and test statistics
permutation_weighted_test <- function(query_df, pathway_genes, background_df,
                                     n_permutations = 10000) {

  # Filter to valid genes and use consistent copy numbers
  valid_query_genes <- intersect(query_df$gene, background_df$gene)
  query_size <- length(valid_query_genes)

  if (query_size == 0) {
    stop("No query genes found in background dataset")
  }

  # Calculate observed test statistic (weighted overlap)
  query_filtered <- query_df[query_df$gene %in% valid_query_genes, ]
  query_with_bg_copies <- merge(query_filtered[, "gene", drop = FALSE],
                               background_df, by = "gene", all.x = TRUE)

  observed_overlap <- sum(query_with_bg_copies$copy_number[
    query_with_bg_copies$gene %in% pathway_genes])

  # Generate null distribution via gene-level permutation
  null_overlaps <- replicate(n_permutations, {
    # Sample genes at gene level (not instance level)
    permuted_genes <- sample(background_df$gene, query_size)
    permuted_df <- background_df[background_df$gene %in% permuted_genes, ]
    sum(permuted_df$copy_number[permuted_df$gene %in% pathway_genes])
  })

  # Calculate empirical p-value
  pvalue <- (sum(null_overlaps >= observed_overlap) + 1) / (n_permutations + 1)

  # Additional statistics
  null_mean <- mean(null_overlaps)
  null_sd <- sd(null_overlaps)
  z_score <- (observed_overlap - null_mean) / null_sd

  return(list(
    pvalue = pvalue,
    observed_overlap = observed_overlap,
    expected_overlap = null_mean,
    fold_enrichment = observed_overlap / null_mean,
    z_score = z_score,
    null_distribution = null_overlaps,
    method = "permutation"
  ))
}

#' Effective sample size corrected weighted test
#'
#' Adjusts the effective sample size to account for clustering effects
#' while maintaining the efficiency of parametric testing.
#'
#' @param query_df Data frame with gene and copy_number columns
#' @param pathway_genes Vector of pathway gene names
#' @param background_df Data frame with all background genes and copy numbers
#' @param correction_method Method for calculating effective sample size
#' @return List with corrected p-value and statistics
effective_sample_corrected_test <- function(query_df, pathway_genes, background_df,
                                           correction_method = "mean_copies") {

  # Calculate standard weighted parameters
  valid_query_genes <- intersect(query_df$gene, background_df$gene)
  query_filtered <- query_df[query_df$gene %in% valid_query_genes, ]
  query_with_bg_copies <- merge(query_filtered[, "gene", drop = FALSE],
                               background_df, by = "gene", all.x = TRUE)

  query_in_pathway <- query_with_bg_copies[query_with_bg_copies$gene %in% pathway_genes, ]
  pathway_in_background <- background_df[background_df$gene %in% pathway_genes, ]

  k_weighted <- sum(query_with_bg_copies$copy_number)
  q_weighted <- sum(query_in_pathway$copy_number)
  m_weighted <- sum(pathway_in_background$copy_number)
  n_weighted <- sum(background_df$copy_number) - m_weighted

  # Calculate effective sample size correction
  if (correction_method == "mean_copies") {
    mean_copy_number <- k_weighted / length(valid_query_genes)
    k_effective <- k_weighted / mean_copy_number
    q_effective <- q_weighted / mean_copy_number
  } else if (correction_method == "variance_adjustment") {
    copy_variance <- var(query_with_bg_copies$copy_number)
    copy_mean <- mean(query_with_bg_copies$copy_number)
    # Design effect for cluster sampling
    design_effect <- 1 + (copy_variance / copy_mean - 1)
    k_effective <- k_weighted / design_effect
    q_effective <- q_weighted / design_effect
  } else {
    stop("Unknown correction method")
  }

  # Apply corrected hypergeometric test
  pvalue_corrected <- phyper(q_effective - 1, m_weighted, n_weighted, k_effective,
                            lower.tail = FALSE)

  # Also calculate uncorrected for comparison
  pvalue_uncorrected <- phyper(q_weighted - 1, m_weighted, n_weighted, k_weighted,
                              lower.tail = FALSE)

  return(list(
    pvalue = pvalue_corrected,
    pvalue_uncorrected = pvalue_uncorrected,
    correction_factor = k_weighted / k_effective,
    observed_overlap = q_weighted,
    effective_overlap = q_effective,
    effective_sample_size = k_effective,
    method = paste0("effective_sample_", correction_method)
  ))
}

# ==============================================================================
# ENHANCED VALIDATION FRAMEWORK
# ==============================================================================

#' Validate null distribution for permutation-based test
#'
#' Tests whether permutation-based approach produces uniform p-values
#' under gene-level null hypothesis sampling.
#'
#' @param n_simulations Number of simulation runs
#' @param n_permutations Number of permutations per test
#' @param background_size Number of background genes
#' @param pathway_size Number of pathway genes
#' @param query_size Number of query genes
#' @param copy_distribution Type of copy number distribution
validate_permutation_null_distribution <- function(n_simulations = 500,
                                                  n_permutations = 1000,
                                                  background_size = 1000,
                                                  pathway_size = 100,
                                                  query_size = 50,
                                                  copy_distribution = "uniform") {

  cat("Validating permutation-based null distribution...\n")
  cat("Parameters: bg =", background_size, "pw =", pathway_size,
      "query =", query_size, "dist =", copy_distribution, "\n")

  # Generate background dataset
  if (copy_distribution == "uniform") {
    copy_numbers <- sample(1:8, background_size, replace = TRUE)
  } else if (copy_distribution == "skewed") {
    copy_numbers <- pmax(1, rpois(background_size, 3))
  } else if (copy_distribution == "realistic") {
    # Based on real genomic copy number distributions
    copy_numbers <- sample(c(rep(1:2, each = background_size * 0.4),
                            rep(3:5, each = background_size * 0.15),
                            rep(6:20, each = background_size * 0.05)),
                          background_size, replace = TRUE)
  }

  background_df <- data.frame(
    gene = paste0("G", 1:background_size),
    copy_number = copy_numbers,
    stringsAsFactors = FALSE
  )

  pathway_genes <- sample(background_df$gene, pathway_size)

  # Run null simulations with gene-level sampling
  pvalues <- replicate(n_simulations, {
    # Sample genes uniformly (gene-level null)
    query_genes <- sample(background_df$gene, query_size)
    query_df <- background_df[background_df$gene %in% query_genes, ]

    # Run permutation test
    result <- permutation_weighted_test(query_df, pathway_genes, background_df,
                                       n_permutations = n_permutations)
    result$pvalue
  })

  # Test uniformity
  ks_result <- ks.test(pvalues, punif)
  type_i_rate <- mean(pvalues < 0.05)

  cat("Results:\n")
  cat("  KS test p-value:", format(ks_result$p.value, digits = 4), "\n")
  cat("  Type I error rate:", round(type_i_rate, 3), "(expected: 0.05)\n")
  cat("  Uniformity test:", ks_result$p.value > 0.05, "\n")
  cat("  Type I control:", abs(type_i_rate - 0.05) < 0.02, "\n")

  return(list(
    pvalues = pvalues,
    ks_pvalue = ks_result$p.value,
    type_i_rate = type_i_rate,
    uniform = ks_result$p.value > 0.05,
    type_i_controlled = abs(type_i_rate - 0.05) < 0.02,
    scenario = paste(copy_distribution, background_size, pathway_size, query_size, sep = "_")
  ))
}

#' Validate effective sample size correction
#'
#' Tests whether effective sample size correction restores uniform p-values.
#'
#' @param correction_method Type of correction to validate
#' @param ... Other parameters passed to validation function
validate_effective_sample_correction <- function(correction_method = "mean_copies",
                                                n_simulations = 500,
                                                background_size = 1000,
                                                pathway_size = 100,
                                                query_size = 50,
                                                copy_distribution = "uniform") {

  cat("Validating effective sample size correction:", correction_method, "\n")
  cat("Parameters: bg =", background_size, "pw =", pathway_size,
      "query =", query_size, "dist =", copy_distribution, "\n")

  # Generate background dataset
  if (copy_distribution == "uniform") {
    copy_numbers <- sample(1:8, background_size, replace = TRUE)
  } else if (copy_distribution == "skewed") {
    copy_numbers <- pmax(1, rpois(background_size, 3))
  } else {
    copy_numbers <- sample(c(rep(1:2, each = background_size * 0.4),
                            rep(3:8, each = background_size * 0.1)),
                          background_size, replace = TRUE)
  }

  background_df <- data.frame(
    gene = paste0("G", 1:background_size),
    copy_number = copy_numbers,
    stringsAsFactors = FALSE
  )

  pathway_genes <- sample(background_df$gene, pathway_size)

  # Run null simulations
  results <- replicate(n_simulations, {
    query_genes <- sample(background_df$gene, query_size)
    query_df <- background_df[background_df$gene %in% query_genes, ]

    result <- effective_sample_corrected_test(query_df, pathway_genes, background_df,
                                             correction_method = correction_method)
    c(corrected = result$pvalue, uncorrected = result$pvalue_uncorrected,
      correction_factor = result$correction_factor)
  }, simplify = TRUE)

  corrected_pvals <- results["corrected", ]
  uncorrected_pvals <- results["uncorrected", ]
  correction_factors <- results["correction_factor", ]

  # Test uniformity for both corrected and uncorrected
  ks_corrected <- ks.test(corrected_pvals, punif)
  ks_uncorrected <- ks.test(uncorrected_pvals, punif)

  type_i_corrected <- mean(corrected_pvals < 0.05)
  type_i_uncorrected <- mean(uncorrected_pvals < 0.05)

  cat("Results:\n")
  cat("  Corrected - KS p-value:", format(ks_corrected$p.value, digits = 4), "\n")
  cat("  Corrected - Type I rate:", round(type_i_corrected, 3), "\n")
  cat("  Uncorrected - KS p-value:", format(ks_uncorrected$p.value, digits = 4), "\n")
  cat("  Uncorrected - Type I rate:", round(type_i_uncorrected, 3), "\n")
  cat("  Mean correction factor:", round(mean(correction_factors), 2), "\n")

  return(list(
    corrected_pvals = corrected_pvals,
    uncorrected_pvals = uncorrected_pvals,
    correction_factors = correction_factors,
    corrected_uniform = ks_corrected$p.value > 0.05,
    corrected_type_i_controlled = abs(type_i_corrected - 0.05) < 0.02,
    improvement = type_i_uncorrected - type_i_corrected
  ))
}

#' Compare power between methods
#'
#' Evaluates detection power for different copy-number weighted approaches
#' using simulated enriched pathways.
#'
#' @param methods Vector of method names to compare
#' @param enrichment_factors Vector of enrichment levels to test
#' @param n_simulations Number of power simulations per scenario
compare_method_power <- function(methods = c("standard", "permutation", "effective_sample"),
                                enrichment_factors = c(1.5, 2, 3, 5),
                                n_simulations = 200,
                                background_size = 1000,
                                pathway_size = 100,
                                query_size = 50) {

  cat("Comparing statistical power between methods...\n")

  # Generate background with realistic copy numbers
  background_df <- data.frame(
    gene = paste0("G", 1:background_size),
    copy_number = sample(c(rep(1:2, each = background_size * 0.4),
                          rep(3:8, each = background_size * 0.1)),
                        background_size, replace = TRUE),
    stringsAsFactors = FALSE
  )

  pathway_genes <- sample(background_df$gene, pathway_size)

  results <- expand.grid(method = methods, enrichment = enrichment_factors,
                        stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(
      power = {
        # Create enriched query set
        n_pathway_in_query <- round(pathway_size * query_size / background_size * enrichment)
        n_pathway_in_query <- min(n_pathway_in_query, min(pathway_size, query_size))
        n_random_in_query <- query_size - n_pathway_in_query

        pvals <- replicate(n_simulations, {
          # Create enriched query
          query_genes <- c(sample(pathway_genes, n_pathway_in_query),
                          sample(setdiff(background_df$gene, pathway_genes), n_random_in_query))
          query_df <- background_df[background_df$gene %in% query_genes, ]

          if (method == "standard") {
            # Standard hypergeometric test
            overlap <- sum(query_genes %in% pathway_genes)
            phyper(overlap - 1, pathway_size, background_size - pathway_size,
                   query_size, lower.tail = FALSE)
          } else if (method == "permutation") {
            result <- permutation_weighted_test(query_df, pathway_genes, background_df,
                                               n_permutations = 1000)
            result$pvalue
          } else if (method == "effective_sample") {
            result <- effective_sample_corrected_test(query_df, pathway_genes, background_df,
                                                     correction_method = "mean_copies")
            result$pvalue
          }
        })

        mean(pvals < 0.05)  # Power = P(reject H0 | H1 true)
      }
    ) %>%
    ungroup()

  # Print results
  cat("\nPower Analysis Results:\n")
  print(results)

  return(results)
}

#' Run comprehensive enhanced validation
#'
#' Executes all enhanced validation components and generates summary report
#'
#' @param save_results Whether to save detailed results to files
#' @return Comprehensive validation results
run_enhanced_validation <- function(save_results = TRUE) {

  cat("================================================================================\n")
  cat("ENHANCED STATISTICAL VALIDATION FRAMEWORK\n")
  cat("Copy-Number Weighted Methods with Statistical Corrections\n")
  cat("================================================================================\n\n")

  start_time <- Sys.time()

  validation_results <- list(
    metadata = list(
      timestamp = start_time,
      framework_version = "enhanced-v1"
    ),
    permutation_validation = list(),
    correction_validation = list(),
    power_comparison = NULL,
    summary = list()
  )

  # 1. Validate permutation-based approach
  cat("PHASE 1: PERMUTATION-BASED METHOD VALIDATION\n")
  cat("============================================\n\n")

  perm_scenarios <- expand.grid(
    copy_dist = c("uniform", "skewed", "realistic"),
    bg_size = c(500, 1000),
    query_size = c(25, 50),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(perm_scenarios)) {
    scenario <- perm_scenarios[i, ]
    cat("Testing scenario:", scenario$copy_dist, scenario$bg_size, scenario$query_size, "\n")

    result <- validate_permutation_null_distribution(
      n_simulations = 300,
      n_permutations = 500,
      background_size = scenario$bg_size,
      pathway_size = scenario$bg_size / 10,
      query_size = scenario$query_size,
      copy_distribution = scenario$copy_dist
    )

    validation_results$permutation_validation[[paste0("scenario_", i)]] <- result
    cat("\n")
  }

  # 2. Validate effective sample size corrections
  cat("PHASE 2: EFFECTIVE SAMPLE SIZE CORRECTION VALIDATION\n")
  cat("==================================================\n\n")

  correction_methods <- c("mean_copies", "variance_adjustment")

  for (method in correction_methods) {
    cat("Testing correction method:", method, "\n")

    result <- validate_effective_sample_correction(
      correction_method = method,
      n_simulations = 300,
      background_size = 1000,
      pathway_size = 100,
      query_size = 50,
      copy_distribution = "uniform"
    )

    validation_results$correction_validation[[method]] <- result
    cat("\n")
  }

  # 3. Compare power between methods
  cat("PHASE 3: POWER COMPARISON ANALYSIS\n")
  cat("=================================\n\n")

  power_results <- compare_method_power(
    methods = c("standard", "permutation", "effective_sample"),
    enrichment_factors = c(1.5, 2.0, 3.0),
    n_simulations = 100
  )

  validation_results$power_comparison <- power_results

  # 4. Generate summary
  perm_pass_rate <- mean(sapply(validation_results$permutation_validation,
                               function(x) x$uniform && x$type_i_controlled))

  correction_pass_rate <- mean(sapply(validation_results$correction_validation,
                                     function(x) x$corrected_uniform && x$corrected_type_i_controlled))

  validation_results$summary <- list(
    total_runtime = as.numeric(difftime(Sys.time(), start_time, units = "mins")),
    permutation_pass_rate = perm_pass_rate,
    correction_pass_rate = correction_pass_rate,
    overall_status = if (perm_pass_rate > 0.8) "PASS" else "REVIEW_NEEDED"
  )

  cat("\n================================================================================\n")
  cat("ENHANCED VALIDATION SUMMARY\n")
  cat("================================================================================\n")
  cat("Permutation method validation pass rate:", round(perm_pass_rate * 100, 1), "%\n")
  cat("Correction method validation pass rate:", round(correction_pass_rate * 100, 1), "%\n")
  cat("Overall status:", validation_results$summary$overall_status, "\n")
  cat("Total runtime:", round(validation_results$summary$total_runtime, 2), "minutes\n")

  if (save_results) {
    save(validation_results, file = "enhanced_validation_results.RData")
    cat("Results saved to: enhanced_validation_results.RData\n")
  }

  return(validation_results)
}

# ==============================================================================
# DEMONSTRATION AND USAGE
# ==============================================================================

if (FALSE) {  # Set to TRUE to run demonstrations

  # Example 1: Basic permutation test
  demo_data <- data.frame(
    gene = paste0("Gene", 1:100),
    copy_number = sample(1:8, 100, replace = TRUE)
  )

  pathway <- sample(demo_data$gene, 20)
  query <- sample(demo_data$gene, 30)
  query_df <- demo_data[demo_data$gene %in% query, ]

  perm_result <- permutation_weighted_test(query_df, pathway, demo_data)
  print(perm_result)

  # Example 2: Effective sample size correction
  corr_result <- effective_sample_corrected_test(query_df, pathway, demo_data)
  print(corr_result)

  # Example 3: Run enhanced validation
  validation_results <- run_enhanced_validation()
}

cat("Enhanced Statistical Validation Framework loaded successfully!\n")
cat("Key functions:\n")
cat("  - permutation_weighted_test(): Statistically valid permutation-based test\n")
cat("  - effective_sample_corrected_test(): Parametric test with clustering correction\n")
cat("  - run_enhanced_validation(): Comprehensive validation of corrected methods\n")