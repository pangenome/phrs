#!/usr/bin/env Rscript

#' Copy-Number-Weighted Hypergeometric Parameter Mapping for ORA
#'
#' This script implements concrete parameter mapping functions from copy-number
#' weighted gene counts to R phyper(q,m,n,k) parameters, based on comprehensive
#' research into mathematical framework and statistical validity.
#'
#' Author: Workgraph Agent (fix-implement-copy task)
#' Date: 2026-04-01
#'
#' Research Foundation:
#' - r_phyper_modifications_research.md
#' - phyper_parameter_modification_analysis.md
#' - mathematical_framework_copy_weighted_ora.md

# Load required libraries
suppressPackageStartupMessages({
  library(utils)
})

# ==============================================================================
# CORE PARAMETER MAPPING FUNCTIONS
# ==============================================================================

#' Calculate Copy-Number-Weighted Hypergeometric Parameters
#'
#' Transforms standard gene counts to copy-number-weighted gene instance counts
#' for use with R's phyper() function. This implements the parameter weighting
#' approach which is mathematically equivalent to instance expansion but
#' computationally superior.
#'
#' @param query_df Data frame with columns: gene_name, copy_number
#' @param pathway_genes Character vector of gene names in the pathway of interest
#' @param background_df Data frame with columns: gene_name, copy_number
#' @param validate_params Logical, whether to validate parameter constraints
#'
#' @return List with components:
#'   - k_weighted: Query instances (sum of copy numbers in query)
#'   - q_weighted: Overlap instances (sum of copy numbers in query ∩ pathway)
#'   - m_weighted: Pathway instances in background
#'   - n_weighted: Non-pathway instances in background
#'   - parameters_standard: Standard (gene-count) parameters for comparison
#'   - fold_enrichment: Copy-weighted fold enrichment ratio
#'   - validation: Parameter validation results
#'
#' @details
#' Parameter Transformation Logic:
#' - Standard: Count unique genes
#' - Weighted: Sum copy numbers for gene instances
#'
#' Mathematical relationships:
#' - k_weighted = Σ(copy_number_i) for all i in query
#' - q_weighted = Σ(copy_number_j) for all j in (query ∩ pathway)
#' - m_weighted = Σ(copy_number_k) for all k in (pathway ∩ background)
#' - n_weighted = Σ(all background copy numbers) - m_weighted
#'
#' @examples
#' # Example with PHR data
#' query <- data.frame(
#'   gene_name = c("OR4F17", "OR4F29", "GENE1"),
#'   copy_number = c(14, 14, 5)
#' )
#' pathway <- c("OR4F17", "OR4F29", "OR4F3", "OR4F5")
#' background <- data.frame(
#'   gene_name = paste0("GENE", 1:1000),
#'   copy_number = rep(2, 1000)
#' )
#' params <- calculate_weighted_phyper_params(query, pathway, background)
#'
calculate_weighted_phyper_params <- function(query_df, pathway_genes, background_df,
                                           validate_params = TRUE) {

  # Input validation
  if (!all(c("gene_name", "copy_number") %in% names(query_df))) {
    stop("query_df must have columns: gene_name, copy_number")
  }
  if (!all(c("gene_name", "copy_number") %in% names(background_df))) {
    stop("background_df must have columns: gene_name, copy_number")
  }
  if (length(pathway_genes) == 0) {
    stop("pathway_genes cannot be empty")
  }

  # Remove any zero-copy genes (annotation errors)
  query_df <- query_df[query_df$copy_number > 0, ]
  background_df <- background_df[background_df$copy_number > 0, ]

  if (nrow(query_df) == 0) {
    stop("No valid genes in query_df after filtering zero copies")
  }
  if (nrow(background_df) == 0) {
    stop("No valid genes in background_df after filtering zero copies")
  }

  # ==========================================================================
  # CALCULATE WEIGHTED PARAMETERS
  # ==========================================================================

  # 1. Query size: k_standard → k_weighted
  k_standard <- length(unique(query_df$gene_name))        # unique genes
  k_weighted <- sum(query_df$copy_number)                 # total instances

  # 2. Overlap: q_standard → q_weighted
  query_pathway_genes <- intersect(query_df$gene_name, pathway_genes)
  q_standard <- length(query_pathway_genes)               # unique overlap genes

  if (length(query_pathway_genes) > 0) {
    query_pathway_df <- query_df[query_df$gene_name %in% pathway_genes, ]
    q_weighted <- sum(query_pathway_df$copy_number)       # overlap instances
  } else {
    q_weighted <- 0
  }

  # 3. Pathway size in background: m_standard → m_weighted
  pathway_in_background <- intersect(background_df$gene_name, pathway_genes)
  m_standard <- length(pathway_in_background)             # unique pathway genes

  if (length(pathway_in_background) > 0) {
    pathway_bg_df <- background_df[background_df$gene_name %in% pathway_genes, ]
    m_weighted <- sum(pathway_bg_df$copy_number)          # pathway instances
  } else {
    m_weighted <- 0
  }

  # 4. Non-pathway background: n_standard → n_weighted
  background_total <- sum(background_df$copy_number)      # total background instances
  n_weighted <- background_total - m_weighted             # non-pathway instances
  n_standard <- nrow(background_df) - m_standard          # non-pathway unique genes

  # ==========================================================================
  # PARAMETER VALIDATION
  # ==========================================================================

  validation_results <- list(passed = TRUE, warnings = character(0), errors = character(0))

  if (validate_params) {
    # Check hypergeometric constraints
    validation_checks <- list(
      "Non-negative parameters" = all(c(k_weighted, q_weighted, m_weighted, n_weighted) >= 0),
      "Integer parameters" = all(c(k_weighted, q_weighted, m_weighted, n_weighted) == floor(c(k_weighted, q_weighted, m_weighted, n_weighted))),
      "Overlap <= query size" = q_weighted <= k_weighted,
      "Overlap <= pathway size" = q_weighted <= m_weighted,
      "Query <= total population" = k_weighted <= (m_weighted + n_weighted)
    )

    failed_checks <- names(validation_checks)[!unlist(validation_checks)]
    if (length(failed_checks) > 0) {
      validation_results$passed <- FALSE
      validation_results$errors <- failed_checks
    }

    # Additional warnings
    if (k_weighted < 10) {
      validation_results$warnings <- c(validation_results$warnings,
                                     "Small query size may reduce statistical power")
    }
    if (m_weighted == 0) {
      validation_results$warnings <- c(validation_results$warnings,
                                     "No pathway genes found in background")
    }
    if (q_weighted == 0) {
      validation_results$warnings <- c(validation_results$warnings,
                                     "No overlap between query and pathway")
    }
  }

  # ==========================================================================
  # CALCULATE ENRICHMENT METRICS
  # ==========================================================================

  # Copy-weighted fold enrichment
  query_pathway_fraction <- if (k_weighted > 0) q_weighted / k_weighted else 0
  background_pathway_fraction <- if ((m_weighted + n_weighted) > 0) {
    m_weighted / (m_weighted + n_weighted)
  } else 0

  fold_enrichment_weighted <- if (background_pathway_fraction > 0) {
    query_pathway_fraction / background_pathway_fraction
  } else Inf

  # Standard fold enrichment for comparison
  fold_enrichment_standard <- if (m_standard > 0 && (m_standard + n_standard) > 0 && k_standard > 0) {
    (q_standard / k_standard) / (m_standard / (m_standard + n_standard))
  } else NA

  # ==========================================================================
  # RETURN RESULTS
  # ==========================================================================

  return(list(
    # Primary weighted parameters for phyper()
    k_weighted = k_weighted,    # query instances
    q_weighted = q_weighted,    # overlap instances
    m_weighted = m_weighted,    # pathway instances in background
    n_weighted = n_weighted,    # non-pathway instances in background

    # Standard parameters for comparison
    parameters_standard = list(
      k_standard = k_standard,
      q_standard = q_standard,
      m_standard = m_standard,
      n_standard = n_standard
    ),

    # Enrichment metrics
    fold_enrichment_weighted = fold_enrichment_weighted,
    fold_enrichment_standard = fold_enrichment_standard,

    # Validation results
    validation = validation_results,

    # Additional metadata
    metadata = list(
      total_background_instances = background_total,
      query_genes_count = nrow(query_df),
      pathway_genes_count = length(pathway_genes),
      background_genes_count = nrow(background_df),
      copy_expansion_factor = k_weighted / k_standard,
      pathway_overlap_genes = query_pathway_genes
    )
  ))
}

#' Run Copy-Number-Weighted Hypergeometric Test
#'
#' Wrapper function that calculates weighted parameters and runs phyper() test
#'
#' @param query_df Data frame with gene_name and copy_number columns
#' @param pathway_genes Character vector of pathway gene names
#' @param background_df Data frame with gene_name and copy_number columns
#' @param alpha Significance level for interpretation
#'
#' @return List with test results including p-value and interpretation
#'
run_weighted_hypergeometric_test <- function(query_df, pathway_genes, background_df,
                                           alpha = 0.05) {

  # Calculate weighted parameters
  params <- calculate_weighted_phyper_params(query_df, pathway_genes, background_df)

  # Check if test is feasible
  if (!params$validation$passed) {
    stop(paste("Parameter validation failed:", paste(params$validation$errors, collapse = "; ")))
  }

  if (params$q_weighted == 0) {
    return(list(
      pvalue = 1.0,
      significant = FALSE,
      method = "copy_weighted_hypergeometric",
      note = "No overlap between query and pathway"
    ))
  }

  # Run hypergeometric test
  # P(X >= q_weighted) where X ~ Hypergeometric(k_weighted, m_weighted, n_weighted)
  pvalue <- phyper(params$q_weighted - 1, params$m_weighted, params$n_weighted,
                   params$k_weighted, lower.tail = FALSE)

  # Calculate expected overlap under null hypothesis
  expected_overlap <- params$k_weighted * params$m_weighted / (params$m_weighted + params$n_weighted)

  return(list(
    pvalue = pvalue,
    pvalue_adj = NA,  # To be filled by multiple testing correction
    significant = pvalue < alpha,
    observed_overlap_weighted = params$q_weighted,
    expected_overlap_weighted = expected_overlap,
    fold_enrichment = params$fold_enrichment_weighted,
    query_size_weighted = params$k_weighted,
    pathway_size_weighted = params$m_weighted,
    background_size_weighted = params$m_weighted + params$n_weighted,
    method = "copy_weighted_hypergeometric",
    parameters = params,
    warnings = params$validation$warnings
  ))
}

#' Verify Mathematical Equivalence with Instance Expansion
#'
#' Validates that parameter weighting produces identical results to instance expansion
#'
#' @param query_df Data frame with gene_name and copy_number
#' @param pathway_genes Character vector of pathway genes
#' @param background_df Data frame with gene_name and copy_number
#' @param tolerance Numerical tolerance for comparison
#'
#' @return List with equivalence test results
#'
verify_equivalence_with_expansion <- function(query_df, pathway_genes, background_df,
                                            tolerance = 1e-12) {

  cat("Verifying mathematical equivalence with instance expansion...\n")

  # Method 1: Parameter weighting approach
  params <- calculate_weighted_phyper_params(query_df, pathway_genes, background_df)
  pval_weighted <- phyper(params$q_weighted - 1, params$m_weighted,
                          params$n_weighted, params$k_weighted, lower.tail = FALSE)

  # Method 2: Instance expansion approach
  query_expanded <- rep(query_df$gene_name, query_df$copy_number)
  background_expanded <- rep(background_df$gene_name, background_df$copy_number)

  q_exp <- sum(query_expanded %in% pathway_genes)
  m_exp <- sum(background_expanded %in% pathway_genes)
  n_exp <- length(background_expanded) - m_exp
  k_exp <- length(query_expanded)

  pval_expansion <- phyper(q_exp - 1, m_exp, n_exp, k_exp, lower.tail = FALSE)

  # Compare parameters
  params_match <- all(c(
    params$q_weighted == q_exp,
    params$m_weighted == m_exp,
    params$n_weighted == n_exp,
    params$k_weighted == k_exp
  ))

  # Compare p-values
  pval_diff <- abs(pval_weighted - pval_expansion)
  pvals_match <- pval_diff < tolerance

  return(list(
    parameters_equivalent = params_match,
    pvalues_equivalent = pvals_match,
    pvalue_difference = pval_diff,
    weighted_pvalue = pval_weighted,
    expansion_pvalue = pval_expansion,
    weighted_params = c(params$q_weighted, params$m_weighted, params$n_weighted, params$k_weighted),
    expansion_params = c(q_exp, m_exp, n_exp, k_exp),
    memory_reduction_factor = length(background_expanded) / nrow(background_df)
  ))
}

# ==============================================================================
# VALIDATION AND EXAMPLE USAGE
# ==============================================================================

#' Load PHR Gene Copy Data
#'
#' Helper function to load and prepare PHR gene copy data for analysis
#'
load_phr_copy_data <- function() {
  if (!file.exists("gene_copy_summary.csv")) {
    stop("PHR gene copy data file 'gene_copy_summary.csv' not found")
  }

  copy_data <- read.csv("gene_copy_summary.csv", stringsAsFactors = FALSE)

  # Focus on protein-coding genes
  protein_genes <- copy_data[copy_data$gene_biotype == "protein_coding", ]

  # Prepare data frame in expected format
  phr_data <- data.frame(
    gene_name = protein_genes$gene_name,
    copy_number = protein_genes$total_copies,
    stringsAsFactors = FALSE
  )

  return(phr_data)
}

#' Create Mock Background Dataset
#'
#' Generate a simulated genome-wide background for testing
#'
create_mock_background <- function(n_genes = 20000, copy_distribution = c(0.7, 0.15, 0.1, 0.03, 0.02)) {
  background <- data.frame(
    gene_name = paste0("BG_GENE", 1:n_genes),
    copy_number = sample(1:5, n_genes, replace = TRUE, prob = copy_distribution),
    stringsAsFactors = FALSE
  )
  return(background)
}

#' Example Analysis: Olfactory Receptor Enrichment
#'
#' Demonstrates the copy-number-weighted ORA approach using PHR data
#'
example_olfactory_analysis <- function() {
  cat("=== Copy-Number-Weighted ORA Example: Olfactory Receptors ===\n\n")

  # Load PHR data (if available)
  if (file.exists("gene_copy_summary.csv")) {
    cat("Loading PHR gene copy data...\n")
    query_data <- load_phr_copy_data()
    cat(sprintf("Loaded %d protein-coding genes with copy data\n", nrow(query_data)))
  } else {
    cat("PHR data not available, using mock query data...\n")
    query_data <- data.frame(
      gene_name = c("OR4F17", "OR4F29", "OR4F3", "GENE1", "GENE2"),
      copy_number = c(14, 14, 10, 5, 3),
      stringsAsFactors = FALSE
    )
  }

  # Define olfactory receptor pathway
  olfactory_genes <- c("OR4F17", "OR4F29", "OR4F3", "OR4F5", "OR1A1", "OR1A2")
  cat(sprintf("Testing pathway with %d olfactory receptor genes\n", length(olfactory_genes)))

  # Create background dataset
  cat("Generating background dataset...\n")
  background_data <- create_mock_background(n_genes = 20000)

  # Add olfactory genes to background with copy numbers that ensure valid constraints
  # For genes that overlap with query, use at least the query copy number
  query_or_genes <- intersect(query_data$gene_name, olfactory_genes)
  or_background_copies <- numeric(length(olfactory_genes))

  for (i in seq_along(olfactory_genes)) {
    gene <- olfactory_genes[i]
    if (gene %in% query_or_genes) {
      query_copies <- query_data$copy_number[query_data$gene_name == gene][1]
      # Background should have at least as many copies as query for overlapping genes
      or_background_copies[i] <- max(query_copies, sample(c(query_copies:(query_copies+5)), 1))
    } else {
      # Non-overlapping pathway genes can have random copies
      or_background_copies[i] <- sample(1:10, 1, prob = c(0.3, 0.2, 0.15, 0.1, 0.08, 0.07, 0.05, 0.03, 0.01, 0.01))
    }
  }

  or_background <- data.frame(
    gene_name = olfactory_genes,
    copy_number = or_background_copies,
    stringsAsFactors = FALSE
  )
  background_data <- rbind(background_data, or_background)

  cat(sprintf("Background: %d genes, %d total instances\n",
              nrow(background_data), sum(background_data$copy_number)))

  # Calculate weighted parameters
  cat("\nCalculating copy-number-weighted parameters...\n")
  params <- calculate_weighted_phyper_params(query_data, olfactory_genes, background_data)

  # Display parameter mapping
  cat("\nParameter Transformation Results:\n")
  cat(sprintf("Query size:      %d genes → %d instances (%.1fx expansion)\n",
              params$parameters_standard$k_standard, params$k_weighted,
              params$k_weighted / params$parameters_standard$k_standard))
  cat(sprintf("Overlap:         %d genes → %d instances\n",
              params$parameters_standard$q_standard, params$q_weighted))
  cat(sprintf("Pathway size:    %d genes → %d instances\n",
              params$parameters_standard$m_standard, params$m_weighted))
  cat(sprintf("Background size: %d genes → %d instances\n",
              params$parameters_standard$m_standard + params$parameters_standard$n_standard,
              params$m_weighted + params$n_weighted))

  # Run hypergeometric test
  cat("\nRunning weighted hypergeometric test...\n")
  test_result <- run_weighted_hypergeometric_test(query_data, olfactory_genes, background_data)

  cat(sprintf("P-value: %.2e\n", test_result$pvalue))
  cat(sprintf("Fold enrichment: %.2f\n", test_result$fold_enrichment))
  cat(sprintf("Significant (α=0.05): %s\n", ifelse(test_result$significant, "YES", "NO")))

  # Verify equivalence with expansion method
  cat("\nVerifying equivalence with instance expansion...\n")
  equivalence <- verify_equivalence_with_expansion(query_data, olfactory_genes, background_data)
  cat(sprintf("Parameters equivalent: %s\n", equivalence$parameters_equivalent))
  cat(sprintf("P-values equivalent: %s (diff = %.2e)\n",
              equivalence$pvalues_equivalent, equivalence$pvalue_difference))
  cat(sprintf("Memory reduction: %.1fx less memory than expansion\n",
              equivalence$memory_reduction_factor))

  return(list(
    parameters = params,
    test_result = test_result,
    equivalence = equivalence
  ))
}

# ==============================================================================
# MAIN EXECUTION (if run as script)
# ==============================================================================

if (!interactive()) {
  cat("Copy-Number-Weighted Hypergeometric Parameter Mapping\n")
  cat("=====================================================\n\n")

  # Run example analysis
  example_results <- example_olfactory_analysis()

  cat("\n=== Example Analysis Complete ===\n")
  cat("Functions available for use:\n")
  cat("- calculate_weighted_phyper_params()\n")
  cat("- run_weighted_hypergeometric_test()\n")
  cat("- verify_equivalence_with_expansion()\n")
}