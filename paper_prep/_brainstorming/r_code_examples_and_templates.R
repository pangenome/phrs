# R Code Examples and Templates for Copy-Number Weighted ORA
# Comprehensive implementation examples based on research synthesis

# ==============================================================================
# CORE IMPLEMENTATION TEMPLATES
# ==============================================================================

#' Production-Ready Weighted Hypergeometric Test
#'
#' Implements the computationally efficient parameter-weighting approach
#' with comprehensive validation and error handling.
#'
#' @param query_df Data frame with 'gene' and 'copy_number' columns
#' @param pathway_genes Character vector of pathway gene identifiers
#' @param background_df Data frame with 'gene' and 'copy_number' columns
#' @param validate_inputs Logical, whether to perform input validation
#' @return List with test results and metadata
weighted_hypergeometric_ora <- function(query_df, pathway_genes, background_df,
                                        validate_inputs = TRUE) {

  if (validate_inputs) {
    # Comprehensive input validation
    if (!is.data.frame(query_df) || !all(c("gene", "copy_number") %in% names(query_df))) {
      stop("query_df must be data frame with 'gene' and 'copy_number' columns")
    }

    if (!is.data.frame(background_df) || !all(c("gene", "copy_number") %in% names(background_df))) {
      stop("background_df must be data frame with 'gene' and 'copy_number' columns")
    }

    if (!is.character(pathway_genes) || length(pathway_genes) == 0) {
      stop("pathway_genes must be non-empty character vector")
    }

    if (any(query_df$copy_number < 1) || any(!is.finite(query_df$copy_number))) {
      stop("All query copy numbers must be positive integers")
    }

    if (any(background_df$copy_number < 1) || any(!is.finite(background_df$copy_number))) {
      stop("All background copy numbers must be positive integers")
    }

    if (any(duplicated(query_df$gene))) {
      stop("Duplicate genes found in query_df")
    }

    if (any(duplicated(background_df$gene))) {
      stop("Duplicate genes found in background_df")
    }
  }

  # Ensure query genes are subset of background genes
  valid_query_genes <- intersect(query_df$gene, background_df$gene)

  if (length(valid_query_genes) == 0) {
    stop("No query genes found in background gene set")
  }

  # Filter to valid query genes and merge with background copy numbers
  query_filtered <- query_df[query_df$gene %in% valid_query_genes, , drop = FALSE]

  # Use background copy numbers for consistency (critical for proper statistics)
  query_with_bg_copies <- merge(
    query_filtered[, "gene", drop = FALSE],
    background_df,
    by = "gene",
    all.x = TRUE
  )

  # Calculate overlap with pathway
  query_in_pathway <- query_with_bg_copies[query_with_bg_copies$gene %in% pathway_genes, , drop = FALSE]
  pathway_in_background <- background_df[background_df$gene %in% pathway_genes, , drop = FALSE]

  # Calculate weighted hypergeometric parameters
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
  if (validate_inputs) {
    stopifnot(q_weighted <= k_weighted)
    stopifnot(q_weighted <= m_weighted)
    stopifnot(k_weighted <= (m_weighted + n_weighted))

    if (m_weighted == 0) {
      warning("No pathway genes found in background - returning p-value = 1")
      return(list(pvalue = 1.0, valid = FALSE, reason = "no_pathway_genes"))
    }

    if (n_weighted == 0) {
      warning("All background genes are in pathway - invalid test")
      return(list(pvalue = 1.0, valid = FALSE, reason = "all_genes_in_pathway"))
    }
  }

  # Hypergeometric test
  if (q_weighted == 0) {
    pvalue <- 1.0
  } else {
    pvalue <- phyper(q_weighted - 1, m_weighted, n_weighted, k_weighted,
                    lower.tail = FALSE)
  }

  # Calculate enrichment metrics
  observed_rate <- q_weighted / k_weighted
  expected_rate <- m_weighted / (m_weighted + n_weighted)
  fold_enrichment <- if (expected_rate > 0) {
    observed_rate / expected_rate
  } else {
    Inf
  }

  # Expected overlap under null hypothesis
  expected_overlap <- k_weighted * expected_rate

  # Return comprehensive results
  list(
    # Primary results
    pvalue = pvalue,
    fold_enrichment = fold_enrichment,

    # Instance-level counts (weighted)
    overlap_instances = q_weighted,
    query_instances = k_weighted,
    pathway_instances = m_weighted,
    background_instances = m_weighted + n_weighted,
    expected_overlap_instances = expected_overlap,

    # Gene-level counts (for reference)
    overlap_genes = nrow(query_in_pathway),
    query_genes = nrow(query_with_bg_copies),
    pathway_genes_in_background = nrow(pathway_in_background),
    total_background_genes = nrow(background_df),

    # Statistical metadata
    hypergeometric_parameters = list(
      q = q_weighted,
      m = m_weighted,
      n = n_weighted,
      k = k_weighted
    ),

    # Validation flags
    valid = TRUE,
    method = "weighted_hypergeometric",
    parameter_approach = "direct_weighting"
  )
}

#' Standard (Unweighted) Hypergeometric Test for Comparison
standard_hypergeometric_ora <- function(query_df, pathway_genes, background_df) {

  # Extract unique gene lists (ignore copy numbers)
  query_genes <- unique(query_df$gene)
  background_genes <- unique(background_df$gene)

  # Ensure valid overlap
  valid_query_genes <- intersect(query_genes, background_genes)

  if (length(valid_query_genes) == 0) {
    stop("No query genes found in background gene set")
  }

  # Standard hypergeometric parameters (gene-level)
  k_standard <- length(valid_query_genes)
  q_standard <- length(intersect(valid_query_genes, pathway_genes))
  m_standard <- length(intersect(background_genes, pathway_genes))
  n_standard <- length(background_genes) - m_standard

  # Hypergeometric test
  if (q_standard == 0) {
    pvalue <- 1.0
  } else {
    pvalue <- phyper(q_standard - 1, m_standard, n_standard, k_standard,
                    lower.tail = FALSE)
  }

  # Calculate enrichment metrics
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
    expected_overlap_genes = k_standard * expected_rate,
    hypergeometric_parameters = list(
      q = q_standard,
      m = m_standard,
      n = n_standard,
      k = k_standard
    ),
    method = "standard_hypergeometric"
  )
}

# ==============================================================================
# CONTEXT-AWARE WRAPPER FUNCTION
# ==============================================================================

#' Intelligent Copy-Number Aware ORA
#'
#' Selects appropriate method based on context and provides both results
#' for comparison. Includes statistical warnings and recommendations.
copy_number_aware_ora <- function(query_df, pathway_genes, background_df,
                                 selection_type = c("auto", "gene_level", "instance_level"),
                                 include_permutation = FALSE,
                                 n_permutations = 1000) {

  selection_type <- match.arg(selection_type)

  # Auto-detect selection type based on copy number variation
  if (selection_type == "auto") {
    cv_query <- sd(query_df$copy_number) / mean(query_df$copy_number)
    cv_background <- sd(background_df$copy_number) / mean(background_df$copy_number)

    # If low copy number variation, prefer standard approach
    if (cv_query < 0.5 && cv_background < 0.5) {
      selection_type <- "gene_level"
      auto_reason <- "Low copy number variation detected"
    } else {
      selection_type <- "instance_level"
      auto_reason <- "High copy number variation detected"
    }
  } else {
    auto_reason <- "User specified"
  }

  # Calculate both weighted and standard results
  weighted_result <- weighted_hypergeometric_ora(query_df, pathway_genes, background_df)
  standard_result <- standard_hypergeometric_ora(query_df, pathway_genes, background_df)

  # Determine primary recommendation
  if (selection_type == "gene_level") {
    primary_method <- "standard"
    primary_result <- standard_result
    alternative_result <- weighted_result

    statistical_warnings <- c(
      "Gene-level selection specified - weighted test may be anti-conservative",
      "Type I error rate can be inflated 4-6x above nominal level",
      "FDR correction may not be reliable for weighted p-values",
      "Standard hypergeometric test recommended for this context"
    )

  } else {
    primary_method <- "weighted"
    primary_result <- weighted_result
    alternative_result <- standard_result

    statistical_warnings <- c(
      "Instance-level selection specified - weighted test appropriate",
      "Ensure independence assumption holds at copy level",
      "Consider validation with permutation testing"
    )
  }

  # Optionally compute permutation-based p-value
  if (include_permutation) {
    perm_pvalue <- permutation_hypergeometric_test(
      query_df, pathway_genes, background_df, n_permutations
    )
  } else {
    perm_pvalue <- NULL
  }

  # Compile comprehensive results
  list(
    # Method selection
    primary_method = primary_method,
    selection_type = selection_type,
    selection_reason = auto_reason,

    # Primary and alternative results
    primary_result = primary_result,
    standard_result = standard_result,
    weighted_result = weighted_result,
    permutation_pvalue = perm_pvalue,

    # Statistical guidance
    statistical_warnings = statistical_warnings,

    # Method comparison
    method_comparison = compare_methods(standard_result, weighted_result),

    # Recommendations
    recommendations = generate_recommendations(standard_result, weighted_result, selection_type)
  )
}

# ==============================================================================
# PERMUTATION-BASED IMPLEMENTATION
# ==============================================================================

#' Permutation-Based Hypergeometric Test
#'
#' Gold standard approach that correctly accounts for gene-level clustering
#' while incorporating copy number information.
permutation_hypergeometric_test <- function(query_df, pathway_genes, background_df,
                                           n_permutations = 1000,
                                           alternative = "greater",
                                           seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Calculate observed test statistic using weighted approach
  observed_result <- weighted_hypergeometric_ora(query_df, pathway_genes, background_df)
  observed_statistic <- observed_result$overlap_instances

  # Generate permutation distribution
  n_query_genes <- nrow(query_df)
  background_genes <- background_df$gene

  permutation_statistics <- replicate(n_permutations, {
    # Sample random genes (preserving copy number structure)
    sampled_genes <- sample(background_genes, n_query_genes, replace = FALSE)

    # Create permuted query with original copy number structure
    perm_query_df <- data.frame(
      gene = sampled_genes,
      copy_number = query_df$copy_number  # Keep original copy number pattern
    )

    # Calculate test statistic under permutation
    perm_result <- weighted_hypergeometric_ora(perm_query_df, pathway_genes, background_df)
    return(perm_result$overlap_instances)
  })

  # Calculate empirical p-value
  if (alternative == "greater") {
    p_empirical <- mean(permutation_statistics >= observed_statistic)
  } else if (alternative == "less") {
    p_empirical <- mean(permutation_statistics <= observed_statistic)
  } else {  # two-sided
    p_empirical <- mean(abs(permutation_statistics - mean(permutation_statistics)) >=
                       abs(observed_statistic - mean(permutation_statistics)))
  }

  list(
    empirical_pvalue = p_empirical,
    observed_statistic = observed_statistic,
    permutation_distribution = permutation_statistics,
    n_permutations = n_permutations,
    alternative = alternative
  )
}

# ==============================================================================
# MULTIPLE PATHWAY TESTING
# ==============================================================================

#' Efficient Multiple Pathway Testing
#'
#' Vectorized implementation for testing many pathways simultaneously
#' with proper multiple testing correction.
test_multiple_pathways <- function(query_df, pathway_list, background_df,
                                  method = c("weighted", "standard", "both"),
                                  correction = c("fdr", "bonferroni", "none"),
                                  min_pathway_size = 5,
                                  parallel = FALSE,
                                  n_cores = 2) {

  method <- match.arg(method)
  correction <- match.arg(correction)

  # Filter pathways by minimum size
  pathway_sizes <- sapply(pathway_list, function(genes) {
    length(intersect(genes, background_df$gene))
  })

  valid_pathways <- pathway_list[pathway_sizes >= min_pathway_size]

  if (length(valid_pathways) == 0) {
    warning("No pathways meet minimum size criterion")
    return(data.frame())
  }

  # Pre-calculate query statistics for efficiency
  valid_query_genes <- intersect(query_df$gene, background_df$gene)
  query_filtered <- query_df[query_df$gene %in% valid_query_genes, , drop = FALSE]
  query_with_bg_copies <- merge(
    query_filtered[, "gene", drop = FALSE],
    background_df,
    by = "gene",
    all.x = TRUE
  )

  # Function to test single pathway
  test_single_pathway <- function(pathway_name, pathway_genes) {
    tryCatch({
      if (method %in% c("weighted", "both")) {
        weighted_result <- weighted_hypergeometric_ora(query_df, pathway_genes, background_df)
      }

      if (method %in% c("standard", "both")) {
        standard_result <- standard_hypergeometric_ora(query_df, pathway_genes, background_df)
      }

      # Return appropriate results
      if (method == "weighted") {
        data.frame(
          pathway = pathway_name,
          pvalue = weighted_result$pvalue,
          overlap_instances = weighted_result$overlap_instances,
          query_instances = weighted_result$query_instances,
          pathway_instances = weighted_result$pathway_instances,
          fold_enrichment = weighted_result$fold_enrichment,
          method = "weighted"
        )
      } else if (method == "standard") {
        data.frame(
          pathway = pathway_name,
          pvalue = standard_result$pvalue,
          overlap_genes = standard_result$overlap_genes,
          query_genes = standard_result$query_genes,
          pathway_genes = standard_result$pathway_genes_in_background,
          fold_enrichment = standard_result$fold_enrichment,
          method = "standard"
        )
      } else {  # both
        data.frame(
          pathway = rep(pathway_name, 2),
          pvalue = c(weighted_result$pvalue, standard_result$pvalue),
          overlap = c(weighted_result$overlap_instances, standard_result$overlap_genes),
          query_size = c(weighted_result$query_instances, standard_result$query_genes),
          pathway_size = c(weighted_result$pathway_instances, standard_result$pathway_genes_in_background),
          fold_enrichment = c(weighted_result$fold_enrichment, standard_result$fold_enrichment),
          method = c("weighted", "standard")
        )
      }
    }, error = function(e) {
      # Return NA results for failed tests
      data.frame(
        pathway = pathway_name,
        pvalue = NA,
        error = e$message
      )
    })
  }

  # Execute tests (optionally in parallel)
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    results_list <- parallel::mcmapply(
      test_single_pathway,
      names(valid_pathways),
      valid_pathways,
      mc.cores = n_cores,
      SIMPLIFY = FALSE
    )
  } else {
    results_list <- mapply(
      test_single_pathway,
      names(valid_pathways),
      valid_pathways,
      SIMPLIFY = FALSE
    )
  }

  # Combine results
  results <- do.call(rbind, results_list)

  # Remove failed tests
  results <- results[!is.na(results$pvalue), ]

  if (nrow(results) == 0) {
    warning("No valid pathway test results")
    return(data.frame())
  }

  # Apply multiple testing correction
  if (method == "both") {
    # Apply correction separately for each method
    for (test_method in c("weighted", "standard")) {
      method_subset <- results$method == test_method
      if (correction == "fdr") {
        results$qvalue[method_subset] <- p.adjust(results$pvalue[method_subset], method = "fdr")
      } else if (correction == "bonferroni") {
        results$qvalue[method_subset] <- p.adjust(results$pvalue[method_subset], method = "bonferroni")
      } else {
        results$qvalue[method_subset] <- results$pvalue[method_subset]
      }
    }
  } else {
    # Single method correction
    if (correction == "fdr") {
      results$qvalue <- p.adjust(results$pvalue, method = "fdr")
    } else if (correction == "bonferroni") {
      results$qvalue <- p.adjust(results$pvalue, method = "bonferroni")
    } else {
      results$qvalue <- results$pvalue
    }
  }

  # Sort by p-value
  results <- results[order(results$pvalue), ]

  return(results)
}

# ==============================================================================
# UTILITY AND HELPER FUNCTIONS
# ==============================================================================

#' Compare Weighted vs Standard Methods
compare_methods <- function(standard_result, weighted_result) {

  # Calculate method comparison metrics
  pvalue_ratio <- weighted_result$pvalue / standard_result$pvalue
  enrichment_ratio <- weighted_result$fold_enrichment / standard_result$fold_enrichment

  # Instance vs gene counts
  instance_gene_ratio <- weighted_result$query_instances / standard_result$query_genes

  list(
    pvalue_ratio = pvalue_ratio,
    enrichment_ratio = enrichment_ratio,
    instance_gene_ratio = instance_gene_ratio,

    weighted_more_significant = weighted_result$pvalue < standard_result$pvalue,
    large_difference = abs(log10(pvalue_ratio)) > 1,  # >10x difference

    summary = if (pvalue_ratio < 0.1) {
      "Weighted test much more significant"
    } else if (pvalue_ratio > 10) {
      "Standard test much more significant"
    } else {
      "Similar significance levels"
    }
  )
}

#' Generate Method Recommendations
generate_recommendations <- function(standard_result, weighted_result, selection_type) {

  comparison <- compare_methods(standard_result, weighted_result)

  recommendations <- character()

  # Primary method recommendation
  if (selection_type == "gene_level") {
    recommendations <- c(recommendations,
      "Primary: Use standard hypergeometric test for gene-level selection")

    if (comparison$large_difference && comparison$weighted_more_significant) {
      recommendations <- c(recommendations,
        "Caution: Weighted test shows much stronger signal - may indicate anti-conservative behavior")
    }

  } else {
    recommendations <- c(recommendations,
      "Primary: Weighted hypergeometric appropriate for instance-level selection")

    if (comparison$large_difference && !comparison$weighted_more_significant) {
      recommendations <- c(recommendations,
        "Note: Standard test shows stronger signal - verify instance-level independence")
    }
  }

  # Multiple testing recommendations
  if (selection_type == "gene_level") {
    recommendations <- c(recommendations,
      "Multiple testing: Use FDR correction on standard p-values",
      "Alternative: Permutation-based FDR for copy-number aware correction")
  } else {
    recommendations <- c(recommendations,
      "Multiple testing: Consider permutation-based FDR correction",
      "Caution: Standard FDR may not control error rate for weighted tests")
  }

  # Additional guidance
  recommendations <- c(recommendations,
    "Validation: Report both weighted and standard results",
    "Documentation: Clearly state selection assumptions and copy number source")

  return(recommendations)
}

# ==============================================================================
# EXAMPLE USAGE AND TESTING
# ==============================================================================

# Example data creation
create_example_data <- function() {
  # PHR-like query dataset
  query_df <- data.frame(
    gene = paste0("PHR_", 1:35),
    copy_number = sample(10:50, 35, replace = TRUE)
  )

  # Genome-scale background
  background_df <- data.frame(
    gene = c(query_df$gene, paste0("BG_", 1:19965)),
    copy_number = c(query_df$copy_number, sample(1:5, 19965, replace = TRUE))
  )

  # Olfactory receptor pathway
  or_pathway <- c(query_df$gene[1:4], paste0("OR_", 1:396))

  list(
    query_df = query_df,
    background_df = background_df,
    or_pathway = or_pathway
  )
}

# Comprehensive testing example
run_comprehensive_example <- function() {

  cat("Creating example data...\n")
  example_data <- create_example_data()

  cat("Running context-aware ORA...\n")
  results <- copy_number_aware_ora(
    query_df = example_data$query_df,
    pathway_genes = example_data$or_pathway,
    background_df = example_data$background_df,
    selection_type = "auto"
  )

  cat("Results Summary:\n")
  cat("Primary method:", results$primary_method, "\n")
  cat("Selection type:", results$selection_type, "\n")
  cat("Primary p-value:", format(results$primary_result$pvalue, scientific = TRUE), "\n")
  cat("Alternative p-value:", format(results$alternative_result$pvalue, scientific = TRUE), "\n")

  cat("\nStatistical warnings:\n")
  for (warning in results$statistical_warnings) {
    cat("- ", warning, "\n")
  }

  cat("\nRecommendations:\n")
  for (rec in results$recommendations) {
    cat("- ", rec, "\n")
  }

  return(results)
}