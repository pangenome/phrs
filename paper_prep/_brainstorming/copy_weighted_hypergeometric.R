# Copy-Number Weighted Hypergeometric Testing for ORA
#
# This implementation provides robust, production-ready functions for performing
# copy-number weighted over-representation analysis using modified hypergeometric
# parameters. Based on comprehensive research showing mathematical equivalence
# between parameter weighting and instance expansion approaches.
#
# Author: Robust R Code Implementation Task
# Date: 2026-04-01
# Version: 1.0

#' Copy-Number Weighted Hypergeometric Test
#'
#' Performs hypergeometric enrichment testing accounting for gene copy numbers.
#' Uses parameter transformation approach which is mathematically equivalent to
#' instance expansion but computationally more efficient.
#'
#' @param query_df Data frame with columns 'gene' and 'copy_number' for query set
#' @param pathway_genes Vector of gene symbols defining the pathway of interest
#' @param background_df Data frame with columns 'gene' and 'copy_number' for background
#' @param validate_inputs Logical, whether to perform comprehensive input validation
#' @param handle_zeros Logical, whether to automatically remove zero-copy genes
#' @param max_copies Integer, maximum allowed copy number (caps extremes)
#' @param min_instances Integer, minimum total instances required for reliable results
#'
#' @return List containing:
#'   - pvalue: Hypergeometric test p-value
#'   - overlap_instances: Copy-weighted overlap count
#'   - query_instances: Total query instances
#'   - pathway_instances: Total pathway instances in background
#'   - background_instances: Total background instances
#'   - fold_enrichment: Fold enrichment ratio
#'   - parameters: Named vector of (q, m, n, k) hypergeometric parameters
#'   - diagnostics: Validation and diagnostic information
#'
#' @examples
#' # Example with PHR-like data
#' query <- data.frame(
#'   gene = c("OR4F17", "OR4F29", "OR4F3", "GENE1", "GENE2"),
#'   copy_number = c(14, 14, 14, 8, 12)
#' )
#'
#' background <- data.frame(
#'   gene = paste0("GENE", 1:1000),
#'   copy_number = rpois(1000, lambda = 2) + 1  # Copy numbers 1-10
#' )
#'
#' pathway <- paste0("OR4F", c(17, 29, 3, 5))
#'
#' result <- weighted_hypergeometric_test(query, pathway, background)
#' print(paste("p-value:", format(result$pvalue, scientific = TRUE)))
#'
#' @export
weighted_hypergeometric_test <- function(query_df,
                                       pathway_genes,
                                       background_df,
                                       validate_inputs = TRUE,
                                       handle_zeros = TRUE,
                                       max_copies = 500,
                                       min_instances = 10) {

  # Store original inputs for diagnostics
  original_query_size <- nrow(query_df)
  original_background_size <- nrow(background_df)

  # Input validation
  if (validate_inputs) {
    validation_result <- validate_hypergeometric_inputs(
      query_df, pathway_genes, background_df,
      max_copies, min_instances
    )
    if (!validation_result$valid) {
      stop(paste("Input validation failed:", validation_result$message))
    }
  }

  # Handle edge cases
  if (handle_zeros) {
    query_df <- filter_zero_copies(query_df)
    background_df <- filter_zero_copies(background_df)
  }

  # Cap extreme copy numbers
  if (!is.null(max_copies) && max_copies > 0) {
    query_df <- cap_extreme_copies(query_df, max_copies)
    background_df <- cap_extreme_copies(background_df, max_copies)
  }

  # Check for empty datasets after filtering
  if (nrow(query_df) == 0) {
    stop("Query dataset is empty after filtering")
  }
  if (nrow(background_df) == 0) {
    stop("Background dataset is empty after filtering")
  }

  # Calculate weighted parameters
  params <- calculate_weighted_parameters(query_df, pathway_genes, background_df)

  # Validate hypergeometric constraints
  constraint_validation <- validate_hypergeometric_constraints(
    params$q_weighted, params$m_weighted, params$n_weighted, params$k_weighted
  )
  if (!constraint_validation$valid) {
    stop(paste("Hypergeometric constraint violation:", constraint_validation$message))
  }

  # Check minimum sample size
  if (!is.null(min_instances) && params$k_weighted < min_instances) {
    warning(paste("Small sample size:", params$k_weighted,
                  "instances. Results may be unreliable."))
  }

  # Perform hypergeometric test
  # Using lower.tail=FALSE for P(X >= q_weighted)
  pvalue <- phyper(params$q_weighted - 1,
                   params$m_weighted,
                   params$n_weighted,
                   params$k_weighted,
                   lower.tail = FALSE)

  # Calculate fold enrichment
  # (observed/expected) = (q_weighted/k_weighted) / (m_weighted/total_weighted)
  total_weighted <- params$m_weighted + params$n_weighted
  expected_fraction <- params$m_weighted / total_weighted
  observed_fraction <- params$q_weighted / params$k_weighted
  fold_enrichment <- observed_fraction / expected_fraction

  # Prepare diagnostics
  diagnostics <- list(
    original_query_genes = original_query_size,
    original_background_genes = original_background_size,
    filtered_query_genes = nrow(query_df),
    filtered_background_genes = nrow(background_df),
    genes_removed = (original_query_size + original_background_size) -
                   (nrow(query_df) + nrow(background_df)),
    pathway_genes_in_query = sum(query_df$gene %in% pathway_genes),
    pathway_genes_in_background = sum(background_df$gene %in% pathway_genes),
    mean_query_copies = mean(query_df$copy_number),
    mean_background_copies = mean(background_df$copy_number),
    copy_model_consistency = validate_copy_model_consistency(query_df, background_df)
  )

  # Return comprehensive results
  return(list(
    pvalue = pvalue,
    overlap_instances = params$q_weighted,
    query_instances = params$k_weighted,
    pathway_instances = params$m_weighted,
    background_instances = total_weighted,
    fold_enrichment = fold_enrichment,
    parameters = c(
      q = params$q_weighted,
      m = params$m_weighted,
      n = params$n_weighted,
      k = params$k_weighted
    ),
    diagnostics = diagnostics,
    method = "copy_weighted_hypergeometric"
  ))
}

#' Calculate Weighted Hypergeometric Parameters
#'
#' Internal function to compute copy-number weighted parameters for phyper()
#'
#' @param query_df Data frame with gene and copy_number columns
#' @param pathway_genes Vector of pathway gene symbols
#' @param background_df Data frame with gene and copy_number columns
#'
#' @return List with q_weighted, m_weighted, n_weighted, k_weighted
#'
#' @noRd
calculate_weighted_parameters <- function(query_df, pathway_genes, background_df) {

  # k_weighted: Total query instances
  k_weighted <- sum(query_df$copy_number)

  # q_weighted: Overlap instances (query genes in pathway)
  query_in_pathway <- query_df[query_df$gene %in% pathway_genes, ]
  q_weighted <- sum(query_in_pathway$copy_number)

  # m_weighted: Pathway instances in background
  pathway_in_background <- background_df[background_df$gene %in% pathway_genes, ]
  m_weighted <- sum(pathway_in_background$copy_number)

  # n_weighted: Non-pathway instances in background
  n_weighted <- sum(background_df$copy_number) - m_weighted

  return(list(
    q_weighted = q_weighted,
    m_weighted = m_weighted,
    n_weighted = n_weighted,
    k_weighted = k_weighted
  ))
}

#' Comprehensive Input Validation
#'
#' Validates inputs for weighted hypergeometric testing
#'
#' @param query_df Query data frame
#' @param pathway_genes Pathway gene vector
#' @param background_df Background data frame
#' @param max_copies Maximum allowed copy number
#' @param min_instances Minimum required instances
#'
#' @return List with valid (logical) and message (character)
#'
#' @noRd
validate_hypergeometric_inputs <- function(query_df, pathway_genes, background_df,
                                          max_copies, min_instances) {

  # Check basic structure
  if (!is.data.frame(query_df)) {
    return(list(valid = FALSE, message = "query_df must be a data frame"))
  }

  if (!is.data.frame(background_df)) {
    return(list(valid = FALSE, message = "background_df must be a data frame"))
  }

  if (!is.vector(pathway_genes)) {
    return(list(valid = FALSE, message = "pathway_genes must be a vector"))
  }

  # Check required columns
  required_cols <- c("gene", "copy_number")

  if (!all(required_cols %in% names(query_df))) {
    missing <- setdiff(required_cols, names(query_df))
    return(list(valid = FALSE,
                message = paste("query_df missing columns:", paste(missing, collapse = ", "))))
  }

  if (!all(required_cols %in% names(background_df))) {
    missing <- setdiff(required_cols, names(background_df))
    return(list(valid = FALSE,
                message = paste("background_df missing columns:", paste(missing, collapse = ", "))))
  }

  # Check data types
  if (!is.numeric(query_df$copy_number)) {
    return(list(valid = FALSE, message = "query_df$copy_number must be numeric"))
  }

  if (!is.numeric(background_df$copy_number)) {
    return(list(valid = FALSE, message = "background_df$copy_number must be numeric"))
  }

  # Check for valid copy numbers
  if (any(query_df$copy_number < 0)) {
    return(list(valid = FALSE, message = "query_df contains negative copy numbers"))
  }

  if (any(background_df$copy_number < 0)) {
    return(list(valid = FALSE, message = "background_df contains negative copy numbers"))
  }

  if (any(!is.finite(query_df$copy_number))) {
    return(list(valid = FALSE, message = "query_df contains non-finite copy numbers"))
  }

  if (any(!is.finite(background_df$copy_number))) {
    return(list(valid = FALSE, message = "background_df contains non-finite copy numbers"))
  }

  # Check for duplicated genes
  if (any(duplicated(query_df$gene))) {
    return(list(valid = FALSE, message = "query_df contains duplicated genes"))
  }

  if (any(duplicated(background_df$gene))) {
    return(list(valid = FALSE, message = "background_df contains duplicated genes"))
  }

  # Check empty inputs
  if (nrow(query_df) == 0) {
    return(list(valid = FALSE, message = "query_df is empty"))
  }

  if (nrow(background_df) == 0) {
    return(list(valid = FALSE, message = "background_df is empty"))
  }

  if (length(pathway_genes) == 0) {
    return(list(valid = FALSE, message = "pathway_genes is empty"))
  }

  # Check reasonable parameter ranges
  if (!is.null(max_copies) && max_copies <= 0) {
    return(list(valid = FALSE, message = "max_copies must be positive"))
  }

  if (!is.null(min_instances) && min_instances < 1) {
    return(list(valid = FALSE, message = "min_instances must be at least 1"))
  }

  # Check if pathway has any genes in background
  pathway_in_bg <- sum(background_df$gene %in% pathway_genes)
  if (pathway_in_bg == 0) {
    return(list(valid = FALSE,
                message = "No pathway genes found in background dataset"))
  }

  # Check if query has sufficient genes
  if (nrow(query_df) < 2) {
    return(list(valid = FALSE,
                message = "Query dataset too small (need at least 2 genes)"))
  }

  return(list(valid = TRUE, message = "All inputs valid"))
}

#' Validate Hypergeometric Distribution Constraints
#'
#' Checks that weighted parameters satisfy hypergeometric constraints
#'
#' @param q_weighted Observed overlap instances
#' @param m_weighted Pathway instances in background
#' @param n_weighted Non-pathway instances in background
#' @param k_weighted Query instances
#'
#' @return List with valid (logical) and message (character)
#'
#' @noRd
validate_hypergeometric_constraints <- function(q_weighted, m_weighted, n_weighted, k_weighted) {

  # Check non-negative integers
  params <- c(q_weighted, m_weighted, n_weighted, k_weighted)
  if (any(params < 0)) {
    return(list(valid = FALSE, message = "Parameters must be non-negative"))
  }

  if (any(params != floor(params))) {
    return(list(valid = FALSE, message = "Parameters must be integers"))
  }

  # Check logical constraints
  if (q_weighted > k_weighted) {
    return(list(valid = FALSE,
                message = paste("Overlap instances (", q_weighted,
                              ") cannot exceed query instances (", k_weighted, ")")))
  }

  if (q_weighted > m_weighted) {
    return(list(valid = FALSE,
                message = paste("Overlap instances (", q_weighted,
                              ") cannot exceed pathway instances (", m_weighted, ")")))
  }

  total_background <- m_weighted + n_weighted
  if (k_weighted > total_background) {
    return(list(valid = FALSE,
                message = paste("Query instances (", k_weighted,
                              ") cannot exceed total background instances (",
                              total_background, ")")))
  }

  # Check feasibility constraints
  min_possible_overlap <- max(0, k_weighted - n_weighted)
  max_possible_overlap <- min(k_weighted, m_weighted)

  if (q_weighted < min_possible_overlap) {
    return(list(valid = FALSE,
                message = paste("Overlap too small: minimum possible is",
                              min_possible_overlap)))
  }

  if (q_weighted > max_possible_overlap) {
    return(list(valid = FALSE,
                message = paste("Overlap too large: maximum possible is",
                              max_possible_overlap)))
  }

  return(list(valid = TRUE, message = "All constraints satisfied"))
}

#' Filter Zero-Copy Genes
#'
#' Removes genes with zero copy numbers and issues warning
#'
#' @param df Data frame with gene and copy_number columns
#'
#' @return Filtered data frame
#'
#' @noRd
filter_zero_copies <- function(df) {
  zero_copy_genes <- df[df$copy_number == 0, "gene"]
  if (length(zero_copy_genes) > 0) {
    warning(paste("Removing", length(zero_copy_genes),
                  "genes with zero copies:",
                  paste(head(zero_copy_genes, 3), collapse = ", "),
                  if (length(zero_copy_genes) > 3) "..." else ""))
    df <- df[df$copy_number > 0, ]
  }
  return(df)
}

#' Cap Extreme Copy Numbers
#'
#' Limits copy numbers to prevent numerical overflow
#'
#' @param df Data frame with gene and copy_number columns
#' @param max_copies Maximum allowed copy number
#'
#' @return Data frame with capped copy numbers
#'
#' @noRd
cap_extreme_copies <- function(df, max_copies) {
  extreme_indices <- df$copy_number > max_copies
  extreme_genes <- df[extreme_indices, ]

  if (nrow(extreme_genes) > 0) {
    warning(paste("Capping", nrow(extreme_genes), "genes at", max_copies,
                  "copies. Genes:",
                  paste(head(extreme_genes$gene, 3), collapse = ", "),
                  if (nrow(extreme_genes) > 3) "..." else ""))
    df$copy_number[extreme_indices] <- max_copies
  }

  return(df)
}

#' Validate Copy Model Consistency
#'
#' Checks correlation between query and background copy numbers for common genes
#'
#' @param query_df Query data frame
#' @param background_df Background data frame
#'
#' @return List with correlation coefficient and warning flags
#'
#' @noRd
validate_copy_model_consistency <- function(query_df, background_df) {

  common_genes <- intersect(query_df$gene, background_df$gene)

  if (length(common_genes) < 5) {
    return(list(
      correlation = NA,
      n_common_genes = length(common_genes),
      warning = "Too few common genes for correlation analysis"
    ))
  }

  # Merge copy numbers for common genes
  query_subset <- query_df[query_df$gene %in% common_genes, ]
  bg_subset <- background_df[background_df$gene %in% common_genes, ]

  merged <- merge(query_subset, bg_subset, by = "gene",
                  suffixes = c("_query", "_background"))

  # Calculate correlation
  correlation <- cor(merged$copy_number_query, merged$copy_number_background)

  # Check for warning conditions
  warning_msg <- NULL
  if (correlation < 0.5) {
    warning_msg <- "Low correlation between query and background copy numbers"
  }

  # Check for extreme differences
  max_query <- max(query_df$copy_number)
  max_bg <- max(background_df$copy_number)

  if (max_query > 10 * max_bg || max_bg > 10 * max_query) {
    extreme_warning <- "Extreme copy number differences between query and background"
    warning_msg <- if (is.null(warning_msg)) extreme_warning else
                   paste(warning_msg, extreme_warning, sep = "; ")
  }

  if (!is.null(warning_msg)) {
    warning(warning_msg)
  }

  return(list(
    correlation = correlation,
    n_common_genes = length(common_genes),
    warning = warning_msg
  ))
}

#' Compare Weighted vs Standard Hypergeometric Results
#'
#' Utility function to compare copy-weighted and standard hypergeometric results
#'
#' @param query_df Query data frame with gene and copy_number columns
#' @param pathway_genes Vector of pathway gene symbols
#' @param background_df Background data frame with gene and copy_number columns
#'
#' @return List comparing both methods
#'
#' @examples
#' # Create example data
#' query <- data.frame(gene = c("A", "B", "C"), copy_number = c(1, 5, 3))
#' pathway <- c("A", "B")
#' background <- data.frame(gene = LETTERS[1:10], copy_number = rep(2, 10))
#'
#' comparison <- compare_weighted_vs_standard(query, pathway, background)
#' print(comparison)
#'
#' @export
compare_weighted_vs_standard <- function(query_df, pathway_genes, background_df) {

  # Standard hypergeometric test
  query_genes <- unique(query_df$gene)
  background_genes <- unique(background_df$gene)

  k_standard <- length(query_genes)
  q_standard <- length(intersect(query_genes, pathway_genes))
  m_standard <- length(intersect(background_genes, pathway_genes))
  n_standard <- length(background_genes) - m_standard

  pval_standard <- phyper(q_standard - 1, m_standard, n_standard, k_standard,
                         lower.tail = FALSE)

  # Weighted hypergeometric test
  weighted_result <- weighted_hypergeometric_test(query_df, pathway_genes,
                                                 background_df, validate_inputs = TRUE)

  # Calculate scaling factors
  scaling_factors <- list(
    query_scaling = weighted_result$query_instances / k_standard,
    pathway_scaling = weighted_result$pathway_instances / m_standard,
    background_scaling = weighted_result$background_instances / length(background_genes),
    overlap_scaling = weighted_result$overlap_instances / q_standard
  )

  return(list(
    standard = list(
      pvalue = pval_standard,
      parameters = c(q = q_standard, m = m_standard, n = n_standard, k = k_standard)
    ),
    weighted = list(
      pvalue = weighted_result$pvalue,
      parameters = weighted_result$parameters
    ),
    comparison = list(
      pvalue_ratio = weighted_result$pvalue / pval_standard,
      log10_pvalue_diff = log10(weighted_result$pvalue) - log10(pval_standard),
      scaling_factors = scaling_factors,
      more_significant = weighted_result$pvalue < pval_standard
    )
  ))
}