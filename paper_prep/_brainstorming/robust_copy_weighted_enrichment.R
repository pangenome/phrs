# Robust Copy-Number Weighted Enrichment Analysis
# Production-ready implementation with comprehensive error handling and validation

#' Copy-Number Weighted Hypergeometric Enrichment Test
#'
#' Performs pathway enrichment analysis incorporating gene copy number weights.
#' This function implements the mathematical formulation from copy_number_weighted_phyper_mathematical_formulation.md
#'
#' @param query_df Data frame with columns 'gene' and 'copy_number'
#' @param pathway_genes Character vector of gene IDs in pathway
#' @param background_df Data frame with columns 'gene' and 'copy_number'
#' @param min_overlap Minimum overlap instances for testing (default: 1)
#' @param validate_inputs Perform input validation (default: TRUE)
#'
#' @return List containing test results and diagnostics
#'
#' @examples
#' query_df <- data.frame(gene = c("A", "B", "C"), copy_number = c(2, 1, 3))
#' pathway_genes <- c("A", "C", "D")
#' background_df <- data.frame(gene = LETTERS[1:10], copy_number = rep(1, 10))
#'
#' result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)
#' print(result$p_value)
weighted_hypergeometric_test <- function(query_df,
                                       pathway_genes,
                                       background_df,
                                       min_overlap = 1,
                                       validate_inputs = TRUE) {

  # Performance timing start
  start_time <- Sys.time()

  # Input validation with comprehensive error handling
  if (validate_inputs) {
    validation_result <- validate_enrichment_inputs(query_df, pathway_genes, background_df)
    if (!validation_result$is_valid) {
      stop("Input validation failed: ", validation_result$error_message)
    }
  }

  tryCatch({
    # Calculate weighted parameters using optimized approach
    params <- calculate_weighted_parameters(query_df, pathway_genes, background_df)

    # Extract parameters
    k_weighted <- params$k_weighted
    q_weighted <- params$q_weighted
    m_weighted <- params$m_weighted
    n_weighted <- params$n_weighted

    # Check minimum overlap requirement
    if (q_weighted < min_overlap) {
      return(list(
        p_value = 1.0,
        fold_enrichment = 0,
        overlap_instances = q_weighted,
        query_instances = k_weighted,
        pathway_instances = m_weighted,
        background_instances = m_weighted + n_weighted,
        status = "insufficient_overlap",
        computation_time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      ))
    }

    # Hypergeometric parameter validation
    param_validation <- validate_weighted_params(q_weighted, m_weighted, n_weighted, k_weighted)
    if (!param_validation$is_valid) {
      stop("Parameter validation failed: ", param_validation$error_message)
    }

    # Perform hypergeometric test with error handling
    p_value <- safe_phyper(q_weighted - 1, m_weighted, n_weighted, k_weighted)

    # Calculate enrichment metrics
    expected_overlap <- (k_weighted * m_weighted) / (m_weighted + n_weighted)
    fold_enrichment <- ifelse(expected_overlap > 0,
                             q_weighted / expected_overlap,
                             Inf)

    # Calculate additional diagnostic metrics
    depletion_p_value <- phyper(q_weighted, m_weighted, n_weighted, k_weighted, lower.tail = TRUE)

    # Calculate confidence intervals for fold enrichment (approximate)
    ci <- calculate_enrichment_ci(q_weighted, k_weighted, m_weighted, n_weighted)

    # Computation time
    end_time <- Sys.time()
    computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Return comprehensive results
    list(
      # Core results
      p_value = p_value,
      fold_enrichment = fold_enrichment,

      # Counts
      overlap_instances = q_weighted,
      query_instances = k_weighted,
      pathway_instances = m_weighted,
      background_instances = m_weighted + n_weighted,

      # Additional metrics
      expected_overlap = expected_overlap,
      depletion_p_value = depletion_p_value,
      enrichment_ci_lower = ci$lower,
      enrichment_ci_upper = ci$upper,

      # Diagnostics
      status = "success",
      computation_time = computation_time,
      parameters_used = list(q = q_weighted, m = m_weighted, n = n_weighted, k = k_weighted)
    )

  }, error = function(e) {
    # Graceful error handling
    list(
      p_value = NA,
      fold_enrichment = NA,
      status = "error",
      error_message = as.character(e),
      computation_time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    )
  })
}

#' Validate inputs for enrichment analysis
#'
#' @param query_df Query gene data frame
#' @param pathway_genes Pathway gene vector
#' @param background_df Background gene data frame
#'
#' @return List with validation result
validate_enrichment_inputs <- function(query_df, pathway_genes, background_df) {

  # Check data frame structures
  if (!is.data.frame(query_df)) {
    return(list(is_valid = FALSE, error_message = "query_df must be a data frame"))
  }

  if (!is.data.frame(background_df)) {
    return(list(is_valid = FALSE, error_message = "background_df must be a data frame"))
  }

  # Check required columns
  required_cols <- c("gene", "copy_number")

  if (!all(required_cols %in% names(query_df))) {
    missing <- setdiff(required_cols, names(query_df))
    return(list(is_valid = FALSE,
                error_message = paste("query_df missing columns:", paste(missing, collapse = ", "))))
  }

  if (!all(required_cols %in% names(background_df))) {
    missing <- setdiff(required_cols, names(background_df))
    return(list(is_valid = FALSE,
                error_message = paste("background_df missing columns:", paste(missing, collapse = ", "))))
  }

  # Check data types
  if (!is.character(query_df$gene) && !is.factor(query_df$gene)) {
    return(list(is_valid = FALSE, error_message = "query_df$gene must be character or factor"))
  }

  if (!is.numeric(query_df$copy_number)) {
    return(list(is_valid = FALSE, error_message = "query_df$copy_number must be numeric"))
  }

  if (!is.character(background_df$gene) && !is.factor(background_df$gene)) {
    return(list(is_valid = FALSE, error_message = "background_df$gene must be character or factor"))
  }

  if (!is.numeric(background_df$copy_number)) {
    return(list(is_valid = FALSE, error_message = "background_df$copy_number must be numeric"))
  }

  # Check pathway genes
  if (!is.character(pathway_genes) && !is.factor(pathway_genes)) {
    return(list(is_valid = FALSE, error_message = "pathway_genes must be character or factor"))
  }

  if (length(pathway_genes) == 0) {
    return(list(is_valid = FALSE, error_message = "pathway_genes cannot be empty"))
  }

  # Check for missing values
  if (any(is.na(query_df$gene))) {
    return(list(is_valid = FALSE, error_message = "query_df$gene contains NA values"))
  }

  if (any(is.na(query_df$copy_number))) {
    return(list(is_valid = FALSE, error_message = "query_df$copy_number contains NA values"))
  }

  if (any(is.na(background_df$gene))) {
    return(list(is_valid = FALSE, error_message = "background_df$gene contains NA values"))
  }

  if (any(is.na(background_df$copy_number))) {
    return(list(is_valid = FALSE, error_message = "background_df$copy_number contains NA values"))
  }

  if (any(is.na(pathway_genes))) {
    return(list(is_valid = FALSE, error_message = "pathway_genes contains NA values"))
  }

  # Check copy number constraints
  if (any(query_df$copy_number < 0)) {
    return(list(is_valid = FALSE, error_message = "query_df$copy_number cannot be negative"))
  }

  if (any(background_df$copy_number < 0)) {
    return(list(is_valid = FALSE, error_message = "background_df$copy_number cannot be negative"))
  }

  if (any(query_df$copy_number != floor(query_df$copy_number))) {
    return(list(is_valid = FALSE, error_message = "query_df$copy_number must be integers"))
  }

  if (any(background_df$copy_number != floor(background_df$copy_number))) {
    return(list(is_valid = FALSE, error_message = "background_df$copy_number must be integers"))
  }

  # Check for empty data
  if (nrow(query_df) == 0) {
    return(list(is_valid = FALSE, error_message = "query_df cannot be empty"))
  }

  if (nrow(background_df) == 0) {
    return(list(is_valid = FALSE, error_message = "background_df cannot be empty"))
  }

  # Check that query genes are in background
  query_not_in_bg <- setdiff(query_df$gene, background_df$gene)
  if (length(query_not_in_bg) > 0) {
    return(list(is_valid = FALSE,
                error_message = paste("Query genes not in background:",
                                     paste(head(query_not_in_bg, 5), collapse = ", "),
                                     ifelse(length(query_not_in_bg) > 5, "...", ""))))
  }

  # Check for duplicate genes
  if (any(duplicated(query_df$gene))) {
    return(list(is_valid = FALSE, error_message = "query_df contains duplicate genes"))
  }

  if (any(duplicated(background_df$gene))) {
    return(list(is_valid = FALSE, error_message = "background_df contains duplicate genes"))
  }

  return(list(is_valid = TRUE, error_message = NULL))
}

#' Calculate weighted hypergeometric parameters
#'
#' @param query_df Query gene data frame
#' @param pathway_genes Pathway gene vector
#' @param background_df Background gene data frame
#'
#' @return List with weighted parameters
calculate_weighted_parameters <- function(query_df, pathway_genes, background_df) {

  # Convert to character for consistent matching
  query_genes <- as.character(query_df$gene)
  pathway_genes <- as.character(pathway_genes)
  background_genes <- as.character(background_df$gene)

  # Calculate k_weighted (total query instances)
  k_weighted <- sum(query_df$copy_number)

  # Calculate q_weighted (overlap instances)
  query_in_pathway <- query_df[query_genes %in% pathway_genes, , drop = FALSE]
  q_weighted <- sum(query_in_pathway$copy_number)

  # Calculate m_weighted (pathway instances in background)
  background_in_pathway <- background_df[background_genes %in% pathway_genes, , drop = FALSE]
  m_weighted <- sum(background_in_pathway$copy_number)

  # Calculate n_weighted (non-pathway instances in background)
  n_weighted <- sum(background_df$copy_number) - m_weighted

  list(
    k_weighted = k_weighted,
    q_weighted = q_weighted,
    m_weighted = m_weighted,
    n_weighted = n_weighted,

    # Additional diagnostics
    query_genes_in_pathway = nrow(query_in_pathway),
    background_genes_in_pathway = nrow(background_in_pathway),
    total_background_genes = nrow(background_df)
  )
}

#' Validate weighted hypergeometric parameters
#'
#' @param q_w Weighted overlap count
#' @param m_w Weighted pathway count
#' @param n_w Weighted non-pathway count
#' @param k_w Weighted query count
#'
#' @return List with validation result
validate_weighted_params <- function(q_w, m_w, n_w, k_w) {

  # Check for numeric input
  params <- c(q_w, m_w, n_w, k_w)
  param_names <- c("q_weighted", "m_weighted", "n_weighted", "k_weighted")

  if (!all(is.numeric(params))) {
    return(list(is_valid = FALSE, error_message = "All parameters must be numeric"))
  }

  # Check for missing values
  if (any(is.na(params))) {
    na_params <- param_names[is.na(params)]
    return(list(is_valid = FALSE,
                error_message = paste("Parameters contain NA:", paste(na_params, collapse = ", "))))
  }

  # Check for infinite values
  if (any(is.infinite(params))) {
    inf_params <- param_names[is.infinite(params)]
    return(list(is_valid = FALSE,
                error_message = paste("Parameters contain infinite values:", paste(inf_params, collapse = ", "))))
  }

  # Check non-negativity
  if (any(params < 0)) {
    neg_params <- param_names[params < 0]
    return(list(is_valid = FALSE,
                error_message = paste("Parameters must be non-negative:", paste(neg_params, collapse = ", "))))
  }

  # Check integer constraint
  if (any(params != floor(params))) {
    non_int_params <- param_names[params != floor(params)]
    return(list(is_valid = FALSE,
                error_message = paste("Parameters must be integers:", paste(non_int_params, collapse = ", "))))
  }

  # Check hypergeometric constraints
  if (q_w > k_w) {
    return(list(is_valid = FALSE,
                error_message = paste("Overlap instances (", q_w, ") exceeds query instances (", k_w, ")")))
  }

  if (q_w > m_w) {
    return(list(is_valid = FALSE,
                error_message = paste("Overlap instances (", q_w, ") exceeds pathway instances (", m_w, ")")))
  }

  if (k_w > (m_w + n_w)) {
    return(list(is_valid = FALSE,
                error_message = paste("Query instances (", k_w, ") exceeds total background instances (", m_w + n_w, ")")))
  }

  # Check feasibility bounds
  min_possible <- max(0, k_w - n_w)
  if (q_w < min_possible) {
    return(list(is_valid = FALSE,
                error_message = paste("Overlap instances (", q_w, ") below feasible minimum (", min_possible, ")")))
  }

  return(list(is_valid = TRUE, error_message = NULL))
}

#' Safe phyper computation with overflow protection
#'
#' @param q Overlap parameter (already decremented)
#' @param m Pathway instances
#' @param n Non-pathway instances
#' @param k Query instances
#'
#' @return p-value or appropriate fallback
safe_phyper <- function(q, m, n, k) {

  # Check for edge cases
  if (q < 0) {
    return(1.0)  # No overlap
  }

  if (m == 0 || k == 0) {
    return(1.0)  # Empty pathway or query
  }

  # Check for extreme parameters that might cause numerical issues
  total_pop <- m + n
  if (total_pop > 1e8) {
    warning("Very large population size may cause numerical instability")
  }

  # Perform calculation with error handling
  tryCatch({
    p_val <- phyper(q, m, n, k, lower.tail = FALSE)

    # Check for numerical issues
    if (is.na(p_val) || is.infinite(p_val)) {
      warning("phyper returned NA or Inf, using fallback approximation")
      return(calculate_normal_approximation(q + 1, m, n, k))
    }

    # Ensure valid probability range
    p_val <- max(0, min(1, p_val))

    return(p_val)

  }, error = function(e) {
    warning("phyper failed: ", e$message, ", using normal approximation")
    return(calculate_normal_approximation(q + 1, m, n, k))
  })
}

#' Normal approximation for hypergeometric when phyper fails
#'
#' @param overlap_obs Observed overlap
#' @param m Pathway instances
#' @param n Non-pathway instances
#' @param k Query instances
#'
#' @return Approximate p-value
calculate_normal_approximation <- function(overlap_obs, m, n, k) {

  N <- m + n
  p <- m / N

  # Hypergeometric mean and variance
  mu <- k * p
  var_hyper <- k * p * (1 - p) * (N - k) / (N - 1)

  if (var_hyper <= 0) {
    return(ifelse(overlap_obs > mu, 0, 1))
  }

  # Continuity correction
  z_score <- (overlap_obs - 0.5 - mu) / sqrt(var_hyper)

  # Return upper tail probability
  p_val <- pnorm(z_score, lower.tail = FALSE)
  return(max(0, min(1, p_val)))
}

#' Calculate approximate confidence interval for fold enrichment
#'
#' @param overlap Observed overlap instances
#' @param query Total query instances
#' @param pathway Pathway instances in background
#' @param non_pathway Non-pathway instances in background
#' @param confidence Confidence level (default 0.95)
#'
#' @return List with lower and upper bounds
calculate_enrichment_ci <- function(overlap, query, pathway, non_pathway, confidence = 0.95) {

  alpha <- 1 - confidence
  z <- qnorm(1 - alpha/2)

  total_bg <- pathway + non_pathway
  expected <- (query * pathway) / total_bg

  if (expected == 0 || overlap == 0) {
    return(list(lower = 0, upper = Inf))
  }

  # Approximate standard error for log fold enrichment
  se_log_fe <- sqrt(1/overlap + 1/expected)

  log_fe <- log(overlap / expected)

  ci_lower <- exp(log_fe - z * se_log_fe)
  ci_upper <- exp(log_fe + z * se_log_fe)

  list(lower = ci_lower, upper = ci_upper)
}

# Performance benchmarking function
#' Benchmark copy-weighted enrichment vs instance expansion
#'
#' @param query_df Query gene data frame
#' @param pathway_genes Pathway genes
#' @param background_df Background gene data frame
#' @param n_reps Number of benchmark repetitions (default 100)
#'
#' @return Benchmark results
benchmark_enrichment_methods <- function(query_df, pathway_genes, background_df, n_reps = 100) {

  # Parameter weighting method timing
  param_times <- replicate(n_reps, {
    start <- Sys.time()
    result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df, validate_inputs = FALSE)
    end <- Sys.time()
    as.numeric(difftime(end, start, units = "secs"))
  })

  # Instance expansion method timing
  expansion_times <- replicate(n_reps, {
    start <- Sys.time()

    # Expand datasets
    query_expanded <- rep(query_df$gene, query_df$copy_number)
    background_expanded <- rep(background_df$gene, background_df$copy_number)

    # Standard phyper
    q_exp <- sum(query_expanded %in% pathway_genes)
    m_exp <- sum(background_expanded %in% pathway_genes)
    n_exp <- length(background_expanded) - m_exp
    k_exp <- length(query_expanded)

    p_val <- phyper(q_exp - 1, m_exp, n_exp, k_exp, lower.tail = FALSE)

    end <- Sys.time()
    as.numeric(difftime(end, start, units = "secs"))
  })

  # Memory usage estimates
  query_instances <- sum(query_df$copy_number)
  background_instances <- sum(background_df$copy_number)
  param_memory <- nrow(query_df) + nrow(background_df)  # Approximate
  expansion_memory <- query_instances + background_instances

  list(
    parameter_method = list(
      mean_time = mean(param_times),
      median_time = median(param_times),
      sd_time = sd(param_times),
      memory_estimate = param_memory
    ),
    expansion_method = list(
      mean_time = mean(expansion_times),
      median_time = median(expansion_times),
      sd_time = sd(expansion_times),
      memory_estimate = expansion_memory
    ),
    speedup_factor = mean(expansion_times) / mean(param_times),
    memory_reduction = 1 - (param_memory / expansion_memory),
    n_repetitions = n_reps
  )
}