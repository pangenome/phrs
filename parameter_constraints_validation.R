# Parameter Constraints Validation for Copy-Number Weighted phyper()
# Implement validation functions for parameter constraints and edge cases

library(tidyverse)
library(yaml)

# Load parameter bounds configuration
load_parameter_bounds <- function(config_file = "copy_weighted_ora_parameter_bounds.yaml") {
  if (file.exists(config_file)) {
    yaml::read_yaml(config_file)
  } else {
    # Default configuration if file not found
    list(
      safety_bounds = list(
        max_total_copies = 1e8,
        max_single_copy = 1e6,
        max_copy_ratio = 1e4,
        min_pathway_size = 3,
        min_query_size = 3,
        min_expected_overlap = 1e-10
      ),
      validation_criteria = list(
        hypergeometric = list(
          non_negativity = TRUE,
          integer_values = TRUE,
          overlap_bounds = TRUE,
          sample_bounds = TRUE,
          feasibility = TRUE
        )
      )
    )
  }
}

# Core parameter constraint validation
validate_hypergeometric_parameters <- function(q_weighted, m_weighted, n_weighted, k_weighted) {
  # Validate hypergeometric parameters for copy-number weighted phyper().
  #
  # Parameters:
  #   q_weighted: overlap instances (gene copies in both query and pathway)
  #   m_weighted: pathway instances in background (total pathway gene copies)
  #   n_weighted: non-pathway instances in background
  #   k_weighted: query instances (total query gene copies)
  #
  # Returns:
  #   List with validation results and detailed constraint checks

  validation_results <- list(
    parameters = list(q = q_weighted, m = m_weighted, n = n_weighted, k = k_weighted),
    constraints_satisfied = TRUE,
    constraint_violations = character(0),
    warnings = character(0),
    edge_cases = character(0)
  )

  # Constraint 1: Non-negativity
  if (any(c(q_weighted, m_weighted, n_weighted, k_weighted) < 0)) {
    validation_results$constraints_satisfied <- FALSE
    validation_results$constraint_violations <- c(
      validation_results$constraint_violations,
      "VIOLATION: All parameters must be non-negative"
    )
  }

  # Special handling for all-zeros edge case
  if (all(c(q_weighted, m_weighted, n_weighted, k_weighted) == 0)) {
    validation_results$edge_cases <- c(validation_results$edge_cases, "all_zeros")
    validation_results$warnings <- c(validation_results$warnings, "All parameters are zero - degenerate case")
    # This is still mathematically valid (0 <= 0 <= min(0,0) is true)
    return(validation_results)
  }

  # Constraint 2: Integer values
  if (!all(c(q_weighted, m_weighted, n_weighted, k_weighted) ==
           as.integer(c(q_weighted, m_weighted, n_weighted, k_weighted)))) {
    validation_results$constraints_satisfied <- FALSE
    validation_results$constraint_violations <- c(
      validation_results$constraint_violations,
      "VIOLATION: All parameters must be integers"
    )
  }

  # Constraint 3: Overlap bounds - q ≤ min(m, k)
  if (q_weighted > min(m_weighted, k_weighted)) {
    validation_results$constraints_satisfied <- FALSE
    validation_results$constraint_violations <- c(
      validation_results$constraint_violations,
      sprintf("VIOLATION: Overlap constraint - q (%d) > min(m=%d, k=%d)",
              q_weighted, m_weighted, k_weighted)
    )
  }

  # Constraint 4: Sample bounds - k ≤ m + n
  if (k_weighted > (m_weighted + n_weighted)) {
    validation_results$constraints_satisfied <- FALSE
    validation_results$constraint_violations <- c(
      validation_results$constraint_violations,
      sprintf("VIOLATION: Sample constraint - k (%d) > m+n (%d)",
              k_weighted, m_weighted + n_weighted)
    )
  }

  # Constraint 5: Feasibility bounds - max(0, k-n) ≤ q ≤ min(k, m)
  min_feasible_q <- max(0, k_weighted - n_weighted)
  max_feasible_q <- min(k_weighted, m_weighted)

  if (q_weighted < min_feasible_q || q_weighted > max_feasible_q) {
    validation_results$constraints_satisfied <- FALSE
    validation_results$constraint_violations <- c(
      validation_results$constraint_violations,
      sprintf("VIOLATION: Feasibility constraint - q (%d) not in [%d, %d]",
              q_weighted, min_feasible_q, max_feasible_q)
    )
  }

  # Edge case detection
  if (q_weighted == 0) {
    validation_results$edge_cases <- c(validation_results$edge_cases, "zero_overlap")
  }

  if (q_weighted == k_weighted && k_weighted > 0) {
    validation_results$edge_cases <- c(validation_results$edge_cases, "complete_query_overlap")
  }

  if (q_weighted == m_weighted && m_weighted > 0) {
    validation_results$edge_cases <- c(validation_results$edge_cases, "complete_pathway_overlap")
  }

  if (k_weighted == m_weighted + n_weighted) {
    validation_results$edge_cases <- c(validation_results$edge_cases, "query_equals_population")
  }

  # Safety bound warnings
  config <- load_parameter_bounds()
  total_instances <- m_weighted + n_weighted

  if (total_instances > config$safety_bounds$max_total_copies) {
    validation_results$warnings <- c(
      validation_results$warnings,
      sprintf("WARNING: Total instances (%d) exceeds safety bound (%d)",
              total_instances, config$safety_bounds$max_total_copies)
    )
  }

  # Check for extreme parameter ratios
  if (max(c(q_weighted, m_weighted, n_weighted, k_weighted)) > 0) {
    max_ratio <- max(c(q_weighted, m_weighted, n_weighted, k_weighted)) /
                 max(1, min(c(q_weighted, m_weighted, n_weighted, k_weighted)))

    if (max_ratio > config$safety_bounds$max_copy_ratio) {
      validation_results$warnings <- c(
        validation_results$warnings,
        sprintf("WARNING: Extreme parameter ratio (%.0f) may cause numerical instability", max_ratio)
      )
    }
  }

  return(validation_results)
}

# Gene-level constraint validation
validate_gene_level_constraints <- function(query_df, pathway_genes, background_df) {
  # Validate constraints at the gene level before parameter calculation.
  #
  # Parameters:
  #   query_df: data frame with gene and copy_number columns
  #   pathway_genes: vector of pathway gene identifiers
  #   background_df: data frame with gene and copy_number columns for background
  #
  # Returns:
  #   List with validation results for gene-level constraints

  validation_results <- list(
    valid = TRUE,
    violations = character(0),
    warnings = character(0),
    corrections_applied = character(0)
  )

  config <- load_parameter_bounds()

  # Check input format
  required_cols <- c("gene", "copy_number")
  if (!all(required_cols %in% names(query_df))) {
    validation_results$valid <- FALSE
    validation_results$violations <- c(
      validation_results$violations,
      "Query dataframe missing required columns: gene, copy_number"
    )
  }

  if (!all(required_cols %in% names(background_df))) {
    validation_results$valid <- FALSE
    validation_results$violations <- c(
      validation_results$violations,
      "Background dataframe missing required columns: gene, copy_number"
    )
  }

  if (!validation_results$valid) {
    return(validation_results)
  }

  # Check copy number constraints
  if (any(query_df$copy_number < 0, na.rm = TRUE)) {
    validation_results$valid <- FALSE
    validation_results$violations <- c(
      validation_results$violations,
      "Negative copy numbers found in query"
    )
  }

  if (any(background_df$copy_number < 0, na.rm = TRUE)) {
    validation_results$valid <- FALSE
    validation_results$violations <- c(
      validation_results$violations,
      "Negative copy numbers found in background"
    )
  }

  # Check for extreme copy numbers
  max_query_copy <- max(query_df$copy_number, na.rm = TRUE)
  max_bg_copy <- max(background_df$copy_number, na.rm = TRUE)

  if (max_query_copy > config$safety_bounds$max_single_copy) {
    validation_results$warnings <- c(
      validation_results$warnings,
      sprintf("Very large copy numbers in query (max: %d)", max_query_copy)
    )
  }

  if (max_bg_copy > config$safety_bounds$max_single_copy) {
    validation_results$warnings <- c(
      validation_results$warnings,
      sprintf("Very large copy numbers in background (max: %d)", max_bg_copy)
    )
  }

  # Check gene ID uniqueness
  if (any(duplicated(query_df$gene))) {
    validation_results$valid <- FALSE
    validation_results$violations <- c(
      validation_results$violations,
      "Duplicate gene identifiers in query"
    )
  }

  if (any(duplicated(background_df$gene))) {
    validation_results$valid <- FALSE
    validation_results$violations <- c(
      validation_results$violations,
      "Duplicate gene identifiers in background"
    )
  }

  # Check query-background relationship
  genes_not_in_bg <- setdiff(query_df$gene, background_df$gene)
  if (length(genes_not_in_bg) > 0) {
    validation_results$warnings <- c(
      validation_results$warnings,
      sprintf("%d query genes not found in background (will be filtered)",
              length(genes_not_in_bg))
    )
  }

  # Check pathway coverage
  pathway_in_bg <- sum(pathway_genes %in% background_df$gene)
  pathway_coverage <- pathway_in_bg / length(pathway_genes)

  if (pathway_coverage < config$safety_bounds$min_background_coverage) {
    validation_results$warnings <- c(
      validation_results$warnings,
      sprintf("Low pathway coverage in background (%.1f%%)", pathway_coverage * 100)
    )
  }

  # Check minimum sizes
  if (length(pathway_genes) < config$safety_bounds$min_pathway_size) {
    validation_results$warnings <- c(
      validation_results$warnings,
      sprintf("Small pathway size (%d genes)", length(pathway_genes))
    )
  }

  if (nrow(query_df) < config$safety_bounds$min_query_size) {
    validation_results$warnings <- c(
      validation_results$warnings,
      sprintf("Small query size (%d genes)", nrow(query_df))
    )
  }

  return(validation_results)
}

# Edge case handler
handle_edge_cases <- function(validation_result, q_weighted, m_weighted, n_weighted, k_weighted) {
  # Handle special edge cases for copy-number weighted hypergeometric test.
  #
  # Returns:
  #   List with appropriate handling for detected edge cases

  edge_case_results <- list(
    requires_special_handling = FALSE,
    special_result = NULL,
    recommendations = character(0)
  )

  config <- load_parameter_bounds()

  for (edge_case in validation_result$edge_cases) {
    edge_case_results$requires_special_handling <- TRUE

    switch(edge_case,
      "zero_overlap" = {
        edge_case_results$special_result <- list(
          pvalue = 1.0,
          fold_enrichment = 0.0,
          overlap_instances = 0,
          query_instances = k_weighted,
          pathway_instances = m_weighted,
          background_instances = m_weighted + n_weighted,
          note = "No overlap between query and pathway"
        )
        edge_case_results$recommendations <- c(
          edge_case_results$recommendations,
          "Consider different pathway or expand gene set"
        )
      },

      "complete_query_overlap" = {
        # Calculate p-value for complete overlap case
        pvalue <- phyper(q_weighted-1, m_weighted, n_weighted, k_weighted, lower.tail = FALSE)
        fold_enrich <- (q_weighted/k_weighted) / (m_weighted/(m_weighted+n_weighted))

        edge_case_results$special_result <- list(
          pvalue = pvalue,
          fold_enrichment = fold_enrich,
          overlap_instances = q_weighted,
          query_instances = k_weighted,
          pathway_instances = m_weighted,
          background_instances = m_weighted + n_weighted,
          note = "Complete overlap - all query genes in pathway"
        )
        edge_case_results$recommendations <- c(
          edge_case_results$recommendations,
          "Verify biological plausibility of complete overlap"
        )
      },

      "complete_pathway_overlap" = {
        # All pathway genes present in query
        pvalue <- phyper(q_weighted-1, m_weighted, n_weighted, k_weighted, lower.tail = FALSE)
        fold_enrich <- (q_weighted/k_weighted) / (m_weighted/(m_weighted+n_weighted))

        edge_case_results$special_result <- list(
          pvalue = pvalue,
          fold_enrichment = fold_enrich,
          overlap_instances = q_weighted,
          query_instances = k_weighted,
          pathway_instances = m_weighted,
          background_instances = m_weighted + n_weighted,
          note = "Complete pathway capture - all pathway instances in query"
        )
        edge_case_results$recommendations <- c(
          edge_case_results$recommendations,
          "Consider larger pathway or more stringent query selection"
        )
      },

      "query_equals_population" = {
        edge_case_results$special_result <- list(
          pvalue = 1.0,
          fold_enrichment = 1.0,
          overlap_instances = q_weighted,
          query_instances = k_weighted,
          pathway_instances = m_weighted,
          background_instances = m_weighted + n_weighted,
          note = "Query contains entire population"
        )
        edge_case_results$recommendations <- c(
          edge_case_results$recommendations,
          "Use different statistical approach - hypergeometric test not appropriate"
        )
      }
    )
  }

  return(edge_case_results)
}

# Constraint violation detection and error handling
detect_constraint_violations <- function(query_df, pathway_genes, background_df) {
  # Comprehensive constraint violation detection with detailed error reporting.
  #
  # Returns:
  #   List with violation details and suggested corrections

  # First validate at gene level
  gene_validation <- validate_gene_level_constraints(query_df, pathway_genes, background_df)

  if (!gene_validation$valid) {
    return(list(
      valid = FALSE,
      level = "gene_level",
      violations = gene_validation$violations,
      warnings = gene_validation$warnings
    ))
  }

  # Calculate parameters with validated inputs
  # Filter query to background genes first
  valid_query_genes <- intersect(query_df$gene, background_df$gene)
  query_filtered <- query_df[query_df$gene %in% valid_query_genes, ]

  if (nrow(query_filtered) == 0) {
    return(list(
      valid = FALSE,
      level = "parameter_level",
      violations = "No query genes found in background after filtering",
      warnings = character(0)
    ))
  }

  # Use background copy numbers for consistency
  query_with_bg_copies <- merge(query_filtered[, "gene", drop = FALSE],
                               background_df, by = "gene", all.x = TRUE)

  # Calculate weighted parameters
  query_in_pathway <- query_with_bg_copies[query_with_bg_copies$gene %in% pathway_genes, ]
  pathway_in_background <- background_df[background_df$gene %in% pathway_genes, ]

  k_weighted <- sum(query_with_bg_copies$copy_number)        # query instances
  q_weighted <- sum(query_in_pathway$copy_number)            # overlap instances
  m_weighted <- sum(pathway_in_background$copy_number)       # pathway instances
  n_weighted <- sum(background_df$copy_number) - m_weighted  # non-pathway instances

  # Validate hypergeometric parameters
  param_validation <- validate_hypergeometric_parameters(q_weighted, m_weighted,
                                                        n_weighted, k_weighted)

  # Combine results
  violation_results <- list(
    valid = gene_validation$valid && param_validation$constraints_satisfied,
    level = ifelse(param_validation$constraints_satisfied, "validated", "parameter_level"),
    gene_violations = gene_validation$violations,
    parameter_violations = param_validation$constraint_violations,
    combined_warnings = c(gene_validation$warnings, param_validation$warnings),
    edge_cases = param_validation$edge_cases,
    parameters = param_validation$parameters
  )

  return(violation_results)
}

# Parameter transformation consistency validation
validate_parameter_transformation_consistency <- function(query_df, pathway_genes, background_df) {
  # Validate that parameter transformations maintain mathematical consistency.
  #
  # Compares weighted parameter approach with instance expansion method.

  # Method 1: Weighted parameter calculation
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

  # Method 2: Instance expansion for comparison
  query_expanded <- rep(query_with_bg_copies$gene, query_with_bg_copies$copy_number)
  background_expanded <- rep(background_df$gene, background_df$copy_number)

  q_expanded <- sum(query_expanded %in% pathway_genes)
  m_expanded <- sum(background_expanded %in% pathway_genes)
  n_expanded <- length(background_expanded) - m_expanded
  k_expanded <- length(query_expanded)

  # Compare parameters
  consistency_results <- list(
    method_1_params = list(q = q_weighted, m = m_weighted, n = n_weighted, k = k_weighted),
    method_2_params = list(q = q_expanded, m = m_expanded, n = n_expanded, k = k_expanded),
    parameters_match = all(c(
      q_weighted == q_expanded,
      m_weighted == m_expanded,
      n_weighted == n_expanded,
      k_weighted == k_expanded
    )),
    consistency_validated = TRUE
  )

  # If parameters don't match, identify discrepancy
  if (!consistency_results$parameters_match) {
    consistency_results$consistency_validated <- FALSE
    consistency_results$discrepancies <- list(
      q_match = q_weighted == q_expanded,
      m_match = m_weighted == m_expanded,
      n_match = n_weighted == n_expanded,
      k_match = k_weighted == k_expanded
    )
  }

  return(consistency_results)
}

# Comprehensive constraint validation test suite
run_constraint_validation_tests <- function() {
  # Run comprehensive test suite for parameter constraint validation.

  cat("=== Parameter Constraint Validation Test Suite ===\n")

  test_results <- list()

  # Test 1: Basic constraint validation
  cat("\nTest 1: Basic hypergeometric constraints\n")
  cat("=========================================\n")

  # Valid parameters
  valid_test <- validate_hypergeometric_parameters(5, 20, 80, 30)
  cat("Valid parameters test:", ifelse(valid_test$constraints_satisfied, "PASS", "FAIL"), "\n")

  # Invalid parameters (q > min(k, m))
  invalid_test1 <- validate_hypergeometric_parameters(25, 20, 80, 30)
  cat("Invalid overlap constraint:", ifelse(!invalid_test1$constraints_satisfied, "PASS", "FAIL"), "\n")

  # Invalid parameters (k > m + n)
  invalid_test2 <- validate_hypergeometric_parameters(5, 20, 30, 60)
  cat("Invalid sample constraint:", ifelse(!invalid_test2$constraints_satisfied, "PASS", "FAIL"), "\n")

  test_results$basic_constraints <- list(
    valid_params = valid_test$constraints_satisfied,
    invalid_overlap = !invalid_test1$constraints_satisfied,
    invalid_sample = !invalid_test2$constraints_satisfied
  )

  # Test 2: Edge case detection
  cat("\nTest 2: Edge case detection\n")
  cat("===========================\n")

  # Zero overlap
  zero_overlap_test <- validate_hypergeometric_parameters(0, 20, 80, 30)
  has_zero_overlap_edge <- "zero_overlap" %in% zero_overlap_test$edge_cases
  cat("Zero overlap detection:", ifelse(has_zero_overlap_edge, "PASS", "FAIL"), "\n")

  # Complete overlap
  complete_overlap_test <- validate_hypergeometric_parameters(20, 20, 80, 20)
  has_complete_overlap_edge <- "complete_query_overlap" %in% complete_overlap_test$edge_cases
  cat("Complete overlap detection:", ifelse(has_complete_overlap_edge, "PASS", "FAIL"), "\n")

  test_results$edge_case_detection <- list(
    zero_overlap = has_zero_overlap_edge,
    complete_overlap = has_complete_overlap_edge
  )

  # Test 3: Gene-level validation
  cat("\nTest 3: Gene-level constraint validation\n")
  cat("=======================================\n")

  # Create test data
  background_df <- data.frame(
    gene = paste0("GENE", 1:100),
    copy_number = sample(1:10, 100, replace = TRUE),
    stringsAsFactors = FALSE
  )

  pathway_genes <- sample(background_df$gene, 20)
  query_df <- background_df[sample(nrow(background_df), 15), ]

  gene_validation <- validate_gene_level_constraints(query_df, pathway_genes, background_df)
  cat("Gene-level validation:", ifelse(gene_validation$valid, "PASS", "FAIL"), "\n")

  test_results$gene_level_validation <- gene_validation$valid

  # Test 4: Parameter transformation consistency
  cat("\nTest 4: Parameter transformation consistency\n")
  cat("===========================================\n")

  consistency_test <- validate_parameter_transformation_consistency(query_df, pathway_genes, background_df)
  cat("Transformation consistency:", ifelse(consistency_test$consistency_validated, "PASS", "FAIL"), "\n")

  test_results$transformation_consistency <- consistency_test$consistency_validated

  # Test 5: Comprehensive constraint violation detection
  cat("\nTest 5: Constraint violation detection\n")
  cat("=====================================\n")

  violation_test <- detect_constraint_violations(query_df, pathway_genes, background_df)
  cat("Violation detection:", ifelse(violation_test$valid, "PASS", "FAIL"), "\n")

  test_results$violation_detection <- violation_test$valid

  # Overall test summary
  all_tests_pass <- all(unlist(test_results))
  cat("\n=== TEST SUMMARY ===\n")
  cat("Overall result:", ifelse(all_tests_pass, "ALL TESTS PASS", "SOME TESTS FAILED"), "\n")

  if (!all_tests_pass) {
    cat("Failed tests:\n")
    failed_tests <- names(test_results)[!unlist(test_results)]
    for (test_name in failed_tests) {
      cat("  -", test_name, "\n")
    }
  }

  return(test_results)
}

# Edge case test suite with extreme scenarios
test_extreme_edge_cases <- function() {
  # Test constraint validation with extreme edge cases and boundary conditions.

  cat("=== Extreme Edge Case Testing ===\n")

  edge_case_results <- list()

  # Test zero copy numbers
  cat("\n1. Zero copy number scenarios\n")
  cat("=============================\n")

  # All zero copies
  zero_test <- validate_hypergeometric_parameters(0, 0, 0, 0)
  has_zero_edge_case <- "all_zeros" %in% zero_test$edge_cases
  cat("All zeros:", ifelse(has_zero_edge_case || zero_test$constraints_satisfied, "HANDLED", "ERROR"), "\n")
  edge_case_results$all_zeros <- has_zero_edge_case || zero_test$constraints_satisfied

  # Test extreme copy number ratios
  cat("\n2. Extreme copy number ratios\n")
  cat("=============================\n")

  extreme_ratio_test <- validate_hypergeometric_parameters(1, 1000000, 100, 50)
  has_ratio_warning <- length(extreme_ratio_test$warnings) > 0
  cat("Extreme ratios warning:", ifelse(has_ratio_warning, "DETECTED", "MISSED"), "\n")
  edge_case_results$extreme_ratios <- has_ratio_warning

  # Test boundary values
  cat("\n3. Boundary value testing\n")
  cat("=========================\n")

  # Exact boundary cases
  boundary_test1 <- validate_hypergeometric_parameters(10, 10, 90, 10)  # q = min(k, m)
  boundary_test2 <- validate_hypergeometric_parameters(50, 50, 50, 100) # k = m + n, q at min feasible

  cat("Boundary case 1 (q = min(k,m)):", ifelse(boundary_test1$constraints_satisfied, "PASS", "FAIL"), "\n")
  cat("Boundary case 2 (k = m+n):", ifelse(boundary_test2$constraints_satisfied, "PASS", "FAIL"), "\n")

  edge_case_results$boundary_cases <- boundary_test1$constraints_satisfied && boundary_test2$constraints_satisfied

  # Test negative values
  cat("\n4. Negative value handling\n")
  cat("==========================\n")

  negative_test <- validate_hypergeometric_parameters(-1, 20, 80, 30)
  cat("Negative values:", ifelse(!negative_test$constraints_satisfied, "REJECTED", "ERROR"), "\n")
  edge_case_results$negative_handling <- !negative_test$constraints_satisfied

  # Test non-integer values
  cat("\n5. Non-integer value handling\n")
  cat("=============================\n")

  non_integer_test <- validate_hypergeometric_parameters(5.5, 20.3, 80.7, 30.1)
  cat("Non-integer values:", ifelse(!non_integer_test$constraints_satisfied, "REJECTED", "ERROR"), "\n")
  edge_case_results$non_integer_handling <- !non_integer_test$constraints_satisfied

  # Overall edge case handling
  all_edge_cases_handled <- all(unlist(edge_case_results))
  cat("\n=== EDGE CASE SUMMARY ===\n")
  cat("All edge cases handled:", ifelse(all_edge_cases_handled, "YES", "NO"), "\n")

  return(edge_case_results)
}

# Main execution
if (!interactive()) {
  cat("PARAMETER CONSTRAINTS VALIDATION\n")
  cat("================================\n")

  # Run test suites
  constraint_tests <- run_constraint_validation_tests()
  edge_case_tests <- test_extreme_edge_cases()

  # Save results
  validation_report <- list(
    timestamp = Sys.time(),
    constraint_validation_tests = constraint_tests,
    edge_case_tests = edge_case_tests,
    overall_status = ifelse(all(unlist(constraint_tests)) && all(unlist(edge_case_tests)),
                           "ALL_TESTS_PASS", "SOME_TESTS_FAILED")
  )

  save(validation_report, file = "parameter_constraints_validation_report.RData")
  cat("\nValidation report saved to: parameter_constraints_validation_report.RData\n")
}