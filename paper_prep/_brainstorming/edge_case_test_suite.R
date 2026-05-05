# Comprehensive Edge Case Test Suite for Copy-Number Weighted phyper() Parameters
# Extended testing for boundary conditions and extreme scenarios

library(tidyverse)
source("parameter_constraints_validation.R")

# Generate edge case datasets
generate_edge_case_datasets <- function() {
  edge_case_data <- list()

  # 1. Zero copy number scenarios
  edge_case_data$zero_copies <- list(
    all_zero = data.frame(gene = paste0("G", 1:10), copy_number = rep(0, 10)),
    mixed_zero = data.frame(gene = paste0("G", 1:10), copy_number = c(rep(0, 5), rep(2, 5))),
    single_zero = data.frame(gene = paste0("G", 1:10), copy_number = c(0, rep(1, 9)))
  )

  # 2. Extreme copy number ratios
  edge_case_data$extreme_ratios <- list(
    high_variance = data.frame(gene = paste0("G", 1:10), copy_number = c(1, 1000000, rep(1, 8))),
    geometric_progression = data.frame(gene = paste0("G", 1:8), copy_number = 2^(0:7)),
    extreme_skew = data.frame(gene = paste0("G", 1:10), copy_number = c(rep(1, 9), 10000))
  )

  # 3. Large-scale datasets
  edge_case_data$large_scale <- list(
    large_uniform = data.frame(gene = paste0("G", 1:100000), copy_number = rep(5, 100000)),
    large_variable = data.frame(gene = paste0("G", 1:50000),
                               copy_number = sample(1:1000, 50000, replace = TRUE))
  )

  # 4. Boundary condition datasets
  edge_case_data$boundary_conditions <- list(
    single_gene_query = data.frame(gene = "G1", copy_number = 10),
    single_gene_pathway = "G1",
    identical_query_pathway = data.frame(gene = paste0("G", 1:5), copy_number = rep(3, 5))
  )

  return(edge_case_data)
}

# Test zero copy number handling
test_zero_copy_scenarios <- function() {
  cat("=== Zero Copy Number Scenario Testing ===\n")

  edge_data <- generate_edge_case_datasets()
  zero_test_results <- list()

  # Test 1: All zero copies
  cat("\n1. All zero copies scenario\n")
  background_zero <- edge_data$zero_copies$all_zero
  pathway_genes <- c("G1", "G2", "G3")
  query_zero <- background_zero[1:3, ]

  tryCatch({
    violation_result <- detect_constraint_violations(query_zero, pathway_genes, background_zero)
    zero_test_results$all_zero <- list(
      valid = violation_result$valid,
      handled_appropriately = !violation_result$valid || length(violation_result$edge_cases) > 0
    )
    cat("   Result: ", ifelse(zero_test_results$all_zero$handled_appropriately, "HANDLED", "ERROR"), "\n")
  }, error = function(e) {
    cat("   Error:", e$message, "\n")
    zero_test_results$all_zero <- list(valid = FALSE, error = e$message)
  })

  # Test 2: Mixed zero/non-zero copies
  cat("\n2. Mixed zero/non-zero scenario\n")
  background_mixed <- edge_data$zero_copies$mixed_zero
  query_mixed <- background_mixed[c(1, 6:8), ]  # Include both zero and non-zero genes

  tryCatch({
    violation_result <- detect_constraint_violations(query_mixed, pathway_genes, background_mixed)
    zero_test_results$mixed_zero <- list(
      valid = violation_result$valid,
      parameters_calculated = !is.null(violation_result$parameters)
    )
    cat("   Result: ", ifelse(violation_result$valid, "VALID", "INVALID"), "\n")
  }, error = function(e) {
    cat("   Error:", e$message, "\n")
    zero_test_results$mixed_zero <- list(valid = FALSE, error = e$message)
  })

  return(zero_test_results)
}

# Test extreme ratio scenarios
test_extreme_ratio_scenarios <- function() {
  cat("\n=== Extreme Copy Number Ratio Testing ===\n")

  edge_data <- generate_edge_case_datasets()
  ratio_test_results <- list()

  # Test high variance scenario
  cat("\n1. High variance copy numbers\n")
  background_extreme <- edge_data$extreme_ratios$high_variance
  pathway_genes <- c("G1", "G2", "G3")
  query_extreme <- background_extreme[1:4, ]

  tryCatch({
    violation_result <- detect_constraint_violations(query_extreme, pathway_genes, background_extreme)
    param_result <- validate_hypergeometric_parameters(
      violation_result$parameters$q,
      violation_result$parameters$m,
      violation_result$parameters$n,
      violation_result$parameters$k
    )

    ratio_test_results$high_variance <- list(
      valid = violation_result$valid,
      has_warnings = length(param_result$warnings) > 0,
      max_ratio = max(background_extreme$copy_number) / min(background_extreme$copy_number[background_extreme$copy_number > 0])
    )
    cat("   Max ratio:", round(ratio_test_results$high_variance$max_ratio, 0), "\n")
    cat("   Warnings generated:", ifelse(ratio_test_results$high_variance$has_warnings, "YES", "NO"), "\n")
  }, error = function(e) {
    cat("   Error:", e$message, "\n")
    ratio_test_results$high_variance <- list(valid = FALSE, error = e$message)
  })

  # Test geometric progression
  cat("\n2. Geometric progression copy numbers\n")
  background_geom <- edge_data$extreme_ratios$geometric_progression
  query_geom <- background_geom[1:4, ]
  pathway_genes_geom <- c("G1", "G2", "G8")  # Include smallest and largest

  tryCatch({
    violation_result <- detect_constraint_violations(query_geom, pathway_genes_geom, background_geom)
    ratio_test_results$geometric <- list(
      valid = violation_result$valid,
      parameter_range = max(c(violation_result$parameters$q, violation_result$parameters$m,
                            violation_result$parameters$n, violation_result$parameters$k))
    )
    cat("   Max parameter value:", ratio_test_results$geometric$parameter_range, "\n")
  }, error = function(e) {
    cat("   Error:", e$message, "\n")
    ratio_test_results$geometric <- list(valid = FALSE, error = e$message)
  })

  return(ratio_test_results)
}

# Test boundary conditions comprehensively
test_comprehensive_boundary_conditions <- function() {
  cat("\n=== Comprehensive Boundary Condition Testing ===\n")

  boundary_test_results <- list()

  # Generate test scenarios systematically
  test_scenarios <- list(
    # Format: c(q, m, n, k, description)
    c(0, 10, 90, 50, "Zero overlap"),
    c(10, 10, 90, 10, "Complete query in pathway"),
    c(10, 10, 90, 100, "Query equals population"),
    c(50, 50, 50, 50, "All parameters equal"),
    c(1, 1000, 1000, 1, "Minimal query, large background"),
    c(999, 1000, 1000, 1000, "Near-complete overlap"),
    c(100, 100, 0, 100, "No non-pathway genes"),
    c(0, 0, 100, 50, "No pathway genes"),
    c(25, 50, 50, 50, "Exact feasibility boundary")
  )

  for (i in seq_along(test_scenarios)) {
    scenario <- test_scenarios[[i]]
    q <- as.integer(scenario[1])
    m <- as.integer(scenario[2])
    n <- as.integer(scenario[3])
    k <- as.integer(scenario[4])
    desc <- scenario[5]

    cat(sprintf("\n%d. %s (q=%d, m=%d, n=%d, k=%d)\n", i, desc, q, m, n, k))

    tryCatch({
      param_result <- validate_hypergeometric_parameters(q, m, n, k)
      edge_case_result <- handle_edge_cases(param_result, q, m, n, k)

      boundary_test_results[[paste0("scenario_", i)]] <- list(
        description = desc,
        parameters = c(q=q, m=m, n=n, k=k),
        constraints_satisfied = param_result$constraints_satisfied,
        edge_cases = param_result$edge_cases,
        requires_special_handling = edge_case_result$requires_special_handling,
        violations = param_result$constraint_violations
      )

      if (param_result$constraints_satisfied) {
        cat("   Status: VALID")
        if (length(param_result$edge_cases) > 0) {
          cat(" (Edge cases:", paste(param_result$edge_cases, collapse=", "), ")")
        }
        cat("\n")
      } else {
        cat("   Status: INVALID -", paste(param_result$constraint_violations, collapse="; "), "\n")
      }

    }, error = function(e) {
      cat("   Error:", e$message, "\n")
      boundary_test_results[[paste0("scenario_", i)]] <- list(
        description = desc,
        error = e$message
      )
    })
  }

  return(boundary_test_results)
}

# Test parameter consistency across transformations
test_parameter_transformation_robustness <- function(n_tests = 50) {
  cat("\n=== Parameter Transformation Robustness Testing ===\n")

  transformation_results <- list()
  consistency_failures <- 0

  for (i in 1:n_tests) {
    # Generate random valid scenario
    n_genes <- sample(50:500, 1)
    background_df <- data.frame(
      gene = paste0("G", 1:n_genes),
      copy_number = sample(1:20, n_genes, replace = TRUE),
      stringsAsFactors = FALSE
    )

    pathway_size <- sample(5:min(50, n_genes/4), 1)
    pathway_genes <- sample(background_df$gene, pathway_size)

    query_size <- sample(5:min(100, n_genes/2), 1)
    query_genes <- sample(background_df$gene, query_size)
    query_df <- background_df[background_df$gene %in% query_genes, ]

    tryCatch({
      consistency_result <- validate_parameter_transformation_consistency(query_df, pathway_genes, background_df)

      if (!consistency_result$consistency_validated) {
        consistency_failures <- consistency_failures + 1
        cat("   Test", i, ": INCONSISTENT\n")
        if (i <= 5) {  # Print details for first few failures
          cat("     Method 1:", paste(names(consistency_result$method_1_params), "=",
                                      consistency_result$method_1_params, collapse=", "), "\n")
          cat("     Method 2:", paste(names(consistency_result$method_2_params), "=",
                                      consistency_result$method_2_params, collapse=", "), "\n")
        }
      }
    }, error = function(e) {
      consistency_failures <- consistency_failures + 1
      if (i <= 5) cat("   Test", i, ": ERROR -", e$message, "\n")
    })
  }

  consistency_rate <- (n_tests - consistency_failures) / n_tests
  cat(sprintf("\nTransformation consistency rate: %.1f%% (%d/%d tests passed)\n",
              consistency_rate * 100, n_tests - consistency_failures, n_tests))

  transformation_results$consistency_rate <- consistency_rate
  transformation_results$failures <- consistency_failures
  transformation_results$total_tests <- n_tests

  return(transformation_results)
}

# Performance testing with large datasets
test_performance_scalability <- function() {
  cat("\n=== Performance and Scalability Testing ===\n")

  performance_results <- list()

  # Test different dataset sizes
  test_sizes <- c(100, 1000, 10000, 50000)

  for (size in test_sizes) {
    cat(sprintf("\nTesting with %d genes...\n", size))

    # Generate large dataset
    background_df <- data.frame(
      gene = paste0("G", 1:size),
      copy_number = sample(1:10, size, replace = TRUE),
      stringsAsFactors = FALSE
    )

    pathway_genes <- sample(background_df$gene, min(size/10, 1000))
    query_genes <- sample(background_df$gene, min(size/5, 2000))
    query_df <- background_df[background_df$gene %in% query_genes, ]

    # Time the validation
    start_time <- Sys.time()

    tryCatch({
      violation_result <- detect_constraint_violations(query_df, pathway_genes, background_df)

      end_time <- Sys.time()
      runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

      performance_results[[paste0("size_", size)]] <- list(
        dataset_size = size,
        runtime_seconds = runtime,
        validation_successful = violation_result$valid,
        memory_efficient = runtime < 30  # Arbitrary threshold
      )

      cat(sprintf("   Runtime: %.2f seconds\n", runtime))
      cat(sprintf("   Validation: %s\n", ifelse(violation_result$valid, "PASS", "FAIL")))

    }, error = function(e) {
      cat("   Error:", e$message, "\n")
      performance_results[[paste0("size_", size)]] <- list(
        dataset_size = size,
        error = e$message
      )
    })
  }

  return(performance_results)
}

# Main comprehensive edge case test runner
run_comprehensive_edge_case_tests <- function() {
  cat("========================================\n")
  cat("COMPREHENSIVE EDGE CASE TEST SUITE\n")
  cat("========================================\n")

  start_time <- Sys.time()

  # Run all test suites
  test_results <- list(
    metadata = list(
      start_time = start_time,
      r_version = R.version.string
    ),
    zero_copy_tests = test_zero_copy_scenarios(),
    extreme_ratio_tests = test_extreme_ratio_scenarios(),
    boundary_condition_tests = test_comprehensive_boundary_conditions(),
    transformation_robustness = test_parameter_transformation_robustness(n_tests = 100),
    performance_tests = test_performance_scalability()
  )

  test_results$metadata$end_time <- Sys.time()
  test_results$metadata$total_runtime <- test_results$metadata$end_time - start_time

  # Generate summary report
  cat("\n========================================\n")
  cat("EDGE CASE TEST SUMMARY\n")
  cat("========================================\n")

  total_boundary_tests <- length(test_results$boundary_condition_tests)
  valid_boundary_tests <- sum(sapply(test_results$boundary_condition_tests, function(x) {
    if ("constraints_satisfied" %in% names(x)) x$constraints_satisfied else FALSE
  }))

  cat(sprintf("Zero copy scenarios: %d tests completed\n",
              length(test_results$zero_copy_tests)))
  cat(sprintf("Extreme ratio scenarios: %d tests completed\n",
              length(test_results$extreme_ratio_tests)))
  cat(sprintf("Boundary conditions: %d/%d tests valid\n",
              valid_boundary_tests, total_boundary_tests))
  cat(sprintf("Transformation consistency: %.1f%% pass rate\n",
              test_results$transformation_robustness$consistency_rate * 100))
  cat(sprintf("Performance tests: %d size scenarios completed\n",
              length(test_results$performance_tests)))

  cat(sprintf("\nTotal runtime: %.1f minutes\n",
              as.numeric(test_results$metadata$total_runtime, units = "mins")))

  # Save results
  save(test_results, file = "comprehensive_edge_case_test_results.RData")
  cat("\nResults saved to: comprehensive_edge_case_test_results.RData\n")

  return(test_results)
}

# Execute comprehensive tests
if (!interactive()) {
  results <- run_comprehensive_edge_case_tests()
}