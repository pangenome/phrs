# Comprehensive Validation Tests for Copy-Number Weighted Enrichment
# Production-ready testing framework with statistical validation

source("robust_copy_weighted_enrichment.R")

#' Comprehensive test suite for weighted enrichment analysis
#'
#' @param test_data_dir Directory containing test datasets (optional)
#' @param verbose Print detailed test results
#'
#' @return Test results summary
run_enrichment_validation_suite <- function(test_data_dir = NULL, verbose = TRUE) {

  if (verbose) cat("Running comprehensive enrichment validation suite...\n")

  results <- list()

  # Test 1: Basic functionality
  if (verbose) cat("Test 1: Basic functionality... ")
  results$basic_functionality <- test_basic_functionality()
  if (verbose) cat(ifelse(results$basic_functionality$passed, "PASS\n", "FAIL\n"))

  # Test 2: Input validation
  if (verbose) cat("Test 2: Input validation... ")
  results$input_validation <- test_input_validation()
  if (verbose) cat(ifelse(results$input_validation$passed, "PASS\n", "FAIL\n"))

  # Test 3: Edge cases
  if (verbose) cat("Test 3: Edge cases... ")
  results$edge_cases <- test_edge_cases()
  if (verbose) cat(ifelse(results$edge_cases$passed, "PASS\n", "FAIL\n"))

  # Test 4: Mathematical equivalence
  if (verbose) cat("Test 4: Mathematical equivalence... ")
  results$mathematical_equivalence <- test_mathematical_equivalence()
  if (verbose) cat(ifelse(results$mathematical_equivalence$passed, "PASS\n", "FAIL\n"))

  # Test 5: Null distribution
  if (verbose) cat("Test 5: Null distribution uniformity... ")
  results$null_distribution <- test_null_distribution()
  if (verbose) cat(ifelse(results$null_distribution$passed, "PASS\n", "FAIL\n"))

  # Test 6: Type I error control
  if (verbose) cat("Test 6: Type I error control... ")
  results$type_i_error <- test_type_i_error_control()
  if (verbose) cat(ifelse(results$type_i_error$passed, "PASS\n", "FAIL\n"))

  # Test 7: Performance benchmarks
  if (verbose) cat("Test 7: Performance benchmarks... ")
  results$performance <- test_performance_benchmarks()
  if (verbose) cat(ifelse(results$performance$passed, "PASS\n", "FAIL\n"))

  # Overall summary
  all_passed <- all(sapply(results, function(x) x$passed))
  results$overall_passed <- all_passed

  if (verbose) {
    cat("\n=== VALIDATION SUMMARY ===\n")
    cat("Overall result:", ifelse(all_passed, "PASS", "FAIL"), "\n")
    cat("Individual tests:\n")
    for (test_name in names(results)) {
      if (test_name != "overall_passed") {
        cat(sprintf("  %s: %s\n", test_name, ifelse(results[[test_name]]$passed, "PASS", "FAIL")))
      }
    }
  }

  return(results)
}

#' Test basic functionality
test_basic_functionality <- function() {

  tryCatch({
    # Create simple test data
    query_df <- data.frame(
      gene = c("GENE1", "GENE2", "GENE3"),
      copy_number = c(2, 1, 3)
    )

    pathway_genes <- c("GENE1", "GENE3", "GENE5")

    background_df <- data.frame(
      gene = paste0("GENE", 1:10),
      copy_number = c(2, 1, 3, 1, 1, 2, 1, 1, 1, 2)
    )

    # Run test
    result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)

    # Basic checks
    checks <- c(
      "has_p_value" = !is.na(result$p_value),
      "p_value_range" = result$p_value >= 0 && result$p_value <= 1,
      "has_fold_enrichment" = !is.na(result$fold_enrichment),
      "has_counts" = all(c("overlap_instances", "query_instances") %in% names(result)),
      "status_success" = result$status == "success"
    )

    list(
      passed = all(checks),
      details = checks,
      result = result
    )

  }, error = function(e) {
    list(passed = FALSE, error = as.character(e))
  })
}

#' Test input validation
test_input_validation <- function() {

  passed_tests <- 0
  total_tests <- 0

  # Valid baseline data
  valid_query <- data.frame(gene = c("A", "B"), copy_number = c(1, 2))
  valid_pathway <- c("A", "C")
  valid_background <- data.frame(gene = c("A", "B", "C"), copy_number = c(1, 2, 1))

  # Test 1: Missing columns
  total_tests <- total_tests + 1
  tryCatch({
    bad_query <- data.frame(gene = c("A", "B"))  # Missing copy_number
    weighted_hypergeometric_test(bad_query, valid_pathway, valid_background)
    # Should not reach here
  }, error = function(e) {
    if (grepl("missing columns", e$message)) passed_tests <- passed_tests + 1
  })

  # Test 2: Wrong data types
  total_tests <- total_tests + 1
  tryCatch({
    bad_query <- data.frame(gene = c("A", "B"), copy_number = c("1", "2"))  # String copy numbers
    weighted_hypergeometric_test(bad_query, valid_pathway, valid_background)
  }, error = function(e) {
    if (grepl("must be numeric", e$message)) passed_tests <- passed_tests + 1
  })

  # Test 3: Negative copy numbers
  total_tests <- total_tests + 1
  tryCatch({
    bad_query <- data.frame(gene = c("A", "B"), copy_number = c(1, -1))
    weighted_hypergeometric_test(bad_query, valid_pathway, valid_background)
  }, error = function(e) {
    if (grepl("cannot be negative", e$message)) passed_tests <- passed_tests + 1
  })

  # Test 4: Non-integer copy numbers
  total_tests <- total_tests + 1
  tryCatch({
    bad_query <- data.frame(gene = c("A", "B"), copy_number = c(1.5, 2))
    weighted_hypergeometric_test(bad_query, valid_pathway, valid_background)
  }, error = function(e) {
    if (grepl("must be integers", e$message)) passed_tests <- passed_tests + 1
  })

  # Test 5: Query genes not in background
  total_tests <- total_tests + 1
  tryCatch({
    bad_query <- data.frame(gene = c("A", "X"), copy_number = c(1, 2))  # X not in background
    weighted_hypergeometric_test(bad_query, valid_pathway, valid_background)
  }, error = function(e) {
    if (grepl("not in background", e$message)) passed_tests <- passed_tests + 1
  })

  # Test 6: Empty pathway
  total_tests <- total_tests + 1
  tryCatch({
    weighted_hypergeometric_test(valid_query, character(0), valid_background)
  }, error = function(e) {
    if (grepl("cannot be empty", e$message)) passed_tests <- passed_tests + 1
  })

  # Test 7: Duplicate genes
  total_tests <- total_tests + 1
  tryCatch({
    bad_query <- data.frame(gene = c("A", "A"), copy_number = c(1, 2))
    weighted_hypergeometric_test(bad_query, valid_pathway, valid_background)
  }, error = function(e) {
    if (grepl("duplicate genes", e$message)) passed_tests <- passed_tests + 1
  })

  list(
    passed = passed_tests == total_tests,
    passed_count = passed_tests,
    total_count = total_tests
  )
}

#' Test edge cases
test_edge_cases <- function() {

  tests_passed <- 0
  total_tests <- 0

  # Test 1: No overlap
  total_tests <- total_tests + 1
  tryCatch({
    query_df <- data.frame(gene = c("A", "B"), copy_number = c(1, 2))
    pathway_genes <- c("C", "D")
    background_df <- data.frame(gene = c("A", "B", "C", "D"), copy_number = c(1, 2, 1, 1))

    result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)

    if (result$overlap_instances == 0 && result$p_value == 1.0) {
      tests_passed <- tests_passed + 1
    }
  }, error = function(e) {
    # Expected to handle gracefully
  })

  # Test 2: Perfect overlap
  total_tests <- total_tests + 1
  tryCatch({
    query_df <- data.frame(gene = c("A", "B"), copy_number = c(1, 2))
    pathway_genes <- c("A", "B", "C")
    background_df <- data.frame(gene = c("A", "B", "C"), copy_number = c(1, 2, 1))

    result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)

    if (result$overlap_instances == result$query_instances && result$p_value >= 0) {
      tests_passed <- tests_passed + 1
    }
  }, error = function(e) {
    # Should handle this case
  })

  # Test 3: Single gene query
  total_tests <- total_tests + 1
  tryCatch({
    query_df <- data.frame(gene = "A", copy_number = 1)
    pathway_genes <- c("A", "B")
    background_df <- data.frame(gene = c("A", "B", "C"), copy_number = c(1, 1, 1))

    result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)

    if (!is.na(result$p_value) && result$status == "success") {
      tests_passed <- tests_passed + 1
    }
  }, error = function(e) {
    # Should handle this case
  })

  # Test 4: All genes in pathway
  total_tests <- total_tests + 1
  tryCatch({
    query_df <- data.frame(gene = c("A", "B"), copy_number = c(1, 2))
    pathway_genes <- c("A", "B")
    background_df <- data.frame(gene = c("A", "B"), copy_number = c(1, 2))

    result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)

    if (!is.na(result$p_value) && result$status == "success") {
      tests_passed <- tests_passed + 1
    }
  }, error = function(e) {
    # Should handle this case
  })

  # Test 5: Very large copy numbers
  total_tests <- total_tests + 1
  tryCatch({
    query_df <- data.frame(gene = c("A", "B"), copy_number = c(1000, 2000))
    pathway_genes <- c("A", "C")
    background_df <- data.frame(gene = c("A", "B", "C"), copy_number = c(1000, 2000, 500))

    result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)

    if (!is.na(result$p_value) && result$status == "success") {
      tests_passed <- tests_passed + 1
    }
  }, error = function(e) {
    # Should handle or warn appropriately
  })

  list(
    passed = tests_passed == total_tests,
    passed_count = tests_passed,
    total_count = total_tests
  )
}

#' Test mathematical equivalence between parameter weighting and instance expansion
test_mathematical_equivalence <- function() {

  n_tests <- 10
  tests_passed <- 0

  for (i in 1:n_tests) {
    tryCatch({
      # Generate random test data
      n_genes <- sample(10:20, 1)
      genes <- paste0("GENE", 1:n_genes)
      copy_numbers <- sample(1:5, n_genes, replace = TRUE)

      background_df <- data.frame(gene = genes, copy_number = copy_numbers)

      # Random query subset
      query_size <- sample(3:8, 1)
      query_genes <- sample(genes, query_size)
      query_df <- background_df[background_df$gene %in% query_genes, ]

      # Random pathway
      pathway_size <- sample(3:8, 1)
      pathway_genes <- sample(genes, pathway_size)

      # Method 1: Parameter weighting
      result_weighted <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)

      # Method 2: Instance expansion
      query_expanded <- rep(query_df$gene, query_df$copy_number)
      background_expanded <- rep(background_df$gene, background_df$copy_number)

      q_exp <- sum(query_expanded %in% pathway_genes)
      m_exp <- sum(background_expanded %in% pathway_genes)
      n_exp <- length(background_expanded) - m_exp
      k_exp <- length(query_expanded)

      if (q_exp > 0) {
        p_val_expanded <- phyper(q_exp - 1, m_exp, n_exp, k_exp, lower.tail = FALSE)

        # Check equivalence (allowing for small numerical differences)
        if (abs(result_weighted$p_value - p_val_expanded) < 1e-10) {
          tests_passed <- tests_passed + 1
        }
      } else {
        # Both should give p-value = 1
        if (result_weighted$p_value == 1.0) {
          tests_passed <- tests_passed + 1
        }
      }

    }, error = function(e) {
      # Count as failed test
    })
  }

  list(
    passed = tests_passed == n_tests,
    passed_count = tests_passed,
    total_count = n_tests
  )
}

#' Test null distribution uniformity
test_null_distribution <- function() {

  # Generate background dataset
  n_genes <- 100
  genes <- paste0("GENE", 1:n_genes)
  copy_numbers <- sample(1:3, n_genes, replace = TRUE)  # Keep copy numbers small for speed
  background_df <- data.frame(gene = genes, copy_number = copy_numbers)

  # Random pathway
  pathway_genes <- sample(genes, 20)

  # Generate null p-values
  n_sim <- 200  # Reduced for speed
  null_pvals <- replicate(n_sim, {
    # Random query with no enrichment bias
    query_size <- sample(10:30, 1)
    query_genes <- sample(genes, query_size)
    query_df <- background_df[background_df$gene %in% query_genes, ]

    result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)
    return(result$p_value)
  })

  # Remove any NA values
  null_pvals <- null_pvals[!is.na(null_pvals)]

  if (length(null_pvals) < 100) {
    return(list(passed = FALSE, reason = "Too few valid p-values"))
  }

  # Test uniformity using Kolmogorov-Smirnov test
  ks_result <- ks.test(null_pvals, punif)

  # Also check basic quantiles
  quantiles <- quantile(null_pvals, c(0.05, 0.5, 0.95))
  quantile_check <- (quantiles[1] < 0.15 && quantiles[1] > 0.0) &&
                   (quantiles[2] < 0.65 && quantiles[2] > 0.35) &&
                   (quantiles[3] < 1.0 && quantiles[3] > 0.85)

  list(
    passed = ks_result$p.value > 0.01 && quantile_check,  # Allow some deviation
    ks_p_value = ks_result$p.value,
    ks_statistic = ks_result$statistic,
    quantiles = quantiles
  )
}

#' Test Type I error control
test_type_i_error_control <- function() {

  # Generate background dataset
  n_genes <- 80
  genes <- paste0("GENE", 1:n_genes)
  copy_numbers <- sample(1:3, n_genes, replace = TRUE)
  background_df <- data.frame(gene = genes, copy_number = copy_numbers)

  # Random pathway
  pathway_genes <- sample(genes, 15)

  # Test Type I error at alpha = 0.05
  alpha <- 0.05
  n_tests <- 500  # Reduced for speed

  false_positives <- sum(replicate(n_tests, {
    # Generate random query (null hypothesis)
    query_size <- sample(8:25, 1)
    query_genes <- sample(genes, query_size)
    query_df <- background_df[background_df$gene %in% query_genes, ]

    result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)
    return(result$p_value < alpha)
  }))

  observed_rate <- false_positives / n_tests
  tolerance <- 0.02  # Allow 2% deviation from nominal alpha

  list(
    passed = abs(observed_rate - alpha) <= tolerance,
    observed_rate = observed_rate,
    expected_rate = alpha,
    false_positives = false_positives,
    total_tests = n_tests,
    tolerance = tolerance
  )
}

#' Test performance benchmarks
test_performance_benchmarks <- function() {

  # Create test dataset
  n_genes <- 50
  genes <- paste0("GENE", 1:n_genes)
  copy_numbers <- sample(1:4, n_genes, replace = TRUE)
  background_df <- data.frame(gene = genes, copy_number = copy_numbers)

  query_genes <- sample(genes, 15)
  query_df <- background_df[background_df$gene %in% query_genes, ]
  pathway_genes <- sample(genes, 12)

  # Run benchmark
  benchmark_result <- benchmark_enrichment_methods(query_df, pathway_genes, background_df, n_reps = 20)

  # Check that parameter method is faster
  speedup <- benchmark_result$speedup_factor
  memory_reduction <- benchmark_result$memory_reduction

  list(
    passed = speedup > 1.0 && memory_reduction > 0,
    speedup_factor = speedup,
    memory_reduction = memory_reduction,
    parameter_time = benchmark_result$parameter_method$mean_time,
    expansion_time = benchmark_result$expansion_method$mean_time
  )
}

#' Generate test report
generate_validation_report <- function(results, output_file = "enrichment_validation_report.txt") {

  report <- c(
    "=== Copy-Number Weighted Enrichment Validation Report ===",
    paste("Generated:", Sys.time()),
    "",
    "OVERALL RESULT:",
    paste("  Status:", ifelse(results$overall_passed, "PASS", "FAIL")),
    "",
    "INDIVIDUAL TEST RESULTS:",
    ""
  )

  for (test_name in names(results)) {
    if (test_name == "overall_passed") next

    test_result <- results[[test_name]]
    status <- ifelse(test_result$passed, "PASS", "FAIL")
    report <- c(report, paste("", test_name, ":", status))

    if (!is.null(test_result$details)) {
      for (detail_name in names(test_result$details)) {
        report <- c(report, paste("    ", detail_name, ":", test_result$details[[detail_name]]))
      }
    }

    if (!is.null(test_result$error)) {
      report <- c(report, paste("    Error:", test_result$error))
    }

    report <- c(report, "")
  }

  # Write report
  cat(paste(report, collapse = "\n"), file = output_file)
  cat("Validation report written to:", output_file, "\n")

  return(report)
}

# Quick validation function for integration testing
quick_validation <- function() {
  cat("Running quick validation...\n")

  # Simple functional test - fix copy number balance
  query_df <- data.frame(gene = c("A", "B", "C"), copy_number = c(1, 1, 1))
  pathway_genes <- c("A", "C", "D")
  background_df <- data.frame(gene = LETTERS[1:10], copy_number = rep(2, 10))  # Higher copy numbers

  result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)

  if (result$status == "success" && !is.na(result$p_value)) {
    cat("Quick validation: PASS\n")
    return(TRUE)
  } else {
    cat("Quick validation: FAIL\n")
    cat("Status:", result$status, "\n")
    if (!is.null(result$error_message)) {
      cat("Error:", result$error_message, "\n")
    }
    return(FALSE)
  }
}