# Comprehensive Statistical Validation Suite for Copy-Number Weighted phyper()
#
# This suite integrates all statistical validation components into a unified framework
# for comprehensive validation of the copy-number weighted hypergeometric test.
#
# Author: AI Assistant
# Task: statistical-validation-framework-3
# Date: 2026-04-01
#
# COMPONENTS INTEGRATED:
# 1. Null Distribution Validation - Tests p-values follow Uniform(0,1) under H0
# 2. Type I Error Rate Validation - Tests false positive control at various alpha levels
# 3. Parameter Constraint Validation - Tests parameter bounds and edge cases
#
# USAGE:
#   source("statistical_validation_suite.R")
#   results <- run_comprehensive_validation()
#   generate_validation_report(results)

library(tidyverse)
library(yaml)

# Source all validation components
source("null_distribution_test.R")
source("type_i_error_validation.R")
source("parameter_constraints_validation.R")
source("edge_case_test_suite.R")
source("constraint_violation_handler.R")
source("debug_weighted_phyper.R")

# ==============================================================================
# COMPREHENSIVE VALIDATION FRAMEWORK
# ==============================================================================

#' Run comprehensive statistical validation suite
#'
#' Executes all validation components and generates integrated results
#'
#' @param validation_config List containing validation parameters (optional)
#' @param verbose Whether to print detailed progress (default: TRUE)
#' @return Comprehensive validation results list
run_comprehensive_validation <- function(validation_config = NULL, verbose = TRUE) {

  # Default configuration if not provided
  if (is.null(validation_config)) {
    validation_config <- list(
      null_distribution = list(
        n_simulations = 2000,
        scenarios = c("uniform", "skewed", "realistic"),
        background_sizes = c(500, 1000, 2000),
        pathway_sizes = c(50, 100, 200),
        query_sizes = c(25, 50, 100)
      ),
      type_i_error = list(
        n_simulations = 2000,
        alpha_levels = c(0.001, 0.01, 0.05, 0.1),
        tolerance = 0.01,
        scenarios = list(
          basic = list(n_genes = 1000, n_pathway = 100, n_query = 50),
          large = list(n_genes = 5000, n_pathway = 500, n_query = 250),
          small = list(n_genes = 200, n_pathway = 20, n_query = 10)
        )
      ),
      parameter_constraints = list(
        edge_case_scenarios = c("zero_copies", "extreme_ratios", "boundary_conditions"),
        constraint_tests = c("parameter_bounds", "feasibility", "consistency"),
        error_recovery_tests = TRUE
      ),
      output = list(
        save_detailed_results = TRUE,
        generate_plots = TRUE,
        export_summary = TRUE
      )
    )
  }

  if (verbose) {
    cat("================================================================================\n")
    cat("COMPREHENSIVE STATISTICAL VALIDATION SUITE\n")
    cat("Copy-Number Weighted Hypergeometric Test Validation\n")
    cat("================================================================================\n")
    cat("Timestamp:", as.character(Sys.time()), "\n\n")
  }

  # Initialize results container
  validation_results <- list(
    metadata = list(
      timestamp = Sys.time(),
      configuration = validation_config,
      r_version = R.version.string,
      package_versions = list()
    ),
    null_distribution = NULL,
    type_i_error = NULL,
    parameter_constraints = NULL,
    summary = NULL,
    overall_status = "PENDING"
  )

  # Track overall validation status
  all_tests_passed <- TRUE
  failed_components <- character(0)

  # ==============================================================================
  # 1. NULL DISTRIBUTION VALIDATION
  # ==============================================================================

  if (verbose) cat("1. Running null distribution validation...\n")

  tryCatch({
    null_results <- list()

    # Test across different scenarios
    for (scenario in validation_config$null_distribution$scenarios) {
      for (bg_size in validation_config$null_distribution$background_sizes) {
        for (pw_size in validation_config$null_distribution$pathway_sizes) {
          for (q_size in validation_config$null_distribution$query_sizes) {

            if (verbose) {
              cat(sprintf("  - Scenario: %s, Background: %d, Pathway: %d, Query: %d\n",
                         scenario, bg_size, pw_size, q_size))
            }

            # Run null distribution test
            test_key <- paste(scenario, bg_size, pw_size, q_size, sep = "_")

            pvalues <- simulate_null_pvalues(
              n_simulations = validation_config$null_distribution$n_simulations,
              background_size = bg_size,
              pathway_size = pw_size,
              query_size = q_size,
              copy_distribution = scenario
            )

            # Validate p-value distribution
            validation <- test_pvalue_uniformity(pvalues)
            null_results[[test_key]] <- list(
              scenario = scenario,
              parameters = list(bg = bg_size, pw = pw_size, q = q_size),
              pvalues = pvalues,
              validation = validation,
              passed = validation$overall_uniform == TRUE
            )

            if (!null_results[[test_key]]$passed) {
              all_tests_passed <- FALSE
              failed_components <- c(failed_components, paste("null_dist", test_key, sep = "_"))
            }
          }
        }
      }
    }

    validation_results$null_distribution <- null_results

    if (verbose) {
      passed_tests <- sum(sapply(null_results, function(x) x$passed))
      total_tests <- length(null_results)
      cat(sprintf("  - Null distribution tests: %d/%d passed\n", passed_tests, total_tests))
    }

  }, error = function(e) {
    cat("ERROR in null distribution validation:", e$message, "\n")
    validation_results$null_distribution <- list(error = e$message)
    all_tests_passed <- FALSE
    failed_components <- c(failed_components, "null_distribution_error")
  })

  # ==============================================================================
  # 2. TYPE I ERROR RATE VALIDATION
  # ==============================================================================

  if (verbose) cat("\n2. Running Type I error rate validation...\n")

  tryCatch({
    type_i_results <- list()

    # Test across different scenarios
    for (scenario_name in names(validation_config$type_i_error$scenarios)) {
      scenario_config <- validation_config$type_i_error$scenarios[[scenario_name]]

      if (verbose) {
        cat(sprintf("  - Scenario: %s (%d genes, %d pathway, %d query)\n",
                   scenario_name, scenario_config$n_genes,
                   scenario_config$n_pathway, scenario_config$n_query))
      }

      # Run Type I error validation
      type_i_validation <- validate_type_i_error_rates(
        n_sims = validation_config$type_i_error$n_simulations,
        alphas = validation_config$type_i_error$alpha_levels,
        n_genes = scenario_config$n_genes,
        n_pathway = scenario_config$n_pathway,
        n_query = scenario_config$n_query,
        tolerance = validation_config$type_i_error$tolerance,
        verbose = FALSE
      )

      type_i_results[[scenario_name]] <- type_i_validation

      if (!type_i_validation$overall_status$all_tests_passed) {
        all_tests_passed <- FALSE
        failed_components <- c(failed_components, paste("type_i_error", scenario_name, sep = "_"))
      }
    }

    validation_results$type_i_error <- type_i_results

    if (verbose) {
      passed_scenarios <- sum(sapply(type_i_results, function(x) x$overall_status$all_tests_passed))
      total_scenarios <- length(type_i_results)
      cat(sprintf("  - Type I error scenarios: %d/%d passed\n", passed_scenarios, total_scenarios))
    }

  }, error = function(e) {
    cat("ERROR in Type I error validation:", e$message, "\n")
    validation_results$type_i_error <- list(error = e$message)
    all_tests_passed <- FALSE
    failed_components <- c(failed_components, "type_i_error_error")
  })

  # ==============================================================================
  # 3. PARAMETER CONSTRAINT VALIDATION
  # ==============================================================================

  if (verbose) cat("\n3. Running parameter constraint validation...\n")

  tryCatch({
    constraint_results <- list()

    # Test parameter bounds validation
    if (verbose) cat("  - Testing parameter bounds...\n")
    bounds_validation <- run_constraint_validation_tests()
    constraint_results$parameter_bounds <- bounds_validation

    if (!bounds_validation$all_tests_passed) {
      all_tests_passed <- FALSE
      failed_components <- c(failed_components, "parameter_bounds")
    }

    # Test edge cases
    if (verbose) cat("  - Testing edge cases...\n")
    edge_case_validation <- run_comprehensive_edge_case_tests()
    constraint_results$edge_cases <- edge_case_validation

    if (!edge_case_validation$all_tests_passed) {
      all_tests_passed <- FALSE
      failed_components <- c(failed_components, "edge_cases")
    }

    # Test constraint violation handling
    if (validation_config$parameter_constraints$error_recovery_tests) {
      if (verbose) cat("  - Testing constraint violation handling...\n")
      violation_validation <- test_constraint_violation_handling()
      constraint_results$violation_handling <- violation_validation

      if (!violation_validation$all_tests_passed) {
        all_tests_passed <- FALSE
        failed_components <- c(failed_components, "violation_handling")
      }
    }

    validation_results$parameter_constraints <- constraint_results

    if (verbose) {
      passed_tests <- sum(sapply(constraint_results, function(x) {
        if ("all_constraints_satisfied" %in% names(x)) return(x$all_constraints_satisfied)
        if ("all_tests_passed" %in% names(x)) return(x$all_tests_passed)
        if ("all_recovery_tests_passed" %in% names(x)) return(x$all_recovery_tests_passed)
        return(FALSE)
      }))
      total_tests <- length(constraint_results)
      cat(sprintf("  - Parameter constraint tests: %d/%d passed\n", passed_tests, total_tests))
    }

  }, error = function(e) {
    cat("ERROR in parameter constraint validation:", e$message, "\n")
    validation_results$parameter_constraints <- list(error = e$message)
    all_tests_passed <- FALSE
    failed_components <- c(failed_components, "parameter_constraints_error")
  })

  # ==============================================================================
  # 4. GENERATE SUMMARY AND FINAL STATUS
  # ==============================================================================

  validation_results$overall_status <- if (all_tests_passed) "PASS" else "FAIL"
  validation_results$failed_components <- failed_components

  # Generate summary statistics
  validation_results$summary <- generate_validation_summary(validation_results)

  if (verbose) {
    cat("\n================================================================================\n")
    cat("VALIDATION SUMMARY\n")
    cat("================================================================================\n")
    cat(sprintf("Overall Status: %s\n", validation_results$overall_status))

    if (length(failed_components) > 0) {
      cat("Failed Components:\n")
      for (comp in failed_components) {
        cat(sprintf("  - %s\n", comp))
      }
    }

    cat("\nComponent Results:\n")
    if (!is.null(validation_results$null_distribution) && !("error" %in% names(validation_results$null_distribution))) {
      null_passed <- sum(sapply(validation_results$null_distribution, function(x) x$passed))
      null_total <- length(validation_results$null_distribution)
      cat(sprintf("  - Null Distribution: %d/%d tests passed\n", null_passed, null_total))
    }

    if (!is.null(validation_results$type_i_error) && !("error" %in% names(validation_results$type_i_error))) {
      type_i_passed <- sum(sapply(validation_results$type_i_error, function(x) x$overall_status$all_tests_passed))
      type_i_total <- length(validation_results$type_i_error)
      cat(sprintf("  - Type I Error: %d/%d scenarios passed\n", type_i_passed, type_i_total))
    }

    if (!is.null(validation_results$parameter_constraints) && !("error" %in% names(validation_results$parameter_constraints))) {
      param_passed <- sum(sapply(validation_results$parameter_constraints, function(x) {
        if ("all_constraints_satisfied" %in% names(x)) return(x$all_constraints_satisfied)
        if ("all_tests_passed" %in% names(x)) return(x$all_tests_passed)
        if ("all_recovery_tests_passed" %in% names(x)) return(x$all_recovery_tests_passed)
        return(FALSE)
      }))
      param_total <- length(validation_results$parameter_constraints)
      cat(sprintf("  - Parameter Constraints: %d/%d tests passed\n", param_passed, param_total))
    }

    cat("================================================================================\n")
  }

  # Save results if requested
  if (validation_config$output$save_detailed_results) {
    save(validation_results, file = "comprehensive_validation_results.RData")
    if (verbose) cat("Detailed results saved to: comprehensive_validation_results.RData\n")
  }

  return(validation_results)
}

# ==============================================================================
# SUMMARY AND REPORTING FUNCTIONS
# ==============================================================================

#' Generate validation summary statistics
#'
#' @param validation_results Results from run_comprehensive_validation()
#' @return Summary statistics list
generate_validation_summary <- function(validation_results) {
  summary <- list(
    timestamp = validation_results$metadata$timestamp,
    overall_status = validation_results$overall_status,
    total_tests = 0,
    passed_tests = 0,
    failed_tests = 0,
    component_summary = list()
  )

  # Null distribution summary
  if (!is.null(validation_results$null_distribution) &&
      !("error" %in% names(validation_results$null_distribution))) {
    null_tests <- length(validation_results$null_distribution)
    null_passed <- sum(sapply(validation_results$null_distribution, function(x) x$passed))

    summary$component_summary$null_distribution <- list(
      total_tests = null_tests,
      passed_tests = null_passed,
      pass_rate = null_passed / null_tests
    )

    summary$total_tests <- summary$total_tests + null_tests
    summary$passed_tests <- summary$passed_tests + null_passed
  }

  # Type I error summary
  if (!is.null(validation_results$type_i_error) &&
      !("error" %in% names(validation_results$type_i_error))) {
    type_i_scenarios <- length(validation_results$type_i_error)
    type_i_passed <- sum(sapply(validation_results$type_i_error, function(x) x$overall_status$all_tests_passed))

    summary$component_summary$type_i_error <- list(
      total_scenarios = type_i_scenarios,
      passed_scenarios = type_i_passed,
      pass_rate = type_i_passed / type_i_scenarios
    )

    summary$total_tests <- summary$total_tests + type_i_scenarios
    summary$passed_tests <- summary$passed_tests + type_i_passed
  }

  # Parameter constraints summary
  if (!is.null(validation_results$parameter_constraints) &&
      !("error" %in% names(validation_results$parameter_constraints))) {
    param_tests <- length(validation_results$parameter_constraints)
    param_passed <- sum(sapply(validation_results$parameter_constraints, function(x) {
      if ("all_tests_passed" %in% names(x)) return(x$all_tests_passed)
      if ("all_constraints_satisfied" %in% names(x)) return(x$all_constraints_satisfied)
      if ("all_recovery_tests_passed" %in% names(x)) return(x$all_recovery_tests_passed)
      return(FALSE)
    }))

    summary$component_summary$parameter_constraints <- list(
      total_tests = param_tests,
      passed_tests = param_passed,
      pass_rate = param_passed / param_tests
    )

    summary$total_tests <- summary$total_tests + param_tests
    summary$passed_tests <- summary$passed_tests + param_passed
  }

  summary$failed_tests <- summary$total_tests - summary$passed_tests
  summary$overall_pass_rate <- if (summary$total_tests > 0) summary$passed_tests / summary$total_tests else 0

  return(summary)
}

#' Generate comprehensive validation report
#'
#' Creates detailed markdown report of validation results
#'
#' @param validation_results Results from run_comprehensive_validation()
#' @param output_file Output file path (default: "statistical_validation_report.md")
#' @return Path to generated report file
generate_validation_report <- function(validation_results, output_file = "statistical_validation_report.md") {

  report_content <- c(
    "# Comprehensive Statistical Validation Report",
    "## Copy-Number Weighted Hypergeometric Test",
    "",
    paste("**Generated:** ", as.character(validation_results$metadata$timestamp)),
    paste("**Overall Status:** ", validation_results$overall_status),
    "",
    "## Executive Summary",
    ""
  )

  # Add executive summary
  summary <- validation_results$summary
  report_content <- c(report_content,
    paste("This report presents the results of comprehensive statistical validation for the copy-number weighted hypergeometric test implementation. A total of", summary$total_tests, "validation tests were executed across three major validation domains:"),
    "",
    "1. **Null Distribution Validation** - Verifying p-values follow Uniform(0,1) under null hypothesis",
    "2. **Type I Error Rate Control** - Validating false positive rates at multiple significance levels",
    "3. **Parameter Constraint Validation** - Testing parameter bounds and edge case handling",
    "",
    paste("**Overall Results:** ", summary$passed_tests, "/", summary$total_tests, " tests passed (", round(summary$overall_pass_rate * 100, 1), "%)"),
    ""
  )

  if (validation_results$overall_status == "FAIL") {
    report_content <- c(report_content,
      "### ⚠️ VALIDATION FAILURES DETECTED",
      "",
      "The following validation components failed:",
      ""
    )

    for (failed_comp in validation_results$failed_components) {
      report_content <- c(report_content, paste("- ", failed_comp))
    }

    report_content <- c(report_content, "")
  }

  # Detailed results sections
  report_content <- c(report_content,
    "## 1. Null Distribution Validation Results",
    ""
  )

  if (!is.null(validation_results$null_distribution) &&
      !("error" %in% names(validation_results$null_distribution))) {

    null_summary <- summary$component_summary$null_distribution
    report_content <- c(report_content,
      paste("**Tests:** ", null_summary$total_tests, " scenarios tested"),
      paste("**Passed:** ", null_summary$passed_tests, " (", round(null_summary$pass_rate * 100, 1), "%)"),
      ""
    )

    # Add scenario details
    for (test_name in names(validation_results$null_distribution)) {
      test_result <- validation_results$null_distribution[[test_name]]
      status_icon <- if (test_result$passed) "✅" else "❌"

      report_content <- c(report_content,
        paste("### Scenario:", test_name, status_icon),
        paste("- Copy Distribution:", test_result$scenario),
        paste("- Parameters: Background =", test_result$parameters$bg,
              ", Pathway =", test_result$parameters$pw,
              ", Query =", test_result$parameters$q),
        paste("- Result:", test_result$validation$overall_conclusion),
        ""
      )
    }
  } else {
    report_content <- c(report_content, "❌ **ERROR:** Null distribution validation failed to execute properly", "")
  }

  # Type I Error section
  report_content <- c(report_content,
    "## 2. Type I Error Rate Validation Results",
    ""
  )

  if (!is.null(validation_results$type_i_error) &&
      !("error" %in% names(validation_results$type_i_error))) {

    type_i_summary <- summary$component_summary$type_i_error
    report_content <- c(report_content,
      paste("**Scenarios:** ", type_i_summary$total_scenarios, " tested"),
      paste("**Passed:** ", type_i_summary$passed_scenarios, " (", round(type_i_summary$pass_rate * 100, 1), "%)"),
      ""
    )

    # Add scenario details
    for (scenario_name in names(validation_results$type_i_error)) {
      scenario_result <- validation_results$type_i_error[[scenario_name]]
      status_icon <- if (scenario_result$overall_status$all_tests_passed) "✅" else "❌"

      report_content <- c(report_content,
        paste("### Scenario:", scenario_name, status_icon),
        paste("- Status:", scenario_result$overall_status$conclusion),
        ""
      )
    }
  } else {
    report_content <- c(report_content, "❌ **ERROR:** Type I error validation failed to execute properly", "")
  }

  # Parameter Constraints section
  report_content <- c(report_content,
    "## 3. Parameter Constraint Validation Results",
    ""
  )

  if (!is.null(validation_results$parameter_constraints) &&
      !("error" %in% names(validation_results$parameter_constraints))) {

    param_summary <- summary$component_summary$parameter_constraints
    report_content <- c(report_content,
      paste("**Test Categories:** ", param_summary$total_tests),
      paste("**Passed:** ", param_summary$passed_tests, " (", round(param_summary$pass_rate * 100, 1), "%)"),
      ""
    )

    # Add component details
    for (test_name in names(validation_results$parameter_constraints)) {
      test_result <- validation_results$parameter_constraints[[test_name]]

      if ("all_tests_passed" %in% names(test_result)) {
        status_icon <- if (test_result$all_tests_passed) "✅" else "❌"
        status_text <- if (test_result$all_tests_passed) "PASS" else "FAIL"
      } else if ("all_constraints_satisfied" %in% names(test_result)) {
        status_icon <- if (test_result$all_constraints_satisfied) "✅" else "❌"
        status_text <- if (test_result$all_constraints_satisfied) "PASS" else "FAIL"
      } else if ("all_recovery_tests_passed" %in% names(test_result)) {
        status_icon <- if (test_result$all_recovery_tests_passed) "✅" else "❌"
        status_text <- if (test_result$all_recovery_tests_passed) "PASS" else "FAIL"
      } else {
        status_icon <- "❓"
        status_text <- "UNKNOWN"
      }

      report_content <- c(report_content,
        paste("### ", str_to_title(gsub("_", " ", test_name)), status_icon),
        paste("- Status:", status_text),
        ""
      )
    }
  } else {
    report_content <- c(report_content, "❌ **ERROR:** Parameter constraint validation failed to execute properly", "")
  }

  # Conclusions and recommendations
  report_content <- c(report_content,
    "## Conclusions and Recommendations",
    ""
  )

  if (validation_results$overall_status == "PASS") {
    report_content <- c(report_content,
      "### ✅ VALIDATION SUCCESSFUL",
      "",
      "All statistical validation tests have passed. The copy-number weighted hypergeometric test implementation demonstrates:",
      "",
      "- **Correct null distribution behavior** - P-values follow expected Uniform(0,1) distribution under null hypothesis",
      "- **Proper Type I error control** - False positive rates are maintained at nominal levels across tested scenarios",
      "- **Robust parameter handling** - Parameter constraints are properly enforced and edge cases handled gracefully",
      "",
      "**Recommendation:** The implementation is statistically valid and ready for production use.",
      ""
    )
  } else {
    report_content <- c(report_content,
      "### ❌ VALIDATION ISSUES DETECTED",
      "",
      "Statistical validation has identified issues that must be addressed:",
      ""
    )

    for (failed_comp in validation_results$failed_components) {
      report_content <- c(report_content, paste("- ", failed_comp))
    }

    report_content <- c(report_content,
      "",
      "**Recommendation:** Address identified issues before using this implementation in production.",
      ""
    )
  }

  # Technical details
  report_content <- c(report_content,
    "## Technical Details",
    "",
    paste("- **R Version:** ", validation_results$metadata$r_version),
    paste("- **Validation Framework Version:** statistical-validation-framework-3"),
    paste("- **Total Validation Tests:** ", summary$total_tests),
    paste("- **Test Configuration:** Comprehensive validation across multiple scenarios"),
    "",
    "### Files Generated",
    "",
    "- `statistical_validation_suite.R` - Integrated validation framework",
    "- `comprehensive_validation_results.RData` - Detailed validation results",
    "- `statistical_validation_report.md` - This validation report",
    "",
    "---",
    "*Report generated by statistical validation framework v3*"
  )

  # Write report to file
  writeLines(report_content, output_file)

  return(output_file)
}

# ==============================================================================
# USAGE DOCUMENTATION AND GUIDELINES
# ==============================================================================

#' Print usage guidelines for the validation framework
#'
#' Displays comprehensive usage instructions and interpretation guidelines
print_validation_usage_guide <- function() {
  cat("================================================================================\n")
  cat("STATISTICAL VALIDATION FRAMEWORK - USAGE GUIDE\n")
  cat("================================================================================\n")
  cat("\n")
  cat("QUICK START:\n")
  cat("  source('statistical_validation_suite.R')\n")
  cat("  results <- run_comprehensive_validation()\n")
  cat("  generate_validation_report(results)\n")
  cat("\n")
  cat("DETAILED USAGE:\n")
  cat("\n")
  cat("1. BASIC VALIDATION:\n")
  cat("   # Run with default parameters (recommended for most users)\n")
  cat("   results <- run_comprehensive_validation()\n")
  cat("\n")
  cat("2. CUSTOM VALIDATION:\n")
  cat("   # Customize validation parameters\n")
  cat("   config <- list(\n")
  cat("     null_distribution = list(n_simulations = 5000),\n")
  cat("     type_i_error = list(alpha_levels = c(0.01, 0.05)),\n")
  cat("     output = list(save_detailed_results = TRUE)\n")
  cat("   )\n")
  cat("   results <- run_comprehensive_validation(config)\n")
  cat("\n")
  cat("3. GENERATE REPORTS:\n")
  cat("   # Create markdown report\n")
  cat("   report_file <- generate_validation_report(results)\n")
  cat("   \n")
  cat("   # View summary\n")
  cat("   print(results$summary)\n")
  cat("\n")
  cat("INTERPRETATION GUIDELINES:\n")
  cat("\n")
  cat("1. OVERALL STATUS:\n")
  cat("   - PASS: All validation tests passed - implementation is statistically valid\n")
  cat("   - FAIL: One or more tests failed - issues must be addressed\n")
  cat("\n")
  cat("2. NULL DISTRIBUTION VALIDATION:\n")
  cat("   - Tests whether p-values follow Uniform(0,1) under null hypothesis\n")
  cat("   - Uses Kolmogorov-Smirnov and Anderson-Darling tests\n")
  cat("   - Failure indicates fundamental statistical issues\n")
  cat("\n")
  cat("3. TYPE I ERROR VALIDATION:\n")
  cat("   - Tests false positive rates at multiple significance levels\n")
  cat("   - Validates error rate control within 1% of nominal levels\n")
  cat("   - Failure indicates inflated or deflated Type I error rates\n")
  cat("\n")
  cat("4. PARAMETER CONSTRAINT VALIDATION:\n")
  cat("   - Tests parameter bounds enforcement\n")
  cat("   - Validates edge case handling and error recovery\n")
  cat("   - Failure indicates robustness issues\n")
  cat("\n")
  cat("TROUBLESHOOTING:\n")
  cat("\n")
  cat("1. NULL DISTRIBUTION FAILURES:\n")
  cat("   - Check weighted parameter calculations\n")
  cat("   - Verify background gene selection\n")
  cat("   - Review hypergeometric parameter mapping\n")
  cat("\n")
  cat("2. TYPE I ERROR FAILURES:\n")
  cat("   - Check null hypothesis implementation\n")
  cat("   - Verify query gene sampling\n")
  cat("   - Review statistical test implementation\n")
  cat("\n")
  cat("3. PARAMETER CONSTRAINT FAILURES:\n")
  cat("   - Check parameter validation logic\n")
  cat("   - Verify edge case handling\n")
  cat("   - Review error recovery mechanisms\n")
  cat("\n")
  cat("================================================================================\n")
}

# Display usage guide when framework is loaded
cat("Statistical Validation Suite loaded successfully!\n")
cat("Run print_validation_usage_guide() for detailed usage instructions.\n")
cat("Quick start: run_comprehensive_validation()\n")