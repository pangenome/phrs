# Constraint Violation Handling Functions for Copy-Number Weighted phyper()
# Robust error handling and recovery strategies for parameter constraint violations

library(tidyverse)
source("parameter_constraints_validation.R")

# Error recovery strategies
get_violation_recovery_strategy <- function(violation_type) {
  # Define recovery strategies for different constraint violations
  recovery_strategies <- list(
    "negative_parameters" = list(
      action = "filter_and_warn",
      description = "Remove genes with negative copy numbers",
      severity = "error"
    ),
    "non_integer_parameters" = list(
      action = "round_and_warn",
      description = "Round copy numbers to nearest integers",
      severity = "warning"
    ),
    "overlap_bounds_violation" = list(
      action = "data_validation_error",
      description = "Fundamental data inconsistency - requires manual review",
      severity = "critical"
    ),
    "sample_bounds_violation" = list(
      action = "data_validation_error",
      description = "Query larger than background - check data integrity",
      severity = "critical"
    ),
    "feasibility_violation" = list(
      action = "parameter_adjustment",
      description = "Adjust parameters to feasible range",
      severity = "warning"
    ),
    "zero_overlap" = list(
      action = "return_null_result",
      description = "No enrichment possible with zero overlap",
      severity = "info"
    ),
    "extreme_values" = list(
      action = "proceed_with_warning",
      description = "Continue with numerical stability warnings",
      severity = "warning"
    )
  )

  return(recovery_strategies[[violation_type]])
}

# Constraint violation handler with recovery
handle_constraint_violation <- function(violation_result, recovery_mode = "strict") {
  # Handle constraint violations with different recovery strategies
  #
  # Parameters:
  #   violation_result: Result from detect_constraint_violations()
  #   recovery_mode: "strict" (fail on any violation), "lenient" (attempt recovery),
  #                  "permissive" (proceed with warnings)
  #
  # Returns:
  #   List with handling result and any corrected data

  handler_result <- list(
    can_proceed = FALSE,
    corrected_data = NULL,
    applied_corrections = character(0),
    warnings = character(0),
    errors = character(0),
    recovery_actions = character(0)
  )

  # If validation passed, no handling needed
  if (violation_result$valid) {
    handler_result$can_proceed <- TRUE
    handler_result$warnings <- c(handler_result$warnings, violation_result$combined_warnings)
    return(handler_result)
  }

  # Process gene-level violations
  if (length(violation_result$gene_violations) > 0) {
    for (violation in violation_result$gene_violations) {
      if (grepl("Negative copy numbers", violation)) {
        if (recovery_mode %in% c("lenient", "permissive")) {
          handler_result$recovery_actions <- c(
            handler_result$recovery_actions,
            "Applied negative copy number filtering"
          )
        } else {
          handler_result$errors <- c(handler_result$errors, violation)
        }
      } else if (grepl("Duplicate gene identifiers", violation)) {
        handler_result$errors <- c(handler_result$errors, violation)
        # Duplicate genes are always critical errors
      } else {
        handler_result$warnings <- c(handler_result$warnings, violation)
      }
    }
  }

  # Process parameter-level violations
  if (length(violation_result$parameter_violations) > 0) {
    for (violation in violation_result$parameter_violations) {
      if (grepl("non-negative", violation)) {
        strategy <- get_violation_recovery_strategy("negative_parameters")
      } else if (grepl("integer", violation)) {
        strategy <- get_violation_recovery_strategy("non_integer_parameters")
      } else if (grepl("Overlap constraint", violation)) {
        strategy <- get_violation_recovery_strategy("overlap_bounds_violation")
      } else if (grepl("Sample constraint", violation)) {
        strategy <- get_violation_recovery_strategy("sample_bounds_violation")
      } else if (grepl("Feasibility constraint", violation)) {
        strategy <- get_violation_recovery_strategy("feasibility_violation")
      } else {
        strategy <- list(action = "unknown", severity = "error")
      }

      # Apply recovery strategy based on mode
      if (strategy$severity == "critical" || recovery_mode == "strict") {
        handler_result$errors <- c(handler_result$errors, violation)
      } else if (strategy$action == "proceed_with_warning" && recovery_mode == "permissive") {
        handler_result$warnings <- c(handler_result$warnings, violation)
        handler_result$can_proceed <- TRUE
      } else {
        handler_result$warnings <- c(handler_result$warnings, paste("Recoverable:", violation))
        handler_result$recovery_actions <- c(handler_result$recovery_actions, strategy$description)
      }
    }
  }

  # Handle edge cases
  if (length(violation_result$edge_cases) > 0) {
    for (edge_case in violation_result$edge_cases) {
      strategy <- get_violation_recovery_strategy("zero_overlap")  # Default for edge cases
      handler_result$warnings <- c(handler_result$warnings, paste("Edge case detected:", edge_case))
      handler_result$can_proceed <- TRUE  # Edge cases are usually handleable
    }
  }

  # Final decision on whether to proceed
  if (length(handler_result$errors) == 0) {
    handler_result$can_proceed <- TRUE
  }

  return(handler_result)
}

# Data correction functions
correct_negative_copy_numbers <- function(df) {
  # Remove or correct negative copy numbers
  original_nrow <- nrow(df)

  # Option 1: Filter out negative values
  corrected_df <- df[df$copy_number >= 0, ]

  correction_info <- list(
    method = "filter_negative",
    original_rows = original_nrow,
    corrected_rows = nrow(corrected_df),
    genes_removed = original_nrow - nrow(corrected_df)
  )

  return(list(data = corrected_df, correction_info = correction_info))
}

correct_non_integer_copy_numbers <- function(df) {
  # Round non-integer copy numbers
  corrected_df <- df
  original_copies <- df$copy_number
  corrected_df$copy_number <- round(df$copy_number)

  correction_info <- list(
    method = "round_to_integer",
    max_adjustment = max(abs(original_copies - corrected_df$copy_number)),
    genes_adjusted = sum(original_copies != corrected_df$copy_number)
  )

  return(list(data = corrected_df, correction_info = correction_info))
}

correct_duplicate_genes <- function(df) {
  # Handle duplicate gene identifiers by aggregating copy numbers
  if (!any(duplicated(df$gene))) {
    return(list(data = df, correction_info = list(method = "no_duplicates")))
  }

  original_nrow <- nrow(df)
  corrected_df <- df %>%
    group_by(gene) %>%
    summarise(copy_number = sum(copy_number), .groups = "drop")

  correction_info <- list(
    method = "aggregate_duplicates",
    original_rows = original_nrow,
    corrected_rows = nrow(corrected_df),
    duplicates_merged = original_nrow - nrow(corrected_df)
  )

  return(list(data = corrected_df, correction_info = correction_info))
}

# Comprehensive data correction pipeline
apply_data_corrections <- function(query_df, background_df, correction_options = list()) {
  # Apply a series of corrections to make data valid for hypergeometric testing

  # Default correction options
  default_options <- list(
    correct_negatives = TRUE,
    correct_non_integers = TRUE,
    correct_duplicates = TRUE,
    filter_query_to_background = TRUE
  )

  options <- modifyList(default_options, correction_options)

  correction_log <- list(
    applied_corrections = character(0),
    correction_details = list()
  )

  corrected_query <- query_df
  corrected_background <- background_df

  # Correction 1: Handle negative copy numbers
  if (options$correct_negatives) {
    if (any(corrected_query$copy_number < 0, na.rm = TRUE)) {
      query_correction <- correct_negative_copy_numbers(corrected_query)
      corrected_query <- query_correction$data
      correction_log$applied_corrections <- c(correction_log$applied_corrections, "query_negative_filter")
      correction_log$correction_details$query_negative <- query_correction$correction_info
    }

    if (any(corrected_background$copy_number < 0, na.rm = TRUE)) {
      bg_correction <- correct_negative_copy_numbers(corrected_background)
      corrected_background <- bg_correction$data
      correction_log$applied_corrections <- c(correction_log$applied_corrections, "background_negative_filter")
      correction_log$correction_details$background_negative <- bg_correction$correction_info
    }
  }

  # Correction 2: Handle non-integer copy numbers
  if (options$correct_non_integers) {
    if (!all(corrected_query$copy_number == as.integer(corrected_query$copy_number), na.rm = TRUE)) {
      query_correction <- correct_non_integer_copy_numbers(corrected_query)
      corrected_query <- query_correction$data
      correction_log$applied_corrections <- c(correction_log$applied_corrections, "query_integer_rounding")
      correction_log$correction_details$query_integer <- query_correction$correction_info
    }

    if (!all(corrected_background$copy_number == as.integer(corrected_background$copy_number), na.rm = TRUE)) {
      bg_correction <- correct_non_integer_copy_numbers(corrected_background)
      corrected_background <- bg_correction$data
      correction_log$applied_corrections <- c(correction_log$applied_corrections, "background_integer_rounding")
      correction_log$correction_details$background_integer <- bg_correction$correction_info
    }
  }

  # Correction 3: Handle duplicate genes
  if (options$correct_duplicates) {
    if (any(duplicated(corrected_query$gene))) {
      query_correction <- correct_duplicate_genes(corrected_query)
      corrected_query <- query_correction$data
      correction_log$applied_corrections <- c(correction_log$applied_corrections, "query_duplicate_merge")
      correction_log$correction_details$query_duplicates <- query_correction$correction_info
    }

    if (any(duplicated(corrected_background$gene))) {
      bg_correction <- correct_duplicate_genes(corrected_background)
      corrected_background <- bg_correction$data
      correction_log$applied_corrections <- c(correction_log$applied_corrections, "background_duplicate_merge")
      correction_log$correction_details$background_duplicates <- bg_correction$correction_info
    }
  }

  # Correction 4: Filter query to background genes
  if (options$filter_query_to_background) {
    original_query_size <- nrow(corrected_query)
    valid_query_genes <- intersect(corrected_query$gene, corrected_background$gene)
    corrected_query <- corrected_query[corrected_query$gene %in% valid_query_genes, ]

    if (nrow(corrected_query) < original_query_size) {
      correction_log$applied_corrections <- c(correction_log$applied_corrections, "query_background_filter")
      correction_log$correction_details$query_filter <- list(
        original_genes = original_query_size,
        filtered_genes = nrow(corrected_query),
        genes_removed = original_query_size - nrow(corrected_query)
      )
    }
  }

  return(list(
    query_df = corrected_query,
    background_df = corrected_background,
    correction_log = correction_log
  ))
}

# Robust constraint validation with error recovery
robust_constraint_validation <- function(query_df, pathway_genes, background_df,
                                       recovery_mode = "lenient",
                                       max_correction_attempts = 3) {
  # Perform constraint validation with automatic error recovery

  validation_attempts <- list()

  for (attempt in 1:max_correction_attempts) {
    cat(sprintf("Validation attempt %d/%d...\n", attempt, max_correction_attempts))

    # Detect violations
    violation_result <- detect_constraint_violations(query_df, pathway_genes, background_df)

    # Handle violations
    handler_result <- handle_constraint_violation(violation_result, recovery_mode)

    validation_attempts[[attempt]] <- list(
      attempt_number = attempt,
      violation_result = violation_result,
      handler_result = handler_result
    )

    # If validation successful, return result
    if (handler_result$can_proceed && violation_result$valid) {
      cat("Validation successful!\n")
      return(list(
        success = TRUE,
        final_violation_result = violation_result,
        final_handler_result = handler_result,
        attempts = validation_attempts,
        corrected_data = NULL
      ))
    }

    # If unrecoverable errors, fail
    if (!handler_result$can_proceed) {
      cat("Unrecoverable validation errors detected.\n")
      return(list(
        success = FALSE,
        final_violation_result = violation_result,
        final_handler_result = handler_result,
        attempts = validation_attempts,
        critical_errors = handler_result$errors
      ))
    }

    # Attempt corrections for next iteration
    if (attempt < max_correction_attempts) {
      cat("Applying data corrections for next attempt...\n")
      correction_result <- apply_data_corrections(query_df, background_df)
      query_df <- correction_result$query_df
      background_df <- correction_result$background_df

      # Update pathway genes to match corrected background
      pathway_genes <- intersect(pathway_genes, background_df$gene)

      validation_attempts[[attempt]]$correction_applied <- correction_result$correction_log
    }
  }

  # Max attempts reached
  cat(sprintf("Maximum correction attempts (%d) reached.\n", max_correction_attempts))
  return(list(
    success = FALSE,
    final_violation_result = violation_result,
    final_handler_result = handler_result,
    attempts = validation_attempts,
    max_attempts_reached = TRUE
  ))
}

# Generate detailed validation report
generate_constraint_validation_report <- function(validation_result, output_file = NULL) {
  # Generate a comprehensive report of constraint validation results

  report_lines <- character(0)

  # Header
  report_lines <- c(report_lines,
    "========================================",
    "CONSTRAINT VALIDATION REPORT",
    "========================================",
    paste("Generated:", Sys.time()),
    paste("Validation Success:", validation_result$success),
    ""
  )

  # Summary
  if (validation_result$success) {
    report_lines <- c(report_lines,
      "VALIDATION SUMMARY: PASSED",
      "All constraints satisfied after validation/correction process.",
      ""
    )
  } else {
    report_lines <- c(report_lines,
      "VALIDATION SUMMARY: FAILED",
      "Critical constraint violations could not be resolved.",
      ""
    )
  }

  # Attempt details
  report_lines <- c(report_lines,
    sprintf("VALIDATION ATTEMPTS: %d", length(validation_result$attempts)),
    ""
  )

  for (i in seq_along(validation_result$attempts)) {
    attempt <- validation_result$attempts[[i]]
    report_lines <- c(report_lines,
      sprintf("Attempt %d:", i),
      sprintf("  Constraints satisfied: %s", attempt$violation_result$valid),
      sprintf("  Can proceed: %s", attempt$handler_result$can_proceed)
    )

    if (length(attempt$violation_result$gene_violations) > 0) {
      report_lines <- c(report_lines, "  Gene-level violations:")
      for (violation in attempt$violation_result$gene_violations) {
        report_lines <- c(report_lines, paste("    -", violation))
      }
    }

    if (length(attempt$violation_result$parameter_violations) > 0) {
      report_lines <- c(report_lines, "  Parameter violations:")
      for (violation in attempt$violation_result$parameter_violations) {
        report_lines <- c(report_lines, paste("    -", violation))
      }
    }

    if (length(attempt$handler_result$recovery_actions) > 0) {
      report_lines <- c(report_lines, "  Recovery actions:")
      for (action in attempt$handler_result$recovery_actions) {
        report_lines <- c(report_lines, paste("    -", action))
      }
    }

    if ("correction_applied" %in% names(attempt)) {
      report_lines <- c(report_lines, "  Corrections applied:")
      for (correction in attempt$correction_applied$applied_corrections) {
        report_lines <- c(report_lines, paste("    -", correction))
      }
    }

    report_lines <- c(report_lines, "")
  }

  # Final parameters if successful
  if (validation_result$success && !is.null(validation_result$final_violation_result$parameters)) {
    params <- validation_result$final_violation_result$parameters
    report_lines <- c(report_lines,
      "FINAL VALIDATED PARAMETERS:",
      sprintf("  q (overlap instances): %d", params$q),
      sprintf("  m (pathway instances): %d", params$m),
      sprintf("  n (non-pathway instances): %d", params$n),
      sprintf("  k (query instances): %d", params$k),
      ""
    )
  }

  # Recommendations
  report_lines <- c(report_lines,
    "RECOMMENDATIONS:",
    "- Review any applied corrections for biological validity",
    "- Consider data quality improvements to reduce future violations",
    "- Monitor extreme parameter values for numerical stability",
    ""
  )

  # Output report
  report_text <- paste(report_lines, collapse = "\n")

  if (!is.null(output_file)) {
    writeLines(report_lines, output_file)
    cat("Validation report saved to:", output_file, "\n")
  } else {
    cat(report_text)
  }

  return(invisible(report_lines))
}

# Testing function for constraint violation handling
test_constraint_violation_handling <- function() {
  cat("=== Constraint Violation Handling Tests ===\n")

  # Test 1: Recoverable violations
  cat("\nTest 1: Recoverable violations (negative copy numbers)\n")
  bad_query <- data.frame(gene = c("G1", "G2", "G3"), copy_number = c(-1, 2, 3))
  background <- data.frame(gene = c("G1", "G2", "G3", "G4"), copy_number = c(1, 2, 3, 4))
  pathway <- c("G1", "G2")

  recovery_result <- robust_constraint_validation(bad_query, pathway, background, recovery_mode = "lenient")
  cat("Recovery success:", recovery_result$success, "\n")

  # Test 2: Critical violations
  cat("\nTest 2: Critical violations (data inconsistency)\n")
  # Create impossible constraint violation - query with gene not in background
  impossible_query <- data.frame(gene = c("NONEXISTENT"), copy_number = c(5))
  background <- data.frame(gene = c("G1", "G2"), copy_number = c(1, 2))

  critical_result <- robust_constraint_validation(impossible_query, "G1", background, recovery_mode = "strict")
  cat("Critical error handling:", !critical_result$success, "\n")

  # Test 3: Generate report
  cat("\nTest 3: Report generation\n")
  report_file <- "test_validation_report.txt"
  generate_constraint_validation_report(recovery_result, report_file)
  cat("Report generated successfully\n")

  return(list(
    recoverable_test = recovery_result$success,
    critical_test = !critical_result$success,
    report_test = file.exists(report_file)
  ))
}

# Execute tests if run directly
if (!interactive()) {
  cat("CONSTRAINT VIOLATION HANDLER TESTING\n")
  cat("====================================\n")

  test_results <- test_constraint_violation_handling()

  overall_success <- all(unlist(test_results))
  cat("\n=== HANDLER TEST SUMMARY ===\n")
  cat("Overall result:", ifelse(overall_success, "ALL TESTS PASS", "SOME TESTS FAILED"), "\n")
}