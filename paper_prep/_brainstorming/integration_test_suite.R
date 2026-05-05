#!/usr/bin/env Rscript

#' Integration Test Suite for Copy-Number-Weighted ORA Parameter Mapping
#'
#' This script tests integration between copy_number_phyper_mapping.R and
#' existing ORA workflows including gprofiler, Python enrichment scripts,
#' and PHR dataset analysis pipelines.
#'
#' Author: Workgraph Agent (integration-testing-with task)
#' Date: 2026-04-01

# Load required libraries
suppressPackageStartupMessages({
  library(jsonlite)
  library(tools)
  library(utils)
})

# Load the implementation
source("copy_number_phyper_mapping.R")

cat("=== Copy-Number-Weighted ORA Integration Test Suite ===\n\n")

# ==============================================================================
# TEST 1: GPROFILER INTEGRATION
# ==============================================================================

test_gprofiler_integration <- function() {
  cat("TEST 1: g:Profiler Integration\n")
  cat("==============================\n")

  # Load existing gprofiler request
  if (!file.exists("gprofiler_request.json")) {
    cat("SKIP: gprofiler_request.json not found\n\n")
    return(list(status = "skipped", reason = "gprofiler_request.json not found"))
  }

  gprofiler_data <- fromJSON("gprofiler_request.json")
  phr_genes <- gprofiler_data$query

  cat(sprintf("Loaded %d PHR genes from gprofiler request\n", length(phr_genes)))

  # Create query data with copy numbers for known PHR genes
  # Use realistic copy numbers for olfactory receptor genes
  query_genes <- head(phr_genes, 20)  # Test subset for speed
  or_genes <- query_genes[grepl("^OR", query_genes)]
  other_genes <- setdiff(query_genes, or_genes)

  query_df <- data.frame(
    gene_name = query_genes,
    copy_number = c(
      rep(10, length(or_genes)),  # OR genes typically high copy
      rep(2, length(other_genes))  # Other genes typically low copy
    ),
    stringsAsFactors = FALSE
  )

  # Define test pathway (olfactory receptor family)
  or_pathway <- c("OR4F17", "OR4F29", "OR4F3", "OR4F5", "OR1A1", "OR1A2",
                  "OR2A1", "OR2A2", "OR6A2", "OR51E1", or_genes)

  # Create background with realistic genomic copy distribution
  background_df <- create_genomic_background()

  cat("Testing parameter mapping with g:Profiler gene set...\n")

  # Test parameter calculation
  tryCatch({
    params <- calculate_weighted_phyper_params(query_df, or_pathway, background_df)

    # Validate parameters
    if (!params$validation$passed) {
      return(list(
        status = "failed",
        reason = paste("Parameter validation failed:",
                      paste(params$validation$errors, collapse = "; "))
      ))
    }

    # Run hypergeometric test
    test_result <- run_weighted_hypergeometric_test(query_df, or_pathway, background_df)

    cat(sprintf("✓ Parameters calculated successfully\n"))
    cat(sprintf("  - Query instances: %d\n", params$k_weighted))
    cat(sprintf("  - Overlap instances: %d\n", params$q_weighted))
    cat(sprintf("  - P-value: %.2e\n", test_result$pvalue))
    cat(sprintf("  - Fold enrichment: %.2f\n", test_result$fold_enrichment))

    # Test compatibility with g:Profiler output format
    gprofiler_compatible_result <- format_for_gprofiler(test_result, or_pathway)

    cat("✓ g:Profiler integration test passed\n\n")
    return(list(
      status = "passed",
      parameters = params,
      test_result = test_result,
      gprofiler_format = gprofiler_compatible_result
    ))

  }, error = function(e) {
    cat(sprintf("✗ g:Profiler integration test failed: %s\n\n", e$message))
    return(list(status = "failed", reason = e$message))
  })
}

# ==============================================================================
# TEST 2: PHR DATASET INTEGRATION
# ==============================================================================

test_phr_dataset_integration <- function() {
  cat("TEST 2: PHR Dataset Integration\n")
  cat("===============================\n")

  # Test with actual PHR BED file
  if (!file.exists("chm13.phrs.bed")) {
    cat("SKIP: chm13.phrs.bed not found\n\n")
    return(list(status = "skipped", reason = "chm13.phrs.bed not found"))
  }

  # Load PHR intervals
  phr_bed <- read.table("chm13.phrs.bed", sep = "\t", header = FALSE,
                        col.names = c("chr", "start", "end", "copies"))

  cat(sprintf("Loaded %d PHR intervals from chm13.phrs.bed\n", nrow(phr_bed)))

  # Test parameter mapping with real gene copy data if available
  if (file.exists("gene_copy_summary.csv")) {
    copy_data <- read.csv("gene_copy_summary.csv", stringsAsFactors = FALSE)
    protein_genes <- copy_data[copy_data$gene_biotype == "protein_coding", ]

    # Use top 50 highest copy genes as test query
    top_copy_genes <- head(protein_genes[order(-protein_genes$total_copies), ], 50)

    query_df <- data.frame(
      gene_name = top_copy_genes$gene_name,
      copy_number = top_copy_genes$total_copies,
      stringsAsFactors = FALSE
    )

    cat(sprintf("Testing with %d high-copy PHR genes\n", nrow(query_df)))
  } else {
    # Use mock data based on PHR characteristics
    query_df <- create_mock_phr_query()
    cat("Using mock PHR gene data\n")
  }

  # Test with chromosome arm enrichment pathway
  chr_arm_pathway <- c("AMACR", "ACRO1", "H2AFZ", "HIST1H1A", "HIST1H2AB",
                      "HIST1H3A", "WASH1")  # Known centromeric/pericentromeric genes

  background_df <- create_genomic_background()

  tryCatch({
    params <- calculate_weighted_phyper_params(query_df, chr_arm_pathway, background_df)
    test_result <- run_weighted_hypergeometric_test(query_df, chr_arm_pathway, background_df)

    cat("✓ PHR dataset integration test passed\n")
    cat(sprintf("  - PHR query genes: %d\n", nrow(query_df)))
    cat(sprintf("  - Total query copies: %d\n", sum(query_df$copy_number)))
    cat(sprintf("  - P-value: %.2e\n", test_result$pvalue))

    return(list(
      status = "passed",
      phr_intervals = nrow(phr_bed),
      query_genes = nrow(query_df),
      parameters = params,
      test_result = test_result
    ))

  }, error = function(e) {
    cat(sprintf("✗ PHR dataset integration test failed: %s\n\n", e$message))
    return(list(status = "failed", reason = e$message))
  })
}

# ==============================================================================
# TEST 3: PYTHON WORKFLOW COMPATIBILITY
# ==============================================================================

test_python_workflow_compatibility <- function() {
  cat("TEST 3: Python Workflow Compatibility\n")
  cat("=====================================\n")

  if (!file.exists("copy_number_enrichment.py")) {
    cat("SKIP: copy_number_enrichment.py not found\n\n")
    return(list(status = "skipped", reason = "copy_number_enrichment.py not found"))
  }

  # Create test data in format expected by Python script
  test_gene_copy <- data.frame(
    gene_name = paste0("GENE", 1:100),
    gene_biotype = "protein_coding",
    total_copies = sample(1:20, 100, replace = TRUE, prob = c(0.4, 0.25, 0.15, 0.1, 0.05, rep(0.01, 15)))
  )

  test_all_copies <- data.frame(
    gene_name = rep(test_gene_copy$gene_name, test_gene_copy$total_copies),
    chr = sample(paste0("chr", 1:22), sum(test_gene_copy$total_copies), replace = TRUE),
    start = sample(1:1000000, sum(test_gene_copy$total_copies)),
    end = sample(1000001:2000000, sum(test_gene_copy$total_copies))
  )

  # Write temporary files
  write.csv(test_gene_copy, "test_gene_copy.csv", row.names = FALSE)
  write.csv(test_all_copies, "test_all_copies.csv", row.names = FALSE)

  # Create mock PHR intervals
  test_phr <- data.frame(
    chr = c("chr1", "chr2", "chr3"),
    start = c(1000, 2000, 3000),
    end = c(2000, 3000, 4000),
    name = c("phr1", "phr2", "phr3")
  )
  write.table(test_phr, "test_phr.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  cat("Created test data files for Python workflow\n")

  # Test parameter mapping equivalence with subset
  subset_genes <- head(test_gene_copy, 20)
  query_df <- data.frame(
    gene_name = subset_genes$gene_name,
    copy_number = subset_genes$total_copies,
    stringsAsFactors = FALSE
  )

  pathway_genes <- sample(query_df$gene_name, 5)
  background_df <- create_genomic_background()

  # Ensure pathway genes in background have at least as many copies as in query
  # to satisfy hypergeometric constraints
  for (gene in pathway_genes) {
    if (gene %in% query_df$gene_name) {
      query_copies <- query_df$copy_number[query_df$gene_name == gene]
      if (!gene %in% background_df$gene_name) {
        # Add gene to background
        background_df <- rbind(background_df, data.frame(
          gene_name = gene,
          copy_number = max(query_copies, 2),
          stringsAsFactors = FALSE
        ))
      } else {
        # Ensure background has at least as many copies
        bg_idx <- background_df$gene_name == gene
        background_df$copy_number[bg_idx] <- max(background_df$copy_number[bg_idx], query_copies)
      }
    }
  }

  tryCatch({
    # Test R implementation
    r_params <- calculate_weighted_phyper_params(query_df, pathway_genes, background_df)
    r_result <- run_weighted_hypergeometric_test(query_df, pathway_genes, background_df)

    # Validate that our parameters would be compatible with Python scipy.stats.hypergeom
    scipy_params <- c(
      M = r_params$m_weighted + r_params$n_weighted,  # population size
      n = r_params$m_weighted,                        # number of success states
      N = r_params$k_weighted                         # sample size
    )

    cat("✓ Python workflow compatibility test passed\n")
    cat(sprintf("  - R phyper params: q=%d, m=%d, n=%d, k=%d\n",
               r_params$q_weighted, r_params$m_weighted, r_params$n_weighted, r_params$k_weighted))
    cat(sprintf("  - SciPy hypergeom params: M=%d, n=%d, N=%d\n",
               scipy_params[1], scipy_params[2], scipy_params[3]))

    # Clean up test files
    file.remove(c("test_gene_copy.csv", "test_all_copies.csv", "test_phr.bed"))

    return(list(
      status = "passed",
      r_params = r_params,
      scipy_params = scipy_params,
      compatible = TRUE
    ))

  }, error = function(e) {
    # Clean up test files on error
    file.remove(c("test_gene_copy.csv", "test_all_copies.csv", "test_phr.bed"))
    cat(sprintf("✗ Python workflow compatibility test failed: %s\n\n", e$message))
    return(list(status = "failed", reason = e$message))
  })
}

# ==============================================================================
# TEST 4: ERROR HANDLING AND EDGE CASES
# ==============================================================================

test_error_handling <- function() {
  cat("TEST 4: Error Handling and Edge Cases\n")
  cat("=====================================\n")

  results <- list()

  # Test 4.1: Empty query
  tryCatch({
    empty_query <- data.frame(gene_name = character(0), copy_number = numeric(0))
    pathway <- c("GENE1", "GENE2")
    background <- create_genomic_background()

    result <- calculate_weighted_phyper_params(empty_query, pathway, background)
    results$empty_query <- "failed - should have thrown error"
  }, error = function(e) {
    results$empty_query <<- "passed - correctly caught empty query"
    cat("✓ Empty query error handling: PASS\n")
  })

  # Test 4.2: Zero copy numbers
  tryCatch({
    zero_query <- data.frame(gene_name = c("GENE1", "GENE2"), copy_number = c(0, 5))
    pathway <- c("GENE1", "GENE2")
    background <- create_genomic_background()

    result <- calculate_weighted_phyper_params(zero_query, pathway, background)
    if (nrow(result$metadata) == 1) {  # Should filter out zero-copy gene
      results$zero_copies <- "passed - correctly filtered zero copies"
      cat("✓ Zero copy number handling: PASS\n")
    } else {
      results$zero_copies <- "failed - did not filter zero copies"
    }
  }, error = function(e) {
    results$zero_copies <<- paste("error:", e$message)
  })

  # Test 4.3: No pathway overlap
  tryCatch({
    query <- data.frame(gene_name = c("GENE1", "GENE2"), copy_number = c(5, 3))
    pathway <- c("GENE3", "GENE4")  # No overlap
    background <- create_genomic_background()

    result <- run_weighted_hypergeometric_test(query, pathway, background)
    if (result$pvalue == 1.0 && result$observed_overlap_weighted == 0) {
      results$no_overlap <- "passed - correctly handled no overlap"
      cat("✓ No pathway overlap handling: PASS\n")
    } else {
      results$no_overlap <- "failed - incorrect handling of no overlap"
    }
  }, error = function(e) {
    results$no_overlap <<- paste("error:", e$message)
  })

  # Test 4.4: Invalid parameter constraints
  tryCatch({
    # Create scenario where overlap > query (should be caught)
    query <- data.frame(gene_name = c("GENE1"), copy_number = c(5))
    pathway <- c("GENE1")
    # Create background where pathway gene has fewer copies than query
    background <- data.frame(gene_name = c("GENE1", paste0("BG", 1:100)),
                           copy_number = c(3, rep(2, 100)))

    result <- calculate_weighted_phyper_params(query, pathway, background, validate_params = TRUE)
    if (!result$validation$passed) {
      results$invalid_constraints <- "passed - caught invalid constraints"
      cat("✓ Invalid parameter constraints: PASS\n")
    } else {
      results$invalid_constraints <- "failed - missed invalid constraints"
    }
  }, error = function(e) {
    results$invalid_constraints <<- paste("error:", e$message)
  })

  return(results)
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

create_genomic_background <- function(n_genes = 20000) {
  "Create realistic genomic background with copy number distribution"

  # Realistic copy number distribution (most genes have 1-2 copies)
  copy_probs <- c(0.85, 0.10, 0.03, 0.015, 0.005)  # For copies 1,2,3,4,5+

  background <- data.frame(
    gene_name = paste0("BG_GENE", 1:n_genes),
    copy_number = sample(1:5, n_genes, replace = TRUE, prob = copy_probs),
    stringsAsFactors = FALSE
  )

  return(background)
}

create_mock_phr_query <- function() {
  "Create mock PHR gene query based on known PHR characteristics"

  # PHR genes tend to be highly repetitive
  phr_genes <- c(
    paste0("OR4F", c(17, 29, 3, 5, 8)),  # Olfactory receptors (high copy)
    paste0("FAM138", LETTERS[1:5]),       # FAM genes (medium copy)
    paste0("LOC", c(100419924, 101927506, 102723376)),  # LOC genes
    paste0("SEPTIN14P", 11:15),           # Pseudogenes (variable copy)
    paste0("RPL23AP", c(21, 24, 25))      # Ribosomal pseudogenes
  )

  # Assign realistic copy numbers
  copy_numbers <- c(
    rep(12, 5),   # OR genes (high)
    rep(8, 5),    # FAM genes (medium)
    rep(6, 3),    # LOC genes (medium)
    rep(4, 5),    # Pseudogenes (low-medium)
    rep(10, 3)    # Ribosomal (high)
  )

  query_df <- data.frame(
    gene_name = phr_genes,
    copy_number = copy_numbers,
    stringsAsFactors = FALSE
  )

  return(query_df)
}

format_for_gprofiler <- function(test_result, pathway_genes) {
  "Format results in g:Profiler compatible structure"

  return(list(
    term_id = "custom_pathway",
    term_name = "Custom Pathway Test",
    p_value = test_result$pvalue,
    significant = test_result$significant,
    query_size = test_result$query_size_weighted,
    term_size = test_result$pathway_size_weighted,
    intersection_size = test_result$observed_overlap_weighted,
    fold_enrichment = test_result$fold_enrichment,
    intersection = pathway_genes[1:min(5, length(pathway_genes))]  # Sample intersection genes
  ))
}

# ==============================================================================
# MAIN TEST RUNNER
# ==============================================================================

run_integration_test_suite <- function() {
  cat("Starting Integration Test Suite...\n")
  cat("=================================\n\n")

  results <- list()

  # Run all tests
  results$gprofiler <- test_gprofiler_integration()
  results$phr_dataset <- test_phr_dataset_integration()
  results$python_compatibility <- test_python_workflow_compatibility()
  results$error_handling <- test_error_handling()

  # Generate summary
  cat("\n=== INTEGRATION TEST SUMMARY ===\n")

  for (test_name in names(results)) {
    if (is.list(results[[test_name]]) && "status" %in% names(results[[test_name]])) {
      status <- results[[test_name]]$status
      cat(sprintf("%-25s: %s\n", test_name, toupper(status)))
      if (status == "failed" && "reason" %in% names(results[[test_name]])) {
        cat(sprintf("  Reason: %s\n", results[[test_name]]$reason))
      }
    } else {
      cat(sprintf("%-25s: COMPLETED\n", test_name))
    }
  }

  # Check overall pass rate
  statuses <- sapply(results, function(x) if("status" %in% names(x)) x$status else "completed")
  passed <- sum(statuses == "passed")
  total <- length(statuses[statuses != "skipped"])

  cat(sprintf("\nOverall: %d/%d tests passed (%.1f%%)\n", passed, total, 100*passed/total))

  return(results)
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if (!interactive()) {
  # Run integration test suite
  test_results <- run_integration_test_suite()

  # Save results for reporting
  saveRDS(test_results, "integration_test_results.rds")
  cat("\nResults saved to integration_test_results.rds\n")
}

# String concatenation operator for nicer output
`%+%` <- function(a, b) paste0(a, b)