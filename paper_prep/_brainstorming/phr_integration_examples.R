# PHR Integration Examples for Copy-Number Weighted Enrichment
# Demonstrates integration with PHR analysis pipeline

source("robust_copy_weighted_enrichment.R")
source("enrichment_validation_tests.R")

#' Example 1: Basic PHR enrichment analysis
#'
#' @param phr_genes Character vector of genes in PHRs
#' @param copy_number_data Data frame with gene copy numbers
#' @param pathway_database Named list of pathway gene sets
#'
#' @return Enrichment results for all pathways
example_phr_basic_enrichment <- function(phr_genes, copy_number_data, pathway_database) {

  cat("Running basic PHR enrichment analysis...\n")

  # Create query data frame from PHR genes
  query_df <- copy_number_data[copy_number_data$gene %in% phr_genes, ]

  if (nrow(query_df) == 0) {
    stop("No PHR genes found in copy number data")
  }

  cat("PHR query set:", nrow(query_df), "genes,", sum(query_df$copy_number), "total instances\n")

  # Run enrichment for each pathway
  results <- list()

  for (pathway_name in names(pathway_database)) {
    pathway_genes <- pathway_database[[pathway_name]]

    tryCatch({
      result <- weighted_hypergeometric_test(query_df, pathway_genes, copy_number_data)

      results[[pathway_name]] <- list(
        pathway = pathway_name,
        p_value = result$p_value,
        fold_enrichment = result$fold_enrichment,
        overlap_instances = result$overlap_instances,
        overlap_genes = sum(query_df$gene %in% pathway_genes),
        pathway_size = length(pathway_genes),
        status = result$status
      )

    }, error = function(e) {
      results[[pathway_name]] <- list(
        pathway = pathway_name,
        p_value = NA,
        fold_enrichment = NA,
        status = "error",
        error_message = as.character(e)
      )
    })
  }

  # Convert to data frame and sort by p-value
  results_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      pathway = x$pathway,
      p_value = ifelse(is.null(x$p_value), NA, x$p_value),
      fold_enrichment = ifelse(is.null(x$fold_enrichment), NA, x$fold_enrichment),
      overlap_instances = ifelse(is.null(x$overlap_instances), NA, x$overlap_instances),
      overlap_genes = ifelse(is.null(x$overlap_genes), NA, x$overlap_genes),
      pathway_size = ifelse(is.null(x$pathway_size), NA, x$pathway_size),
      status = x$status,
      stringsAsFactors = FALSE
    )
  }))

  # Sort by p-value
  results_df <- results_df[order(results_df$p_value, na.last = TRUE), ]

  # Add FDR correction
  valid_pvals <- !is.na(results_df$p_value)
  results_df$fdr <- NA
  if (sum(valid_pvals) > 0) {
    results_df$fdr[valid_pvals] <- p.adjust(results_df$p_value[valid_pvals], method = "fdr")
  }

  cat("Enrichment analysis complete. Found", sum(results_df$p_value < 0.05, na.rm = TRUE), "nominally significant pathways.\n")

  return(results_df)
}

#' Example 2: Compare copy-weighted vs standard enrichment
#'
#' @param phr_genes Character vector of genes in PHRs
#' @param copy_number_data Data frame with gene copy numbers
#' @param pathway_genes Character vector of pathway genes
#'
#' @return Comparison results
example_compare_enrichment_methods <- function(phr_genes, copy_number_data, pathway_genes) {

  cat("Comparing copy-weighted vs standard enrichment...\n")

  # Prepare data
  query_df <- copy_number_data[copy_number_data$gene %in% phr_genes, ]

  # Copy-weighted enrichment
  weighted_result <- weighted_hypergeometric_test(query_df, pathway_genes, copy_number_data)

  # Standard enrichment (treat all genes as copy number 1)
  standard_query_df <- data.frame(
    gene = unique(query_df$gene),
    copy_number = 1
  )
  standard_background_df <- data.frame(
    gene = unique(copy_number_data$gene),
    copy_number = 1
  )

  standard_result <- weighted_hypergeometric_test(standard_query_df, pathway_genes, standard_background_df)

  # Create comparison
  comparison <- data.frame(
    method = c("Copy-weighted", "Standard"),
    p_value = c(weighted_result$p_value, standard_result$p_value),
    fold_enrichment = c(weighted_result$fold_enrichment, standard_result$fold_enrichment),
    overlap_count = c(weighted_result$overlap_instances, standard_result$overlap_instances),
    query_count = c(weighted_result$query_instances, standard_result$query_instances),
    stringsAsFactors = FALSE
  )

  cat("Results:\n")
  print(comparison)

  # Calculate effect of copy weighting
  p_value_ratio <- weighted_result$p_value / standard_result$p_value
  enrichment_ratio <- weighted_result$fold_enrichment / standard_result$fold_enrichment

  cat("\nCopy weighting effects:\n")
  cat("  P-value ratio (weighted/standard):", round(p_value_ratio, 4), "\n")
  cat("  Fold enrichment ratio (weighted/standard):", round(enrichment_ratio, 4), "\n")

  return(list(
    comparison = comparison,
    p_value_ratio = p_value_ratio,
    enrichment_ratio = enrichment_ratio,
    weighted_result = weighted_result,
    standard_result = standard_result
  ))
}

#' Example 3: PHR-specific pathway analysis with copy number stratification
#'
#' @param phr_genes Character vector of genes in PHRs
#' @param copy_number_data Data frame with gene copy numbers
#' @param pathway_genes Character vector of pathway genes
#'
#' @return Stratified analysis results
example_copy_number_stratification <- function(phr_genes, copy_number_data, pathway_genes) {

  cat("Running copy number stratification analysis...\n")

  # Prepare PHR query data
  query_df <- copy_number_data[copy_number_data$gene %in% phr_genes, ]

  # Stratify by copy number levels
  copy_strata <- list(
    "low_copy" = query_df[query_df$copy_number == 1, ],
    "medium_copy" = query_df[query_df$copy_number == 2, ],
    "high_copy" = query_df[query_df$copy_number >= 3, ]
  )

  results <- list()

  for (stratum_name in names(copy_strata)) {
    stratum_df <- copy_strata[[stratum_name]]

    if (nrow(stratum_df) > 0) {
      cat("Testing", stratum_name, "stratum:", nrow(stratum_df), "genes\n")

      result <- weighted_hypergeometric_test(stratum_df, pathway_genes, copy_number_data)

      results[[stratum_name]] <- list(
        stratum = stratum_name,
        gene_count = nrow(stratum_df),
        instance_count = sum(stratum_df$copy_number),
        p_value = result$p_value,
        fold_enrichment = result$fold_enrichment,
        overlap_instances = result$overlap_instances
      )
    } else {
      results[[stratum_name]] <- list(
        stratum = stratum_name,
        gene_count = 0,
        instance_count = 0,
        p_value = NA,
        fold_enrichment = NA,
        overlap_instances = 0
      )
    }
  }

  # Create results table
  results_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      stratum = x$stratum,
      gene_count = x$gene_count,
      instance_count = x$instance_count,
      p_value = ifelse(is.null(x$p_value), NA, x$p_value),
      fold_enrichment = ifelse(is.null(x$fold_enrichment), NA, x$fold_enrichment),
      overlap_instances = x$overlap_instances,
      stringsAsFactors = FALSE
    )
  }))

  cat("\nStratification results:\n")
  print(results_df)

  return(results_df)
}

#' Example 4: Performance demonstration with PHR data
#'
#' @param phr_genes Character vector of genes in PHRs
#' @param copy_number_data Data frame with gene copy numbers
#' @param pathway_genes Character vector of pathway genes
#'
#' @return Performance comparison results
example_performance_demonstration <- function(phr_genes, copy_number_data, pathway_genes) {

  cat("Running performance demonstration...\n")

  # Prepare PHR query data
  query_df <- copy_number_data[copy_number_data$gene %in% phr_genes, ]

  # Run benchmark
  benchmark_result <- benchmark_enrichment_methods(query_df, pathway_genes, copy_number_data, n_reps = 50)

  cat("Performance Results:\n")
  cat("  Parameter method mean time:", round(benchmark_result$parameter_method$mean_time * 1000, 2), "ms\n")
  cat("  Expansion method mean time:", round(benchmark_result$expansion_method$mean_time * 1000, 2), "ms\n")
  cat("  Speedup factor:", round(benchmark_result$speedup_factor, 2), "x\n")
  cat("  Memory reduction:", round(benchmark_result$memory_reduction * 100, 1), "%\n")

  # Memory usage
  cat("\nMemory Usage:\n")
  cat("  Parameter method estimate:", benchmark_result$parameter_method$memory_estimate, "objects\n")
  cat("  Expansion method estimate:", benchmark_result$expansion_method$memory_estimate, "objects\n")

  return(benchmark_result)
}

#' Example 5: Complete PHR analysis workflow
#'
#' @param phr_bed_file Path to BED file with PHR coordinates
#' @param gene_copy_file Path to CSV file with gene copy numbers
#' @param pathway_file Path to GMT file with pathway definitions
#'
#' @return Complete analysis results
example_complete_phr_workflow <- function(phr_bed_file = "chm13.phrs.bed",
                                        gene_copy_file = "gene_copy_summary.csv",
                                        pathway_file = NULL) {

  cat("Running complete PHR enrichment workflow...\n")

  # Step 1: Load PHR data
  if (file.exists(phr_bed_file)) {
    cat("Loading PHR regions from", phr_bed_file, "...\n")
    # Note: This would require additional processing to extract genes from regions
    # For demonstration, using mock data
  }

  # Step 2: Load gene copy data
  if (file.exists(gene_copy_file)) {
    cat("Loading gene copy data from", gene_copy_file, "...\n")
    copy_data <- tryCatch({
      read.csv(gene_copy_file, stringsAsFactors = FALSE)
    }, error = function(e) {
      cat("Could not load gene copy data, using mock data\n")
      create_mock_copy_data()
    })
  } else {
    cat("Gene copy file not found, using mock data\n")
    copy_data <- create_mock_copy_data()
  }

  # Step 3: Define pathways (mock for demonstration)
  pathways <- create_mock_pathways()

  # Step 4: Define PHR genes (mock for demonstration)
  phr_genes <- sample(copy_data$gene, min(30, nrow(copy_data)))

  # Step 5: Run validation
  cat("Running validation tests...\n")
  validation_passed <- quick_validation()

  if (!validation_passed) {
    warning("Validation tests failed - results may be unreliable")
  }

  # Step 6: Run enrichment analysis
  cat("Running enrichment analysis...\n")
  enrichment_results <- example_phr_basic_enrichment(phr_genes, copy_data, pathways)

  # Step 7: Performance benchmark
  if (nrow(enrichment_results) > 0) {
    top_pathway <- pathways[[enrichment_results$pathway[1]]]
    performance_results <- example_performance_demonstration(phr_genes, copy_data, top_pathway)
  }

  # Step 8: Generate summary report
  summary_report <- list(
    phr_gene_count = length(phr_genes),
    total_copy_instances = sum(copy_data[copy_data$gene %in% phr_genes, "copy_number"]),
    pathways_tested = nrow(enrichment_results),
    significant_pathways = sum(enrichment_results$p_value < 0.05, na.rm = TRUE),
    validation_passed = validation_passed,
    top_result = enrichment_results[1, ]
  )

  cat("\n=== WORKFLOW SUMMARY ===\n")
  cat("PHR genes:", summary_report$phr_gene_count, "\n")
  cat("Total copy instances:", summary_report$total_copy_instances, "\n")
  cat("Pathways tested:", summary_report$pathways_tested, "\n")
  cat("Significant pathways (p < 0.05):", summary_report$significant_pathways, "\n")
  cat("Validation passed:", summary_report$validation_passed, "\n")

  if (summary_report$significant_pathways > 0) {
    cat("\nTop result:\n")
    top <- summary_report$top_result
    cat("  Pathway:", top$pathway, "\n")
    cat("  P-value:", scientific_notation(top$p_value), "\n")
    cat("  Fold enrichment:", round(top$fold_enrichment, 2), "\n")
  }

  return(list(
    summary = summary_report,
    enrichment_results = enrichment_results,
    copy_data = copy_data,
    pathways = pathways
  ))
}

# Helper functions for mock data generation
create_mock_copy_data <- function() {
  n_genes <- 100
  genes <- paste0("GENE", sprintf("%04d", 1:n_genes))
  copy_numbers <- sample(1:5, n_genes, replace = TRUE, prob = c(0.4, 0.3, 0.15, 0.1, 0.05))

  data.frame(
    gene = genes,
    copy_number = copy_numbers,
    stringsAsFactors = FALSE
  )
}

create_mock_pathways <- function() {
  all_genes <- paste0("GENE", sprintf("%04d", 1:100))

  list(
    "DNA_REPAIR" = sample(all_genes, 15),
    "CELL_CYCLE" = sample(all_genes, 20),
    "APOPTOSIS" = sample(all_genes, 12),
    "METABOLISM" = sample(all_genes, 25),
    "SIGNALING" = sample(all_genes, 18)
  )
}

# Helper function for scientific notation
scientific_notation <- function(x, digits = 2) {
  if (is.na(x)) return("NA")
  if (x == 0) return("0")
  if (x >= 0.01) return(round(x, 4))
  return(formatC(x, format = "e", digits = digits))
}

#' Quick integration test
test_phr_integration <- function() {
  cat("Testing PHR integration...\n")

  # Create test data
  copy_data <- create_mock_copy_data()
  pathways <- create_mock_pathways()
  phr_genes <- sample(copy_data$gene, 20)

  # Test basic functionality
  tryCatch({
    results <- example_phr_basic_enrichment(phr_genes, copy_data, pathways)

    if (is.data.frame(results) && nrow(results) > 0) {
      cat("PHR integration test: PASS\n")
      return(TRUE)
    } else {
      cat("PHR integration test: FAIL - no results\n")
      return(FALSE)
    }

  }, error = function(e) {
    cat("PHR integration test: FAIL -", e$message, "\n")
    return(FALSE)
  })
}