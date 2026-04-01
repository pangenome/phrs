#!/usr/bin/env Rscript

# Computational Benchmarking of Copy-Number-Weighted Hypergeometric Testing
# Benchmarks different approaches to implementing weighted ORA using R phyper()

# Load required libraries
library(microbenchmark)
library(pryr)
library(ggplot2)
library(dplyr)
library(tibble)

# Set seed for reproducibility
set.seed(42)

# ==============================================================================
# APPROACH 1: DIRECT WEIGHTED PHYPER() WITH CALCULATED PARAMETERS
# ==============================================================================

weighted_phyper_direct <- function(query_df, pathway_genes, background_df) {
  # Calculate weighted parameters directly from copy numbers
  k_weighted <- sum(query_df$copy_number)                    # query instances

  query_in_pathway <- query_df[query_df$gene %in% pathway_genes, ]
  q_weighted <- sum(query_in_pathway$copy_number)            # overlap instances

  pathway_in_background <- background_df[background_df$gene %in% pathway_genes, ]
  m_weighted <- sum(pathway_in_background$copy_number)       # pathway instances

  n_weighted <- sum(background_df$copy_number) - m_weighted  # non-pathway instances

  # Parameter validation
  stopifnot(q_weighted <= k_weighted)
  stopifnot(q_weighted <= m_weighted)
  stopifnot(k_weighted <= (m_weighted + n_weighted))

  # Hypergeometric test
  pvalue <- phyper(q_weighted-1, m_weighted, n_weighted, k_weighted,
                   lower.tail = FALSE)

  list(
    pvalue = pvalue,
    params = list(q = q_weighted, m = m_weighted, n = n_weighted, k = k_weighted),
    method = "direct_weighted"
  )
}

# ==============================================================================
# APPROACH 2: INSTANCE EXPANSION + STANDARD PHYPER()
# ==============================================================================

expansion_phyper <- function(query_df, pathway_genes, background_df) {
  # Expand datasets by copy numbers
  query_expanded <- rep(query_df$gene, query_df$copy_number)
  background_expanded <- rep(background_df$gene, background_df$copy_number)

  # Calculate standard hypergeometric parameters on expanded data
  k_exp <- length(query_expanded)                            # expanded query size
  q_exp <- sum(query_expanded %in% pathway_genes)           # expanded overlap
  m_exp <- sum(background_expanded %in% pathway_genes)      # expanded pathway size
  n_exp <- length(background_expanded) - m_exp              # expanded non-pathway

  # Standard hypergeometric test
  pvalue <- phyper(q_exp-1, m_exp, n_exp, k_exp, lower.tail = FALSE)

  list(
    pvalue = pvalue,
    params = list(q = q_exp, m = m_exp, n = n_exp, k = k_exp),
    method = "instance_expansion"
  )
}

# ==============================================================================
# APPROACH 3: CUSTOM HYPERGEOMETRIC IMPLEMENTATION
# ==============================================================================

# Custom hypergeometric implementation with direct weighting
custom_hypergeometric <- function(query_df, pathway_genes, background_df) {
  # Calculate weighted parameters
  k_weighted <- sum(query_df$copy_number)
  query_in_pathway <- query_df[query_df$gene %in% pathway_genes, ]
  q_weighted <- sum(query_in_pathway$copy_number)
  pathway_in_background <- background_df[background_df$gene %in% pathway_genes, ]
  m_weighted <- sum(pathway_in_background$copy_number)
  n_weighted <- sum(background_df$copy_number) - m_weighted

  # Custom hypergeometric probability calculation
  # P(X >= q) = sum(dhyper(i, m, n, k)) for i = q to min(k, m)
  max_overlap <- min(k_weighted, m_weighted)

  if (q_weighted > max_overlap) {
    pvalue <- 0  # Impossible overlap
  } else {
    # Calculate probability mass for each possible overlap >= q_weighted
    probabilities <- sapply(q_weighted:max_overlap, function(i) {
      dhyper(i, m_weighted, n_weighted, k_weighted)
    })
    pvalue <- sum(probabilities)
  }

  list(
    pvalue = pvalue,
    params = list(q = q_weighted, m = m_weighted, n = n_weighted, k = k_weighted),
    method = "custom_hypergeometric"
  )
}

# ==============================================================================
# APPROACH 4: VECTORIZED MULTIPLE PATHWAY TESTING
# ==============================================================================

vectorized_multiple_pathways <- function(query_df, pathway_list, background_df) {
  # Pre-calculate query totals once
  k_weighted <- sum(query_df$copy_number)
  background_total <- sum(background_df$copy_number)

  # Vectorized calculation for multiple pathways
  results <- lapply(names(pathway_list), function(pathway_name) {
    pathway_genes <- pathway_list[[pathway_name]]

    # Calculate weighted parameters
    query_in_pathway <- query_df[query_df$gene %in% pathway_genes, ]
    q_weighted <- sum(query_in_pathway$copy_number)

    pathway_in_background <- background_df[background_df$gene %in% pathway_genes, ]
    m_weighted <- sum(pathway_in_background$copy_number)
    n_weighted <- background_total - m_weighted

    # Skip invalid parameter combinations
    if (q_weighted > k_weighted || q_weighted > m_weighted ||
        k_weighted > (m_weighted + n_weighted)) {
      return(list(pathway = pathway_name, pvalue = NA, method = "vectorized"))
    }

    pvalue <- phyper(q_weighted-1, m_weighted, n_weighted, k_weighted,
                     lower.tail = FALSE)

    list(
      pathway = pathway_name,
      pvalue = pvalue,
      params = list(q = q_weighted, m = m_weighted, n = n_weighted, k = k_weighted),
      method = "vectorized"
    )
  })

  return(results)
}

# ==============================================================================
# DATA GENERATION FUNCTIONS
# ==============================================================================

generate_copy_number_data <- function(n_genes, mean_copies = 3, max_copies = 50) {
  # Generate realistic copy number distribution
  # Use gamma distribution to model copy number variation
  copy_numbers <- round(rgamma(n_genes, shape = 2, scale = mean_copies/2))
  copy_numbers[copy_numbers == 0] <- 1  # Ensure minimum 1 copy
  copy_numbers[copy_numbers > max_copies] <- max_copies  # Cap extreme values

  data.frame(
    gene = paste0("GENE_", sprintf("%06d", 1:n_genes)),
    copy_number = copy_numbers,
    stringsAsFactors = FALSE
  )
}

generate_pathway_sets <- function(background_genes, n_pathways = 100,
                                pathway_size_range = c(10, 500)) {
  # Generate realistic pathway gene sets
  pathway_list <- list()

  for (i in 1:n_pathways) {
    pathway_size <- sample(pathway_size_range[1]:pathway_size_range[2], 1)
    pathway_genes <- sample(background_genes, pathway_size)
    pathway_list[[paste0("PATHWAY_", sprintf("%03d", i))]] <- pathway_genes
  }

  return(pathway_list)
}

create_benchmark_datasets <- function() {
  # Create datasets of different sizes for benchmarking

  # Small dataset (PHR-scale)
  small_bg <- generate_copy_number_data(1000, mean_copies = 2)
  small_query <- small_bg[sample(nrow(small_bg), 35), ]
  small_pathways <- generate_pathway_sets(small_bg$gene, n_pathways = 50,
                                        pathway_size_range = c(10, 100))

  # Medium dataset
  medium_bg <- generate_copy_number_data(5000, mean_copies = 3)
  medium_query <- medium_bg[sample(nrow(medium_bg), 500), ]
  medium_pathways <- generate_pathway_sets(medium_bg$gene, n_pathways = 200,
                                         pathway_size_range = c(20, 300))

  # Large dataset
  large_bg <- generate_copy_number_data(20000, mean_copies = 4)
  large_query <- large_bg[sample(nrow(large_bg), 2000), ]
  large_pathways <- generate_pathway_sets(large_bg$gene, n_pathways = 1000,
                                        pathway_size_range = c(50, 1000))

  list(
    small = list(background = small_bg, query = small_query, pathways = small_pathways),
    medium = list(background = medium_bg, query = medium_query, pathways = medium_pathways),
    large = list(background = large_bg, query = large_query, pathways = large_pathways)
  )
}

# ==============================================================================
# BENCHMARKING FUNCTIONS
# ==============================================================================

benchmark_single_pathway <- function(query_df, pathway_genes, background_df) {
  # Benchmark all approaches on a single pathway
  cat("Benchmarking dataset with", nrow(query_df), "query genes and",
      nrow(background_df), "background genes\n")

  # Memory measurement function
  measure_memory <- function(expr) {
    gc() # Clean up before measurement
    mem_before <- pryr::mem_used()
    result <- expr
    mem_after <- pryr::mem_used()
    mem_used <- as.numeric(mem_after - mem_before)
    list(result = result, memory_bytes = mem_used)
  }

  # Benchmark each approach
  cat("Running microbenchmark...\n")
  timing_results <- microbenchmark(
    direct = weighted_phyper_direct(query_df, pathway_genes, background_df),
    expansion = expansion_phyper(query_df, pathway_genes, background_df),
    custom = custom_hypergeometric(query_df, pathway_genes, background_df),
    times = 100,
    unit = "ms"
  )

  # Memory measurements (single run each due to overhead)
  cat("Measuring memory usage...\n")
  mem_direct <- measure_memory(
    weighted_phyper_direct(query_df, pathway_genes, background_df)
  )
  mem_expansion <- measure_memory(
    expansion_phyper(query_df, pathway_genes, background_df)
  )
  mem_custom <- measure_memory(
    custom_hypergeometric(query_df, pathway_genes, background_df)
  )

  # Verify mathematical equivalence
  cat("Verifying mathematical equivalence...\n")
  direct_result <- weighted_phyper_direct(query_df, pathway_genes, background_df)
  expansion_result <- expansion_phyper(query_df, pathway_genes, background_df)
  custom_result <- custom_hypergeometric(query_df, pathway_genes, background_df)

  pvalue_tolerance <- 1e-12
  params_equivalent <- all(
    direct_result$params$q == expansion_result$params$q,
    direct_result$params$m == expansion_result$params$m,
    direct_result$params$n == expansion_result$params$n,
    direct_result$params$k == expansion_result$params$k
  )

  pvalues_equivalent <- abs(direct_result$pvalue - expansion_result$pvalue) < pvalue_tolerance
  custom_equivalent <- abs(direct_result$pvalue - custom_result$pvalue) < pvalue_tolerance

  # Compile results
  list(
    timing = timing_results,
    memory = list(
      direct = mem_direct$memory_bytes,
      expansion = mem_expansion$memory_bytes,
      custom = mem_custom$memory_bytes
    ),
    equivalence = list(
      parameters_match = params_equivalent,
      direct_expansion_pvals_match = pvalues_equivalent,
      direct_custom_pvals_match = custom_equivalent,
      pvalue_differences = list(
        direct_vs_expansion = abs(direct_result$pvalue - expansion_result$pvalue),
        direct_vs_custom = abs(direct_result$pvalue - custom_result$pvalue)
      )
    ),
    dataset_info = list(
      query_genes = nrow(query_df),
      background_genes = nrow(background_df),
      total_query_instances = sum(query_df$copy_number),
      total_background_instances = sum(background_df$copy_number)
    )
  )
}

benchmark_multiple_pathways <- function(query_df, pathway_list, background_df, n_pathways = 10) {
  # Benchmark vectorized approach vs sequential single-pathway tests
  cat("Benchmarking multiple pathway testing with", n_pathways, "pathways\n")

  # Select subset of pathways for testing
  test_pathways <- pathway_list[1:min(n_pathways, length(pathway_list))]

  # Sequential approach (baseline)
  sequential_time <- system.time({
    sequential_results <- lapply(test_pathways, function(pathway_genes) {
      weighted_phyper_direct(query_df, pathway_genes, background_df)
    })
  })

  # Vectorized approach
  vectorized_time <- system.time({
    vectorized_results <- vectorized_multiple_pathways(query_df, test_pathways, background_df)
  })

  # Extract p-values for comparison
  sequential_pvals <- sapply(sequential_results, function(x) x$pvalue)
  vectorized_pvals <- sapply(vectorized_results, function(x) x$pvalue)

  # Check equivalence
  pval_differences <- abs(sequential_pvals - vectorized_pvals)
  max_difference <- max(pval_differences, na.rm = TRUE)

  list(
    timing = list(
      sequential_elapsed = sequential_time[["elapsed"]],
      vectorized_elapsed = vectorized_time[["elapsed"]],
      speedup = sequential_time[["elapsed"]] / vectorized_time[["elapsed"]]
    ),
    equivalence = list(
      max_pvalue_difference = max_difference,
      pvalues_equivalent = max_difference < 1e-12
    ),
    n_pathways = n_pathways
  )
}

# ==============================================================================
# COMPREHENSIVE BENCHMARK EXECUTION
# ==============================================================================

run_comprehensive_benchmark <- function() {
  cat("=== R phyper() Computational Efficiency Benchmark ===\n\n")

  # Generate test datasets
  cat("Generating benchmark datasets...\n")
  datasets <- create_benchmark_datasets()

  # Initialize results storage
  all_results <- list()

  # Benchmark each dataset size
  for (size_name in names(datasets)) {
    cat("\n=== Benchmarking", toupper(size_name), "dataset ===\n")

    data <- datasets[[size_name]]

    # Single pathway benchmark (first pathway)
    first_pathway <- data$pathways[[1]]
    single_results <- benchmark_single_pathway(
      data$query, first_pathway, data$background
    )

    # Multiple pathway benchmark
    multi_results <- benchmark_multiple_pathways(
      data$query, data$pathways, data$background, n_pathways = 10
    )

    all_results[[size_name]] <- list(
      single_pathway = single_results,
      multiple_pathways = multi_results
    )
  }

  return(all_results)
}

# ==============================================================================
# RESULTS ANALYSIS AND REPORTING
# ==============================================================================

analyze_benchmark_results <- function(results) {
  cat("\n=== Benchmark Results Analysis ===\n\n")

  # Create summary table
  summary_data <- data.frame()

  for (size_name in names(results)) {
    single_res <- results[[size_name]]$single_pathway
    multi_res <- results[[size_name]]$multiple_pathways

    # Extract timing statistics
    timing_stats <- summary(single_res$timing)

    # Add row to summary
    summary_data <- rbind(summary_data, data.frame(
      Dataset = size_name,
      Query_Genes = single_res$dataset_info$query_genes,
      Background_Genes = single_res$dataset_info$background_genes,
      Query_Instances = single_res$dataset_info$total_query_instances,
      Background_Instances = single_res$dataset_info$total_background_instances,

      # Timing (median milliseconds)
      Direct_Time_ms = median(timing_stats$time[timing_stats$expr == "direct"]) / 1e6,
      Expansion_Time_ms = median(timing_stats$time[timing_stats$expr == "expansion"]) / 1e6,
      Custom_Time_ms = median(timing_stats$time[timing_stats$expr == "custom"]) / 1e6,

      # Memory usage (bytes)
      Direct_Memory_bytes = single_res$memory$direct,
      Expansion_Memory_bytes = single_res$memory$expansion,
      Custom_Memory_bytes = single_res$memory$custom,

      # Performance ratios
      Time_Speedup_Direct_vs_Expansion =
        median(timing_stats$time[timing_stats$expr == "expansion"]) /
        median(timing_stats$time[timing_stats$expr == "direct"]),
      Memory_Reduction_Direct_vs_Expansion =
        single_res$memory$expansion / single_res$memory$direct,

      # Multiple pathway speedup
      Multi_Pathway_Speedup = multi_res$timing$speedup,

      # Equivalence checks
      Parameters_Equivalent = single_res$equivalence$parameters_match,
      PValues_Equivalent = single_res$equivalence$direct_expansion_pvals_match,

      stringsAsFactors = FALSE
    ))
  }

  return(summary_data)
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

if (!interactive()) {
  # Run benchmarks when script is executed directly
  cat("Starting comprehensive phyper() computational benchmark...\n")

  # Ensure required packages are available
  required_packages <- c("microbenchmark", "pryr", "ggplot2", "dplyr", "tibble")
  missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

  if (length(missing_packages) > 0) {
    cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages, repos = "https://cran.r-project.org")
  }

  # Load libraries
  for (pkg in required_packages) {
    library(pkg, character.only = TRUE)
  }

  # Run benchmark
  results <- run_comprehensive_benchmark()

  # Analyze results
  summary_table <- analyze_benchmark_results(results)

  # Save results
  saveRDS(results, "phyper_benchmark_results.rds")
  write.csv(summary_table, "phyper_benchmark_summary.csv", row.names = FALSE)

  # Display summary
  cat("\n=== PERFORMANCE SUMMARY ===\n")
  print(summary_table)

  cat("\nDetailed results saved to:\n")
  cat("- phyper_benchmark_results.rds (full results)\n")
  cat("- phyper_benchmark_summary.csv (summary table)\n")
}