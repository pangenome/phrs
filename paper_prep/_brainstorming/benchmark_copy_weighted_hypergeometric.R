# Performance Benchmarking for Copy-Number Weighted Hypergeometric Testing
#
# Comprehensive benchmarking suite comparing parameter weighting vs instance
# expansion approaches across different data scales and copy number distributions.
#
# Author: Robust R Code Implementation Task
# Date: 2026-04-01
# Version: 1.0

# Source the implementation
source("copy_weighted_hypergeometric.R")

# Load required libraries
if (!require(microbenchmark, quietly = TRUE)) {
  cat("Installing microbenchmark package for precise timing...\n")
  install.packages("microbenchmark", quiet = TRUE)
  library(microbenchmark)
}

#' Instance Expansion Implementation (for comparison)
#'
#' Traditional approach that expands datasets by copy number
instance_expansion_test <- function(query_df, pathway_genes, background_df) {

  # Expand query and background by copy numbers
  query_expanded <- rep(query_df$gene, query_df$copy_number)
  background_expanded <- rep(background_df$gene, background_df$copy_number)

  # Calculate standard hypergeometric parameters on expanded data
  q_exp <- sum(query_expanded %in% pathway_genes)
  m_exp <- sum(background_expanded %in% pathway_genes)
  n_exp <- length(background_expanded) - m_exp
  k_exp <- length(query_expanded)

  # Hypergeometric test
  pvalue <- phyper(q_exp - 1, m_exp, n_exp, k_exp, lower.tail = FALSE)

  return(list(
    pvalue = pvalue,
    overlap_instances = q_exp,
    query_instances = k_exp,
    pathway_instances = m_exp,
    background_instances = m_exp + n_exp
  ))
}

#' Generate Test Datasets
#'
#' Creates synthetic datasets of varying scales for benchmarking
generate_test_data <- function(scale = "small") {

  scales <- list(
    small = list(
      query_genes = 10,
      background_genes = 1000,
      pathway_genes = 50,
      copy_lambda = 2
    ),
    medium = list(
      query_genes = 35,
      background_genes = 5000,
      pathway_genes = 200,
      copy_lambda = 5
    ),
    large = list(
      query_genes = 100,
      background_genes = 20000,
      pathway_genes = 1000,
      copy_lambda = 10
    ),
    phr_like = list(
      query_genes = 35,
      background_genes = 20000,
      pathway_genes = 400,
      copy_lambda = 30  # High copy numbers like PHRs
    )
  )

  if (!scale %in% names(scales)) {
    stop("Scale must be one of:", paste(names(scales), collapse = ", "))
  }

  params <- scales[[scale]]

  set.seed(42)  # For reproducible benchmarks

  # Generate query dataset
  query_genes <- paste0("QUERY", 1:params$query_genes)
  query_copies <- rpois(params$query_genes, lambda = params$copy_lambda) + 1

  query_df <- data.frame(
    gene = query_genes,
    copy_number = query_copies
  )

  # Generate background dataset
  background_genes <- paste0("GENE", 1:params$background_genes)
  background_copies <- rpois(params$background_genes, lambda = params$copy_lambda) + 1

  background_df <- data.frame(
    gene = background_genes,
    copy_number = background_copies
  )

  # Pathway genes (some overlap with query)
  pathway_genes <- c(
    paste0("QUERY", 1:min(5, params$query_genes)),  # Some query genes in pathway
    paste0("GENE", 1:params$pathway_genes)          # Additional pathway genes
  )

  return(list(
    query_df = query_df,
    background_df = background_df,
    pathway_genes = pathway_genes,
    scale_info = params
  ))
}

#' Memory Usage Estimation
#'
#' Estimates memory usage for different approaches
estimate_memory_usage <- function(query_df, background_df) {

  # Calculate instance counts
  query_instances <- sum(query_df$copy_number)
  background_instances <- sum(background_df$copy_number)

  # Instance expansion memory (character vectors)
  # Assuming average gene name is 8 characters (64 bits) + overhead
  expansion_memory_bytes <- (query_instances + background_instances) * 64

  # Parameter weighting memory (original dataframes)
  # gene name + copy number per row
  weighting_memory_bytes <- (nrow(query_df) + nrow(background_df)) * (64 + 8)

  return(list(
    expansion_mb = expansion_memory_bytes / (1024^2),
    weighting_mb = weighting_memory_bytes / (1024^2),
    memory_reduction_factor = expansion_memory_bytes / weighting_memory_bytes,
    query_instances = query_instances,
    background_instances = background_instances
  ))
}

#' Single Scale Benchmark
#'
#' Benchmarks both approaches on a single data scale
benchmark_single_scale <- function(scale, n_iterations = 10) {

  cat(paste("Benchmarking", scale, "scale with", n_iterations, "iterations...\n"))

  # Generate test data
  test_data <- generate_test_data(scale)

  # Memory usage estimation
  memory_info <- estimate_memory_usage(test_data$query_df, test_data$background_df)

  cat(paste("  Dataset:", nrow(test_data$query_df), "query genes,",
            nrow(test_data$background_df), "background genes\n"))
  cat(paste("  Instances:", memory_info$query_instances, "query,",
            memory_info$background_instances, "background\n"))
  cat(paste("  Memory reduction:", round(memory_info$memory_reduction_factor, 1), "x\n"))

  # Benchmark both approaches
  timing_results <- microbenchmark(
    parameter_weighting = {
      result <- weighted_hypergeometric_test(
        test_data$query_df, test_data$pathway_genes, test_data$background_df,
        validate_inputs = FALSE  # Skip validation for pure timing
      )
    },
    instance_expansion = {
      result <- instance_expansion_test(
        test_data$query_df, test_data$pathway_genes, test_data$background_df
      )
    },
    times = n_iterations,
    unit = "ms"
  )

  # Verify both methods give same results
  weighted_result <- weighted_hypergeometric_test(
    test_data$query_df, test_data$pathway_genes, test_data$background_df
  )

  expansion_result <- instance_expansion_test(
    test_data$query_df, test_data$pathway_genes, test_data$background_df
  )

  pvalue_diff <- abs(weighted_result$pvalue - expansion_result$pvalue)
  results_match <- pvalue_diff < 1e-12

  # Summarize timing
  timing_summary <- summary(timing_results)
  weighting_median <- timing_summary$median[timing_summary$expr == "parameter_weighting"]
  expansion_median <- timing_summary$median[timing_summary$expr == "instance_expansion"]
  speedup_factor <- expansion_median / weighting_median

  cat(paste("  Parameter weighting median:", round(weighting_median, 2), "ms\n"))
  cat(paste("  Instance expansion median:", round(expansion_median, 2), "ms\n"))
  cat(paste("  Speedup factor:", round(speedup_factor, 1), "x\n"))
  cat(paste("  Results equivalent:", results_match, "\n\n"))

  return(list(
    scale = scale,
    timing_results = timing_results,
    timing_summary = timing_summary,
    memory_info = memory_info,
    speedup_factor = speedup_factor,
    results_match = results_match,
    pvalue_difference = pvalue_diff,
    dataset_info = test_data$scale_info
  ))
}

#' Comprehensive Benchmark Suite
#'
#' Runs benchmarks across all scales and summarizes results
run_comprehensive_benchmark <- function(n_iterations = 10) {

  cat("==========================================\n")
  cat("Copy-Weighted Hypergeometric Benchmarking\n")
  cat("==========================================\n\n")

  scales <- c("small", "medium", "large", "phr_like")
  benchmark_results <- list()

  # Run benchmarks for each scale
  for (scale in scales) {
    benchmark_results[[scale]] <- benchmark_single_scale(scale, n_iterations)
  }

  # Create summary table
  cat("SUMMARY TABLE\n")
  cat("=============\n")
  cat(sprintf("%-12s %-8s %-8s %-10s %-10s %-12s %-10s\n",
              "Scale", "QGenes", "BGenes", "QInst", "BGInst", "MemReduction", "Speedup"))
  cat(sprintf("%-12s %-8s %-8s %-10s %-10s %-12s %-10s\n",
              "--------", "------", "------", "-----", "-----", "-----------", "-------"))

  for (scale in scales) {
    result <- benchmark_results[[scale]]
    cat(sprintf("%-12s %-8d %-8d %-10d %-10d %-12.1fx %-10.1fx\n",
                scale,
                result$dataset_info$query_genes,
                result$dataset_info$background_genes,
                result$memory_info$query_instances,
                result$memory_info$background_instances,
                result$memory_info$memory_reduction_factor,
                result$speedup_factor))
  }

  cat("\n")

  # Performance analysis
  cat("PERFORMANCE ANALYSIS\n")
  cat("====================\n")

  speedups <- sapply(benchmark_results, function(x) x$speedup_factor)
  memory_reductions <- sapply(benchmark_results, function(x) x$memory_info$memory_reduction_factor)

  cat(paste("Average speedup:", round(mean(speedups), 1), "x\n"))
  cat(paste("Average memory reduction:", round(mean(memory_reductions), 1), "x\n"))

  # Check if speedup improves with scale
  if (speedups["large"] > speedups["small"]) {
    cat("✓ Speedup increases with dataset size (good scalability)\n")
  } else {
    cat("! Speedup does not improve with scale\n")
  }

  # Check result consistency
  all_match <- all(sapply(benchmark_results, function(x) x$results_match))
  if (all_match) {
    cat("✓ All methods produce equivalent results\n")
  } else {
    cat("✗ Some methods produce different results\n")
  }

  cat("\n")

  # Recommendations
  cat("RECOMMENDATIONS\n")
  cat("===============\n")

  min_speedup_scale <- names(speedups)[which.min(speedups)]
  max_speedup_scale <- names(speedups)[which.max(speedups)]

  if (min(speedups) >= 2) {
    cat("✓ Parameter weighting consistently faster across all scales\n")
  } else {
    cat(paste("! Parameter weighting may be slower for", min_speedup_scale, "scale\n"))
  }

  if (max(memory_reductions) > 10) {
    cat("✓ Significant memory savings, especially for high copy number data\n")
  }

  if (speedups["phr_like"] > speedups["medium"]) {
    cat("✓ Particularly beneficial for PHR-like high copy number datasets\n")
  }

  cat("\n==========================================\n")

  return(list(
    individual_results = benchmark_results,
    summary_stats = list(
      mean_speedup = mean(speedups),
      mean_memory_reduction = mean(memory_reductions),
      all_results_match = all_match
    )
  ))
}

#' Copy Number Distribution Impact Analysis
#'
#' Tests how different copy number distributions affect performance
analyze_copy_distribution_impact <- function() {

  cat("Analyzing impact of copy number distributions...\n\n")

  # Fixed dataset size, varying copy number distributions
  n_genes_query <- 50
  n_genes_background <- 1000

  distributions <- list(
    uniform_low = list(min = 1, max = 3, name = "Uniform 1-3"),
    uniform_medium = list(min = 1, max = 10, name = "Uniform 1-10"),
    poisson_low = list(lambda = 2, name = "Poisson(λ=2)"),
    poisson_high = list(lambda = 20, name = "Poisson(λ=20)"),
    extreme = list(extreme = TRUE, name = "Extreme (1-100)")
  )

  results <- list()

  for (dist_name in names(distributions)) {
    dist_params <- distributions[[dist_name]]

    set.seed(42)

    # Generate copy numbers based on distribution
    if ("extreme" %in% names(dist_params)) {
      # Mix of single copies and high copies
      query_copies <- c(rep(1, 45), rep(50, 5))
      background_copies <- c(rep(1, 900), rep(50, 100))
    } else if ("lambda" %in% names(dist_params)) {
      query_copies <- rpois(n_genes_query, dist_params$lambda) + 1
      background_copies <- rpois(n_genes_background, dist_params$lambda) + 1
    } else {
      query_copies <- sample(dist_params$min:dist_params$max, n_genes_query, replace = TRUE)
      background_copies <- sample(dist_params$min:dist_params$max, n_genes_background, replace = TRUE)
    }

    # Create datasets
    query_df <- data.frame(
      gene = paste0("Q", 1:n_genes_query),
      copy_number = query_copies
    )

    background_df <- data.frame(
      gene = paste0("G", 1:n_genes_background),
      copy_number = background_copies
    )

    pathway_genes <- paste0("Q", 1:10)  # First 10 query genes in pathway

    # Timing test
    timing <- microbenchmark(
      parameter_weighting = weighted_hypergeometric_test(
        query_df, pathway_genes, background_df, validate_inputs = FALSE
      ),
      instance_expansion = instance_expansion_test(
        query_df, pathway_genes, background_df
      ),
      times = 5,
      unit = "ms"
    )

    # Memory analysis
    memory_info <- estimate_memory_usage(query_df, background_df)

    # Summarize
    timing_summary <- summary(timing)
    weighting_time <- timing_summary$median[timing_summary$expr == "parameter_weighting"]
    expansion_time <- timing_summary$median[timing_summary$expr == "instance_expansion"]
    speedup <- expansion_time / weighting_time

    results[[dist_name]] <- list(
      distribution = dist_params$name,
      total_instances = sum(query_copies) + sum(background_copies),
      mean_copy_number = mean(c(query_copies, background_copies)),
      max_copy_number = max(c(query_copies, background_copies)),
      speedup = speedup,
      memory_reduction = memory_info$memory_reduction_factor,
      weighting_time_ms = weighting_time,
      expansion_time_ms = expansion_time
    )

    cat(sprintf("%-15s: %4.1fx speedup, %4.1fx memory reduction, %d total instances\n",
                dist_params$name, speedup, memory_info$memory_reduction_factor,
                sum(query_copies) + sum(background_copies)))
  }

  cat("\nConclusions:\n")
  speedups <- sapply(results, function(x) x$speedup)
  highest_speedup <- names(speedups)[which.max(speedups)]

  if (max(speedups) > 2 * min(speedups)) {
    cat(paste("- Copy number distribution significantly affects performance\n"))
    cat(paste("- Highest speedup with", results[[highest_speedup]]$distribution, "\n"))
  } else {
    cat("- Performance relatively stable across copy number distributions\n")
  }

  return(results)
}

#' Validation Overhead Analysis
#'
#' Measures the cost of input validation and error checking
analyze_validation_overhead <- function() {

  cat("Analyzing validation overhead...\n\n")

  # Use medium-scale test data
  test_data <- generate_test_data("medium")

  # Benchmark with and without validation
  timing_results <- microbenchmark(
    with_full_validation = weighted_hypergeometric_test(
      test_data$query_df, test_data$pathway_genes, test_data$background_df,
      validate_inputs = TRUE, handle_zeros = TRUE, max_copies = 500
    ),
    minimal_validation = weighted_hypergeometric_test(
      test_data$query_df, test_data$pathway_genes, test_data$background_df,
      validate_inputs = FALSE, handle_zeros = FALSE, max_copies = NULL
    ),
    times = 20,
    unit = "ms"
  )

  timing_summary <- summary(timing_results)
  full_validation_time <- timing_summary$median[timing_summary$expr == "with_full_validation"]
  minimal_validation_time <- timing_summary$median[timing_summary$expr == "minimal_validation"]
  validation_overhead <- full_validation_time - minimal_validation_time

  cat(sprintf("Full validation:    %6.2f ms\n", full_validation_time))
  cat(sprintf("Minimal validation: %6.2f ms\n", minimal_validation_time))
  cat(sprintf("Validation overhead:%6.2f ms (%.1f%%)\n",
              validation_overhead,
              100 * validation_overhead / full_validation_time))

  if (validation_overhead / full_validation_time < 0.2) {
    cat("✓ Validation overhead is reasonable (<20% of total time)\n")
  } else {
    cat("! High validation overhead - consider optimizing for production use\n")
  }

  return(list(
    full_validation_time_ms = full_validation_time,
    minimal_validation_time_ms = minimal_validation_time,
    validation_overhead_ms = validation_overhead,
    validation_overhead_percent = 100 * validation_overhead / full_validation_time
  ))
}

#' Main Benchmark Runner
#'
#' Runs all benchmarking analyses
main_benchmark <- function() {
  cat("Starting comprehensive benchmark suite...\n\n")

  # Main performance benchmark
  main_results <- run_comprehensive_benchmark(n_iterations = 10)

  # Copy number distribution analysis
  dist_results <- analyze_copy_distribution_impact()

  # Validation overhead analysis
  validation_results <- analyze_validation_overhead()

  # Final summary
  cat("\nFINAL RECOMMENDATIONS\n")
  cat("=====================\n")

  avg_speedup <- main_results$summary_stats$mean_speedup
  avg_memory_reduction <- main_results$summary_stats$mean_memory_reduction

  if (avg_speedup >= 3 && avg_memory_reduction >= 10) {
    cat("✓ EXCELLENT: Parameter weighting strongly recommended\n")
  } else if (avg_speedup >= 1.5 && avg_memory_reduction >= 5) {
    cat("✓ GOOD: Parameter weighting recommended\n")
  } else {
    cat("? MIXED: Consider dataset characteristics\n")
  }

  cat(paste("- Average", round(avg_speedup, 1), "x faster\n"))
  cat(paste("- Average", round(avg_memory_reduction, 1), "x less memory\n"))

  if (validation_results$validation_overhead_percent < 20) {
    cat("- Validation overhead acceptable for production\n")
  } else {
    cat("- Consider disabling validation for performance-critical applications\n")
  }

  return(list(
    main_benchmark = main_results,
    distribution_analysis = dist_results,
    validation_analysis = validation_results
  ))
}

# Run benchmark if this script is executed directly
if (sys.nframe() == 0) {
  benchmark_results <- main_benchmark()
}