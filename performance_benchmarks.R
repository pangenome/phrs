# Comprehensive Performance Benchmarks for Copy-Number Weighted Enrichment
# Production-ready performance testing and optimization validation

source("robust_copy_weighted_enrichment.R")

#' Comprehensive benchmark suite
#'
#' @param scale_factors Vector of scaling factors for dataset size
#' @param output_file File path for benchmark results
#' @param detailed_output Include detailed timing breakdowns
#'
#' @return Benchmark results
run_comprehensive_benchmarks <- function(scale_factors = c(1, 2, 5, 10),
                                       output_file = "performance_benchmarks.csv",
                                       detailed_output = TRUE) {

  cat("Running comprehensive performance benchmarks...\n")

  benchmark_results <- list()

  for (scale in scale_factors) {
    cat("Testing scale factor", scale, "...\n")

    # Generate scaled test dataset
    base_genes <- 50
    n_genes <- base_genes * scale
    genes <- paste0("GENE", sprintf("%04d", 1:n_genes))

    # Generate realistic copy number distribution
    copy_numbers <- generate_realistic_copy_distribution(n_genes)
    background_df <- data.frame(gene = genes, copy_number = copy_numbers)

    # Query and pathway sizes scale with dataset
    query_size <- min(round(15 * sqrt(scale)), n_genes)
    pathway_size <- min(round(12 * sqrt(scale)), n_genes)

    query_genes <- sample(genes, query_size)
    query_df <- background_df[background_df$gene %in% query_genes, ]
    pathway_genes <- sample(genes, pathway_size)

    # Run detailed benchmarks
    scale_results <- run_single_scale_benchmark(query_df, pathway_genes, background_df, scale)
    benchmark_results[[paste0("scale_", scale)]] <- scale_results

    if (detailed_output) {
      cat("  Scale", scale, "results:\n")
      cat("    Parameter method:", round(scale_results$parameter_mean_time * 1000, 2), "ms\n")
      cat("    Expansion method:", round(scale_results$expansion_mean_time * 1000, 2), "ms\n")
      cat("    Speedup:", round(scale_results$speedup_factor, 2), "x\n")
      cat("    Memory reduction:", round(scale_results$memory_reduction * 100, 1), "%\n\n")
    }
  }

  # Create summary data frame
  summary_df <- create_benchmark_summary(benchmark_results, scale_factors)

  # Write results
  write.csv(summary_df, output_file, row.names = FALSE)
  cat("Benchmark results written to:", output_file, "\n")

  # Performance analysis
  analysis <- analyze_performance_trends(summary_df)

  return(list(
    detailed_results = benchmark_results,
    summary = summary_df,
    analysis = analysis
  ))
}

#' Run benchmark for a single scale
run_single_scale_benchmark <- function(query_df, pathway_genes, background_df, scale, n_reps = 30) {

  # Benchmark parameter method
  param_times <- replicate(n_reps, {
    start_time <- Sys.time()
    result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df, validate_inputs = FALSE)
    end_time <- Sys.time()
    as.numeric(difftime(end_time, start_time, units = "secs"))
  })

  # Benchmark expansion method
  expansion_times <- replicate(n_reps, {
    start_time <- Sys.time()

    # Instance expansion
    query_expanded <- rep(query_df$gene, query_df$copy_number)
    background_expanded <- rep(background_df$gene, background_df$copy_number)

    # Standard calculation
    q_exp <- sum(query_expanded %in% pathway_genes)
    m_exp <- sum(background_expanded %in% pathway_genes)
    n_exp <- length(background_expanded) - m_exp
    k_exp <- length(query_expanded)

    # phyper call
    if (q_exp > 0) {
      p_val <- phyper(q_exp - 1, m_exp, n_exp, k_exp, lower.tail = FALSE)
    } else {
      p_val <- 1.0
    }

    end_time <- Sys.time()
    as.numeric(difftime(end_time, start_time, units = "secs"))
  })

  # Memory estimates
  query_instances <- sum(query_df$copy_number)
  background_instances <- sum(background_df$copy_number)

  param_memory <- nrow(query_df) + nrow(background_df)
  expansion_memory <- query_instances + background_instances

  # Calculate statistics
  list(
    scale_factor = scale,
    n_query_genes = nrow(query_df),
    n_background_genes = nrow(background_df),
    query_instances = query_instances,
    background_instances = background_instances,

    parameter_mean_time = mean(param_times),
    parameter_median_time = median(param_times),
    parameter_sd_time = sd(param_times),

    expansion_mean_time = mean(expansion_times),
    expansion_median_time = median(expansion_times),
    expansion_sd_time = sd(expansion_times),

    speedup_factor = mean(expansion_times) / mean(param_times),
    memory_reduction = 1 - (param_memory / expansion_memory),

    param_memory_estimate = param_memory,
    expansion_memory_estimate = expansion_memory,

    n_repetitions = n_reps
  )
}

#' Generate realistic copy number distribution
generate_realistic_copy_distribution <- function(n_genes) {
  # Based on observed genomic copy number distributions
  # Most genes have 1-2 copies, some have higher
  copy_probs <- c(0.6, 0.25, 0.1, 0.04, 0.01)  # For copies 1-5
  sample(1:5, n_genes, replace = TRUE, prob = copy_probs)
}

#' Create benchmark summary data frame
create_benchmark_summary <- function(benchmark_results, scale_factors) {

  summary_rows <- lapply(names(benchmark_results), function(scale_name) {
    result <- benchmark_results[[scale_name]]

    data.frame(
      scale_factor = result$scale_factor,
      n_query_genes = result$n_query_genes,
      n_background_genes = result$n_background_genes,
      query_instances = result$query_instances,
      background_instances = result$background_instances,
      parameter_mean_time_ms = result$parameter_mean_time * 1000,
      expansion_mean_time_ms = result$expansion_mean_time * 1000,
      speedup_factor = result$speedup_factor,
      memory_reduction_pct = result$memory_reduction * 100,
      param_memory = result$param_memory_estimate,
      expansion_memory = result$expansion_memory_estimate,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summary_rows)
}

#' Analyze performance trends
analyze_performance_trends <- function(summary_df) {

  # Scaling analysis
  scale_correlation <- cor(summary_df$scale_factor, summary_df$speedup_factor, use = "complete.obs")

  # Memory efficiency
  avg_memory_reduction <- mean(summary_df$memory_reduction_pct, na.rm = TRUE)

  # Time complexity analysis
  param_time_growth <- analyze_time_growth(summary_df$scale_factor, summary_df$parameter_mean_time_ms)
  expansion_time_growth <- analyze_time_growth(summary_df$scale_factor, summary_df$expansion_mean_time_ms)

  # Overall assessment
  performance_grade <- assess_performance(summary_df)

  list(
    scale_correlation = scale_correlation,
    avg_memory_reduction = avg_memory_reduction,
    param_time_complexity = param_time_growth,
    expansion_time_complexity = expansion_time_growth,
    performance_grade = performance_grade,
    min_speedup = min(summary_df$speedup_factor, na.rm = TRUE),
    max_speedup = max(summary_df$speedup_factor, na.rm = TRUE),
    consistent_improvement = all(summary_df$speedup_factor > 1.0, na.rm = TRUE)
  )
}

#' Analyze time growth patterns
analyze_time_growth <- function(scale_factors, times) {

  if (length(unique(scale_factors)) < 3) {
    return("insufficient_data")
  }

  # Fit linear and quadratic models
  linear_fit <- lm(times ~ scale_factors)
  quad_fit <- lm(times ~ scale_factors + I(scale_factors^2))

  # Compare fits
  linear_r2 <- summary(linear_fit)$r.squared
  quad_r2 <- summary(quad_fit)$r.squared

  if (quad_r2 - linear_r2 > 0.1) {
    return("quadratic")
  } else if (linear_r2 > 0.8) {
    return("linear")
  } else {
    return("sublinear")
  }
}

#' Assess overall performance
assess_performance <- function(summary_df) {

  # Criteria for performance grades
  min_speedup <- min(summary_df$speedup_factor, na.rm = TRUE)
  avg_speedup <- mean(summary_df$speedup_factor, na.rm = TRUE)
  min_memory_reduction <- min(summary_df$memory_reduction_pct, na.rm = TRUE)

  if (min_speedup > 2.0 && min_memory_reduction > 50) {
    return("excellent")
  } else if (min_speedup > 1.5 && min_memory_reduction > 30) {
    return("good")
  } else if (min_speedup > 1.1 && min_memory_reduction > 10) {
    return("adequate")
  } else {
    return("poor")
  }
}

#' Memory profiling benchmark
benchmark_memory_usage <- function(max_scale = 5) {

  cat("Running memory usage benchmark...\n")

  memory_results <- list()

  for (scale in 1:max_scale) {
    # Create test data
    n_genes <- 50 * scale
    genes <- paste0("GENE", sprintf("%04d", 1:n_genes))
    copy_numbers <- generate_realistic_copy_distribution(n_genes)
    background_df <- data.frame(gene = genes, copy_number = copy_numbers)

    query_size <- min(15 * scale, n_genes)
    query_genes <- sample(genes, query_size)
    query_df <- background_df[background_df$gene %in% query_genes, ]
    pathway_genes <- sample(genes, min(12 * scale, n_genes))

    # Memory measurement (approximate)
    param_objects <- c(
      object.size(query_df),
      object.size(background_df),
      object.size(pathway_genes)
    )

    # Expansion method memory
    query_expanded <- rep(query_df$gene, query_df$copy_number)
    background_expanded <- rep(background_df$gene, background_df$copy_number)

    expansion_objects <- c(
      object.size(query_expanded),
      object.size(background_expanded),
      object.size(pathway_genes)
    )

    memory_results[[scale]] <- list(
      scale = scale,
      param_memory_bytes = sum(param_objects),
      expansion_memory_bytes = sum(expansion_objects),
      memory_ratio = sum(expansion_objects) / sum(param_objects),
      query_instances = length(query_expanded),
      background_instances = length(background_expanded)
    )

    # Clean up large objects
    rm(query_expanded, background_expanded)
  }

  # Convert to data frame
  memory_df <- do.call(rbind, lapply(memory_results, function(x) {
    data.frame(
      scale = x$scale,
      param_memory_mb = x$param_memory_bytes / (1024^2),
      expansion_memory_mb = x$expansion_memory_bytes / (1024^2),
      memory_ratio = x$memory_ratio,
      query_instances = x$query_instances,
      background_instances = x$background_instances
    )
  }))

  cat("Memory usage results:\n")
  print(memory_df)

  return(memory_df)
}

#' Stress test with extreme parameters
stress_test_extreme_parameters <- function() {

  cat("Running stress tests with extreme parameters...\n")

  stress_tests <- list(
    "large_copy_numbers" = list(
      query_df = data.frame(gene = c("A", "B"), copy_number = c(1000, 2000)),
      background_df = data.frame(gene = c("A", "B", "C"), copy_number = c(1000, 2000, 500)),
      pathway_genes = c("A", "C")
    ),

    "many_genes_low_copy" = list(
      query_df = data.frame(gene = paste0("G", 1:1000), copy_number = rep(1, 1000)),
      background_df = data.frame(gene = paste0("G", 1:5000), copy_number = rep(1, 5000)),
      pathway_genes = paste0("G", 1:500)
    ),

    "unbalanced_copy_distribution" = list(
      query_df = data.frame(gene = c("A", "B", "C"), copy_number = c(1, 100, 1)),
      background_df = data.frame(gene = c("A", "B", "C", "D"), copy_number = c(1, 100, 1, 1)),
      pathway_genes = c("B", "D")
    )
  )

  stress_results <- list()

  for (test_name in names(stress_tests)) {
    cat("Running stress test:", test_name, "... ")

    test_data <- stress_tests[[test_name]]

    tryCatch({
      start_time <- Sys.time()
      result <- weighted_hypergeometric_test(
        test_data$query_df,
        test_data$pathway_genes,
        test_data$background_df
      )
      end_time <- Sys.time()

      computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

      stress_results[[test_name]] <- list(
        status = "success",
        computation_time = computation_time,
        p_value = result$p_value,
        parameters_valid = !is.na(result$p_value)
      )

      cat("PASS (", round(computation_time * 1000, 2), "ms)\n")

    }, error = function(e) {
      stress_results[[test_name]] <- list(
        status = "failed",
        error = as.character(e)
      )

      cat("FAIL:", e$message, "\n")
    })
  }

  return(stress_results)
}

#' Generate comprehensive performance report
generate_performance_report <- function(benchmark_results, output_file = "performance_report.txt") {

  report <- c(
    "=== Copy-Number Weighted Enrichment Performance Report ===",
    paste("Generated:", Sys.time()),
    "",
    "BENCHMARK SUMMARY:",
    ""
  )

  # Add summary statistics
  summary_df <- benchmark_results$summary
  analysis <- benchmark_results$analysis

  report <- c(report,
    paste("Performance Grade:", toupper(analysis$performance_grade)),
    paste("Consistent Improvement:", analysis$consistent_improvement),
    paste("Average Memory Reduction:", round(analysis$avg_memory_reduction, 1), "%"),
    paste("Speedup Range:", round(analysis$min_speedup, 2), "x to", round(analysis$max_speedup, 2), "x"),
    ""
  )

  # Add detailed results table
  report <- c(report, "DETAILED RESULTS:")
  report <- c(report, "Scale\tParam_Time(ms)\tExpansion_Time(ms)\tSpeedup\tMemory_Reduction(%)")

  for (i in 1:nrow(summary_df)) {
    row <- summary_df[i, ]
    report <- c(report, paste(
      row$scale_factor,
      round(row$parameter_mean_time_ms, 2),
      round(row$expansion_mean_time_ms, 2),
      round(row$speedup_factor, 2),
      round(row$memory_reduction_pct, 1),
      sep = "\t"
    ))
  }

  # Write report
  cat(paste(report, collapse = "\n"), file = output_file)
  cat("Performance report written to:", output_file, "\n")

  return(report)
}

# Quick performance validation
quick_performance_check <- function() {
  cat("Running quick performance check...\n")

  # Simple functional test - just verify the parameter method works
  query_df <- data.frame(gene = c("A", "B"), copy_number = c(2, 3))
  background_df <- data.frame(gene = c("A", "B", "C"), copy_number = c(2, 3, 1))
  pathway_genes <- c("A", "C")

  # Test parameter method
  tryCatch({
    result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)

    if (result$status == "success" && !is.na(result$p_value)) {
      cat("Quick performance check: PASS\n")
      cat("  Parameter method functional: YES\n")
      cat("  P-value computed:", round(result$p_value, 4), "\n")
      return(TRUE)
    } else {
      cat("Quick performance check: FAIL - method not functional\n")
      return(FALSE)
    }

  }, error = function(e) {
    cat("Quick performance check: FAIL - error:", e$message, "\n")
    return(FALSE)
  })
}