#!/usr/bin/env Rscript

# Verification script for null-distribution-validation task
# Tests that:
# 1. R functions test p-value uniformity under null hypothesis
# 2. Kolmogorov-Smirnov tests pass (for actual uniform p-values)
# 3. Null simulation functions work correctly

source("null_distribution_test.R", local = TRUE)

# Test 1: Verify p-value uniformity testing function exists and works
cat("Testing p-value uniformity functions...\n")
set.seed(42)
uniform_pvals <- runif(100)
test_result <- test_pvalue_uniformity(uniform_pvals)

if (!test_result$overall_uniform) {
  cat("ERROR: Uniformity test failed for truly uniform p-values\n")
  quit(status = 1)
}

# Test 2: Verify Kolmogorov-Smirnov test passes for uniform data
cat("Testing Kolmogorov-Smirnov functionality...\n")
if (!test_result$ks_test$uniform) {
  cat("ERROR: KS test failed for uniform p-values\n")
  quit(status = 1)
}

# Test 3: Verify null simulation functions work
cat("Testing null simulation functions...\n")
sim_pvals <- simulate_null_pvalues(n_simulations = 10, background_size = 50,
                                   pathway_size = 10, query_size = 5)

if (length(sim_pvals) != 10 || any(is.na(sim_pvals)) ||
    any(sim_pvals < 0) || any(sim_pvals > 1)) {
  cat("ERROR: Null simulation functions not working correctly\n")
  quit(status = 1)
}

# Test 4: Verify comprehensive validation function exists
cat("Testing comprehensive validation framework...\n")
tryCatch({
  # Test with minimal scenario to verify it works
  minimal_scenarios <- list(
    list(name = "test", background_size = 50, pathway_size = 5,
         query_size = 10, copy_distribution = "uniform")
  )
  result <- validate_null_distribution_comprehensive(minimal_scenarios, n_simulations = 10)
  if (length(result) != 1) {
    stop("Comprehensive validation returned unexpected results")
  }
}, error = function(e) {
  cat("ERROR: Comprehensive validation framework failed:", e$message, "\n")
  quit(status = 1)
})

cat("\n✓ All verification checks passed!\n")
cat("✓ R functions test p-value uniformity under null hypothesis\n")
cat("✓ Kolmogorov-Smirnov tests pass for uniform data\n")
cat("✓ Null simulation functions work correctly\n")
cat("✓ Comprehensive validation framework functional\n")
quit(status = 0)