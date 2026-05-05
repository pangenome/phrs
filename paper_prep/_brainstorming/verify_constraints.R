#!/usr/bin/env Rscript
# Quick verification script for parameter constraint validation
# Returns exit code 0 if all validations pass, 1 if any fail

library(tidyverse, quietly = TRUE)
source("parameter_constraints_validation.R")

# Run basic validation tests
cat("Verifying parameter constraint validation system...\n")

# Test 1: Basic constraint validation
result1 <- validate_hypergeometric_parameters(5, 20, 80, 30)
test1_pass <- result1$constraints_satisfied

# Test 2: Edge case detection
result2 <- validate_hypergeometric_parameters(0, 20, 80, 30)
test2_pass <- "zero_overlap" %in% result2$edge_cases

# Test 3: Constraint violation detection
background_df <- data.frame(gene = paste0("G", 1:10), copy_number = rep(2, 10))
query_df <- background_df[1:3, ]
result3 <- detect_constraint_violations(query_df, c("G1", "G2"), background_df)
test3_pass <- result3$valid

all_pass <- test1_pass && test2_pass && test3_pass

cat("Parameter bounds validation: ", ifelse(test1_pass, "PASS", "FAIL"), "\n")
cat("Edge case handling: ", ifelse(test2_pass, "PASS", "FAIL"), "\n")
cat("Constraint violation detection: ", ifelse(test3_pass, "PASS", "FAIL"), "\n")
cat("Overall verification: ", ifelse(all_pass, "PASS", "FAIL"), "\n")

if (all_pass) {
  quit(status = 0)
} else {
  quit(status = 1)
}