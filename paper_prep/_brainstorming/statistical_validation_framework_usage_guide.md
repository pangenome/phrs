# Statistical Validation Framework Usage Guide

## Copy-Number Weighted Hypergeometric Test Validation

**Framework Version:** statistical-validation-framework-3  
**Created:** 2026-04-01  
**Task:** statistical-validation-framework-3  

---

## Overview

The Statistical Validation Framework provides comprehensive validation for the copy-number weighted hypergeometric test implementation. It integrates three critical validation domains:

1. **Null Distribution Validation** - Verifies p-values follow Uniform(0,1) under null hypothesis
2. **Type I Error Rate Control** - Validates false positive rates at multiple significance levels  
3. **Parameter Constraint Validation** - Tests parameter bounds and edge case handling

## Quick Start

```r
# Load the validation framework
source("statistical_validation_suite.R")

# Run comprehensive validation with default parameters
results <- run_comprehensive_validation()

# Generate detailed report
generate_validation_report(results)

# View summary
print(results$summary)
```

## Detailed Usage

### 1. Basic Validation

```r
# Run with default parameters (recommended for most users)
results <- run_comprehensive_validation()
```

**Default Configuration:**
- 2000 simulations for null distribution testing
- Tests uniform, skewed, and realistic copy number scenarios
- Multiple background/pathway/query size combinations
- Type I error testing at α = 0.001, 0.01, 0.05, 0.1
- Comprehensive parameter constraint validation

### 2. Custom Validation Configuration

```r
# Create custom configuration
custom_config <- list(
  null_distribution = list(
    n_simulations = 5000,
    scenarios = c("uniform", "realistic"),
    background_sizes = c(1000, 2000),
    pathway_sizes = c(100, 200),
    query_sizes = c(50, 100)
  ),
  type_i_error = list(
    n_simulations = 3000,
    alpha_levels = c(0.01, 0.05),
    tolerance = 0.01,
    scenarios = list(
      basic = list(n_genes = 1000, n_pathway = 100, n_query = 50),
      large = list(n_genes = 5000, n_pathway = 500, n_query = 250)
    )
  ),
  parameter_constraints = list(
    edge_case_scenarios = c("zero_copies", "extreme_ratios"),
    constraint_tests = c("parameter_bounds", "feasibility"),
    error_recovery_tests = TRUE
  ),
  output = list(
    save_detailed_results = TRUE,
    generate_plots = TRUE,
    export_summary = TRUE
  )
)

# Run with custom configuration
results <- run_comprehensive_validation(custom_config)
```

### 3. Individual Component Testing

```r
# Test only null distribution validation
source("null_distribution_test.R")
pvalues <- simulate_null_pvalues(n_simulations = 1000)
null_validation <- validate_pvalue_distribution(pvalues)

# Test only Type I error rates  
source("type_i_error_validation.R")
type_i_results <- validate_type_i_error_rates(n_sims = 1000)

# Test only parameter constraints
source("parameter_constraints_validation.R")
constraint_results <- validate_comprehensive_parameter_bounds()
```

### 4. Report Generation

```r
# Generate detailed markdown report
report_file <- generate_validation_report(results, "my_validation_report.md")

# Generate summary statistics
summary_stats <- generate_validation_summary(results)

# Print usage guidelines
print_validation_usage_guide()
```

## Interpretation Guidelines

### Overall Status

- **PASS** ✅ - All validation tests passed; implementation is statistically valid and ready for production use
- **FAIL** ❌ - One or more tests failed; issues must be addressed before production deployment

### Null Distribution Validation

**Purpose:** Validates that p-values from the weighted hypergeometric test follow the expected Uniform(0,1) distribution under the null hypothesis.

**Tests Performed:**
- Kolmogorov-Smirnov test for uniformity
- Anderson-Darling test for uniformity (if nortest package available)
- Visual Q-Q plot inspection
- Distribution moment checks

**Interpretation:**
- **PASS:** P-values properly follow uniform distribution - statistical foundation is sound
- **FAIL:** P-values deviate from uniform distribution - indicates fundamental statistical issues

**Common Failure Causes:**
- Incorrect weighted parameter calculations
- Improper background gene selection
- Bugs in hypergeometric parameter mapping
- Copy number aggregation errors

### Type I Error Rate Validation

**Purpose:** Validates that false positive rates are maintained at nominal levels across different significance thresholds and experimental scenarios.

**Tests Performed:**
- False positive rate testing at α = 0.001, 0.01, 0.05, 0.1
- Multiple experimental scenario validation
- Tolerance testing (error rates within 1% of nominal levels)
- Copy number distribution robustness testing

**Interpretation:**
- **PASS:** Error rates are controlled within acceptable bounds
- **FAIL:** Error rates are inflated (liberal) or deflated (conservative)

**Common Failure Causes:**
- Incorrect null hypothesis implementation  
- Biased query gene sampling
- Statistical test implementation errors
- Copy number bias in gene selection

### Parameter Constraint Validation

**Purpose:** Validates that parameter bounds are properly enforced and edge cases are handled gracefully.

**Tests Performed:**
- Parameter bound enforcement: q_weighted ≤ min(m_weighted, k_weighted)
- Edge case handling: zero copies, extreme ratios, boundary conditions
- Constraint violation detection and recovery
- Error handling robustness

**Interpretation:**
- **PASS:** Parameters are properly validated and constrained
- **FAIL:** Parameter validation is insufficient or error handling is inadequate

**Common Failure Causes:**
- Missing parameter validation logic
- Inadequate edge case handling
- Poor error recovery mechanisms
- Inconsistent constraint enforcement

## Configuration Reference

### Null Distribution Parameters

```r
null_distribution = list(
  n_simulations = 2000,                          # Number of simulation replicates
  scenarios = c("uniform", "skewed", "realistic"), # Copy number distributions
  background_sizes = c(500, 1000, 2000),         # Background gene counts
  pathway_sizes = c(50, 100, 200),               # Pathway gene counts  
  query_sizes = c(25, 50, 100)                   # Query gene counts
)
```

### Type I Error Parameters

```r
type_i_error = list(
  n_simulations = 2000,                          # Number of simulation replicates
  alpha_levels = c(0.001, 0.01, 0.05, 0.1),     # Significance levels to test
  tolerance = 0.01,                              # Acceptable deviation (1%)
  scenarios = list(                              # Test scenarios
    basic = list(n_genes = 1000, n_pathway = 100, n_query = 50),
    large = list(n_genes = 5000, n_pathway = 500, n_query = 250),
    small = list(n_genes = 200, n_pathway = 20, n_query = 10)
  )
)
```

### Parameter Constraint Parameters

```r
parameter_constraints = list(
  edge_case_scenarios = c(                       # Edge cases to test
    "zero_copies", 
    "extreme_ratios", 
    "boundary_conditions"
  ),
  constraint_tests = c(                          # Constraint types
    "parameter_bounds", 
    "feasibility", 
    "consistency"
  ),
  error_recovery_tests = TRUE                    # Test error handling
)
```

### Output Parameters

```r
output = list(
  save_detailed_results = TRUE,                  # Save .RData file
  generate_plots = TRUE,                         # Generate diagnostic plots
  export_summary = TRUE                          # Export summary statistics
)
```

## Troubleshooting

### Common Issues and Solutions

#### 1. Package Dependencies Missing

**Error:** Package 'nortest' not available
```r
# Install missing packages
install.packages("nortest")
```

#### 2. Source File Not Found

**Error:** cannot open file 'debug_weighted_phyper.R'
```r
# Ensure all required files are in the working directory:
# - debug_weighted_phyper.R
# - null_distribution_test.R  
# - type_i_error_validation.R
# - parameter_constraints_validation.R
# - edge_case_test_suite.R
# - constraint_violation_handler.R
```

#### 3. Memory Issues with Large Simulations

**Error:** Cannot allocate vector of size X
```r
# Reduce simulation counts
config$null_distribution$n_simulations <- 1000
config$type_i_error$n_simulations <- 1000
```

#### 4. Null Distribution Test Failures

**Investigation Steps:**
1. Check parameter calculation functions
2. Verify background gene selection logic  
3. Review hypergeometric test implementation
4. Examine copy number aggregation

#### 5. Type I Error Rate Failures

**Investigation Steps:**
1. Verify null hypothesis implementation
2. Check for bias in gene sampling
3. Review statistical test parameters
4. Examine error rate calculations

#### 6. Parameter Constraint Failures  

**Investigation Steps:**
1. Check parameter validation functions
2. Review edge case handling logic
3. Verify error recovery mechanisms
4. Examine constraint enforcement

### Validation Performance

**Expected Runtime:**
- Basic validation (default): 10-15 minutes
- Comprehensive validation (extended): 30-45 minutes  
- Quick validation (reduced simulations): 2-5 minutes

**Memory Requirements:**
- Basic validation: 1-2 GB RAM
- Extended validation: 2-4 GB RAM
- Large-scale validation: 4-8 GB RAM

## Output Files

The validation framework generates several output files:

| File | Description |
|------|-------------|
| `statistical_validation_suite.R` | Integrated validation framework |
| `comprehensive_validation_results.RData` | Detailed validation results data |
| `statistical_validation_report.md` | Comprehensive validation report |
| `statistical_validation_framework_usage_guide.md` | This usage guide |

## Best Practices

### 1. Regular Validation

- Run comprehensive validation after any changes to core implementation
- Use quick validation for development iterations
- Archive validation results for reproducibility

### 2. Configuration Management

- Save custom configurations for reproducible testing
- Document any deviations from default parameters
- Version control validation configurations

### 3. Result Interpretation

- Always review the detailed report, not just overall status
- Investigate any unexpected failures thoroughly
- Document and track validation issues over time

### 4. Production Deployment

- **CRITICAL:** Only deploy implementations that achieve PASS status
- Archive validation reports with deployment records
- Establish validation as a mandatory pre-deployment step

## Framework Architecture

### Component Integration

```
statistical_validation_suite.R (main framework)
├── null_distribution_test.R (null distribution validation)
├── type_i_error_validation.R (Type I error control)  
├── parameter_constraints_validation.R (constraint validation)
├── edge_case_test_suite.R (edge case testing)
├── constraint_violation_handler.R (error handling)
└── debug_weighted_phyper.R (core implementation)
```

### Validation Workflow

1. **Initialization** - Load components and configuration
2. **Null Distribution Testing** - Validate p-value distributions  
3. **Type I Error Testing** - Validate error rate control
4. **Parameter Constraint Testing** - Validate bounds and edge cases
5. **Summary Generation** - Aggregate results and generate reports
6. **Report Creation** - Generate detailed markdown documentation

## Support and Maintenance

### Version History

- **v3 (statistical-validation-framework-3)** - Integrated comprehensive framework
- **v2** - Individual component validation
- **v1** - Basic validation implementation  

### Contributing

When extending or modifying the validation framework:

1. Maintain backward compatibility with existing configurations
2. Add comprehensive documentation for new validation components
3. Include example usage in the usage guide
4. Validate new components with existing test cases
5. Update version history and changelog

### Contact

For issues with the statistical validation framework:

1. Check troubleshooting section above
2. Review detailed error messages in validation output
3. Consult framework architecture for component interactions
4. Document issues for future framework improvements

---

*Statistical Validation Framework v3 - Comprehensive validation for copy-number weighted hypergeometric testing*