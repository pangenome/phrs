# Comprehensive Statistical Validation Report: Copy-Number Weighted Hypergeometric Test

**Generated:** 2026-04-01  
**Task:** statistical-validation-framework  
**Status:** CRITICAL ISSUES IDENTIFIED

## Executive Summary

Comprehensive statistical validation of the copy-number weighted hypergeometric test has identified **fundamental statistical issues** that require immediate attention. While the mathematical implementation is correct, the approach violates independence assumptions when applied to gene-level sampling scenarios, leading to systematic Type I error inflation.

### Key Findings

- ✅ **Mathematical Equivalence**: Parameter weighting approach is mathematically equivalent to instance expansion
- ✅ **Parameter Constraints**: All constraint validation and edge case handling works correctly  
- ❌ **Null Distribution**: Non-uniform p-values under null hypothesis across all tested scenarios
- ❌ **Type I Error Control**: Systematic inflation of false positive rates (2-5x expected levels)
- ❌ **Statistical Validity**: Approach violates hypergeometric independence assumptions

### Critical Issue

**The weighted hypergeometric test produces anti-conservative p-values when genes are sampled at the gene level but tested at the instance level, violating the independence assumptions required for valid hypergeometric testing.**

## Detailed Validation Results

### 1. Null Distribution Validation

**Status: FAILED** - All 81 scenarios tested showed non-uniform p-values

**Test Summary:**
- Scenarios tested: 81 (3 copy distributions × 3 background sizes × 3 pathway sizes × 3 query sizes)
- Scenarios passing uniformity: **0** (0%)
- Consistent pattern: KS test p-value ≈ 0 across all scenarios

**Representative Results:**
```
Scenario: uniform_500_50_30
- Background: 500 genes, Pathway: 50 genes, Query: 30 genes
- Type I error rate: 0.222 (expected: 0.05)
- KS test p-value: 0 (uniformity rejected)

Scenario: realistic_1000_100_50  
- Background: 1000 genes, Pathway: 100 genes, Query: 50 genes
- Type I error rate: 0.196 (expected: 0.05)
- KS test p-value: 0 (uniformity rejected)
```

### 2. Type I Error Rate Analysis

**Pattern Identified:** Type I error inflation increases with:
1. Higher average copy numbers
2. Greater copy number variance
3. Larger effective sample sizes

**Copy Number Effects:**
```
Copy Distribution    Type I Rate    Expected
All CN=1            0.051          0.05    ✓ (standard case)
All CN=5            0.089          0.05    ✗ (1.8x inflation)
CN 1-8 uniform      0.222          0.05    ✗ (4.4x inflation)
CN 1-20 variable    0.312          0.05    ✗ (6.2x inflation)
```

### 3. Mathematical Equivalence Validation

**Status: PASSED** - Parameter weighting method is mathematically equivalent to instance expansion

**Test Results:**
- Equivalence tests: 100/100 passed
- Parameter matching rate: 100%
- Maximum p-value difference: 0 (numerical precision)

### 4. Parameter Constraint Validation

**Status: PASSED** - All constraint validation mechanisms work correctly

**Tests Passed:**
- ✅ Basic hypergeometric constraints
- ✅ Edge case detection and handling
- ✅ Boundary condition validation  
- ✅ Parameter transformation consistency
- ✅ Constraint violation recovery

## Root Cause Analysis

### The Independence Violation Problem

The hypergeometric distribution assumes that each "draw" from the population is independent. However, in the copy-number weighted scenario:

1. **Gene-Level Sampling**: Researchers typically select genes (not instances) for analysis
2. **Instance-Level Testing**: The weighted test counts all copies of selected genes
3. **Clustered Selection**: When a gene is selected, ALL its copies are selected together

This creates a **cluster sampling** effect that violates independence assumptions:

```
Standard Model:        Gene1, Gene2, Gene3, ...    (independent draws)
Weighted Reality:      Gene1×14, Gene2×3, Gene3×8, ... (clustered draws)
```

### Why This Causes Anti-Conservative P-Values

1. **Effective Sample Size Inflation**: Copy number weighting increases k without proportional increase in variance
2. **Reduced Sampling Variance**: Clustered selection reduces the random variation expected under the null
3. **Hypergeometric Assumption Violation**: The test assumes random sampling of instances, but receives clustered instances

## Recommendations

### 1. When to Use Copy-Number Weighted Testing

**Appropriate Scenarios:**
- Instance-level sampling (e.g., RNA-seq reads, protein molecules)
- Known enrichment testing (positive controls)
- Comparative analysis (relative enrichment between conditions)
- Dosage effect hypothesis testing

**Inappropriate Scenarios:**
- Gene-level discovery analysis
- Null hypothesis testing with gene-level sampling
- Significance testing for pathway enrichment from gene lists

### 2. Statistical Corrections

#### Option A: Effective Sample Size Correction
```r
# Adjust k to account for clustering
k_effective <- k_weighted / mean(copy_numbers)
```

#### Option B: Permutation-Based Testing
```r
# Generate null distribution via gene-level permutation
null_pvals <- replicate(n_perms, {
  permuted_genes <- sample(background_genes, query_size)
  # Calculate weighted test statistic
})
```

#### Option C: Hybrid Approach
```r
# Use copy weighting for effect size, standard test for significance
effect_size <- weighted_fold_enrichment
significance <- standard_hypergeometric_pvalue
```

### 3. Validation Framework Enhancements

#### Required Validation Components

1. **Gene-Level Null Validation**
   - Test uniformity under gene-level sampling
   - Verify Type I error control
   - Cross-validate with permutation tests

2. **Instance-Level Null Validation**  
   - Test uniformity under instance-level sampling
   - Validate theoretical predictions
   - Benchmark against known enriched sets

3. **Power Analysis Framework**
   - Compare detection power vs. standard methods
   - Analyze power vs. copy number variance
   - Evaluate false discovery rate control

4. **Robustness Testing**
   - Test across different copy number distributions
   - Validate parameter constraint handling
   - Stress test with extreme scenarios

## Proposed Solutions

### 1. Enhanced Validation Suite

```r
run_enhanced_validation <- function() {
  # Gene-level null validation with corrections
  gene_level_results <- validate_gene_level_null()
  
  # Instance-level validation 
  instance_level_results <- validate_instance_level_null()
  
  # Permutation-based validation
  permutation_results <- validate_permutation_approach()
  
  # Comparative power analysis
  power_analysis <- compare_detection_power()
  
  return(comprehensive_results)
}
```

### 2. Corrected Implementation

```r
corrected_weighted_test <- function(query_df, pathway_genes, background_df, 
                                   method = "permutation") {
  if (method == "permutation") {
    return(permutation_based_test(...))
  } else if (method == "effective_sample") {
    return(effective_sample_correction(...))
  } else {
    warning("Using uncorrected method - results may be anti-conservative")
    return(standard_weighted_test(...))
  }
}
```

## Implementation Timeline

### Phase 1: Critical Fixes (Immediate)
1. Document statistical limitations in existing code
2. Implement permutation-based alternative
3. Add warnings for inappropriate use cases
4. Create corrected validation framework

### Phase 2: Enhanced Methods (1-2 weeks)
1. Implement effective sample size corrections
2. Develop hybrid significance/effect size approaches
3. Create comprehensive power analysis tools
4. Validate corrections with simulation studies

### Phase 3: Production Ready (2-4 weeks)
1. Integrate all validation components
2. Create user guidance and decision trees
3. Develop method selection algorithms
4. Comprehensive documentation and examples

## Conclusion

The statistical validation has successfully identified critical issues with the copy-number weighted hypergeometric approach that must be addressed before production use. While the mathematical implementation is correct, the violation of independence assumptions makes the approach inappropriate for standard gene-level discovery analyses.

The validation framework itself is working correctly - it has detected real statistical problems that require methodological solutions, not just implementation fixes. This demonstrates the importance of comprehensive statistical validation for novel analytical approaches.

**Next Steps:**
1. Implement corrected methods (permutation-based testing)
2. Enhance validation framework with appropriate null models  
3. Provide clear guidance on when weighted approaches are statistically valid
4. Validate corrected approaches against established benchmarks

This validation exercise has prevented deployment of a statistically flawed method and provides a roadmap for developing robust alternatives that maintain the biological insights of copy-number weighting while ensuring statistical validity.