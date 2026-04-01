# Mathematical Verification Report: Copy-Number Parameter Mapping

**Task:** mathematical-verification-of  
**Date:** 2026-04-01  
**Status:** In Progress  

## Executive Summary

This report provides comprehensive mathematical verification of the copy-number-weighted hypergeometric parameter mapping implementation in `copy_number_phyper_mapping.R`. The analysis reveals both strengths and critical issues that require attention.

## 1. Mathematical Foundation Verification

### 1.1 Parameter Transformation Logic

The implementation correctly implements the theoretical framework:

**Standard → Weighted Parameter Mapping:**
- `k_standard → k_weighted`: ✓ Correctly sums copy numbers for query instances
- `q_standard → q_weighted`: ✓ Correctly sums copy numbers for overlap instances  
- `m_standard → m_weighted`: ✓ Correctly sums copy numbers for pathway instances in background
- `n_standard → n_weighted`: ✓ Correctly calculates non-pathway instances as total - pathway

**Verification:** Manual calculations confirm the transformation logic is mathematically sound.

### 1.2 Mathematical Equivalence with Instance Expansion

**Test Results:**
- ✓ Parameter equivalence: TRUE
- ✓ P-value equivalence: TRUE (difference = 0.00e+00)
- ✓ Memory efficiency: 1.8x reduction factor

**Conclusion:** The parameter weighting approach is mathematically equivalent to instance expansion while being computationally superior.

## 2. Hypergeometric Constraint Validation

### 2.1 Critical Issue Identified

**CONSTRAINT VIOLATION:** The implementation can produce parameters that violate hypergeometric constraints.

**Specific Case Found:**
```
Query instances (k_weighted):    284
Overlap instances (q_weighted):  72  
Pathway instances (m_weighted):  9
Background instances (n_weighted): 31107
```

**Violation:** q_weighted (72) > m_weighted (9)

This violates the fundamental hypergeometric constraint: overlap ≤ pathway size.

### 2.2 Root Cause Analysis

The constraint violation occurs when:
1. Query genes have high copy numbers
2. Background pathway genes have low copy numbers  
3. Overlap genes exist with copy numbers exceeding background pathway totals

**Example Scenario:**
- Query gene "OR4F17" has 14 copies
- Background "OR4F17" has only 2 copies  
- Overlap instances = 14 (from query)
- Pathway instances = 2 (from background)
- Result: 14 > 2 → constraint violation

### 2.3 Theoretical Implications

This reveals a fundamental issue with the parameter mapping approach:
- Query copy numbers can exceed background copy numbers for the same gene
- The hypergeometric model assumes sampling from a finite population
- When query "sample" exceeds population size, the model breaks down

## 3. Edge Case Testing Results

### 3.1 Basic Parameter Mapping
✓ **PASS**: Simple cases work correctly
- k_weighted: 24 instances (expected 24) ✓
- q_weighted: 6 instances (expected 6) ✓
- All constraints satisfied ✓

### 3.2 Zero Overlap Cases  
✓ **PASS**: No overlap handled correctly
- q_weighted = 0 ✓
- p-value = 1.0 ✓

### 3.3 Complete Overlap Cases
✓ **PASS**: Full overlap handled correctly
- q_weighted = k_weighted ✓

### 3.4 High Copy Number Cases
⚠️ **ISSUE**: Can trigger constraint violations
- Single gene with 50 copies works in isolation
- Fails when background has fewer copies

### 3.5 Data Cleaning
✓ **PASS**: Zero copy genes properly filtered

## 4. Statistical Properties Assessment

### 4.1 P-value Calculation
- ✓ Uses correct phyper() parameterization
- ✓ Lower.tail=FALSE for P(X ≥ q) correctly implemented
- ✓ Handles q-1 adjustment properly

### 4.2 Fold Enrichment Calculation
- ✓ Copy-weighted fractions calculated correctly
- ✓ Handles zero denominators appropriately  

### 4.3 Parameter Validation System
- ✓ Checks non-negative parameters
- ✓ Checks integer constraints
- ✓ Checks overlap ≤ query constraint
- ⚠️ **FAILS** on overlap ≤ pathway constraint in edge cases
- ✓ Checks query ≤ total population constraint

## 5. Critical Issues Summary

### 5.1 Constraint Violation Problem
**Severity:** HIGH  
**Impact:** Makes hypergeometric test invalid  
**Frequency:** Occurs with PHR data and high copy number scenarios

### 5.2 Recommended Solutions

**Option 1: Background Copy Number Adjustment**
```r
# Ensure background copy numbers ≥ query copy numbers for overlapping genes
for (gene in overlap_genes) {
  bg_copies <- background_df$copy_number[background_df$gene_name == gene]
  query_copies <- query_df$copy_number[query_df$gene_name == gene]
  if (bg_copies < query_copies) {
    background_df$copy_number[background_df$gene_name == gene] <- query_copies
    warning(paste("Adjusted background copy number for", gene))
  }
}
```

**Option 2: Query Copy Number Capping**
```r
# Cap query copy numbers to background maximums
for (gene in overlap_genes) {
  bg_copies <- background_df$copy_number[background_df$gene_name == gene]
  if (query_df$copy_number[query_df$gene_name == gene] > bg_copies) {
    query_df$copy_number[query_df$gene_name == gene] <- bg_copies
    warning(paste("Capped query copy number for", gene))
  }
}
```

**Option 3: Error Detection and Handling**
```r
# Detect and report constraint violations explicitly
if (q_weighted > m_weighted) {
  stop(paste("Constraint violation: overlap instances (", q_weighted, 
             ") > pathway instances (", m_weighted, 
             "). Check copy number consistency between query and background."))
}
```

## 6. Mathematical Correctness Assessment

### 6.1 Core Algorithm: ✅ CORRECT
- Parameter transformation logic is mathematically sound
- Equivalence with instance expansion verified
- Statistical test implementation correct

### 6.2 Constraint Handling: ❌ INCOMPLETE
- Missing validation for q ≤ m constraint in realistic scenarios
- No handling of copy number inconsistencies between query and background
- Can produce invalid hypergeometric parameters

### 6.3 Edge Case Robustness: ⚠️ PARTIAL
- Handles most edge cases correctly
- Fails on copy number inconsistency scenarios
- Missing comprehensive constraint enforcement

## 7. Verification Required Items

Per task specification, verifying:

### ✅ Mathematical constraints verified
- ✅ q≤k constraint: Properly validated
- ❌ q≤m constraint: **VIOLATION FOUND** 
- ✅ k≤m+n constraint: Properly validated

### ✅ Mathematical equivalence with instance expansion
- ✅ Parameters match exactly
- ✅ P-values identical (difference = 0.00e+00)
- ✅ Memory efficiency demonstrated

### ⚠️ Edge cases tested
- ✅ Zero overlap: Handled correctly
- ✅ Single gene pathways: Basic cases work  
- ❌ Extreme copy numbers: **CONSTRAINT VIOLATIONS**

### ✅ Statistical properties preservation
- ✅ P-value calculation correct
- ✅ Fold enrichment calculation correct
- ✅ Parameter validation framework present

## 8. Recommendations

### 8.1 Immediate Actions Required
1. **Fix constraint violation**: Implement copy number consistency checks
2. **Add pre-flight validation**: Check query vs background copy number compatibility  
3. **Improve error messages**: Provide clear guidance when constraints fail

### 8.2 Implementation Improvements
1. Add `validate_copy_consistency()` function
2. Add `adjust_background_copies()` option
3. Enhance parameter validation with detailed constraint checking

### 8.3 Documentation Updates
1. Document copy number consistency requirements
2. Add troubleshooting section for constraint violations
3. Provide guidance on background dataset preparation

## 9. Conclusion

**Mathematical Verification Status: PARTIAL PASS**

The core mathematical framework is **correct and validated**, but the implementation has a **critical constraint handling issue** that prevents use with real datasets like PHR data. The parameter transformation logic is sound and mathematically equivalent to instance expansion, but the system fails when query copy numbers exceed background copy numbers for the same genes.

**Recommendation:** Address constraint violation issue before production use. The mathematical foundation is solid, but operational robustness needs improvement.