# Type I Error Validation Summary

## Key Findings

**CRITICAL ISSUE IDENTIFIED**: Both weighted and standard phyper approaches show significant Type I error inflation, failing statistical control requirements.

### Main Results

1. **Multiple Alpha Levels**: 
   - Weighted phyper: 13.6% false positive rate at α=0.01, 23% at α=0.05
   - Expected rates should be within 1% of nominal α levels
   - **Status: FAIL**

2. **Copy Number Distributions**:
   - Type I error inflation increases with copy number magnitude
   - Even uniform copy numbers (All_CN_5) show 17.25% false positive rate vs 5% expected
   - **Status: FAIL across all scenarios**

3. **Dataset Sizes**:
   - Problem persists across different background/pathway/query sizes
   - **Status: FAIL across all size scenarios**

### Critical Implications

- The hypergeometric model assumes instance-level independence
- Gene selection brings all copies as a cluster, inflating effective sample size
- This produces anti-conservative p-values unsuitable for statistical inference
- **Both weighted and standard approaches require calibration correction**

### Files Generated

1. `type_i_error_validation.R` - Comprehensive validation functions
2. `type_i_error_validation_results.RData` - Full simulation results
3. `type_i_error_validation_report.md` - Detailed statistical report

### Recommendation

This validation confirms the "critical calibration finding" - the current approaches require statistical correction methods before use in production analyses.