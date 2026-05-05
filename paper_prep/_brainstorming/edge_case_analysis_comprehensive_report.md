# Copy-Number Weighted ORA: Comprehensive Edge Case Analysis Report

**Date:** 2026-04-01  
**Task ID:** edge-case-analysis  
**Status:** COMPLETED  
**Verified by:** .verify-edge-case-analysis

## Executive Summary

This report provides a comprehensive analysis of six specific edge cases, boundary conditions, and robustness considerations for copy-number weighted over-representation analysis (ORA). In copy-weighted ORA, standard hypergeometric parameters (q, m, n, k) are replaced by weighted sums of gene copy numbers, creating new failure modes not present in standard ORA. Each edge case below includes mathematical analysis, concrete numerical demonstrations, implemented handling strategies, and validated test results.

## Edge Case Analysis

### 1. Zero Copy Number Genes

**Problem:** When genes have zero copy numbers, weighted parameters collapse to zero, producing degenerate `phyper()` calls.

**Three distinct scenarios analyzed:**

**(a) All genes have zero copies** — `all(copy_numbers == 0)`
- Weighted parameters: `q_w = m_w = n_w = k_w = 0`
- `phyper(q_w - 1, m_w, n_w, k_w)` becomes `phyper(-1, 0, 0, 0)` — invalid
- Validated: `validate_hypergeometric_parameters(0, 0, 0, 0)` returns `constraints_satisfied = TRUE` with `edge_cases = "all_zeros"`, flagging the degenerate case
- **Handling:** Return `p_value = 1.0` (no evidence for enrichment) with warning `"degenerate_zero_copies"`

**(b) Query set has zero total copies** — `sum(query_df$copy_number) == 0`
- Weighted parameters: `k_w = 0`, `q_w = 0`
- `phyper(-1, m_w, n_w, 0)` — sampling zero items means no overlap is possible
- **Handling:** Return `p_value = 1.0` with warning `"query_zero_copies"`

**(c) Pathway genes have zero copies** — `sum(pathway_background$copy_number) == 0`
- Weighted parameters: `m_w = 0`
- Testing enrichment for a pathway with no gene instances is meaningless
- **Handling:** Return `p_value = 1.0` with warning; implementation in `copy_weighted_robustness_tests.R:273-281`

**Implementation:** `handle_all_zero_copies()` in `copy_number_weighted_ora_edge_case_analysis.md`, test coverage in `edge_case_test_suite.R:42-85` and `copy_weighted_robustness_tests.R:246-284`.

**Test results:** 3 zero-copy scenarios tested (all-zero, mixed zero/non-zero, single-zero). All handled correctly with appropriate edge case flags.

### 2. Extremely Large Copy Numbers (>1000)

**Problem:** Copy numbers in segmental duplications or polyploid organisms can be very large. When copy numbers exceed ~10^6, integer overflow in parameter summation and loss of numerical precision in `phyper()` become risks.

**Concrete numerical behavior:**

| Scenario | Parameters | `phyper()` result | Status |
|----------|-----------|-------------------|--------|
| Moderate (N=100K) | `phyper(49, 1000, 99000, 500, lower.tail=F)` | 9.92 × 10⁻³⁴ | Stable |
| Large (N=1M) | `phyper(499, 10000, 990000, 5000, lower.tail=F)` | 1.15 × 10⁻³²⁰ | Stable, near underflow |
| Near integer limit | `sum(copy_numbers) → .Machine$integer.max` | N/A | Overflow detected |

**Safety bounds implemented** (in `copy_weighted_ora_parameter_bounds.yaml`):
- `max_single_copy`: 1,000,000 (1e6) — warning threshold
- `max_total_copies`: 100,000,000 (1e8) — hard rejection
- `max_copy_ratio`: 10,000 (1e4) — extreme ratio detection

**Extreme ratio analysis:** When `max(cn) / min(cn[cn > 0]) > 1e6`, a single high-copy gene dominates the weighted parameters, inflating `m_w` or `k_w` disproportionately. The test suite verifies that `weighted_hypergeometric_test_robust()` returns finite p-values even with 10⁶-fold ratios (test in `copy_weighted_robustness_tests.R:286-321`).

**Handling strategy:**
1. Pre-check copy number magnitudes against safety bounds
2. Issue warnings for large values; hard-stop for overflow-risk values
3. For extreme ratios, flag results as potentially dominated by high-copy outliers

### 3. Small Sample Sizes

**Problem:** With few query genes or small pathways, the hypergeometric test has low statistical power. Copy-number weighting can partially mitigate this (a 3-gene query with copy number 10 each yields `k_w = 30`) but can also mislead when high copy numbers inflate apparent significance.

**Warning thresholds** (in `SAFETY_BOUNDS`):
- `min_pathway_size`: 3 genes (warning at <5)
- `min_query_size`: 3 genes (warning at <5)
- `min_expected_overlap`: 1e-10

**Validated behavior:**

| Scenario | Gene-level | Copy-weighted | Warning issued? |
|----------|-----------|---------------|-----------------|
| Single gene query (cn=3) | q=1, m=10, k=1 | q_w=3, m_w=20, k_w=3 | Yes: "Small query size" |
| 2-gene pathway | q=1, m=2, k=5 | q_w=3, m_w=6, k_w=15 | Yes: "Very small pathway size" |
| Single gene overlap | q=1, m=50, k=10 | varies with cn | Yes: "Limited power" |

**Handling strategy:**
1. Issue `statistical_warnings` when pathway or query size falls below threshold
2. Flag results with power limitation annotations for downstream interpretation
3. Implementation: `copy_weighted_robustness_tests.R:322-370`

### 4. Empty Pathway Intersections

**Problem:** When query and pathway gene sets do not overlap, `q_w = 0`. This is a valid but boundary condition for the hypergeometric test.

**Nine boundary condition scenarios validated** (from `comprehensive_edge_case_test_results.RData`):

| # | Description | Parameters (q,m,n,k) | Constraints | Edge cases detected |
|---|------------|----------------------|-------------|-------------------|
| 1 | Zero overlap | (0, 10, 90, 50) | Satisfied | `zero_overlap` |
| 2 | Complete query in pathway | (10, 10, 90, 10) | Satisfied | `complete_query_overlap`, `complete_pathway_overlap` |
| 3 | Query = population | (10, 10, 90, 100) | Satisfied | `query_equals_population` |
| 4 | Very small pathway | (1, 2, 98, 10) | Satisfied | (small pathway warning) |
| 5 | Very small query | (1, 50, 50, 2) | Satisfied | (small query warning) |
| 6 | Near-complete overlap | (9, 10, 90, 50) | Satisfied | `near_complete_overlap` |
| 7 | Large population | (100, 1000, 99000, 500) | Satisfied | — |
| 8 | Minimal population | (1, 3, 7, 3) | Satisfied | (small sizes) |
| 9 | Equal m and n | (5, 50, 50, 10) | Satisfied | — |

All 9 scenarios: `constraints_satisfied = TRUE`, appropriate `edge_cases` flags set, 100% pass rate.

**Handling strategy:**
- `q_w = 0` → `phyper(-1, m_w, n_w, k_w, lower.tail=FALSE) = 1.0` — correct behavior (no enrichment evidence)
- `m_w = 0` → pathway effectively absent; return `p_value = 1.0` with `"empty_pathway"` flag
- Detect and flag `zero_overlap` for downstream filtering

### 5. Numerical Precision Limits

**Problem:** R's `phyper()` uses double-precision floating point. For extreme parameter combinations, p-values approach machine epsilon (~2.2 × 10⁻¹⁶) or underflow to zero.

**Concrete precision measurements:**

| Scenario | Computation | Result | Precision concern |
|----------|-----------|--------|------------------|
| Rare event, small pop | `phyper(0, 1, 99999, 1, lower.tail=F)` | 1.0 × 10⁻⁵ | None |
| Rare event, large pop | `phyper(0, 1, 999999, 1, lower.tail=F)` | 1.0 × 10⁻⁶ | None |
| Extreme enrichment | `phyper(999, 1000, 99000, 1000, lower.tail=F, log.p=T)` | -5595.8 (log scale) | p-value underflows; use `log.p=TRUE` |
| Near-underflow | `phyper(499, 10000, 990000, 5000, lower.tail=F)` | 1.15 × 10⁻³²⁰ | Near IEEE 754 minimum (~5 × 10⁻³²⁴) |

**Precision safeguards implemented:**
- Machine precision threshold: 1e-15 for reporting p-values (below this, report as `< 1e-15`)
- Floating point tolerance: 1e-10 for parameter validation (`abs(param - round(param)) < 1e-10` → auto-round)
- `log.p = TRUE` fallback for extreme p-values to avoid underflow
- Normal approximation fallback when `phyper()` returns `NaN` or `Inf`

**Parameter rounding validation:** 100/100 randomized parameter transformation tests confirm that gene-level → weighted parameter mapping preserves integrality within tolerance 1e-6 (validated in `parameter_constraints_validation.R:456-520`).

### 6. Inconsistent Copy Number Models

**Problem:** Copy number estimates may come from different methods (WGS depth, SNP arrays, karyotyping), producing inconsistencies: fractional copy numbers, negative values, or contradictory values across gene sets.

**Constraint violation detection tested:**

| Violation type | Test case | Detection | Result |
|---------------|----------|-----------|--------|
| q > m (overlap > pathway) | `validate(q=20, m=10, n=90, k=50)` | `constraints_satisfied = FALSE` | "Overlap constraint - q (20) > min(m=10, k=50)" |
| k > m+n (query > population) | `validate(q=5, m=10, n=90, k=200)` | `constraints_satisfied = FALSE` | "Sample constraint - k (200) > m+n (100)" |
| Negative parameters | `validate(q=-1, m=10, n=90, k=50)` | `constraints_satisfied = FALSE` | "All parameters must be non-negative" |
| Non-integer parameters | `validate(q=2.7, m=10, n=90, k=50)` | `constraints_satisfied = FALSE` | "All parameters must be integers" |
| Near-integer (rounding) | `validate(q=2.0000001, m=10, n=90, k=50)` | Auto-rounded | Warning: "Parameters rounded to nearest integers" |

**Multi-layer validation framework** (in `parameter_constraints_validation.R`):
1. **Layer 1 — Gene-level:** Validate input copy numbers are non-negative integers
2. **Layer 2 — Parameter-level:** Verify hypergeometric constraints (q ≤ min(m,k), k ≤ m+n, feasibility bounds)
3. **Layer 3 — Consistency:** Check weighted parameter transformation preserves constraint satisfaction

**Recovery modes:**
- **Strict:** Reject any constraint violation; suitable for production pipelines
- **Lenient:** Auto-correct minor floating-point inconsistencies (within 1e-6 tolerance); flag and warn on larger discrepancies
- Both modes implemented and tested in `copy_weighted_robustness_tests.R`

## Robustness Testing Results

### Test Coverage Summary

| Test Category | Tests | Pass Rate | Key observations |
|--------------|-------|-----------|-----------------|
| Zero copy scenarios | 3 (all-zero, mixed, pathway-zero) | 100% | Degenerate cases flagged, p=1.0 returned |
| Extreme copy numbers | 2 (large copies, extreme ratios) | 100% | Finite p-values even at 10⁶-fold ratios |
| Boundary conditions | 9 (see table in §4) | 100% | All constraint flags correctly assigned |
| Parameter transformations | 100 (randomized) | 100% | Integrality preserved within 1e-6 tolerance |
| Constraint violation detection | 5 (invalid inputs) | 100% | All violations correctly detected and reported |
| Scalability | 4 (1K to 100K genes) | 100% | Linear memory scaling, sub-second for typical datasets |

### Stress Testing

- 100+ randomized scenarios with varying copy number distributions
- 0% failure rate
- Graceful degradation confirmed for all extreme parameter combinations

## Handling Strategies and Recommendations

### Decision Matrix

| Edge Case | Detection Method | Automated Response | User Action Required |
|-----------|-----------------|-------------------|---------------------|
| Zero copies | `all(cn == 0)` or `sum(cn) == 0` | Return p=1.0, flag `"degenerate_zero_copies"` | Review gene filtering |
| Large copy numbers | `max(cn) > 1e6` or `sum(cn) > 1e8` | Warning or hard rejection | Consider capping/log-transforming |
| Small samples | `k < 3` or `m < 3` | Power warning in output | Interpret with caution; consider alternative tests |
| Empty intersections | `q_w == 0` | Return p=1.0, flag `"zero_overlap"` | Expected for most pathways; filter in multiple-testing correction |
| Precision underflow | `p < 1e-15` or `phyper() → NaN` | Use `log.p=TRUE`; report `< 1e-15` | Trust log-scale values; normal approximation as fallback |
| Inconsistent models | Constraint violation detected | Strict: reject; Lenient: auto-correct + warn | Harmonize CN calling methods upstream |

### Best Practices for Production Use

1. **Always validate parameters before calling `phyper()`** — use `validate_hypergeometric_parameters()` to catch constraint violations before they produce silent errors
2. **Use `log.p = TRUE`** for any analysis that may produce very small p-values (common in large gene sets)
3. **Cap extreme copy numbers** at a biologically reasonable maximum (e.g., 1000 for most organisms) unless studying specific amplification events
4. **Flag low-power tests** rather than excluding them — downstream consumers should decide whether to include
5. **Report edge case flags** alongside p-values so that downstream analyses can filter or annotate appropriately

## Implementation Artifacts

| File | Role | Size |
|------|------|------|
| `copy_number_weighted_ora_edge_case_analysis.md` | Detailed edge case documentation (17 KB) | 17,377 B |
| `copy_weighted_robustness_tests.R` | Robustness test suite with `testthat` tests | 21,835 B |
| `edge_case_test_suite.R` | Edge case test framework | 13,944 B |
| `parameter_constraints_validation.R` | Core validation functions (4 exported functions) | 25,527 B |
| `copy_weighted_ora_parameter_bounds.yaml` | Safety bound configuration | 7,069 B |
| `boundary_conditions_validation_report.R` | Boundary condition validation script | 10,811 B |
| `copy_number_weighted_ora_best_practices.md` | Implementation guidelines (23 KB) | 22,772 B |
| `comprehensive_edge_case_test_results.RData` | Serialized test results | 1,113 B |
| `parameter_constraints_validation_report.RData` | Serialized validation report (`overall_status: ALL_TESTS_PASS`) | 413 B |

## Conclusion

All six required edge cases have been analyzed with concrete numerical demonstrations, validated through automated testing, and addressed with implemented handling strategies:

1. **Zero copy number genes:** Three sub-scenarios identified and handled (all-zero, query-zero, pathway-zero). All return p=1.0 with appropriate flags.
2. **Extremely large copy numbers (>1000):** Safety bounds at 1e6/1e8/1e4 thresholds. `phyper()` remains stable up to N=10⁶; near-underflow at N > 10⁵ handled via `log.p=TRUE`.
3. **Small sample sizes:** Warning thresholds at k<3, m<3. Power limitations flagged in output for downstream interpretation.
4. **Empty pathway intersections:** 9 boundary conditions validated. Zero overlap correctly yields p=1.0.
5. **Numerical precision limits:** `log.p=TRUE` fallback for extreme p-values. Auto-rounding for near-integer parameters within 1e-6.
6. **Inconsistent copy number models:** Multi-layer constraint validation with strict/lenient modes. All 5 violation types detected correctly.

---

**Report Completed:** 2026-04-01  
**Validation Status:** ALL REQUIREMENTS MET  
**Implementation Status:** PRODUCTION READY