# Copy-Number Weighted ORA Performance and Accuracy Analysis

**Date:** 2026-04-01  
**Task:** performance-and-accuracy  
**Objective:** Comprehensive analysis of performance benchmarks and accuracy validation for copy-number weighted Over-Representation Analysis

## Executive Summary

This report synthesizes performance benchmarks and accuracy validation results for copy-number weighted hypergeometric testing (weighted ORA) versus standard ORA approaches. The analysis reveals a critical trade-off: **substantial computational performance gains with severely compromised statistical accuracy**.

### Key Findings

- **Performance**: Copy-weighted approach shows 2.4-3.7× speedup for medium-large datasets with sub-linear scaling
- **Accuracy**: Weighted approach has **inflated Type I error** (4-6× above nominal) and **compromised FDR control** 
- **Mathematical Equivalence**: Perfect equivalence between parameter weighting and instance expansion methods
- **Scalability**: Direct weighted approach essential for large-scale analyses (>50K instances)
- **Recommendation**: Use weighted approach only when instance-level independence is justified

## Performance Analysis

### Runtime Performance Benchmarks

Comprehensive benchmarks across three dataset scales reveal clear performance advantages for the copy-weighted parameter approach:

| Dataset Scale | Direct Time (ms) | Expansion Time (ms) | Speedup Factor | Memory Reduction |
|---------------|------------------|---------------------|----------------|------------------|
| **Small** (PHR-scale) | 0.110 | 0.073 | 0.7× (expansion faster) | Modest |
| **Medium** | 0.214 | 0.510 | **2.4× faster** | ~34% |
| **Large** | 0.648 | 2.428 | **3.7× faster** | >50% |

#### Performance Insights

1. **Dataset Size Threshold**: Crossover point around 1,000-2,000 total instances where direct weighted parameters become more efficient
2. **Scaling Characteristics**:
   - **Direct approach**: Sub-linear scaling (Time ∝ n^0.48) - excellent scalability
   - **Instance expansion**: Near-linear scaling (Time ∝ n^0.95) - poor scalability
3. **Memory Efficiency**: Direct approach requires O(unique_genes) vs O(total_instances) for expansion

### Computational Complexity

The fundamental algorithmic differences create distinct complexity profiles:

```
Direct Weighted Approach:
- Time Complexity: O(n^0.48) - sub-linear, excellent scaling
- Memory Complexity: O(unique genes)
- Bottlenecks: Parameter validation, background statistics calculation

Instance Expansion Approach:  
- Time Complexity: O(n^0.95) - near-linear, poor scaling
- Memory Complexity: O(total instances)  
- Bottlenecks: Vector expansion via rep(), memory allocation overhead
```

### Scalability Analysis

Performance advantage increases dramatically with dataset size:

- **Small datasets (<2K instances)**: Either approach acceptable, slight preference for expansion
- **Medium datasets (2K-50K)**: Direct approach strongly preferred (2-4× speedup)  
- **Large datasets (>50K)**: Direct approach essential (>3× speedup, memory constraints)

### Multiple Pathway Testing

Vectorized approach provides consistent but modest speedup:
- Small dataset: 1.07× speedup
- Medium dataset: 1.05× speedup
- Diminishing returns due to small per-pathway computation overhead

## Accuracy Analysis

### Statistical Validity Assessment

Comprehensive validation across 8 test categories reveals **serious statistical limitations** for the copy-weighted approach under gene-level sampling:

#### Type I Error Control

| Alpha Level | Weighted Rate | Standard Rate | Weighted Controlled | Standard Controlled |
|-------------|---------------|---------------|-------------------|-------------------|
| 0.01 | 0.142 | 0.004 | ❌ **14× inflation** | ✅ Conservative |
| 0.05 | 0.225 | 0.018 | ❌ **4.5× inflation** | ✅ Controlled |
| 0.10 | 0.263 | 0.071 | ❌ **2.6× inflation** | ✅ Controlled |

**Critical Finding**: Weighted phyper() does not control Type I error at any conventional significance level.

#### Copy Number Magnitude Effect

Type I error inflation correlates directly with copy number magnitude:

| Copy Number Scenario | Type I (weighted) | Type I (standard) | Mean CN | Inflation Factor |
|---------------------|------------------|------------------|---------|------------------|
| All CN=1 | 0.019 | 0.016 | 1.0 | 1.2× |
| All CN=5 | 0.153 | 0.022 | 5.0 | 7.0× |
| CN 1-3 (uniform) | 0.104 | 0.019 | 2.0 | 5.5× |
| CN 1-20 (uniform) | 0.287 | 0.017 | 10.5 | 16.9× |

#### Null Distribution Properties

Validation of null p-value distributions:

- **Weighted phyper**: Non-uniform distribution (KS p-value ~0), anti-conservative
- **Standard phyper**: Conservative but well-controlled (Type I < α)
- **Mathematical equivalence**: 100% parameter matching between weighted and expansion approaches

### False Discovery Rate Analysis

Testing 100 pathways with 10 truly enriched (BH correction at FDR < 0.05):

| Method | Mean Discoveries | True Positives | False Positives | Actual FDR | Target FDR |
|--------|-----------------|----------------|-----------------|------------|------------|
| **Weighted** | 22.4 | 7.4 | 15.0 | **66.1%** | 5% |
| **Standard** | 0.3 | 0.3 | 0.1 | 2.6% | 5% |

**Critical Issue**: BH correction fails to control FDR for weighted phyper() due to anti-conservative null p-values.

### Power Analysis

ROC analysis with controlled false positive rates:

| False Positive Rate | True Positive Rate (Weighted) | True Positive Rate (Standard) |
|--------------------|------------------------------|------------------------------|
| 0.01 | 0.626 | **0.740** |
| 0.05 | 0.818 | **0.850** |
| 0.10 | 0.920 | **0.930** |

**Finding**: Standard phyper() achieves higher true positive rates at every controlled false positive rate, demonstrating better calibrated sensitivity.

## Root Cause Analysis: Statistical Anti-Conservative Behavior

### Theoretical Explanation

The anti-conservative behavior has a clear mathematical basis:

1. **Model Assumption**: Hypergeometric test models draws as independent samples from gene *instances*
2. **Actual Sampling**: In ORA, genes are selected as units, creating **clustering** of instances
3. **Clustering Effect**: Instances of same gene are perfectly correlated, violating independence assumption
4. **Statistical Impact**: Clustering inflates effective sample size, reducing variance and producing systematically small p-values

This is analogous to the "design effect" in survey statistics where cluster sampling requires inflation corrections.

### Copy Number Scaling Effect

Even with uniform copy numbers, scaling all parameters by constant c produces:
- Hypergeometric distribution with smaller relative variance than original
- Systematically deflated p-values
- Anti-conservative behavior proportional to copy number magnitude

## Performance vs Accuracy Trade-offs

### Trade-off Matrix

| Aspect | Copy-Weighted Approach | Standard Approach |
|--------|----------------------|------------------|
| **Computational Speed** | **2.4-3.7× faster** | Baseline |
| **Memory Usage** | **34-97% reduction** | Baseline |
| **Scalability** | **Sub-linear (n^0.48)** | Linear+ (n^0.95) |
| **Type I Error Control** | ❌ **Severely inflated** | ✅ Well-controlled |
| **FDR Control** | ❌ **Ineffective** | ✅ Reliable |
| **Power at Controlled FPR** | ❌ Lower | ✅ Higher |
| **Mathematical Validity** | ✅ Exact equivalence | ✅ Standard |

### Decision Framework

**Use Copy-Weighted Approach When:**
- Instance-level independence assumption holds
- Individual gene copies can independently be in query set
- Selection is truly at instance level (e.g., genomic region overlaps)
- Computational performance is critical and statistical inflation is acceptable

**Use Standard Approach When:**
- Genes are selected as units (most biological analyses)
- Type I error control is critical
- FDR guarantees are needed
- Statistical validity outweighs computational efficiency

## Computational Complexity Documentation

### Algorithm Complexity Analysis

#### Direct Weighted Parameter Method
```
Time Complexity: O(G log G + P) where G = genes, P = pathways
- Gene filtering and merging: O(G log G)
- Parameter calculation: O(G)
- phyper() call: O(1)
- Multiple pathways: O(P × G)

Space Complexity: O(G)
- Gene dataframes: O(G)
- No instance expansion required
```

#### Instance Expansion Method
```
Time Complexity: O(I + P) where I = total instances, P = pathways
- Vector expansion: O(I)
- Set intersections: O(I)
- phyper() call: O(1)

Space Complexity: O(I)
- Expanded gene vectors: O(I)
- Memory scales with copy number magnitude
```

#### Scaling Comparison
For datasets with mean copy number c:
- Direct approach: O(G) independent of copy numbers
- Expansion approach: O(c × G) scales linearly with copy numbers

### Performance Bottlenecks

**Primary Bottlenecks Identified:**

1. **Instance Expansion**: `rep()` function becomes expensive with large copy numbers
2. **Memory Allocation**: Large expanded vectors require significant overhead
3. **Custom Implementation**: R-level loops create substantial computational overhead

**Optimization Opportunities:**
1. Parallel processing for multiple pathways
2. Background statistics caching for repeated queries
3. C++ implementation for very large-scale applications

## Recommendations

### Production Implementation Strategy

**Tiered Approach Based on Use Case:**

1. **High-Accuracy Applications** (recommended default):
   ```r
   # Use standard phyper() for statistical validity
   standard_ora <- function(query_genes, pathway_genes, background_genes) {
     # Standard hypergeometric without copy weighting
   }
   ```

2. **Performance-Critical Applications** (with statistical caveats):
   ```r
   # Use weighted approach with explicit warnings
   weighted_ora_with_warnings <- function(query_df, pathway_genes, background_df) {
     warning("Weighted ORA has inflated Type I error - use with caution")
     return(weighted_hypergeometric_test(query_df, pathway_genes, background_df))
   }
   ```

3. **Hybrid Validation Approach**:
   ```r
   # Provide both results for comparison
   comprehensive_ora <- function(query_df, pathway_genes, background_df) {
     standard_result <- standard_ora(query_df$gene, pathway_genes, background_df$gene)
     weighted_result <- weighted_hypergeometric_test(query_df, pathway_genes, background_df)
     
     return(list(
       standard = standard_result,
       weighted = weighted_result,
       agreement = standard_result$significant && weighted_result$significant
     ))
   }
   ```

### Quality Assurance Guidelines

1. **Development Phase**:
   - Use instance expansion for validation and testing
   - Implement comprehensive parameter validation
   - Include edge case handling for extreme copy numbers

2. **Production Phase**:
   - Default to standard approach for statistical validity
   - Use weighted approach only with explicit justification
   - Always report both weighted and standard results when feasible
   - Include warnings about Type I error inflation

3. **Statistical Mitigation**:
   - Implement permutation-based p-values for copy-weighted analyses
   - Apply conservative corrections to compensate for inflation
   - Use empirical null distributions based on data structure

## Conclusion

The copy-number weighted ORA approach presents a **classic performance-accuracy trade-off**:

### Performance Advantages
- **2.4-3.7× speedup** for medium-large datasets
- **Sub-linear scaling** essential for large-scale analyses
- **Substantial memory reduction** (34-97%)
- **Perfect mathematical equivalence** with instance expansion

### Statistical Limitations
- **Severely inflated Type I error** (4-6× above nominal)
- **Compromised FDR control** (66% actual vs 5% target)
- **Anti-conservative behavior** scaling with copy number magnitude
- **Lower power at controlled false positive rates**

### Final Recommendation

**Use standard ORA as the default approach** for most biological analyses where genes are selected as units. Reserve copy-weighted ORA for specialized cases where:

1. Instance-level independence is justified
2. Computational performance is critical
3. Statistical inflation is acceptable or can be mitigated
4. Results are validated through permutation or other robust methods

The substantial performance gains do not justify the severe statistical limitations for most enrichment analyses in genomics and systems biology.

---

## Technical Specifications

**Benchmark Environment:**
- Platform: Linux 4.19.0-27-amd64
- R Version: 4.3.0+
- Validation: 2000+ simulations per test condition
- Reproducibility: Fixed random seeds throughout

**Test Coverage:**
- 8 statistical validation categories  
- 3 performance benchmark scales
- Multiple copy number distributions
- Various pathway and query sizes
- Both null and enriched scenarios