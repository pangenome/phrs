# Copy-Number Weighted ORA: Comprehensive Methodology Synthesis

**Date:** 2026-04-01  
**Document Type:** Methodology Synthesis  
**Scope:** Integration of mathematical framework, statistical validation, and practical implementation

## Executive Summary

This document synthesizes the comprehensive research into copy-number weighted over-representation analysis (ORA), integrating mathematical formulations, statistical validation findings, and practical implementation considerations. The research reveals both the theoretical elegance and significant practical limitations of copy-weighted hypergeometric tests.

### Key Findings

1. **Mathematical Framework**: Copy-weighted ORA represents a well-defined mathematical transformation of standard hypergeometric parameters, scaling gene counts to gene instance counts based on copy numbers.

2. **Critical Statistical Limitation**: Under gene-level sampling (the typical ORA scenario), weighted hypergeometric tests exhibit anti-conservative behavior with Type I error rates 4-6× above nominal levels due to clustering effects.

3. **Implementation Complexity**: While mathematically straightforward, copy-weighted methods require sophisticated statistical frameworks, permutation-based approaches, and careful validation to achieve reliable results.

4. **Practical Recommendations**: Standard unweighted ORA remains the preferred approach for most applications, with weighted methods reserved for specific use cases where instance-level independence assumptions are justified.

## Mathematical Foundation

### Theoretical Framework

The copy-weighted hypergeometric distribution transforms standard ORA parameters through gene copy number scaling:

**Standard Parameters:**
- N: Total background genes
- K: Genes in pathway  
- n: Query genes
- k: Overlap genes

**Weighted Transformation:**
```
N_w = Σ(copy_number_i) over all background genes
K_w = Σ(copy_number_j) over pathway genes  
n_w = Σ(copy_number_k) over query genes
k_w = Σ(copy_number_l) over overlap genes
```

**Modified Distribution:**
```
P(X = k_w | N_w, K_w, n_w) = C(K_w, k_w) × C(N_w-K_w, n_w-k_w) / C(N_w, n_w)
```

### Parameter Scaling Relationships

The transformation involves scaling factors that capture the relative copy number distributions:

```
S_bg = mean background copy number = N_w / N
S_path = mean pathway copy number = K_w / K  
S_query = mean query copy number = n_w / n
S_overlap = mean overlap copy number = k_w / k

Effective scaling = (S_path × S_query) / (S_bg × S_overlap)
```

### Mathematical Validation

Comprehensive testing (200 test cases) confirmed 100% mathematical equivalence between the weighted parameter approach and full instance expansion, validating the mathematical correctness of the transformation.

## Statistical Properties and Validation

### Null Distribution Behavior

**Critical Finding**: Under gene-level sampling, copy-weighted hypergeometric tests produce anti-conservative p-value distributions:

| Method | Type I Error (α=0.05) | Mean p-value | KS Test Result |
|--------|----------------------|--------------|----------------|
| Weighted phyper | 0.217 | 0.558 | Non-uniform |
| Standard phyper | 0.024 | 0.599 | Conservative |

### Root Cause Analysis

The anti-conservative behavior stems from a fundamental mismatch between model assumptions and actual sampling:

1. **Model Assumption**: The hypergeometric test assumes independent sampling of gene instances
2. **Actual Sampling**: Genes are selected as units, bringing all copies together (clustering)
3. **Effect**: Clustering inflates effective sample size, reducing variance and producing systematically small p-values

This is analogous to the "design effect" in survey statistics, where cluster sampling violates independence assumptions.

### Copy Number Magnitude Effects

Type I error inflation correlates strongly with copy number magnitude:

| Copy Number Scenario | Mean CN | Type I Error (Weighted) | Type I Error (Standard) |
|---------------------|---------|-------------------------|-------------------------|
| All CN=1 | 1.0 | 0.019 | 0.016 |
| CN 1-3 | 2.0 | 0.104 | 0.019 |
| CN 1-8 | 4.7 | 0.198 | 0.021 |
| CN 1-20 | 10.5 | 0.287 | 0.017 |

### False Discovery Rate Control

**Critical Limitation**: Benjamini-Hochberg FDR correction fails with weighted p-values:

| Method | Target FDR | Actual FDR | Mean Discoveries | False Positives |
|--------|------------|------------|------------------|-----------------|
| Weighted | 0.05 | **0.661** | 22.4 | 15.0 |
| Standard | 0.05 | 0.026 | 0.3 | 0.1 |

The 66% actual FDR occurs because BH correction assumes uniform or conservative null p-values, which weighted tests violate.

### Power Analysis

While weighted tests show apparent power advantages at low enrichment levels, this reflects their anti-conservative nature rather than genuine statistical power:

| Enrichment Factor | Power (Weighted) | Power (Standard) | Interpretation |
|-------------------|------------------|------------------|----------------|
| 1.0 (null) | 0.210 | 0.024 | Weighted inflated |
| 1.5 | 0.480 | 0.170 | Weighted artifact |
| 3.0 | 0.954 | 0.860 | Both high power |

ROC analysis at controlled false positive rates demonstrates that standard phyper achieves higher true positive rates at every FPR level, confirming better-calibrated performance.

## Advanced Statistical Methods

### Permutation-Based Approaches

When copy number effects are important but gene-level selection applies, permutation methods provide the gold standard:

#### Gene-Level Permutation
```R
# Shuffle gene labels while preserving copy number structure
permutation_test_gene_level <- function(query_genes, pathways, copy_numbers, n_perm = 10000) {
  null_distribution <- replicate(n_perm, {
    shuffled_genes <- sample(names(copy_numbers), size = length(query_genes))
    calculate_enrichment_stats(shuffled_genes, pathways, copy_numbers)
  })
  
  # Calculate empirical p-values
  return(calculate_empirical_pvalues(observed_stats, null_distribution))
}
```

#### Copy-Weighted Permutation
```R
# Sample genes proportional to their copy numbers
permutation_test_weighted <- function(query_genes, pathways, copy_numbers, n_perm = 10000) {
  total_query_copies <- sum(copy_numbers[query_genes])
  
  null_distribution <- replicate(n_perm, {
    null_genes <- sample(names(copy_numbers), 
                        size = total_query_copies,
                        replace = TRUE,
                        prob = copy_numbers / sum(copy_numbers))
    calculate_pathway_overlaps(null_genes, pathways, copy_numbers)
  })
  
  return(calculate_empirical_pvalues(observed_overlaps, null_distribution))
}
```

### Multiple Testing Correction Adaptations

#### Adaptive FDR Methods
For copy-weighted data, standard FDR methods require modification:

1. **Empirical Null Estimation**: Account for non-uniform null p-value distributions
2. **Stratified Correction**: Apply correction within copy number strata
3. **Weight-Adjusted FDR**: Modify correction factors based on copy number heterogeneity

#### Permutation-Based Multiple Testing
- **Max-T Procedure**: Control family-wise error rate using permutation maxima
- **Step-Down FDR**: Apply Benjamini-Hochberg procedure with permutation-derived null distributions

## Computational Considerations

### Numerical Precision Challenges

Copy-weighted parameters can become very large:
- Standard N ≤ 50,000
- Weighted N_w can exceed 2,000,000
- Large integers may cause combinatorial overflow

**Solutions:**
1. **High-Precision Arithmetic**: Use logarithmic calculations
2. **Normal Approximation**: For extreme parameters (N_w > 1,000,000)
3. **Stirling's Approximation**: For very large combinatorial terms

### Memory-Efficient Implementation

```R
# Sparse matrix representations for pathway data
create_sparse_pathway_matrix <- function(pathways, genes, copy_numbers) {
  library(Matrix)
  
  # Create sparse binary pathway matrix
  pathway_matrix <- sparseMatrix(i = pathway_indices, j = gene_indices, x = 1)
  
  # Apply copy number weights
  copy_weights <- copy_numbers[genes]
  weighted_matrix <- pathway_matrix %*% Diagonal(x = copy_weights)
  
  return(weighted_matrix)
}
```

### Performance Optimization

1. **Caching**: Pre-compute background statistics
2. **Vectorization**: Batch process multiple pathways
3. **Parallel Processing**: Distribute permutation computations
4. **Early Termination**: Stop permutations when sufficient precision achieved

## Implementation Decision Framework

### Method Selection Guidelines

**Use Standard (Unweighted) Hypergeometric When:**
- Genes are selected as units (typical ORA scenario)
- Type I error control is critical
- FDR guarantees are required
- Copy numbers are used for annotation only

**Use Weighted Hypergeometric When:**
- Instance-level selection genuinely applies
- Individual genomic copies can independently belong to query set
- Exploratory analysis where sensitivity matters more than specificity
- **Always pair with standard results for comparison**

**Use Permutation-Based Methods When:**
- Copy number effects are important AND gene-level selection applies
- Highest statistical rigor is required
- Computational resources permit extensive validation

### Quality Control Framework

#### Pre-Analysis Validation
1. **Copy Number Distribution Assessment**
   - Check for extreme outliers (>1000 copies)
   - Validate copy number coefficient of variation
   - Assess correlation between copy numbers and pathway membership

2. **Background Gene Set Validation**
   - Ensure completeness and appropriate scope
   - Check for systematic biases in gene selection
   - Validate copy number consistency

#### Analysis Validation
1. **Statistical Property Testing**
   - Null distribution validation via simulation
   - Type I error rate assessment
   - Power analysis with known controls

2. **Method Comparison**
   - Run both standard and weighted methods
   - Compare results against biological expectations
   - Assess concordance between methods

#### Post-Analysis Validation
1. **Biological Plausibility Assessment**
   - Evaluate enriched pathways for biological relevance
   - Check consistency with known biology
   - Assess copy number effects on pathway interpretation

2. **Sensitivity Analysis**
   - Test robustness to copy number estimation errors
   - Evaluate stability across different background definitions
   - Assess impact of extreme copy number exclusions

## Pathway-Specific Considerations

### High-Copy Pathways
**Examples**: Olfactory receptors, immunoglobulins, histone genes
- **Effect**: Increased statistical power in weighted analysis
- **Consideration**: May dominate enrichment signals
- **Recommendation**: Monitor for over-interpretation of high-copy pathway enrichments

### Low-Copy Pathways
**Examples**: Transcription factors, essential genes
- **Effect**: Decreased power compared to standard ORA
- **Consideration**: May miss true enrichments in weighted analysis
- **Recommendation**: Use standard ORA as primary method for low-copy pathways

### Mixed-Copy Pathways
**Examples**: Metabolic pathways with diverse gene families
- **Effect**: Complex, unpredictable effects
- **Consideration**: Weighted and standard results may differ substantially
- **Recommendation**: Report both results and discuss biological implications

## Reporting Standards and Best Practices

### Required Documentation
1. **Method Specification**
   - State whether gene-level or instance-level selection applies
   - Document copy number estimation methodology
   - Specify statistical approach and parameters used

2. **Statistical Reporting**
   - Report both standard and weighted results when different
   - Document multiple testing correction method
   - Acknowledge anti-conservative behavior if using weighted tests

3. **Validation Documentation**
   - Include null distribution validation results
   - Report positive/negative control performance
   - Document sensitivity analysis findings

### Recommended Analysis Pipeline

#### Stage 1: Standard ORA
```R
# Primary analysis using established methods
standard_results <- lapply(pathways, function(pathway) {
  fisher_test(query_genes, pathway$genes, background_genes)
})

# Apply FDR correction
standard_results$p_adj <- p.adjust(standard_results$p_value, method = "fdr")
```

#### Stage 2: Copy-Weighted Analysis (If Justified)
```R
# Secondary analysis with copy number weighting
weighted_results <- lapply(pathways, function(pathway) {
  weighted_hypergeometric_test(query_df, pathway$genes, background_df)
})

# Apply appropriate multiple testing correction
weighted_results$p_adj <- apply_weighted_correction(weighted_results$p_value)
```

#### Stage 3: Permutation Validation (If Required)
```R
# Validation analysis using permutation methods
permutation_results <- lapply(pathways, function(pathway) {
  permutation_test_gene_level(query_genes, pathway$genes, copy_numbers)
})
```

#### Stage 4: Result Integration and Interpretation
```R
# Combine results for comparative analysis
integrated_results <- merge_analysis_results(standard_results, 
                                            weighted_results, 
                                            permutation_results)

# Generate interpretation framework
interpretation <- interpret_multi_method_results(integrated_results)
```

## Use Case Applications

### PHR (Pseudohomologous Region) Analysis
- **Context**: Extreme copy number variation (up to 672× for some genes)
- **Recommendation**: Use permutation methods as primary approach
- **Rationale**: Standard methods miss copy number effects; weighted methods are anti-conservative

### Genome-Wide Expression Analysis
- **Context**: Moderate copy number variation, large gene sets
- **Recommendation**: Standard ORA with copy number as annotation
- **Rationale**: Gene-level selection applies; copy effects typically modest

### Targeted Gene Panel Analysis
- **Context**: Small gene sets, variable copy numbers
- **Recommendation**: Report both methods with extensive validation
- **Rationale**: Small sample sizes make statistical properties critical

### Clinical Diagnostic Applications
- **Context**: High accuracy requirements, established thresholds needed
- **Recommendation**: Conservative permutation-based approaches
- **Rationale**: Type I error control is paramount in clinical settings

## Future Directions and Research Needs

### Methodological Improvements
1. **Calibrated Weighted Tests**: Develop correction factors for clustering effects
2. **Adaptive Methods**: Create algorithms that automatically select appropriate approaches
3. **Bayesian Frameworks**: Incorporate copy number uncertainty into analysis

### Computational Advances
1. **Scalable Permutation**: Develop efficient algorithms for large-scale permutation testing
2. **GPU Acceleration**: Leverage parallel computing for intensive calculations
3. **Memory Optimization**: Create memory-efficient implementations for large datasets

### Biological Integration
1. **Copy Number Context**: Better integration of copy number mechanisms into interpretation
2. **Pathway Weighting**: Develop methods to weight pathways based on their copy number sensitivity
3. **Multi-Modal Analysis**: Integrate copy number effects with other genomic features

## Conclusions and Recommendations

### Summary of Key Findings

1. **Mathematical Validity**: Copy-weighted hypergeometric tests are mathematically well-defined and correctly implement the intended transformation from gene counts to gene instance counts.

2. **Statistical Limitations**: Under typical gene-level sampling scenarios, weighted tests exhibit severe anti-conservative behavior (Type I error 4-6× nominal) due to clustering effects that violate independence assumptions.

3. **FDR Control Failure**: Standard multiple testing corrections fail with weighted p-values, leading to actual FDR rates of ~66% versus target 5%.

4. **Method Performance**: At controlled false positive rates, standard hypergeometric tests consistently outperform weighted tests in terms of true positive rates.

### Primary Recommendations

1. **Default to Standard ORA**: Use unweighted hypergeometric tests as the primary method for pathway enrichment analysis.

2. **Limited Weighted Applications**: Reserve copy-weighted methods for scenarios where instance-level independence is genuinely justified.

3. **Permutation Gold Standard**: When copy number effects are important and gene-level selection applies, use permutation-based approaches.

4. **Comparative Reporting**: Always report both standard and weighted results when copy numbers vary significantly.

5. **Comprehensive Validation**: Implement thorough statistical validation including null distribution testing, positive/negative controls, and sensitivity analysis.

### Implementation Priority

1. **High Priority**: Implement robust standard ORA with copy number annotation
2. **Medium Priority**: Develop permutation-based methods for copy-sensitive analyses
3. **Low Priority**: Implement weighted parametric methods with extensive caveats

### Final Assessment

Copy-number weighted ORA represents an intellectually appealing approach to incorporating gene dosage effects into pathway analysis. However, the research demonstrates that mathematical elegance does not guarantee statistical validity under realistic sampling conditions. The anti-conservative behavior under gene-level sampling, combined with the failure of standard multiple testing corrections, makes weighted hypergeometric tests inappropriate for most practical applications.

The standard (unweighted) hypergeometric test, supplemented with copy number information for interpretation and annotation, remains the most statistically sound approach for pathway enrichment analysis. When copy number effects are genuinely important to the biological question, permutation-based methods provide the necessary statistical rigor, albeit at increased computational cost.

This synthesis represents a comprehensive evaluation of copy-weighted ORA methodology, integrating theoretical, statistical, and practical considerations to provide evidence-based guidance for pathway enrichment analysis in the presence of gene copy number variation.

## Source Artifacts

This synthesis integrates findings from:
- `mathematical_framework_copy_weighted_ora.md` - Mathematical foundations and theoretical framework
- `advanced_statistical_methods_copy_weighted_ora.md` - Permutation methods and advanced approaches  
- `weighted_phyper_statistical_validation_report.md` - Comprehensive statistical validation results
- `statistical_best_practices_weighted_ora.md` - Implementation guidelines and best practices
- `statistical_framework_implementation_guide.md` - Practical implementation decision framework
- Multiple R implementation files and validation scripts

---

*Methodology synthesis completed 2026-04-01 by methodology-synthesis-documentation task. Integrates comprehensive research findings into unified methodological framework.*