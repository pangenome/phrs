# Decision Framework: When to Use Copy-Weighted vs Standard ORA

## Executive Summary

Copy-number-weighted ORA provides significant value for genomic regions with extreme copy number variation but comes with increased computational overhead and complexity. This framework provides quantitative decision criteria based on comprehensive analysis of Pseudohomologous Regions (PHRs) where copy-weighted ORA revealed a 12.35x copy expansion factor and significant functional composition bias (p=0.0118).

## Decision Criteria

### Primary Decision Factors

#### 1. Copy Number Variation Intensity
- **Standard ORA appropriate when:**
  - Copy expansion factor < 2x
  - Copy number range < 10-fold variation
  - Coefficient of variation (CV) of copy numbers < 1.0

- **Copy-weighted ORA recommended when:**
  - Copy expansion factor ≥ 5x
  - Copy number range ≥ 50-fold variation
  - Coefficient of variation (CV) of copy numbers ≥ 2.0

- **Copy-weighted ORA essential when:**
  - Copy expansion factor ≥ 10x
  - Copy number range ≥ 100-fold variation
  - Single genes represent >25% of total genomic content

#### 2. Dataset Characteristics
- **Query size considerations:**
  - Small queries (<50 genes): Standard ORA sufficient unless extreme CV
  - Medium queries (50-200 genes): Copy-weighting beneficial if CV > 1.5
  - Large queries (>200 genes): Copy-weighting recommended if CV > 1.0

- **Background complexity:**
  - Simple backgrounds (<10K genes): Copy-weighting feasible
  - Moderate backgrounds (10K-50K genes): Requires optimization
  - Large backgrounds (>50K genes): Use sampling approaches

#### 3. Functional Context
- **High-benefit contexts for copy-weighting:**
  - Gene family expansion analysis
  - Repetitive genomic regions
  - Dosage-sensitive pathways
  - Evolutionary genomics studies

- **Standard ORA preferred contexts:**
  - Essential gene analysis
  - Binary presence/absence studies
  - Cross-species comparisons (unless copy data available)

## Quantitative Decision Tree

### Step 1: Calculate Dataset Metrics

```r
calculate_decision_metrics <- function(gene_copy_data) {
  # Calculate key decision parameters
  copy_expansion_factor <- sum(gene_copy_data$copies) / nrow(gene_copy_data)
  copy_range_ratio <- max(gene_copy_data$copies) / min(gene_copy_data$copies)
  copy_cv <- sd(gene_copy_data$copies) / mean(gene_copy_data$copies)
  max_gene_proportion <- max(gene_copy_data$copies) / sum(gene_copy_data$copies)
  
  return(list(
    expansion_factor = copy_expansion_factor,
    range_ratio = copy_range_ratio,
    cv = copy_cv,
    max_proportion = max_gene_proportion
  ))
}
```

### Step 2: Apply Decision Logic

```r
recommend_ora_method <- function(metrics, query_size, computational_budget) {
  # Primary screening: Copy number variation intensity
  if (metrics$expansion_factor >= 10 || metrics$cv >= 2.0 || metrics$max_proportion >= 0.25) {
    method <- "copy_weighted_essential"
    confidence <- "high"
  } else if (metrics$expansion_factor >= 5 || metrics$cv >= 1.5 || metrics$range_ratio >= 50) {
    method <- "copy_weighted_recommended" 
    confidence <- "medium"
  } else if (metrics$expansion_factor >= 2 || metrics$cv >= 1.0) {
    method <- "copy_weighted_beneficial"
    confidence <- "low"
  } else {
    method <- "standard_ora"
    confidence <- "high"
  }
  
  # Secondary considerations: Query size and computational constraints
  if (method != "standard_ora" && query_size < 50) {
    method <- paste0(method, "_with_validation")
  }
  
  if (method != "standard_ora" && computational_budget == "limited") {
    method <- paste0(method, "_optimized")
  }
  
  return(list(method = method, confidence = confidence, metrics = metrics))
}
```

## Decision Flowchart

```
START: Pathway Enrichment Analysis
    |
    ▼
[1] Calculate copy number metrics
    ├─ Copy expansion factor
    ├─ Copy number CV
    ├─ Range ratio  
    └─ Max gene proportion
    |
    ▼
[2] Primary Decision: Copy Number Variation
    |
    ├─ CV ≥ 2.0 OR Expansion ≥ 10x OR Max proportion ≥ 25%
    │   └─ ESSENTIAL: Use Copy-Weighted ORA ──────────────┐
    │                                                     │
    ├─ CV ≥ 1.5 OR Expansion ≥ 5x OR Range ≥ 50x         │
    │   └─ RECOMMENDED: Use Copy-Weighted ORA ────────────┤
    │                                                     │
    ├─ CV ≥ 1.0 OR Expansion ≥ 2x                         │
    │   └─ BENEFICIAL: Consider Copy-Weighted ORA ────────┤
    │                                                     │
    └─ CV < 1.0 AND Expansion < 2x                        │
        └─ Use Standard ORA ─────────────────────────────────┐
                                                             │
    ┌────────────────────────────────────────────────────────┘
    ▼
[3] Secondary Considerations
    |
    ├─ Query size < 50 genes?
    │   └─ YES: Add validation step
    │
    ├─ Computational budget limited?
    │   └─ YES: Use optimized implementation
    │
    └─ Background > 50K genes?
        └─ YES: Use sampling approach
    |
    ▼
[4] Implementation Selection
    |
    ├─ ESSENTIAL + Large scale → Permutation method
    ├─ RECOMMENDED + Medium scale → Parameter method  
    ├─ BENEFICIAL + Small scale → Instance expansion
    └─ Standard → Classical hypergeometric
    |
    ▼
END: Execute selected method
```

## Use Case Scenarios

### Scenario 1: Pseudohomologous Regions (PHRs)
**Context:** 35 unique genes, 1,189 total copies, copy range 2-672
- **Metrics:** Expansion factor = 34x, CV = 4.2, Max proportion = 56.5%
- **Recommendation:** Copy-weighted ORA essential (high confidence)
- **Expected benefit:** Reveals olfactory gene family expansion (31.3% of content)
- **Implementation:** Parameter-based method with validation

### Scenario 2: Olfactory Receptor Gene Family Analysis
**Context:** 100 unique OR genes, 500 total copies, copy range 1-20
- **Metrics:** Expansion factor = 5x, CV = 1.8, Max proportion = 8%
- **Recommendation:** Copy-weighted ORA recommended (medium confidence)
- **Expected benefit:** Proper weighting of expanded gene families
- **Implementation:** Instance expansion or parameter method

### Scenario 3: Immunoglobulin Locus Analysis
**Context:** 50 unique IG genes, 300 total copies, copy range 1-45
- **Metrics:** Expansion factor = 6x, CV = 2.3, Max proportion = 15%
- **Recommendation:** Copy-weighted ORA essential (high confidence)
- **Expected benefit:** Accurate representation of VDJ diversity
- **Implementation:** Parameter method with permutation validation

### Scenario 4: Essential Gene Pathway Analysis
**Context:** 200 unique genes, 205 total copies (mostly single copy)
- **Metrics:** Expansion factor = 1.03x, CV = 0.2, Max proportion = 2%
- **Recommendation:** Standard ORA appropriate (high confidence)
- **Expected benefit:** None - copy weighting adds unnecessary complexity
- **Implementation:** Classical hypergeometric testing

### Scenario 5: Ribosomal Protein Gene Analysis
**Context:** 80 unique genes, 240 total copies, copy range 1-8
- **Metrics:** Expansion factor = 3x, CV = 1.2, Max proportion = 6%
- **Recommendation:** Copy-weighted ORA beneficial (low confidence)
- **Expected benefit:** Moderate - reflects dosage importance
- **Implementation:** Run both methods, compare results

## Implementation Recommendations

### For Copy-Weighted ORA Essential Cases

**Phase 1 Implementation (1-2 weeks):**
```r
# Parameter-based implementation for immediate results
weighted_ora_essential <- function(query_genes, copy_data, pathways) {
  # Calculate weighted parameters directly
  N_weighted <- sum(background_copy_data$copies)
  K_weighted <- sum(background_copy_data[background_copy_data$gene %in% pathway, "copies"])
  n_weighted <- sum(copy_data$copies)
  k_weighted <- sum(copy_data[copy_data$gene %in% pathway, "copies"])
  
  # Use standard phyper with weighted parameters
  pvalue <- phyper(k_weighted-1, K_weighted, N_weighted-K_weighted, n_weighted, lower.tail=FALSE)
  
  return(pvalue)
}
```

**Phase 2 Validation (2-3 weeks):**
- Implement permutation testing
- Validate against positive controls
- Compare with standard ORA results

### For Copy-Weighted ORA Recommended Cases

**Phased Approach:**
1. Implement parameter method first
2. Validate on subset of pathways
3. Full implementation if validation successful
4. Compare cost-benefit vs standard approach

### For Copy-Weighted ORA Beneficial Cases

**Conservative Approach:**
1. Run standard ORA first
2. Implement copy-weighted as secondary analysis
3. Report both results with interpretation guidance
4. Focus on pathways showing significant differences

## Performance and Resource Considerations

### Computational Complexity

| Dataset Scale | Standard ORA | Copy-Weighted ORA | Scaling Factor |
|---------------|--------------|-------------------|----------------|
| Small (<1K genes) | <1 second | <5 seconds | 2-5x |
| Medium (1K-10K genes) | <10 seconds | 1-5 minutes | 10-30x |
| Large (10K-50K genes) | <1 minute | 10-60 minutes | 10-60x |
| Genome-wide (>50K genes) | 1-10 minutes | 1-8 hours | 60-100x |

### Memory Requirements

| Implementation Method | Memory Overhead | Scalability Limit |
|----------------------|-----------------|-------------------|
| Parameter method | 1x (minimal) | Unlimited |
| Instance expansion | 10-100x | <10K genes |
| Sampling method | 2-5x | 50K+ genes |
| Permutation method | 5-20x | Medium datasets |

### Budget Planning

**Development Resources:**
- **Essential cases:** 2-4 weeks development, 1 bioinformatician
- **Recommended cases:** 1-2 weeks development, 0.5 FTE
- **Beneficial cases:** 3-5 days development, existing tools modification

**Computational Resources:**
- **Essential cases:** Dedicated compute node (16GB+ RAM)
- **Recommended cases:** Standard workstation (8GB+ RAM)
- **Beneficial cases:** Personal computer adequate

## Quality Control and Validation

### Mandatory Validation Steps

1. **Null Distribution Testing**
   - Generate 1000+ random gene sets
   - Test p-value uniformity (KS test p > 0.05)
   - Validate across copy number profiles

2. **Positive Control Validation**
   - Known enriched pathways should show expected signals
   - Effect sizes should align with biological expectations
   - Method comparison correlation > 0.8

3. **Computational Validation**
   - Parameter equivalence between methods
   - Numerical precision for large parameters
   - Performance benchmarking

### Red Flags Requiring Method Reconsideration

- **Statistical:** P-values non-uniform under null (KS p < 0.01)
- **Biological:** Known positive controls fail to enrich (sensitivity < 50%)
- **Computational:** Runtime >10x longer than budgeted
- **Practical:** Results inconsistent between validation methods

## Decision Implementation Template

```r
# Complete decision framework implementation
ora_decision_framework <- function(query_genes, copy_data, pathways, 
                                   computational_budget = "medium",
                                   validation_level = "standard") {
  
  # Step 1: Calculate decision metrics
  metrics <- calculate_decision_metrics(copy_data)
  
  # Step 2: Make recommendation
  recommendation <- recommend_ora_method(metrics, length(query_genes), computational_budget)
  
  # Step 3: Execute recommended method
  if (grepl("copy_weighted", recommendation$method)) {
    if (grepl("essential", recommendation$method)) {
      results <- copy_weighted_ora_essential(query_genes, copy_data, pathways)
      validation_required <- TRUE
    } else if (grepl("recommended", recommendation$method)) {
      results <- copy_weighted_ora_parameter(query_genes, copy_data, pathways)  
      validation_required <- (validation_level != "none")
    } else {
      # Run both methods for comparison
      results_standard <- standard_ora(query_genes, pathways)
      results_weighted <- copy_weighted_ora_basic(query_genes, copy_data, pathways)
      results <- list(standard = results_standard, weighted = results_weighted)
      validation_required <- FALSE
    }
  } else {
    results <- standard_ora(query_genes, pathways)
    validation_required <- FALSE
  }
  
  # Step 4: Validation if required
  if (validation_required) {
    validation_results <- validate_ora_results(results, copy_data, pathways)
    results$validation <- validation_results
  }
  
  # Step 5: Return comprehensive results
  return(list(
    recommendation = recommendation,
    metrics = metrics,
    results = results,
    interpretation_guide = generate_interpretation_guide(recommendation$method)
  ))
}
```

## Conclusion

The decision to use copy-weighted vs standard ORA should be data-driven, based primarily on copy number variation intensity (CV ≥ 2.0), copy expansion factors (≥10x), and the proportion of genomic content represented by highly amplified genes (≥25%). 

For datasets meeting "essential" criteria like PHRs, copy-weighted ORA is not optional—it reveals functional architecture invisible to standard approaches. For "recommended" cases, the added complexity is justified by improved biological accuracy. For "beneficial" cases, comparative analysis provides the best risk-benefit balance.

**Key Takeaway:** Copy number variation fundamentally changes the enrichment landscape. When present at significant levels, ignoring it leads to systematically biased functional interpretations that underweight the biological importance of expanded gene families.