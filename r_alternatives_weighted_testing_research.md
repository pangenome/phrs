# R Alternatives to phyper() for Weighted Hypergeometric Testing: Comprehensive Research

## Executive Summary

This research evaluates alternative R implementations to `phyper()` for copy-number weighted hypergeometric testing in Over-Representation Analysis (ORA). While parameter-modified `phyper()` provides an efficient solution, several alternative approaches exist in both base R and specialized packages that offer different trade-offs in terms of flexibility, performance, and statistical features.

**Key Findings:**
1. **Base R alternatives** (dhyper, fisher.test) provide mathematical equivalence but limited weighted support
2. **BioConductor packages** (GOstats, clusterProfiler) offer enrichment workflows but lack direct weighted implementations 
3. **Specialized CRAN packages** provide some weighted statistical methods but not specifically for hypergeometric testing
4. **Permutation-based approaches** offer greatest flexibility for custom weighting schemes

## Research Methodology

### Investigation Approach
- Systematic review of R documentation and package repositories
- Analysis of existing implementations and their APIs
- Evaluation against weighted testing requirements
- Performance and usability assessments

### Evaluation Criteria
1. **Weighted Support**: Native support for instance-based or weighted testing
2. **Computational Efficiency**: Memory and runtime performance characteristics  
3. **Statistical Robustness**: Adherence to proper statistical theory
4. **API Usability**: Integration ease and development workflow
5. **Community Support**: Package maintenance and ecosystem integration

---

## Base R Hypergeometric Functions

### 1. dhyper() - Density Function

**Function Signature:**
```r
dhyper(x, m, n, k, log = FALSE)
```

**Purpose:** Computes probability mass function values for hypergeometric distribution.

**Weighted Testing Application:**
```r
# Calculate exact probability of observed weighted overlap
weighted_probability <- dhyper(q_weighted, m_weighted, n_weighted, k_weighted)

# Can be used for exact p-value calculation via summation
exact_pvalue <- sum(dhyper(q_weighted:min(m_weighted, k_weighted), 
                          m_weighted, n_weighted, k_weighted))
```

**Evaluation:**
- ✅ **Weighted Support**: Full compatibility with weighted parameters
- ✅ **Mathematical Equivalence**: Identical to phyper() results when summed appropriately
- ⚠️ **Computational Efficiency**: Slower for p-value calculation (requires summation loop)
- ✅ **Statistical Robustness**: Direct hypergeometric theory implementation
- ⚠️ **API Usability**: Requires manual p-value computation

**Use Case Recommendation:** Best for exact probability calculations or when full probability distribution is needed.

### 2. qhyper() - Quantile Function

**Function Signature:**
```r
qhyper(p, m, n, k, lower.tail = TRUE, log.p = FALSE)
```

**Purpose:** Computes quantiles of hypergeometric distribution.

**Weighted Testing Application:**
```r
# Find critical value for significance threshold
alpha <- 0.05
critical_value <- qhyper(1 - alpha, m_weighted, n_weighted, k_weighted, 
                        lower.tail = FALSE)

# Compare observed overlap to critical threshold
is_significant <- q_weighted >= critical_value
```

**Evaluation:**
- ✅ **Weighted Support**: Compatible with weighted parameters
- ✅ **Computational Efficiency**: Direct quantile computation
- ✅ **Statistical Robustness**: Proper hypergeometric quantiles
- ⚠️ **API Usability**: Indirect approach for p-value testing
- ⚡ **Special Use**: Valuable for power analysis and threshold determination

**Use Case Recommendation:** Ideal for determining significance thresholds or power analysis.

### 3. rhyper() - Random Number Generation

**Function Signature:**
```r
rhyper(nn, m, n, k)
```

**Purpose:** Generates random samples from hypergeometric distribution.

**Weighted Testing Application:**
```r
# Generate null distribution for empirical p-value calculation
n_simulations <- 10000
null_samples <- rhyper(n_simulations, m_weighted, n_weighted, k_weighted)

# Calculate empirical p-value
empirical_pvalue <- mean(null_samples >= q_weighted)

# Bootstrap confidence intervals
bootstrap_ci <- quantile(null_samples, c(0.025, 0.975))
```

**Evaluation:**
- ✅ **Weighted Support**: Generates weighted hypergeometric samples
- ⚠️ **Computational Efficiency**: Slower than analytical methods but useful for validation
- ✅ **Statistical Robustness**: Exact sampling from weighted distribution
- ✅ **API Usability**: Simple interface for simulation-based testing
- ⚡ **Special Value**: Essential for method validation and null distribution visualization

**Use Case Recommendation:** Critical for validating weighted implementations and generating confidence intervals.

---

## Fisher's Exact Test Family

### 4. fisher.test() - Standard Implementation

**Function Signature:**
```r
fisher.test(x, y = NULL, workspace = 200000, hybrid = FALSE, 
           control = list(), or = 1, alternative = "two.sided",
           conf.int = TRUE, conf.level = 0.95, simulate.p.value = FALSE, B = 2000)
```

**Purpose:** Performs Fisher's exact test for independence in contingency tables.

**Weighted Testing Application:**
```r
# Convert hypergeometric parameters to 2x2 contingency table
contingency_matrix <- matrix(c(
  q_weighted,                           # overlap instances
  k_weighted - q_weighted,              # query non-overlap instances  
  m_weighted - q_weighted,              # pathway non-query instances
  n_weighted - (k_weighted - q_weighted) # background remainder
), nrow = 2, byrow = TRUE)

result <- fisher.test(contingency_matrix, alternative = "greater")
```

**Evaluation:**
- ✅ **Weighted Support**: Accepts weighted contingency tables
- ✅ **Mathematical Equivalence**: Identical p-values to phyper() for one-sided tests
- ⚠️ **Computational Efficiency**: Network algorithm can be slower for large tables
- ✅ **Statistical Robustness**: Well-established exact test
- ✅ **API Usability**: Rich output including confidence intervals and odds ratios
- ⚡ **Added Value**: Provides effect size estimates (odds ratio) with confidence intervals

**Use Case Recommendation:** Excellent when effect size estimates and confidence intervals are needed.

### 5. Weighted Fisher Variants

**Package Options:**
- `exactRankTests::fisher.exact()` - Enhanced implementation
- `GeneNet::fisher.test.4d()` - 4-dimensional Fisher test
- `coin::fisher_test()` - Permutation-based Fisher test

**Example with coin Package:**
```r
library(coin)

# Create contingency table from weighted parameters
tab <- table(factor(c(rep("query", k_weighted), rep("background", sum(m_weighted, n_weighted) - k_weighted))),
             factor(c(rep("pathway", q_weighted), rep("other", k_weighted - q_weighted),
                     rep("pathway", m_weighted - q_weighted), rep("other", n_weighted - (k_weighted - q_weighted)))))

result <- fisher_test(tab, alternative = "greater", distribution = "exact")
```

**Evaluation:**
- ✅ **Weighted Support**: Enhanced flexibility for weighted tables
- ⚠️ **Computational Efficiency**: Variable, depends on implementation
- ✅ **Statistical Robustness**: Multiple exact algorithms available
- ⚠️ **API Usability**: Package-specific interfaces vary
- ⚡ **Special Features**: Some offer permutation-based alternatives

---

## BioConductor Enrichment Analysis Packages

### 6. GOstats Package

**Primary Functions:**
```r
library(GOstats)
# hyperGTest() - Hypergeometric testing for GO terms
# fisherTest() - Fisher's exact test wrapper
```

**Current Implementation:**
```r
# Standard GOstats workflow (unweighted)
params <- new("GOHyperGParams",
              geneIds = query_genes,
              universeGeneIds = background_genes, 
              annotation = annotation_package,
              ontology = "BP",
              pvalueCutoff = 0.05,
              conditional = FALSE,
              testDirection = "over")

result <- hyperGTest(params)
```

**Weighted Modification Potential:**
```r
# Custom weighted implementation using GOstats infrastructure
weighted_hypergeometric_gostats <- function(query_df, background_df, go_annotation) {
  # Extract GO terms for genes
  query_go <- merge(query_df, go_annotation, by = "gene")
  background_go <- merge(background_df, go_annotation, by = "gene")
  
  # Group by GO term and calculate weighted parameters
  results <- ddply(query_go, .(go_term), function(term_data) {
    pathway_genes <- unique(term_data$gene)
    
    # Calculate weighted parameters (using approach from phyper research)
    k_weighted <- sum(query_df$copy_number)
    q_weighted <- sum(query_df[query_df$gene %in% pathway_genes, "copy_number"])
    m_weighted <- sum(background_df[background_df$gene %in% pathway_genes, "copy_number"])
    n_weighted <- sum(background_df$copy_number) - m_weighted
    
    # Weighted hypergeometric test
    pvalue <- phyper(q_weighted-1, m_weighted, n_weighted, k_weighted, lower.tail = FALSE)
    
    return(data.frame(go_term = term_data$go_term[1],
                     pvalue = pvalue,
                     overlap_instances = q_weighted,
                     pathway_instances = m_weighted))
  })
  
  return(results)
}
```

**Evaluation:**
- ❌ **Native Weighted Support**: No built-in weighted testing
- ✅ **Extension Potential**: Framework allows custom implementations
- ✅ **Statistical Robustness**: Well-tested hypergeometric implementations
- ✅ **API Usability**: Clean object-oriented interface
- ✅ **Ecosystem Integration**: Excellent Bioconductor integration

**Use Case Recommendation:** Best for GO/pathway analysis when extended with custom weighted implementations.

### 7. clusterProfiler Package

**Primary Functions:**
```r
library(clusterProfiler)
# enrichGO() - GO enrichment analysis
# enrichKEGG() - KEGG pathway enrichment
# GSEA() - Gene Set Enrichment Analysis
```

**Current Implementation:**
```r
# Standard clusterProfiler workflow
ego <- enrichGO(gene = query_genes,
                universe = background_genes,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)
```

**Weighted Modification Potential:**
```r
# Custom weighted enrichment using clusterProfiler infrastructure
weighted_enrich_go <- function(query_df, background_df, orgdb) {
  # Extract GO annotations
  go_annotations <- AnnotationDbi::select(orgdb, keys = background_df$gene, 
                                         columns = c("GOALL"), keytype = "SYMBOL")
  
  # Calculate weighted enrichment for each GO term
  go_results <- go_annotations %>%
    group_by(GOALL) %>%
    do({
      pathway_genes <- unique(.$SYMBOL)
      
      # Weighted hypergeometric parameters
      k_weighted <- sum(query_df$copy_number)
      q_weighted <- sum(query_df[query_df$gene %in% pathway_genes, "copy_number"])
      m_weighted <- sum(background_df[background_df$gene %in% pathway_genes, "copy_number"])
      n_weighted <- sum(background_df$copy_number) - m_weighted
      
      # Test and format results
      pvalue <- phyper(q_weighted-1, m_weighted, n_weighted, k_weighted, lower.tail = FALSE)
      
      data.frame(
        ID = .$GOALL[1],
        pvalue = pvalue,
        p.adjust = p.adjust(pvalue, method = "BH"),  # Apply after all terms computed
        GeneRatio = paste(q_weighted, k_weighted, sep = "/"),
        BgRatio = paste(m_weighted, m_weighted + n_weighted, sep = "/")
      )
    })
  
  return(go_results)
}
```

**Evaluation:**
- ❌ **Native Weighted Support**: No built-in weighted testing
- ✅ **Extension Potential**: Modular design enables custom implementations
- ✅ **Visualization**: Excellent plotting and visualization tools
- ✅ **Statistical Features**: Multiple test correction, effect sizes
- ✅ **API Usability**: Modern tidyverse-compatible interface

**Use Case Recommendation:** Ideal for comprehensive pathway analysis workflows when extended with weighted methods.

---

## Specialized CRAN Packages

### 8. Weighted Statistical Testing Packages

#### a) `survey` Package - Weighted Analysis

**Functions:**
```r
library(survey)
# svychisq() - Weighted chi-square tests
# svyttest() - Weighted t-tests
```

**Application to Hypergeometric Testing:**
```r
# Convert hypergeometric to weighted chi-square framework
weighted_contingency_test <- function(query_df, pathway_genes, background_df) {
  # Create expanded survey design
  query_expanded <- query_df[rep(seq_len(nrow(query_df)), query_df$copy_number), ]
  background_expanded <- background_df[rep(seq_len(nrow(background_df)), background_df$copy_number), ]
  
  # Create contingency data
  all_data <- rbind(
    data.frame(gene = query_expanded$gene, group = "query"),
    data.frame(gene = background_expanded$gene, group = "background")
  )
  
  all_data$in_pathway <- all_data$gene %in% pathway_genes
  
  # Weighted chi-square test
  design <- svydesign(ids = ~1, data = all_data)
  result <- svychisq(~group + in_pathway, design = design)
  
  return(result)
}
```

**Evaluation:**
- ⚠️ **Weighted Support**: Designed for survey weights, not copy number weights
- ❌ **Hypergeometric Specificity**: Not designed for hypergeometric testing
- ⚠️ **Computational Efficiency**: Requires data expansion
- ✅ **Statistical Robustness**: Well-tested weighted statistical methods
- ❌ **Direct Applicability**: Requires significant adaptation

#### b) `Hmisc` Package - Statistical Functions

**Relevant Functions:**
```r
library(Hmisc)
# wtd.chi.sq() - Weighted chi-square test
# wtd.t.test() - Weighted t-test
```

**Application:**
```r
# Weighted chi-square for independence testing
weighted_chi_square_test <- function(query_df, pathway_genes, background_df) {
  # Create weighted contingency table
  all_genes <- union(query_df$gene, background_df$gene)
  
  contingency_data <- data.frame(
    gene = all_genes,
    query_weight = ifelse(all_genes %in% query_df$gene, 
                         query_df$copy_number[match(all_genes, query_df$gene)], 0),
    background_weight = background_df$copy_number[match(all_genes, background_df$gene)],
    in_pathway = all_genes %in% pathway_genes
  )
  
  # Weighted chi-square test
  result <- wtd.chi.sq(contingency_data$in_pathway, 
                      contingency_data$query_weight > 0,
                      weight = contingency_data$background_weight)
  
  return(result)
}
```

**Evaluation:**
- ⚠️ **Weighted Support**: Limited weighted contingency table support
- ❌ **Hypergeometric Theory**: Uses chi-square approximation, not exact hypergeometric
- ✅ **Computational Efficiency**: Reasonably efficient
- ⚠️ **Statistical Robustness**: Approximation rather than exact test
- ❌ **Direct Applicability**: Different statistical framework

### 9. Permutation-Based Approaches

#### a) `coin` Package - Conditional Inference

**Functions:**
```r
library(coin)
# independence_test() - General permutation test framework
# chisq_test() - Permutation-based chi-square
```

**Weighted Hypergeometric Implementation:**
```r
permutation_weighted_hypergeometric <- function(query_df, pathway_genes, background_df, 
                                               n_permutations = 10000) {
  
  # Observed test statistic (weighted overlap)
  observed_overlap <- sum(query_df[query_df$gene %in% pathway_genes, "copy_number"])
  
  # Permutation testing
  permuted_overlaps <- replicate(n_permutations, {
    # Randomly sample query-sized set from background (respecting copy numbers)
    background_expanded <- rep(background_df$gene, background_df$copy_number)
    query_size <- sum(query_df$copy_number)
    
    permuted_query <- sample(background_expanded, size = query_size, replace = FALSE)
    permuted_overlap <- sum(table(permuted_query[permuted_query %in% pathway_genes]))
    
    return(permuted_overlap)
  })
  
  # Calculate empirical p-value
  pvalue <- mean(permuted_overlaps >= observed_overlap)
  
  return(list(
    pvalue = pvalue,
    observed = observed_overlap,
    permutation_distribution = permuted_overlaps
  ))
}
```

**Evaluation:**
- ✅ **Weighted Support**: Fully flexible for any weighting scheme
- ⚠️ **Computational Efficiency**: Slower than analytical methods
- ✅ **Statistical Robustness**: Exact permutation-based inference
- ✅ **API Usability**: Clear permutation framework
- ⚡ **Maximum Flexibility**: Handles any complex weighting or dependency structure

#### b) Custom Bootstrap/Permutation Implementation

**Resampling-Based Weighted Testing:**
```r
bootstrap_weighted_hypergeometric <- function(query_df, pathway_genes, background_df,
                                             method = c("permutation", "bootstrap"),
                                             n_resamples = 10000) {
  
  method <- match.arg(method)
  
  observed_overlap <- sum(query_df[query_df$gene %in% pathway_genes, "copy_number"])
  query_size <- sum(query_df$copy_number)
  
  if (method == "permutation") {
    # Permutation: resample query from background without replacement
    null_distribution <- replicate(n_resamples, {
      background_expanded <- rep(background_df$gene, background_df$copy_number)
      sampled_genes <- sample(background_expanded, size = query_size, replace = FALSE)
      sum(table(sampled_genes[sampled_genes %in% pathway_genes]))
    })
  } else {
    # Bootstrap: resample with replacement
    null_distribution <- replicate(n_resamples, {
      background_expanded <- rep(background_df$gene, background_df$copy_number)
      sampled_genes <- sample(background_expanded, size = query_size, replace = TRUE)
      sum(table(sampled_genes[sampled_genes %in% pathway_genes]))
    })
  }
  
  # Calculate p-value and confidence intervals
  pvalue <- mean(null_distribution >= observed_overlap)
  ci_lower <- quantile(null_distribution, 0.025)
  ci_upper <- quantile(null_distribution, 0.975)
  
  return(list(
    pvalue = pvalue,
    observed = observed_overlap,
    expected = mean(null_distribution),
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    null_distribution = null_distribution
  ))
}
```

**Evaluation:**
- ✅ **Weighted Support**: Complete flexibility for weighting schemes
- ❌ **Computational Efficiency**: Computationally intensive
- ✅ **Statistical Robustness**: Non-parametric, assumption-free
- ✅ **Confidence Intervals**: Natural CI computation from resampling
- ⚡ **Diagnostic Value**: Full null distribution available for examination

---

## Performance and Efficiency Analysis

### Computational Complexity Comparison

| Method | Time Complexity | Memory Complexity | Suitability |
|--------|----------------|-------------------|-------------|
| phyper() modified | O(1) | O(unique_genes) | ⭐⭐⭐⭐⭐ |
| dhyper() summation | O(range) | O(unique_genes) | ⭐⭐⭐⭐ |
| fisher.test() | O(min(m,n,k)!) | O(unique_genes) | ⭐⭐⭐ |
| Permutation tests | O(n_perm × k) | O(total_instances) | ⭐⭐ |
| GOstats/clusterProfiler | O(n_pathways) | O(unique_genes) | ⭐⭐⭐⭐ |

### Memory Usage Benchmarks

**Test Data:** PHR dataset (35 genes, 1189 instances)

```r
benchmark_memory_usage <- function() {
  methods <- list(
    phyper_modified = function() phyper(56-1, 800, 1999200, 1189, lower.tail=FALSE),
    dhyper_summation = function() sum(dhyper(56:min(800,1189), 800, 1999200, 1189)),
    fisher_exact = function() fisher.test(matrix(c(56, 1133, 744, 1998456), nrow=2))$p.value,
    permutation = function() permutation_weighted_hypergeometric(query_df, pathway_genes, background_df, 1000)$pvalue
  )
  
  memory_usage <- sapply(methods, function(method) {
    gc()
    mem_before <- sum(gc()[,2])
    result <- method()
    mem_after <- sum(gc()[,2])
    return(mem_after - mem_before)
  })
  
  return(memory_usage)
}
```

**Expected Results:**
- phyper() modified: ~0.1 MB
- dhyper() summation: ~0.2 MB  
- fisher.test(): ~0.5 MB
- Permutation: ~50 MB

### Runtime Performance Benchmarks

```r
library(microbenchmark)

runtime_comparison <- microbenchmark(
  phyper_modified = phyper(56-1, 800, 1999200, 1189, lower.tail=FALSE),
  dhyper_sum = sum(dhyper(56:min(800,1189), 800, 1999200, 1189)),
  fisher_exact = fisher.test(matrix(c(56, 1133, 744, 1998456), nrow=2))$p.value,
  times = 1000
)
```

**Expected Performance Ranking:**
1. phyper() modified: ~0.1 ms
2. dhyper() summation: ~2 ms
3. fisher.test(): ~10 ms
4. Permutation (1000 reps): ~100 ms

---

## Statistical Robustness Assessment

### Null Distribution Validation

**Test Framework:**
```r
validate_method_null_distribution <- function(test_method, n_simulations = 1000) {
  null_pvalues <- replicate(n_simulations, {
    # Generate null query (no enrichment)
    null_query <- sample(background_df$gene, size = nrow(query_df), replace = FALSE)
    null_query_df <- data.frame(
      gene = null_query,
      copy_number = sample(query_df$copy_number)  # Permute copy numbers
    )
    
    # Run test method
    test_method(null_query_df, pathway_genes, background_df)
  })
  
  # Test uniformity of p-values under null
  ks_result <- ks.test(null_pvalues, punif)
  
  return(list(
    pvalues = null_pvalues,
    ks_pvalue = ks_result$p.value,
    uniform = ks_result$p.value > 0.05,
    type_i_error_005 = mean(null_pvalues < 0.05),
    type_i_error_001 = mean(null_pvalues < 0.01)
  ))
}
```

**Expected Results for Well-Calibrated Methods:**
- KS test p-value > 0.05 (uniformity not rejected)
- Type I error ≈ α (e.g., ~0.05 for α=0.05)

### Power Analysis Comparison

**Test Framework:**
```r
compare_statistical_power <- function(methods, effect_sizes, n_reps = 100) {
  power_results <- expand.grid(
    method = names(methods),
    effect_size = effect_sizes,
    power = NA
  )
  
  for (i in seq_len(nrow(power_results))) {
    method_name <- as.character(power_results$method[i])
    effect_size <- power_results$effect_size[i]
    
    # Generate queries with specified effect size
    significant_tests <- replicate(n_reps, {
      enriched_query <- generate_enriched_query(effect_size)
      pvalue <- methods[[method_name]](enriched_query, pathway_genes, background_df)
      return(pvalue < 0.05)
    })
    
    power_results$power[i] <- mean(significant_tests)
  }
  
  return(power_results)
}
```

---

## API Usability and Integration

### Function Interface Comparison

**Standardized API Design:**
```r
# Proposed unified interface for all methods
weighted_enrichment_test <- function(query_df, pathway_genes, background_df, 
                                   method = c("phyper", "dhyper", "fisher", "permutation"),
                                   n_permutations = 10000,
                                   alternative = c("greater", "less", "two.sided"),
                                   conf.level = 0.95) {
  
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  
  # Input validation
  validate_input_data(query_df, pathway_genes, background_df)
  
  # Method dispatch
  result <- switch(method,
    phyper = weighted_phyper_test(query_df, pathway_genes, background_df),
    dhyper = weighted_dhyper_test(query_df, pathway_genes, background_df),
    fisher = weighted_fisher_test(query_df, pathway_genes, background_df),
    permutation = weighted_permutation_test(query_df, pathway_genes, background_df, n_permutations)
  )
  
  # Standardized output format
  structure(
    list(
      method = method,
      pvalue = result$pvalue,
      statistic = result$statistic,
      alternative = alternative,
      conf.int = result$conf.int,
      estimate = result$estimate,
      data.name = deparse(substitute(query_df))
    ),
    class = "weighted_enrichment_test"
  )
}
```

### Integration with Existing Workflows

**BioConductor Integration:**
```r
# Extend existing enrichment objects
setClass("WeightedEnrichmentResult",
         contains = "GeneSetCollection",
         slots = c(
           pvalues = "numeric",
           weighted_statistics = "list",
           copy_number_model = "data.frame"
         ))

# Method for existing generic functions
setMethod("summary", "WeightedEnrichmentResult",
          function(object) {
            cat("Weighted Enrichment Analysis Results\n")
            cat("Copy-number model:", nrow(object@copy_number_model), "genes\n")
            cat("Significant pathways (p < 0.05):", sum(object@pvalues < 0.05), "\n")
          })
```

---

## Gap Analysis and Missing Functionality

### Current Limitations

1. **Native Weighted Support**:
   - No R packages provide native weighted hypergeometric testing
   - Existing packages require manual parameter calculation or data expansion

2. **Computational Optimization**:
   - Most alternatives less efficient than parameter-modified phyper()
   - Limited vectorization for multiple pathway testing

3. **Statistical Features**:
   - Few packages provide confidence intervals for weighted enrichment
   - Limited support for effect size estimation in weighted context

4. **User Experience**:
   - No standardized interface across weighted testing approaches
   - Limited documentation for weighted applications

### Recommended Developments

#### 1. Weighted Hypergeometric Package

**Package Structure:**
```r
# Proposed R package: weightedHypergeometric

#' Weighted hypergeometric testing
#' @export
whyper <- function(q, m, n, k, weights = NULL, lower.tail = TRUE) {
  if (is.null(weights)) {
    return(phyper(q, m, n, k, lower.tail = lower.tail))
  }
  
  # Implement weighted parameter calculation
  weighted_params <- calculate_weighted_parameters(q, m, n, k, weights)
  return(phyper(weighted_params$q, weighted_params$m, 
                weighted_params$n, weighted_params$k, lower.tail = lower.tail))
}

#' Weighted enrichment testing for pathways
#' @export
weighted_enrichment <- function(query_df, pathways_list, background_df, 
                               method = "hypergeometric", correction = "BH") {
  # Implementation with multiple pathway support
}
```

#### 2. BioConductor Extension Package

**Package Structure:**
```r
# Proposed BioConductor package: weightedGSEA

#' S4 class for weighted gene sets
setClass("WeightedGeneSet", 
         contains = "GeneSet",
         slots = c(copy_weights = "numeric"))

#' Weighted over-representation analysis
#' @export
weightedORA <- function(geneList, geneSet, universe, weights) {
  # Implementation using efficient weighted hypergeometric testing
}
```

---

## Integration Strategy Recommendations

### 1. Short-term: Extend Existing Packages

**Strategy:** Contribute weighted functionality to established packages

**Priority Targets:**
1. **clusterProfiler**: Add weighted ORA functions
2. **GOstats**: Extend hyperGTest for weighted testing  
3. **fgsea**: Add weighted GSEA algorithms

**Implementation Approach:**
```r
# Example contribution to clusterProfiler
enrichGO_weighted <- function(gene, weights, universe, universe_weights, ...) {
  # Create weighted gene/universe dataframes
  query_df <- data.frame(gene = gene, copy_number = weights)
  background_df <- data.frame(gene = universe, copy_number = universe_weights)
  
  # Use existing clusterProfiler infrastructure with weighted parameters
  # ...
}
```

### 2. Medium-term: Dedicated Weighted Package

**Package Goals:**
- Unified interface for weighted statistical testing
- Multiple algorithm implementations (exact, approximate, permutation)
- Comprehensive validation and diagnostic tools
- Integration with major enrichment analysis workflows

**Core Functions:**
```r
# Core weighted testing functions
weighted.hypergeometric.test()
weighted.fisher.test() 
weighted.permutation.test()

# Enrichment analysis workflows
weighted.pathway.enrichment()
weighted.go.enrichment()
weighted.kegg.enrichment()

# Diagnostic and validation tools
validate.weighted.method()
power.analysis.weighted()
null.distribution.check()
```

### 3. Long-term: Statistical Method Developments

**Research Priorities:**
1. **Exact algorithms** for large-scale weighted testing
2. **Approximate methods** with theoretical guarantees
3. **Bayesian approaches** incorporating copy number uncertainty
4. **Machine learning integration** for copy number modeling

---

## Conclusions and Recommendations

### Primary Recommendations

1. **For immediate use**: Continue with parameter-modified `phyper()` approach
   - Most efficient and mathematically sound
   - Easy integration with existing workflows
   - Well-validated statistical properties

2. **For enhanced functionality**: Develop wrapper using `fisher.test()`
   - Provides confidence intervals and effect size estimates
   - Maintains exact statistical properties
   - Compatible with existing R statistical frameworks

3. **For flexible research**: Implement permutation-based testing
   - Maximum flexibility for complex weighting schemes
   - Natural confidence interval computation
   - Valuable for method validation and diagnostics

### Method Selection Guide

**Use `phyper()` modified when:**
- High-throughput analysis (many pathways/gene sets)
- Computational efficiency is critical
- Standard hypergeometric assumptions are met
- Integration with existing ORA workflows needed

**Use `fisher.test()` when:**
- Effect size estimates required
- Confidence intervals needed
- Small number of tests (computational cost acceptable)
- Compatibility with contingency table frameworks desired

**Use permutation testing when:**
- Complex dependency structures exist
- Non-standard weighting schemes required
- Method validation and diagnostics needed
- Computational time is not limiting

**Use specialized packages when:**
- Comprehensive pathway analysis workflows needed
- Rich visualization and reporting required
- Integration with annotation databases essential
- Multiple testing correction and meta-analysis planned

### Implementation Priorities

1. **Immediate**: Standardized wrapper functions for common use cases
2. **Short-term**: Contribution to existing BioConductor packages  
3. **Medium-term**: Dedicated weighted enrichment analysis package
4. **Long-term**: Advanced statistical method development

This comprehensive evaluation provides the foundation for informed selection and implementation of weighted hypergeometric testing approaches in R, supporting the continued development of copy-number-aware enrichment analysis methodologies.