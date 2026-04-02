# Copy-Number-Weighted ORA Investigation

## Executive Summary

Investigation of custom over-representation analysis (ORA) that weights by gene copy number rather than unique gene symbols, applied to PHR (Pseudohomologous Region) gene data.

## Data Context

**Current Dataset:**
- 35 unique gene symbols
- 1,189 total gene instances 
- Copy numbers range: 2-672 (extreme: MIR8078 with 672 copies)
- Most protein-coding genes: 2-39 copies
- Standard ORA treats each unique symbol equally (ignores copy number)

## Key Investigation Areas

### 1. Custom Background Construction

**Challenge:** Standard ORA backgrounds assume each gene appears once. Copy-number-weighted ORA requires backgrounds where genes appear proportional to their copy counts.

**Approaches:**

#### A) Instance-Based Background
- Represent each gene copy as separate entry in background
- Example: DUX4 (18 copies) appears 18 times in background
- **Pros:** Direct representation of biological reality
- **Cons:** Massive background inflation (genome-wide would have ~millions of entries)

#### B) Weighted Sampling Background  
- Sample genes proportional to copy number distribution
- Maintain manageable background size while preserving copy-number proportions
- **Implementation:** Bootstrap sampling with copy-number weights

#### C) Stratified Background by Copy Number
- Create separate backgrounds for different copy-number strata
- Test enrichment within each stratum independently
- **Advantage:** Controls for copy-number bias effects

### 2. Statistical Framework Modifications

**Core Issue:** Standard hypergeometric test assumes sampling without replacement from finite population with known composition.

#### Modified Hypergeometric Parameters:
- **Standard:** phyper(q, m, n, k)
  - q = genes in intersection
  - m = genes in pathway (background)  
  - n = genes not in pathway (background)
  - k = total genes in query

- **Copy-Weighted:** Need to redefine based on:
  - q = weighted intersection score
  - m = weighted pathway gene count
  - n = weighted non-pathway count
  - k = weighted query gene count

#### Alternative Statistical Approaches:
1. **Weighted Fisher's Exact Test**
2. **Permutation-Based Tests** with copy-aware shuffling
3. **Negative Binomial Models** accounting for copy-number overdispersion

### 3. R Implementation Strategy

#### Option A: Modified phyper() Parameters
```r
# Weighted query size
k_weighted <- sum(query_genes$copy_number)

# Weighted pathway composition
m_weighted <- sum(pathway_genes$copy_number) 
n_weighted <- sum(background_genes$copy_number) - m_weighted

# Test statistic: weighted overlap
q_weighted <- sum(intersect_genes$copy_number)

# Test
pval <- phyper(q_weighted-1, m_weighted, n_weighted, k_weighted, lower.tail=FALSE)
```

#### Option B: Instance-Based Approach  
```r
# Expand datasets by copy number
query_expanded <- rep(query_genes$symbol, query_genes$copy_number)
background_expanded <- rep(background_genes$symbol, background_genes$copy_number)

# Standard hypergeometric test on expanded sets
fisher.test(contingency_table)
```

### 4. Expected Changes from Standard ORA

#### Statistical Power Changes:
- **High-copy genes:** Dramatically increased weight in analysis
- **Single-copy genes:** Proportionally reduced influence
- **Overall:** Potential for different significantly enriched pathways

#### Interpretation Changes:
- **Standard ORA:** "X% of pathway genes are present"  
- **Copy-weighted ORA:** "X% of pathway gene instances are present"

#### Example Impact (PHR data):
- Standard: MIR8078 = 1 vote (2.9% of 35 genes)
- Weighted: MIR8078 = 672 votes (56.5% of 1,189 instances)

### 5. Implementation Complexity Assessment

#### Low Complexity (1-2 days):
- Modified phyper() with weighted parameters
- Instance expansion approach for small datasets

#### Medium Complexity (1-2 weeks):
- Custom permutation testing framework
- Weighted background construction algorithms
- Statistical validation against known controls

#### High Complexity (1-2 months):
- Tool integration with g:Profiler/similar platforms
- Scalable genome-wide background construction
- Comprehensive statistical framework with multiple test corrections

### 6. Validation Strategy

#### Controls Needed:
1. **Known pathways:** Test on pathways with well-characterized copy-number effects
2. **Null simulations:** Permutation tests with copy-aware shuffling  
3. **Comparison benchmarks:** Standard ORA results as baseline

#### Success Metrics:
- Statistically valid p-values under null
- Biologically meaningful enrichment differences
- Computational efficiency for genome-wide analysis

## Recommendations

### Immediate Implementation (Phase 1):
1. **Instance-based approach** for current PHR dataset (manageable scale)
2. **Modified phyper()** parameters with copy-number weights
3. **Direct comparison** to standard ORA results to quantify differences

### Advanced Development (Phase 2):
1. **Permutation testing framework** for robust statistical validation
2. **Scalable background construction** for genome-wide analysis
3. **Tool integration** for broader accessibility

### Critical Considerations:
- **Extreme outliers** (MIR8078: 672 copies) may dominate analysis
- **Biological interpretation:** Whether copy number should weight equally across gene types
- **Multiple testing:** How copy-number weighting affects FDR correction

## Expected Runtime Estimates

- **Phase 1 Implementation:** 2-5 days
- **Statistical validation:** 1-2 weeks  
- **Genome-wide scaling:** 2-4 weeks
- **Full tool development:** 2-3 months

## Next Steps

1. Implement instance-based weighted ORA for PHR dataset
2. Compare results with standard g:Profiler analysis
3. Develop permutation testing validation
4. Scale to genome-wide background construction

