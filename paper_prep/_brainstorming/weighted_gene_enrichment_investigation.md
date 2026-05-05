# Weighted Gene Enrichment Methods Investigation

## Executive Summary

This investigation evaluates weighted gene enrichment approaches for analyzing Pseudohomologous Region (PHR) regions where genes have variable copy numbers (2-39 copies per gene in the current dataset). Traditional binary enrichment methods treat all genes equally, potentially undervaluing highly duplicated genes that may have greater functional impact.

**Key Findings:**
- **Best Method**: GSEA with pre-ranked lists using copy numbers as ranking weights
- **Available Tools**: fgsea (R), GSEA desktop/web, GSEApy (Python)
- **Implementation Complexity**: Moderate (requires background construction)
- **Runtime Estimate**: <5 minutes for 1,189 gene copies
- **Alternative**: Weighted hypergeometric tests with custom implementations

---

## Method 1: Gene Set Enrichment Analysis (GSEA) with Pre-ranked Lists

### Overview
GSEA with pre-ranked gene lists is the most established method for copy-number weighted enrichment analysis. Instead of using gene expression fold-changes for ranking, gene copy numbers serve as the ranking metric.

### Implementation Approach
```
1. Create ranked gene list: Gene_Name, Copy_Count (descending order)
2. Run GSEA against pathway databases (GO, KEGG, etc.)
3. Use permutation testing to assess significance
4. Report normalized enrichment scores (NES) and FDR q-values
```

### Tools Available

#### 1. fgsea (R Package)
- **Pros**: Fast algorithm, well-maintained, integrates with Bioconductor
- **Cons**: Requires R environment setup
- **Runtime**: ~30 seconds for 1,189 genes vs 15,000 pathway gene sets
- **Installation**: `install.packages("fgsea")`

```r
library(fgsea)
library(msigdbr)

# Example implementation
gene_ranks <- setNames(copy_counts$total_copies, copy_counts$gene_name)
pathways <- msigdbr(species = "Homo sapiens", category = "C5") # GO terms

fgsea_results <- fgsea(pathways = pathways, 
                       stats = gene_ranks,
                       minSize = 15,
                       maxSize = 500)
```

#### 2. GSEApy (Python)
- **Pros**: Python ecosystem, command line interface
- **Cons**: Can be slower than fgsea
- **Runtime**: ~2-5 minutes for similar analysis
- **Installation**: `pip install gseapy`

```python
import gseapy as gp

# Pre-ranked GSEA
pre_res = gp.prerank(rnk=gene_ranks,  # pandas Series
                     gene_sets='GO_Biological_Process_2023',
                     outdir='output',
                     permutation_num=1000)
```

#### 3. GSEA Desktop/Web
- **Pros**: User-friendly interface, well-established
- **Cons**: Java dependency, slower for large analyses
- **Runtime**: ~5-10 minutes with web interface

### Background Construction Requirements
- **Gene Universe**: All protein-coding genes in human genome (~20,000)
- **Copy Number Assignment**: 
  - PHR genes: Use actual copy counts (2-39)
  - Non-PHR genes: Assign count of 1 (single copy assumed)
- **Ranking**: Sort by copy count (descending) to prioritize high-copy genes

### Statistical Considerations
- **Permutation Strategy**: Permute gene rankings while preserving pathway structure
- **Multiple Testing**: FDR correction across all tested pathways
- **Effect Size**: Normalized Enrichment Score (NES) accounts for pathway size

---

## Method 2: Weighted Hypergeometric Tests

### Overview
Extension of Fisher's exact test where genes contribute weighted probabilities based on copy number rather than binary presence/absence.

### Mathematical Framework
Traditional hypergeometric test:
```
P(X ≥ k) where:
- N = total genes in universe
- K = genes in pathway 
- n = genes in query set
- k = overlap between query and pathway
```

Weighted version:
```
P(W ≥ w) where:
- W = sum of weights for overlapping genes
- w = observed weighted overlap
- Weights = copy counts
```

### Implementation Options

#### 1. BioConductor GOstats with Custom Weights
- **Status**: Requires custom modification
- **Complexity**: High (need to modify internal functions)
- **Maintenance**: Risk of breaking with updates

#### 2. Custom R Implementation
```r
weighted_hypergeometric <- function(query_genes, query_weights, 
                                   pathway_genes, background_weights) {
  # Calculate observed weighted overlap
  overlap_genes <- intersect(query_genes, pathway_genes)
  observed_weight <- sum(query_weights[overlap_genes])
  
  # Permutation-based p-value
  null_dist <- replicate(10000, {
    perm_genes <- sample(names(background_weights), length(query_genes),
                        prob = background_weights, replace = FALSE)
    perm_overlap <- intersect(perm_genes, pathway_genes)
    sum(background_weights[perm_overlap])
  })
  
  p_value <- mean(null_dist >= observed_weight)
  return(p_value)
}
```

#### 3. Python Implementation with scipy.stats
- **Approach**: Use weighted sampling for permutation tests
- **Runtime**: ~10-30 seconds per pathway (slower than GSEA)

### Limitations
- **No Standard Implementation**: Requires custom coding
- **Computational Cost**: Higher than GSEA for large pathway collections
- **Statistical Properties**: Less well-characterized than GSEA

---

## Method 3: Background Weighting Strategies

### Genome-wide Copy Number Background

#### Option A: Assume Single Copy for Non-PHR Genes
- **Rationale**: Most genes are single copy; PHR duplications are exceptional
- **Implementation**: Set background weight = 1 for all non-PHR genes
- **Pros**: Conservative, computationally simple
- **Cons**: May miss subtle copy number effects

#### Option B: Incorporate Known Duplication Data
- **Data Sources**: 
  - Segmental Duplications Database
  - DGV (Database of Genomic Variants)  
  - Copy Number Variation studies
- **Implementation**: Use empirical copy counts from population data
- **Pros**: More accurate background
- **Cons**: Data complexity, population variation

#### Option C: Model-based Background
- **Approach**: Estimate copy counts based on genomic features
- **Features**: Gene length, chromosome location, sequence composition
- **Pros**: Systematic approach
- **Cons**: Model uncertainty, computational complexity

### Recommended Strategy
Use **Option A** initially (single copy background) for simplicity and interpretability. This creates a conservative test where only PHR-specific duplications contribute to enrichment.

---

## Method 4: Alternative Weighted Approaches

### 1. Weighted Gene Ontology Analysis (wGSA)
- **Tool**: Not currently maintained
- **Method**: Direct weighting of GO term enrichment
- **Status**: Deprecated, use GSEA instead

### 2. GREAT (Genomic Regions Enrichment of Annotations Tool)
- **Application**: Regional enrichment with distance weighting
- **Relevance**: Could weight genes by distance from PHR centers
- **Implementation**: Web tool or rGREAT R package

### 3. Camera (Competitive Gene Set Test)
- **Package**: edgeR/limma
- **Application**: Originally for RNA-seq, adaptable to copy counts
- **Advantage**: Accounts for inter-gene correlation

```r
library(limma)
# Adapt for copy counts
design <- model.matrix(~copy_counts$is_phr_gene)
camera_res <- camera(copy_matrix, pathway_indices, design)
```

---

## Tool Availability and Implementation Complexity

### Complexity Rankings (1-5 scale, 5 = most complex)

| Method | Setup | Implementation | Maintenance | Total |
|--------|-------|---------------|-------------|-------|
| fgsea GSEA | 2 | 2 | 1 | **5** |
| GSEApy | 2 | 2 | 2 | **6** |
| Weighted Hypergeometric | 3 | 4 | 4 | **11** |
| Camera | 2 | 3 | 2 | **7** |
| Custom R/Python | 4 | 5 | 5 | **14** |

### Resource Requirements
- **Memory**: <1GB for current dataset size
- **CPU**: Single-core sufficient for GSEA methods
- **Dependencies**: R/Python environment with bioinformatics packages
- **Data Storage**: ~100MB for pathway databases

---

## Runtime Estimates for 1,189 Gene Copies

### GSEA Methods
- **fgsea**: 30-60 seconds (15,000 pathway tests)
- **GSEApy**: 2-5 minutes (same tests)
- **GSEA Web**: 5-10 minutes (manual interface overhead)

### Hypergeometric Methods  
- **Per pathway**: 10-30 seconds with 10,000 permutations
- **Full analysis**: 2-6 hours for 15,000 pathways
- **Optimization**: Pre-compute permutation null distributions

### Bottlenecks
1. **Pathway Database Size**: Scales linearly with number of gene sets
2. **Permutation Count**: Higher permutations = better p-value precision
3. **Background Size**: Larger gene universe = slower permutation sampling

---

## Validation Requirements

### Method Validation (Addressing Task Requirements)

✅ **At least 2 weighted enrichment methods evaluated**:
1. GSEA with pre-ranked copy counts (recommended)
2. Weighted hypergeometric tests (custom implementation)

✅ **Implementation requirements documented**:
- Software dependencies specified
- Code examples provided
- Runtime estimates given

✅ **Background weighting approach specified**:
- Single copy assumption for non-PHR genes
- Alternative strategies outlined
- Recommended conservative approach

### Additional Validation Steps
1. **Positive Controls**: Test with known duplicated pathways (e.g., olfactory receptors)
2. **Negative Controls**: Verify no enrichment with random gene sets
3. **Sensitivity Analysis**: Compare results across different copy count thresholds
4. **Comparison with Binary**: Quantify differences from traditional enrichment

---

## Recommendations

### Primary Recommendation: GSEA with Pre-ranked Lists
1. **Tool**: fgsea R package for performance, GSEApy for Python integration
2. **Ranking**: Use copy counts directly (2-39 for PHR genes, 1 for background)
3. **Pathways**: Start with GO Biological Process and Molecular Function
4. **Parameters**: min_size=15, max_size=500, 1000 permutations

### Implementation Timeline
- **Week 1**: Setup fgsea environment, prepare ranked gene lists
- **Week 2**: Run initial analysis, validate against current binary results  
- **Week 3**: Refine parameters, expand to additional pathway databases
- **Week 4**: Documentation and comparison reporting

### Future Extensions
1. **Multi-level Analysis**: Incorporate both copy count and expression data
2. **Pathway Topology**: Weight based on pathway position/centrality
3. **Cross-species**: Compare enrichment patterns across species with different duplication patterns
4. **Dynamic Analysis**: Track copy number evolution and functional shifts

---

## Code Templates and Quick Start

### fgsea Implementation Template
```r
# Install packages if needed
if (!require("fgsea")) install.packages("fgsea")
if (!require("msigdbr")) install.packages("msigdbr")

# Load libraries
library(fgsea)
library(msigdbr)
library(dplyr)

# Prepare gene rankings (replace with actual data)
gene_copy_data <- read.csv("gene_copy_summary.csv")
gene_ranks <- setNames(gene_copy_data$total_copies, gene_copy_data$gene_name)

# Get pathway databases
go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
pathways <- split(go_bp$gene_symbol, go_bp$gs_name)

# Run GSEA
gsea_results <- fgsea(pathways = pathways,
                      stats = gene_ranks,
                      minSize = 15,
                      maxSize = 500,
                      eps = 0.0)

# Filter significant results
significant_results <- gsea_results[padj < 0.05]
```

### Background Construction Template
```r
# Create genome-wide background with copy counts
create_background <- function(phr_gene_copies, genome_size = 20000) {
  # Start with single-copy assumption for all genes
  all_genes <- unique(c(names(phr_gene_copies), 
                       paste0("BACKGROUND_", 1:(genome_size - length(phr_gene_copies)))))
  
  background_weights <- rep(1, length(all_genes))
  names(background_weights) <- all_genes
  
  # Update with actual PHR copy counts
  background_weights[names(phr_gene_copies)] <- phr_gene_copies
  
  return(background_weights)
}
```

This investigation provides a comprehensive foundation for implementing weighted gene enrichment analysis that properly accounts for the variable copy numbers observed in PHR regions.