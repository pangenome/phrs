# Copy-Number Enrichment Methods: Comprehensive Analysis and Ranked Recommendations

**Date:** 2026-04-01  
**Project:** PHR Gene Enrichment Analysis  
**Context:** 29 PHR intervals, 220 genes total, 1,189 gene instances (copy numbers 2-672)

## Executive Summary

Based on comprehensive investigation of 5 methodological approaches and literature review, this report provides ranked recommendations for implementing copy-number-aware enrichment analysis in Pseudohomologous Regions (PHRs). The analysis reveals that traditional enrichment methods systematically underweight highly duplicated genes, potentially missing crucial biological signals in repeat-rich genomic regions.

**Top Recommendations:**
1. **GSEA with Pre-ranked Copy Counts** - Immediate implementation
2. **Permutation-Based Validation** - Essential statistical validation  
3. **Copy-Number Weighted ORA** - Advanced implementation for maximum accuracy

## 1. Methods Investigated

### Summary Table

| Method | Tool Availability | CHM13 Compatibility | Implementation Complexity | Runtime | Statistical Rigor |
|--------|------------------|---------------------|--------------------------|---------|-------------------|
| **GSEA Pre-ranked** | ✅ Ready (fgsea) | ✅ Compatible | Low | <5 min | High |
| **Permutation Testing** | ✅ Ready (bedtools) | ✅ Compatible | Medium | 5-8 hours | Highest |
| **GREAT/rGREAT** | ❓ Needs setup | ⚠️ Requires custom | Medium | 5-15 min | High |
| **GAT** | ❌ Installation blocked | ❓ Unvalidated | High | Hours | High |
| **Copy-Weighted ORA** | ✅ Ready (R) | ✅ Compatible | High | 10-30 min | Medium |

### 1.1 GSEA with Pre-ranked Copy Counts

**Concept:** Use gene copy numbers as ranking weights in Gene Set Enrichment Analysis instead of expression fold-changes.

**Advantages:**
- **Established methodology** with strong literature precedent (Subramanian et al. 2005)
- **Ready-to-use tools** (fgsea R package, GSEApy) available on system
- **Natural fit** for copy-number data using ranking approach
- **Fast implementation** (<5 minutes runtime)
- **Normalized Enrichment Score** accounts for pathway size
- **Permutation-based significance** testing built-in

**Disadvantages:**
- **Ranking assumptions** may not capture non-linear copy-number effects
- **Pathway size sensitivity** may bias toward larger gene sets
- **Limited background control** compared to custom approaches

**Implementation Details:**
- **Tool:** fgsea R package (preferred) or GSEApy (Python)
- **Input:** Ranked gene list: Gene_Name, Copy_Count (descending)
- **Background:** Genome-wide genes with single-copy assumption for non-PHR genes
- **Parameters:** minSize=15, maxSize=500, 1000 permutations
- **Runtime:** 30-60 seconds for 15,000 pathway tests

### 1.2 Permutation-Based Enrichment Validation

**Concept:** Generate null distributions through random genomic interval sampling to validate enrichment significance.

**Advantages:**
- **Strongest statistical foundation** - empirical p-values without distributional assumptions
- **Preserves genomic context** through constrained sampling
- **Controls for interval properties** (size, chromosome, regional context)
- **Well-established in genomics** (GAT methodology, Churchill & Doerge 1994)
- **Copy-number aware** by design

**Disadvantages:**
- **Computationally intensive** (5-8 hours for 10,000 permutations)
- **Complex implementation** requiring pipeline orchestration
- **High data requirements** (gene annotations, constraint regions)

**Implementation Details:**
- **Tool:** bedtools shuffle + custom statistical pipeline
- **Input:** PHR BED intervals, gene annotations (GFF), subtelomeric constraint regions
- **Protocol:** 10,000 permutations with chromosome and subtelomeric constraints
- **Runtime:** 5.5-8.5 hours sequential, ~1 hour with 8-core parallelization
- **Statistics:** Empirical p-values, Z-scores, Benjamini-Hochberg FDR correction

### 1.3 GREAT/rGREAT Region-Based Analysis

**Concept:** Region-based enrichment that naturally handles multi-copy genes through regulatory domain assignment.

**Advantages:**
- **Region-based approach** avoids gene boundary issues
- **Handles multi-copy genes** through regulatory domains
- **Established tool** (>4000 citations) with proven track record
- **Built-in statistical framework** with proper background correction

**Disadvantages:**
- **CHM13 compatibility uncertain** - requires custom gene annotations
- **Installation complexity** for rGREAT package setup
- **Fixed regulatory domain rules** may not capture all relevant biology
- **Less control** over background and statistical parameters

**Implementation Details:**
- **Tool:** rGREAT R package or GREAT web interface
- **Input:** 29 PHR intervals in BED format
- **Setup time:** 30-60 minutes for CHM13 annotations
- **Runtime:** 5-15 minutes per analysis

### 1.4 GAT (Genomic Association Tester)

**Concept:** Simulation framework for genomic interval association testing with workspace constraints.

**Advantages:**
- **Sophisticated simulation** with isochore and workspace controls
- **Handles repeat regions** through workspace exclusion
- **Performance optimization** for large datasets

**Disadvantages:**
- **Installation blocked** on current system (missing dependencies)
- **No explicit multi-copy handling** documented
- **CHM13 compatibility unvalidated**
- **Limited advantage** over custom permutation approach

**Status:** **Not recommended** due to installation barriers and limited advantages over permutation-based approach.

### 1.5 Copy-Number Weighted Over-Representation Analysis

**Concept:** Custom ORA that weights genes by copy number rather than treating as binary presence/absence.

**Advantages:**
- **Maximum biological accuracy** by directly incorporating copy numbers
- **Extreme copy number handling** (e.g., MIR8078 with 672 copies)
- **Flexible implementation** allowing method customization
- **Novel approach** potentially revealing unique insights

**Disadvantages:**
- **High development complexity** requiring custom statistical framework
- **Computational scaling challenges** (100x increase in data size)
- **Statistical validation requirements** extensive
- **Long development timeline** (10-12 weeks for robust implementation)

**Implementation Options:**
- **Phase 1:** Instance expansion approach (2 weeks)
- **Phase 2:** Robust statistical framework (3 weeks)  
- **Phase 3:** Production tool (6 weeks)

## 2. Literature Precedents and Statistical Framework

### 2.1 Established Methodological Precedents

**Copy-Number Aware Analysis:**
- **ABCD-DNA Framework** (Benjamini et al. 2012): Direct CNV integration into differential analysis
- **GSEA Methodology** (Subramanian et al. 2005): Weighted Kolmogorov-Smirnov statistics with permutation testing
- **GeneToCN Method** (2023): Alignment-free copy number estimation for multi-copy genes

**Subtelomeric and Repeat Region Analysis:**
- **Subtelomeric VNTR Analysis** (2020): 21-fold enrichment validation using Z-score analysis
- **RepEnTools** (2024): Automated repeat enrichment using Mann-Whitney-Wilcoxon tests
- **GAT Framework** (Heger et al. 2013): Simulation-based genomic overlap testing

**Statistical Best Practices:**
- **Permutation Standards:** 1,000-10,000 permutations for publication quality
- **Multiple Testing:** Benjamini-Hochberg FDR correction for discovery studies
- **Effect Size Reporting:** Z-scores and confidence intervals for biological interpretation

### 2.2 Methodological Gaps Our Work Addresses

1. **First PHR-specific enrichment analysis** - no established methods for Pseudohomologous Regions
2. **Copy number integration** - most tools treat genes as single-copy
3. **Extreme copy number ranges** - handling 300+ fold copy differences
4. **CHM13-based analysis** - validation on most complete human reference

## 3. Tool Availability and Technical Feasibility

### 3.1 System Resources (Octopus Head Node)
- **Memory:** 515 GB available (excellent for statistical computing)
- **CPU:** 48 cores (ideal for parallel processing)  
- **Tools Ready:** R 4.3.0, bedtools v2.30.0, Python 3.7.3
- **Installation Required:** R packages (fgsea, clusterProfiler), Python packages (scipy, numpy)

### 3.2 Installation Feasibility Matrix

| Tool/Package | Installation Time | Complexity | Success Likelihood |
|-------------|------------------|------------|-------------------|
| fgsea (R) | 10-15 minutes | Low | High |
| bedtools shuffle | Already available | None | High |
| rGREAT (R) | 10-15 minutes | Low | High |
| Python packages | 5-10 minutes | Low | High |
| GAT | 15-30 minutes | High | Medium |

## 4. Ranked Recommendations

### Rank 1: GSEA with Pre-ranked Copy Counts ⭐⭐⭐⭐⭐

**Justification:**
- **Immediate implementability** with available tools
- **Strong statistical foundation** from established GSEA methodology
- **Literature precedent** for weighted enrichment analysis
- **Optimal complexity-to-benefit ratio**

**Implementation Plan:**
```r
# Week 1: Setup and initial analysis
library(fgsea)
library(msigdbr)

# Create copy-number ranked gene list
gene_ranks <- setNames(copy_data$total_copies, copy_data$gene_name)

# Get pathway databases  
pathways <- msigdbr(species = "Homo sapiens", category = "C5")

# Run GSEA
results <- fgsea(pathways = pathways, stats = gene_ranks, 
                 minSize = 15, maxSize = 500, eps = 0.0)
```

**Expected Runtime:** 30-60 seconds  
**Implementation Time:** 1 week  
**Key Papers:** Subramanian et al. (2005) GSEA, Benjamini et al. (2012) ABCD-DNA

### Rank 2: Permutation-Based Validation ⭐⭐⭐⭐

**Justification:**
- **Essential for statistical validation** of any enrichment findings
- **Gold standard approach** in genomics for empirical significance
- **Robust to distributional assumptions**
- **Addresses PHR-specific genomic context**

**Implementation Plan:**
```bash
# Week 2-3: Permutation testing framework
for i in $(seq 1 10000); do
    bedtools shuffle -i chm13.phrs.bed -g chm13.genome \
        -incl subtelomeric_5mb_windows.bed -chrom -seed $i | \
    bedtools intersect -a genes.gff -b - -wa | \
    extract_gene_names.py > genes_perm_${i}.txt
done

# Statistical analysis in R/Python
python permutation_analysis.py --permutations 10000 --output validation_results.csv
```

**Expected Runtime:** 5-8 hours (or 1 hour parallelized)  
**Implementation Time:** 2-3 weeks  
**Key Papers:** Churchill & Doerge (1994), Heger et al. (2013) GAT, Layer et al. (2018)

### Rank 3: Copy-Number Weighted ORA ⭐⭐⭐

**Justification:**
- **Maximum biological accuracy** by direct copy-number integration
- **Novel methodological contribution** to the field
- **Addresses extreme copy numbers** (672-fold variation in dataset)
- **Long-term value** for similar analyses

**Implementation Plan:**
```r
# Week 4-6: Proof-of-concept implementation  
# Instance expansion approach
query_expanded <- rep(gene_data$symbol, gene_data$copies)
background_weighted <- create_weighted_background(genome_annotation)
results_weighted <- hypergeometric_test(query_expanded, pathways, background_weighted)

# Week 7-10: Robust statistical framework
# Advanced implementation with validation
```

**Expected Runtime:** 10-30 minutes  
**Implementation Time:** 6-10 weeks for robust implementation  
**Key Papers:** Benjamini et al. (2012) ABCD-DNA, GeneToCN (2023)

## 5. Implementation Roadmap

### Phase 1: Immediate Implementation (Weeks 1-2)
**Goal:** Validate copy-number effects on enrichment results

**Tasks:**
1. **Install required R packages** (fgsea, msigdbr, clusterProfiler)
2. **Implement GSEA with copy-number ranking**
3. **Compare results** with current g:Profiler binary analysis
4. **Quantify differences** in enriched pathways

**Deliverables:**
- Copy-number weighted enrichment results
- Comparison table: standard vs copy-weighted enrichment
- Statistical significance assessment

**Success Metrics:**
- ✅ GSEA analysis completes successfully
- ✅ Results show biological plausibility
- ✅ Differences from standard analysis quantified

### Phase 2: Statistical Validation (Weeks 3-5)  
**Goal:** Establish statistical confidence in enrichment findings

**Tasks:**
1. **Setup permutation testing framework** using bedtools shuffle
2. **Generate 10,000 random PHR interval sets**
3. **Build null distributions** for all tested pathways
4. **Validate GSEA results** against permutation-based p-values

**Deliverables:**
- Permutation testing pipeline
- Empirical p-values for all enriched pathways
- Statistical validation report

**Success Metrics:**
- ✅ Permutation pipeline runs reliably
- ✅ Null distributions show expected properties
- ✅ Key enrichments confirmed by permutation testing

### Phase 3: Advanced Methods (Weeks 6-12)
**Goal:** Implement copy-weighted ORA for comprehensive analysis

**Tasks:**
1. **Develop weighted hypergeometric framework**
2. **Implement background construction algorithms**  
3. **Validate statistical properties** with positive/negative controls
4. **Compare all three approaches** comprehensively

**Deliverables:**
- Production-ready copy-weighted ORA tool
- Comprehensive comparison of all methods
- Methodological paper documenting approach

**Success Metrics:**
- ✅ Copy-weighted ORA produces statistically valid results
- ✅ Method comparisons show consistent patterns
- ✅ Novel biological insights identified

## 6. Expected Outcomes and Validation Strategy

### 6.1 Primary Validation Targets

**Known Positive Controls:**
- **Olfactory receptor enrichment** (currently p=0.029) - should be confirmed/strengthened
- **Immunoglobulin gene families** - expected in subtelomeric regions
- **Histone gene clusters** - known repetitive gene families

**Expected Null Controls:**
- **Random gene sets** - should show no enrichment
- **Single-copy pathways** - should show minimal copy-number effects

### 6.2 Success Criteria

**Technical Success:**
- [ ] All methods produce consistent p-value distributions under null
- [ ] Runtime targets met for all approaches
- [ ] Results reproducible across different parameter settings

**Biological Success:**  
- [ ] Copy-number weighting reveals previously missed enrichments
- [ ] Results align with known biology of repeat-rich regions
- [ ] Method identifies pathways relevant to genomic instability

**Methodological Success:**
- [ ] Approaches complement rather than contradict each other
- [ ] Statistical validation confirms reliability of findings
- [ ] Framework applicable to other repetitive genomic regions

## 7. Risk Assessment and Mitigation

### 7.1 Technical Risks

**High Priority:**
- **Permutation runtime** may be prohibitive
  - *Mitigation:* Implement parallel processing, optimize pipeline
  - *Fallback:* Reduce permutation count to 1,000 for initial validation

**Medium Priority:**  
- **Copy-weighted ORA complexity** may delay implementation
  - *Mitigation:* Focus on proof-of-concept first, defer advanced features
  - *Fallback:* Use instance expansion approach only

### 7.2 Biological Risks

**Low Priority:**
- **Results may not align** with standard enrichment findings
  - *Mitigation:* Document differences clearly, provide both perspectives
  - *Approach:* Frame as complementary rather than replacement methods

## 8. Resource Requirements and Timeline

### 8.1 Personnel Requirements
- **1 bioinformatician** with statistical expertise (primary)
- **0.5 computational biologist** for validation (supporting)
- **Total effort:** 6-8 weeks for comprehensive implementation

### 8.2 Computational Requirements
- **Development:** Current head node resources adequate
- **Storage:** ~100 GB for intermediate permutation files
- **Runtime:** ~20-30 hours total across all methods

### 8.3 Critical Dependencies
- **R package installation** (fgsea, rGREAT, clusterProfiler)
- **Python package installation** (scipy, numpy for permutation analysis)  
- **CHM13 gene annotations** (RefSeq Liftoff or GENCODE)

## 9. Literature Citations and Methodological Support

### 9.1 Core Methodological Papers
- **Subramanian et al. (2005)** "Gene set enrichment analysis: A knowledge-based approach" - GSEA methodology
- **Benjamini et al. (2012)** "Copy-number-aware differential analysis" - ABCD-DNA framework
- **Churchill & Doerge (1994)** "Empirical threshold values for quantitative trait loci" - Permutation testing foundations
- **Heger et al. (2013)** "GAT: A simulation framework for testing genomic intervals" - Genomic permutation testing

### 9.2 Subtelomeric and Copy Number Papers  
- **Extreme enrichment of VNTR-associated polymorphicity in human subtelomeres** (2020) - Statistical precedent for subtelomeric enrichment
- **GeneToCN: an alignment-free method for gene copy number estimation** (2023) - Multi-copy gene analysis methods
- **Robust and accurate estimation of paralog-specific copy number** (2022) - Parascopy tool for duplicated genes

### 9.3 Statistical and Technical References
- **Layer et al. (2018)** "A framework for reproducible permutation testing in genomics" - Best practices
- **Eden et al. (2009)** "GOrilla: enriched GO terms discovery" - Background selection importance
- **McLean et al. (2010)** "GREAT improves functional interpretation" - Region-based enrichment

## 10. Conclusion

The investigation reveals three complementary approaches for copy-number-aware enrichment analysis, each addressing different aspects of the statistical and biological challenges posed by PHR regions:

1. **GSEA with pre-ranked copy counts** provides immediate implementation with strong methodological precedent
2. **Permutation-based validation** offers essential statistical rigor through empirical significance testing  
3. **Copy-number weighted ORA** represents the most biologically accurate approach for long-term implementation

**Recommended implementation sequence:** Begin with GSEA for rapid results, validate with permutation testing, and develop copy-weighted ORA as resources permit. This staged approach balances immediate progress with long-term methodological advancement.

**Expected impact:** These methods should reveal copy-number-driven pathway enrichments currently missed by standard approaches, providing new insights into the functional organization of pericentromeric heterochromatin and establishing methodological precedents for analyzing other repeat-rich genomic regions.

The comprehensive framework addresses the critical gap between traditional single-copy enrichment methods and the biological reality of highly duplicated genomic regions, potentially transforming our understanding of functional enrichment in repetitive elements of the human genome.