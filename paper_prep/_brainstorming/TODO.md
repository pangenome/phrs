# TODO: PHR-specific gene enrichment analysis

## Context

Angela Gyamfi (Heather Mefford's student, Memphis) ran GSEA on genes within
**1 Mb** of PHR boundaries. Results show strong OR gene enrichment (146-fold,
z=18.0) and innate immune gene enrichment, with histone/keratin/MHC class II
depletion. But the 1 Mb window is ~10x wider than the median PHR (105 kb),
so the analysis captures the *neighborhood*, not just the PHRs themselves.

We want to redo this with genes **within the actual PHR intervals** only.

## What we have

### PHR coordinates on CHM13

The file `CHM13-HG002.sub-telo-phrs.bed` contains PHR boundaries for both
CHM13 and HG002. The CHM13 rows (lines 1-37) use PanSN naming:

```
CHM13#0#chr1    2706    297706    chr1,chr2,...
CHM13#0#chr1    248369277    248384276    chr1,chr2,...
...
```

Column 4 lists which chromosomes share sequence at that arm. There are 37
CHM13 PHR intervals covering the 41 arms with inter-chromosomal signal
(some chromosomes have both p and q arm entries).

**To make a standard BED**: strip the `CHM13#0#` prefix, giving coordinates
directly on CHM13v2.0 chromosome names.

### Gene annotations on CHM13

`chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz` — RefSeq gene models lifted to
CHM13v2.0. This is already in the repo.

Alternatively, Andrea's pipeline used CHM13 Liftoff GFF3 files at:
```
/moosefs/guarracino/HPRCv2/PHR_III/hprc_annotations/*.gff3.gz
```
The CHM13 reference annotation is among these 465 files. This is what was
used for the per-community gene enrichment in section 9 of the report.

For GO enrichment, we should use GENCODE on CHM13 if available, or the
RefSeq Liftoff already in the repo.

### Angela's existing gene catalog

`PHR_Subtelomeric Regions_Summary_March 2026.xlsx` has 511 gene entries
(327 unique) within CHM13 subtelomeric regions, mapped to the 15 Leiden
communities. Sheets:
- Overview: biotype counts (65 protein-coding, 101 lncRNA, 130 pseudogene, ...)
- Protein Coding Genes: 65 unique with disease annotations
- lncRNAs: 101 unique with function annotations
- Novel T2T: copies resolved by T2T but missing from GRCh38
- Gene Clustering by Chromosome: per-arm gene lists
- PCA+Gene Clustering: genes mapped to 15 Leiden communities (C1-C15)

### Angela's GSEA results (1 Mb window)

- `Figure1_GSEA_BP_vertical.pdf` — GO Biological Process
- `Figure_GSEA_MF_vertical.pdf` — GO Molecular Function
- `PHR_enrichment_summary.xlsx` — top-line results (ORA + permutation)
- `PHR_enrichment_all_results.xlsx` — all 89 significant GSEA terms

## Task: PHR-only GO enrichment

### Step 1: Extract CHM13 PHR BED

From `CHM13-HG002.sub-telo-phrs.bed`, extract CHM13 rows and strip PanSN:

```bash
grep '^CHM13#0#' CHM13-HG002.sub-telo-phrs.bed \
  | sed 's/CHM13#0#//' \
  > chm13.phrs.bed
```

This gives ~37 intervals (the actual PHR regions on CHM13, not 1 Mb windows).

### Step 2: Get CHM13 gene annotations

Option A — use RefSeq Liftoff already in repo:
```bash
zcat chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz \
  | awk '$3 == "gene"' \
  | bedtools intersect -a - -b chm13.phrs.bed -wa \
  > phrs.genes.gff3
```

Option B — use GENCODE CHM13 annotation (download if needed):
```
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/
GRCh38_mapping/gencode.v46.annotation.gff3.gz
```
But we need CHM13 coords, not GRCh38. The CAT/Liftoff annotations from
HPRC are the right choice. Get the CHM13 one from Andrea's cluster:
```
/moosefs/guarracino/HPRCv2/PHR_III/hprc_annotations/
```

### Step 3: Run GO enrichment

With the gene list from step 2, run ORA (over-representation analysis)
using clusterProfiler in R or gprofiler2:

```r
library(clusterProfiler)
library(org.Hs.eg.db)

# PHR genes (from bedtools intersect)
phr_genes <- read_gene_list("phrs.genes.txt")

# Background: all genes in the genome (or all genes in subtelomeric 500kb)
# Using all genes as background is conservative

enrichGO(gene = phr_genes,
         OrgDb = org.Hs.eg.db,
         ont = "BP",
         pAdjustMethod = "BH",
         pvalueCutoff = 0.05)
```

Key decision: **what is the background set?**
- All human genes → most conservative, tests "are PHR genes enriched vs genome"
- Genes in 500 kb subtelomeric flanks → tests "are PHR genes enriched vs the broader subtelomeric region" (this is what Angela's 1 Mb GSEA approximated)
- Genes in 1 Mb of chromosome ends → closest to Angela's original analysis

Recommendation: run with **all human genes** as background first (simplest,
most defensible), then optionally with subtelomeric background.

### Step 4: Compare to Angela's results

The PHR-only analysis should:
- **Sharpen** the OR gene signal (OR genes are within PHRs, not just nearby)
- **Lose** the negative enrichments for histones/keratins/MHC (those are in
  chromosome interiors, not near PHRs — they appeared in the 1 Mb GSEA
  because the window was wide enough to be informative about what's NOT
  near chromosome ends generally)
- **Confirm** Angela's community-specific gene mapping in the xlsx

### Step 5: Reconcile with Andrea's report (section 9)

Andrea's report already has per-community gene enrichment (section 9) using
174K gene annotations across 233 samples. Key findings:
- 374 unique genes across 39 arms
- Predominantly pseudogenes (28-86%) and ncRNA
- OR4F family in 7 communities
- DUX4L in C1, MTCO in C7, SHOX in C15
- RPL23AP45 spans 10 communities (most widespread hub gene)

The CHM13-only GO enrichment should be consistent with this but provide
the standard GO framework (BP, MF, KEGG) for the paper.

## Files to get from Andrea's cluster

```
# The main PHR TSV (per-haplotype, all 18,827 flanks)
/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv

# CHM13 gene annotation (Liftoff GFF3)
/moosefs/guarracino/HPRCv2/PHR_III/hprc_annotations/  (the CHM13 file)

# Community assignments
# (embedded in Andrea's R scripts / RDS files — need the arm-level
#  and sequence-level Leiden community labels)
```

## Notes

- The `chm13-annotations.bed` file in the repo has centromere, PAR, XTR,
  PHR-sex, and PHR-acro coordinates — the latter two are from the 2023
  acrocentric paper. The `CHM13-HG002.sub-telo-phrs.bed` is the genome-wide
  PHR file from the current analysis.
- Angela's xlsx maps genes to communities already — we may be able to use
  her gene list directly for GO enrichment without re-running the bedtools
  intersection, but should verify coordinates match.
- The RefSeq Liftoff in the repo (`chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz`)
  is from August 2024 — should be current enough.
