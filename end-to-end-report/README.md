# Subtelomeric Analysis Report

## Overview

**What it does.** Identifies and characterizes inter-chromosomal subtelomeric sequence sharing across 233 HPRCv2 human samples, validates the findings with 3D genome data from multiple technologies and organisms.

**How.** Extract the terminal 500 kb from each chromosome arm across 465 haplotypes (contigs ≥ 1 Mb). Align all-vs-all at ≥95% identity to find which arms share sequence across different chromosomes. Build a pangenome graph from these shared regions, then cluster arms into communities by graph similarity. Annotate communities with genes and repeats. Validate by testing whether same-community arms sit closer in 3D nuclear space — using Hi-C (6 human samples), Pore-C and CiFi (HG002), Dip-C single-cell (16 GM12878 cells), sperm single-cell (20 cells), RPE-1 (3 datasets), and mouse meiotic Hi-C (4 stages).

**Key metrics.** 233 samples (465 haplotypes), 48 chromosome arms (41 with inter-chromosomal signal), 15,668 PHR sequences, 15 arm-level communities, 50 sequence-level communities. 3D validation across 4 technologies, 2 cell types (lymphoblastoid cell lines [LCL], retinal pigment epithelium [RPE-1]), sperm, and mouse meiosis.

**Result.** Chromosome arms that share subtelomeric sequence cluster into discrete communities. These communities are reflected in 3D nuclear proximity — arms in the same community are physically closer than arms in different communities. This holds in bulk Hi-C, single-cell Dip-C, haploid sperm, and mouse meiotic cells.

**Conclusion.** Subtelomeric regions form a structured system of inter-chromosomal sharing, shaped by recurrent ectopic exchange and reflected in 3D nuclear organization across cell types and species.

---

## Table of Contents

### Part I: Subtelomeric Sequence Analysis

- [Pipeline](report/01_pipeline.md) — Contig classification, flank extraction, all-vs-all alignment, inter-chromosomal region detection, pangenome graph, community detection (arm-level 15 communities, sequence-level 50 communities)
- [Annotation](report/02_annotation.md) — TAR1 prevalence and positional distribution, internal (TTAGGG)n islands (boundary enrichment, length, motif composition, cross-arm status), terminal telomere tract length
- [Gene enrichment](report/03_gene_enrichment.md) — Biotype composition, olfactory receptors, Ambrosini duplicon blocks, subterminal families, community-specific genes, hub genes
- [Within-arm heterogeneity](report/04_heterogeneity.md) — Allele vs paralog distance, silhouette analysis, cross-arm affinity, population structure, type discordance, region length, gene replacement scoring, two-domain model test

### Part II: 3D Genome Validation

- [Hi-C/Pore-C validation](report/05_hic_validation.md) — Community-based B/W ratio, flanking analysis, Mantel test, ARI, per-arm-pair correlation, community-free sequence-level, multi-resolution (5kb–100kb), no-acrocentric control, O/E normalization, RPE-1 cell-type validation, per-community enrichment
- [Dip-C + sperm](report/06_dipc_validation.md) — GM12878 single-cell 3D (16 cells), community-free arm-level, supplementary communities, sperm single-cell (20 cells)
- [Integrated interpretation](report/07_integrated.md) — Convergent evidence table, flanking paradox, meiotic bouquet, D4Z4-CTCF-lamin model, nucleolar association, causal feedback loop, testable predictions

### Part III: Cross-Species and Cross-System Validation

- [Mouse T2T](report/08_mouse.md) — Mouse subtelomeric pipeline (B6 + CAST), community detection, Hi-C validation across 4 meiotic stages (1Mb/2Mb/4Mb windows), window-size optimization (1Mb–33Mb), flanking Hi-C, community-free sequence-level
- [RPE-1 self-vs-self](report/09_rpe1_self.md) — RPE-1 pipeline (self-discovered communities), Hi-C validation, comparison with HPRC communities, flanking control

### Part IV: Pedigree Analysis

- [Pedigree recombination](report/14_pedigree_recombination.md) — Subtelomeric recombination in two pedigrees (WashU 3-generation, CEPH1463 4-generation). Inter-chromosomal exchange in PHR regions, NAHR/gene conversion events, cross-assembler validation (hifiasm vs verkko), transmission across generations

### Part V: Summary and Context

- [Limitations](report/10_limitations.md) — Sample composition, methodological limitations, 3D validation limitations
- [Summary of key findings](report/11_summary.md) — 10 key biological findings with metrics
- [Literature and novelty](report/12_literature.md) — Confirmed literature claims, 27 novel contributions, testable predictions from existing data
- [Appendix + References](report/13_appendix.md) — External tools, references
