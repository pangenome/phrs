# Extended Data Figure 4 — Gene enrichment, copy-weighted GO, pseudogene gradient

**a, GO:BP top terms (PHR-only gene set, n=23 protein-coding genes).** Top 10 GO Biological Process terms ranked by adjusted p-value (gprofiler ORA). Spliceosome / snRNP terms (blue) dominate the list (top hit: "formation of quadruple SL/U4/U5/U6 snRNP", p_adj = 1.45e-3); intersection / term-size labels shown to the right of each bar. The original GSEA at `Figure1_GSEA_BP_vertical.pdf` used a 1 Mb window around PHR boundaries (~10× wider than the 105 kb median PHR), so neighborhood genes inflate the hit list — a PHR-only re-run is flagged as a Gap in `WORK_DECOMPOSITION.md` and is reproduced here. Source: `phr_GO_BP_enrichment.csv`.

**b, Copy-weighted vs deduplicated GO enrichment.** Copy-weighted fold enrichment for the 5 GO terms whose significance changes when each gene copy (rather than each unique symbol) is counted. "Sensory perception of smell" reaches 598× over expectation (copy_pvalue ≈ 0; dedup p = 0.040). Source: `improved_copy_weighted_vs_deduplicated_comparison.csv`.

**c, High-copy gene families (top 15, coding + pseudogene).** Copy count across PHRs for the 15 highest-copy protein-coding and pseudogene families. OR4F variants, IL9R/IL9RP1, DUX4, FRG2/FRG2B all carry 16–20 copies. Excludes lncRNA / miRNA (e.g. MIR8078 = 672). Source: `gene_copy_summary.csv`.

**d, OR4F pseudogenisation gradient by arm.** Pseudogene fraction of OR4F annotations per arm (n = 5,023 entries across 16 arms). Range 11.1% (chr7p) to 99.8% (chr15q); population mean 62.1% (red dashed line). Source: `or4f_pseudogene_fraction.csv`.

(Word count ≈ 199.)
