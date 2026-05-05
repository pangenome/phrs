# SURVEY: Tabular Data Inventory (Nature MS prep)

**Scope.** All `*.csv / *.tsv / *.RData / *.rds / *.parquet / *.bed / *.gff[3] / *.txt` files under repo root that look like tables. Subdirectories `scripts/`, `paper_prep/`, `end-to-end-report/` contain no tabular data files — every table sits at the project root.

**Total cataloged: 46 / 46 files (100%).** All inspected directly except for binary R serializations (`*.RData`, `*.rds`), where contents are inferred from filename + companion CSV exports.

**Use-tier legend.**
- `MAIN` — directly drives a main-text figure / table
- `EXT` — strong candidate for an Extended Data table or panel
- `SI` — Supplementary Information / supporting table
- `SUP` — supporting input list (gene IDs etc.) — methods footnote at most
- `SKIP` — text dump / debug / superseded; not for the paper

---

## 1. PHR coordinate tracks and gene annotations (core figure inputs)

| PATH | SHAPE | TOPIC | USE | NOTES |
|---|---|---|---|---|
| `chm13.phrs.bed` | 37 × 4 | PHR call set on CHM13 (all chromosomes) | **MAIN** | Canonical PHR intervals; cols = chrom, start, end, sharing-pattern (comma-sep chrom list). Drives every PHR-track panel. |
| `chm13.phrs.no_acro.bed` | 29 × 4 | PHR call set excluding acrocentric short arms | **MAIN** | Acro-excluded variant used for the conservative gene-enrichment analyses. |
| `CHM13-HG002.sub-telo-phrs.bed` | 112 × 4 | Per-haplotype PHR intervals (CHM13 + HG002 hap0/hap1) | **EXT** | 3× the row count vs the CHM13-only call set → contains HG002 phased PHRs. Good for an Extended Data "PHRs across 3 haplotypes" panel. |
| `chm13-annotations.bed` | 95 × 4 | CHM13 reference annotations (centromere, gap features) | **SUP** | Background reference track used for ideogram overlays. |
| `phrs.genes.gff3` | 412 × 9 (Liftoff GFF3) | All genes lifted into PHR intervals (incl. acrocentrics) | **EXT** | Source for "gene census in PHRs" Extended Data table; rich attribute column (biotype, copy_num_ID, sequence_ID). |
| `phrs.no_acro.genes.gff3` | 361 × 9 | Same, acrocentrics excluded | **EXT** | Companion to the no_acro analyses. |
| `phrs.no_acro.coding_genes.gff3` | 23 × 9 | Protein-coding-only, no acrocentrics | **MAIN** | The 23-gene shortlist that anchors the GO/KEGG figure. |

## 2. Per-copy gene catalogs (origin of every "gene copies in PHRs" figure)

| PATH | SHAPE | TOPIC | USE | NOTES |
|---|---|---|---|---|
| `all_gene_copies_by_arm.csv` | 1,189 × 13 | One row per gene copy in PHRs, with arm/biotype/sharing-pattern | **MAIN** | The "rows-per-copy" master table; supports per-arm bar charts and per-biotype breakdowns. |
| `genome_wide_gene_copies.csv` | 58,230 × 5 | Genome-wide copy-number background per gene | **EXT** | Background distribution for copy-weighted enrichment; also a clean SI table. |
| `comprehensive_copy_background.csv` | 58,230 × 8 | Same per-gene rows annotated with `in_phrs`/`in_genome` flags + chromosome list | **SI** | Superset of `genome_wide_gene_copies.csv`; pick one for SI, not both. |
| `gene_copy_summary.csv` | 35 × 6 | Top high-copy gene families with arm + community lists | **MAIN** | Compact "DUX4, BAGE2, …" headline table. Strong Extended Data candidate as-is. |
| `enriched_genes_detailed_map.csv` | 72 × 14 | Per-gene functional cluster + community + brief function | **MAIN** | Drives the annotated network/cluster figure (snRNP, olfactory, etc.). |

## 3. Functional enrichment — standard ORA (gprofiler outputs)

| PATH | SHAPE | TOPIC | USE | NOTES |
|---|---|---|---|---|
| `phr_GO_BP_enrichment.csv` | 35 × 9 | gprofiler GO:BP enrichment, all PHR genes | **MAIN** | Lead enrichment table. Top hit = "formation of quadruple SL/U4/U5/U6 snRNP" (p_adj=1.45e-3). |
| `phr_GO_MF_enrichment.csv` | 6 × 9 | GO:MF enrichment, all PHR genes | **EXT** | Top hit = "U4 snRNA binding" (9.1e-5). |
| `phr_KEGG_enrichment.csv` | 0 rows + comment line | KEGG enrichment | **SKIP** | Header-only file; comment "No significant KEGG pathways found". Cite as a negative result in methods, do not publish as a table. |
| `phr_no_acro_GO_BP_enrichment.csv` | 35 × 9 | GO:BP, acrocentrics excluded | **EXT** | Confirms GO:BP signal is robust to acro removal. |
| `phr_no_acro_GO_MF_enrichment.csv` | 6 × 9 | GO:MF, acrocentrics excluded | **EXT** | Companion to above. |
| `phr_coding_only_GO_BP_enrichment.csv` | 7 × 9 | GO:BP restricted to 23 coding genes | **EXT** | Olfactory + sensory chemical stimulus signal. |
| `phr_coding_only_GO_MF_enrichment.csv` | 9 × 9 | GO:MF restricted to coding genes | **EXT** | Olfactory receptor activity tops the list. |

## 4. Copy-weighted / corrected enrichment (custom analysis)

| PATH | SHAPE | TOPIC | USE | NOTES |
|---|---|---|---|---|
| `phr_copy_weighted_enrichment.csv` | 5 × 11 | First-pass copy-weighted GO enrichment (small) | **SI** | Largely superseded by the `improved_*` variants below. |
| `improved_copy_weighted_enrichment.csv` | 6 × 13 | Final copy-weighted GO enrichment with fold/p_adj | **MAIN** | Carries the "olfactory signal is real after copy correction" claim (fold=598 over expectation). |
| `copy_weighted_permutation_results.csv` | 6 × 4 | Empirical null permutation p-values | **EXT** | Robustness control for copy-weighted enrichment. |
| `copy_weighted_functional_analysis.csv` | 3 × 4 | High-level category counts (e.g., olfactory_receptor) | **SI** | Tiny rollup; mainly for inline narrative. |
| `copy_weighted_vs_deduplicated_comparison.csv` | 2 × 8 | Per-term comparison: dedup p-value vs copy p-value | **EXT** | Use as a small Extended Data table that shows method matters. |
| `improved_copy_weighted_vs_deduplicated_comparison.csv` | 5 × 11 | Same idea, more terms + fold-change ratios | **EXT** | Pick this over the smaller `copy_weighted_vs_deduplicated_comparison.csv` for the figure. |
| `copy_number_vs_standard_ora_comparison.csv` | 4 × 7 | Method-vs-method qualitative summary | **EXT** | Already in the right shape for an Extended Data "methods comparison" panel. |
| `comparison_table.csv` | 18 × 9 | PHR-only vs Angela-1Mb GO/category comparison with rank & interpretation | **MAIN** | Reads almost like a finished Extended Data table; minimal cleanup needed. |

## 5. Method validation & benchmarks (R serializations + summary CSVs)

| PATH | SHAPE | TOPIC | USE | NOTES |
|---|---|---|---|---|
| `phyper_benchmark_summary.csv` | 3 × 14 | phyper() benchmark: timing + p-value match across small/medium/large | **EXT** | Compact benchmark table for an Extended Data "performance" panel. |
| `performance_summary.csv` | 4 × 8 | Same flavour: speedup / memory ratio across datasets | **EXT** | Companion to the phyper benchmark; merge into one table. |
| `phyper_benchmark_detailed_results.rds` | binary R serialisation, ~1.6 KB | Detailed benchmark records | **SI** | Underlying data for `phyper_benchmark_summary.csv`; ship as supplementary file, not as a printed table. |
| `benchmark_results.RData` | binary R serialisation, ~1.7 KB | Generic benchmark results object | **SI** | Likely the workspace dump matching `performance_summary.csv`; SI-only. |
| `comprehensive_validation_results.RData` | binary R serialisation, ~700 KB | Full validation suite (the largest binary table here) | **SI** | Source for `null_distribution_validation_report.txt`; ship as supplementary. |
| `type_i_error_validation_results.RData` | binary R serialisation, ~25 KB | Type-I error sweep | **EXT** | Worth re-exporting one summary CSV (alpha vs empirical FPR) for an Extended Data calibration panel. |
| `error_handling_validation_results.rds` | binary R serialisation, ~1.5 KB | Edge-case / error-path validation | **SI** | Methods/SI only. |
| `integration_test_results.rds` | binary R serialisation, ~1 KB | Pipeline integration smoke results | **SKIP** | CI smoke output; not a paper artifact. |

## 6. Gene / Entrez ID lists (gprofiler inputs and gene-name dumps)

| PATH | SHAPE | TOPIC | USE | NOTES |
|---|---|---|---|---|
| `phrs.gene_names.txt` | 412 × 1 | All gene symbols inside PHRs (per-copy, with duplicates) | **SUP** | Direct gprofiler input. Methods footnote. |
| `phrs.unique_gene_names.txt` | 245 × 1 | Deduplicated gene symbols | **SUP** | Standard-ORA input list. |
| `phrs.entrez_ids.txt` | 113 × 1 | Entrez IDs (per-copy) | **SUP** | Methods footnote. |
| `phrs.unique_entrez_ids.txt` | 69 × 1 | Deduplicated Entrez IDs | **SUP** | Methods footnote. |
| `phrs.no_acro.gene_names.txt` | 220 × 1 | Acrocentric-excluded gene symbols | **SUP** | Methods footnote. |
| `phrs.no_acro.coding_gene_names.txt` | 23 × 1 | The 23-gene coding shortlist | **SI** | Worth printing in full as a small inline SI list (matches `phrs.no_acro.coding_genes.gff3`). |
| `gene_list_for_gprofiler_no_acro.txt` | 220 × 1 | Same content as `phrs.no_acro.gene_names.txt` | **SKIP** | Duplicate of the no_acro gene-names list — reference one, drop the other. |

## 7. Text dumps / human-readable reports (not tables)

| PATH | SHAPE | TOPIC | USE | NOTES |
|---|---|---|---|---|
| `phr_GO_BP_dotplot.txt` | 122 lines | Pretty-printed GO:BP top-20 (text dotplot) | **SKIP** | Same data as `phr_GO_BP_enrichment.csv`; cite the CSV. |
| `phr_GO_MF_dotplot.txt` | 38 lines | Pretty-printed GO:MF | **SKIP** | Same as above for MF. |
| `phr_enrichment_summary.txt` | 28 lines | Plain-text rollup of the enrichment runs | **SKIP** | Narrative summary; redundant with the CSVs. |
| `null_distribution_validation_report.txt` | 173 lines | Validation report (header "COPY-NUMBER WEIGHTED PHYPER() NULL DISTRIBUTION VALIDATION REPORT", 2026-04-01) | **SI** | Useful as a methods appendix; could lift one summary table out of it. |

---

## Opportunities — strong Extended Data / SI candidates

These are tables whose contents would slot directly into the manuscript with minimal reshaping:

1. **PHR call set across CHM13 + HG002 haplotypes** — `CHM13-HG002.sub-telo-phrs.bed` (112 rows) is currently underused. An Extended Data table listing every PHR per haplotype with sharing pattern would (a) make the call set fully reproducible and (b) let reviewers see CHM13 ↔ HG002 concordance at a glance.
2. **Method-comparison "what changed" panel** — `comparison_table.csv` (PHR-only vs Angela-1Mb, with rank + interpretation column) and `improved_copy_weighted_vs_deduplicated_comparison.csv` together form a natural Extended Data figure: one column for the dedup p-value, one for the copy-weighted p-value, fold-change ratio, and a written interpretation. Almost no work needed.
3. **The 23-gene coding shortlist as a printed table** — `phrs.no_acro.coding_genes.gff3` + `phrs.no_acro.coding_gene_names.txt`: print all 23 with chromosome, copy count (from `gene_copy_summary.csv`), biotype, and brief function (from `enriched_genes_detailed_map.csv`). This is the headline gene set; printing it in full kills any "which 23 genes?" reviewer question.
4. **High-copy gene families** — `gene_copy_summary.csv` (35 rows incl. DUX4 ×18, BAGE2, …) is already in the exact shape of an Extended Data table. Add a "biotype" colour column in the figure and ship.
5. **Statistical robustness summary** — combine `phyper_benchmark_summary.csv` + `performance_summary.csv` + a one-row summary lifted from `null_distribution_validation_report.txt` and `type_i_error_validation_results.RData` (re-export needed) into a single Extended Data "method validation" table covering speed, memory, type-I error, and null calibration.
6. **Gene census per chromosome arm** — `all_gene_copies_by_arm.csv` (1,189 rows) is too large to print, but a derived per-arm summary (counts by arm × biotype) would make a clean Extended Data heatmap and is one `groupby` away.
7. **Olfactory receptor copy table** — `enriched_genes_detailed_map.csv` filtered to the olfactory cluster is small and self-explanatory; pulling it out as its own SI table would reinforce the OR-cluster claim that's central to the GO results.

## Notes for the synthesis architect

- **CHM13-only vs no_acro vs coding_only is the analysis tree** — make sure each figure caption states which gene set it uses; the three tracks above are *not* interchangeable and the rows-per-table differ accordingly.
- **`comprehensive_copy_background.csv` is a strict superset of `genome_wide_gene_copies.csv`** (8 cols vs 5, same row count). Drop one before the SI freeze.
- **`gene_list_for_gprofiler_no_acro.txt` ≡ `phrs.no_acro.gene_names.txt`** (both 220 rows, identical first lines). Same comment.
- **The KEGG file is empty by design** — `phr_KEGG_enrichment.csv` only has a header + a "no significant pathways" comment line. Cite as a negative result; do not table.
- **R serializations need re-export.** None of the `*.RData` / `*.rds` files were inspected directly; if a downstream task wants to print, e.g., the type-I error sweep, an `R --no-save` round trip is required to write a CSV.
