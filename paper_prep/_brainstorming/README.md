# Off-target materials from prior manuscript session

The files in this directory were produced during a prior synthesis run that hill-climbed on copy-number-weighted ORA / hypergeometric-test / OR4F-gene-family minutiae. They do **not** anchor on the canonical paper described in `paper_prep/synthesis/ABSTRACT.md` ("Concerted evolution and unorthodox recombination of human subtelomeres" — implicit pangenome graph methodology, 465-assembly HPRC v2 dataset, extended interchromosomal subtelomeric homology, NJ-tree cladistic analysis, PCA + community detection, Hi-C 3D-proximity hypothesis).

Retained here for possible salvage. **Not part of the canonical manuscript pipeline.**

For canonical material, see `paper_prep/synthesis/ABSTRACT.md` (the anchoring abstract). The figures under `paper_prep/figures/` have NOT been moved; their alignment with the abstract is still being audited (see task `audit-canonical-materials`). The render-pipeline artifacts that previously lived in `synthesis/` (`MANUSCRIPT_DRAFT.{md,typ}`, `STATS_AUDIT.md`, `render.log`, etc.) have all been relocated here as well — see "Second relocation batch" below.

**Date of relocation:** 2026-05-05
**Source:** repository top-level (158 files moved from `/moosefs/erikg/phrs/`)
**Move method:** `git mv` (history preserved for tracked files)

---

## Files moved, by category

### A. Copy-number-weighted ORA / hypergeometric / phyper methodology

Theory, implementation, validation, benchmarks, edge-case suites, and reports for the copy-number-weighted hypergeometric framework that was hill-climbed in the prior session.

- `base_r_phyper_benchmark.R`
- `benchmark_copy_weighted_hypergeometric.R`
- `benchmark_results.RData`
- `boundary_conditions_validation_report.R`
- `build_genome_wide_copy_background.py`
- `comprehensive_copy_background.csv`
- `comprehensive_parameter_mapping_validation_report.md`
- `comprehensive_r_phyper_synthesis_report.md`
- `comprehensive_statistical_validation.R`
- `comprehensive_statistical_validation_report.md`
- `comprehensive_validation_results.RData`
- `constraint_validation_tests.R`
- `constraint_violation_handler.R`
- `copy_aware_findings_summary.md`
- `copy_number_enrichment_methods_report.md`
- `copy_number_parameter_mapping_documentation.md`
- `copy_number_phyper_mapping.R`
- `copy_number_vs_standard_ora_comparison.csv`
- `copy_number_weighted_ora_best_practices_final.md`
- `copy_number_weighted_ora_best_practices_guide.md`
- `copy_number_weighted_ora_investigation.md`
- `copy_number_weighted_ora_methodology_synthesis.md`
- `copy_number_weighted_ora_performance_analysis.md`
- `copy_number_weighted_phyper_mathematical_formulation.md`
- `copy_weighted_enrichment.R`
- `copy_weighted_functional_analysis.csv`
- `copy_weighted_go_enrichment.R`
- `copy_weighted_hypergeometric.R`
- `copy_weighted_hypergeometric_documentation.md`
- `copy_weighted_ora_parameter_bounds.yaml`
- `copy_weighted_permutation_results.csv`
- `copy_weighted_vs_deduplicated_comparison.csv`
- `demo_copy_weighted_hypergeometric.R`
- `edge_case_analysis_comprehensive_report.md`
- `edge_case_test_results.md`
- `edge_case_test_suite.R`
- `enhanced_statistical_validation_framework.R`
- `enrichment_validation_tests.R`
- `error_handling_validation.R`
- `error_handling_validation_results.rds`
- `final_mathematical_validation_summary.md`
- `final_statistical_validation_report.md`
- `improved_copy_weighted_enrichment.R`
- `improved_copy_weighted_enrichment.csv`
- `improved_copy_weighted_vs_deduplicated_comparison.csv`
- `install_benchmark_packages.R`
- `integration_test_results.rds`
- `integration_test_suite.R`
- `integration_testing_final_report.md`
- `map_copy_number_research_completion_summary.md`
- `mathematical_formulation_copy_number_parameter_mapping.md`
- `mathematical_verification_report.md`
- `multi_copy_or_significance.md`
- `null_distribution_test.R`
- `null_distribution_validation_report.txt`
- `parameter_constraints_validation.R`
- `performance_benchmarks.R`
- `performance_summary.csv`
- `phyper_benchmark_detailed_results.rds`
- `phyper_benchmark_summary.csv`
- `phyper_computational_benchmark.R`
- `phyper_computational_benchmark_report.md`
- `phyper_parameter_modification_analysis.md`
- `protein_coding_enrichment_report.md`
- `r_alternatives_comparative_analysis.md`
- `r_alternatives_weighted_testing_research.md`
- `r_code_examples_and_templates.R`
- `r_implementation_guidelines.md`
- `r_phyper_modifications_research.md`
- `robust_copy_weighted_enrichment.R`
- `scalability_analysis_report.md`
- `simplified_performance_benchmark.R`
- `statistical_best_practices_weighted_ora.md`
- `statistical_validation_conclusions_and_recommendations.md`
- `statistical_validation_final_suite.R`
- `statistical_validation_framework_usage_guide.md`
- `statistical_validation_report.md`
- `statistical_validation_suite.R`
- `terminology_validation_report.md`
- `test_copy_weighted_hypergeometric.R`
- `test_null_validation_functions.R`
- `type_i_error_validation.R`
- `type_i_error_validation_report.md`
- `type_i_error_validation_results.RData`
- `type_i_error_validation_summary.md`
- `validation_examples_copy_number_mapping.R`
- `validation_report.md`
- `verification_system_bug_report.md`
- `verify_constraints.R`
- `verify_mathematical_validation.sh`
- `verify_null_distribution_validation.R`
- `verify_weighted_phyper_equivalence.R`
- `weighted_gene_enrichment_investigation.md`
- `weighted_phyper_mathematical_verification.md`
- `weighted_phyper_statistical_validation_report.md`
- `weighted_phyper_verification_summary.md`
- `workflow_compatibility_report.md`
- `workgraph_failure_modes_feedback.md`

### B. Gene-family-specific deep research (OR4F, DUX4/FRG2, miRNA, OR biology)

Literature deep-dives focused on gene families flagged by the off-target enrichment analysis. The canonical paper does not depend on per-family minutiae beyond the abstract's mention of the 4q–10q DUX4 clade.

- `olfactory_receptor_research.md`
- `or4f_functional_characterization.md`
- `deep_research_olfactory_receptors.md`
- `deep_research_synthesis.md`
- `deep_research_dux4_frg2.md`
- `deep_research_tubb8.md`
- `mirna_silencing_research_report.md`
- `phr_or_biology_connection.md`

### C. TUBB8 paralog deep research

- `tubb8_copy_literature.md`
- `tubb8_section1_literature_search.md`
- `tubb8_t2t_context.md`
- `tubb8b_copy_literature_section.md`
- `tubb8b_paralog.md`

### D. g:Profiler / GO / KEGG enrichment runs and outputs

API requests, results, parsers, dotplots, and per-arm GO enrichment tables.

- `compare_enrichment_simple.py`
- `comparison_table.csv`
- `gene_list_for_gprofiler_no_acro.txt`
- `gprofiler_request_coding_only.json`
- `gprofiler_request_no_acro.json`
- `gprofiler_results.json`
- `gprofiler_results_coding_only.json`
- `gprofiler_results_no_acro.json`
- `parse_coding_enrichment_results.py`
- `parse_gprofiler_results.py`
- `parse_gprofiler_results_no_acro.py`
- `phr_GO_BP_dotplot.txt`
- `phr_GO_BP_enrichment.csv`
- `phr_GO_MF_dotplot.txt`
- `phr_GO_MF_enrichment.csv`
- `phr_KEGG_enrichment.csv`
- `phr_coding_only_GO_BP_enrichment.csv`
- `phr_coding_only_GO_MF_enrichment.csv`
- `phr_copy_weighted_enrichment.R`
- `phr_copy_weighted_enrichment.csv`
- `phr_enrichment_summary.txt`
- `phr_gene_enrichment_report.md`
- `phr_gene_enrichment_synthesis.md`
- `phr_integration_examples.R`
- `phr_no_acro_GO_BP_enrichment.csv`
- `phr_no_acro_GO_MF_enrichment.csv`
- `run_coding_enrichment.sh`
- `run_coding_genes_enrichment.py`
- `run_gprofiler.sh`
- `create_plots.py`

### E. Off-target synthesis docs and supporting tables

Top-level synthesis documents from the prior session (executive summaries, decision frameworks, fact-checks, gene-copy tables) that are anchored on enrichment analysis rather than on the canonical implicit-pangenome-graph + 3D-organisation story.

- `EXECUTIVE_SUMMARY.md`
- `TODO.md`
- `andrea_phr_reconciliation.md`
- `decision_framework_section.md`
- `implementation_roadmap_section.md`
- `fact_check_report.md`
- `enriched_genes_per_arm.md`
- `enriched_genes_detailed_map.csv`
- `gene_copy_summary.csv`
- `genome_wide_gene_copies.csv`
- `all_gene_copies_by_arm.csv`

### F. GSEA figures and Excel summary outputs

- `Figure1_GSEA_BP_vertical.pdf`
- `Figure_GSEA_MF_vertical.pdf`
- `PHR_enrichment_summary.xlsx`
- `PHR_enrichment_all_results.xlsx`

### G. Wrapper scripts for off-target validation

- `R` (wrapper invoking `verify_null_distribution_validation.R`)
- `functions` (validation script for R null hypothesis tests)

---

## Verification (2026-05-05)

- 158 files moved (155 from category A–G inventory, plus the three wrappers).
- Pre/post SHA-256 of `paper_prep/synthesis/ABSTRACT.md`: `0771b2cd7b975bc209b8912c136fe28af510c7501cebe441c5e083ffa04c07af` (unchanged).
- `paper_prep/synthesis/REFERENCES.bib` left in place.
- `paper_prep/figures/` untouched.
- Top-level subtelomere data assets (`chm13.phrs.bed`, `chm13.phrs.no_acro.bed`, `CHM13-HG002.sub-telo-phrs.bed`, `chm13-annotations.bed`, `chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz`, `hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz`) left in place.
- Top-level canonical figure assets (`identity_heatmap_chr*.pdf`, `p_*.{pdf,png}`) left in place.
- Top-level canonical/uncertain-but-on-topic docs (`subtelomeric_analysis_report.{md,pdf}`, `framing_synthesis.md`, `plot-impg-coverage.R`, `phrs.*` gene reference files, `PHR_Subtelomeric Regions_Summary_March 2026.xlsx`, `T2T_Subtelomeric_Gene_Summary (1).xlsx`) left in place per the "leave when uncertain" rule.
- `paper_prep/synthesis/MANUSCRIPT_DRAFT.{md,typ}`, `MANUSCRIPT_SKELETON.md`, `ARCHITECT_TASK_BRIEF.md`, `WORK_DECOMPOSITION.md`, `ACCEPTANCE_CHECKLIST.md`, `STATS_AUDIT.md`, `SCRIPT_INVENTORY.md`, `TALK_OUTLINE_15MIN.md`, `VERSIONS.md`, `CAPTIONS.md`, `NOVEL_CONTRIBUTIONS.tsv`, `LIMITATIONS_X_FINDINGS.tsv` were inspected and found to anchor on the canonical abstract (41 arms → 15 communities, 232 individuals × 465 near-complete assemblies, Hi-C/Pore-C/CiFi/Dip-C/sperm/mouse meiotic, pedigree, out-of-Africa). They are kept in `paper_prep/synthesis/`.

---

## Second relocation batch (lead-author-directed)

**Date:** 2026-05-05 (later same day)
**Method:** `git mv` (history preserved)
**Reason:** The lead author (Erik Garrison) reviewed the rendered `MANUSCRIPT_DRAFT.pdf` and rejected it as "minutiae synthesis" and "from another dimension" relative to the canonical abstract. Despite the parking agent's judgment above that the synthesis docs reference canonical topics (41 arms, 465 near-complete assemblies, Hi-C, etc.), the lead author determined that these materials, as exposition, *do not deliver the actual paper described in `ABSTRACT.md`* — even when the topical surface looks right, the framing, priorities, and HPRC-v2-companion positioning are wrong. **Lead-author judgment overrides agent classification.** All synthesis-dir docs from the prior session (except `ABSTRACT.md` and `REFERENCES.bib`) were therefore relocated here.

Files relocated in this batch (from `paper_prep/synthesis/` to here):

- `ACCEPTANCE_CHECKLIST.md`
- `ARCHITECT_TASK_BRIEF.md`
- `CAPTIONS.md`
- `LIMITATIONS_X_FINDINGS.tsv`
- `BRAINSTORM.md` (formerly `MANUSCRIPT_DRAFT.md` — the off-target draft itself; renamed to clarify that this was a brainstorm of the prior session, NOT a draft of the canonical Nature companion)
- `BRAINSTORM.typ` (formerly `MANUSCRIPT_DRAFT.typ` — typst conversion of the brainstorm)
- `MANUSCRIPT_SKELETON.md`
- `NOVEL_CONTRIBUTIONS.tsv`
- `SCRIPT_INVENTORY.md`
- `STATS_AUDIT.md`
- `stats_audit/` (subdirectory: f7501 FDR, mantel multires)
- `TALK_OUTLINE_15MIN.md`
- `VERSIONS.md`
- `WORK_DECOMPOSITION.md`
- `pandoc_convert.log`
- `render.log`

After this batch, `paper_prep/synthesis/` contains only:
- `ABSTRACT.md` — canonical anchor (committed in same change as this relocation)
- `REFERENCES.bib` — citations database (will be reused by the canonical pipeline)

**Salvage potential:** non-zero. The topical references in these docs (e.g., the 41-arms / 15-communities analysis, Hi-C / Pore-C cross-checks, mantel multires correlation tests, per-arm-per-superpop Fisher tests) may correspond to legitimate analyses underlying the abstract's claims. The audit task `audit-canonical-materials` will catalog them; whether any specific text or analysis is salvaged into the canonical manuscript is a downstream decision.
