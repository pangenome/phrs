# v7 Gene Enrichment Report-Backed Assets

Task: `review-zoom-v7-gene-enrichment-report-backed`

Scope: slide-ready text and tables for the review zoom v7 gene-enrichment/GO section. No final deck source was edited.

## Bottom Line

The v7 deck can say that HPRCv2 subtelomeric communities carry recurrent gene-family architecture: OR4F/olfactory receptor blocks, C1 D4Z4/DUX4L, the C5 DDX11L/WASH/FAM138/IQSEC3 module, C15 PAR1 genes such as SHOX/GTPBP6/P2RY8, C7 acrocentric MTCO pseudogenes, and a pseudogene/ncRNA-rich duplicon backbone.

The v7 deck should not say that BH-significant GO enrichment proves any of those biology points. The canonical HPRCv2 community-family Fisher screen has 116 rows and 0 BH-significant rows. The strongest rows have BH q = 0.07118; C3 OR and C7 MTCO are near-miss/candidate presence patterns, not definitive enriched classes.

Use wording such as:

- "recurrent gene-family architecture"
- "copy-aware candidate signal"
- "community-arm support"
- "report-backed presence pattern"

Avoid wording such as:

- "BH-significant GO enrichment proves..."
- "validated copy-weighted ORA enrichment"
- "definitive enriched GO class"
- "OR4F GO enrichment from the current canonical Fisher screen"

## Canonical Sources

Primary report text:

- `subtelomeric_analysis_report.md:477-600`, section 9.
- `end-to-end-report/report/03_gene_enrichment.md:1-130`, the report-extracted version of the same section.

Primary HPRCv2 tables listed by both reports:

- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_summary_table.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_enrichment.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_families.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_specific_genes.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/cross_community_genes.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_enrichment_fisher.tsv`

Historical/context GO tables:

- `paper_prep/_brainstorming/phr_GO_BP_enrichment.csv`
- `paper_prep/_brainstorming/phr_GO_MF_enrichment.csv`
- `paper_prep/_brainstorming/phr_coding_only_GO_BP_enrichment.csv`
- `paper_prep/_brainstorming/phr_coding_only_GO_MF_enrichment.csv`
- `paper_prep/_brainstorming/phr_copy_weighted_enrichment.csv`
- `paper_prep/_brainstorming/improved_copy_weighted_enrichment.csv`

## Claim Boundaries

### What Can Be Claimed

1. The canonical analysis grouped gene annotations by 15 arm-level communities and counted shared gene instances across arms. It reports 374 unique genes across 39 arms, 576 shared gene instances within communities, 93 community-specific genes, and 216 genes found in 2 or more communities. Sources: `subtelomeric_analysis_report.md:481-483`; `end-to-end-report/report/03_gene_enrichment.md:5-7`.

2. Subtelomeric gene content is mostly pseudogene/ncRNA-rich, with PAR1 as the protein-coding outlier. C15 has 32.1% protein-coding content; most other populated communities are roughly 4-9% protein-coding. Sources: `subtelomeric_analysis_report.md:487-507`; `community_summary_table.tsv` rows 2, 4, 6, 8, and 16 for C1, C3, C5, C7, and C15 examples.

3. OR4F/olfactory receptor architecture is a recurrent presence pattern. The report says 10 OR4F family genes occur across 7 communities, with OR4F5 and OR4F8P the most widespread at 14 arms each. C3 has an OR family row across all 6 C3 arms, but its Fisher row is not BH-significant. Sources: `subtelomeric_analysis_report.md:509-515`; `community_gene_families.tsv` row 73; `community_gene_enrichment.tsv` rows 152, 153, 170, and 198; `community_enrichment_fisher.tsv` row 4.

4. C1 chr4_q/chr10_q has a D4Z4/DUX4L signature. The report maps Ambrosini block 7 to C1 and describes 22 C1-specific DUX4L pseudogenes. Sources: `subtelomeric_analysis_report.md:531` and `:573`; `community_specific_genes.tsv` rows 5-26; `community_summary_table.tsv` row 2.

5. C5 supports a DDX11L/WASH/FAM138/IQSEC3 module. The report maps Ambrosini family C to C5 and notes DDX11L family members across 3-4 arms; IQSEC3 is detected in C5 on chr12_p with 453 samples. Sources: `subtelomeric_analysis_report.md:535`, `:551`, and `:562`; `community_gene_families.tsv` rows 91, 93, and 94; `community_gene_enrichment.tsv` row 270.

6. C15 PAR1 is a protein-coding outlier with biologically interpretable genes. The report lists 10 C15-specific genes, including SHOX, GTPBP6, and P2RY8, and reports 32.1% protein-coding content. Sources: `subtelomeric_analysis_report.md:505`, `:507`, and `:574`; `community_specific_genes.tsv` rows 33-42; `community_gene_enrichment.tsv` rows 1134, 1139, and 1142.

7. C7 acrocentric p-arms carry a candidate MTCO pseudogene signal. The report lists C7-specific MTCO genes and describes acrocentric/rDNA-adjacent context. The MTCO family row spans all 5 C7 arms, but the Fisher row is a BH near-miss. Sources: `subtelomeric_analysis_report.md:572`; `community_gene_families.tsv` row 103; `community_specific_genes.tsv` rows 60-64, 73-75, and 82; `community_enrichment_fisher.tsv` row 5.

8. The cross-community backbone is mostly duplicated pseudogene/ncRNA content. RPL23AP45 spans 10 communities/21 arms, SEPTIN14P22 spans 9/22, DDX11L16 spans 9/20, FAM138D spans 9/17, and WASH6P spans 8/21. Sources: `subtelomeric_analysis_report.md:580-598`; `cross_community_genes.tsv` rows 2, 4, 6, 10, and 11.

### What Cannot Be Claimed

1. Do not claim a BH-significant HPRCv2 gene-family enrichment result. The canonical Fisher table has 116 rows, 0 `significant=TRUE` rows, and minimum BH q = 0.07118. Sources: `community_enrichment_fisher.tsv` rows 2-5 and full-row count; report conclusion at `subtelomeric_analysis_report.md:598`.

2. Do not claim the historical standard GO/ORA tables are copy-number-aware. They collapse repeated subtelomeric copies to unique gene symbols or to a small coding-only gene list. They can be shown as context or method contrast only.

3. Do not claim the copy-weighted ORA p-values or fold enrichments are final statistics. The older method reports striking fold enrichments in `improved_copy_weighted_enrichment.csv`, but the methodology files flag anti-conservative behavior and failed FDR control. Sources: `paper_prep/_brainstorming/copy_number_weighted_ora_methodology_synthesis.md:13-19`, `:67-104`, and `:388-418`; `paper_prep/_brainstorming/statistical_validation_conclusions_and_recommendations.md:50-56` and `:86-94`.

4. Do not call the exploratory 598x smell/olfactory copy-weighted rows an OR4F validation. The improved table's smell/olfactory rows list IL9R/IL9RP genes, not OR4F genes. Source: `paper_prep/_brainstorming/improved_copy_weighted_enrichment.csv` rows 2-3.

## Files In This Bundle

- `gene_family_function_slide_table.tsv`: compact slide-ready table of report-backed gene-family/function signals, claim scope, and source rows.
- `gene_family_function_slide_table.md`: markdown rendering of the same table for direct slide drafting.
- `go_function_context_table.tsv`: compact table of historical GO/ORA context rows with explicit labels.
- `go_function_context_table.md`: markdown rendering of the GO/method contrast table.
- `slide_text_blocks.md`: 3 slide-ready text blocks with speaker-note hooks.
- `gene_enrichment_report_backed_summary.svg`: optional compact visual summary panel.
- `SLIDE_PATCH.md`: recommended v7 deck patch and speaker-note language for the fan-in renderer.
- `asset_manifest.tsv`: manifest of the bundle contents.

## Recommended Use

For v7, replace or tighten the current gene-enrichment section around slides 14m-14c. A clean version is two slides:

1. "Gene enrichment: report-backed biology, conservative statistics" using `gene_family_function_slide_table.md` or the SVG summary panel.
2. "GO/ORA context: useful contrast, not final proof" using `go_function_context_table.md`.

The existing browser backup panels for D4Z4 and OR4F can remain as backup visuals if slide count allows, but their notes should inherit the same caveat: they illustrate recurrent architecture and copy-aware support, not BH-significant GO proof.
