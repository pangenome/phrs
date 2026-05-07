# Review Zoom v5 Gene Enrichment Inventory

Task: `review-zoom-v5-copy-number-enrichment-inventory`

Scope: audit the copy-number-aware and gene-family enrichment sources that could feed a new review-zoom v5 section. This directory is an inventory only; no deck source was edited.

## Bottom Line

Use the HPRCv2 arm-level community enrichment TSVs as the canonical source for review-zoom v5 gene-family claims:

- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_summary_table.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_enrichment.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_families.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_enrichment_fisher.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_specific_genes.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/cross_community_genes.tsv`
- `/moosefs/guarracino/HPRCv2/scripts/community/analyze-community-enrichments.R`

Use `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv` as the canonical numeric source for the OR4F pseudogenization gradient.

Treat `paper_prep/_brainstorming/*copy*`, `*weighted*`, and older `*enrichment*` artifacts as parked exploratory work unless a downstream task re-runs and re-validates them. They are useful for understanding why copy number matters, but they should not override the HPRCv2 community enrichment tables in the v5 deck.

## Canonical Versus Exploratory

The current manuscript/report gene-enrichment section is an HPRCv2 community analysis, not a genome-wide copy-weighted ORA section. Its statistics are arm/community-level summaries: gene presence by community, gene-family rollups, community-specific genes, cross-community hub genes, and Fisher tests over community x gene-family tables. `paper_prep/surveys/SURVEY_03_gene_enrichment.md` lines 83-100 inventory these arm-level TSVs, and lines 125-136 describe the deployed method.

The old copy-weighted ORA stream is different. It uses CHM13/no-acro PHR gene lists, per-copy gene tables, and a genome-wide copy background to compare standard/deduplicated ORA against copy-weighted hypergeometric tests. `paper_prep/surveys/SURVEY_DATA_inventory.md` lines 50-61 cataloged those tables as attractive figure inputs, but later review-zoom work parked the artifacts under `paper_prep/_brainstorming/` as noncanonical relative to the current paper pipeline (`slides/v2-review-zoom/_revision_assets/14_gene_enrichment_or4f/README.md` lines 27-31).

The methodological caveat is substantial. `paper_prep/_brainstorming/copy_number_weighted_ora_methodology_synthesis.md` lines 13-19 says weighted ORA is mathematically defined but anti-conservative under gene-level sampling, and lines 67-104 show inflated Type I error and failed FDR control. `paper_prep/_brainstorming/statistical_validation_conclusions_and_recommendations.md` lines 48-56 and 84-91 explicitly warn against deploying the current weighted hypergeometric implementation for production use. One concrete data-quality warning: `paper_prep/_brainstorming/improved_copy_weighted_enrichment.R` lines 55-67 maps sensory smell and olfactory receptor activity to IL9R/IL9RP genes, while the prose surveys often describe the same 598x signal as OR4F-driven. Do not render that as a validated OR4F GO result without reanalysis.

## Copy-Number Awareness

There are three distinct meanings of "copy-number-aware" in these sources:

1. HPRCv2 community enrichment tables are copy/support-aware at the pangenome/community level. They count arms, samples, annotations, and gene-family occurrences across 465 haplotype/reference sequences. They are canonical for slide claims about OR4F presence, DUX4L, RPL23A/DDX11L/WASH/FAM138 duplicon backbones, MTCO, SHOX/GTPBP6, and TAR1 context.
2. Standard/deduplicated ORA tables (`phr_GO_*`, `phr_coding_only_*`) collapse genes to unique symbols or coding lists. They are not copy-number-aware, but they are useful as contrast rows.
3. Copy-weighted ORA tables (`improved_copy_weighted_*`, `gene_copy_summary.csv`, `comprehensive_copy_background.csv`) explicitly weight by copies against a genome-wide background. These are the only true genome-wide copy-number-weighted ORA outputs, but they are parked/exploratory and methodologically caveated.

## PHR-Interval Specificity Caveat

Separate arm/community membership from called/rendered PHR intervals. The HPRCv2 community tables assign arms to communities across the pangenome. The CHM13 reference interval extract `chm13.phrs.bed` has 37 called arm intervals, not every community arm. Among the arm-level HPRCv2 communities, the CHM13 extract lacks called intervals for:

- C5: `chr6_p`
- C7: `chr13_p`
- C14: `chrY_q`
- C15: `chrY_p`

Therefore, a signal can be valid as an arm/community signal without being valid as a CHM13 rendered-interval claim for every arm in that community. Downstream slides should label these as "community arms" or "pangenome HPR sequence support" unless the plotted interval is explicitly present in `chm13.phrs.bed`.

## Slide-Worthy Signals

Highest-confidence signals for v5:

- OR4F pseudogenization gradient: 5,023 OR4F annotations across 16 arms; 3,117 pseudogene and 1,906 coding annotations; 62.1 percent pseudogene overall; 11.1 percent at chr7p to 99.8 percent at chr15q.
- OR family community signal: 10 OR4F genes across 7 communities; OR4F5 and OR4F8P each on 14 arms. This is a qualitative presence pattern; C3 OR Fisher q is about 0.071, not significant.
- D4Z4/DUX4L: C1 chr4q/chr10q has 26 community-specific genes, including 22 DUX4L pseudogenes, and 86.4 percent pseudogene content.
- Duplicon backbone: RPL23AP45 spans 10 communities/21 arms; SEPTIN14P22 spans 9/22; DDX11L16 spans 9/20; FAM138D spans 9/17; WASH6P spans 8/21.
- C5 module: DDX11L/WASH/FAM138 family rows plus IQSEC3 in C5. Caveat that C5 arm `chr6_p` is a community member but lacks a called CHM13 interval in `chm13.phrs.bed`.
- C15 PAR1: 10 community-specific genes including SHOX and GTPBP6; 32.1 percent protein-coding content. Caveat that `chrY_p` lacks a called CHM13 interval row.
- C7 acrocentric MTCO: 48 C7-specific genes, including 9 MTCO mitochondrial pseudogenes; Fisher MTCO p is small but q is about 0.071. Caveat that `chr13_p` lacks a called CHM13 interval row.
- TAR1: useful as repeat context, not a gene-enrichment result. C2 has the highest mean TAR1 count per sequence (2.51); overall TAR1 has 21,544 entries with 66.9 percent within 10 kb of the telomere.

Important non-results/caveats:

- HPRCv2 arm-level Fisher tests: 116 community-family tests, no BH-significant result.
- C4 chr7q/chr12q has zero gene annotations despite being a community/minimal-PHR signal.
- Copy-weighted GO fold-enrichment values (598x, 928x, 309x, 825x) are exploratory until revalidated with a calibrated/permutation method and corrected GO mappings.

## Files In This Directory

- `source_inventory.tsv`: source-by-source audit with status, scope, use/reject decision, and caveat.
- `candidate_enrichment_signals.tsv`: ranked signal list with source, statistic, support, community/arms, PHR-interval specificity, copy-awareness, and caveat fields.
- `slide_recommendations.md`: recommended 3-slide v5 section plus optional backup material.

## Recommended Next Action

For the downstream figure task, start with `candidate_enrichment_signals.tsv` ranks 1-4 and render from the HPRCv2 arm-level tables plus the OR4F pseudogene-fraction table. Keep the copy-weighted ORA rows as a backup/methods-caveat panel unless a downstream agent re-runs a statistically calibrated copy-number-aware analysis.
