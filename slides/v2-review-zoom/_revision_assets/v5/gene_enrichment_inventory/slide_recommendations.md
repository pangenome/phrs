# Slide Recommendations

Recommended v5 section: 3 slides, plus one optional backup/methods slide if the presenter wants to discuss copy-weighted ORA.

## Slide 1: OR4F Is The Clean Copy-Rich Visual

Primary source:

- `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_enrichment.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_enrichment_fisher.tsv`

Talk claim:

OR4F gives the clearest visual entry point into copy-rich subtelomeric biology: 5,023 OR4F annotations across 16 arms, with a 62.1 percent overall pseudogene fraction and an 11.1 percent to 99.8 percent arm gradient.

Recommended visual:

- Reuse or update the existing OR4F pseudogenization gradient.
- Add one compact annotation line: 10 OR4F genes across 7 HPRCv2 communities; OR4F5 and OR4F8P each on 14 arms.
- Add one caveat in small text: OR family community presence is qualitative; C3 OR q is about 0.071, so no BH-significant family enrichment claim.

Why this slide:

It is the strongest signal with a canonical numeric plot source and already has review-zoom asset support. It also lets v5 introduce copy-rich gene biology without first explaining the entire community-enrichment pipeline.

## Slide 2: Duplicon Backbone And Diagnostic Communities

Primary source:

- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/cross_community_genes.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_families.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_specific_genes.tsv`

Talk claim:

The HPRCv2 communities recover a repeated subtelomeric duplicon backbone rather than isolated one-off genes: RPL23AP45 spans 10 communities/21 arms, SEPTIN14P22 spans 9/22, DDX11L16 spans 9/20, FAM138D spans 9/17, and WASH6P spans 8/21.

Recommended visual:

- A horizontal hub-gene bar or compact table with gene, biotype, communities, arms.
- Add two callout chips:
  - C1 chr4_q/chr10_q: 22 DUX4L pseudogenes specific to C1, D4Z4/FSHD-relevant.
  - C5 chr12_p/chr9_p/chr20_q/chr6_p: DDX11L/WASH/FAM138 module plus IQSEC3 gene-level anchor.

Required caveat:

For C5, `chr6_p` is a community member but lacks a called CHM13 PHR interval in `chm13.phrs.bed`. Label this as a community-arm signal unless a downstream figure explicitly renders only arms with called intervals.

Why this slide:

It explains the biology behind the OR4F slide: the same subtelomeric regions are dominated by recurrent duplicated families, not just olfactory receptors. It is also more robust than the parked copy-weighted GO outputs because it comes directly from the HPRCv2 community tables.

## Slide 3: Boundaries, Exceptions, And Non-Results

Primary source:

- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_summary_table.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_specific_genes.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_enrichment_fisher.tsv`
- `./chm13.phrs.bed`

Talk claim:

The community gene story has useful exceptions, and those exceptions keep the interpretation honest: C15/PAR1 is protein-coding rich with SHOX and GTPBP6, C7/acros carries an MTCO pseudogene signal, C4 has zero gene annotations despite community membership, and no family-community Fisher test survives BH correction.

Recommended visual:

- Three small rows or tiles:
  - C15 PAR1: 32.1 percent protein-coding; 10 specific genes including SHOX and GTPBP6.
  - C7 acrocentric p-arms: 48 specific genes including 9 MTCO mitochondrial pseudogenes; Fisher q about 0.071.
  - C4 chr7_q/chr12_q: zero genes, useful as the "community does not equal gene claim" caveat.
- A footer line: 116 family-community Fisher tests; no BH-significant row.

Required caveat:

The CHM13 called-interval extract does not cover every community arm. Missing community-arm rows include C5 `chr6_p`, C7 `chr13_p`, C14 `chrY_q`, and C15 `chrY_p`. Do not force those arms into CHM13 rendered-interval claims.

Why this slide:

This prevents overclaiming. It also gives downstream narrative authors a clean place to state that v5 is explaining qualitative gene-family architecture, not claiming a newly significant GO/ORA result from the HPRCv2 community Fisher tests.

## Optional Backup: Copy-Weighted ORA Needs Reanalysis

Primary source:

- `paper_prep/_brainstorming/improved_copy_weighted_enrichment.csv`
- `paper_prep/_brainstorming/improved_copy_weighted_vs_deduplicated_comparison.csv`
- `paper_prep/_brainstorming/gene_copy_summary.csv`
- `paper_prep/_brainstorming/copy_number_weighted_ora_methodology_synthesis.md`
- `paper_prep/_brainstorming/statistical_validation_conclusions_and_recommendations.md`

Talk claim:

There is a recoverable copy-weighted ORA workstream, but it should be shown only as an exploratory motivation unless re-run. It reports large fold enrichments, including 309x GTP binding/GTPase for GTPBP6/IQSEC3 and 928x transcription regulation for DUX4/FRG2/FRG2B, but the current weighted hypergeometric implementation is statistically caveated and the smell/olfactory mapping in the improved script is not reliable enough for an OR4F claim.

Recommended use:

- Keep out of the main 3-slide section unless the downstream figure task revalidates the method.
- If included, label it "exploratory copy-weighted ORA" and pair it with the HPRCv2 canonical gene support for GTPBP6/IQSEC3.
- Do not say "copy-weighted ORA proves OR4F olfactory enrichment" from the current parked tables. The canonical OR4F slide should instead use the HPRCv2 OR4F pseudogenization gradient and community presence metrics.

## Build Order For Downstream Tasks

1. Render Slide 1 first from the OR4F gradient plus HPRCv2 OR rows.
2. Render Slide 2 from `cross_community_genes.tsv` and `community_gene_families.tsv`.
3. Render Slide 3 as a compact caveat/exceptions slide using `community_summary_table.tsv`, `community_specific_genes.tsv`, `community_enrichment_fisher.tsv`, and `chm13.phrs.bed`.
4. Only then decide whether the optional copy-weighted ORA backup is worth reanalysis.
