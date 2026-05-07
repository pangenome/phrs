# Review Zoom v5 Gene Enrichment Narrative

Task: `review-zoom-v5-copy-number-enrichment-narrative`

Scope: speaker-usable narrative text for a review-zoom v5 section on
copy-number-aware gene enrichment. This directory is a documentation handoff
for downstream rendering; no deck source was edited here.

## Recommended Framing

Use the HPRCv2 arm-level community enrichment tables and the OR4F arm-gradient
table as the production sources for v5. The section should explain
copy-number-aware interpretation in plain language: standard ORA collapses each
gene symbol to one hit, while PHR/subtelomeric biology is dominated by repeated
gene instances, duplicated families, pseudogenes, and repeat-adjacent
annotations.

The older genome-wide copy-weighted ORA workstream is useful as motivation and
as a possible backup, but it should not be presented as validated production
evidence without reanalysis. Prior validation found the weighted
hypergeometric implementation anti-conservative under gene-level sampling, with
failed null calibration and FDR concerns. Use it only as an exploratory
contrast unless a calibrated permutation or empirical-null analysis is run.

## Canonical Source Decisions

Primary v5 evidence:

- `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_summary_table.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_enrichment.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_families.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/cross_community_genes.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_specific_genes.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_enrichment_fisher.tsv`

Repo-local inventory source:

- `slides/v2-review-zoom/_revision_assets/v5/gene_enrichment_inventory/README.md`
- `slides/v2-review-zoom/_revision_assets/v5/gene_enrichment_inventory/candidate_enrichment_signals.tsv`
- `slides/v2-review-zoom/_revision_assets/v5/gene_enrichment_inventory/slide_recommendations.md`

Exploratory ORA caveat sources:

- `paper_prep/_brainstorming/improved_copy_weighted_enrichment.csv`
- `paper_prep/_brainstorming/improved_copy_weighted_vs_deduplicated_comparison.csv`
- `paper_prep/_brainstorming/copy_number_weighted_ora_methodology_synthesis.md`
- `paper_prep/_brainstorming/statistical_validation_conclusions_and_recommendations.md`

## Required Interval Caveat

Community-arm assignment is not identical to a called CHM13 PHR interval in
`chm13.phrs.bed`. The v5 text should say "community arms", "HPRCv2 arm-level
support", or "pangenome HPR sequence support" unless the rendered interval is
explicitly present in the CHM13 BED file.

Community-assigned arms without called CHM13 PHR interval rows include:

- C5 `chr6_p`
- C7 `chr13_p`
- C14 `chrY_q`
- C15 `chrY_p`

## Files

- `slide_text.md`: four-slide talk-section proposal with concise titles,
  subtitles, and bullets.
- `method_box.md`: slide-ready plain-language explanation of copy-number-aware
  enrichment.
- `speaker_notes.md`: speaker notes, biological interpretation, and caveats.

## Bottom Line

The strongest v5 narrative is not "we discovered significant GO enrichment."
It is: subtelomeric PHRs are copy-rich, and the HPRCv2 community tables show
that repeated gene-family architecture is a real part of the signal. OR4F is
the clean visual entry point; DDX11L/WASH/FAM138/RPL23A define the recurring
duplicon backbone; IQSEC3/GTPBP6 explain why GTP-related terms are plausible
but still require validated copy-weighted statistics.
