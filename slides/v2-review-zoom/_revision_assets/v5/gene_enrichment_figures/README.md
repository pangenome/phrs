# Review Zoom v5 Copy-Number-Aware Gene Enrichment Figures

Task: `review-zoom-v5-copy-number-enrichment-figures`

Scope: slide-ready revision assets for the review-zoom v5 gene enrichment section. No deck source was edited.

## Assets

Run the generator from the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v5/gene_enrichment_figures/make_gene_enrichment_figures.R
```

Expected rendered assets:

| Asset | Purpose |
|---|---|
| `copy_aware_method_comparison.png` / `.pdf` | Compact method/comparison panel explaining why v5 uses HPRCv2 copy-aware support rather than standard or parked weighted ORA claims. |
| `ranked_copy_aware_gene_signals.png` / `.pdf` | Ranked candidate-signal panel for OR4F, DUX4L, duplicon backbone, C5, GTPBP6/IQSEC3, C15 PAR1, and C7 MTCO signals. |
| `community_family_signal_map.png` / `.pdf` | Community/family map linking selected signals to HPRCv2 communities and arms, with explicit CHM13 called-interval caveats. |

Supporting TSVs generated and used directly by the figure script:

| TSV | Role |
|---|---|
| `method_comparison_support.tsv` | Method/comparison lanes, source scope, and caveats. |
| `ranked_signal_support.tsv` | Ranked signal rows, plotted support-arm metric, key slide metric, interval note, caveat, and source paths. |
| `community_family_map_support.tsv` | Tile-level support table for the community/family map. |
| `community_interval_status.tsv` | Community-arm membership checked against `chm13.phrs.bed`; identifies community arms missing called CHM13 PHR intervals. |

## Canonical Sources Used

The generator reads these source tables directly:

| Path | Use |
|---|---|
| `slides/v2-review-zoom/_revision_assets/v5/gene_enrichment_inventory/source_inventory.tsv` | Inventory status and source classification context. |
| `slides/v2-review-zoom/_revision_assets/v5/gene_enrichment_inventory/candidate_enrichment_signals.tsv` | Ranked candidate signals, source paths, and caveats used in `ranked_signal_support.tsv`. |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv` | Canonical OR4F arm-level pseudogenization gradient: 5,023 annotations across 16 arms. |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_summary_table.tsv` | HPRCv2 arm-level community summaries, community arms, gene counts, biotype percentages, and C4/C15/C7 context. |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_enrichment.tsv` | Gene-by-community support rows for OR4F genes, IQSEC3, GTPBP6, SHOX/PAR1 anchors, and samples-in-community metrics. |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_families.tsv` | Community family support rows for OR, RPL, SEPTIN, DDX11L, WASH, FAM138, and related tiles. |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_enrichment_fisher.tsv` | Fisher/BH caveat: 116 community-family rows, 0 BH-significant rows; C3 OR and C7 MTCO q about 0.071. |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_specific_genes.tsv` | Community-specific DUX4L, MTCO, and C15 signal rows. |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/cross_community_genes.tsv` | Cross-community hub-gene support for RPL23AP45, SEPTIN14P22, DDX11L16, FAM138D, and WASH6P. |
| `chm13.phrs.bed` | Called CHM13 PHR interval gate used to label community arms that should not be overclaimed as called reference intervals. |

## Rejected Or Stale Sources

These were intentionally not used as final plotted truth:

| Path or source group | Decision |
|---|---|
| `paper_prep/_brainstorming/improved_copy_weighted_enrichment.csv` and `paper_prep/_brainstorming/improved_copy_weighted_vs_deduplicated_comparison.csv` | Parked exploratory weighted-ORA outputs. The method panel mentions this workstream only as context because upstream validation warns that the current weighted hypergeometric implementation is anti-conservative and not production-ready. |
| `paper_prep/_brainstorming/phr_GO_*.csv` and `paper_prep/_brainstorming/phr_coding_only_GO_*.csv` | Historical standard/deduplicated ORA outputs. They are not copy-number-aware and are used only as method contrast. |
| `paper_prep/_brainstorming/gene_copy_summary.csv` and `paper_prep/_brainstorming/all_gene_copies_by_arm.csv` | Useful historical copy-count context, but not the canonical HPRCv2 arm-level community source for this v5 figure bundle. |
| `/moosefs/guarracino/HPRCv2/PHR_III/sequence_level/enrichment/*` | Context only. These use the sequence-level partition, while this section narrates the 15-community HPRCv2 arm-level analysis. |
| `paper_prep/figures/ed4/figure_ed4.R` | Not reused. The inventory notes stale hard-coded worktree paths and mixed canonical/parked source decisions. |
| `slides/v2-review-zoom/_typst/assets/s14_or4f.png` and older OR4F crops | Not used as sources. They are deck/cache images rather than canonical numeric tables. |

## Interpretation Boundaries

The figures use labels such as "copy-aware support" and "copy-number-aware candidate signal" because the canonical HPRCv2 tables are support-aware arm/community summaries, not genome-wide copy-weighted ORA statistics. The ranked plot is ordered by the upstream inventory rank; bar length is a support-arm metric and should not be read as a q-value or final statistical effect size.

Interval-specific wording is intentionally restricted. The community tables describe HPRCv2 arm/community signals across pangenome PHR annotations. A community can contain arms that lack called CHM13 PHR rows in `chm13.phrs.bed`; the figure bundle calls out:

| Community | Missing called CHM13 PHR arm |
|---|---|
| C5 | `chr6_p` |
| C7 | `chr13_p` |
| C14 | `chrY_q` |
| C15 | `chrY_p` |

Use "community arm", "HPRCv2 arm-level support", or "pangenome HPR sequence support" for those signals unless a downstream panel explicitly renders only called CHM13 intervals.
