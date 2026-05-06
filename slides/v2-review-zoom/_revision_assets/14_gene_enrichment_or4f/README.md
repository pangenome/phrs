# Review Zoom 14: OR4F / Gene-Family Enrichment Inventory

Task: `review-zoom-14-gene-enrichment-or4f`
Date: 2026-05-06 UTC
Scope: revision assets only. No deck source was edited.

## Recommended Source

Use the regenerated candidate in this directory for the zoom review set:

- `or4f_gene_family_signal.png` / `.pdf`: slide-ready candidate visual.
- `or4f_gene_family_signal_table.md`: compact text-table version.
- `make_or4f_gene_family_signal.R`: reproducible local generator.

The recommended canonical source for the visible OR4F plot is:

`/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv`

The recommended canonical source for community gene-family context is the
arm-level HPRCv2 enrichment table set:

- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_enrichment.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_families.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_enrichment_fisher.tsv`
- `/moosefs/guarracino/HPRCv2/scripts/community/analyze-community-enrichments.R`

Use the old copy-number-aware ORA outputs only as historical context. They are
recoverable and useful for explaining why Erik previously emphasized OR4F and
other copy-rich families, but commit `c50b99c` explicitly parked them under
`paper_prep/_brainstorming/` as noncanonical relative to the current paper
pipeline.

## Candidate Claim

OR4F is the clean slide visual: 5,023 annotations across 16 subtelomeric arms
form an 11.1% to 99.8% pseudogenization gradient. The broader gene-family
analysis supports the same biology as a qualitative PHR-content signal:
OR-family rows are present across communities, DUX4L marks the C1/D4Z4
community, and RPL23A/SEPTIN/DDX11L families form a copy-rich subtelomeric
duplicon backbone. Do not say the OR family is BH-significant at the
community-family Fisher level; the best OR row is C3 with q = 0.071.

## Output Files

| File | Purpose |
|---|---|
| `or4f_gene_family_signal.png` | 16:9 slide-ready PNG: OR4F pseudogenization gradient plus compact gene-family signal cards. |
| `or4f_gene_family_signal.pdf` | Vector/PDF counterpart for review or future deck integration. |
| `or4f_gene_family_signal_table.md` | Slide-ready table and one-line claim. |
| `or4f_gene_family_signal_table.tsv` | Machine-readable backing table generated from source data. |
| `make_or4f_gene_family_signal.R` | Generator. Reads HPRCv2 OR4F/community TSVs and the parked Erik copy-summary table. |

## Provenance Timeline

| Commit | Date UTC | Task / agent | What matters here | Classification |
|---|---:|---|---|---|
| `340e61d` | 2026-03-31 20:52 | initial PHR enrichment data import | Added Angela-style enrichment workbooks, GSEA PDFs, CHM13 interval/annotation inputs, and `subtelomeric_analysis_report.md`. | Historical input. |
| `0d4c53f` | 2026-03-31 21:16 | `step-3-run-go` | Added repo-root g:Profiler GO CSVs and parser scripts. | Historical ORA output. |
| `847ec37` | 2026-03-31 22:06 | `rerun-enrichment-excluding` | Added no-acrocentric GO rerun outputs. | Historical ORA output. |
| `724ec49` | 2026-04-01 12:54 | `go-enrichment-on`, agent-61 | Added protein-coding-only GO outputs and found 23 coding genes, including 4 OR4F coding genes. | Historical ORA output. |
| `6dcec23` | 2026-04-02 01:34 | `implement-copy-number`, agent-73 | Added `gene_copy_summary.csv`, `all_gene_copies_by_arm.csv`, `enriched_genes_per_arm.md`, and copy-number-aware enrichment scripts/tables. | Recoverable but now parked. |
| `ccc152f` | 2026-04-02 19:22 | `executive-summary-combined`, agent-576 | Added the comprehensive `EXECUTIVE_SUMMARY.md` with the OR4F/DUX4/IL9R copy-family narrative. | Historical synthesis. |
| `e245ba2` | 2026-04-02 23:56 | `reformat-gene-family`, agent-585 | Reformatted the gene-family table in `EXECUTIVE_SUMMARY.md`. | Historical synthesis. |
| `c50b99c` | 2026-05-05 23:26 | `park-off-target`, agent-744 | Moved the enrichment hill-climb artifacts into `paper_prep/_brainstorming/` and labeled them not canonical for the current manuscript. | Stale-path warning. |
| `438fb02` / `88954d7` | 2026-05-05 05:41 / 05:42 | `figure-ed3-ed4-annotation-genes`, agent-697 | Built ED4, including high-copy gene families and OR4F gradient. The committed ED4 script has a stale absolute `REPO_ROOT` into agent-697. | Useful render, script needs path audit before reuse. |
| `863e756` / `6cf052c` | 2026-05-06 00:49 / 00:50 | `bog-v2-slide-14`, agent-774 | Added `slides/v2/slide_14_gene_biology.R/.md/.pdf/.png`; the R script reads the HPRCv2 OR4F source directly. | Best committed slide-level source. |
| `cb973ab` / `fd9a250` | 2026-05-06 02:16 / 02:17 | `build-bog-v2-3`, agent-878 | Created zoom crops `slide_14a_dux4.png`, `slide_14b_or4f.png`, `slide_14c_tar1.png`. Crop recipe not committed. | Stale/intermediate crops. |
| `10bee88` / `4862ec7` | 2026-05-06 19:21 / 19:24 | `render-bog-annotated-zoom`, agent-951 | Copied slide 14 assets into the current review-zoom deck. `s14_or4f.png` is an exact copy of the agent-878 crop. | Current deck cache, not source. |
| `a2ac7d1` | 2026-05-06 20:40 | `review-zoom-git-provenance-audit` | Established the current deck lineage and warned that slide 14 crops are non-reproducible. | Required upstream audit. |
| `cda3b92` / `6d1136b` | 2026-05-06 20:56 | `review-zoom-14-gene-background`, agent-972 | Produced slidelet guidance that preserves the OR4F gradient and documents the same crop risk. | Related revision guidance. |

## Source Inventory

| Path | Contents | Provenance / status | Recommendation |
|---|---|---|---|
| `slides/v2-review-zoom/_revision_assets/git_provenance/README.md` | Upstream deck asset audit. Slide 14 rows identify current OR4F as a copied agent-878 crop and recommend rebuilding from source. | Required first-read dependency; commit `a2ac7d1`. | Cite for deck lineage. |
| `slides/v2/slide_14_gene_biology.R` | Three-panel slide 14 generator for DUX4, OR4F, TAR1. Reads HPRCv2 off-tree tables directly. | Added by `863e756`, merged as `6cf052c`. | Best committed slide-level source for the existing OR4F gradient. |
| `slides/v2/slide_14_gene_biology.md` | Existing slide 14 bullets, metrics, speaker notes, and source paths. | Added by `863e756`, merged as `6cf052c`. | Use for slide copy and denominator sanity checks. |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv` | OR4F per-arm total, pseudogene count, coding count, pseudogene fraction. 17 lines; mtime 2026-03-19 01:13 UTC. | Off-tree HPRCv2 analysis product; local HPRCv2 directory is not a git checkout. | Canonical numeric source for OR4F gradient. |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_enrichment.tsv` | Gene-by-community presence, arm counts, samples, biotypes. 1,160 lines; mtime 2026-03-28 21:57 UTC. | Output of HPRCv2 community enrichment pipeline. | Canonical for OR4F genes across communities and hub-gene arm/community counts. |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_families.tsv` | Community-level gene-family rollups including OR, RPL, SEPTIN, DDX11L, WASH, MTCO. 131 lines; mtime 2026-03-28 21:57 UTC. | Output of HPRCv2 community enrichment pipeline. | Canonical qualitative gene-family context. |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_enrichment_fisher.tsv` | 116 community-family Fisher tests with BH correction. All `significant = FALSE`. | Output of HPRCv2 community enrichment pipeline. | Use for caveat: OR presence is qualitative, not BH-significant. |
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze-community-enrichments.R` | Script producing the HPRCv2 community enrichment TSVs. | Off-tree script; mtime 2026-03-28 20:11 UTC. | Method anchor for canonical tables. |
| `subtelomeric_analysis_report.md` | Section 9 summarizes community gene enrichment, OR4F, DUX4L, hub genes, and Fisher caveat. | Added by `4adebed`; not part of the parked `_brainstorming` directory. | Best repo-local prose source for the HPRCv2 enrichment narrative. |
| `paper_prep/surveys/SURVEY_03_gene_enrichment.md` | Survey of `end-to-end-report/report/03_gene_enrichment.md`, including source TSV inventory and figure gaps. | Added by `3db7023`. | Use as a compact repo-local guide to HPRCv2 section 9. |
| `paper_prep/figures/ed4/caption.md` | ED4 caption: high-copy gene families and OR4F gradient. | Added by `88954d7`. | Useful wording, but confirm source paths. |
| `paper_prep/figures/ed4/figure_ed4.R` | ED4 generator. Includes a stale absolute `REPO_ROOT <- "/moosefs/erikg/phrs/.wg-worktrees/agent-697"`. | Added by `88954d7`. | Do not reuse without path repair. |
| `paper_prep/_brainstorming/gene_copy_summary.csv` | Recovered Erik copy-number table used here for contextual copy counts: OR4F = 72, DUX4/FRG2/FRG2B = 54, IL9R/IL9RP = 58. | Added by `6dcec23`, moved by `c50b99c` into parked `_brainstorming`. | Context only, not recommended as canonical slide source. |
| `paper_prep/_brainstorming/EXECUTIVE_SUMMARY.md` | Comprehensive older enrichment summary and compact gene-family catalog. | Added by `ccc152f`, reformatted by `e245ba2`, moved by `c50b99c`. | Recoverable historical synthesis; do not treat as canonical without caveat. |
| `slides/v2-review-zoom/_typst/assets/s14_or4f.png` | Current visible review-zoom OR4F asset. | Added by `10bee88`, merged as `4862ec7`; exact blob match to agent-878 crop `slide_14b_or4f.png`. | Avoid as a source; it is a crop cache. |

## Locked Metrics

| Metric | Value | Source |
|---|---:|---|
| OR4F annotations | 5,023 total | `or4f_pseudogene_fraction.csv` |
| OR4F pseudogene / coding split | 3,117 pseudogene / 1,906 coding | `or4f_pseudogene_fraction.csv` |
| OR4F pseudogene fraction range | 11.1% at chr7p to 99.8% at chr15q | `or4f_pseudogene_fraction.csv` |
| OR4F overall pseudogene fraction | 62.1% | `or4f_pseudogene_fraction.csv` |
| OR4F genes across communities | 10 OR4F genes in 7 communities | `community_gene_enrichment.tsv` |
| Most widespread OR4F genes | OR4F5 and OR4F8P on 14 arms each | `community_gene_enrichment.tsv`; also `subtelomeric_analysis_report.md` section 9.2 |
| Community-family Fisher result | 116 tests; no BH-significant family-community enrichments; C3 OR q = 0.071 | `community_enrichment_fisher.tsv` |
| D4Z4 comparator | C1 chr4q/chr10q; 22 DUX4L pseudogenes specific to C1 | `subtelomeric_analysis_report.md` section 9.5 |
| Duplicon backbone comparator | RPL23AP45 10 communities / 21 arms; SEPTIN14P22 9 / 22; DDX11L16 9 / 20 | `community_gene_enrichment.tsv`; `subtelomeric_analysis_report.md` section 9.6 |

## Regeneration

From repo root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/14_gene_enrichment_or4f/make_or4f_gene_family_signal.R
```

Expected outputs:

- `or4f_gene_family_signal.png` at 3200 x 1800 pixels.
- `or4f_gene_family_signal.pdf`.
- `or4f_gene_family_signal_table.md`.
- `or4f_gene_family_signal_table.tsv`.

The script intentionally keeps the old copy-number table as a labeled context
row rather than making it the plot source. This preserves Erik's recovered
enrichment work without letting parked intermediates override the HPRCv2 source
tables.
