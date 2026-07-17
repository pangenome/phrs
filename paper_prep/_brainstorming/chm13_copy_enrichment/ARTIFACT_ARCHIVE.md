# Chm13 Copy Enrichment Artifact Archive

This directory was reduced during repository size surgery on 2026-07-17.
The full pre-rewrite repository history, including all large intermediate
ontology-enrichment TSV/TSV.GZ ledgers and failed exploratory result layers,
is preserved in the private backup repository:

`git@github.com:ekg/phrs-backup.git`

Backup `main` was synchronized through the archive-manifest commits immediately
before surgery.

The cleaned `main` keeps the paper-facing V7 copy-number-aware ontology layer,
community attribution, functional-component heatmap artifacts and the compact
reports/manifests needed to understand the result. Bulky intermediate layers
removed from `main` history include:

- `analysis_ready/`
- `calibration_results/`
- `engine_terms/`
- `family_aware_v2/`
- `ontology_v3/`
- `ontology_v4/`
- `ontology_v5/`
- large V6 physical-copy/source/term ledgers
- `outputs/`
- `results/`
- `slurm_logs/`
- `sources/`
- `tests/`

These removals are repository-hygiene changes only. They do not alter the V7
final result tables, validation reports, community attribution outputs or
functional-component heatmap review artifacts retained on `main`.
