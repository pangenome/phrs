# Dip-C validation assets for review zoom v6

Task: `review-zoom-v6-dipc-validation-slides`

This directory stages the single-cell 3D validation assets added to the v6
review-zoom deck. The intent is to reuse already-rendered Dip-C/sperm outputs
where they exist, and only generate light summary plots from existing TSVs.
No Dip-C, hickit, or 3D modelling pipeline was re-run.

## Reproducible commands

Run from this directory:

```bash
./prepare_dipc_validation_assets.sh
```

The script:

- copies existing rendered PDFs into `source_pdfs/`;
- converts those PDFs to slide-ready PNGs in `pdf_pngs/` with
  `guix shell ghostscript poppler -- pdftoppm -r 220 -png -singlefile`;
- runs `make_dipc_validation_summary_plots.R` to create two summary plots in
  `plots/` from the TSV sources listed below;
- writes `source_manifest.tsv`, `conversion_log.tsv`,
  `plots/wb_negative_control_summary.tsv`, and
  `plots/community_free_rho_summary.tsv`.

## Staged rendered PDFs

| Asset | Original source | v6 use |
| --- | --- | --- |
| `gm12878_mantel_scatter` | `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_mantel_scatter.pdf` | Converted to PNG and used on slide `11a`. |
| `gm12878_radial_community` | `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_radial_community.pdf` | Converted to PNG and used on slide `11a`. |
| `sperm_all20_mantel_scatter` | `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_mantel_scatter.pdf` | Converted to PNG and used on slide `11a`. |
| `sperm_all20_radial_community` | `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_radial_community.pdf` | Converted to PNG and used on slide `11a`. |
| `sperm_all20_by_arm_type_arm` | `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/overlay_plots/sperm_all20.by_arm-type.arm.pdf` | Staged/converted for provenance; not placed in the concise v6 section. |
| `sperm_all20_by_arm_type_per_cell` | `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/overlay_plots/sperm_all20.by_arm-type.per-cell.pdf` | Staged/converted for provenance; not placed in the concise v6 section. |

## Summary TSV sources

The v6 summary plots use only existing TSV outputs:

| Plot | Inputs |
| --- | --- |
| `plots/wb_negative_control.png` | `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_summary.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_per_cell.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_per_community_per_cell.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_summary.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_per_cell.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_per_community_per_cell.tsv`. |
| `plots/community_free_rho_distribution.png` | `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_community_free_per_cell.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_community_free_arm.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_community_free_per_cell.tsv`; `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_community_free_arm.tsv`. |

## Metric checks

The slide text intentionally follows `end-to-end-report/report/06_dipc_validation.md`
and `paper_prep/surveys/SURVEY_06_dipc_validation.md`:

- GM12878 Dip-C: W/B = 0.931, 6.9% closer; Fisher p = 2.4e-05; Mantel rho =
  0.296, p = 0.002.
- Sperm scHi-C: W/B = 0.401, 60% closer; Fisher p = 3.9e-51; Mantel rho =
  0.202, p = 0.023.
- Negative control: non-sharing `S_all` arms are farther apart. GM12878
  `S_all` W/B = 1.106 with 0/16 cells below 1. Sperm `S_all` W/B = 1.397
  with 1/20 cells below 1.
- Community-free per-cell rho: GM12878 median rho = 0.093 with 15/16 cells
  positive; sperm median rho = 0.029 with 15/20 cells positive. The sperm
  pooled arm-level rho is shown as a caveat, not a headline positive result.

PBMC is intentionally not included in the v6 slide section because the report
states that PBMC community-free analysis is unavailable.
