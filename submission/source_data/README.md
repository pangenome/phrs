# Source Data

One spreadsheet per figure, one sheet per panel, as the journal asks.

| file | sheets | written by |
|---|---|---|
| `Source_Data_Fig5.xlsx` | Fig5a (community ratios), Fig5b (cumulative divergence), Fig5c (per-arm ratio curves) and Fig5c_median, Fig5d (candidate graph dice inside and outside the tract) and Fig5d_windows (per 2 kb window) | `paper_prep/figures/fig5_concerted/figure_fig5_concerted.py` |
| `Source_Data_ED4.xlsx` | ED4: every alignment record of the dot plots (FastGA, identity = matches / block length) | `submission/scripts/make_pedigree_source_data.py` |
| `Source_Data_ED7.xlsx` | ED7B: 2 kb windows of panel B (released and corrected assemblies); ED7C: the 51 drawn non-homologous runs | `submission/scripts/make_pedigree_source_data.py` |

ED7A is the pedigree drawing and has no data. The grey same-chromosome chains of ED7C come from
`data/fig5_corrected_v11_PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz`.
Source Data for the other figures are not assembled yet.
