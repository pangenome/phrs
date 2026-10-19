# ED3 - subtelomeric module maps

Manuscript label `fig:ed3`; asset `submission/fig/ExtendedDataFigures/ED_Fig3_module_maps.pdf`.

This figure lives in `ed3_module_maps/` rather than `ed3/` because the legacy
`paper_prep/figures/ed3/` directory holds an unrelated, superseded figure
(TAR1 prevalence / internal TTAGGG islands / telomere length) that is no longer referenced by
`paper.tex`. Do not merge the two directories: the numbering here refers to the manuscript's
Extended Data ordering, not to the legacy directory names.

Regenerate: `Rscript paper_prep/figures/ed3_module_maps/figure_ed3_module_maps.R`
(input `module_tiling.tsv`, produced by `scripts/cladistics/arm_module_architecture.py`).
