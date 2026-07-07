# Manuscript Promotion Decision

The maternal homologous-context ribbon views are folded into the main Fig. 5
manuscript figure as panels C and D, rather than being maintained as a separate
extended-data figure.

Current manuscript integration:

- The manuscript-facing combined PNG is
  `submission/fig/MainFigures/Fig5_whole_genome_recombination.png`.
- Panel B is the paternal PAN027-versus-PAN011 scan.
- Panel C is `PAN027mat_vs_PAN010_joint`.
- Panel D is `PAN028mat_vs_PAN027_joint`.
- The individual maternal PNGs remain copied under
  `submission/fig/ExtendedDataFigures` as source assets but are not referenced
  by `submission/paper.tex`.

The combined figure is assembled by
`paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft/scripts/compose_main_fig5.py`,
which crops the repeated per-panel legends and appends one shared legend.
