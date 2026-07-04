# Manuscript Promotion Decision

Decision: do not expand `submission/paper.tex` or copy these maternal PNGs into
`submission/fig/ExtendedDataFigures` in this commit.

Rationale:

- The active submission currently defines a single Extended Data figure,
  `fig:ed1`, for sequence-to-3D contact replicates.
- The repository instructions say only `fig:ed1` and `fig:fig1` through
  `fig:fig5` exist in the current submission draft and warn not to re-add
  `fig:ed2` or later references without an explicit figure-scheme change.
- These maternal whole-genome ribbons are draft Extended Data companions to
  Fig. 5, generated under `paper_prep/_brainstorming`, not yet part of the
  submission figure plan.
- The current assets are suitable for manual promotion after author decision:
  they include SVG, PDF, PNG, run TSVs, summaries, merge audits, conversion
  status, a draft caption, and method notes.

Manual promotion path, if approved later:

1. Decide whether these two maternal comparisons should be a new Extended Data
   figure or folded into a revised Fig. 5 supplement.
2. Copy the selected PNG or PDF assets from
   `paper_prep/_brainstorming/fig5_extended_maternal_whole_genome_ribbons` into
   `submission/fig/ExtendedDataFigures` with final names.
3. Add a new figure environment and label in `submission/paper.tex`, updating
   all cross-references consistently.
4. Rebuild the manuscript from `submission/` with `make` or `bash compile.sh`.

Current candidate figure assets:

- `PAN027mat_vs_PAN010_joint/PAN027mat_vs_PAN010_joint.whole_genome_ribbon.png`
- `PAN027mat_vs_PAN010_joint/PAN027mat_vs_PAN010_joint.whole_genome_homologous_context_ribbon.png`
- `PAN028mat_vs_PAN027_joint/PAN028mat_vs_PAN027_joint.whole_genome_ribbon.png`
- `PAN028mat_vs_PAN027_joint/PAN028mat_vs_PAN027_joint.whole_genome_homologous_context_ribbon.png`

The PDF and SVG siblings are present beside each PNG.
