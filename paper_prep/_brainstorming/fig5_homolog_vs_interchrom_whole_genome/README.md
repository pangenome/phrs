# Fig5 homolog-vs-interchromosomal whole-genome view

This directory contains whole-genome visualizations built from the corrected
IMPG class-winner summaries in
`paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/summaries/`
for `PAN027pat_vs_PAN011_joint`.

The canonical plot for this task is
`fig5_homolog_vs_interchrom_whole_genome.pdf`. It shows each query chromosome at
actual CHM13 length. Colored tracts are merged 2 kb query windows where the best
interchromosomal IMPG hit has higher estimated identity than the best
same-chromosome/homolog hit. Color encodes the winning target chromosome.
Separate stacked panels show the `1:1`, `4:4`, and `10:10` SweepGA/FastGA
mapping bases from the corrected summaries.

Context overlays:

- PAR1/PAR2 intervals from `data/chm13-annotations.bed` on chrX/chrY.
- Acrocentric PHR intervals from `data/chm13-annotations.bed` on
  chr13/14/15/21/22.
- PAN027 and PAN028 chr9q-to-chr3q candidate windows from the Fig5 scaffold-jump
  sensitivity config.
- Centromeres from `data/chm13-annotations.bed`.

Run from the repository root:

```bash
Rscript paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome/scripts/plot_homolog_vs_interchrom_whole_genome.R
```

or regenerate and validate all task outputs:

```bash
bash paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome/validate_outputs.sh
```

Canonical outputs:

- `fig5_homolog_vs_interchrom_whole_genome.pdf`
- `fig5_homolog_vs_interchrom_whole_genome.png`
- `fig5_homolog_vs_interchrom_whole_genome.svg`
- `context_intervals.tsv`
- `top_pair_summary_for_plot.tsv`

This directory also preserves the independently generated 10:10-only upstream
asset set:

- `fig5_homolog_vs_interchrom_whole_genome.10to10.png`
- `fig5_homolog_vs_interchrom_whole_genome.10to10.svg`
- `scripts/make_whole_genome_plot.sh`
- `scripts/plot_whole_genome_homolog_vs_interchrom.R`
