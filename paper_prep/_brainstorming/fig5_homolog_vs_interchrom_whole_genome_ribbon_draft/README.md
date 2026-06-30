# Fig5 whole-genome 10:10 ribbon draft

Draft whole-genome ribbon view for the corrected `PAN027pat_vs_PAN011_joint`
SweepGA/F32 10:10 IMPG class-winner scan.

This uses the same source data as the Fig5 zoom/ribbon and the corrected
whole-genome homolog-vs-interchrom overview:

`paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs/PAN027pat_vs_PAN011_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz`

Design:

- top genome: PAN027 paternal haplotype child query
- lower genomes: PAN011 father donor haplotype 1 and haplotype 2
- chromosomes are concatenated in chromosome order with actual chromosome-length
  scaling within each genome track
- ribbons show adjacent 2 kb query windows where the best interchromosomal IMPG
  match beats the best same-chromosome/homologous match
- drawn ribbons are filtered to runs of at least 10 kb and mean interchromosomal
  identity of at least 0.95

Build:

```bash
bash paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft/scripts/make_whole_genome_ribbon_draft.sh
```

Outputs:

- `fig5_homolog_vs_interchrom_whole_genome_ribbon_draft.svg`
- `fig5_homolog_vs_interchrom_whole_genome_ribbon_draft.pdf`
- `fig5_homolog_vs_interchrom_whole_genome_ribbon_draft.png`
- `fig5_homologous_recombination_context_ribbon_draft.svg`
- `fig5_homologous_recombination_context_ribbon_draft.pdf`
- `fig5_homologous_recombination_context_ribbon_draft.png`
- `whole_genome_ribbon_runs.tsv`
- `whole_genome_ribbon_summary.tsv`
- `whole_genome_homologous_context_runs.tsv`
- `whole_genome_homologous_context_summary.tsv`

The homologous-context variant uses the same geometry but adds the paired
same-chromosome competitor from each interchrom-over-same IMPG window as a
light-gray ribbon layer. The colored interchromosomal/non-homologous ribbons are
preserved on top.
