# Fig5 homolog-vs-interchrom whole-genome view

Whole-genome visualization of the corrected IMPG class-winner scan for
`PAN027pat_vs_PAN011_joint`.

- input tracts: `../fig5_pre_impg_depth_filtered_similarity/summaries/homolog_vs_interchrom_top_tracts.tsv`
- input summaries: `../fig5_pre_impg_depth_filtered_similarity/summaries/homolog_vs_interchrom_overall.tsv`
- chromosome lengths: PAN027 paternal query FASTA `.fai`
- plotted signal: query windows where the best interchromosomal IMPG match
  beats the best same-chromosome/homologous match

Build:

```bash
bash paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome/scripts/make_whole_genome_plot.sh
```

