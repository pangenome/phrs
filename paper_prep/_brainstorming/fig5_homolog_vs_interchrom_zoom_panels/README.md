# Fig5 homolog-vs-interchrom zoom panels

Zoomed subtelomeric panels for the corrected PAN027 paternal homolog-vs-
interchrom IMPG class-winner scan.

- input: `fig5_pre_impg_depth_filtered_similarity` 10:10 class-winner TSV
- window size: 2 kb query windows
- zoom size: 500 kb from the relevant telomere
- plotted signal: windows where the best interchromosomal IMPG match beats the
  best same-chromosome/homolog match
- panels: chrX p/PAR1, chr22 p, chr5 q, chr15 p, chr9 q

Build:

```bash
bash paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/scripts/make_zoom_panels.sh
```

