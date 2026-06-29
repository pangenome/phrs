# Fig5 homolog-vs-interchrom zoom panels

Zoomed subtelomeric panels for the corrected PAN027 paternal homolog-vs-
interchrom IMPG class-winner scan. The query is PAN027#2, the paternal hap2
assembly of PAN027, compared against PAN011#1+PAN011#2 as a combined father
target.

- input: `fig5_pre_impg_depth_filtered_similarity` 10:10 class-winner TSV
- window size: 2 kb query windows
- zoom size: 500 kb from the relevant telomere
- orientation: p-arm telomeres are shown at the left edge; q-arm telomeres are
  shown at the right edge with actual q-arm coordinate endpoints annotated
- break glyphs mark the omitted chromosome-body side of each zoom
- plotted signal: windows where the best interchromosomal IMPG match in the
  PAN011 father target beats the best same-chromosome/homolog match
- target summaries include the PAN011 target haplotype label (`h1`/`h2`) where
  available, for example `chr3 h2`
- panels: chrX p/PAR1, chr22 p, chr5 q, chr15 p, chr9 q

Build:

```bash
bash paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/scripts/make_zoom_panels.sh
```
