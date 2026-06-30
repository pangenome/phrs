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
- scale: per-row coordinate endpoints are drawn on the tracks; a single 100 kb
  scale bar replaces the full local-position axis
- layout: panels are stacked in one centered column, with row-specific p/q
  orientation indicated by the telomere-side signal and break glyph
- retained panels: chrXp/PAR1 control, chr5q, chr9q; chr22p and chr15p
  acrocentric rows are intentionally omitted from this view
- colors: best interchromosomal targets are collapsed to chrY, chr1, chr3, and
  other target
- PHR span: the shaded box behind each row marks the range covered by the
  plotted interchromosome-over-homolog windows in that 500 kb subtelomeric view
- break glyphs mark the omitted chromosome-body side of each zoom
- plotted signal: windows where the best interchromosomal IMPG match in the
  PAN011 father target beats the best same-chromosome/homolog match
- panels: chrX p/PAR1, chr5 q, chr9 q

Build:

```bash
bash paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/scripts/make_zoom_panels.sh
python3 paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/scripts/validate_centered_layout.py
```

The build script regenerates `zoom_window_segments.tsv` and
`zoom_panel_summary.tsv` when the upstream class-winner gzip is available. In
repo-only worktrees where that upstream output is absent, it reuses the
committed TSV snapshots and regenerates only the PDF/PNG/SVG layout artifacts.
