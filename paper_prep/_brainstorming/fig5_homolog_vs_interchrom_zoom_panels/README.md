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
- chromosome-end window: the gray bar spans the first/last 500 kb of the
  PAN027#2 chromosome from the alignment FASTA; the dark cap marks the physical
  chromosome start/end in that displayed window
- layout: panels are stacked in one centered column, with row-specific p/q
  orientation indicated by the telomere-side signal and break glyph; output is
  compacted to a 13.2 x 3.6 inch plotting area
- retained panels: chrXp/PAR1 control, chr5q, chr9q; chr22p and chr15p
  acrocentric rows are intentionally omitted from this view
- colors: best interchromosomal targets are collapsed to chrY, chr1, chr3, and
  other target
- PHR intervals: the thin bracket above each row is drawn from
  `zoom_phr_intervals.tsv`, which is regenerated from the WashU
  population-derived PHR table
  `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/all-vs-all.1Mb.p95.id95.len.tsv`.
  The table gives offsets within the 500 kb telomere-trimmed sequence; the
  script converts those offsets back to PAN027#2 full-chromosome coordinates
  using the sequence start embedded in the `seq` field, then clips/projects the
  intervals into the displayed PAN027#2 chromosome-coordinate zoom window.
  Bracket labels report the projected full-chromosome PHR/PAR coordinates, so
  terminal sequence outside the population-defined PHR remains visible when the
  PHR/PAR does not extend exactly to the chromosome end.
- break glyphs mark the omitted chromosome-body side of each zoom
- plotted signal: windows where the best interchromosomal IMPG match in the
  PAN011 father target beats the best same-chromosome/homolog match
- panels: chrX p/PAR1, chr5 q, chr9 q

Build:

```bash
bash paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/scripts/make_zoom_panels.sh
python3 paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/scripts/validate_centered_layout.py
```

The build script regenerates `zoom_window_segments.tsv`,
`zoom_panel_summary.tsv`, and `zoom_phr_intervals.tsv` when the upstream
class-winner gzip is available. In repo-only worktrees where that upstream
output is absent, it reuses the committed window/summary TSV snapshots, refreshes
`zoom_phr_intervals.tsv` directly from the WashU PHR table, and regenerates the
PDF/PNG/SVG layout artifacts.
