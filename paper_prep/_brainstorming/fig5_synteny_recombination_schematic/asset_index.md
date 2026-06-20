# Fig5 Synteny Schematic Review Asset Index

This index lists the inspectable assets relevant to the Fig5 synteny recombination schematic review pack. The schematic pack is SVG/table-only in this worktree; PDF conversion was unavailable during prototype generation.

## Schematic Prototype Pack

- `fig5_synteny_recombination_full.svg` - Full chromosome/arm-context schematic for the PAR1 positive control and two strict-path chr9q/chr3q PHR candidates, with neutral chromosome-length context and source-to-child flow ribbons.
- `fig5_synteny_recombination_focus.svg` - Focused 500 kb native-window schematic for the same three events, using a consistent physical scale and labeled flow segments.
- `event_manifest.tsv` - Event-level manifest selecting exactly the PAR1 positive control, PAN027 chr9q PHR candidate, and corrected PAN028 strict chr9q PHR candidate, with coordinate provenance and donor/side-fragment summaries.
- `selected_segments.tsv` - Segment-level strict primary-path table used for drawing geometry, including native query intervals, recovered target-side intervals, event roles, identity/Jaccard, and community annotations when exactly joined.
- `coordinate_provenance.md` - Audit note documenting that coordinates are native sample assembly window coordinates, not CHM13-projected coordinates.
- `pdf_conversion_status.txt` - Conversion status showing that `fig5_synteny_recombination_full.pdf` and `fig5_synteny_recombination_focus.pdf` were not generated because no local SVG-to-PDF converter was available.

## Comparison Assets

- `../fig5_sweepga_1to1_redraw/fig5_sweepga_1to1_redraw.pdf` - Existing strict `nb=1` sweepGA 1:1 compact redraw for comparison against the schematic's event selection and color-track interpretation.
- `../fig5_sweepga_1to1_redraw/fig5_sweepga_1to1_redraw.svg` - SVG version of the same strict redraw, useful for inspecting labels and geometry when PDF rendering is inconvenient.
- `../fig5_sweepga_1to1_redraw/conservative_segments.tsv` - Source strict segment table feeding the schematic manifest and selected segment table; this is the geometry source to compare against.
- `../fig5_sweepga_1to1_redraw/summary_counts.tsv` - Filter/count audit for native n4, native `nb=1`, and conservative sweepGA 1:1 stages.
- `../fig5_sweepga_1to1_redraw/phr_intervals.tsv` - PHR overlay intervals used by the strict redraw, useful only as context for the candidate tracks.
- `../fig5_sweepga_1to1_redraw/validation_report.tsv` - Redraw validation table for query-length, coordinate, `nb`, and rendered-rectangle checks.
- `../fig5_par1_phr_candidate_panels/fig5_par1_phr_candidate_panels.pdf` - Superseded compact PAR1/PHR candidate panel; useful only as a warning comparison because panel C is misleading under the strict manifest.
- `../fig5_par1_phr_candidate_panels/fig5_par1_phr_candidate_panels.svg` - SVG version of the superseded compact panel; do not use panel C as the schematic target.
- `../fig5_par1_phr_candidate_panels/panel_event_summary.tsv` - Superseded panel summary showing old panel C as PAN028 chr3q with chr7p, chr16q hap1, and chr20q hap2 side fragments, not the corrected PAN028 chr9q event.

## Superseded Panel C Caveat

The older compact asset `../fig5_par1_phr_candidate_panels/fig5_par1_phr_candidate_panels.{pdf,svg}` is superseded for panel C. Under strict `nb=1` sweepGA 1:1 data, the old PAN028 chr3q panel does not map chr3q to chr9q; it shows side fragments from chr7p, chr16q hap1, and chr20q hap2. The corrected schematic event is `PAN028_chr9q_chr3q_PHR_candidate`: PAN028 chr9q query with chr3q primary donor and chr15q plus chr16q side fragments, as recorded in `event_manifest.tsv` and `selected_segments.tsv`.
