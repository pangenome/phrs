# Fig5 Synteny Schematic Visual Review

## Short Review

1. **Full context plus focused windows:** Yes. `fig5_synteny_recombination_full.svg` gives chromosome/arm context for all three selected events, while `fig5_synteny_recombination_focus.svg` zooms into the 500 kb native windows with 100 kb scale bars and labeled segment flows. The pair should be reviewed together; the full view is context, the focus view is the interpretable recombination-window view.

2. **Three strict primary-path events:** Mostly yes, with one important correction relative to older assets. The schematic manifest contains exactly the PAR1 positive control, PAN027 chr9q PHR candidate, and PAN028 chr9q PHR candidate. The second PHR event checked here is `PAN028_chr9q_chr3q_PHR_candidate`, not the older compact PAN028 chr3q panel: `event_manifest.tsv` reports PAN028 maternal hap1 chr9q (`chr9:134380985-134880985`) with 22 strict support rows, 499901 bp total support, 449356 bp same-chromosome chr9q context, 34172 bp chr3q primary donor, and chr15q/chr16q side fragments of 15166 bp and 1207 bp. `selected_segments.tsv` agrees with this role breakdown and has no chr7p or chr20q rows for the corrected PAN028 event.

3. **Coordinate clarity:** The prototype is explicit that coordinates are native sample assembly window coordinates and 0-based half-open intervals, not CHM13 projections. This is documented in `coordinate_provenance.md`, repeated in `README.md`, and visible in the focus SVG labels. The manuscript-facing caption or panel label must preserve this, because the chromosome labels alone could otherwise be mistaken for CHM13/reference coordinates.

4. **Flow/ribbon labels:** The flow labels are clearer than the older color-track interpretation. The legend and segment labels distinguish same-chromosome context, PAR1 positive-control flow, primary PHR donor, side fragment, and low-confidence tail. The focus view is especially useful because the labels sit on the native 500 kb windows; the full view is less self-sufficient because many ribbons are necessarily compressed.

5. **Misleading title/label/scale issues to fix before manuscript integration:** Do not present the full view as cytoband-accurate or as a precise whole-chromosome scale drawing. It uses `data/chm13.chrom.sizes` for approximate chromosome-length context and neutral schematic bands, and some native sample coordinates extend beyond CHM13 size lines. The title "strict-path synteny/recombination schematic" is acceptable for brainstorming, but a manuscript panel should avoid implying a validated clean crossover; the autosomal rows remain candidate terminal PHR exchange-compatible patches with same-chromosome context and side fragments. Also keep the "native assembly coordinates, not CHM13" caveat in the figure text or caption.

6. **Recommendation:** Use the focus view as the primary candidate for manuscript integration, optionally paired with a simplified full-context inset if space allows. Keep the combined full/focus pair in brainstorming until the coordinate caveat and schematic-scale caveat are carried into the caption. Do not reuse the older compact panel C.

## Superseded Comparison Asset

`../fig5_par1_phr_candidate_panels/fig5_par1_phr_candidate_panels.{pdf,svg}` is superseded and misleading for panel C. Under strict `nb=1` sweepGA 1:1, the old PAN028 chr3q panel does not show a chr3q-to-chr9q event; `panel_event_summary.tsv` shows chr7p, chr16q hap1, and chr20q hap2 side fragments. The corrected schematic event is `PAN028_chr9q_chr3q_PHR_candidate`, with PAN028 chr9q as the query, chr3q as the primary donor, and chr15q plus chr16q side fragments according to `event_manifest.tsv`.
