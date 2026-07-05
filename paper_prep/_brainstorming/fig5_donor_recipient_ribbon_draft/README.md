# Fig5 donor-recipient ribbon

Ribbon-style view for the three PAN027 paternal Fig5 recipient windows.
This complements, rather than replaces,
`paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/`.

Inputs (all in-repo now — the figure is self-contained, no moosefs needed):

- recipient/window calls:
  `paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/zoom_window_segments.tsv`
- projected WashU population PHR/PAR intervals (the black bars):
  `paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/zoom_phr_intervals.tsv`
- WashU donor-population PHR length table:
  `data/fig5_washu.all-vs-all.1Mb.p95.id95.len.tsv`
- PAN011 (joint father target) chromosome lengths:
  `data/fig5_PAN027pat_vs_PAN011_joint.target.fa.fai`

The last two were vendored from moosefs
(`/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/all-vs-all.1Mb.p95.id95.len.tsv`
and the `pedigree_whole_genome_wfmash_p95_updated_bin` target `.fai`); override
the vendored copies with the `FIG5_PHR_TABLE` / `FIG5_TARGET_FAI` env vars to
point back at the moosefs sources.

Outputs:

- `fig5_donor_recipient_ribbon_draft.svg`
- `fig5_donor_recipient_ribbon_draft.pdf`
- `fig5_donor_recipient_ribbon_draft.png`
- `donor_recipient_runs.tsv`
- `caption.md`
- `methods.md`
- `COMMUNITY_LINKAGE_NOTE.md`

Geometry:

- Each top recipient track is the same PAN027#2 500 kb subtelomeric window used
  in the homolog-vs-interchrom zoom panel.
- Ribbons link each drawn recipient run to the PAN011 donor interval from the
  IMPG class-winning row (`inter_group_a`).
- Donor rows are displayed as unique 500 kb PAN011 windows containing one or
  more winning donor intervals. They are p-tip, q-tip, or local 500 kb windows
  depending on where the winning donor interval lies.
- Gray triangles mark the end of a displayed 500 kb window where the chromosome
  continues outside the plotted interval.
- Black bars show the population-defined PHR/PAR intervals on recipient and
  donor tracks when available.
- In the Xp/Yp panel, PAR1 is drawn across the full visible 500 kb p-tip
  window. The underlying PHR-call table only marks a 120 kb high-sharing call in
  this family, but that is not the full PAR1 annotation.
- Adjacent 2 kb windows are merged into a run when they share the same donor
  sequence/haplotype and remain adjacent in both recipient and donor space.
- Dominant donors for the three panels (`chrY`, `chr1`, `chr3`) are always
  drawn; other donor runs are drawn when they have at least 4 kb support. The
  `chr9_q -> chr7` row and redundant `chr5_q -> chr1 h1` donor row are
  suppressed for the clean display and flagged in `donor_recipient_runs.tsv`.

Regenerate:

```bash
# Needs: python3 (standard library only) and, for the PDF/PNG, rsvg-convert
# (librsvg, available via guix). Without rsvg-convert the script still writes
# the SVG and warns "SVG only: no rsvg-convert found."
bash paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/scripts/make_donor_recipient_ribbons.sh
```

Writes `fig5_donor_recipient_ribbon_draft.{svg,pdf,png}` and
`donor_recipient_runs.tsv` into this directory. Deterministic: the same inputs
reproduce byte-identical SVG/PNG.

Vendor the result into the manuscript:

```bash
cp paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/fig5_donor_recipient_ribbon_draft.pdf \
   submission/fig/MainFigures/Fig5_pedigree_untangle.pdf
```

- `caption.md` records the Fig. 5 caption language used in `submission/paper.tex`.
- `methods.md` records the matching direct-alignment workflow summary.
