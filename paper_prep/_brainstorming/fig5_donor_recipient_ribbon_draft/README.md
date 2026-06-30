# Fig5 donor-recipient ribbon draft

Draft ribbon-style view for the three PAN027 paternal Fig5 recipient windows.
This complements, rather than replaces,
`paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/`.

Inputs:

- recipient/window calls:
  `paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/zoom_window_segments.tsv`
- projected WashU population PHR/PAR intervals:
  `paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/zoom_phr_intervals.tsv`
- PAN011 target chromosome lengths:
  `/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027pat_vs_PAN011_joint.target.fa.fai`

Outputs:

- `fig5_donor_recipient_ribbon_draft.svg`
- `fig5_donor_recipient_ribbon_draft.pdf`
- `fig5_donor_recipient_ribbon_draft.png`
- `donor_recipient_runs.tsv`

Geometry:

- Each top recipient track is the same PAN027#2 500 kb subtelomeric window used
  in the homolog-vs-interchrom zoom panel.
- Ribbons link each drawn recipient run to the PAN011 donor interval from the
  IMPG class-winning row (`inter_group_a`).
- Donor rows are displayed as unique 500 kb PAN011 windows containing one or
  more winning donor intervals. They are p-tip, q-tip, or local 500 kb windows
  depending on where the winning donor interval lies.
- Adjacent 2 kb windows are merged into a run when they share the same donor
  sequence/haplotype and remain adjacent in both recipient and donor space.
- Dominant donors for the three panels (`chrY`, `chr1`, `chr3`) are always
  drawn; other donor runs are drawn when they have at least 4 kb support.
  Remaining 2 kb singleton "other" hits are summarized in text below the panel.
  The per-run TSV records both the run and the donor 500 kb window used for
  plotting.

Build:

```bash
bash paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/scripts/make_donor_recipient_ribbons.sh
```
