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
- community/Jaccard check:
  `data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`,
  `paper_prep/manuscript_revision/C0_continuum/arm_pair_similarity_long.tsv`,
  and `paper_prep/manuscript_revision/C0_continuum/community_similarity_summary.tsv`

Outputs:

- `fig5_donor_recipient_ribbon_draft.svg`
- `fig5_donor_recipient_ribbon_draft.pdf`
- `fig5_donor_recipient_ribbon_draft.png`
- `donor_recipient_runs.tsv`
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

Build:

```bash
bash paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/scripts/make_donor_recipient_ribbons.sh
```
