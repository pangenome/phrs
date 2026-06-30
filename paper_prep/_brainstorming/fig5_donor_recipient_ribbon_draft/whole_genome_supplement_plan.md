# Fig5 Whole-Genome Supplement Plan

## Recommendation

Use a supplementary whole-genome context figure now, but keep it
framed as context for the focused Fig. 5 donor-recipient ribbons rather than as
a new genome-wide recombination claim.

The recommended primary artifact is the existing PAN027 paternal
length-scaled whole-genome track:

- `paper_prep/_brainstorming/fig5_whole_genome_length_scaled_tracks/fig5_whole_genome_length_scaled_tracks.PAN027_hap2_paternal_vs_PAN011.pdf`
- `paper_prep/_brainstorming/fig5_whole_genome_length_scaled_tracks/fig5_whole_genome_length_scaled_tracks.PAN027_hap2_paternal_vs_PAN011.png`
- `paper_prep/_brainstorming/fig5_whole_genome_length_scaled_tracks/fig5_whole_genome_length_scaled_tracks.PAN027_hap2_paternal_vs_PAN011.svg`

This is preferable to adding only a manuscript note because it directly answers
the cherry-picking concern: the same child-recipient / father-donor comparison
used in Fig. 5 is shown in chromosome-scale context. It also matches the revised
visual grammar better than the earlier matrix-style whole-genome draft because
it uses length-scaled chromosome tracks, actual query coordinates, and a
PAN027 paternal <- PAN011 father page label. The focused Fig. 5 panel should
remain the main-text display because only the PAR1 positive control plus the two
autosomal PHR candidates are interpretable at 2 kb resolution with donor
interval ribbons.

Pair the figure with the compact source table generated here:

- `paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/pan027_paternal_high_conf_interchrom_winning_tracts.tsv`

The table is a conservative audit table from the existing 10:10 IMPG
class-winner tract summary. It keeps rows with:

- `basis == 10:10`
- tract length `bp >= 10000`
- `mean_inter_identity >= 0.95`

It labels the Fig. 5 contexts (`Fig5_PAR1_positive_control`,
`Fig5_chr5q_chr1p_candidate`, `Fig5_chr9q_chr3q_candidate`) and leaves other
rows as `acrocentric_context` or `whole_genome_context`. These labels are for
review/audit only and should not be used to make additional recombination
claims without separate follow-up.

## Inputs Inspected

Focused Fig. 5 ribbon and zoom inputs:

- `paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/fig5_donor_recipient_ribbon_draft.svg`
- `paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/fig5_donor_recipient_ribbon_draft.pdf`
- `paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/fig5_donor_recipient_ribbon_draft.png`
- `paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/donor_recipient_runs.tsv`
- `paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/zoom_window_segments.tsv`
- `paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/zoom_panel_summary.tsv`
- `paper_prep/_brainstorming/fig5_homolog_vs_interchrom_zoom_panels/zoom_phr_intervals.tsv`

Existing whole-genome / genome-context artifacts:

- `paper_prep/_brainstorming/fig5_whole_genome_length_scaled_tracks/README.md`
- `paper_prep/_brainstorming/fig5_whole_genome_length_scaled_tracks/length_scaled_track_manifest.tsv`
- `paper_prep/_brainstorming/fig5_whole_genome_length_scaled_tracks/length_scaled_track_segments.tsv`
- `paper_prep/_brainstorming/fig5_whole_genome_length_scaled_tracks/length_scaled_track_summary.tsv`
- `paper_prep/_brainstorming/fig5_untangle_whole_genome_overview/fig5_untangle_whole_genome_overview.pdf`
- `paper_prep/_brainstorming/fig5_untangle_whole_genome_overview/fig5_untangle_whole_genome_overview.png`
- `paper_prep/_brainstorming/fig5_untangle_whole_genome_overview/fig5_untangle_whole_genome_overview.svg`
- `paper_prep/_brainstorming/fig5_untangle_whole_genome_overview/untangle_whole_genome_segments.tsv`
- `paper_prep/_brainstorming/fig5_untangle_whole_genome_overview/untangle_whole_genome_summary.tsv`
- `paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome/fig5_homolog_vs_interchrom_whole_genome.pdf`
- `paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome/fig5_homolog_vs_interchrom_whole_genome.10to10.png`
- `paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome/fig5_homolog_vs_interchrom_whole_genome.10to10.svg`

Existing class-winner source tables:

- `paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/summaries/homolog_vs_interchrom_overall.tsv`
- `paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/summaries/homolog_vs_interchrom_pair_summary.tsv`
- `paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/summaries/homolog_vs_interchrom_top_tracts.tsv`

The length-scaled track manifest records the external source manifests used by
the existing plots:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency32/summaries/query_grid_chop_filter_manifest.tsv`
- `/moosefs/erikg/phrs/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/query_grid_filter_manifest.tsv`
- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_untangle_whole_genome_overview/untangle_whole_genome_segments.tsv`

## Rationale

The existing 10:10 corrected class-winner scan already shows the expected PAR1
positive control and the two autosomal windows used in Fig. 5. The overall
summary reports:

- `10:10`: 177,998 windows scanned, 174,942 with both same/homolog and
  interchromosomal matches, 1,778 windows where the best interchromosomal match
  beats the best same/homologous match, totaling 3.556 Mb.
- The two largest PAR1 tracts are `chrX -> chrY` at `chrX:14000-106000`
  (92 kb) and `chrX:110000-156000` (46 kb).
- The selected autosomal candidate tracts are represented in the high-confidence
  support table as `chr5 -> chr1` at `chr5:182052000-182080000` (28 kb) plus a
  second `chr5 -> chr1` tract at `chr5:182008000-182022000` (14 kb), and
  `chr9 -> chr3` at `chr9:136168000-136188000` (20 kb).

The scan also contains many acrocentric-context tracts and smaller scattered
signals. Those rows are useful for demonstrating that Fig. 5 was selected from
a genome-wide scan, but they should not be promoted into new biological claims
in the current main submission. The main text should preserve the existing
framing: PAR1 as the positive control, followed by two putative autosomal PHR
exchanges.

## Use In Submission

Recommended language for the figure/table callout:

> To place the focused Fig. 5 windows in genome-wide context, we also inspected
> the complete PAN027 paternal haplotype versus PAN011 father class-winner scan
> (Supplementary Fig. X; Supplementary Table Y). The genome-wide view recovers
> the expected Xp/Yp PAR1 signal and the two autosomal PHR candidate windows
> highlighted in Fig. 5, while additional interchromosomal-winning tracts are
> retained as context rather than interpreted as recombination candidates.

Do not add a new main-text statement such as "genome-wide paternal PHR
recombination is widespread in PAN027" from these artifacts alone.

## Later Work

If the submission needs a publication-polished whole-genome supplement, the next
task should redraw the PAN027 paternal page from
`fig5_whole_genome_length_scaled_tracks` with the exact Fig. 5 ribbon palette
and labels. That should still use the existing 2 kb class-winner and filtered
PAF outputs listed above; it should not rerun whole-genome alignment or IMPG.
