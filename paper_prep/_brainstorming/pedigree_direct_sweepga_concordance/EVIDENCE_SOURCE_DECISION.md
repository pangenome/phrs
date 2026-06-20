# Fig5 pedigree evidence-source decision

Date: 2026-06-20

## Decision

Use a hybrid evidence policy for Fig5/pedigree. Direct haplotype-to-parent haplotype `sweepGA` is cleaner than the graph-derived path for the PAR1 positive control and for the large, interpretable blocks of the two chr9q/chr3q autosomal candidates, but it is not clean enough to replace graph/`odgi untangle` globally. The graph-derived `event_manifest.tsv` and `selected_segments.tsv` should remain the primary event definition and provenance source; direct `sweepGA` should become the primary visual/alignment source only for segments that are concordant with those graph calls.

Overall recommendation: graph/untangle remains primary for event boundaries and fragment interpretation; direct `sweepGA` becomes primary support/control evidence for clean agreeing blocks and should be shown as the cleaner alignment validation layer.

## Evidence reviewed

- Direct raw many:many/no-scaffold PAFs: `raw_paf/*.sweepga_many_many_j0.paf.gz`
- Direct filtered PAFs: `filtered_paf/*.four_many_noscaffold.paf.gz`, with `one_one`, `one_many`, `two_many`, and `simple_i95_l1k_q80` retained as controls.
- Direct-vs-graph concordance: `summaries/direct_vs_graph_concordance.tsv`
- Direct file summaries: `summaries/paf_file_summary.tsv`
- Direct plots: `plots/direct_sweepga_concordance_focused.{svg,pdf}`, `plots/direct_sweepga_full_genome_overview.{svg,pdf}`, and `plots/direct_sweepga_concordance_review.svg`
- Graph event definition: `../fig5_synteny_recombination_schematic/event_manifest.tsv`
- Graph selected strict segments: `../fig5_synteny_recombination_schematic/selected_segments.tsv`
- Prior strict sweepGA/untangle redraw: `../fig5_sweepga_1to1_redraw/conservative_segments.tsv`

The direct concordance table contains 38 selected-segment comparisons: 28 agree, 9 are discordant to another target arm/haplotype, and 1 is partial/inconclusive. Counts alone are not sufficient; the discordance is concentrated in small side fragments, donor slivers, and one low-confidence tail, while the large same-chromosome context and major chr3q donor intervals usually agree.

## Event decisions

| Event | Direct concordance | Clarity of donor/product | Ambiguity or multimapping | Agreement with graph | Coordinate interpretability | Visualization suitability | Evidence-source decision |
|---|---:|---|---|---|---|---|---|
| `PAR1_XY_positive_control` | 6 agree / 0 discordant; 499,675 bp agreeing | 5/5. Alternating PAN011 chrX/chrY parental haplotype blocks recover the expected male PAR structure. | 5/5. No discordant selected rows. | 5/5. Direct calls reproduce all selected graph intervals. | 5/5. Native PAN027 chrX query coordinates and PAN011 chrX/chrY target windows are straightforward. | 5/5. Direct PAF is suitable for both full-window and focus visualization. | **Direct sweepGA primary** for this control; graph files remain provenance. |
| `PAN027_chr9q_chr3q_PHR_candidate` | 8 agree / 2 discordant; 498,300 bp agreeing and 1,699 bp discordant | 4/5. The large chr9q context and chr3q donor rows are clear; chr15q side fragment is also directly recovered. | 3/5. chr16q side fragment and tiny chr20q tail are not recovered as expected; direct best hits stay on chr9p/chr9q. | 4/5 overall. Agreement covers nearly the full 500 kb candidate, but not every side/tail fragment. | 4/5. Direct PAF coordinates are native and interpretable; graph selected segments still provide better fragment labels and patch annotations. | 4/5. Direct works well for full-genome/focus panels for chr9q and chr3q; graph should annotate discordant side/tail fragments. | **Hybrid: direct sweepGA primary for concordant chr9q/chr3q and chr15q blocks; graph primary for chr16q side fragment and chr20q low-confidence tail.** |
| `PAN028_chr9q_chr3q_PHR_candidate` | 14 agree / 7 discordant / 1 partial; 449,226 bp agreeing, 32,344 bp discordant, 947 bp partial | 3/5. Broad chr9q context is clear, and several chr3q donor rows agree, but donor slivers split across PAN027 haplotypes are less clean. | 2/5. Discordance affects chr16q/chr15q side fragments, several small chr3q donor rows, and the terminal chr3q row; direct best hits often remain on chr9q or go to chr19p. | 3/5. The broad graph candidate is supported, but segment-level replacement would lose important graph-derived fragment calls. | 3/5. Native coordinates are usable, but the source relationship is harder to explain from direct PAF alone because short graph-defined fragments are absorbed by longer direct alignments. | 3/5. Direct is suitable for full-window context and selected clean blocks; focus visualization needs graph segment labels/provenance to avoid overstating certainty. | **Hybrid, leaning graph primary: direct sweepGA supports the event and can show clean large blocks, but graph/untangle remains primary for event interpretation and ambiguous fragments.** |
| `PAN027_maternal_vs_PAN010_control/comparison` | No selected graph event rows; raw and filtered PAFs exist for both PAN010 haplotypes. | Not scoreable at event level. | Not scoreable at event level. | Not scoreable because `selected_segments.tsv` has no corresponding event. | File-level native coordinates are available. | Useful as a file-level control only, not a Fig5 event panel source yet. | **Support/control only unless a graph-selected candidate is added later.** |

## Segment-level policy

Use `summaries/direct_vs_graph_concordance.tsv` as the decision table:

- Rows labeled `direct_sweepga_can_be_primary_for_segment` may drive direct alignment geometry in review or figure-redraw work.
- Rows labeled `keep_graph_primary_pending_manual_review` should keep `../fig5_synteny_recombination_schematic/selected_segments.tsv` as the source of truth for that fragment.
- Rows labeled `inconclusive_keep_graph_primary` should not be promoted to direct-primary evidence.

For autosomal Fig5, do not collapse the event into direct-only PAF geometry. The direct many:many/four:many alignments often expose the dominant longer homologous block cleanly, but the graph-derived strict path is still better for naming short side fragments, preserving patch-pattern annotations, and explaining why a terminal sliver is a candidate fragment rather than background subtelomeric multimapping.

## Recommended Fig5 data-source stack

Primary event definition:

- `../fig5_synteny_recombination_schematic/event_manifest.tsv`
- `../fig5_synteny_recombination_schematic/selected_segments.tsv`

Primary direct evidence for clean rows:

- `summaries/direct_vs_graph_concordance.tsv`
- `filtered_paf/PAN027pat_vs_PAN011_hap1.four_many_noscaffold.paf.gz`
- `filtered_paf/PAN027pat_vs_PAN011_hap2.four_many_noscaffold.paf.gz`
- `filtered_paf/PAN028mat_vs_PAN027_hap1.four_many_noscaffold.paf.gz`
- `filtered_paf/PAN028mat_vs_PAN027_hap2.four_many_noscaffold.paf.gz`

Support/control direct evidence:

- `raw_paf/*.sweepga_many_many_j0.paf.gz`
- `filtered_paf/*.one_one_noscaffold.paf.gz`
- `filtered_paf/*.one_many_noscaffold.paf.gz`
- `filtered_paf/*.two_many_noscaffold.paf.gz`
- `filtered_paf/*.simple_i95_l1k_q80.paf.gz`
- `filtered_paf/PAN027mat_vs_PAN010_hap*.paf.gz` as maternal comparison controls with no current selected event.

Prior strict sweepGA/untangle provenance:

- `../fig5_sweepga_1to1_redraw/conservative_segments.tsv`
- `/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_sweepga_1to1_noscaffold/<pair>.e50000.m1000.j0.8.n4.native_nb1.sweepga_1to1_noscaffold.paf`
- `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv` for community/status labels only.

## Notes for downstream figure work

1. For the PAR1 control, a direct-only panel is defensible because every selected segment agrees and the donor/product relationship is visually simple.
2. For PAN027, direct `four_many_noscaffold` PAFs should be used to draw the main chr9q context and chr3q donor evidence. Keep graph-derived labels for chr16q and the chr20q tail, or mark them as graph-only/low-confidence rather than implying direct confirmation.
3. For PAN028, direct PAFs are best used as validation overlays or simplified block evidence. A direct-primary panel would under-represent several graph-defined donor/side fragments; keep the graph-selected segment table as the backbone.
4. The direct package summaries currently contain absolute paths from the upstream worktree that generated the package. For durable scripts, resolve files relative to this directory and use the relative paths listed in `asset_manifest.tsv`.
