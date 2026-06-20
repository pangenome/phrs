# Coordinate Provenance Audit

## Inputs audited

The schematic manifest uses the strict primary-path segment table only for drawing geometry:

- `paper_prep/_brainstorming/fig5_sweepga_1to1_redraw/conservative_segments.tsv`

That table is documented by `paper_prep/_brainstorming/fig5_sweepga_1to1_redraw/README.md` as native `odgi untangle` PAF rows filtered to `nb:i:1`, then passed through `sweepGA --num-mappings 1:1 --scaffold-jump 0`. The selected target-side coordinates in `selected_segments.tsv` were recovered from the already-existing strict PAF intermediates at:

- `/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_sweepga_1to1_noscaffold/<pair>.e50000.m1000.j0.8.n4.native_nb1.sweepga_1to1_noscaffold.paf`

`patches.tsv` is used only as optional annotation for community/status/pattern labels when an exact query interval, target arm, and target haplotype match is present:

- `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv`

It is not used to choose events or to draw segment geometry.

## Native assembly coordinates, not CHM13 projection

The plotted paths are named as sample assembly windows, for example:

- `PAN027#2#chr9.paternal:135704825-136204824_chr9_qarm`
- `PAN028#1#chr9.haplotype1:134380985-134880984_chr9_qarm`
- `PAN027#2#chrX.paternal:12265-512264_chrX_parm`

These names encode the sample, haplotype, chromosome label, haplotype label, and a 500 kb source interval. The local plotted segment interval is a 0-based half-open offset inside that 500 kb window. For display, the manifest projects local offsets onto the native source window by adding the path-name start coordinate. For example, local `[446944,472441)` in `PAN027#2#chr9.paternal:135704825-136204824_chr9_qarm` becomes native `chr9:136151769-136177266` in the PAN027 paternal assembly window.

No CHM13 liftover, CHM13 BED, or CHM13 projection table is used by this audit. The chromosome names in the path strings are labels attached to the sample assembly paths/windows; they should not be interpreted as CHM13 coordinates.

The path-name end coordinate is inclusive in the source label because `end - start + 1 = 500000`. The manifest therefore reports source windows as 0-based half-open by using `start` and `end + 1`, while preserving the original path string in every row. Segment native intervals are reported as 0-based half-open coordinates.

## T2T/reference status and limitations

The available files support a window-level statement, not a whole-chromosome reference statement. The strict redraw README says the source PAFs come from native `odgi untangle` over WashU pedigree paths and the plotted query denominator is recovered from the path-name interval. The selected sequence names are 500 kb subtelomeric windows labeled by sample, haplotype, chromosome, and arm.

From these files alone, we can infer:

- The schematic should draw 500 kb extracted subtelomeric windows anchored/labeled to chromosome arms.
- Query and target coordinates in the tables are native sample assembly window coordinates parsed from those path names and strict PAF offsets.
- Target-side exact intervals are recoverable from the existing strict PAF rows for the selected conservative segments and are included in `selected_segments.tsv`.

From these files alone, we cannot infer:

- That every involved source assembly is complete T2T across the full chromosome.
- That the displayed coordinates are CHM13-projected coordinates.
- That a candidate autosomal PHR event is a validated clean crossover. The chr9q/chr3q examples remain candidate terminal PHR exchange patches with same-chromosome context and side fragments.

## Event selection outcome

The selected review-facing events are exactly the requested strict-path events:

1. `PAR1_XY_positive_control`: PAN027 paternal chrXp query in `PAN027_vs_PAN011`, with chrYp PAR1 donor blocks totaling 150873 bp.
2. `PAN027_chr9q_chr3q_PHR_candidate`: PAN027 paternal chr9q query in `PAN027_vs_PAN011`, with chr3q primary donor segments totaling 45290 bp, plus chr15q/chr16q side fragments and a tiny chr20q low-confidence tail.
3. `PAN028_chr9q_chr3q_PHR_candidate`: PAN028 maternal chr9q query in `PAN028_vs_PAN027`, with chr3q primary donor segments totaling 34172 bp and a chr15q side fragment. This is the strict chr9q event, not the earlier misleading PAN028 chr3q panel.
