# Fig5 Synteny Recombination Schematic Manifest

This directory is a data-audit handoff for a future Fig5 recombination/synteny schematic. It does not contain a final figure and does not modify `submission/`.

## Files

- `event_manifest.tsv` - one row per selected review-facing event, including event class, transmission, native source windows, involved haplotypes, primary donor arms, side fragments, and recommended schematic tracks.
- `selected_segments.tsv` - one row per strict conservative segment from the selected query windows. It includes local query offsets, projected native query coordinates, target source windows, recovered target-side intervals from the strict PAF, identity/Jaccard, optional community annotations, and drawing role.
- `coordinate_provenance.md` - coordinate-system and T2T/window audit. Read this before drawing labels.

## Event Set

The manifest intentionally selects exactly three strict primary-path events:

1. `PAR1_XY_positive_control` - `PAN027_vs_PAN011`, child/query `PAN027#2#chrX.paternal:12265-512264_chrX_parm`, chrYp PAR1 positive-control blocks.
2. `PAN027_chr9q_chr3q_PHR_candidate` - `PAN027_vs_PAN011`, child/query `PAN027#2#chr9.paternal:135704825-136204824_chr9_qarm`, terminal chr3q primary donor segments with chr15q/chr16q side fragments and a tiny chr20q low-confidence tail.
3. `PAN028_chr9q_chr3q_PHR_candidate` - `PAN028_vs_PAN027`, child/query `PAN028#1#chr9.haplotype1:134380985-134880984_chr9_qarm`, the preferred strict chr9q event with chr3q donor segments and a chr15q side fragment.

The third event deliberately replaces the earlier PAN028 chr3q review panel. The strict PAN028 chr9q path is the review-facing candidate requested here; the older chr3q panel should not be used as a substitute for this schematic.

## Drawing Guidance

Use `selected_segments.tsv` for geometry. Draw the child/recombinant haplotype as the destination track, then draw source tracks for same-chromosome context and non-homologous donor blocks. Prefer block/flow/spline encodings over a many-color interval stack.

Use the `event_role` column to group segments:

- `same-chromosome context` - homologous/source context, useful as a parental chromosome track but not the recombination highlight.
- `PAR positive control` - chrYp blocks in the male PAR1 sanity-check event.
- `primary donor` - main chr3q donor blocks for autosomal candidate examples.
- `side fragment` - secondary non-homologous blocks that should be shown explicitly when they require another source chromosome.
- `low-confidence tail` - tiny/low-score fragment; show as a caveat if retained.

Coordinates are native sample assembly window coordinates, not CHM13. Segment intervals are 0-based half-open. Target-side exact intervals are recovered from the existing strict PAF rows where possible; do not invent CHM13 coordinates or whole-chromosome T2T context from this manifest.

## Inputs Not To Use For Geometry

Do not select or draw geometry from permissive multimap/nth-best rows. `patches.tsv` may provide community/status labels only when joined exactly to a strict segment.
