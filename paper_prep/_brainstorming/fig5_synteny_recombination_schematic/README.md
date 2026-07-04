# Fig5 Synteny Recombination Schematic Prototype

This directory contains a manifest-backed prototype for a manuscript-facing Fig5
synteny/recombination schematic. It does not modify `submission/` or replace the
manuscript Fig5.

## Files

- `event_manifest.tsv` - one row per selected review-facing event, including event class, transmission, native source windows, involved haplotypes, primary donor arms, side fragments, and recommended schematic tracks.
- `selected_segments.tsv` - one row per strict conservative segment from the selected query windows. It includes local query offsets, projected native query coordinates, target source windows, recovered target-side intervals from the strict PAF, identity/Jaccard, optional community annotations, and drawing role.
- `coordinate_provenance.md` - coordinate-system and T2T/window audit. Read this before drawing labels.
- `../fig5_washu_recombination_comparison.md` - citation/status note for the
  comparison with the WashU pedigree recombination annotation and the current
  boundary between direct preprint citation and personal communication.
- `plot_synteny_recombination_schematic.py` - standard-library Python SVG renderer.
- `fig5_synteny_recombination_full.svg` - full schematic prototype: three event rows drawn as source/product/source native 0-500 kb terminal-window views. PAR1 shows chrX source, child chrX product, and chrY source; autosomal candidates show chr9q context, child chr9q product, and chr3q donor. chr15q/chr16q/chr20q side mappings are caveat markers only.
- `fig5_synteny_recombination_full.pdf` - PDF rendering of the full source/product/source terminal-window prototype, generated with Guix `librsvg` / `rsvg-convert 2.54.5`.
- `fig5_synteny_recombination_focus.svg` - event-window prototype: the same three event rows zoomed to the native 500 kb assembly windows with a consistent physical scale and 100 kb scale bars.
- `fig5_synteny_recombination_focus.pdf` - PDF rendering of the event-window prototype, generated with Guix `librsvg` / `rsvg-convert 2.54.5`.
- `pdf_conversion_status.txt` - records whether local SVG-to-PDF conversion was available during regeneration.

For these static SVG schematics, the preferred converter is Guix `librsvg`
(`rsvg-convert`). Inkscape is a heavier fallback if a future SVG needs features
outside librsvg's rendering support.

## Inspect / Regenerate

Open the SVGs directly in a browser or vector editor:

```bash
xdg-open paper_prep/_brainstorming/fig5_synteny_recombination_schematic/fig5_synteny_recombination_full.svg
xdg-open paper_prep/_brainstorming/fig5_synteny_recombination_schematic/fig5_synteny_recombination_focus.svg
```

Regenerate from the repository root:

```bash
python3 paper_prep/_brainstorming/fig5_synteny_recombination_schematic/plot_synteny_recombination_schematic.py
```

Convert the SVGs to PDFs with Guix:

```bash
guix shell librsvg -- bash -lc '
  set -euo pipefail
  cd /moosefs/erikg/phrs
  RSVG="$(command -v rsvg-convert || find "$GUIX_ENVIRONMENT/bin" -name rsvg-convert -print -quit)"
  "$RSVG" -f pdf -o paper_prep/_brainstorming/fig5_synteny_recombination_schematic/fig5_synteny_recombination_full.pdf paper_prep/_brainstorming/fig5_synteny_recombination_schematic/fig5_synteny_recombination_full.svg
  "$RSVG" -f pdf -o paper_prep/_brainstorming/fig5_synteny_recombination_schematic/fig5_synteny_recombination_focus.pdf paper_prep/_brainstorming/fig5_synteny_recombination_schematic/fig5_synteny_recombination_focus.svg
'
```

The script reads only:

- `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/event_manifest.tsv`
- `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/selected_segments.tsv`

No heavy alignment, PAF filtering, or event reselection is performed.

## Coordinate Source and Unbanded Tracks

Segment geometry and all displayed interval labels come from native sample
assembly coordinates in `selected_segments.tsv` and `event_manifest.tsv`.
Coordinates are 0-based half-open native assembly windows parsed from path names
and strict primary-path PAF provenance. They are not CHM13-projected
coordinates.

The full schematic deliberately does not draw whole-chromosome-length tracks.
Each evidence row is a native terminal assembly window with a local 0-500 kb
axis. This keeps terminal donor patches visible at their actual local offsets
instead of shrinking them onto chromosome-scale ideograms. The terminal-window
tracks are deliberately plain and unbanded because these are extracted native
assembly windows, not cytogenetic ideograms or G-banded chromosomes.

## Event Set

The manifest intentionally selects exactly three strict primary-path events:

1. `PAR1_XY_positive_control` - `PAN027_vs_PAN011`, child/query `PAN027#2#chrX.paternal:12265-512264_chrX_parm`, chrYp PAR1 positive-control blocks.
2. `PAN027_chr9q_chr3q_PHR_candidate` - `PAN027_vs_PAN011`, child/query `PAN027#2#chr9.paternal:135704825-136204824_chr9_qarm`, terminal chr3q primary donor segments with chr15q/chr16q side fragments and a tiny chr20q low-confidence tail.
3. `PAN028_chr9q_chr3q_PHR_candidate` - `PAN028_vs_PAN027`, child/query `PAN028#1#chr9.haplotype1:134380985-134880984_chr9_qarm`, the preferred strict chr9q event with chr3q donor segments and a chr15q side fragment.

The third event deliberately replaces the earlier PAN028 chr3q review panel. The strict PAN028 chr9q path is the review-facing candidate requested here; the older chr3q panel should not be used as a substitute for this schematic.

## Drawing Guidance

Use `selected_segments.tsv` for geometry. Draw the child/recombinant haplotype as the middle destination track, with one source/donor window above and one source/donor window below. For PAR1, bracket the child chrX window with chrX and chrY PAR1 source windows. For autosomal candidates, bracket the child chr9q product with chr9q same-chromosome context above and chr3q primary donor below. Prefer block/flow/spline encodings over a many-color interval stack.

Use the `event_role` column to group segments:

- `same-chromosome context` - homologous/source context, useful as a parental chromosome track but not the recombination highlight.
- `PAR positive control` - chrYp blocks in the male PAR1 sanity-check event.
- `primary donor` - main chr3q donor blocks for autosomal candidate examples.
- `side fragment` - secondary non-homologous blocks such as chr15q/chr16q; show as caveat markers/labels on the product, not as equivalent parent/source tracks.
- `low-confidence tail` - tiny/low-score fragment such as chr20q; show as a caveat marker if retained.

Coordinates are native sample assembly window coordinates, not CHM13. Segment intervals are 0-based half-open. Target-side exact intervals are recovered from the existing strict PAF rows where possible; do not invent CHM13 coordinates or whole-chromosome T2T context from this manifest.

## Inputs Not To Use For Geometry

Do not select or draw geometry from permissive multimap/nth-best rows. `patches.tsv` may provide community/status labels only when joined exactly to a strict segment.

## Interpretation Boundary

The PAR1 row is a positive control for the drawing and strict-path handoff. The
two autosomal rows are PHR candidates compatible with terminal exchange and
include same-chromosome context plus side fragments where present. These
schematics do not claim event-level validation or a clean crossover.

For the chr9q/chr3q candidate, the manuscript may cite the public WashU
pedigree preprint for assembly provenance, but the statement that the relevant
PAN011 to PAN027 paternal lineage has no ordinary chr3 crossover remains based
on M. Cechova and T. Marschall personal communication until their recombination
annotation is public. See `../fig5_washu_recombination_comparison.md`.
