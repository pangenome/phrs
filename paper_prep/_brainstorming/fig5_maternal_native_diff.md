# Fig. 5 PAN027 Maternal Native Untangle Visual Diff

Scope: PAN027 maternal haplotype only (`PAN027#1`, inherited from PAN010),
drawn against PAN010 target haplotypes for the Fig. 5-style maternal panel.
This is a visualization/provenance comparison, not a recombination mechanism
analysis. The manuscript asset
`submission/fig/MainFigures/Fig5_pedigree_untangle.pdf` was not edited.

## Inputs and Semantics

The exact original plotting recipe for the committed Fig. 5 PDF was not found
in the repo. The PDF identifies as an R graphics output, and the project notes
state that current submission main-figure assets are slide-derived placeholders
that should eventually be regenerated from figure sources. I therefore made a
standardized same-grammar comparison from the old/source and native data.

Important orientation check: the old file named
`PAN027_vs_PAN010.e50000.m1000.bed.gz` is not child-query only. Streaming the
file shows all four WashU samples as query paths against PAN010 reference paths:
`PAN010`, `PAN011`, `PAN027`, and `PAN028` all appear in the query column, while
the reference column is PAN010. For this Fig. 5 maternal comparison, the old
source view was filtered to `PAN027#1` query paths only. The native Slurm BEDPE
and PAF outputs are already scoped to PAN027 query paths against PAN010 targets.

Source rows loaded for the comparison:

| View | Rows loaded after PAN027#1/PAN010 filtering |
|---|---:|
| old BED first-best (`nth.best=1`) | 3,426 |
| native BEDPE n1 | 3,295 |
| native BEDPE n4 | 9,810 |
| native PAF n1 | 3,295 |
| native PAF n4 | 9,810 |
| sweepGA PAF 1:many | 5,526 |
| sweepGA PAF 2:many | 8,026 |
| sweepGA PAF 4:many | 9,810 |

## Artifacts

Generated under `paper_prep/_brainstorming/fig5_maternal_native_diff/`:

- `fig5_maternal_native_diff.svg` - side-by-side visual diff. Rows are
  `PAN027#1` query arms; x-axis is each 500 kb flank; color encodes PAN010
  target arm; sublanes show nth-best ranks; red outlines mark binned intervals
  where a regenerated view lacks the old first-best target arm.
- `source_manifest.tsv` - exact old/native/sweepGA inputs and loaded row counts.
- `segments_compact.tsv` - compact coalesced segments used for the SVG.
- `arm_comparison.tsv` - per-query-arm summary aligned on 1 kb query bins.

Reproducible script:
`scripts/pedigree/plot_fig5_maternal_native_diff.py`.

Large native and sweepGA PAF/BEDPE intermediates are external only and were not
copied into git. The script reads:

- old BED:
  `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/PAN027_vs_PAN010.e50000.m1000.bed.gz`
- native BEDPE:
  `/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/PAN027_vs_PAN010.e50000.m1000.j0.8.n4.bedpe.gz`
- native PAF:
  `/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/PAN027_vs_PAN010.e50000.m1000.j0.8.n4.paf.gz`
- sweepGA PAN027/PAN010 external outputs generated cheaply from the
  uncompressed native PAF in
  `/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_sweepga_test/`.

## Decision

Does the plotted maternal untangle view change visually?

Yes, if the regenerated panel shows native n-best/multimap content. No, not in
the substantive first-best schematic sense. The native BEDPE n1 and native PAF
n1 views agree with each other and are almost identical to the old
`PAN027#1`-filtered first-best view. On a 1 kb binned comparison, native PAF n1
has only 16 kb where the old first-best target arm is absent, affecting three
query arms: `chr21p` (10 kb), `chr13p` (3 kb), and `chr18p` (3 kb). This is
visually minor relative to the 48-arm, 500 kb-per-arm Fig. 5 panel.

The visible change appears when the native n4 or sweepGA-filtered PAF views are
drawn as multimap-aware views. Those views add alternate equivalent targets
rather than replacing the old first-best target. In the native PAF n4 and
native BEDPE n4 views, the old first-best target arm remains present in every
1 kb bin where the old view has a call (`old_n1_target_absent_bp_binned = 0`),
but alternate target arms cover 4.439 Mb of multi-target bins and 21.834 Mb of
rank-2-plus-present bins across the PAN027 maternal query arms. The added
content is concentrated in known ambiguous subtelomeric families, especially
acrocentric p-arm rows:

| Query arm | Multi-target bp in native PAF n4 | Target-arm set shown |
|---|---:|---|
| `chr21p` | 498,000 | `chr13p,chr14p,chr15p,chr21p,chr22p` |
| `chr22p` | 496,000 | `chr13p,chr14p,chr15p,chr21p,chr22p` |
| `chr13p` | 465,000 | `chr13p,chr14p,chr15p,chr21p,chr22p` |
| `chr14p` | 335,000 | `chr13p,chr14p,chr15p,chr21p,chr22p` |
| `chr15p` | 310,000 | `chr13p,chr14p,chr15p,chr21p,chr22p` |

Other subtelomeric multimapping appears in rows such as `chr16p`, `chr9q`,
`chr3q`, `chr5q`, and `chr7p`, where n4 adds target arms including `chr16p`,
`chr16q`, `chr9q`, `chr11p`, `chr6p`, `chr6q`, `chr19p`, `chr5q`, `chr7p`,
and `chr3q`. These are additional equivalent or near-equivalent donor-arm
assignments in the visual grammar, not a change to the first-best story.

The sweepGA filters mostly preserve that interpretation. `sweepGA 4:many` is
identical in row count to native PAF n4 here. `sweepGA 2:many` keeps 8,026 of
9,810 native PAF rows and has only 1 kb where the old first-best target is
absent. `sweepGA 1:many` keeps 5,526 rows and visibly reduces, but does not
remove, the alternate-target overlay; it also has only 1 kb where the old
first-best target is absent.

## Figure Defensibility

The current Fig. 5 maternal panel remains defensible as a schematic/first-best
view of PAN027 maternal (`PAN027#1`) against PAN010, provided it is described
as a first-best/source rendering and not as an exhaustive accounting of
equivalent subtelomeric mappings. The manuscript asset itself is unchanged by
this task.

A later figure-replacement task is warranted if Fig. 5 is expected to carry
provenance-clean native odgi output or explicitly show ambiguous multimapping.
The replacement is not needed to rescue the current first-best visual claim;
it is needed only to make the figure source current and to decide whether the
published panel should show alternate target ribbons or keep them out of the
main schematic.

## Recommended Future Source

Future regenerated Fig. 5 should use the native odgi PAF as the source of
record:

`/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/PAN027_vs_PAN010.e50000.m1000.j0.8.n4.paf.gz`

For a first-best schematic matching the current panel, filter this PAF to
`PAN027#1` query paths, PAN010 target paths, and `nb:i:1`. For a
multimap-aware replacement, use the same native PAF with a documented
sweepGA-filtered overlay, preferably the external `2:many` result as a compact
middle ground:

`/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_sweepga_test/PAN027_vs_PAN010.e50000.m1000.j0.8.n4.uncompressed.2_many.paf`

The next task should regenerate the full Fig. 5 figure from native odgi PAF for
all intended pedigree panels, compare first-best-only versus multimap-overlay
designs, and make an author-facing choice about whether alternate target arms
belong in the main figure or in a supplemental/provenance panel. That task
should replace assets only after the visual design is chosen; this task does
not replace `submission/fig/MainFigures/Fig5_pedigree_untangle.pdf`.
