# Fig5 joint-parent direct sweepGA schematics

Sibling schematic comparison for the corrected joint-parent direct sweepGA
experiment in `../pedigree_direct_sweepga_joint_parent/`.

The segment tables use the same column schema as
`../fig5_synteny_recombination_sweepga_1to1/selected_segments.sweepga_1to1.tsv`
so the outputs can be compared directly with the original Fig5 schematic and
the earlier separate-haplotype direct-sweepGA sibling.

Coordinates are local 0-based half-open offsets inside each 500 kb telomeric
source window from:

`/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/washu.1Mb.telo_500kb_trimmed.fa.gz`

They are not whole-genome alignment coordinates and are not CHM13 projections.

Regenerate from the repository root after the joint PAF package has been built:

```bash
python3 paper_prep/_brainstorming/fig5_synteny_recombination_joint_parent/build_selected_segments_from_joint_parent_paf.py --filter-id one_one_noscaffold
python3 paper_prep/_brainstorming/fig5_synteny_recombination_joint_parent/build_selected_segments_from_joint_parent_paf.py --filter-id four_many_noscaffold
python3 paper_prep/_brainstorming/fig5_synteny_recombination_joint_parent/build_selected_segments_from_joint_parent_paf.py --filter-id many_many_noscaffold
python3 paper_prep/_brainstorming/fig5_synteny_recombination_joint_parent/plot_synteny_recombination_joint_parent.py --filter-id one_one_noscaffold
python3 paper_prep/_brainstorming/fig5_synteny_recombination_joint_parent/plot_synteny_recombination_joint_parent.py --filter-id four_many_noscaffold
```

Rendered required siblings:

- `fig5_synteny_recombination_joint_parent_1to1_full.pdf`
- `fig5_synteny_recombination_joint_parent_4many_full.pdf`

Interpretation: joint `1:1` is a diagnostic control and removes all rendered
chr3q donor rows for these Fig5 events. Joint `4:many` is the preferred visible
comparison because it matches raw many:many for the selected events while still
using a bounded mapping filter; it preserves PAR1 and recovers the PAN028 chr3q
donor segment, but the direct joint PAF does not recover the PAN027 chr3q donor
seen in the original untangle-derived schematic.
