# Fig5 raw FASTA f16 chop/filter sensitivity panels

This directory corrects the earlier `fig5_raw_fasta_sweepga_f16_chopped_panels`
draft. The earlier draft used whole-genome raw FASTA-derived SweepGA/FastGA
PAFs, but it extracted the query windows before the final SweepGA 1:1 filtering
step. This directory instead uses whole-genome chopped PAFs, applies the final
SweepGA filtering genome-wide, and only then subsets the three display windows.

Inputs:

- Source package:
  `/moosefs/erikg/phrs/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16`
- Raw layer: SweepGA/FastGA with FastGA k-mer frequency 16, many:many,
  scaffold jump 0.
- Chop lengths: 2 kb, 5 kb, 10 kb.
- Filter modes:
  - no merge; ANI (`--scaffold-jump 0 --scoring ani`)
  - no merge; log-length ANI (`--scaffold-jump 0 --scoring log-length-ani`)
  - 50 kb merge; ANI (`--scaffold-jump 50k --scoring ani`)
  - 50 kb merge; log-length ANI (`--scaffold-jump 50k --scoring log-length-ani`)

The existing no-merge ANI full-genome filtered PAFs are reused from the source
package. The other filter modes are generated into ignored `filtered_paf/` files
by the Slurm array script in `scripts/run_filter_mode_array.sbatch`.

Outputs:

- `fig5_raw_fasta_sweepga_f16_chop_filter_sensitivity_panels.pdf`
- `fig5_raw_fasta_sweepga_f16_chop_filter_sensitivity_panels.svg`
- `chop_filter_panel_segments.tsv`
- `chop_filter_panel_summary.tsv`
- `chop_filter_panel_manifest.tsv`

The plot uses absolute query chromosome coordinates on every row. The grey
background is only the 500 kb display window; it is not the filtering scope.
