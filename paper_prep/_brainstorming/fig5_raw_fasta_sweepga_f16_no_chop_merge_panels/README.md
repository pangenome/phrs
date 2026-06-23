# Fig5 raw FASTA f16 no-chop merge panels

This directory tests the alternative requested after the chop/filter
sensitivity panel: do not split the PAF into fixed-size chunks before the final
SweepGA filtering step. Instead, start from raw whole-genome f16 many:many PAFs,
apply SweepGA 1:1 filtering directly to those raw alignment records, and compare
no-merge vs 50 kb scaffold-merging modes.

Inputs:

- Raw whole-genome f16 many:many PAFs from
  `/moosefs/erikg/phrs/.wg-worktrees/agent-2649/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/raw_paf/`
- Display windows in `config/panel_windows.tsv`.

Filter modes:

- raw no chop; no merge; ANI:
  `sweepga --num-mappings 1:1 --scaffold-jump 0 --scoring ani`
- raw no chop; 50 kb merge; ANI:
  `sweepga --num-mappings 1:1 --scaffold-jump 50k --scoring ani`
- raw no chop; 50 kb merge; log-length ANI:
  `sweepga --num-mappings 1:1 --scaffold-jump 50k --scoring log-length-ani`

The heavy filtered PAF outputs are ignored under `filtered_paf/`; their paths
and SHA-256 checksums are recorded in `raw_merge_panel_manifest.tsv`. The plot
subsets to the configured windows only after whole-genome filtering and uses
absolute query chromosome coordinates.

Observed result:

- The raw many:many f16 layer contains chr3 support for both candidate windows
  (`raw_frequency16_chr3_chr9_support.tsv`): about 262 kb query-union chr3
  support in PAN027 and PAN028, alongside full-window chr9 support.
- Direct no-chop 1:1 filtering loses PAN027 chr3 under all tested no-chop
  modes.
- Direct no-chop 1:1 filtering keeps a single PAN028 chr3 block of 24.3 kb
  with ANI scoring, but loses chr3 under the regular 50 kb merge +
  log-length-ANI scoring mode.
- This differs from the 2 kb chopped genome-wide filter panel, where chr3
  survives in both candidates.
