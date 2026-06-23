# Fig5 raw FASTA f16 chop/filter sensitivity panels

This directory contains corrected Fig5 candidate panels from raw FASTA
SweepGA/FastGA `-f16` whole-genome PAFs that were chopped before filtering.

The filtering order is:

1. Start from whole-genome raw f16 PAFs in
   `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16`.
2. Use whole-genome chopped PAFs at 2 kb, 5 kb, and 10 kb.
3. Apply sweepGA filtering to the whole-genome chopped PAFs.
4. Subset the two Fig5 chr9q-to-chr3q candidate windows only for plotting and
   summary extraction.

The committed panel uses absolute query chromosome coordinates rather than
local 0-500 kb window coordinates. Heavy filtered PAF intermediates are ignored
under `filtered_paf/`; the committed `filter_manifest.tsv` records their paths,
hashes, scoring mode, and scaffold setting.

## Filter modes

- `no_merge_ani`: `sweepga --num-mappings 1:1 --scaffold-jump 0 --scoring ani`
- `no_merge_log_length_ani`: `sweepga --num-mappings 1:1 --scaffold-jump 0 --scoring log-length-ani`
- `scaffold50k_ani`: `sweepga --num-mappings 1:1 --scaffold-jump 50k --scoring ani`
- `scaffold50k_log_length_ani`: `sweepga --num-mappings 1:1 --scaffold-jump 50k --scoring log-length-ani`

## Chr3 survival summary

Both Fig5 chr9q-to-chr3q candidate windows retain chr3 support in every tested
combination of chop length and filter mode:

- PAN027 paternal chr9q candidate: chr3 survives at 2 kb, 5 kb, and 10 kb in
  all four modes. Union-covered chr3 support ranges from 10,000 bp to 66,277 bp.
- PAN028 maternal chr9q candidate: chr3 survives at 2 kb, 5 kb, and 10 kb in
  all four modes. Union-covered chr3 support ranges from 10,000 bp to 43,121 bp.

Primary outputs:

- `fig5_raw_fasta_sweepga_f16_chop_filter_sensitivity_panels.pdf`
- `fig5_raw_fasta_sweepga_f16_chop_filter_sensitivity_panels.svg`
- `chop_filter_panel_summary.tsv`
- `chop_filter_panel_segments.tsv`
- `chop_filter_target_breakdown.tsv`
- `filter_manifest.tsv`
