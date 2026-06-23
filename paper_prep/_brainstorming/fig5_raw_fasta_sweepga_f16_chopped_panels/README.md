# Fig5 raw FASTA SweepGA/FastGA f16 chopped candidate panels

Superseded: use
`paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_chop_filter_sensitivity_panels/`
for the corrected whole-genome chopped-then-filtered comparison panels.

This directory contains the earlier compact Fig5 evidence panel generated
from raw FASTA-derived whole-genome SweepGA/FastGA `--fastga-frequency 16`
PAFs, not from `odgi untangle` output.

Default evidence layer:

1. raw whole-genome f16 many:many PAFs from
   `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16`
2. Slurm-side extraction of the configured evidence windows from those
   whole-genome raw PAFs
3. validated `pafchop-rs` query-axis chopping at 2 kb, overlap 0
4. `sweepga --num-mappings 1:1 --scaffold-jump 0 --scoring ani` filtering
5. compact extraction of the PAR1 control plus the PAN027/PAN028 chr9q to chr3q
   candidate windows

Heavy PAF chopping, filtering, scanning, and figure rendering are run through
Slurm by:

```bash
paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_chopped_panels/scripts/submit_raw_fasta_chopped_panel.sh
```

Tracked deliverables:

- `fig5_raw_fasta_sweepga_f16_chopped_panels.pdf`
- `fig5_raw_fasta_sweepga_f16_chopped_panels.svg`
- `raw_fasta_chopped_panel_segments.tsv`
- `raw_fasta_chopped_panel_summary.tsv`
- `slurm_jobs.tsv`

Important caveat: this older draft performs the final 1:1 ANI filtering after
query-window extraction, so its final filtering competition is not genome-wide.
The superseding directory starts from whole-genome chopped PAFs, filters
genome-wide, and subsets windows only for display.

Large intermediate raw-window, 2 kb chopped, and 1:1 ANI-filtered PAFs are
intentionally ignored under `work/` and `evidence_paf/`. Their absolute paths
and SHA-256 checksums are recorded in `slurm_jobs.tsv`, along with the raw f16
Slurm job IDs and the panel-generation Slurm job ID.
