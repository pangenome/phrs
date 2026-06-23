# Fig5 raw FASTA SweepGA f16 multiway panels

This package renders the Fig5 inspection view for the whole-genome raw FASTA
SweepGA/FastGA frequency-16 evidence.  The source-of-truth multiway input is the
completed whole-genome raw many:many PAF layer from:

`/moosefs/erikg/phrs/.wg-worktrees/agent-2649/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/raw_paf/*.sweepga_frequency16_many_many_j0.paf.gz`

The existing chopped many:many, four:many and one:one PAFs from that same
package are included only as comparison layers.  They are not used to define the
raw multiway support.

## Outputs

- `fig5_raw_fasta_sweepga_f16_multiway_panels.pdf`, `.svg` and `.png`: full
  multiway panel, with raw many:many support shown above the comparison layers
  for each event.
- `preview_png/fig5_raw_fasta_sweepga_f16_multiway_panels.png`: PNG preview.
- `multiway_candidate_support.tsv`: row-level support for the PAR1 control,
  PAN027 chr9q->chr3q window and PAN028 chr9q->chr3q window across raw and
  comparison layers.
- `PAR1_XY_positive_control.raw_support.tsv`,
  `PAN027_chr9q_chr3q_PHR_candidate.raw_support.tsv`,
  `PAN028_chr9q_chr3q_PHR_candidate.raw_support.tsv`: raw many:many row-level
  support tables for each displayed window.
- `multiway_candidate_summary.tsv`: per-window and per-layer row counts,
  row-multiplicity, query-union coverage and target-chromosome support.
- `raw_chr3_chr9_other_support_summary.tsv`: compact pre-1:1 summary of chr3,
  chr9 and other target support.
- `input_manifest.tsv`: source paths, checksums and raw f16 Slurm provenance.

## Regeneration

The extractor streams over gzipped PAFs and writes only compact TSV summaries;
it does not vendor heavy PAF slices.  On Slurm:

```bash
cd paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_multiway_panels
sbatch scripts/run_multiway_panel.sbatch
```

For a local validation run on a node where streaming the gzipped inputs is
acceptable:

```bash
python3 scripts/extract_multiway_support.py --panel-dir .
Rscript scripts/plot_multiway_panels.R .
bash scripts/validate_outputs.sh
```

The raw whole-genome PAFs are 2.2-5.0 GB each, so regeneration should normally
be done with the Slurm wrapper or on a compute node.
