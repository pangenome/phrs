# Fig5 Length-Scaled Whole-Genome Tracks

This directory renders a more readable whole-genome view for the Fig5 pedigree
alignment evidence.

The key geometry change is that each chromosome is drawn as its own horizontal
row, with the row length proportional to the native query chromosome length.
The PDF has one page per queried child haplotype/transmission:

- `PAN027_hap1_maternal_vs_PAN010`
- `PAN027_hap2_paternal_vs_PAN011`
- `PAN028_hap1_maternal_vs_PAN027`

Each page contains three comparable evidence layers:

- strict untangle primary-path geometry, without chopping;
- updated `wfmash -p95` whole-genome alignment, 2 kb query-grid chopped, then
  SweepGA `1:1` ANI filtering;
- SweepGA/FastGA whole-genome alignment with `--fastga-frequency 32`, 2 kb
  query-grid chopped, then SweepGA `1:1` ANI filtering.

For the two filtered alignment layers, the display uses the precomputed 1 Mb
whole-genome support table generated from the 2 kb query-grid filtered PAFs.
If a display bin contains retained interchromosomal support, the bin is colored
and labeled by the strongest interchromosomal target chromosome; otherwise it is
shown as same-chromosome gray or no-support background. This intentionally
prevents same-chromosome support from visually hiding weaker interchromosomal
target support in the whole-genome view.

Run from this directory:

```bash
bash run_length_scaled_tracks.sh
```

Primary output:

- `fig5_whole_genome_length_scaled_tracks.pdf`

Per-transmission outputs:

- `fig5_whole_genome_length_scaled_tracks.PAN027_hap1_maternal_vs_PAN010.{pdf,png,svg}`
- `fig5_whole_genome_length_scaled_tracks.PAN027_hap2_paternal_vs_PAN011.{pdf,png,svg}`
- `fig5_whole_genome_length_scaled_tracks.PAN028_hap1_maternal_vs_PAN027.{pdf,png,svg}`

Tables:

- `length_scaled_track_segments.tsv`
- `length_scaled_track_chromosomes.tsv`
- `length_scaled_track_summary.tsv`
- `length_scaled_track_manifest.tsv`
