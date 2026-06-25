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

For the two filtered alignment layers, the display now scans the filtered PAFs
directly and plots exact retained interchromosomal query intervals. It no longer
uses the earlier 1 Mb binned support table for PAF geometry. Same-chromosome PAF
records are not plotted one-by-one; the gray chromosome row is only a coordinate
baseline. This keeps tiny 2 kb retained chr3 hits from being visually expanded
into 1 Mb blocks.

Each page also includes a candidate-window zoom band. The whole-genome rows keep
chromosome-scale geometry, while the zoom band shows the configured PAR/chr9q
candidate windows on a local scale so exact 2 kb PAF records are visible.

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

In `length_scaled_track_segments.tsv`, PAF-derived rows have
`source_row_type=filtered_paf_record_exact_query_interval`, and
`segment_start`/`segment_end` are the original PAF query start/end coordinates.
