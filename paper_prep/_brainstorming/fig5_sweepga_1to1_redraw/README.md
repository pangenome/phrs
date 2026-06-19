# Fig5 sweepGA 1:1 inspection redraw

This directory contains an author-facing inspection redraw of the pedigree/Fig. 5-style untangle panel.
It is exploratory output only and does not update `submission/paper.tex` or `submission/fig/MainFigures/Fig5_pedigree_untangle.pdf`.

## Inputs and filtering

For each WashU transmission, the script starts from the existing native odgi untangle m1000 n4 PAF in `/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/`.
The conservative plotted view first filters the native PAF to rows whose optional tags contain `nb:i:1`, then runs:

```bash
/moosefs/erikg/sweepga/target/release/sweepga --num-mappings 1:1 --scaffold-jump 0 --output-file <conservative.paf> <native_nb1.paf>
```

This differs from plain sweepGA 1:1 because plain sweepGA is run directly on the native n4 PAF and can retain `nb:i:2`, `nb:i:3`, or `nb:i:4` alternates when those are equivalent under its filter. The conservative redraw removes those alternates before sweepGA, so the plotted PAF contains only `nb:i:1` rows.

It also differs from native first-best alone because sweepGA applies an additional reciprocal 1:1 plane-sweep filter with scaffolding disabled (`--scaffold-jump 0`).

## Row counts

| Pair | Native n4 | Native nb1 | Plain sweepGA n4 1:1 no scaffold | Conservative nb1 -> sweepGA 1:1 no scaffold | Plain nb values | Conservative nb values |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| `PAN027_vs_PAN010` | 9810 | 3295 | 3551 | 2737 | 1,2,3,4 | 1 |
| `PAN027_vs_PAN011` | 9938 | 3331 | 3622 | 2835 | 1,2,3,4 | 1 |
| `PAN028_vs_PAN027` | 8742 | 2948 | 3164 | 2270 | 1,2,3,4 | 1 |

## Files

- `fig5_sweepga_1to1_redraw.svg` and `fig5_sweepga_1to1_redraw.pdf`: compact author-facing redraw.
- `summary_counts.tsv`: row counts, `nb` values, inter-chromosomal row counts, source paths, and exact command/filter per stage.
- `conservative_segments.tsv`: compact plotted segment table after coalescing adjacent same-target rows.
- `validation_report.tsv`: coordinate, `nb`, query-length, and SVG/PDF rectangle-bound checks from the last regeneration.

## Coordinate correction

The superseded first redraw incorrectly used PAF column 2 (`fields[1]`) as the full plotting denominator. In these odgi PAFs that field can be 41,248 even when `qstart`/`qend` span hundreds of kilobases, which created invalid coalesced rows and off-panel SVG rectangles.

The corrected redraw derives `query_length` from the child/query path name interval, for example `PAN027#1#chr3.maternal:9503-509502_chr3_parm` gives 500,000 bp. The raw PAF column 2 value is retained only as `paf_query_length` for audit. Rows with invalid child/query coordinates are dropped before coalescing, and rendering validates that all rectangles stay within their p/q-arm track bounds.

Regenerate from the repository root with:

```bash
python3 scripts/pedigree/plot_fig5_sweepga_1to1_redraw.py
```
