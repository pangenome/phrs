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

Regenerate from the repository root with:

```bash
python3 scripts/pedigree/plot_fig5_sweepga_1to1_redraw.py
```
