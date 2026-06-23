# Findings

Final SweepGA filtering was run on Slurm job `1705985` using
`/home/erikg/.cargo/bin/sweepga` version `0.1.1`
(`a0d7ac0c3312080d67de96d85cdcad9ce0c5a7e523897109b7f598c186ab85a6`).
The matrix contains 128 filtered cells: two event-bearing source PAFs times
four scaffold jumps, four minimum-alignment-length settings, two scoring modes,
and two mapping modes. The raw `many:many` f16 PAF for each comparison is also
summarized as an unfiltered baseline.

The readout in `candidate_window_summary.tsv` uses absolute query chromosome
coordinates for all three windows. Each event has 65 rows: one raw baseline plus
64 filtered cells.

## Candidate-window summary

| event | raw chr3 union bp | filtered cells with chr3 | interpretation |
|---|---:|---:|---|
| PAR1 positive control | 2,111 | 3/64 | Expected X/Y PAR1 support is present in all cells; small chr3 overlap is off-target noise. |
| PAN027 chr9q->chr3q | 261,767 | 32/64 | Raw multiway support is large; final filtering preserves chr3 in half the cells, usually at 131-180 kb union bp when present. |
| PAN028 chr9q->chr3q | 261,731 | 48/64 | Raw multiway support is large; final filtering preserves chr3 in three quarters of cells, with 24-140 kb union bp depending on scoring/mapping mode. |

Across the scaffold-jump settings, `10k`, `20k`, and `50k` behave similarly for
these candidate windows at fixed scoring/mapping/min-length. The strongest
ambiguation is driven by scoring and mapping mode: `log-length-ani` plus
`1:1` often redirects the candidate windows to chr9/self-nearby support, while
`ani` and/or `4:many` preserve more chr3 rows. Minimum alignment length changes
the chr3 union size but does not by itself erase the chr3 signal in the cells
where chr3 survives.

`--scaffold-mass` was held at the SweepGA default `10k`; `--overlap` was left
at the SweepGA default. No scaffold-mass sensitivity was added because the
candidate-window chr3 survival pattern did not show a discrete change around
the tested scaffold-jump values.
