# Fig. 5 follow-up transmission check

## Question

The manuscript Fig. 5 candidate interrogates the paternal haplotype of PAN027
against PAN011. This follow-up applies the same strict-path procedure to
PAN028, the child of PAN027, and asks whether the implicated chromosome ends are
visible in the transmitted PAN028 maternal haplotype.

## Procedure

Input segments are the existing strict `nb=1` plus sweepGA `1:1` no-scaffold
rows in `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/selected_segments.tsv`.
No alignment is recomputed here. The script extracts the chr9q -> chr3q
candidate for `PAN027_vs_PAN011` and the matching `PAN028_vs_PAN027` chr9q
candidate, then summarizes the same 500 kb local query coordinate system used by
the Fig. 5 schematic prototype.

## Result

The result is clean enough for a candidate Fig. 5 update. PAN028's maternal
chr9q window retains the same implicated chromosome-end classes seen in the
PAN027 paternal haplotype: chr9q context plus chr3q primary-donor sequence and
chr15q/chr16q side fragments. In PAN028, 335,908 bp of same-chromosome
chr9q context maps directly to PAN027's paternal hap2, and the two diagnostic
side fragments also map to PAN027 paternal hap2 (1,207 bp chr16q
and 15,166 bp chr15q). The child has 34,172 bp of chr3q
primary-donor sequence overall, split between PAN027 maternal and paternal chr3q
sources, so the chr3q end is present in the transmitted haplotype but is not a
single intact paternal-hap2-only block.

The 493 bp chr20q low-confidence tail in PAN027 is not recovered in PAN028.
That is desirable for the figure update: the candidate should focus on chr9q,
chr3q, chr15q, and chr16q, and omit chr20q from the interpreted event model.

## Summary Table

| Event | Same chr9q context | chr3q primary | chr16q side | chr15q side | chr20q low-conf | Weighted identity |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| PAN027 paternal hap2 <- PAN011 | 445,737 | 45,290 | 1,207 | 7,273 | 493 | 99.9195 |
| PAN028 maternal hap1 <- PAN027 | 449,356 | 34,172 | 1,207 | 15,166 | 0 | 99.8415 |

## Candidate Figure Update

Use `fig5_followup_transmission_check.svg` as a candidate companion/update panel:
two aligned rows show the originally interrogated PAN027 paternal haplotype and
the PAN028 maternal follow-up. Black-outlined blocks are inter-chromosomal
mappings; gray-outlined blocks are same-chromosome chr9q context.

## Outputs

- `fig5_followup_transmission_check.svg`
- `transmission_event_summary.tsv`
- `pan028_source_breakdown.tsv`
- `local_interval_comparison.tsv`

## Interpretation Boundary

This is a transmission consistency check, not a new de novo-event proof. It
supports the expected transmission of the implicated chromosome-end pattern into
PAN028, with the caveat that chr3q support in PAN028 is distributed across both
PAN027 haplotypes.
