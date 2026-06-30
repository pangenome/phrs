# Community/Jaccard check for donor-recipient ribbon draft

Question checked: are the PAN027 paternal recipient/donor exchanges shown in the
ribbon draft expected from the population PHR community structure, or are they
off-community donor choices?

Inputs checked:

- `data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
- `paper_prep/manuscript_revision/C0_continuum/arm_pair_similarity_long.tsv`
- `paper_prep/manuscript_revision/C0_continuum/community_similarity_summary.tsv`

Result:

| Ribbon exchange | Arm community | Arm-pair similarity | Community context |
| --- | --- | ---: | --- |
| `chrX_p -> chrY_p` | `C15/C15` | 0.7114 | Expected PAR link; the only/peak pair in C15. |
| `chr5_q -> chr1_p` | `C11/C11` | 0.4606 | Strong expected C11 link; second-ranked C11 pair after `5q-6q`. |
| `chr9_q -> chr3_q` | `C3/C3` | 0.4885 | Strong expected C3 link; second-ranked C3 pair after `19p-3q`. |

Interpretation:

The primary exchanges highlighted in the ribbon draft sit within established
arm-level PHR communities and correspond to strong arm-pair Jaccard/similarity
links. This supports interpreting them as exchanges among known high-sharing
PHR relationships rather than random off-community donor matches.

Secondary observations:

- The small `chr9_q -> chr7_p` signal is also same-community (`C3/C3`,
  similarity 0.4028), so it is biologically plausible but visually secondary.
- The tiny `chr5_q -> chr7_p` singleton is cross-community (`C11/C3`,
  similarity 0.2794) and is not part of the clean highlighted exchange.
- The second `chr1` donor haplotype in the `chr5_q` panel is not a different
  community relationship; it is still the same `5q/1p` C11 relationship and is
  redundant for the draft display.

Caveat:

This is an arm/community-level check against the PHR Jaccard analysis. It does
not by itself prove exact haplotype-level tract identity; that remains supported
by the IMPG/SweepGA-derived winning-window evidence used to draw the ribbons.
