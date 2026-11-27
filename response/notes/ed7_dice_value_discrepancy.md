# Open discrepancy: the released chr9q haplotype scores 0.9918 or 0.994 in the tract

The manuscript states, in three places (Results, the Fig. 5 discussion of the released assembly, and the
Extended Data Fig. 7 caption), that against the daughter's sequence the **released** chr9q haplotype
reaches IMPG dice **0.9918** in the called tract.

Recomputing the same quantity from the competitive similarity run
(`/tmp/flagger_v11/mgall/all_sim_poa.tsv`: all candidate ends in one graph, sweepga many-to-many, FastGA
frequency 5000, `impg similarity --gfa-engine poa`, 2 kb bins; SELF.p window 136,009,149-136,209,149,
tract = window offsets 158,852-178,852 = chr9 136,168,000-136,188,000):

| statistic | value |
|---|---|
| median dice over the 110 in-tract bins | **0.9943** |
| mean dice (also bp-weighted) | 0.9935 |
| median with top and bottom 10% trimmed | 0.9945 |
| median block identity (gap-compressed, same run) | 0.9940 |
| 20 kb bins instead of 2 kb | 0.9943 (single bin) |

None reproduces 0.9918. The difference is 0.25 percentage points and changes no conclusion - the released
haplotype remains clearly below the corrected one (1.0000) and below the father's chr3q (0.9990), which is
the whole point of the panel - but two values for one quantity must not circulate.

**To resolve:** the setup behind 0.9918 (tract definition, bin size, candidate set, and whether it is dice
or block identity). If it came from a re-run after her correction pipeline, that run should be cited as the
provenance; if it is a transcription of a different statistic, the three manuscript mentions should move to
the recomputed 0.994.

**Until then** the manuscript is internally consistent at 0.9918 (three mentions, no competing value in the
text), so this is a provenance question, not an inconsistency in the paper.
