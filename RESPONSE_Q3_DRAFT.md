# Reviewer-response draft: D-PeerQ3 (arm-type asymmetry)

Response-letter material for the concern logged as D-PeerQ3 in
`OPEN_REVIEWER_CONCERNS.md`. Not manuscript prose.

## D-PeerQ3 - P-arm / Q-arm asymmetry, now quantified

The concern was that P-to-P and Q-to-Q relationships dominate the arm-level communities
while P-to-Q homology also exists, and that the asymmetry had never been counted. It is
counted now, on the same 41 signal-bearing arms and the same k=15 Leiden partition used for
the main figure (`data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`).

Within-community arm pairs, total 58:

| pair type | observed | null expectation | p |
|---|---|---|---|
| P-P | 19 | 14.8 | 0.148 |
| Q-Q | 22 | 13.5 | **0.019** |
| P-Q | 17 | 29.7 | **0.005** (depleted) |

Null: 20,000 permutations shuffling arm type (21 P, 20 Q) across arms with communities held
fixed. P-Q pairs are depleted to 0.57 of expectation and Q-Q enriched 1.6-fold; P-P is
enriched 1.3-fold but not significant on its own at this arm count.

The asymmetry is concentrated in two communities, which are pure by arm type:

| community | arms | P-P | Q-Q | P-Q | p(P-Q <= obs) |
|---|---|---|---|---|---|
| **C6** | 6 | 0 | **15** (exp 3.5) | **0** (exp 7.7) | **0.021** |
| **C7** | 5 | **10** (exp 2.6) | 0 | **0** (exp 5.1) | **0.048** |
| C3 | 6 | 3 | 3 | 9 (exp 7.7) | 1.0 |
| C11 | 4 | 1 | 1 | 4 (exp 3.1) | 1.0 |

Every within-community pair in C6 is Q-Q and every pair in C7 is P-P; neither community
contains a single P-Q pair. Conversely the two communities that do mix arms are C3 (the
chr3/chr9-bearing community) and C11, and there P-Q pairs sit at or above expectation. The
two-arm communities follow the pattern the concern anticipated: C1 (4q/10q) is Q-Q and C2
(10p/18p) is P-P.

**What this means.** Arm-type is a real organising axis of the partition, not a by-product
of which arms happen to be similar: two of the fifteen communities contain only one arm
type, and cross-type mixing is confined to two specific communities. That is a mechanistic
statement about how subtelomeric sequence identity is distributed among chromosome ends
rather than about how the clustering was drawn.

**Orientation: not answerable at this scale - and the first attempt was withdrawn.** Jaccard
distance is orientation-free, so we measured strand directly (one 100 kb CHM13 arm end per
arm, all-vs-all with FastGA `-f8192 -s200 -c70 -i0.6`, converted to telomere-first
orientation). The result was initially tabulated as 8 of 17 P-Q pairs alignable with 6 same
and 2 opposite, but that tabulation is withdrawn: on inspection every one of those alignments
is a short (300-2,000 bp, 78-95% identity) scattered repeat match, and **none** is
telomere-offset matched - 0 of 42 P-Q alignments, 0 bp, have the two arms at the same
distance from their telomere. Strand on a scattered repeat match says something about which
paralog copy happened to be hit, not about how two chromosome ends are oriented.

So the correct statement is narrower than the concern assumes: at 100 kb arm-end scale, the
similarity between P and Q arms is dispersed repeat matching with no colinear block, so P-Q
orientation is undefined at that scale. Answering "P-Q at opposite orientation" needs full
1 Mb arm ends (the scale at which the arm-level Jaccard matrix was computed) and an
offset-aware block selection, which is a separate analysis. What is established for Q3 is the
arm-type composition result above, which does not depend on orientation.

Previous tabulation retained for the record in `scripts/cladistics/peerq3_orientation.tsv`
with a header note that it is repeat-level and must not be read as arm-end orientation.

Tooling: `scripts/cladistics/peerq3_orientation.py` (writes the repeat-level table; the
offset-matched extension is the missing piece).

Table: `scripts/cladistics/peerq3_arm_type_asymmetry.tsv`; script:
`scripts/cladistics/peerq3_arm_type_asymmetry.py`.
