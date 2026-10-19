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

**Orientation: measured at the scale that answers it, and the answer is 'same, almost
always'.** The right data are the HPRCv2 all-vs-all alignments of the telomeric 500 kb flanks
(`all-vs-all.1Mb.p95`, one PAF per CHM13 flank, 46 arm ends). Orientation is scored per
alignment block - no colinearity or offset matching required, because the question is about
the modules, not about whole-arm colinearity - with flanks re-oriented telomere-first
(p-arm telomere at low coordinate, q-arm telomere at high, so strand '+' means same
orientation for P-P and Q-Q pairs and strand '-' for P-Q pairs). Homolog-vs-homolog pairs are
excluded, since they are trivially same-orientation.

Blocks of >= 5 kb, pairs with >= 20 kb aligned, 764 arm pairs, 9.6 Gbp:

| pair type | arm pairs | same telomere-relative orientation |
|---|---|---|
| all | 764 | **98.6%** |
| P-P | 152 | 97.8% |
| Q-Q | 228 | **100.0%** |
| P-Q | 384 | 98.4% |

The per-pair distribution is not an average over a mixed bag: the **median pair is 100%**
same-orientation, the 10th percentile is also 100%, and only 7 of 764 pairs fall below 90%.
So the shared modules are reassorted in ORDER but are inserted in a consistent
telomere-to-centromere orientation wherever they land, including between p and q arms, whose
coordinate systems run in opposite directions.

The exceptions are informative rather than noise - they are the known inverted-repeat
regions: chr8p versus chr1p (0.0% same, 31.8 Mb) and versus chr2q (5.7%), with chr4q/chr14p
and chr20p/chr10q the only other pairs below 90%. The classical direct-repeat pairs are
100% same-orientation: chr4q/chr10q, chr10p/chr18p, chr15p/chr21p, and chr1p/chr5q.

**Mechanistic reading.** A direct (same-orientation) bias is what repeat-mediated exchange
predicts: gene conversion and NAHR are strongly biased toward direct repeats, while exchange
between inverted repeats yields dicentric products and is selected against. Order-scrambled
but orientation-consistent modules are therefore the expected signature of repeated
duplication plus conversion in these regions, not of arbitrary reshuffling.

An earlier pass on 100 kb single windows produced a differently-signed answer; that pass was
unreliable (all P-Q alignments there were 300-2,000 bp scattered repeat matches, none
telomere-offset matched at the arm-flank scale) and is retained only as a cautionary note in
`peerq3_orientation.tsv`.

Scripts: `scripts/cladistics/arm_module_orientation.py` (the result above, writes
`arm_module_orientation.tsv`) and `scripts/cladistics/peerq3_orientation.py` (the superseded
100 kb pass, retained with a warning header).

Table: `scripts/cladistics/peerq3_arm_type_asymmetry.tsv`; script:
`scripts/cladistics/peerq3_arm_type_asymmetry.py`.
