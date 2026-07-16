# V6 physical-copy ontology run protocol

Status: frozen before target-term statistics were computed.

## Estimand and frozen hypotheses

The estimand for hypothesis `(collection, relation, term_id)` is the sum of
`physical_copy_cn` over coordinate-distinct CHM13 copies carrying that exact
edge in the independently gated V6 map.  `collection` is one of `GO_BP`,
`GO_MF`, `GO_CC`, or `Reactome`; `relation` is `direct` or `ancestor`.
Direct and ancestor rows are separate hypotheses even when their `term_id` is
the same.  The four multiplicity families combine both declared relation
layers within collection.  No name, family, source identifier, or term label
is substituted for physical-copy CN.

`FROZEN_HYPOTHESES.tsv.gz` is built only from
`GENOMEWIDE_SOURCE_ASSIGNMENTS.tsv.gz`, `SOURCE_TERMS.tsv.gz`, and the
target-blind CHM13 cytoband arm boundaries.  The freeze step refuses columns
whose names begin with `phr_` or otherwise encode target membership.  Terms
are not filtered by PHR burden, PHR presence, p-value, or biological name.
The PHR BED and `PHYSICAL_COPY_TERM_EDGES.tsv.gz` may be opened only after the
catalog and its checksum manifest have been written.

## Spatial null and assignments

The randomization unit is the complete rigid PHR block.  Overlapping or
abutting intervals on an arm are one block; widths, component order, internal
gaps, arm, and terminal-distance stratum never change.  Blocks are translated
uniformly over all valid integer starts on the same arm.  The prespecified
terminal strata are `[0,0.5 Mb)`, `[0.5,1 Mb)`, `[1,2 Mb)`, `[2,5 Mb)`, and
`[5 Mb, arm end)`.  The mask is empty apart from arm bounds.  All annotations,
copy clusters, tandem organization, source assignments, term edges, and CN
weights remain fixed at their genome coordinates.  The same joint placement
is used for every hypothesis.  Midpoint assignment is primary; any positive
base overlap on the identical placements is the boundary sensitivity.

## Multiplicity and empirical uncertainty

All tests are one-sided enrichment tests with ties counted as exceedances and
plus-one empirical p-values.  BH and BY are calculated separately within the
four collection families.  Single-step Westfall--Young maxT is calculated
within each collection and across all four collections as a global safeguard,
using pooled standardized observed/null statistics.  Every empirical p-value
has a two-sided 95% Clopper--Pearson interval.  Sequential decisions use the
more conservative two-sided 98.75% interval (alpha 0.0125 per look).

## Staged Monte Carlo stopping rule

The sufficient initial screen is 99,999 common, valid joint placements for all
hypotheses and both assignments.  A primary hypothesis is unresolved when a
98.75% BH or global-maxT interval straddles 0.05.  Low raw exceedance counts and
point estimates in `[0.04,0.06]` are reported as precision diagnostics but do
not by themselves override a resolved multiplicity decision.  A positive
observed count with an all-zero initial null is structurally unresolved and
must be audited through the cap.

Selective extension is permitted only for hypotheses whose collection and
global maxT decisions are already resolved at the complete initial screen.
Their raw/BH counts may be extended on the continuing PCG64DXSM stream at
249,999, 599,999, and 999,999 total placements.  Resolved candidates stop at
the first checkpoint.  MaxT values remain explicitly tied to the complete
99,999-hypothesis screen; a maxT-unresolved hypothesis is reported unresolved
rather than being used to trigger a reflexive million-permutation rerun of
every term.  Candidate lists, reasons, checkpoints, and stopping outcomes are
immutable artifacts once emitted.

Primary evidence requires both BH upper confidence bound at or below 0.05 and
global maxT upper confidence bound at or below 0.05.  Families and source
genes enter only contributor metadata after inference and cannot create a
hypothesis or discovery.

## Reproducibility constants

- Primary master seed: `2026071301`.
- NumPy generator: `PCG64DXSM`, child spawn key `[0]`.
- Initial permutations: `99999`.
- Extension checkpoints: `249999`, `599999`, `999999`.
- Initial count array type: unsigned 16-bit after a hard check that every
  possible physical-CN burden is at most 65,535.
- Required independent gate: `PRE_RUN_V6_GATE.json` status `PASS`,
  `enrichment_authorized=true`, and exact agreement with
  `PRE_RUN_V6_AUDITED_RELEASE.sha256.tsv` and `INPUT_MANIFEST.tsv`.

The manuscript is not an input and is not modified by this run.
