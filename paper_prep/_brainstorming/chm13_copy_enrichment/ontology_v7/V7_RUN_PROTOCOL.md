# V7 whole-genome non-PHR physical-copy ontology run protocol

Status: frozen before V7 term counts, p-values, adjusted values, or summaries
were computed.

## Scientific estimand and inherited frozen objects

For every V6-frozen exact hypothesis
`T = (collection, relation, term_id)`, V7 tests whether the CHM13 PHR
compartment carries excess annotation-bearing physical-copy burden relative to
the complete ontology-eligible CHM13 genome outside PHRs. This is a complete
finite-population copy-burden enrichment analysis. It is not a gene-list ORA,
placement test, randomization, permutation test, or sampling analysis.

V7 consumes, without altering or refreezing:

- `../ontology_v6/GENOMEWIDE_SOURCE_MAP.tsv.gz`, the gated coordinate-distinct
  CHM13 physical-copy/source map;
- `../ontology_v6/PHYSICAL_COPY_TERM_EDGES.tsv.gz`, the exact direct and
  ancestor GO/Reactome copy-term edges;
- `../ontology_v6/FROZEN_HYPOTHESES.tsv.gz` and
  `../ontology_v6/FROZEN_HYPOTHESES.json`, the target-blind exact hypothesis
  catalog and freeze manifest; and
- `data/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`, used only after
  inference to label PHR contributor arms in a descriptive community summary.

Each coordinate-distinct physical copy has `physical_copy_cn = 1` in the V6
release and contributes one unit to each exact term it carries. Copies are
never collapsed by gene name, functional source, source gene, gene or locus
family, duplicated segment, ontology provenance, or community. Direct and
ancestor rows remain separate hypotheses even when they share a term ID.

## Eligible universe and paired assignments

The frozen ontology-eligible universe is exactly the set of V6 physical-copy
rows with at least one V6 frozen direct-or-ancestor edge. Eligibility is
therefore target-independent and identical for midpoint and overlap analyses.
Let its complete physical-copy count be `N`.

For each assignment separately:

- primary midpoint membership is `phr_midpoint_cn = 1`;
- paired any-overlap sensitivity membership is `phr_any_overlap_cn = 1`;
- the PHR count is `K`; and
- the non-PHR background is the exact set complement, of size `N - K`.

No PHR copy may occur in the corresponding non-PHR background. Every eligible
copy must occur exactly once in one side of the partition. The entire non-PHR
set is used; it is neither sampled nor spatially matched.

## Prespecified term table and exact inference

For every frozen hypothesis `T`, V7 computes coordinate-copy counts

- `a_T`: PHR copies carrying `T`;
- `b_T = K - a_T`: PHR eligible copies not carrying `T`;
- `c_T`: non-PHR copies carrying `T`;
- `d_T = (N - K) - c_T`: non-PHR eligible copies not carrying `T`;
- `M_T = a_T + c_T`: whole-genome term burden.

The one-sided enrichment probability is the complete finite-population upper
tail

`P[X >= a_T | X ~ Hypergeometric(N, M_T, K)]`.

The hypergeometric distribution is only the exact calculation for this fixed,
fully enumerated physical-copy population; it does not imply that copies are
biological replicates or that a generic gene-list analysis was performed.

Reported effects are PHR and non-PHR copy fractions, their unsmoothed fold
enrichment with explicit `inf`/`NA` zero states, a Haldane--Anscombe-smoothed
fold and odds ratio (0.5 in all four cells), the unsmoothed odds ratio with
explicit zero states, `a_T - c_T`, and `a_T - K M_T / N`.

## Multiplicity and fixed primary support rule

For each assignment, BH and BY are applied within each V6 collection family
(`GO_BP`, `GO_MF`, `GO_CC`, `Reactome`) over all direct plus ancestor rows.
Holm and Bonferroni are additionally applied across all V6-frozen hypotheses.
No term is omitted because its PHR burden is zero.

The prespecified primary support rule is:

1. assignment is midpoint;
2. within-collection BH `q <= 0.05`; and
3. global Holm-adjusted `p <= 0.05`.

BY and global Bonferroni are conservative context and do not replace that
rule. The overlap analysis is a paired boundary sensitivity, not a second
route to a primary claim.

## Prohibited primary code paths and output ordering

V7 primary inference has no RNG, sampling, permutation, interval placement,
rigid-block translation, same-arm control, terminal-stratum control,
subtelomeric-matched control, adjacent control, or random spatial-control code
path. The implementation must fail validation if prohibited inference tokens
or imports appear in its executable analysis source.

The execution order is fixed:

1. validate V6 checksums/schema and construct the eligible-copy partition;
2. enumerate all frozen term counts and complete finite-population tails;
3. adjust the complete frozen p-value families and fix term-level decisions;
4. emit the V7-specific exact PHR contributor ledger and mapping/cohort audits;
5. generate source/community summaries from already inferred term rows; and
6. validate partitions, counts, hypotheses, contributor semantics, ordering,
   and the prohibited-code-path gate before writing the final report.

Communities and source labels are descriptive post-inference groupings. They
do not create, merge, select, filter, or test hypotheses.

## Interpretation boundary

V7 answers which exact ontology terms have excess annotation-bearing physical
copy burden in the CHM13 PHR compartment relative to the rest of the eligible
CHM13 genome. It does not test transcription, expression, translation,
protein activity, dosage effect, retained pseudogene function, biological
independence of nearby copies, or population prevalence.
