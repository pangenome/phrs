# Manuscript candidate text: CHM13 V6 physical-copy ontology analysis

**Draft status.** This file supplies candidate Results, Methods, and
limitations language for later manual integration. It does not edit
`submission/paper.tex`. V6 replaces all V1--V5 inferential conclusions for
this estimand. Exact values, contributor ledgers, terminology constraints, and
the independently PASSed evidence index are in
[`FINAL_V6_REPORT.md`](FINAL_V6_REPORT.md).

## Candidate Results text

### Physical-copy ontology analysis identifies specific functional-term enrichments

We tested whether CHM13 PHRs contain excess physical-copy burden for
target-blind GO and Reactome terms. Each coordinate-anchored, source-mapped
copy was retained as one unit: N term-bearing copies contributed N even when
they shared a name, functional source, or post-inference family label. Among
31,235 direct or ancestor term hypotheses, 143 exact term rows met the
prespecified primary rule requiring both the collection-wise BH Monte Carlo
upper bound and the all-ontology single-step maxT upper bound to be at most
0.05; 31,092 were certified non-passing and none remained Monte Carlo
unresolved. The supported rows comprised 88 GO biological-process, 30 GO
molecular-function, 22 GO cellular-component, and 3 Reactome hypotheses.

Supported terms formed several biologically coherent but post-inference
display groups. Endosomal trafficking, exocyst, secretion, and actin
organization were represented by separately tested exact terms. Direct WASH
complex (`GO:0071203`), exocyst (`GO:0000145`), and Arp2/3
complex-mediated actin nucleation (`GO:0034314`) each had 10 physical copies
against a null mean of 0.497 and median of 0 (95% placement interval 0--2), a
median difference of 10, and a smoothed burden ratio of 10.53. WASH complex
had BH q=4.87 x 10^-4 (sequential bounds 4.64 x 10^-4 to 2.96 x 10^-3), BY
q=0.00412, and collection and global maxT p=1 x 10^-5 (sequential upper bound
6.08 x 10^-5). Direct exocytosis (`GO:0006887`) had 11 copies versus null
mean 1.053 and median 1 (0--3), ratio 7.41, BH q=0.00133, and global maxT
p=0.00658; direct endosomal transport (`GO:0016197`) had 12 versus mean 2.880
and median 3 (2--5), ratio 3.70, q=0.00133, and global maxT p=0.01626.

These burdens remain physical-copy counts. The WASH-complex, exocyst, and
Arp2/3 rows share ten coordinate-distinct contributors pointing to the WASHC1
functional source: WASHC1 itself plus LOC100996442, WASH3P, WASH4P, two
WASH5P copies, WASH6P, WASH7P, and two WASH8P copies. One WASHC1 source thus
contributed 10, not one. Exocytosis adds one VAMP7 copy, and endosomal
transport adds two SNX18-source copies. The label “WASH/endosomal-actin
system” summarizes these exact term rows after inference; it is not a tested
WASH-family hypothesis. Pseudogene/source mappings record auditable
annotation-bearing coordinates, not expression or activity.

DNA and chromosome-associated GO terms and three Reactome UPR-related rows
shared a second exact contributor pattern. Direct DNA helicase activity
(`GO:0003678`) had 10 copies versus a null mean of 0.298 and median 0 (0--2),
a median difference of 10 and smoothed ratio of 13.16. Its BH q was 0.00110
(sequential bounds 0.00100--0.00667), BY q was 0.01024, and its collection and
global maxT p-values were 1 x 10^-5 (sequential upper bound 6.08 x 10^-5).
Reactome direct XBP1(S) activates chaperone genes (`R-HSA-381038`) had the
same 10-copy burden and null distribution, BH q=0.00104, BY q=0.00896, and
global maxT p=1 x 10^-5. Its ancestors IRE1alpha activates chaperones
(`R-HSA-381070`) and Unfolded Protein Response (`R-HSA-381119`) were separate
certified rows. Each of these selected rows receives one unit from each of ten
DDX11L coordinate copies mapped by directed evidence to DDX11. The shared
source and name group do not collapse the burden to one and do not establish
expression or functional equivalence among the copies.

Other specific supported rows had distinct contributors and effect sizes.
Direct rRNA binding (`GO:0019843`) had 11 RPL23A-source physical copies versus
a null mean of 0.518 and median 0 (0--2), a ratio of 11.29, BH q=0.00110, and
global maxT p=1 x 10^-5. The broader ancestor sexual reproduction
(`GO:0019953`) had 33 copies versus a null mean of 11.513 and median 11
(6--17), a median difference of 22 and ratio of 2.79; its contributors were 19
SEPTIN14-source, 11 RPL23A-source, 2 TRPC6-source, and 1 TUBB8 copy, and its
global maxT p was 0.04071 (sequential bounds 0.03916--0.04230). Ancestor
ribonucleoside triphosphate phosphatase activity (`GO:0017111`) had 32 copies
versus mean 8.867 and median 9 (3--15), difference 23 and ratio 3.47; it had
19 SEPTIN14-, 10 DDX11-, and three tubulin-source contributors and global
maxT p=0.04071. These source/name concentrations explain term-level burdens
after inference; no RPL23A, SEPTIN14, TUBB8, or other gene-family enrichment
was tested.

The five prespecified name cohorts illustrate why term inference and family
description must remain separate. DUX4/DUX4L had 107 genome copies, 68 PHR
midpoint copies, and 65 ontology-eligible PHR copies. All 65 contributed to
Reactome Zygotic genome activation (`R-HSA-9819196`), rather than collapsing
to one DUX4 source, but the exact term was non-passing (observed 65, null mean
14.045 and median 7 [0--54], p=0.00751, BH q=0.575, global maxT p=0.988).
DDX11L had 10 of 10 PHR copies eligible and all 10 contributed to the
certified DNA-helicase row. TUBB8 had 2 of 7 PHR copies eligible; those two
contributed to the 28-copy organelle-organization row, which was non-passing
under global maxT despite BH q=0.00133. OR4F had 4 of 11 PHR copies eligible;
those four were part of a 7-copy Reactome olfactory-receptor-expression row
that was non-passing (q=1, global maxT p=0.99992). The name-defined WASH cohort
had 9 of 9 PHR copies eligible; those nine plus LOC100996442 produced the
certified 10-copy WASH-complex burden. These are exact term results with named
contributors, not DUX, DDX11L, TUBB8, OR4F, or WASH family tests.

Source coverage was incomplete and was reported explicitly. Of 61,312 CHM13
physical copies, 31,966 (52.14%) were ontology eligible: 20,405 through an
exact-self source and 11,561 through an explicit directed human `Related
functional gene` record. In the 402-copy primary PHR-midpoint set, 187
(46.52%) were eligible and 215 were not: 200 unsupported, 8 ambiguous, 5
type-only, and 2 unresolved records failed closed without emitting an
ontology edge. Thus absence from an enriched term is not evidence that an
unmapped copy or biological function is absent from a PHR.

### V6 supersedes earlier copy-ontology conclusions

For this estimand, V6 replaces all V1--V5 inferential conclusions. No earlier
term selection, burden, p-value, adjusted value, uncertainty classification,
family-level claim, or system-level conclusion is carried into the V6 result.
The V6 claim set consists only of independently validated exact term rows in
the frozen 31,235-hypothesis catalog. System labels and contributor families
are post-inference descriptions and cannot create, rescue, merge, or demote a
term-level inference.

## Candidate Methods text

### Genome-wide source mapping and physical-copy ontology matrix

We represented the CHM13v2.0 RefSeq/Liftoff gene annotation as 61,312
coordinate-anchored physical records. Separately identified GFF rows were
retained even when coordinates and strand coincided, preserving 12 records
beyond coordinate-plus-strand uniqueness. Every row had positive
`physical_copy_cn`; the genome-wide sum was 61,312. Before joining PHR
coordinates, we froze a source assignment and evidence row for every physical
record. An exact-self route admitted a copy when its own exact registry
identity carried GO or Reactome terms. A second route admitted an explicit,
directed human `Related functional gene` record when it uniquely named a
term-bearing functional source. Ambiguous, type-only, unresolved, and
unsupported records emitted no functional source. Names, aliases, presumed
pseudogene parents, sequence similarity, and gene- or sequence-family labels
were not annotation routes.

We reconstructed direct GO assertions from the frozen human gene2go release
and direct Reactome assignments from Reactome v96, then propagated only their
recorded ontology/pathway ancestors. Direct and ancestor edges remained
separately labeled. The frozen matrix contained 31,966 eligible physical
copies, 626,048 direct term-edge CN, and 2,929,709 direct-plus-ancestor edge
CN. Exact-self mapping accounted for 20,405 copies and explicit-related
mapping for 11,561. The target-blind catalog contained 19,181 GO biological-
process, 6,261 GO molecular-function, 2,677 GO cellular-component, and 3,116
Reactome direct-or-ancestor hypotheses. The PHR BED was opened only after the
source matrix and hypothesis catalog were frozen.

### Physical-copy burden and exact contributors

For each exact `(collection, relation, term_id)` hypothesis, we summed
`physical_copy_cn` over term-bearing records whose midpoint lay inside a CHM13
PHR. Any positive base overlap was retained as a paired boundary sensitivity.
Coordinate-distinct physical records were never deduplicated by gene name,
stable identifier, functional source, ontology source assertion, biotype,
family, or tandem organization. Thus N source-mapped copies contributed N.
Duplicate database evidence did not create duplicate physical copies: each
copy-term edge was represented once with its provenance. We retained exact
copy IDs, coordinates, CN, source assignments, directed evidence record IDs,
term distances, and midpoint/overlap membership for every observed term
contributor.

The statistic estimates regional annotation-bearing physical-copy burden. It
does not measure RNA or protein abundance, molecular activity, functional
dosage, or equivalence between a related copy and its annotation-bearing
source. Gene-name, source, and family summaries were joined only after every
term statistic and inferential value had been fixed.

### Spatial randomization and multiplicity

The primary null translated each of 37 complete rigid PHR blocks uniformly
among valid integer starts on its observed chromosome arm and within its
prespecified terminal-distance stratum: 0--0.5 Mb, 0.5--1 Mb, 1--2 Mb, 2--5
Mb, or 5 Mb to the arm end. Block widths, component order, gaps, annotation
clusters, source assignments, ontology edges, and physical-CN weights were
held fixed. We used 99,999 common valid placements generated by PCG64DXSM with
master seed 2026071301 and spawn key `[0]`; midpoint and any-overlap burdens
were recounted on the identical placement coordinates.

All tests were one-sided enrichment tests with ties counted as exceedances and
plus-one empirical p-values. We calculated two-sided 95% and sequential 98.75%
plus-one-transformed Clopper--Pearson intervals. BH and BY were applied within
each of the four ontology collections over direct plus ancestor rows.
Single-step Westfall--Young maxT values were computed within each collection
and across all 31,235 terms from the complete common null vector. Primary
support required both the BH sequential upper endpoint and the global-maxT
sequential upper endpoint to be at most 0.05. The complete 99,999-placement
screen resolved all primary decisions, leaving no selective-extension
candidate and no Monte Carlo-unresolved row.

### Post-inference summaries and independent validation

We used shared exact contributors and ontology labels to organize supported
rows for presentation only after inference. A display system did not merge
terms, reduce multiplicity, receive a p-value, or constitute a family test.
Likewise, functional-source and gene-name concentrations identified physical
contributors but were not substitute burden statistics.

An independent validator imported none of the production mapping, sampling,
inference, or prior validation modules. It rebuilt the genome-wide matrix from
the raw GFF, frozen per-copy evidence, raw gene2go, GO OBO, Reactome relations
and all-level assignments, cytobands, and PHR BED; regenerated all 99,999
placements; reconstructed representative null columns; and recalculated
observed burdens, exact contributors, empirical uncertainty, BH/BY, and
collection/global maxT arithmetic. All 46 required checks passed, including
61,312 assignment/evidence bijections, 31,235 hypotheses, 62,470 result rows,
17,579 exact any-overlap contributors, and computational-provenance hashes.

## Candidate limitations text

This analysis concerns annotation-bearing physical-copy burden in regional
intervals of one CHM13 assembly; it does not estimate human population or
pangenome prevalence. Mapping is deliberately incomplete. Only 187 of 402
primary PHR-midpoint copies had an admissible term-bearing functional source,
whereas 215 failed closed as unsupported, ambiguous, type-only, or unresolved.
These records remain physical copies, and their missing ontology edges must
not be interpreted as absence of expression, molecular function, or pathway
biology.

An explicit directed `Related functional gene` assertion provides auditable
annotation provenance but does not prove functional equivalence, expression,
translation, or activity of a related copy. The presence of a pseudogene or
noncoding copy in a term burden similarly records a coordinate and source edge
rather than a functional assay. The observed statistic is distinct from RNA
or protein abundance, dosage, chromosome-contact frequency, or causal
coordination.

GO and Reactome are hierarchical, and many of the 143 supported term rows
share identical physical contributors. Consequently, the result is 143
supported exact term rows, not 143 independent biological systems. Semantic
system labels and driver groups were formed after inference and have no
system- or family-level p-value. Broader ancestor terms can contain mixed
contributors, so a named cohort's contribution must not be substituted for
the total term burden. Finally, raw and BH values alone are insufficient when
the global maxT safeguard fails, and no contributor pattern can rescue a
non-passing exact term.

## Minimum reporting guardrails

- State that the supported objects are 143 exact
  `(collection, relation, term_id)` rows; do not call them 143 independent
  systems.
- State plainly that N source-mapped physical copies contribute N. Never
  replace physical CN with unique gene, source, or family counts.
- Keep direct and ancestor ontology rows distinct. Do not merge their burdens
  or adjusted values.
- Identify “WASH/endosomal-actin,” “DDX11/helicase-UPR,” and other system
  labels as post-inference summaries without their own tests.
- Identify WASHC1-, DDX11-, RPL23A-, SEPTIN14-, DUX-, TUBB8-, and OR4F-related
  groups as contributor descriptions, not family enrichment.
- Distinguish functional annotation-bearing copy burden from expression,
  translation, activity, dosage, functional equivalence, and causal biology.
- Report coverage: 187/402 primary PHR-midpoint copies are eligible and
  215/402 lack an admissible source. Never present missing annotation as
  biological absence.
- Keep DUX ZGA, TUBB8 organelle organization, and OR4F receptor expression as
  exact `CERTIFIED_NONPASS` examples; their descriptive burdens do not create
  family hypotheses.
- Do not restore or blend any V1--V5 inferential conclusion for this estimand.
