# CHM13 physical-copy enrichment: final author review

**Review status:** the analysis implementation, calibration, frozen inputs, and
completed permutation runs are auditable; no biological enrichment claim is
ready for the manuscript. The current result is best summarized as a validated
method with unresolved biological candidates. This report does not modify or
recommend automatic modification of `submission/paper.tex`.

The compact rerun guide is [REPRODUCIBILITY_README.md](REPRODUCIBILITY_README.md).
The versioned machine audit is
[RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json), and the exact checksums of the
source-sized result products are in
[SOURCE_OUTPUT_CHECKSUMS.tsv](results/v1/SOURCE_OUTPUT_CHECKSUMS.tsv).

## Executive decision

The final run used the frozen set of 37 CHM13 PHR intervals, a universe of
61,312 coordinate-distinct CHM13 gene loci, 402 midpoint-assigned target loci,
and 412 loci under the any-overlap sensitivity definition
([PREFLIGHT.json](results/v1/config/PREFLIGHT.json)). Each of the same-arm,
terminal cross-arm, and adjacent-region nulls completed 99,999 valid joint
region-set permutations, for 299,997 permutations in total and a minimum
plus-one Monte Carlo p-value of 0.00001
([RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json)).

The prespecified term-level numerical screen was:

> Primary midpoint assignment under the same-arm, terminal-stratum-preserving
> spatial null; within-collection BH q <= 0.05 **and** global maxT p <= 0.05.

Forty-one statistic-by-term rows, representing 36 unique terms, meet that
point-estimate screen; all 41 rows are marked `extension_required=1`, all 36
unique terms are `single_driver_sensitive=1`, and 8 of the 35 composition terms
that meet the screen are sensitive to the choice of spatial background
([RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json)). Thus there are **no resolved,
robust term-enrichment conclusions** in this run. The 41 rows are candidates for
extension and diagnosis, not discoveries
([primary_evidence.tsv](results/v1/final_tables/primary_evidence.tsv)).

The prespecified all-locus burden test also points to an excess: 402 midpoint
loci versus a same-arm null median of 246 and central 95% interval of 205--321,
with a smoothed burden ratio of 1.595 and 0/99,999 exceedances. However, its
empirical p=0.00001 is also marked for extension to 999,999 permutations, so it
is not yet a resolved confirmatory statement
([all_burden_results.tsv](results/v1/final_tables/all_burden_results.tsv)).

### Author-facing recommendation

- **Manuscript-ready now:** the descriptive CHM13 inventory, the physical-copy
  estimand, the three spatial backgrounds, the target-blind annotation and term
  filtering procedure, and the fact that the implementation passed its frozen
  engineering calibration. These claims should cite the analysis as methods or
  resource characterization, not as proof of biological enrichment.

- **Not manuscript-ready:** enrichment of total physical-copy burden, WASH or
  endosomal functions, helicase/DDX11L, DUX4/ZGA, interleukin signaling,
  pseudogenes, lncRNAs, miRNAs, or any other individual term. The limiting
  conditions are the mandatory Monte Carlo extension and copy/family influence;
  eight composition candidates additionally fail a background sensitivity.

- **Never confirmatory in this package:** the deduplicated-symbol ORA comparator,
  overlap assignment, terminal and adjacent backgrounds, leave-one-driver
  analyses, and narrative grouping of related terms. These are explicitly
  sensitivity or exploratory objects under the statistical contract
  ([STATISTICAL_SPEC.md](STATISTICAL_SPEC.md)).

## What is being estimated

The observation unit is a **physical CHM13 GFF gene locus**, not a unique gene
symbol. Two annotations with the same symbol at different coordinates are two
observations. A locus can contribute once to each applicable frozen term, but a
duplicate annotation edge cannot multiply the locus. This preserves the
quantity the biological question asks about: how many physical annotated copies
fall inside a set of PHR regions.

The analysis separates three term-level estimands:

- **Copy burden** counts physical loci bearing a term. It can rise because the
  regions contain more genes overall, because the term occupies a larger share
  of genes, or both.

- **Composition** is the term-bearing physical-copy count divided by the number
  of target loci annotated in that collection. It asks whether a term occupies
  an unusual share of the locally annotated copies. The reported composition
  odds ratio is a continuity-corrected descriptive effect size; it is not a
  model-based odds ratio.

- **Arm breadth** counts target arms with at least one term-bearing physical
  copy. It distinguishes a dispersed pattern from a dense single-arm array.

The all-locus burden test is separate from the term families. The exact
definitions, midpoint and overlap rules, empirical tails, and effect-size
formulae are frozen in [STATISTICAL_SPEC.md](STATISTICAL_SPEC.md); their
implementation is in [COPY_engine.py](COPY_engine.py).

This design does not claim duplicate copies are independent. The complete PHR
interval set is randomized against fixed genome annotation, so neighboring
copies, tandem arrays, duplicated blocks, and correlated term labels travel
together in the null distribution. That is why the historical weighted or
instance-expanded hypergeometric test is not the primary analysis.

## Spatial nulls and decision rules

### Primary same-arm null

Each PHR block is translated on its observed chromosome arm while retaining its
width, internal gaps, arm, and one of the frozen telomere-distance strata. Joint
placements preserve the same-arm non-overlap constraint. The primary candidate
space includes the observed placement, as required for exchangeability. Gene
density, GC, repeat content, and realized gene count are not matching variables.

For every highlighted term below, “primary” means this same-arm null plus
midpoint assignment. The highlighted threshold is always within-collection BH
q <= 0.05 **and** global maxT p <= 0.05. BH is calculated separately for each
collection-by-statistic family; global maxT spans all primary midpoint term
hypotheses across collections and statistics. Family sizes are 4,732 GO BP,
1,648 GO MF, 1,044 GO CC, 1,038 HGNC-group, 2,166 Reactome, and 19 biotype
hypotheses per statistic, totaling 10,647 term hypotheses
([multiplicity_families.tsv](results/v1/final_tables/multiplicity_families.tsv)).

The screen above is necessary but not sufficient for release. The contract also
requires Monte Carlo confidence to preserve the decision. A result with fewer
than 100 exceedances, an adjusted value near its threshold, or a confidence
bound capable of changing the decision must be extended without restarting the
RNG stream. All provisional rows triggered that rule
([primary_evidence.tsv](results/v1/final_tables/primary_evidence.tsv)).

### Terminal cross-arm sensitivity

This background allows movement to another arm in a frozen class defined by arm
orientation, acrocentric status, and autosome/sex-chromosome status, while
retaining block geometry and telomere-distance stratum. Its collection-wise BH
q <= 0.05 is a robustness check only; it cannot create a primary result.

### Adjacent-region sensitivity

This background translates each block within the immediately proximal 5-Mb
same-arm annulus, retaining geometry and non-overlap. Its collection-wise BH
q <= 0.05 is also a robustness check only. The deterministic nearest proximal
placement is an effect comparator and has no p-value. Definitions and
non-estimability rules for both sensitivity backgrounds are in
[STATISTICAL_SPEC.md](STATISTICAL_SPEC.md).

### Boundary sensitivity

Midpoint and any-positive-base overlap assignments were recounted from the same
placements. None of the 41 provisional primary rows has
`boundary_sensitive=1`, a favorable but non-rescuing diagnostic
([primary_evidence.tsv](results/v1/final_tables/primary_evidence.tsv)). The
overlap assignment is not included in the primary global maxT family and cannot
be used to declare a discovery.

## Annotation universe and coverage

The physical universe is built only from CHM13v2 RefSeq/Liftoff `gene` features.
The frozen PHR BED is zero-based and half-open; the GFF is converted from
one-based closed coordinates. The source paths, byte counts, and input SHA-256
values are recorded in
[PROVENANCE.tsv](analysis_ready/PROVENANCE.tsv), and all prepared-object hashes
are in [MANIFEST.sha256](analysis_ready/MANIFEST.sha256).

Biotype annotation is complete by construction for all 61,312 physical loci.
Stable identifiers map for 60,411 loci (98.530%); GO covers 20,363 loci
(33.212%), HGNC groups cover 26,913 (43.895%), and Reactome covers 11,365
(18.536%)
([coverage_by_biotype.tsv](outputs/coverage_by_biotype.tsv)). Coverage is highly
biotype-dependent: protein-coding loci have 94.302% GO, 77.824% HGNC-group, and
56.737% Reactome coverage, whereas generic pseudogenes have 0.730%, 8.166%, and
0.000%, respectively
([coverage_by_biotype.tsv](outputs/coverage_by_biotype.tsv)).

This distinction is central to interpretation. The biotype analysis can describe
pseudogene and noncoding copies comprehensively, but GO, HGNC, and especially
Reactome cannot be treated as comprehensive functional annotation of those
classes. “No pathway enrichment” is therefore not evidence that pseudogene or
noncoding copies lack biological roles.

Term maps were frozen without using target membership. Eligibility required at
least 5 genome-wide physical loci on at least 2 arms
([PREFLIGHT.json](results/v1/config/PREFLIGHT.json)). GO assertions with `NOT`
qualifiers were excluded, obsolete identifiers were checked against the frozen
ontology, and direct GO annotations were retained without automatic ancestor
propagation. The copy-to-term adapter preserves both coordinate-derived
`copy_id` and engine `locus_id`; its exact collection checksums are in
[engine_terms/MANIFEST.json](engine_terms/MANIFEST.json).

## Calibration and run audit

The independent calibration passed every frozen engineering gate. Under the
complete synthetic regional null at nominal 0.05, observed rejection rates were
0.059 for total burden, 0.0286 for term-copy burden, 0.0424 for composition,
0.0112 for breadth, 0.048 for any BH rejection, and 0.021 for any global-maxT
rejection
([combined.json](calibration_results/job-20260713T162803Z/combined.json)). The
positive failure control rejected 10,536 of 40,000 clustered-null draws
(0.2634; 95% interval 0.2591--0.2677), confirming that the suite detects the
anti-conservative instance-expanded hypergeometric approach
([combined.json](calibration_results/job-20260713T162803Z/combined.json)).

Calibration validates mechanics under the simulated conditions; it does not
validate a biological result or prove that a spatial exchangeability assumption
is true for every real annotation landscape. The complete independent review is
[VALIDATION_REPORT.md](VALIDATION_REPORT.md), while the machine-readable
calibration authority is
[combined.json](calibration_results/job-20260713T162803Z/combined.json).

For the final data, every mode passed the batch, placement, candidate-space,
acceptance, cached-array, and zero-denominator gates. The primary minimum
candidate count was 250,000, adjacent minimum candidate count was 4,500,001,
and the minimum acceptance rate was 1.0 in all modes
([RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json)). The audit covers 11,730
transient files and 19,414,754,100 bytes, with per-file digests retained in
[transient_file_checksums.tsv.gz](results/v1/transient_file_checksums.tsv.gz)
and totals in [RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json).

The initial full permutation jobs completed their inferential products but were
cancelled when the original leave-one scan proved non-scaling. The exact cached
subtraction finalizer was checked byte-for-byte against the reference diagnostic
implementation before use. The successful finalizer array and aggregation job
completed with zero exit codes; all attempts and retained scheduler accounting
are machine-readable in
[SLURM_COMPLETION.json](results/v1/SLURM_COMPLETION.json),
[slurm_accounting.tsv](results/v1/slurm_accounting.tsv), and
[FINALIZER_EQUIVALENCE.json](results/v1/FINALIZER_EQUIVALENCE.json).

## Primary results

### Total physical-copy burden

The following table keeps the null and threshold next to every highlighted
burden result. The all-locus test is the single prespecified burden test, so its
threshold is empirical p <= 0.05 with a Monte Carlo confidence decision; it is
not part of term-level BH or global maxT.

| Assignment/background | Observed loci | Null median (central 95%) | Smoothed ratio | Empirical tail | Null and threshold | Decision | Source |
| --- | ---: | ---: | ---: | ---: | --- | --- | --- |
| Midpoint, primary | 402 | 246 (205--321) | 1.595 | 0/99,999; p=0.00001 | Same-arm, terminal-stratum spatial null; p <= 0.05 plus resolved MC confidence | Directionally positive; extend to 999,999 | [machine row](results/v1/final_tables/all_burden_results.tsv) |
| Midpoint, terminal sensitivity | 402 | 201 (158--271) | 1.964 | 0/99,999; p=0.00001 | Terminal-matched cross-arm spatial null; sensitivity p <= 0.05 plus resolved MC confidence | Same direction; extend to 999,999 | [machine row](results/v1/final_tables/all_burden_results.tsv) |
| Midpoint, adjacent sensitivity | 402 | 193 (127--263) | 2.078 | 0/99,999; p=0.00001 | Same-arm proximal-annulus spatial null; sensitivity p <= 0.05 plus resolved MC confidence | Same direction; extend to 999,999 | [machine row](results/v1/final_tables/all_burden_results.tsv) |

The burden pattern is the broadest candidate because it is not tied to one
functional annotation source and has the same direction under every background.
It is still not a final inferential claim: all three tails sit at the Monte Carlo
floor with zero exceedances and are explicitly marked for extension
([all_burden_results.tsv](results/v1/final_tables/all_burden_results.tsv)).

### Representative term-level candidates

These rows illustrate the main candidate stories without pretending they are
independent findings. Every row uses **primary same-arm midpoint; BH q <= 0.05
and global maxT p <= 0.05**, and every row is unresolved because
`extension_required=1`.

| Candidate | Observed and effect | Multiplicity result | Driver/background diagnosis | Null and threshold | Source |
| --- | --- | --- | --- | --- | --- |
| Wiskott-Aldrich Syndrome protein family, HGNC_GROUP:14, copy burden and breadth | 7 copies on 7 arms; null median 0, central 95% 0--2; count ratio 8.492 | raw p=0.00001; BH q=0.00346; global maxT p=0.01216 | Both sensitivity backgrounds retain BH support; 6/7 copies map back to the same HGNC family and leave-largest-family-out leaves 1 arm | Primary same-arm midpoint; BH q <= 0.05 and global maxT p <= 0.05 | [evidence](results/v1/final_tables/primary_evidence.tsv), [all nulls](results/v1/final_tables/all_term_results.tsv.gz), [drivers](results/v1/final_tables/primary_driver_summary.tsv), [leave-one](results/v1/final_tables/primary_leave_one_sensitivity.tsv) |
| WASH complex, GO:0071203, composition | 3/23 GO-CC-annotated copies on 3 arms; null median 0; descriptive composition OR 22.366, log2 effect 4.483 | raw p=0.00001; BH q=0.000522; global maxT p=0.00001 | Terminal q=0.000614 and adjacent q=0.001740; 2/3 copies are in HGNC_GROUP:14; removing them leaves 1 arm | Primary same-arm midpoint; BH q <= 0.05 and global maxT p <= 0.05 | [evidence](results/v1/final_tables/primary_evidence.tsv), [cross-null](results/v1/final_tables/composition_cross_null.tsv.gz), [drivers](results/v1/final_tables/primary_driver_group_counts.tsv), [leave-one](results/v1/final_tables/primary_leave_one_sensitivity.tsv) |
| DNA helicase activity, GO:0003678, composition | 2/21 GO-MF-annotated copies, both displayed as DDX11L8 on 2 arms; null median 0; log2 effect 3.955 | raw p=0.00001; BH q=0.001268; global maxT p=0.00001 | Terminal q=0.001831 and adjacent q=0.003296; either physical copy removal reduces breadth to 1 | Primary same-arm midpoint; BH q <= 0.05 and global maxT p <= 0.05 | [evidence](results/v1/final_tables/primary_evidence.tsv), [cross-null](results/v1/final_tables/composition_cross_null.tsv.gz), [drivers](results/v1/final_tables/primary_driver_group_counts.tsv), [leave-one](results/v1/final_tables/primary_leave_one_sensitivity.tsv) |
| Zygotic genome activation, R-HSA-9819196, composition | 1/10 Reactome-annotated copies, DUX4 on chr4_q; null median 0; log2 effect 3.487 | raw p=0.00001; BH q=0.000333; global maxT p=0.00192 | Terminal q=0.000433 and adjacent q=0.001160; removing DUX4 leaves zero observed copies | Primary same-arm midpoint; BH q <= 0.05 and global maxT p <= 0.05 | [evidence](results/v1/final_tables/primary_evidence.tsv), [cross-null](results/v1/final_tables/composition_cross_null.tsv.gz), [drivers](results/v1/final_tables/primary_driver_group_counts.tsv), [leave-one](results/v1/final_tables/primary_leave_one_sensitivity.tsv) |
| Interleukin-9 signaling, R-HSA-8985947, composition | 1/10 Reactome-annotated copies, IL9R on chrX_q; null median 0; log2 effect 3.487 | raw p=0.00001; BH q=0.000333; global maxT p=0.00422 | Terminal q=0.000433 and adjacent q=0.000471; removing the sole copy erases the signal | Primary same-arm midpoint; BH q <= 0.05 and global maxT p <= 0.05 | [evidence](results/v1/final_tables/primary_evidence.tsv), [cross-null](results/v1/final_tables/composition_cross_null.tsv.gz), [drivers](results/v1/final_tables/primary_driver_group_counts.tsv), [leave-one](results/v1/final_tables/primary_leave_one_sensitivity.tsv) |

The coherent recurrence of WASH/endosomal labels is biologically legible, but
it is not independent replication: overlapping ontology terms annotate the same
small set of physical loci. DDX11L8, DUX4, and IL9R candidates are similarly
annotations of one or two physical-copy drivers. The proper reporting unit is
therefore the driver module, not the number of ontology labels.

## Robustness across spatial backgrounds

Among the 35 provisional primary composition terms, 27 retain direction and
within-collection BH q <= 0.05 under both the terminal cross-arm and adjacent
backgrounds, while 8 fail at least one such sensitivity
([RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json)). The plot
[null_sensitivity_composition.svg](results/v1/plots/null_sensitivity_composition.svg)
is a diagnostic visualization; all exact values and flags are in
[composition_cross_null.tsv.gz](results/v1/final_tables/composition_cross_null.tsv.gz).

The eight failures are listed explicitly below. Every row first met the primary
same-arm midpoint screen of BH q <= 0.05 and global maxT p <= 0.05. Sensitivity
support is judged at within-collection BH q <= 0.05 under the named background;
sensitivities have no primary global maxT value.

| Collection / term | Primary BH q; global maxT p | Terminal BH q | Adjacent BH q | Failed sensitivity | Source |
| --- | ---: | ---: | ---: | --- | --- |
| GO BP GO:0001501, skeletal system development | 0.001392; 0.01445 | 0.004579 | 0.196954 | Adjacent | [machine row](results/v1/final_tables/composition_cross_null.tsv.gz) |
| GO BP GO:0010507, negative regulation of autophagy | 0.001392; 0.00240 | 0.001753 | 1.000000 | Adjacent | [machine row](results/v1/final_tables/composition_cross_null.tsv.gz) |
| GO BP GO:0016064, immunoglobulin-mediated immune response | 0.001392; 0.04234 | 1.000000 | 0.003786 | Terminal | [machine row](results/v1/final_tables/composition_cross_null.tsv.gz) |
| GO BP GO:0019221, cytokine-mediated signaling pathway | 0.001392; 0.04234 | 0.001753 | 1.000000 | Adjacent | [machine row](results/v1/final_tables/composition_cross_null.tsv.gz) |
| GO CC GO:0000145, exocyst | 0.000522; 0.00073 | 0.110200 | 0.001740 | Terminal | [machine row](results/v1/final_tables/composition_cross_null.tsv.gz) |
| GO CC GO:0005770, late endosome | 0.001816; 0.03605 | 0.007693 | 0.790540 | Adjacent | [machine row](results/v1/final_tables/composition_cross_null.tsv.gz) |
| GO MF GO:0004896, cytokine receptor activity | 0.001268; 0.02053 | 0.001831 | 1.000000 | Adjacent | [machine row](results/v1/final_tables/composition_cross_null.tsv.gz) |
| Reactome R-HSA-451927, Interleukin-2 family signaling | 0.000333; 0.00422 | 0.000433 | 0.061202 | Adjacent | [machine row](results/v1/final_tables/composition_cross_null.tsv.gz) |

All eight retain a positive displayed effect under the sensitivity backgrounds;
their failure is loss of FDR support rather than direction reversal
([composition_cross_null.tsv.gz](results/v1/final_tables/composition_cross_null.tsv.gz)).
That still makes them background-dependent and unsuitable for a confirmatory
claim.

## Pseudogene and noncoding findings

Biotype is the one collection with complete physical-locus coverage, so these
results are more interpretable as locus inventory than GO or pathway results.
They still use the same spatial permutation framework and multiplicity rules.

| Biotype observation | Primary same-arm midpoint result | Threshold beside result | Interpretation | Source |
| --- | --- | --- | --- | --- |
| Generic pseudogene | 204/402 target loci (50.746%); burden null median 77 and count ratio 2.548; burden BH q=0.0000633 but global maxT p=0.38145; composition BH q=0.00019 but global maxT p=0.98574; breadth BH q=0.000095 but global maxT p=0.84833 | Primary term rule: BH q <= 0.05 **and** global maxT p <= 0.05 | Strong descriptive copy excess, but fails global maxT for every statistic and is not a primary call | [machine rows](results/v1/final_tables/all_term_results.tsv.gz) |
| Transcribed pseudogene | 18/402 target loci (4.478%); burden null median 5 and count ratio 3.091; burden BH q=0.0000633 but global maxT p=0.34061; composition BH q=0.051015 and global maxT p=0.99862 | Primary term rule: BH q <= 0.05 **and** global maxT p <= 0.05 | Burden is exploratory; composition fails BH and global maxT | [machine rows](results/v1/final_tables/all_term_results.tsv.gz) |
| lncRNA | 104/402 target loci (25.871%); burden null median 69 and count ratio 1.495; burden BH q=0.0000633 but global maxT p=0.73634; composition log2 effect -0.161 and BH q=1.0 | Primary term rule: BH q <= 0.05 **and** global maxT p <= 0.05 | More physical lncRNA copies accompany the larger total locus burden; there is no evidence for a larger lncRNA share | [machine rows](results/v1/final_tables/all_term_results.tsv.gz) |
| miRNA | 51/402 target loci (12.687%); burden null median 11 and count ratio 3.119; burden BH q=0.031303 but global maxT p=0.99867; composition BH q=0.65189 | Primary term rule: BH q <= 0.05 **and** global maxT p <= 0.05 | Exploratory burden pattern only; no globally adjusted or composition result | [machine rows](results/v1/final_tables/all_term_results.tsv.gz) |
| rRNA | 0/402 target loci despite 982 genome loci; one-sided enrichment p=1.0 | Primary enrichment test; p <= 0.05 and applicable multiplicity criteria | No enrichment claim; absence may reflect the target coordinates and the CHM13 gene-feature inventory, not a general statement about subtelomeric rDNA | [machine rows](results/v1/final_tables/all_term_results.tsv.gz), [coverage](outputs/coverage_by_biotype.tsv) |

The pseudogene composition signal is stable at within-biotype BH q=0.00019
under the primary, terminal, and adjacent backgrounds, but sensitivity support
cannot substitute for its failed primary global maxT test
([all_term_results.tsv.gz](results/v1/final_tables/all_term_results.tsv.gz)). The
appropriate current wording is descriptive: “half of the midpoint-assigned
CHM13 PHR gene features are annotated as generic pseudogenes,” not “PHRs are
enriched for pseudogenes.”

## Copy-driver diagnostics

The target-blind family map assigns each eligible locus to its smallest
genome-wide frozen HGNC group, with lexical tie-breaking. It maps 26,913 loci;
34,399 unmapped loci receive distinct singleton groups, and target membership
was not used
([PREFLIGHT.json](results/v1/config/PREFLIGHT.json)). No qualifying
target-independent sequence-identity cluster map was available. The engine's
per-locus identity fallback is retained only for schema completeness and **must
not** be described as sequence-cluster evidence
([RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json)).

Every one of the 36 unique provisional terms is marked single-driver-sensitive
under the frozen rule: removal of one family or copy reverses direction, reduces
breadth to one, moves a raw p-value above 0.05, or halves the absolute effect
([RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json)). The following grouping
explicitly lists all provisional terms by their largest family/copy driver. It
is a diagnostic synthesis of
[primary_evidence.tsv](results/v1/final_tables/primary_evidence.tsv),
[primary_driver_group_counts.tsv](results/v1/final_tables/primary_driver_group_counts.tsv),
and [primary_driver_summary.tsv](results/v1/final_tables/primary_driver_summary.tsv).

| Largest driver | Provisional unique terms attached to that driver | Diagnosis |
| --- | --- | --- |
| **HGNC_GROUP:14, duplicated WASH/Wiskott-Aldrich family** | protein targeting to lysosome; exocytosis; endosome organization; endosomal transport; endocytic recycling; Arp2/3 complex-mediated actin nucleation; regulation of Arp2/3 complex-mediated actin nucleation; retrograde transport, endosome to Golgi; early endosome to late endosome transport; alpha-tubulin binding; WASH complex; early endosome membrane; recycling endosome; recycling endosome membrane; Wiskott-Aldrich Syndrome protein family | **Explicit duplicated-family-dominated set.** The family-level result has 7 observed copies, 6 assigned to HGNC_GROUP:14, and breadth falls from 7 arms to 1 after removing that family; many ontology terms have 2 of 3 observed copies in the same family. [Machine evidence](results/v1/final_tables/primary_driver_group_counts.tsv), [summary](results/v1/final_tables/primary_driver_summary.tsv), [leave-one](results/v1/final_tables/primary_leave_one_sensitivity.tsv) |
| HGNC_GROUP:1331, WASHC1-centered single-copy driver | negative regulation of autophagy; positive regulation of pseudopodium assembly; regulation of protein ubiquitination; phosphatidylinositol 3-kinase inhibitor activity; exocyst; late endosome; autophagosome; centriole | Mostly one WASHC1-bearing physical locus, sometimes paired with a second locus; leave-one diagnoses copy dependence. [Machine evidence](results/v1/final_tables/primary_driver_group_counts.tsv), [summary](results/v1/final_tables/primary_driver_summary.tsv) |
| HGNC_GROUP:716, IL9R-centered single-copy driver | immunoglobulin-mediated immune response; cytokine-mediated signaling pathway; interleukin-9-mediated signaling pathway; cytokine receptor activity; Interleukin-2 family signaling; Interleukin-9 signaling | One IL9R physical locus drives each label; several also fail a spatial-background sensitivity. [Machine evidence](results/v1/final_tables/primary_driver_group_counts.tsv), [cross-null](results/v1/final_tables/composition_cross_null.tsv.gz) |
| DDX11L8 coordinate copies, represented as locus-singleton family fallbacks | DNA metabolic process; DNA helicase activity; helicase activity; hydrolase activity acting on acid anhydrides in phosphorus-containing anhydrides | Two coordinate-distinct loci with the same displayed symbol occur on two arms; removing either reduces breadth to one. This is duplicated-symbol influence, not validated sequence-cluster evidence. [Machine evidence](results/v1/final_tables/primary_driver_group_counts.tsv), [leave-one](results/v1/final_tables/primary_leave_one_sensitivity.tsv) |
| DUX4 singleton | negative regulation of G0 to G1 transition; Zygotic genome activation | One chr4_q DUX4 locus; removal erases the observed term count. [Machine evidence](results/v1/final_tables/primary_driver_group_counts.tsv), [leave-one](results/v1/final_tables/primary_leave_one_sensitivity.tsv) |
| HGNC_GROUP:715, SHOX singleton | skeletal system development | One chrX_p SHOX locus; removal erases the observed term count, and the adjacent-background BH q is 0.196954. [Machine evidence](results/v1/final_tables/primary_driver_group_counts.tsv), [cross-null](results/v1/final_tables/composition_cross_null.tsv.gz) |

The duplicated WASH family is the clearest case requested for explicit warning:
ontology redundancy makes one physical family appear as many terms. Those terms
should be shown as one candidate module if the analysis is extended, not as a
long list of independent enrichments.

## Comparator analysis

The conventional comparator deduplicates the universe and query by symbol and
performs ordinary over-representation analysis. It tests a different estimand,
discards spatial correlation and physical-copy multiplicity, and has
`primary_inference=0` on every row. Ten of 10,647 comparator rows have BH
q <= 0.05, but none can confirm, refute, or select a copy-aware regional result
([RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json),
[COMPARATOR_deduplicated_symbol_ORA_not_copy_aware.tsv.gz](results/v1/final_tables/COMPARATOR_deduplicated_symbol_ORA_not_copy_aware.tsv.gz)).

The historical instance-expanded hypergeometric method is not provided as an
inferential comparator because calibration deliberately demonstrates its Type-I
inflation under clustered copies. It is a failure control only.

## Confirmatory versus exploratory ledger

| Object | Status | Permitted claim |
| --- | --- | --- |
| Frozen target and physical-copy inventory | Confirmatory descriptive object | Exact CHM13 interval/locus counts and annotation coverage, with no population extrapolation |
| Primary same-arm midpoint all-locus burden | Prespecified confirmatory test, unresolved | Report as a candidate excess pending stream-preserving extension; do not use “significantly enriched” |
| Primary same-arm midpoint terms passing BH and global maxT point estimates | Prespecified confirmatory candidates, unresolved | Report as candidates only; all require extension and all are driver-sensitive |
| Terminal cross-arm and adjacent-region results | Prespecified sensitivity | State whether a primary candidate persists or fails; never create a discovery |
| Any-overlap assignment | Prespecified sensitivity | Boundary robustness only |
| Family/copy removal | Diagnostic | Identify concentration and fragile terms; never rescue or create significance |
| Pseudogene/noncoding patterns failing global maxT | Exploratory/descriptive | Report target inventory and effect estimates with failed threshold stated beside them |
| Deduplicated-symbol ORA | Comparator only | Illustrate how the estimand changes after symbol deduplication |
| Biological grouping of ontology labels | Exploratory synthesis | Use for follow-up prioritization, not multiplicity-adjusted evidence counting |

## Limitations

1. **Single-reference scope.** The analysis describes physical loci in CHM13v2
   and its frozen PHR coordinates. It does not estimate prevalence or copy-number
   variation across HPRCv2 individuals.

2. **Conditional exchangeability.** The primary p-values require observed block
   locations to be exchangeable with valid translations after conditioning on
   arm, geometry, mask, and telomere-distance stratum. Calibration validates the
   sampler under synthetic conditions, not the biological truth of this
   assumption.

3. **Annotation incompleteness.** GO and pathway coverage is sparse for
   pseudogene and noncoding classes. Functional absence cannot be inferred from
   an unannotated physical copy.

4. **Direct annotation semantics.** GO assignments are direct frozen assertions;
   automatic ancestor propagation is not used. Related ontology labels overlap
   heavily and should not be narrated as independent biological replications.

5. **Family map is diagnostic, not sequence homology.** HGNC groups cover only a
   subset of loci, unmapped loci become singletons, and no acceptable independent
   identity-cluster map was available.

6. **Monte Carlo precision is unfinished at candidate tails.** The completed
   runs meet the planned initial depth, but the candidates sit at or near the
   resolution where the statistical contract requires continuation to unused RNG
   draws.

7. **One-sided scope.** The inferential package tests enrichment. Zero counts or
   apparent deficits are not confirmatory depletion results.

8. **Physical copies are not expression or function.** The estimand is annotated
   coordinate-level gene features. It does not establish transcription,
   protein activity, dosage, or a mechanistic role in PHR biology.

## Claims matrix for manuscript review

| Proposed wording | Ready? | Reason / required action |
| --- | --- | --- |
| “The frozen CHM13 PHR set contains 402 midpoint-assigned physical gene loci.” | **Yes, descriptive** | Exact coordinate recount is independently reproduced. [Machine input](analysis_ready/chm13_phr_gene_midpoint.tsv) |
| “The permutation framework preserves physical-copy and local spatial correlation.” | **Yes, methods** | This is the implemented estimand and randomization unit, validated by exact sampler and recount tests. [Engine](COPY_engine.py), [calibration](calibration_results/job-20260713T162803Z/combined.json) |
| “CHM13 PHRs are enriched for physical gene copies.” | **No** | Primary, terminal, and adjacent burdens are directionally positive but all have 0/99,999 exceedances and require continuation to 999,999 draws. [Machine rows](results/v1/final_tables/all_burden_results.tsv) |
| “WASH/endosomal functions are enriched in PHRs.” | **No** | Provisional point-estimate rule passes, but the module is dominated by duplicated WASH-family copies and every term requires extension. [Evidence](results/v1/final_tables/primary_evidence.tsv), [drivers](results/v1/final_tables/primary_driver_group_counts.tsv) |
| “DDX11L helicase functions are enriched across PHR arms.” | **No** | Two coordinate copies of one displayed symbol drive the terms; removal of either reduces breadth to one, and extension remains required. [Drivers](results/v1/final_tables/primary_driver_group_counts.tsv), [leave-one](results/v1/final_tables/primary_leave_one_sensitivity.tsv) |
| “DUX4 links PHRs to zygotic genome activation.” | **No** | The result is one DUX4 locus and vanishes when it is removed. [Drivers](results/v1/final_tables/primary_driver_group_counts.tsv), [leave-one](results/v1/final_tables/primary_leave_one_sensitivity.tsv) |
| “PHRs are enriched for pseudogenes.” | **No** | The descriptive count is large and sensitivity BH values are small, but every primary biotype statistic fails global maxT. [Machine rows](results/v1/final_tables/all_term_results.tsv.gz) |
| “Generic pseudogenes comprise 204 of 402 midpoint-assigned CHM13 PHR loci.” | **Yes, descriptive with scope** | Biotype is complete for the physical-locus universe; do not convert this inventory into an enrichment or pangenome claim. [Machine rows](results/v1/final_tables/all_term_results.tsv.gz), [coverage](outputs/coverage_by_biotype.tsv) |
| “Symbol-level ORA confirms the copy-aware analysis.” | **No** | The comparator tests a different, non-spatial estimand and has `primary_inference=0`. [Comparator](results/v1/final_tables/COMPARATOR_deduplicated_symbol_ORA_not_copy_aware.tsv.gz) |

## Required next analysis before any enrichment claim

1. Continue all three saved RNG streams to at least the prespecified extension
   depth, without regenerating prior placements. Reaggregate and require the
   Monte Carlo confidence-bound decision checks to pass.

2. Re-evaluate the all-locus burden and all term candidates from the complete
   extended outputs. Do not extend only selected terms; placements and complete
   multiplicity families must remain joint.

3. Treat the duplicated WASH family as a single biological candidate module in
   narrative review, while retaining the full term families for statistical
   correction.

4. If a target-independent sequence-identity clustering becomes available,
   freeze it before inspecting revised results and rerun the diagnostic only. Do
   not retrofit clusters selected from target membership.

5. Consider replication in HPRCv2 physical assemblies as a separate study. A
   CHM13 result, even if resolved, should not be generalized to population
   prevalence without that analysis.

## Artifact index

### Design and implementation

- [STATISTICAL_SPEC.md](STATISTICAL_SPEC.md) — normative estimands, nulls,
  multiplicity, diagnostics, calibration gates, and extension policy.
- [COPY_engine.py](COPY_engine.py) — spatial permutation engine.
- [build_inputs.py](build_inputs.py) — CHM13 physical-locus and interval input
  builder.
- [build_term_maps.py](build_term_maps.py) — frozen source-to-copy annotation
  builder.
- [prepare_engine_terms.py](prepare_engine_terms.py) — coordinate-keyed adapter
  from `copy_id` to `locus_id`.
- [validate_engine_run.py](validate_engine_run.py) — black-box coordinate recount
  and p/BH/maxT validator.
- [calibration_suite.py](calibration_suite.py) and
  [combine_calibration.py](combine_calibration.py) — synthetic calibration and
  aggregate gate evaluation.

### Frozen inputs

- [data/chm13.phrs.bed](../../../data/chm13.phrs.bed) — target PHR intervals.
- [data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz](../../../data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz)
  — physical CHM13 gene features.
- [analysis_ready/MANIFEST.sha256](analysis_ready/MANIFEST.sha256) — prepared
  input hashes, sizes, and row counts.
- [sources/SOURCE_MANIFEST.tsv](sources/SOURCE_MANIFEST.tsv) — external source
  releases, retrieval metadata, and checksums.
- [engine_terms/MANIFEST.json](engine_terms/MANIFEST.json) — final per-collection
  term-map counts and checksums.
- [results/v1/config/final.args](results/v1/config/final.args) — exact final
  collection and family-map arguments.
- [results/v1/config/PREFLIGHT.json](results/v1/config/PREFLIGHT.json) — frozen
  preflight and design parameters.

### Final results and diagnostics

- [all_burden_results.tsv](results/v1/final_tables/all_burden_results.tsv) — all
  locus-burden rows across backgrounds and assignments.
- [all_term_results.tsv.gz](results/v1/final_tables/all_term_results.tsv.gz) —
  complete term results for every mode, assignment, collection, and statistic.
- [primary_evidence.tsv](results/v1/final_tables/primary_evidence.tsv) — rows
  meeting the provisional primary point-estimate screen.
- [multiplicity_families.tsv](results/v1/final_tables/multiplicity_families.tsv)
  — family sizes and threshold counts.
- [composition_cross_null.tsv.gz](results/v1/final_tables/composition_cross_null.tsv.gz)
  — primary, terminal, and adjacent composition comparison.
- [primary_driver_summary.tsv](results/v1/final_tables/primary_driver_summary.tsv),
  [primary_driver_group_counts.tsv](results/v1/final_tables/primary_driver_group_counts.tsv),
  and [primary_leave_one_sensitivity.tsv](results/v1/final_tables/primary_leave_one_sensitivity.tsv)
  — copy/family concentration and paired removal diagnostics.
- [null_sensitivity_composition.svg](results/v1/plots/null_sensitivity_composition.svg)
  — diagnostic plot, not a discovery figure.
- [RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json),
  [SOURCE_OUTPUT_CHECKSUMS.tsv](results/v1/SOURCE_OUTPUT_CHECKSUMS.tsv), and
  [SLURM_COMPLETION.json](results/v1/SLURM_COMPLETION.json) — final machine gates,
  product integrity, and scheduler provenance.
