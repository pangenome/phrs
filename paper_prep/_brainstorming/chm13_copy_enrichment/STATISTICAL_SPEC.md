# Statistical specification: copy-number-aware regional enrichment in CHM13 PHRs

- **Status:** reviewed analysis contract
- **Scope:** CHM13 v2.0 physical gene loci and the 37 intervals in
  `data/chm13.phrs.bed`
- **Primary inference:** region-set randomization on the reference genome
- **Not in scope:** manuscript text, HPRCv2 projected annotations, and biological
  interpretation

This document is normative for the preparation, permutation-engine, analysis,
and independent-validation tasks in `chm13_copy_enrichment`. An implementation
must either conform to every item marked **MUST** or record a versioned protocol
deviation before inspecting enrichment results. A result from a nonconforming
run is diagnostic only.

## 1. Decision and rationale

The inferential question is whether a *set of CHM13 genomic regions* contains an
unusual burden or composition of annotated physical loci relative to exchangeable
region sets. Therefore:

1. Each physical CHM13 GFF3 gene locus is an observation. Copies at different
   coordinates remain different observations even when their symbols, aliases,
   products, or term memberships are identical.
2. The randomization unit is a complete PHR interval (or a rigid same-arm block
   of intervals), never a locus, gene symbol, term label, or individual copy.
3. Each randomized interval retains its width and chromosome-arm stratum. Moving
   continuous intervals keeps nearby loci, tandem arrays, duplicated blocks, and
   overlapping term labels together. This is the mechanism that represents
   spatial and copy correlation under the null.
4. Total locus burden and functional composition are separate estimands. Term
   counts alone mix the two and MUST NOT be described as composition enrichment.

This design does not assume that duplicate copies are independent Bernoulli
draws. It deliberately retains their genomic clustering.

## 2. Frozen inputs and provenance

Every run MUST record SHA-256 checksums, byte sizes, modification dates, parser
version, Git commit, command line, and output-schema version for:

- `data/chm13.phrs.bed`, interpreted as zero-based, half-open BED intervals;
- `data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz`, the sole source of physical
  CHM13 gene loci;
- the CHM13 sequence dictionary/chromosome sizes and the frozen chromosome-arm
  boundaries used by input preparation;
- each frozen locus-to-term map and each frozen family/identity-cluster map;
- any exclusion mask used to construct placement candidates; and
- this specification and the implementation used to execute it.

The prepared target MUST contain 37 interval rows before any union operation.
The interval-set checksum and total union length MUST be reported. Overlapping
target rows, if any, are retained as provenance but coalesced for locus counting
so that a physical locus is not counted twice merely because two target rows
overlap.

The primary locus assignment is by gene-feature midpoint: locus \(i\) belongs to
region set \(R\) when its genomic midpoint lies in the half-open union of \(R\).
The midpoint of a feature with half-open coordinates \([s_i,e_i)\) is
\(s_i + \lfloor(e_i-s_i)/2\rfloor\). The prespecified assignment sensitivity is
any positive-base overlap between \([s_i,e_i)\) and the union of \(R\). Boundary,
strand, pseudogene, ncRNA, duplicate-symbol, and multi-copy cases MUST be tested
by input preparation.

No locus may be propagated to a homologous chromosome, inferred from a gene
symbol, or imported from an HPRCv2 assembly. Missing term annotation means
"unannotated for that collection," not zero copies and not evidence against the
term.

## 3. Analysis objects and notation

Let:

- \(G=\{1,\ldots,N\}\) be all eligible physical CHM13 gene loci, with stable
  `locus_id`, coordinates, arm, biotype, symbol, and optional family and
  identity-cluster identifiers;
- \(R_0=\{r_1,\ldots,r_J\}\), \(J=37\), be the observed PHR interval set;
- \(R_b\), \(b=1,\ldots,B\), be valid randomized region sets;
- \(I_i(R)\) indicate primary midpoint assignment of locus \(i\) to the union
  of \(R\);
- \(A_{it}\) indicate membership of locus \(i\) in frozen term \(t\). A locus
  is counted at most once for a term, even if multiple annotation paths map it
  to that term; and
- \(a(i)\), \(f(i)\), and \(h(i)\) identify locus \(i\)'s chromosome arm,
  prespecified gene family, and prespecified sequence-identity cluster.

Term maps MUST be frozen without reference to whether a locus lies in \(R_0\).
The annotation collections (for example, GO BP, GO MF, GO CC, or a pathway
database) remain separate collections. Ancestor propagation, evidence-code
filters, obsolete-ID replacement, synonym handling, and database versions MUST
be fixed in the term-map report.

Terms may be filtered only by target-blind genome-wide rules. The default
testability rule is at least five annotated physical loci on at least two CHM13
chromosome arms in \(G\). A collection may declare a stricter rule before the
target is joined. Terms cannot be filtered by their observed PHR count,
unadjusted p-value, effect direction, or apparent biological interest.

## 4. Estimands and statistics

All primary tests are one-sided enrichment tests. Prespecified depletion analyses,
if executed, are exploratory and constitute distinct multiplicity families.

### 4.1 Total physical-copy burden

The genome-wide total-burden estimand is the excess number of physical gene loci
covered by the complete PHR region set:

\[
  L(R)=\sum_{i\in G} I_i(R), \qquad
  \Delta_L=L(R_0)-\operatorname{median}_{b}L(R_b).
\]

Report \(L(R_0)\), the null median and central 95% randomization interval,
\(\Delta_L\), and the smoothed burden ratio

\[
  \mathrm{BR}=\frac{L(R_0)+0.5}
                    {\operatorname{mean}_b L(R_b)+0.5}.
\]

The raw integer \(L\) is the test statistic. Gene biotypes may be examined as a
prespecified secondary burden family, but the all-locus burden is the sole
genome-wide burden test.

For every retained functional term, its total physical-copy burden is

\[
  C_t(R)=\sum_{i\in G} I_i(R)A_{it}.
\]

The integer \(C_t\) is also a primary test statistic. It asks whether PHRs cover
more physical copies carrying term \(t\), allowing both higher local gene burden
and different term composition to contribute. Report \(C_t(R_0)\), its null
median and central 95% interval, the count difference, and the smoothed count
ratio \((C_t(R_0)+0.5)/(\operatorname{mean}_b C_t(R_b)+0.5)\). Matching
candidate regions on realized gene count would condition away both \(L\) and part
of \(C_t\), and is prohibited in the primary analysis.

### 4.2 Term composition

For collection \(c\), let \(E_c\subseteq G\) be loci eligible for its
composition denominator: loci with at least one retained annotation in that
collection. Define

\[
  D_c(R)=\sum_{i\in E_c} I_i(R).
\]

For retained term \(t\), \(A_{it}=0\) outside \(E_c\), so the \(C_t\) from
Section 4.1 can equivalently be written

\[
  C_t(R)=\sum_{i\in E_c} I_i(R)A_{it},
\]

and the composition statistic is

\[
  Q_t(R)=
  \begin{cases}
    C_t(R)/D_c(R), & D_c(R)>0;\\
    0,             & D_c(R)=0.
  \end{cases}
\]

Thus total annotated-copy burden can vary across placements, while \(Q_t\)
asks what share belongs to term \(t\). The empirical test compares \(Q_t(R_0)\)
with the same statistic for randomized region sets; no binomial independence is
assumed.

Report the observed physical-copy count, denominator, observed proportion, null
median proportion, and a descriptive Jeffreys-smoothed composition odds ratio:

\[
  \mathrm{COR}_t=
  \frac{(C_t(R_0)+0.5)/(D_c(R_0)-C_t(R_0)+0.5)}
       {(\widetilde C_t+0.5)/(\widetilde D_c-\widetilde C_t+0.5)},
\]

where \(\widetilde C_t\) and \(\widetilde D_c\) are null medians. This ratio is
an effect-size summary, not a model-based odds ratio and not its own test.
Define the composition effect \(E_t=\log_2(\mathrm{COR}_t)\). Raw-count ratios
may be reported as burden-sensitive descriptions but MUST be labeled separately
from composition.

### 4.3 Arm breadth

Term breadth is the number of target arms on which at least one physical copy
with term \(t\) is covered:

\[
  W_t(R)=\sum_{a\in\mathcal A_0}
       \mathbf 1\!\left\{\sum_{i:a(i)=a} I_i(R)A_{it}>0\right\},
\]

where \(\mathcal A_0\) is the fixed set of arms represented by the observed
37 interval templates. The empirical test uses \(W_t\). Report observed breadth,
null median and 95% interval, and the arms contributing to the observed count.
Breadth distinguishes a term distributed over arms from one driven by many
copies in a single local array.

## 5. Placement units

The preparation step MUST construct a `placement_block` table. On each arm,
overlapping or abutting target intervals form a block. A block stores the ordered
interval widths and intervening gaps in telomere-oriented coordinates. The entire
vector is translated rigidly during randomization; widths and internal gaps never
change. With one target interval on an arm, that interval is its own block.

At every placement, loci are assigned anew from coordinates. The engine MUST NOT
sample a number of genes and then choose loci, expand symbols into copies, shuffle
term labels, or independently move members of a duplicated array. Continuous
placement automatically preserves the local co-occurrence structure of whatever
duplicated block or neighboring gene cluster is sampled.

Candidate windows may overlap target PHRs, including the original placement.
Excluding the observed placement would break symmetry. Distinct blocks placed on
the same arm may not overlap one another unless the corresponding observed blocks
overlap; joint rejection or exact sampling without replacement is required.

## 6. Primary null: same-arm, terminally constrained region randomization

### 6.1 Candidate space

Coordinates are converted to distance from the telomeric end of the assigned
arm. Candidate placements MUST:

- remain on the same named CHM13 chromosome arm as their template;
- retain the exact interval widths and internal block gaps;
- remain entirely within the arm and outside the frozen exclusion mask;
- retain the template's terminal-distance stratum; and
- satisfy the joint non-overlap rule for multiple blocks on an arm.

Terminal-distance strata are fixed, in base pairs from the telomeric arm end, as
midpoint distance \([0,0.5\text{ Mb})\), \([0.5,1\text{ Mb})\),
\([1,2\text{ Mb})\), \([2,5\text{ Mb})\), and \([5\text{ Mb},\text{arm end})\).
These cut points are target-blind and MUST NOT be widened after observing term
counts. A rigid block is assigned using its span midpoint. A candidate start is
valid when the translated block's midpoint is in the same stratum as the
observed template and every component interval satisfies all other constraints.

If an exclusion mask is used, candidate placement additionally must have the
same excluded-base fraction as the template to within one percentage point.
If no biological exclusion mask is justified, the mask is empty apart from arm
bounds, and the run records that fact. Gene density, term counts, GC content,
and repeat content are not matching variables in the primary null.

### 6.2 Sampling algorithm

For each Monte Carlo replicate:

1. Process arms in canonical chromosome order. Within an arm, process placement
   blocks in stable `block_id` order.
2. Draw uniformly from all valid integer translations of the block in its
   candidate space. An exact candidate list or an exactly uniform integer-start
   sampler with documented rejection is acceptable. A fixed grid is not
   acceptable unless all valid starts are on that grid by definition and the
   biological estimand is explicitly changed.
3. If two same-arm blocks collide, redraw the complete set for that arm. Do not
   favor whichever block happened to be sampled first.
4. Combine the arm placements into one joint 37-template region set \(R_b\),
   coalesce its union only for overlap assignment, and recompute all locus,
   term, breadth, family, and identity-cluster counts from genomic coordinates.
5. Use the same joint \(R_b\) for every term and statistic. This retains
   dependence among terms and enables maxT control.

The candidate count per block, number of rejected draws, and acceptance rate per
arm MUST be retained. Fewer than 100 distinct valid translations for any primary
block makes that block's primary randomization non-estimable under this contract;
the analysis stops rather than widening the stratum post hoc.

### 6.3 Exchangeability assumption

The sharp primary null is that, conditional on the fixed arm, block geometry,
terminal-distance stratum, arm bounds, and exclusion-mask overlap, the observed
block placement is exchangeable with every valid integer translation. Genomic
annotations and their spatial correlation are fixed. The null does not require
independent loci, independent duplicate copies, independent term memberships,
or similar gene density among candidates.

Uniform placement within the declared candidate space must be scientifically
defensible. A rejection of this assumption by the placement-bias checks in
Section 13 invalidates the primary p-values.

## 7. Prespecified sensitivity backgrounds

Sensitivity analyses reuse the frozen inputs, statistics, term filters, and
observed target. They do not create new confirmatory discoveries. Their purpose
is to expose dependence on the conditioning set.

### 7.1 Matched terminal, cross-arm background

This background relaxes exact-arm conditioning while preserving terminal
context. Arms are partitioned before analysis by p/q orientation,
acrocentric/non-acrocentric status, and autosome/sex-chromosome status. Within
each partition, interval-block templates are assigned without replacement to a
uniformly sampled arm and then uniformly translated within the same fixed
terminal-distance stratum used by the primary null. Widths, internal gaps, masks,
and joint non-overlap remain fixed.

If a partition has too few eligible arms for without-replacement assignment, its
templates remain on their observed arms and the loss of cross-arm sensitivity is
reported. No arm may receive more blocks than allowed by its feasible candidate
space. This analysis asks whether results persist when arm-specific annotation
landscapes are not conditioned away.

### 7.2 Same-arm adjacent-region background

For each observed block, define the proximal annulus from the block's non-telomeric
edge out to 5 Mb farther from the telomere, truncated at the centromere-side arm
boundary. Candidate placements are exact-width rigid translations wholly inside
that annulus and the frozen mask. Sample one valid placement per block uniformly,
with the same joint same-arm non-overlap rule. The original PHR placement is not
part of this annulus.

In addition to Monte Carlo results, report the deterministic immediately proximal
comparator: the first valid equal-width placement after the observed block. It is
an effect-size comparator only and receives no p-value. If fewer than 100 annulus
placements exist, report the comparator and mark the adjacent empirical analysis
non-estimable; do not extend the annulus after looking at results.

The terminal cross-arm and adjacent-region backgrounds test different assumptions
and MUST be reported separately, never pooled into one null distribution.

## 8. Empirical p-values and Monte Carlo precision

For an enrichment statistic \(T\in\{L,C_t,Q_t,W_t\}\), the unadjusted empirical
p-value is

\[
  \widehat p_T=\frac{1+\sum_{b=1}^{B}\mathbf 1\{T(R_b)\geq T(R_0)\}}
                      {B+1}.
\]

Ties use `>=`. No zero p-values are permitted. The primary same-arm run and each
estimable sensitivity background require at least **99,999 valid joint region-set
permutations**; rejected placements do not count toward \(B\). A pilot may use
19,999 permutations for runtime planning but is not reportable inference.

Extend a run, preserving its RNG stream, to at least 999,999 valid permutations
when any of the following holds:

- a raw, BH-adjusted, or maxT-adjusted value is within 20% of its reporting
  threshold;
- the observed exceedance count is below 100;
- the attainable p-value floor is too coarse for the planned family size; or
- the Monte Carlo interval would change a pass/fail decision.

For every p-value, report the exceedance count, \(B\), and a two-sided 95%
Clopper--Pearson interval for the underlying exceedance probability. A result can
pass a significance threshold only when the interval's upper bound is below that
threshold; an indeterminate result triggers more permutations, not optimistic
rounding. For maxT, construct the interval from the maximum-statistic exceedance
count. For a BH decision, apply BH to the vector of upper confidence bounds as a
conservative Monte Carlo decision check in addition to reporting BH on the point
estimates.

## 9. Zero counts, degeneracy, and missing annotations

- If the observed enrichment count is zero, its one-sided enrichment p-value is
  1 because every randomized count is at least zero.
- If all null and observed values are identical, set raw p, q, and maxT p to 1,
  report zero null variance, and label the hypothesis non-informative.
- If the null count is always zero but the observed count is positive, the p-value
  is initially \(1/(B+1)\). The run MUST be extended to 999,999 and candidate
  generation audited before reporting; this pattern can reveal a placement bug or
  support mismatch.
- If \(D_c(R)=0\), set \(Q_t(R)=0\) for the one-sided statistic and retain the
  replicate. Do not drop it. Report how often this occurs; a rate above 1% for a
  collection is a calibration failure for composition inference in that
  collection.
- Effect-size denominators use the stated 0.5 continuity correction. It affects
  display only, never the randomization statistic or p-value.
- Loci missing from a particular annotation collection are excluded only from
  that collection's composition denominator. They remain in total burden.
- A locus with multiple terms contributes once to every applicable term. This
  induces correlation that is preserved by using common region-set permutations.

## 10. Multiplicity families and decision rules

Multiplicity families are declared before target statistics are inspected.

1. **Primary burden family:** the all-locus \(L\) test is one prespecified test.
   Secondary biotype-specific burden tests form one separate family.
2. **Primary term-copy-burden families:** one family per annotation collection,
   containing every target-blind testable term under the primary same-arm null
   and midpoint assignment, using \(C_t\).
3. **Primary composition families:** one family per annotation collection,
   containing every target-blind testable term under the primary same-arm null and
   midpoint assignment.
4. **Primary breadth families:** one family per annotation collection, containing
   the same testable terms but using \(W_t\).
5. **Global primary family:** all testable primary term-copy-burden, composition,
   and breadth hypotheses across all collections. It excludes the single
   all-locus burden test, assignment sensitivities, background sensitivities, and
   driver diagnostics.
6. **Sensitivity families:** any p-values produced for a sensitivity background
   are adjusted within collection and statistic exactly as above but are labeled
   sensitivity q-values. They cannot replace a failed primary test.
7. **Exploratory depletion:** if run, it forms direction-specific families of the
   same structure and is never combined with enrichment after seeing the sign.

Apply Benjamini--Hochberg separately within every collection-by-statistic family
to its empirical p-values and report all terms in each family, not only selected
rows. BH is the within-collection FDR summary, contingent on successful null-FDR
calibration in Section 13. If BH exceeds its calibration tolerance, report BY
q-values as the dependence-robust fallback and suppress BH discovery labels.

Also calculate single-step Westfall--Young maxT adjusted p-values. For every
hypothesis \(t\), pool observed and randomized statistics, then define

\[
  Z_{bt}=\{T_t(R_b)-\bar T_t\}/s_t, \quad b=0,\ldots,B,
\]

where \(\bar T_t\) and \(s_t\) use all \(B+1\) values symmetrically. If
\(s_t=0\), set every \(Z_{bt}=0\) when all values are equal; if only the observed
value exceeds a constant randomized value, set \(Z_{0t}=+\infty\) and
\(Z_{bt}=0\) for \(b>0\). For family \(F\), let
\(M_b=\max_{u\in F}Z_{bu}\), and compute

\[
  \widehat p^{\max T}_t=
  \frac{1+\sum_{b=1}^{B}\mathbf 1\{M_b\ge Z_{0t}\}}{B+1}.
\]

Report collection-by-statistic maxT values and a global-primary maxT value. The
latter controls the probability of any false rejection under the complete joint
randomization null. Strong FWER language requires the additional subset-pivotality
assumption and MUST NOT be used unless independent validation justifies it.
Primary term evidence requires both the relevant within-collection FDR
\(q\le0.05\) and global maxT \(p\le0.05\); the two quantities are still reported
separately.

## 11. Copy-driver diagnostics

Diagnostics are mandatory for every term with primary copy-burden, composition,
or breadth BH \(q\le0.10\), or global maxT \(p\le0.10\), and may be tabulated
for all terms. They do not form a route to significance.

For observed term copies, report:

- counts by physical `locus_id`, symbol, chromosome arm, prespecified family, and
  identity cluster;
- the largest family fraction and largest identity-cluster fraction;
- the Herfindahl concentration and effective number
  \(n_{\mathrm{eff}}=1/\sum_g(C_{tg}/C_t)^2\) for both family and identity groups;
- breadth after removing the largest family and largest identity cluster; and
- whether multiple annotations describe the same physical locus.

Family and identity clusters MUST be frozen independently of target membership.
Gene symbol is a diagnostic grouping, not the primary observation unit. Missing
family/cluster identifiers become singleton groups rather than being discarded.
When \(C_t=0\), concentration and effective-number fields are `NA`, not zero.

## 12. Leave-one-family and identity-cluster sensitivity

For each diagnostic term, perform:

1. a leave-largest-family-out analysis;
2. a leave-largest-identity-cluster-out analysis; and
3. when there are at most 20 contributing groups, exhaustive leave-one-group-out
   analyses; otherwise, leave out groups in descending contribution order until
   at least 80% of observed term copies has been covered.

Removal applies to the same loci in the observed set and in every saved randomized
set. Recompute denominators, \(Q_t\), \(W_t\), empirical p-values, and effect
sizes without resampling placements. This preserves paired Monte Carlo precision.
These are influence analyses, not new discovery families.

Flag a term as **single-driver-sensitive** if removing any one family or identity
cluster reverses the effect direction, reduces arm breadth to one, moves the raw
empirical p above 0.05, or reduces \(|E_t|\) by at least 50%. A flag describes
concentration; it neither invalidates the physical-copy estimand nor warrants a
biological conclusion.

Repeat the full primary analysis under any-overlap locus assignment as a separate
assignment sensitivity, using identical placement coordinates and recalculated
overlaps. Terms whose direction changes or whose absolute log2 composition effect
changes by more than \(0.5\max(|E_t|,0.25)\) receive a **boundary-sensitive**
flag.

## 13. Calibration simulations and hard failure criteria

Calibration uses synthetic or pseudo-observed region sets only. It must not use
claimed biological enrichments as truth.

### 13.1 Exact sampler tests

On synthetic arms small enough to enumerate every valid translation, compare at
least 100,000 sampled placements with the exact distribution. Test uniform start
probabilities, both arm orientations, boundary inclusion, masks, rigid multi-
interval blocks, and same-arm collision handling. The chi-square or exact
multinomial test must have p >= 0.001, total-variation distance must be <=0.01,
and every valid start must be reachable. Widths, gaps, arms, and strata must match
their templates in 100% of draws.

On the real candidate spaces, compare observed sampling frequencies with known
candidate probabilities by arm and terminal-distance sub-bin. Standardized
frequency residuals must show no systematic edge trend; no absolute residual may
exceed 5 after multiplicity review. Report candidate counts, start-coordinate
quantiles, rejection rates, and endpoint frequencies. These checks specifically
detect interval-placement bias.

### 13.2 Regional sharp-null calibration

Generate at least 1,000 pseudo-observed region sets from each null sampler. For
each pseudo-observation, calculate all-locus burden, term-copy burden,
composition, breadth, BH, and maxT p-values against at least 9,999 exchangeable
region sets, using permutation-rank reuse when symmetry is maintained. Assess
rejection rates at 0.10, 0.05, and 0.01 with 95% binomial intervals, overall and
stratified by term prevalence, copy-cluster size, arm, interval width, and
terminal-distance sub-bin.

The calibration includes both complete-null term collections and mixed
collections with a prespecified 10% of terms receiving planted regional effects.
Use the latter to estimate BH false discovery proportion without treating the
planted terms as null. Effect locations and sizes are fixed before simulation and
span low, medium, and high copy-cluster sizes.

A calibration family fails if, at any threshold \(\alpha\):

- the lower 95% binomial confidence limit for its Type I error exceeds \(\alpha\);
- its point estimate exceeds \(\alpha+0.02\) (or \(1.5\alpha\), whichever is
  larger);
- the realized FDR for BH at q=0.05 exceeds 0.07 or has a lower 95% interval
  above 0.05; or
- the complete-null probability of any global maxT rejection at 0.05 exceeds
  0.07 or has a lower 95% interval above 0.05.

The engine must stop inference for a failed family. More permutations can resolve
Monte Carlo imprecision but cannot repair structural inflation.

### 13.3 Clustered-copy and terminal-gradient stress tests

Construct null synthetic annotations with all of the following: singleton loci,
tandem copies, long duplicated blocks containing several neighboring genes,
heavy-tailed family sizes, correlated term memberships, arm-specific gene density,
and strong but null terminal-distance gradients. Draw pseudo-targets as interval
sets using the declared null. The region-randomization test must meet the Type I
criteria above without fragmenting clusters.

As a required positive failure control, apply the historical instance-expanded
weighted hypergeometric test to gene- or region-selected clustered copies. The
calibration suite is itself inadequate if this deliberately mismatched method
does not show detectable Type I inflation in at least one high-copy scenario
(lower 95% binomial limit above 0.05 at nominal 0.05). Also include a planted
regional enrichment scenario to verify nonzero power; this is an engineering
test, not a biological result.

### 13.4 Additional stop conditions

No inferential result may be released if any of the following occurs:

- the input audit cannot reproduce 37 target rows or finds synthetic propagated
  loci, duplicate `locus_id` values, incompatible contig names, or unexplained
  coordinate/count differences;
- any primary placement block has fewer than 100 valid translations;
- fewer than 99,999 valid joint permutations complete;
- same-arm placements change arms, widths, internal gaps, or terminal strata;
- direct coordinate recounts disagree with cached locus counts;
- composition has \(D_c(R_b)=0\) in more than 1% of replicates;
- a sampler acceptance rate is below 20%, unless an exactly uniform enumerated
  candidate sampler replaces rejection sampling;
- rerunning the same saved placements changes any statistic or adjusted p-value;
- independent validation cannot reproduce raw counts, exceedances, BH families,
  or maxT maxima; or
- a reporting decision is unresolved by its Monte Carlo confidence interval at
  the maximum planned \(B\).

Failures are reported as failures, not as absence of enrichment.

## 14. Why weighted/expanded hypergeometric inference is not primary

The historical copy-weighted method replaces unique genes by copy counts and
uses

\[
  P\{X\ge q\},\quad
  X\sim\operatorname{Hypergeometric}(m,n,k),
\]

where \(m,n,k,q\) are totals of expanded gene instances. Modifying `phyper()`
parameters is algebraically equivalent to explicitly expanding each symbol into
its copies. That equivalence proves only that the two implementations calculate
the same hypergeometric tail. It does **not** prove that either tail is the null
distribution for a CHM13 regional-selection process.

The hypergeometric model samples individual instances without replacement. PHR
membership instead arises from genomic interval placement. A single interval can
select an entire tandem family or duplicated multi-gene block, so neighboring
copies have correlated inclusion indicators. Treating them as independent
instances understates the sampling variance and can make p-values anti-
conservative. The archived validation report documents this failure under
clustered gene-level sampling, including Type I error near 0.22 at nominal 0.05
and failed BH control in its simulations. Those numerical results are historical
diagnostics, not validation of the present CHM13 design.

Weighted hypergeometric inference would be valid only for a different experiment
in which physical copies themselves were exchangeable draws from a finite
instance population. That is not the stated experiment. Consequently neither
weighted parameters nor instance expansion may supply primary or sensitivity
p-values here, and no additional multiplicity correction can repair their wrong
sampling unit.

A conventional deduplicated-symbol ORA may be included only in a table headed
**Comparator: deduplicated-symbol ORA; not the copy-aware regional test**. Its
universe is frozen unique symbols, its query is unique symbols with at least one
target locus, and its ordinary hypergeometric/BH results answer whether symbols
are over-represented under symbol sampling. It discards physical copy burden and
spatial correlation, cannot confirm the primary analysis, and must not be used
to select terms for primary reporting. The weighted/expanded historical test may
appear only as a calibration failure control, not as an inferential comparator.

## 15. Reproducibility and random seeds

Use a counter-based or splittable generator whose algorithm and library version
are recorded. The reference implementation uses NumPy `PCG64DXSM` with
`SeedSequence` and these decimal master seeds:

| Purpose | Master seed |
|---|---:|
| Primary same-arm null | `2026071301` |
| Terminal cross-arm sensitivity | `2026071302` |
| Adjacent-region sensitivity | `2026071303` |
| Midpoint/overlap assignment checks | `2026071304` |
| Exact-sampler calibration | `2026071390` |
| Regional Type I/FDR/maxT calibration | `2026071391` |
| Clustered-copy/gradient stress tests | `2026071392` |
| Planted-effect engineering test | `2026071393` |

Parallel workers receive child streams by deterministic `SeedSequence.spawn` in
stable replicate-batch order, not by worker completion order or process ID. Save
the master seed, child spawn key, batch boundaries, generator name/version, and
the placement coordinates for every replicate. Extending \(B\) continues unused
child streams and never restarts or overwrites earlier draws.

Two executions with the same inputs, implementation version, thread count, and
seed manifest MUST produce byte-identical placement and result tables. A second
implementation used for independent validation may use a different RNG, but it
must reproduce observed counts exactly and agree with null summaries and p-values
within their prespecified Monte Carlo uncertainty.

## 16. Required outputs

The engine and final run MUST emit machine-readable, schema-versioned tables for:

- input checksums, reference coordinates, arm bounds, masks, and analysis options;
- target intervals and rigid placement blocks;
- candidate counts and placement-QC summaries;
- every randomized block coordinate or a lossless, auditable equivalent;
- per-replicate all-locus burden, term-copy burden, denominator, composition,
  breadth, family, and identity-cluster statistics;
- complete term results, including filtered/non-testable terms and reasons;
- empirical exceedances, \(B\), Monte Carlo intervals, BH q-values, collection
  maxT values, and global maxT values;
- burden, any-overlap, terminal, adjacent, family/identity leave-out, and
  deduplicated-symbol comparator results in separately labeled fields/files;
- calibration rejection rates and every hard-gate outcome; and
- software versions, Git commits, commands, RNG metadata, and elapsed time.

The human-readable report must describe methods and diagnostics without selecting
biological stories. It must state the number of hypotheses in every multiplicity
family, including terms with p=1, and distinguish `not testable`, `not significant`,
and `analysis invalid`.

## 17. Reconciliation with earlier drafts and archived work

At the reviewed baseline commit `0d9ce46`, neither
`paper_prep/_brainstorming/chm13_copy_enrichment/DESIGN.md` nor
`paper_prep/_brainstorming/chm13_copy_enrichment/TOOLCHAIN.md` exists in any
reachable branch history, and the target directory itself is absent. Therefore
there was no extant draft text to merge line by line. This specification is the
controlling design and explicitly supersedes any unmerged or subsequently
recovered `DESIGN.md`/`TOOLCHAIN.md` where they conflict.

The following archived materials were reviewed for failure history:

- `paper_prep/_brainstorming/README.md`, which marks the historical copy-weighted
  ORA corpus as off-target and noncanonical;
- `weighted_phyper_statistical_validation_report.md`, which demonstrates that
  algebraic instance-expansion equivalence does not ensure calibrated inference;
- `statistical_best_practices_weighted_ora.md` and
  `statistical_validation_conclusions_and_recommendations.md`, which warn against
  production use under clustered selection; and
- `copy_number_weighted_phyper_mathematical_formulation.md`, which documents the
  parameter mapping but not a region-selection null.

This contract supersedes their recommendations to shuffle genes or use
copy-weighted hypergeometric p-values. It retains only their useful distinction
between algebraic implementation equivalence and statistical validity. The new
toolchain must implement coordinate-based region-set randomization, shared joint
permutations, empirical tails, declared multiplicity families, and calibration
gates as specified here.

## 18. Implementation acceptance checklist

- [ ] The primary unit is a stable physical CHM13 locus; duplicate symbols remain
  separate rows.
- [ ] Exactly 37 target interval templates are accounted for, with union counting
  preventing accidental double counts.
- [ ] The primary sampler moves complete rigid intervals/blocks on the same arm,
  with exact widths, terminal strata, masks, and local gene correlation retained.
- [ ] Terminal cross-arm and adjacent-region backgrounds are separate sensitivity
  analyses.
- [ ] All-locus and term-specific physical-copy burden, term composition, and arm
  breadth use the distinct estimands above.
- [ ] Zero counts, unannotated loci, ties, null degeneracy, and non-estimable
  candidate spaces follow the explicit rules above.
- [ ] Every reportable null has at least 99,999 valid joint permutations, fixed
  seeds, plus-one p-values, exceedance counts, and Monte Carlo intervals.
- [ ] BH families, collection maxT, and global-primary maxT are complete and
  reproducible.
- [ ] Copy-driver, leave-one-family, leave-one-identity, and midpoint/overlap
  diagnostics reuse the same placements.
- [ ] Type I, FDR, maxT, clustered-copy, terminal-gradient, and interval-placement
  calibration gates pass; the legacy weighted method fails its positive failure
  control in a clustered-copy scenario.
- [ ] Deduplicated-symbol ORA, if present, is visibly labeled as a comparator and
  cannot drive term selection.
- [ ] No manuscript file is read as an analysis input or edited as an output.
