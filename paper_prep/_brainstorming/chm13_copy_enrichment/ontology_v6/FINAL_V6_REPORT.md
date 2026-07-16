# Final synthesis of CHM13 V6 physical-copy ontology enrichment

## Release decision

**V6 supports 143 exact GO/Reactome term rows as regionally enriched in CHM13
PHRs.** Support applies to the frozen `(collection, relation, term_id)` row,
not to a gene family or to a biological system assembled after seeing the
results. The primary midpoint analysis classified 143 of 31,235 hypotheses as
`CERTIFIED_PASS`, 31,092 as `CERTIFIED_NONPASS`, and none as
`MC_UNRESOLVED`. The supported rows comprise 25 direct and 63 ancestor GO
biological-process terms, 20 direct and 10 ancestor GO molecular-function
terms, 13 direct and 9 ancestor GO cellular-component terms, and 1 direct and
2 ancestor Reactome terms.

The independent validator returned an unqualified **PASS on 46/46 checks**,
including physical-copy semantics, source mapping, ontology closure, observed
burdens, null placements, Monte Carlo arithmetic, multiplicity, exact
contributors, and provenance. It reconstructed all 61,312 raw physical
copies, all 1,686,727 source-term rows, all 31,235 hypotheses and 62,470
midpoint/overlap result rows, and all 17,579 any-overlap contributor rows. See
[`V6_VALIDATION_REPORT.md`](V6_VALIDATION_REPORT.md) and the machine verdict
[`V6_VALIDATION.json`](V6_VALIDATION.json).

**This report retires every V1--V5 inferential conclusion for this estimand.**
No earlier positive, negative, unresolved, family-level, or system-level
conclusion should be combined with V6. V6's raw evidence inputs include a
hash-frozen remediated evidence file produced during earlier work, but no
V1--V5 term selection, burden, p-value, adjusted value, decision, or narrative
is carried forward.

## What the V6 result means

V6 tests an **annotation-bearing physical-copy burden**. For exact term
`T`, the observed midpoint statistic is

\[
S_T = \sum_i \mathrm{physical\_copy\_cn}_i
      I(i\text{ has the frozen edge to }T)
      I(i\text{ has its midpoint in a CHM13 PHR}).
\]

Every real V6 assignment row has `physical_copy_cn=1`. Therefore, **N
source-mapped coordinate copies carrying a term contribute N**, even if all N
copies point to one functional source, share a gene-name prefix, or would be
placed in one gene or sequence family. Source, unique gene, and family counts
are contributor metadata; they are not substitutes for `S_T`. The independent
adversarial fixture makes the distinction explicit: seven CN=1 coordinates
plus one CN=3 coordinate on one source/gene/family produce burden 10, whereas
each forbidden unique-source, unique-gene, and unique-family collapse produces
1.

In plain terms: N source-mapped copies contribute N.

This is not an expression, activity, dosage, or functional-equivalence test.
A source-mapped pseudogene copy contributes a coordinate and an auditable
ontology edge; it does not thereby demonstrate transcription, translation,
protein activity, pathway participation, or retained ancestral function.

## Exact supported term rows and effect sizes

The inferential claims are the exact term rows in
[`TERM_RESULTS.tsv.gz`](results/release/TERM_RESULTS.tsv.gz). Direct and
ancestor relations are separate hypotheses. The labels below retain that
distinction; no terms have been merged into a system-level test.

`Observed / genome` gives PHR-midpoint and genome-wide physical CN. The null
column is mean / median [2.5th--97.5th placement quantiles]. `Difference` is
observed minus the null median. The released burden ratio is the prespecified
smoothed effect `(observed + 0.5) / (null mean + 0.5)`, not a unique-gene or
family ratio. Contributor counts sum physical copies.

| Exact tested row | Observed / genome physical CN | Null mean / median [95% placement interval] | Difference | Smoothed burden ratio | Exact functional-source contributions |
|---|---:|---:|---:|---:|---|
| GO CC direct exocyst (`GO:0000145`) | 10 / 30 | 0.4969349693496935 / 0 [0--2] | 10 | 10.532281766431778 | WASHC1 source: 10 copies |
| GO CC direct WASH complex (`GO:0071203`) | 10 / 31 | 0.4969349693496935 / 0 [0--2] | 10 | 10.532281766431778 | WASHC1 source: 10 copies |
| GO BP direct Arp2/3 complex-mediated actin nucleation (`GO:0034314`) | 10 / 45 | 0.4969349693496935 / 0 [0--2] | 10 | 10.532281766431778 | WASHC1 source: 10 copies |
| GO BP direct exocytosis (`GO:0006887`) | 11 / 110 | 1.052770527705277 / 1 [0--3] | 10 | 7.406116869692901 | WASHC1: 10; VAMP7: 1 |
| GO BP direct endosomal transport (`GO:0016197`) | 12 / 179 | 2.8800488004880047 / 3 [2--5] | 9 | 3.698171457818974 | WASHC1: 10; SNX18: 2 |
| GO MF direct DNA helicase activity (`GO:0003678`) | 10 / 68 | 0.2978829788297883 / 0 [0--2] | 10 | 13.15982453391822 | DDX11 source: 10 copies |
| Reactome direct XBP1(S) activates chaperone genes (`R-HSA-381038`) | 10 / 75 | 0.2978829788297883 / 0 [0--2] | 10 | 13.15982453391822 | DDX11 source: 10 copies |
| GO MF direct rRNA binding (`GO:0019843`) | 11 / 240 | 0.5182951829518295 / 0 [0--2] | 11 | 11.293385447099782 | RPL23A source: 11 copies |
| GO BP ancestor sexual reproduction (`GO:0019953`) | 33 / 2,292 | 11.512745127451275 / 11 [6--17] | 22 | 2.788704800158167 | SEPTIN14: 19; RPL23A: 11; TRPC6: 2; TUBB8: 1 |
| GO MF ancestor ribonucleoside triphosphate phosphatase activity (`GO:0017111`) | 32 / 1,396 | 8.867448674486745 / 9 [3--15] | 23 | 3.4694612299843444 | SEPTIN14: 19; DDX11: 10; TUBB: 1; TUBB8: 1; TUBB8B: 1 |
| GO CC ancestor vesicle (`GO:0031982`) | 60 / 7,378 | 26.70570705707057 / 27 [20--34] | 33 | 2.223798112399232 | DDX11: 10; FABP5: 1; KRT18: 1; RARRES2: 2; RPL23A: 11; SEPTIN14: 19; SNX18: 2; TUBB: 1; TUBB8: 1; TUBB8B: 1; VAMP7: 1; WASHC1: 10 |

All 143 certified rows had 0 exceedances among 99,999 common placements. Thus
their plus-one empirical p-value is `0.00001`, with plus-one-transformed 95%
Clopper--Pearson interval `[0.00001, 0.00004688811415007221]` and sequential
98.75% interval `[0.00001, 0.00006075045029073397]`. Their observed burdens
range from 10 to 60 copies and their smoothed burden ratios range from
2.1113463602870506 to 13.15982453391822. These common raw-tail values do not
make the terms interchangeable: collection adjustment and maxT depend on the
complete null vector and on each row's standardized effect.

The table below gives all adjusted quantities needed to interpret the selected
rows. Brackets are sequential 98.75% Monte Carlo intervals. BH and BY were
computed within each ontology collection over direct plus ancestor rows;
collection and global columns are single-step maxT safeguards. The support rule
required both the BH upper bound and the global-maxT upper bound to be at most
0.05.

| Exact tested row | BH q [sequential bounds] | BY q | Collection maxT p [sequential bounds] | Global maxT p [sequential bounds] | Primary status |
|---|---:|---:|---:|---:|---|
| GO CC direct exocyst (`GO:0000145`) | 0.0004867272727272728 [0.00046444547474362003--0.0029568900986962697] | 0.004122509168057452 | 0.00001 [0.00001--0.00006075045029073397] | 0.00001 [0.00001--0.00006075045029073397] | `CERTIFIED_PASS` |
| GO CC direct WASH complex (`GO:0071203`) | 0.0004867272727272728 [0.00046444547474362003--0.0029568900986962697] | 0.004122509168057452 | 0.00001 [0.00001--0.00006075045029073397] | 0.00001 [0.00001--0.00006075045029073397] | `CERTIFIED_PASS` |
| GO BP direct Arp2/3 complex-mediated actin nucleation (`GO:0034314`) | 0.0013320138888888889 [0.0012953864056011006--0.00809204435435117] | 0.013904782717817914 | 0.00001 [0.00001--0.00006075045029073397] | 0.00001 [0.00001--0.00006075045029073397] | `CERTIFIED_PASS` |
| GO BP direct exocytosis (`GO:0006887`) | 0.0013320138888888889 [0.0012953864056011006--0.00809204435435117] | 0.013904782717817914 | 0.00631 [0.005702269484217936--0.006962574139855258] | 0.00658 [0.005959092930561787--0.007245716247357593] | `CERTIFIED_PASS` |
| GO BP direct endosomal transport (`GO:0016197`) | 0.0013320138888888889 [0.0012953864056011006--0.00809204435435117] | 0.013904782717817914 | 0.01626 [0.015278122834583864--0.017285639785549308] | 0.01626 [0.015278122834583864--0.017285639785549308] | `CERTIFIED_PASS` |
| GO MF direct DNA helicase activity (`GO:0003678`) | 0.001098421052631579 [0.001000040324738357--0.006672957355619042] | 0.010236614962315746 | 0.00001 [0.00001--0.00006075045029073397] | 0.00001 [0.00001--0.00006075045029073397] | `CERTIFIED_PASS` |
| Reactome direct XBP1(S) activates chaperone genes (`R-HSA-381038`) | 0.0010386666666666667 [0.0009501624586484753--0.006309946770197568] | 0.008955053211090762 | 0.00001 [0.00001--0.00006075045029073397] | 0.00001 [0.00001--0.00006075045029073397] | `CERTIFIED_PASS` |
| GO MF direct rRNA binding (`GO:0019843`) | 0.001098421052631579 [0.001000040324738357--0.006672957355619042] | 0.010236614962315746 | 0.00001 [0.00001--0.00006075045029073397] | 0.00001 [0.00001--0.00006075045029073397] | `CERTIFIED_PASS` |
| GO BP ancestor sexual reproduction (`GO:0019953`) | 0.0013320138888888889 [0.0012953864056011006--0.00809204435435117] | 0.013904782717817914 | 0.04061 [0.03906683120484941--0.04219462666710295] | 0.04071 [0.03916498950314897--0.04229645913000627] | `CERTIFIED_PASS` |
| GO MF ancestor ribonucleoside triphosphate phosphatase activity (`GO:0017111`) | 0.001098421052631579 [0.001000040324738357--0.006672957355619042] | 0.010236614962315746 | 0.02966 [0.028336437105545634--0.0310260394181524] | 0.04071 [0.03916498950314897--0.04229645913000627] | `CERTIFIED_PASS` |
| GO CC ancestor vesicle (`GO:0031982`) | 0.0004867272727272728 [0.00046444547474362003--0.0029568900986962697] | 0.004122509168057452 | 0.01975 [0.018667894338789202--0.02087552561163051] | 0.03616 [0.034701530514595555--0.03766033948239383] | `CERTIFIED_PASS` |

## Functional systems are post-inference summaries

The following organization is for biological display only. It was constructed
after every exact term burden and decision was fixed. It does not merge GO
nodes, combine collections, reduce multiplicity, or create a system-level
p-value.

1. **Endosomal trafficking, exocyst, secretion, and actin organization.**
   Forty-one supported rows have exactly the same ten WASHC1-source physical
   contributors; specific examples include exocyst, WASH complex, and
   Arp2/3-complex-mediated actin nucleation. Other separately tested rows add
   exact contributors: exocytosis adds one VAMP7 copy, and endosomal transport
   adds two SNX18-source copies. “WASH/endosomal-actin system” is therefore a
   post-inference summary of exact terms, not a tested WASH-family enrichment.

2. **DNA helicase, replication/chromosome organization, and the three
   Reactome UPR rows.** Sixty-five supported rows have exactly the same ten
   DDX11-source physical contributors. GO DNA helicase activity and Reactome
   XBP1(S) activates chaperone genes are two specific direct rows; Reactome
   IRE1alpha activates chaperones (`R-HSA-381070`) and Unfolded Protein
   Response (`R-HSA-381119`) are separate ancestor rows. This shared driver
   signature does not turn the 65 correlated rows into 65 independent
   biological systems, and it does not demonstrate that the DDX11L copies are
   expressed or active in either process.

3. **Ribosomal/nucleolar annotations.** Eleven RPL23A-source physical copies
   contribute 11 to each of the separately tested direct GO MF rows rRNA
   binding, cadherin binding, and TORC2 complex binding. Related supported GO
   BP terms include regulation of nucleolar large-rRNA transcription, but
   their exact contributors must be checked term by term rather than inferred
   from this display label.

4. **Broader reproduction, nucleotide-hydrolysis, cortex, complex, and
   vesicle ancestors.** These rows have larger mixed-source burdens. For
   example, sexual reproduction has 33 copies from four functional sources,
   ribonucleoside triphosphate phosphatase activity has 32 from five sources,
   and vesicle has 60 from 12 sources. SEPTIN14- and RPL23A-source copies are
   strong contributors to some of these rows, but no SEPTIN14, RPL23A, TUBB8,
   or other family enrichment was tested.

## Exact physical-copy examples and contributor ledgers

The examples below show why source or family collapse is forbidden. Every
listed record has `physical_copy_cn=1` in the validated
[`EXACT_TERM_CONTRIBUTORS.tsv.gz`](results/release/EXACT_TERM_CONTRIBUTORS.tsv.gz).
Coordinates in the `copy_id` strings are one-based closed intervals.

### Ten WASHC1-source copies contribute ten

The direct exocyst, WASH-complex, and Arp2/3-nucleation rows each receive one
unit from each of these ten records:

- `CHM13v2.0|chr11:5962-21167|-|LOC100996442`
- `CHM13v2.0|chr12:7478-23077|-|WASH8P`
- `CHM13v2.0|chr15:99730108-99745901|+|WASH3P`
- `CHM13v2.0|chr16:7007-11013|-|WASH4P`
- `CHM13v2.0|chr19:11934-21953|-|WASH5P`
- `CHM13v2.0|chr20:80190-92223|-|WASH5P_1`
- `CHM13v2.0|chr20:66187777-66202791|+|WASH7P`
- `CHM13v2.0|chr3:201079785-201097323|+|WASH8P_1`
- `CHM13v2.0|chr9:8314-24625|-|WASHC1`
- `CHM13v2.0|chrX:154249567-154252407|+|WASH6P`

All ten point to functional source `NCBIGene:100287171` (WASHC1); WASHC1
itself uses the exact-self route and the other nine use an explicit directed
human `Related functional gene` record. One source is still ten physical
copies. The name-defined WASH cohort contains nine of the ten because
`LOC100996442` does not begin with `WASH`; therefore the audited named-cohort
burden is 9 while the exact WASH-complex term burden is 10.

### Ten DDX11-source copies contribute ten

DNA helicase activity and the three Reactome UPR-related rows each receive one
unit from each of these ten records:

- `CHM13v2.0|chr11:3434-5924|+|DDX11L17`
- `CHM13v2.0|chr12:5295-7527|+|DDX11L8`
- `CHM13v2.0|chr15:99745854-99748389|-|DDX11L9`
- `CHM13v2.0|chr16:4481-7016|+|DDX11L10`
- `CHM13v2.0|chr19:4460-6999|+|DDX11L1`
- `CHM13v2.0|chr20:72708-75239|+|DDX11L9_1`
- `CHM13v2.0|chr20:66202741-66205279|-|DDX11L5_1`
- `CHM13v2.0|chr3:201097274-201099490|-|DDX11L8_1`
- `CHM13v2.0|chr9:5788-8326|+|DDX11L5`
- `CHM13v2.0|chrX:154252395-154254921|-|DDX11L16`

All ten point through explicit directed records to one functional source,
`NCBIGene:1663` (DDX11). Collapsing that source or the DDX11L name group would
change the tested burden from 10 to the forbidden value 1.

### Eleven RPL23A-source copies contribute eleven

The direct rRNA-binding row receives one unit from RPL23AP25 at
`chr1:248375192-248375720`, RPL23AP60 at
`chr10:134745871-134746398`, RPL23AP24 at
`chr16:96318486-96318958`, RPL23AP87 at
`chr17:84250806-84264715`, RPL23AP79 at
`chr19:61696473-61697001`, RPL23AP88 at
`chr2:242685066-242685590`, RPL23AP21 at `chr20:11293-11818`,
RPL23AP4 at `chr21:45077112-45077637`, RPL23AP84 at
`chr4:193562943-193563471`, RPL23AP45 at
`chr5:182034678-182035205`, and RPL23AP47 at
`chr9:150605300-150605772`. All eleven point to `NCBIGene:6147` (RPL23A)
through explicit directed records. One functional source therefore contributes
11 physical copies, not one source count.

## Named cohorts: exact burdens are not family hypotheses

The five prespecified name cohorts are entry audits and contributor summaries,
not ontology hypotheses. “Unresolved” in this table means that a PHR physical
copy did not emit an admissible functional source; it does not mean the copy is
biologically absent or inactive. The selected term is an exact ontology row,
and its total burden can include contributors outside the name cohort.

| Name cohort | Genome physical CN | PHR midpoint physical CN | Ontology-eligible PHR CN | PHR CN without an admissible source | Exact selected row | Total observed term CN (named-cohort contribution) | Null mean / median [95%] | BH q; global maxT p | Status |
|---|---:|---:|---:|---:|---|---:|---:|---:|---|
| DUX4/DUX4L | 107 | 68 | 65 | 3 | Reactome direct Zygotic genome activation (`R-HSA-9819196`) | 65 (65) | 14.044820448204481 / 7 [0--54] | 0.575373023255814; 0.98807 | `CERTIFIED_NONPASS` |
| DDX11L | 12 | 10 | 10 | 0 | GO MF direct DNA helicase activity (`GO:0003678`) | 10 (10) | 0.2978829788297883 / 0 [0--2] | 0.001098421052631579; 0.00001 | `CERTIFIED_PASS` |
| TUBB8 | 16 | 7 | 2 | 5 | GO BP ancestor organelle organization (`GO:0006996`) | 28 (2) | 14.00429004290043 / 14 [9--20] | 0.0013320138888888889; 0.59552 | `CERTIFIED_NONPASS` |
| OR4F | 15 | 11 | 4 | 7 | Reactome direct Expression and translocation of olfactory receptors (`R-HSA-9752946`) | 7 (4) | 3.2308523085230854 / 3 [0--7] | 1.0; 0.99992 | `CERTIFIED_NONPASS` |
| WASH | 18 | 9 | 9 | 0 | GO CC direct WASH complex (`GO:0071203`) | 10 (9) | 0.4969349693496935 / 0 [0--2] | 0.0004867272727272728; 0.00001 | `CERTIFIED_PASS` |

The DUX term has 750/99,999 raw exceedances, plus-one p=0.00751, 95% interval
[0.006984582162182748, 0.008064278254190906], and sequential interval
[0.006845701839263833, 0.008218992946982594]. Its smoothed burden ratio is
4.50332131862692, but its adjusted safeguards do not pass. The OR4F term has
3,424/99,999 exceedances, p=0.03425, 95% interval
[0.0331314589032472, 0.035395660083812806], sequential interval
[0.032829672112938864, 0.03571237546393364], and smoothed ratio
2.010264513249786; it also does not pass. The TUBB8 example shows the converse
adjustment pattern: its smoothed burden ratio is 1.9649358855692631 and its raw
p and BH q are small, but its global maxT upper bound is
0.599397102696665, so it is non-passing. These rows remain exact term-level
results and cannot be promoted to DUX, TUBB8, or OR4F family enrichment.

The three Reactome UPR rows happen to share the ten DDX11L contributor copies,
and DDX11L has 185 source-term rows that carry all ten copies. DUX has 83
source-term rows that carry all 65 eligible copies. TUBB8's two eligible
copies have distinct exact sources and no source term carries both; OR4F's
four eligible copies have four exact sources and no source term carries all
four. WASH has 134 WASHC1-source term rows that carry all nine name-defined
copies. These contributor patterns explain term burdens after inference; none
is a tested family-level result.

## Mapping coverage and unresolved physical copies

Annotation absence is never interpreted as biological absence. V6 retained
all physical records and failed closed when it could not establish an exact
term-bearing source.

| Scope | Total physical CN | Exact-self eligible CN | Explicit-related eligible CN | Total ontology-eligible CN | Ontology-ineligible CN | Eligible fraction |
|---|---:|---:|---:|---:|---:|---:|
| Genome | 61,312 | 20,405 | 11,561 | 31,966 | 29,346 | 52.137% |
| PHR midpoint, primary | 402 | 20 | 167 | 187 | 215 | 46.517% |
| PHR any overlap | 412 | 23 | 170 | 193 | 219 | 46.845% |

At the primary PHR midpoint, the 215 copies without a functional source
comprise 200 `UNSUPPORTED_FAIL_CLOSED`, 8 `AMBIGUOUS_FAIL_CLOSED`, 5
`TYPE_ONLY`, and 2 `UNRESOLVED` records. Genome-wide, the corresponding counts
are 28,407, 45, 890, and 4. These rows remain in coverage denominators but emit
no ontology edge. The route and biotype breakdown is retained in
[`MAPPING_COVERAGE.tsv`](results/release/MAPPING_COVERAGE.tsv).

The 15 name-cohort PHR copies without an admissible source are also retained
explicitly:

- DUX4/DUX4L: DUX4L28 and DUX4L29 are `UNRESOLVED` because no frozen
  admissible copy-specific source relation exists; DBET is `TYPE_ONLY` because
  alias-only term transfer was rejected.
- TUBB8: TUBB8P7 and four TUBB8P8 records (`TUBB8P8`, `TUBB8P8_1`,
  `TUBB8P8_2`, `TUBB8P8_3`) are `AMBIGUOUS_FAIL_CLOSED` because the sequence
  best hit competes with nomenclature.
- OR4F: OR4F4 and OR4F4_1 are `AMBIGUOUS_FAIL_CLOSED` because of a GFF/HGNC
  biotype conflict; OR4F8BP is `AMBIGUOUS_FAIL_CLOSED` after resolved sequence
  competition; OR4F2P, OR4F28P, OR4F8P, and OR4F7P are `TYPE_ONLY` because of
  near-tied copy-specific candidates.

Their exclusion means only that V6 has no admissible edge for them. It is not
evidence that DUX, TUBB8, OR4F, or any other biology is absent from the PHRs.

## Spatial null, uncertainty, and sensitivity

The primary null translated 37 complete rigid PHR blocks uniformly among
valid integer starts on the same arm and within the same prespecified
terminal-distance stratum. It preserved block widths, internal geometry,
copy clusters, source/term edges, and physical-CN weights. All hypotheses used
the same 99,999 PCG64DXSM joint placements (seed `2026071301`, spawn key
`[0]`). Candidate spaces ranged from 250,000 to 497,500 starts. Midpoint was
primary; any positive base overlap on the identical placements was the paired
boundary sensitivity.

The any-overlap sensitivity remained directionally elevated for the displayed
terms, but it is not a second promotion rule. For example, WASH complex and
exocyst each had 10 copies against null mean 0.6440164401644016 and median 0
[0--2], while DNA helicase activity had 10 against 0.32677326773267734 and 0
[0--2]. The complete overlap rows and their sensitivity-only maxT values are
in `TERM_RESULTS.tsv.gz` under `assignment=overlap` and
`mc_status=SENSITIVITY_COMPLETE`.

The complete 99,999-placement screen resolved the prespecified sequential
decision for every primary row. No selective extension candidate existed and
no million-placement reflexive rerun was triggered. Zero unresolved rows means
there is no licensed “near-significant” category in V6.

## Limitations and prohibited interpretations

- V6 is a regional analysis of one CHM13 assembly. It does not estimate
  population prevalence, pangenome copy-number frequency, or individual-level
  enrichment.
- Only 187 of 402 primary PHR-midpoint copies are ontology eligible. The 215
  ineligible copies are retained physical observations; missing source mapping
  is not negative evidence about biological function.
- A directed `Related functional gene` record licenses a frozen annotation
  edge for this analysis. It does not establish functional equivalence between
  the coordinate copy and source, nor expression or activity of the copy.
- GO/Reactome terms are nested and many supported rows share identical
  contributors. There are 143 supported term rows, not 143 demonstrably
  independent functional systems.
- “WASH/endosomal-actin,” “DDX11/helicase-UPR,” “ribosomal/nucleolar,” and the
  broader display labels above are post-inference summaries without their own
  test statistic or adjusted p-value.
- WASHC1-, DDX11-, RPL23A-, SEPTIN14-, DUX-, TUBB8-, and OR4F-associated
  contributors are driver descriptions, not gene-family enrichment results.
- Physical annotation burden is distinct from transcription, translation,
  protein dosage, molecular activity, chromosome contact, and causal
  coordination.
- A small raw or BH value cannot override the global maxT safeguard, as the
  TUBB8 organelle-organization example shows. Conversely, contributor
  inspection cannot rescue a non-passing term or demote a certified term.

## Reproducibility and claim-source index

Only artifacts covered by the independently PASSed V6 release are cited as
evidence here:

- Independent verdict and reconstructed checks:
  [`V6_VALIDATION_REPORT.md`](V6_VALIDATION_REPORT.md) and
  [`V6_VALIDATION.json`](V6_VALIDATION.json).
- Full midpoint and overlap statistics, effect sizes, uncertainty, BH/BY, and
  maxT values:
  [`results/release/TERM_RESULTS.tsv.gz`](results/release/TERM_RESULTS.tsv.gz).
- Exact copy IDs, coordinates, CN, term relation, source assignments, evidence
  record IDs, and midpoint/overlap membership:
  [`results/release/EXACT_TERM_CONTRIBUTORS.tsv.gz`](results/release/EXACT_TERM_CONTRIBUTORS.tsv.gz).
- Mapping coverage by scope, route, and biotype:
  [`results/release/MAPPING_COVERAGE.tsv`](results/release/MAPPING_COVERAGE.tsv).
- Named-cohort entry audit:
  [`results/release/NAMED_COHORT_INFERENCE_AUDIT.tsv`](results/release/NAMED_COHORT_INFERENCE_AUDIT.tsv)
  together with the independently reconstructed source-map audit
  [`NAMED_COHORT_AUDIT.tsv`](NAMED_COHORT_AUDIT.tsv).
- Frozen hypothesis families and sizes:
  [`results/release/MULTIPLICITY_FAMILIES.tsv`](results/release/MULTIPLICITY_FAMILIES.tsv).
- Placement geometry and randomization provenance:
  [`V6_RUN_PROTOCOL.md`](V6_RUN_PROTOCOL.md),
  [`results/release/NULL_PLACEMENT_MANIFEST.tsv`](results/release/NULL_PLACEMENT_MANIFEST.tsv),
  and
  [`results/release/COMPUTATIONAL_PROVENANCE.json`](results/release/COMPUTATIONAL_PROVENANCE.json).
- Exact released-file hashes:
  [`results/release/RELEASE_SHA256.tsv`](results/release/RELEASE_SHA256.tsv).

The candidate manuscript prose is in
[`MANUSCRIPT_CANDIDATE_TEXT.md`](MANUSCRIPT_CANDIDATE_TEXT.md). It is a handoff
artifact only; `submission/paper.tex` was not edited.
