# V6 independent validation report

**Strict result: PASS.** This is not a qualified pass. The copy-number, mapping, ontology, geometry, Monte Carlo, multiplicity, and provenance gates all had to pass; a failure in physical-copy semantics would force the machine result to `FAIL`.

## Scope and independence

The validator in `independent_validation/validate_v6_independent.py` imports none of the production mapping, sampler, inference, or prior validation modules. It rebuilt the matrix from the raw CHM13 GFF, frozen per-copy evidence, raw NCBI gene2go, GO OBO, Reactome pathway relations/all-level assignments, cytobands, and the PHR BED. NumPy was used only to reproduce the specified PCG64DXSM stream and vectorize interval recounts; Clopper--Pearson bounds were recomputed independently with R `qbeta`.

## Fail-closed source and copy gates

All 61,312 raw gene features form an exact bijection with 61,312 source assignments and 61,312 evidence rows. All coordinate-coincident records were retained (12 extra physical rows beyond coordinate+strand uniqueness), every row has positive physical CN, and summed CN is 61,312. The audit reconstructed every assignment and evidence-row digest.

Exactly 31,966 physical copies are ontology-eligible: 20,405 exact-self mappings and 11,561 mappings backed by a directed human `Related functional gene` record. All 45 ambiguous and 28,407 unsupported rows emit no functional source. The 168 remediated pilot records match their frozen upstream evidence and hashes. No name, alias, family, or sequence-only route entered ontology inference.

## Direct/ancestor terms and physical-CN recount

The raw ontology reconstruction matched all 1,686,727 `SOURCE_TERMS` rows field-for-field, including direct/ancestor relation, minimum distance, direct leaves, evidence codes, qualifiers, record identifiers, release labels, and ontology digests. Joining those terms to coordinate-distinct copies reconstructed 2,929,709 genome edge-CN (626,048 direct), 16,763 midpoint edge-CN, and 17,579 any-overlap edge-CN.

Every one of the 31,235 frozen hypotheses matched its independently rebuilt genome CN, arm count, source count, target-blind order, and multiplicity family. Both released observed columns (62,470 rows) equal sums of physical CN; all 17,579 exact contributor rows were reconstructed from raw coordinates, source evidence, and ontology closure.

### Named cohorts

| Cohort | Genome physical CN | PHR midpoint CN | Inference-eligible PHR CN |
|---|---:|---:|---:|
| DUX4_DUX4L | 107 | 68 | 65 |
| DDX11L | 12 | 10 | 10 |
| TUBB8 | 16 | 7 | 2 |
| OR4F | 15 | 11 | 4 |
| WASH | 18 | 9 | 9 |

These are physical-copy burdens, not unique-gene or family counts. The adversarial fixture placed seven CN=1 coordinates plus one CN=3 coordinate on one source/gene/family: the required burden is 10, while every source/gene/family collapse incorrectly returns 1. The validator rejected all three collapsed estimands.

## Null placements and representative retained counts

All 99,999 joint placements were regenerated from seed 2026071301 and spawn key `[0]`. All 100/100 deterministic gzip batch hashes and all 100/100 canonical-coordinate hashes matched. Each of 37 rigid blocks remained on its source arm in the same terminal-distance stratum; component widths and CN/source/term edges remained fixed.

The transient full count matrices are not committed, so the audit regenerated exact columns for leading, negative, and every named-system representative. Released mean, quantiles, maximum, exceedances, and standardized observed value matched for midpoint and any-overlap:

| Role/cohort | Term | Midpoint observed | Midpoint null mean | Midpoint exceedances | Overlap status |
|---|---|---:|---:|---:|---|
| leading_exocyst | `direct|GO:0000145` | 10 | 0.496935 | 0/99,999 | PASS |
| negative_metal_binding | `ancestor|GO:0046872` | 1 | 9.091521 | 99,998/99,999 | PASS |
| DUX4_DUX4L | `direct|R-HSA-9819196` | 65 | 14.044820 | 750/99,999 | PASS |
| DDX11L | `direct|GO:0003678` | 10 | 0.297883 | 0/99,999 | PASS |
| TUBB8 | `ancestor|GO:0006996` | 28 | 14.004290 | 0/99,999 | PASS |
| OR4F | `direct|R-HSA-9752946` | 7 | 3.230852 | 3,424/99,999 | PASS |
| WASH | `direct|GO:0071203` | 10 | 0.496935 | 0/99,999 | PASS |

## Empirical inference and multiplicity

For every result row, the audit recalculated `(exceedances + 1)/(permutations + 1)` for raw, collection-maxT, and global-maxT p-values; two-sided 95% and sequential 98.75% plus-one-transformed Clopper--Pearson intervals; BH and BY within each collection across direct plus ancestor rows; BH-adjusted sequential endpoints; count differences; and burden ratios. MaxT exceedance nesting and monotonicity in standardized observed burden also hold.

The primary classification independently returns 143 `CERTIFIED_PASS`, 31,092 `CERTIFIED_NONPASS`, and 0 `MC_UNRESOLVED` rows. No selective extension candidate exists; the complete 99,999-placement screen is sufficient under the frozen stopping rule.

## Claim classification

### Supported

There are 143 supported term-level CHM13 regional physical-copy enrichments. Support attaches only to the exact `(collection, relation, term)` rows whose sequential BH and global-maxT upper bounds are at most 0.05. A supported term may have DDX11L- or WASH-associated contributors, but this does not create a family hypothesis.

### Unresolved

None. Zero primary rows have a sequential overall decision that straddles 0.05.

### Descriptive

The named burdens 65/10/2/4/9 are validated coordinate counts. DUX zygotic-genome-activation and OR4F receptor-expression examples are descriptive because their exact rows fail complete multiplicity/maxT safeguards. TUBB8 is likewise not promoted to a family-level inference. DDX11L helicase and WASH-complex term rows are certified, but only at their frozen term keys.

### Unsupported

- Any unique-gene, unique-source, or gene-family collapse as the tested statistic.
- A DUX, DDX11L, TUBB8, OR4F, or WASH family-level hypothesis not present in the frozen catalog.
- Causal functional coordination, expression, or chromosome-contact mechanisms.
- Extrapolation from CHM13 regional copies to population prevalence or the human pangenome.
- Terms or biological systems selected after inspecting PHR membership.

## Computational provenance

Live `sacct` reconfirmed jobs 1761733, 1761734, and 1761735 as `COMPLETED` with exit code `0:0`; recorded stdout/stderr byte counts and SHA-256 digests match. The initial arrays resolve to engine `2c57923dbe1201bc472f71395e35f84442a943ddb6ba1ac5fa8d7db3825dc37a` at Git commit `63dd2d751a055b860226754158fab220e07663ba`. The released finalization resolves to engine `ecf1c26205d6f28d2f87aa5939b7eadec149e4e26db50ca8f27d30f0168eec82` at Git commit `689a51db0cc4179fdca664a7c9a5e10734792303`. All frozen raw inputs and all released artifacts match their checksum manifests.

The 200 transient count-array payloads are intentionally absent from the checkout. Their inventory, shapes, assignments, statistics, and hashes are retained; this audit therefore reconstructed representative columns from the fully regenerated placement stream instead of claiming to re-hash absent payloads. This does not qualify the strict result because every placement batch, every observed term, all inferential arithmetic, and the prespecified representative null columns passed independently.

## Machine verdict

`V6_VALIDATION.json`: `status=PASS`, `strict_pass=true`, `copy_number_semantics_pass=true`, `checks=46/46`.
