# V6 pre-run independent gate

## Decision

**PASS.** The exact frozen ontology_v6 physical-copy map is authorized as an input to the downstream enrichment task. The authorization is digest-bound: any input or audited release byte change requires this gate to be rerun. No enrichment was run here.

The gate independently recovered **61,312** raw physical GFF3 loci, **402** PHR-midpoint memberships, **412** any-overlap memberships, **1,686,727** source-term closure rows, and **2,929,709** physical-copy term edges. It did not import `build_genomewide_source_map.py` or `check_genomewide_source_map.py`.

## Evidence-before-arithmetic contract

Source assignment was adjudicated before the PHR BED was opened and before copy-term arithmetic. An ontology source is admissible only when it is a term-bearing exact self identifier or the target of the exact frozen directed human NCBI Gene `Related functional gene` pair for that copy's own identifier. Names, aliases, family membership, reverse-only relations, and one-way sequence similarity do not emit terms. Ambiguous and termless-self rows fail closed.

## Copy-number retention

Source assertions are deduplicated only within a source. Expansion is one row per physical copy and source term. Each burden was independently recomputed as the sum of `physical_copy_cn` (and separately `phr_midpoint_cn`), and contributor-row counts were required to equal those sums. Repeated sources were never reduced to one contributor.

| Cohort | Genome physical copies | PHR physical copies | Reviewed PHR ontology contributors / term burden |
|---|---:|---:|---:|
| DUX4_DUX4L | 107 | 68 | 65 |
| DDX11L | 12 | 10 | 10 |
| TUBB8 | 16 | 7 | 2 |
| OR4F | 15 | 11 | 4 |
| WASH | 18 | 9 | 9 |

DUX4 contributes 65 to every one of its direct and ancestor terms; DDX11 contributes 10 and WASHC1 contributes 9. TUBB8's two reviewed PHR contributors and OR4F's four contributors retain their separate exact functional sources; their copies are not collapsed or combined through a family union.

## Mapping coverage

| Scope | Physical CN | Own exact identity CN | Ontology contributor CN | Ineligible CN | Direct edge CN | Closure edge CN |
|---|---:|---:|---:|---:|---:|---:|
| genome | 61312 | 60430 | 31966 | 29346 | 626048 | 2929709 |
| phr_midpoint | 402 | 402 | 187 | 215 | 3633 | 16763 |
| phr_any_overlap | 412 | 412 | 193 | 219 | 3831 | 17579 |

[`PRE_RUN_V6_MAPPING_COVERAGE.tsv`](PRE_RUN_V6_MAPPING_COVERAGE.tsv) reports genome, PHR-midpoint, and PHR-any-overlap coverage as totals and as complete evidence-route, biotype, and route-by-biotype partitions.

## Machine gate

- [`PRE_RUN_V6_GATE.json`](PRE_RUN_V6_GATE.json) is the strict PASS/BLOCK decision and authorization record.
- [`PRE_RUN_V6_GATE_EVIDENCE.tsv`](PRE_RUN_V6_GATE_EVIDENCE.tsv) records every hard check, observed value, expected value, and evidence source.
- [`PRE_RUN_V6_AUDITED_RELEASE.sha256.tsv`](PRE_RUN_V6_AUDITED_RELEASE.sha256.tsv) binds PASS to the exact production artifacts audited.
- [`independent_gate_genomewide_source_map.py`](independent_gate_genomewide_source_map.py) is the standalone reconstruction; its tests include explicit prohibited-route and copy-collapse controls.

Downstream enrichment must require `status == PASS`, `enrichment_authorized == true`, and unchanged audited digests. A BLOCK result or any digest drift denies inference.

Validation status: **PASS** (118 recorded hard checks).
