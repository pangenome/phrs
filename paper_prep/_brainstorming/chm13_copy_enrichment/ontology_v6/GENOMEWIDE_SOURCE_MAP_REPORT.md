# Genome-wide CHM13 copy-to-ontology source map

## Outcome

V6 reconciles **61,312 coordinate-distinct CHM13v2.0 GFF3 gene records** one-for-one. Every row has `physical_copy_cn=1`; coordinate-coincident records remain separate. Locus annotation identity, functional source identity, and ontology-term eligibility are separate fields.

The target-blind assignment and copy-source evidence tables were written and SHA-256 frozen before the PHR BED was opened. The PHR join is therefore an audit/weight layer, not evidence used to select a source.

No enrichment, hypothesis test, p-value, or biological significance claim is produced by this release.

## Source contract

An exact raw stable identifier may carry only its own frozen GO/Reactome terms. A non-self copy may inherit only through an exact directed human NCBI Gene `Related functional gene` row, or through a copy-specific disposition already frozen by the remediated pilot. Family roots, names, aliases, and one-way sequence similarity never authorize inheritance. Ambiguous and unsupported copies remain in every coverage denominator and emit no term edges.

## Coverage

| Scope | Physical CN | Own annotation identity CN | Ontology contributors | Ineligible CN | Copy-term edges (closure) |
|---|---:|---:|---:|---:|---:|
| genome | 61312 | 60430 | 31966 | 29346 | 2929709 |
| phr_midpoint | 402 | 402 | 187 | 215 | 16763 |
| phr_any_overlap | 412 | 412 | 193 | 219 | 17579 |

## Named-cohort audit

The cohort definitions are evaluated only after the source freeze.

| Cohort | Genome copies | PHR copies | PHR term contributors | Sources | Status |
|---|---:|---:|---:|---|---|
| DUX4_DUX4L | 107 | 68 | 65 | DUX4 | PASS |
| DDX11L | 12 | 10 | 10 | DDX11 | PASS |
| TUBB8 | 16 | 7 | 2 | TUBB8|TUBB8B | PASS |
| OR4F | 15 | 11 | 4 | OR4F17|OR4F29|OR4F3|OR4F5 | PASS |
| WASH | 18 | 9 | 9 | WASHC1 | PASS |

DUX4/DUX4L has 68 PHR physical copies and 65 defensible DUX4-source contributors. DDX11L is 10/10; TUBB8 is 2/7; OR4F is 4/11; and WASH is 9/9. `NAMED_COHORT_TERM_BURDENS.tsv.gz` exposes the exact direct and ancestor term burden behind these contributor counts.

## Copy-number semantics

`PHYSICAL_COPY_TERM_EDGES.tsv.gz` retains `physical_copy_cn`, `phr_midpoint_cn`, and `phr_any_overlap_cn` on every edge. Source-level assertions are deduplicated once in `SOURCE_TERMS.tsv.gz`, but propagation creates one edge for every eligible physical copy. Thus N copies assigned to a source contribute N to each of that source's direct terms and to each true-path ancestor.

## Source dispositions

- `AMBIGUOUS_FAIL_CLOSED`: 45 physical copies
- `EXACT_SELF`: 20,405 physical copies
- `EXPLICIT_RELATED_FUNCTIONAL_GENE`: 11,561 physical copies
- `TYPE_ONLY`: 890 physical copies
- `UNRESOLVED`: 4 physical copies
- `UNSUPPORTED_FAIL_CLOSED`: 28,407 physical copies

## Reproduction and release files

Run from the repository root:

```bash
python3 paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v6/build_genomewide_source_map.py
python3 paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v6/check_genomewide_source_map.py
python3 -m unittest paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v6/test_genomewide_source_map.py
```

`INPUT_MANIFEST.tsv` pins every upstream byte stream. `SOURCE_ASSIGNMENT_FREEZE.sha256.tsv` records the pre-target map/evidence digest. `OUTPUT_MANIFEST.sha256.tsv` pins the release tables, code, report, and validation JSON. All large release tables are deterministic gzip files with gzip modification time zero.

Validation status: **PASS** (21 checks).
