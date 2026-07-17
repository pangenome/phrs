# Final V7 whole-genome non-PHR physical-copy ontology report

## Release verdict

**V7 validation status: PASS (18/18 checks).**

The corrected V7 analysis tested all **31,235** frozen exact direct/ancestor
GO and Reactome hypotheses against the complete ontology-eligible CHM13
non-PHR genome. The primary midpoint universe contains **187 PHR copies**
and **31,779 non-PHR copies** (N = 31,966). The paired any-overlap
sensitivity contains 193 PHR and 31,773 non-PHR copies. No background was
sampled, no PHR interval was moved, and no spatial control defines V7 inference.

The prespecified primary rule (within-collection BH q <= 0.05 and global
Holm-adjusted p <= 0.05) supports **209 exact midpoint term rows**. The paired
overlap sensitivity supports **211 rows**. Primary supports by collection are
`{"GO_BP": 120, "GO_CC": 37, "GO_MF": 45, "Reactome": 7}`; by relation they are `{"ancestor": 119, "direct": 90}`. All 29779 midpoint rows with zero PHR burden
were retained and tested.

## What was counted and tested

For each frozen `(collection, relation, term_id)` row, the result table gives
`a_T`, `b_T`, `c_T`, and `d_T` over coordinate-distinct physical copies.
The p-value is the exact upper tail for observing at least `a_T` term-bearing
copies among K PHR copies in the complete fixed population of N eligible
copies containing M_T term-bearing copies. This hypergeometric calculation is
a complete finite-population copy-burden test, not generic weak gene-list ORA.

Every inherited V6 copy has CN=1. Multiple DUX4/DUX4L, DDX11L, WASH,
RPL23A, OR4F, and TUBB8 coordinates therefore contribute one unit each even
when they share a functional source. The named-cohort audit explicitly shows
the larger coordinate-copy counts beside the forbidden source-collapsed counts.

## Leading exact supported rows

| Collection | Relation | Exact term | a / K | c / non-PHR | Fold | Exact p | BH q | Global Holm |
|---|---|---|---:|---:|---:|---:|---:|---:|
| GO_BP | direct | negative regulation of G0 to G1 transition (`GO:0070317`) | 65 / 187 | 14 / 31779 | 789.0126050420168 | 2.87e-136 | 5.5e-132 | 8.96e-132 |
| GO_BP | ancestor | regulation of G0 to G1 transition (`GO:0070316`) | 65 / 187 | 18 / 31779 | 613.6764705882354 | 1.7e-133 | 1.63e-129 | 5.3e-129 |
| Reactome | direct | Zygotic genome activation (ZGA) (`R-HSA-9819196`) | 65 / 187 | 50 / 31779 | 220.9235294117647 | 2.35e-118 | 7.31e-115 | 7.33e-114 |
| Reactome | ancestor | Maternal to zygotic transition (MZT) (`R-HSA-9816359`) | 65 / 187 | 271 / 31779 | 40.76079878445843 | 2.27e-81 | 3.54e-78 | 7.09e-77 |
| GO_MF | direct | RNA polymerase II transcription regulatory region sequence-specific DNA binding (`GO:0000977`) | 74 / 187 | 512 / 31779 | 24.561810661764707 | 8.45e-79 | 5.29e-75 | 2.64e-74 |
| GO_MF | ancestor | DNA-binding transcription activator activity (`GO:0001216`) | 74 / 187 | 553 / 31779 | 22.7407722582704 | 1.51e-76 | 4.73e-73 | 4.72e-72 |
| GO_MF | direct | DNA-binding transcription activator activity, RNA polymerase II-specific (`GO:0001228`) | 74 / 187 | 570 / 31779 | 22.062538699690403 | 1.16e-75 | 2.43e-72 | 3.64e-71 |
| GO_CC | direct | nuclear membrane (`GO:0031965`) | 65 / 187 | 349 / 31779 | 31.650935445811562 | 4.95e-75 | 1.32e-71 | 1.54e-70 |
| GO_CC | ancestor | nuclear envelope (`GO:0005635`) | 66 / 187 | 386 / 31779 | 29.057299603779338 | 4.7e-74 | 6.28e-71 | 1.47e-69 |
| GO_BP | ancestor | negative regulation of cell cycle process (`GO:0010948`) | 65 / 187 | 439 / 31779 | 25.162133190406003 | 3.38e-69 | 2.16e-65 | 1.06e-64 |
| GO_BP | direct | apoptotic process (`GO:0006915`) | 65 / 187 | 453 / 31779 | 24.384495520062327 | 2.15e-68 | 1.03e-64 | 6.71e-64 |
| GO_BP | ancestor | negative regulation of cell cycle (`GO:0045786`) | 65 / 187 | 459 / 31779 | 24.06574394463668 | 4.67e-68 | 1.79e-64 | 1.46e-63 |
| GO_MF | direct | transcription cis-regulatory region binding (`GO:0000976`) | 65 / 187 | 487 / 31779 | 22.682087208600073 | 1.54e-66 | 2.41e-63 | 4.8e-62 |
| GO_BP | ancestor | positive regulation of DNA-templated transcription (`GO:0045893`) | 84 / 187 | 1349 / 31779 | 10.58195613308333 | 2.68e-62 | 8.56e-59 | 8.36e-58 |
| GO_MF | ancestor | DNA-binding transcription factor activity (`GO:0003700`) | 75 / 187 | 974 / 31779 | 13.085819543423119 | 3.49e-61 | 4.37e-58 | 1.09e-56 |
| GO_BP | direct | negative regulation of cell population proliferation (`GO:0008285`) | 65 / 187 | 607 / 31779 | 18.19798430080434 | 7.16e-61 | 1.96e-57 | 2.23e-56 |
| GO_BP | ancestor | regulation of cell cycle process (`GO:0010564`) | 75 / 187 | 1099 / 31779 | 11.597441524380454 | 1.39e-57 | 3.34e-54 | 4.34e-53 |
| GO_MF | direct | sequence-specific double-stranded DNA binding (`GO:1990837`) | 66 / 187 | 790 / 31779 | 14.197617274758006 | 2.28e-55 | 2.38e-52 | 7.13e-51 |
| GO_MF | direct | DNA-binding transcription factor activity, RNA polymerase II-specific (`GO:0000981`) | 83 / 187 | 1625 / 31779 | 8.680072398190045 | 6.09e-55 | 5.45e-52 | 1.9e-50 |
| GO_BP | ancestor | regulation of cell cycle (`GO:0051726`) | 75 / 187 | 1296 / 31779 | 9.834558823529411 | 1.1e-52 | 2.34e-49 | 3.43e-48 |

The table above is ordered by global Holm value and copy excess; the complete
results, including nonsupports and zero-PHR rows, are in `TERM_RESULTS.tsv.gz`.

## Validation gates

- **inherited_v6_frozen_objects:** `PASS`.
- **all_frozen_hypotheses_tested_twice:** `PASS`.
- **paired_assignments_complete:** `PASS`.
- **exact_frozen_hypothesis_keys:** `PASS`.
- **zero_phr_burdens_retained:** `PASS`.
- **all_term_2x2_tables_partition:** `PASS`.
- **eligible_universe_constant:** `PASS`.
- **phr_nonphr_exact_complement:** `PASS`.
- **no_phr_contributor_marked_background:** `PASS`.
- **exact_contributor_recount_matches_every_a_T:** `PASS`.
- **coordinate_copies_not_source_collapsed:** `PASS`.
- **all_named_cohort_audits_pass:** `PASS`.
- **no_stochastic_or_placement_runtime_path:** `PASS`.
- **multiplicity_complete:** `PASS`.
- **prespecified_primary_rule_recomputed:** `PASS`.
- **community_summary_post_inference_only:** `PASS`.
- **post_inference_community_sums_match_supported_terms:** `PASS`.
- **contributor_ledger_inherited_and_statused:** `PASS`.

The exact contributor ledger reuses the independently validated V6 any-overlap
contributor rows and adds V7 midpoint/overlap partition and decision statuses.
`MAPPING_COVERAGE.tsv` proves the eligible partition and empty intersection.
`COMMUNITY_TERM_SUMMARY.tsv` was generated only after all exact term decisions
were fixed; it is a descriptive source/community display and created no test.

## Interpretation boundary

V7 answers which exact ontology terms have excess annotation-bearing physical
copy burden in CHM13 PHRs relative to the rest of the eligible CHM13 genome.
It does not test expression, protein activity, dosage effect, retained function
of pseudogenes, biological independence of adjacent copies, or population prevalence.

This V7 result supersedes V6's placement-null inference for the biological
question specified here; the V6 map and hypothesis catalog remain the frozen inputs.
