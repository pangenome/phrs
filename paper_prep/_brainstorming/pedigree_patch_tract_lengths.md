# WashU Pedigree Patch Tract-Length Summary

Primary input table:
`/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv`

Manuscript aggregate table cross-check:
`/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/all_pedigrees_patches.tsv`

Length field used: `patch_size` (assembly-derived child-flank interval length
from the WashU `odgi untangle` recombination patch table).

High-quality filter reproduced from `submission/paper.tex`: interchromosomal
patches with `min_score >= 0.8` and `500 <= patch_size <= 100000`.

Resolution audit:

- The committed WashU untangle BEDs are
  `PAN027_vs_PAN010.e50000.m1000.bed.gz`,
  `PAN027_vs_PAN011.e50000.m1000.bed.gz` and
  `PAN028_vs_PAN027.e50000.m1000.bed.gz`.
- The original `odgi untangle` command used `-e 50000 -m 1000`.
- The WashU PGGB graph was induced with wfmash `segment-length: 1000`; the
  PGGB log shows `wfmash -s 1000 -l 3000`. This is a graph-construction
  seed scale, not a hard lower bound on alignment or tract lengths.
- The observed minimum high-quality `patch_size` is 1001 bp.
  Therefore the current high-confidence patch table, which was produced from
  `odgi untangle -m 1000`, is effectively left-truncated at about 1 kb for
  tract-length comparison.
- A lower-merge rerun from the same graph is feasible and was run with
  `odgi untangle -e 50000 -m 0 -j 0 -n 1` for the three assayed transmissions.
  The companion parser merges adjacent best-hit intervals assigned to the same
  donor chrarm/haplotype and summarizes interchromosomal merged intervals with
  `min_score >= 0.8` and no 500 bp lower-size cutoff. This is a graph/untangle-
  resolved interval analysis, not the exact high-quality m1000 patch table used
  for the 494/538 within-community statistic.

Primary denominator: N = 538 high-quality WashU interchromosomal
candidate patches.

Aggregate table cross-check: `ds == "WashU"` contains 538 rows.

Overall length distribution: min 1001 bp; Q1 1136
bp; median 1517 bp; Q3 2981 bp; Q90
8247.4 bp; Q95 15448.3 bp; max 77700 bp.

Comparison ranges:

- Short primate NCO/gene-conversion tract means cited in the manuscript:
  22--95 bp.
- Longer primate CO-associated tract means cited in the manuscript:
  318--688 bp.
- "Near" windows are descriptive two-fold windows around those cited ranges:
  11--190 bp and 159--1376 bp. They are not event validation criteria.

Current-table compatibility counts:

- In 22--95 bp conversion-like range:
  0/538
  (0.0%; below current-table
  resolution).
- Near conversion-like two-fold window:
  0/538
  (0.0%; below
  current-table resolution).
- In 318--688 bp CO-associated range:
  0/538
  (0.0%; below current-table
  resolution).
- Near CO-associated two-fold window:
  228/538
  (42.4%; overlaps the
  observed lower edge but is still descriptive only).

Interpretation: the high-quality WashU candidate patch intervals in the current
manuscript tables are not resolved in the 22--95 bp or 318--688 bp cited tract
ranges. The apparent zero counts in those bins are therefore a resolution limit
of the existing graph/untangle/patch table, not evidence that such tract lengths
are biologically absent. This analysis supports only a current-table
compatibility/proportion statement about assembly-derived patch sizes; it does
not validate event-level conversion or crossover mechanisms.

## Summary Table

Full tabular output: `scripts/pedigree/patch_tract_length_summary.tsv`.
Lower-merge tabular output, when generated:
`scripts/pedigree/patch_tract_lower_merge_summary.tsv`.

| group_type | group_value | n | min | q25 | median | q75 | q95 | max | conv 22-95 | CO 318-688 | near CO 159-1376 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| overall | all_high_quality_interchromosomal | 538 | 1001 | 1136 | 1517 | 2981 | 15448.3 | 77700 | 0 | 0 | 228 |
| pattern | acros_like | 262 | 1001 | 1136 | 1552 | 3963 | 20224.7 | 77700 | 0 | 0 | 104 |
| pattern | gene_conversion_like | 136 | 1001 | 1136 | 1455.5 | 2461 | 6305.8 | 15448 | 0 | 0 | 54 |
| pattern | sandwich_same_hap | 120 | 1001 | 1136 | 1258 | 2589 | 7348.8 | 28349 | 0 | 0 | 62 |
| pattern | crossover_like | 18 | 1002 | 1146 | 1447 | 2696.2 | 19350 | 27969 | 0 | 0 | 8 |
| pattern | complex | 2 | 6288 | 6973 | 7658 | 8343 | 8891 | 9028 | 0 | 0 | 0 |
| community_status | within_community | 494 | 1001 | 1136 | 1449 | 2690 | 15249.2 | 77700 | 0 | 0 | 218 |
| community_status | cross_community | 44 | 1001 | 1876.8 | 3074.5 | 7100.5 | 19740.2 | 26882 | 0 | 0 | 10 |
| overlaps_phr | True | 506 | 1001 | 1136 | 1550 | 3131.2 | 15445.8 | 77700 | 0 | 0 | 203 |
| overlaps_phr | False | 32 | 1001 | 1001 | 1001 | 1231 | 13781.6 | 28339 | 0 | 0 | 25 |
| has_phr | True | 538 | 1001 | 1136 | 1517 | 2981 | 15448.3 | 77700 | 0 | 0 | 228 |
| transmission | PAN028 maternal (hap1) vs PAN027 (mother) | 310 | 1001 | 1136 | 1724.5 | 4138 | 21061.6 | 77700 | 0 | 0 | 120 |
| transmission | PAN027 maternal (hap1) vs PAN010 (mother) | 167 | 1011 | 1136 | 1447 | 2159 | 4845.9 | 39881 | 0 | 0 | 73 |
| transmission | PAN027 paternal (hap2) vs PAN011 (father) | 61 | 1006 | 1136 | 1207 | 2589 | 7273 | 25497 | 0 | 0 | 35 |

## Lower-Merge Untangle Summary

Input lower-merge BEDs are generated by:
`bash scripts/pedigree/run_patch_tract_lower_merge.sh`

The generated intermediates are:
`paper_prep/_brainstorming/pedigree_patch_tract_lower_untangle/*.e50000.m0.n1.bed.gz`
and are intentionally not committed because they are large.

These rows are adjacent best-hit intervals merged by donor chrarm/haplotype
from `odgi untangle -e 50000 -m 0 -j 0 -n 1`, filtered to interchromosomal
merged intervals with `min_score >= 0.8` and length <= 100 kb. They are a
resolution check on the graph/untangle process, not the high-quality m1000 patch
denominator used for the manuscript's 494/538 within-community statistic.

Lower-merge overall: N = 37400; min 1 bp;
Q1 3 bp; median 17 bp; Q3
64 bp; max 55098 bp.

- In 22--95 bp conversion-like range:
  11600/37400
  (31.0%).
- In 318--688 bp CO-associated range:
  927/37400
  (2.5%).
- Near CO-associated two-fold window:
  3445/37400
  (9.2%).

| group_type | group_value | n | min | q25 | median | q75 | q95 | max | conv 22-95 | CO 318-688 | near CO 159-1376 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| lower_merge_overall | m0_n1_merged_interchromosomal | 37400 | 1 | 3 | 17 | 64 | 272.1 | 55098 | 11600 | 927 | 3445 |
| lower_merge_pattern | interchr_merged | 28974 | 1 | 3 | 15 | 59 | 252 | 55098 | 8494 | 610 | 2442 |
| lower_merge_pattern | sandwich_same_hap | 5095 | 1 | 8 | 28 | 71 | 250.3 | 27410 | 2202 | 165 | 504 |
| lower_merge_pattern | gene_conversion_like | 2537 | 1 | 4 | 21 | 107 | 837.2 | 17829 | 576 | 140 | 390 |
| lower_merge_pattern | crossover_like | 794 | 1 | 9 | 26 | 74.8 | 379.1 | 1305 | 328 | 12 | 109 |
| lower_merge_transmission | PAN028 maternal (hap1) vs PAN027 (mother) | 14054 | 1 | 3 | 19 | 69 | 335 | 55098 | 4393 | 365 | 1451 |
| lower_merge_transmission | PAN027 paternal (hap2) vs PAN011 (father) | 11963 | 1 | 5 | 26 | 74 | 252.9 | 18440 | 4248 | 337 | 1292 |
| lower_merge_transmission | PAN027 maternal (hap1) vs PAN010 (mother) | 11383 | 1 | 2 | 11 | 38 | 237 | 39471 | 2959 | 225 | 702 |
