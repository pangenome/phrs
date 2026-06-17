# C0a Arm-Level Continuum Characterization

Date: 2026-06-17

## Scope

This is the compute-light C0a arm-level analysis. It uses the 41 x 41 arm-level Jaccard distance matrix and the arm-level Leiden k=15 assignments only. Similarity is computed as `1 - Jaccard distance` for off-diagonal arm pairs; the nonzero matrix diagonal is ignored. This report does not traverse or summarize the 15,668 x 15,668 sequence-level evidence, so it should be treated as arm-level support for the continuum framing rather than sequence-level proof.

## Inputs

- Distance matrix: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`
- Leiden assignments: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
- Matrix arms: 41
- Off-diagonal arm pairs: 820

## Arm-Level Similarity Distributions

| category | n | mean | median | q05 | q95 | max |
|---|---:|---:|---:|---:|---:|---:|
| all_off_diagonal | 820 | 0.0685 | 0.0097 | 0.0000 | 0.3567 | 0.9847 |
| within_leiden_community | 58 | 0.4506 | 0.4312 | 0.2572 | 0.7993 | 0.9847 |
| between_leiden_communities | 762 | 0.0394 | 0.0078 | 0.0000 | 0.1921 | 0.4155 |

The arm-level matrix shows a two-tier pattern: within-Leiden arm pairs are much more similar on average than between-community pairs, but the off-diagonal background is not a hard zero class. That supports continuum language at arm level: named communities sit on a broad background of lower, variable similarity rather than forming perfectly isolated blocks.

## Named Systems

| system | arms | pair_count | mean | median | min | max | peak_pair |
|---|---|---:|---:|---:|---:|---:|---|
| C1_D4Z4_DUX4 | 4q, 10q | 1 | 0.5830 | 0.5830 | 0.5830 | 0.5830 | 10q-4q (0.5830) |
| C2_10p_18p_TUBB8B | 10p, 18p | 1 | 0.6202 | 0.6202 | 0.6202 | 0.6202 | 10p-18p (0.6202) |
| C6_q_arm_sextet | 1q, 13q, 17q, 19q, 21q, 22q | 15 | 0.5003 | 0.4338 | 0.1633 | 0.8169 | 21q-22q (0.8169) |
| C7_acrocentric_p | 13p, 14p, 15p, 21p, 22p | 10 | 0.4428 | 0.4383 | 0.4000 | 0.4794 | 15p-22p (0.4794) |
| C11_OR4F_core_5q_6q | 5q, 6q | 1 | 0.6027 | 0.6027 | 0.6027 | 0.6027 | 5q-6q (0.6027) |
| C11_OR4F_full | 1p, 5q, 6q, 8p | 6 | 0.4058 | 0.3988 | 0.2520 | 0.6027 | 5q-6q (0.6027) |
| C14_Xq_Yq | Xq, Yq | 1 | 0.9847 | 0.9847 | 0.9847 | 0.9847 | Xq-Yq (0.9847) |
| C15_Xp_Yp | Xp, Yp | 1 | 0.7114 | 0.7114 | 0.7114 | 0.7114 | Xp-Yp (0.7114) |

## C6 / q-arm Sextet Density

C6 arms tested: 1q, 13q, 17q, 19q, 21q, 22q.

| region | n | mean | median | q05 | q95 | max | fold_vs_C6_within |
|---|---:|---:|---:|---:|---:|---:|---:|
| C6_within_sextet | 15 | 0.5003 | 0.4338 | 0.2297 | 0.8142 | 0.8169 | 1.0000 |
| C6_to_non_C6 | 210 | 0.0318 | 0.0165 | 0.0000 | 0.1129 | 0.2188 | 15.7084 |
| non_C6_pairs_only | 595 | 0.0706 | 0.0073 | 0.0000 | 0.3679 | 0.9847 | 7.0892 |
| off_diagonal_excluding_C6_within | 805 | 0.0605 | 0.0090 | 0.0000 | 0.3158 | 0.9847 | 8.2735 |
| between_leiden_communities | 762 | 0.0394 | 0.0078 | 0.0000 | 0.1921 | 0.4155 | 12.6893 |

Exact six-arm set test over all 4496388 possible six-arm subsets: C6 observed mean = 0.5003, null mean = 0.0685, null 95th percentile = 0.1414, percentile = 100.00%, exact greater-tail p = 4.448e-07. Wilcoxon greater-tail p versus all off-diagonal pairs excluding C6-within = 2.786e-10; versus between-community pairs = 3.652e-11.

Interpretation: the q-arm sextet block is visibly and quantitatively denser than the off-diagonal background in the 41-arm matrix. Because C6 was defined from this same arm-level similarity matrix, the exact set and Wilcoxon tests are descriptive diagnostics for heatmap-density decision-making, not independent discovery p-values. The result supports showing C6 as a dense local region while avoiding language that treats the sextet as a closed, sequence-level clade without the separate sequence-level evidence.

## Outputs

- `arm_assignments_used.tsv`: normalized arm labels and Leiden assignments used in this run.
- `arm_pair_similarity_long.tsv`: all 820 off-diagonal arm-pair similarities.
- `similarity_distribution_summary.tsv`: all / within-community / between-community similarity distributions.
- `community_similarity_summary.tsv`: per-community arm-level density and peak pairs.
- `named_system_peak_similarities.tsv`: named-system means and peak similarities.
- `c6_neighborhood_density.tsv`: C6 block, C6-neighbor, and background density summaries.
- `c6_within_pair_similarities.tsv`: the 15 C6 sextet arm pairs, sorted by similarity.
- `c6_exact_set_test.tsv`: exact six-arm-set density diagnostic plus Wilcoxon comparisons.
- `arm_level_similarity_diagnostic.png` and `.pdf`: community-ordered heatmap and similarity distribution diagnostic.
- `compute_arm_level_continuum.R`: reproducible script for all outputs.
