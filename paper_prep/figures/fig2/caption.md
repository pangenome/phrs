# Figure 2 — Within-community heterogeneity, the two-domain model, and population history

**Inter-chromosomal exchange leaves a quantitative population signature.**
**(a) Allele vs. paralog distance.** Per-community Wilcoxon paired test
across nine multi-arm communities (5,946 pairs total;
`heterogeneity/allele_vs_paralog_distance.tsv`, SURVEY_04 §1.1). Bars show
the fraction of pairs where the cross-paralog comparison is closer than the
allelic. Allele is closer in 8/9 (overall p < 1e-300); C7 reverses
(70.5 % paralog-closer, p = 2.0e-7).
**(b) Two-domain model.** Top: per-arm Spearman ρ between number of
chromosome contributors and distance from the telomere (39/48 arms with a
negative gradient; `plots/two_domain_test.tsv`).
Bottom: binned mean number of chromosome contributors vs. distance for the
FISH-era anchor arms chr4p, chr4q, chr22q (`plots/two_domain_binned_means.tsv`),
with piecewise-regression breakpoints (39/41 arms favour piecewise over
linear; `plots/two_domain_changepoint.tsv`). 16/19 arms with internal-TTAGGG
blocks have an ITS within 50 kb of the breakpoint (`plots/its_breakpoint_coloc.tsv`).
**(c) Cross-arm superpop enrichment + Hudson Fst.** Left: −log₁₀ BH-adjusted
Fisher q for cross-arm vs. self-arm superpopulation composition, community
× arm (10/19 significant at q < 0.05;
`heterogeneity/cross_arm_superpop_enrichment.tsv`). Right: Hudson pairwise Fst
on cross-arm fields (`heterogeneity/fst_superpop_matrix.tsv`); AFR vs. non-AFR
Fst 0.10–0.15.
**(d) Out-of-Africa tree.** UPGMA dendrogram from the cross-arm Fst matrix
(9 arms × 5 superpopulations). AFR splits first; AMR/EAS/EUR/SAS form a tight
non-AFR clade.
