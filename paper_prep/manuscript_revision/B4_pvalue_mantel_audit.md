# B4 Pointwise P-value and Mantel Reporting Audit

Date: 2026-06-17

Scope: mechanical-to-judgment audit of pointwise Spearman p-values and
astronomical p-values in the active manuscript, active submission figure labels,
and still-present paper-prep captions/figure generators that could be reused by
an integrator. No new p-values were computed. Cached data and B0 source inventory
are sufficient for the recommendations below; no Slurm script is needed.

## Reading Guide

Reporting classes:

- **effect size**: keep the point estimate and sample count as descriptive; do
  not use the nominal pointwise p-value as primary inference.
- **Mantel p**: if an inferential p-value is needed for sequence-contact matrix
  concordance, use a row+column permutation Mantel/permutation statistic already
  present in cached global-test outputs, with an explicit permutation floor.
- **caveated**: retain only with the design caveat stated nearby.
- **remove**: drop the p-value or move it out of headline text/labels.

The central problem is dependence. PHR pairs, arm pairs, and matrix cells share
rows, columns, arms, genomes, community labels, and read-placement histories. A
nominal `cor.test` p-value over all dots treats those dots as independent. It is
useful as a software diagnostic, but it should not carry the manuscript's
inferential weight.

## Primary Audit Table

| Location | Current statistic | Independence problem | Recommended reporting class | Exact manuscript / figure-label change |
|---|---|---|---|---|
| `submission/paper.tex:219-222` Results, human pointwise contact paragraph | Human pointwise Spearman effects: HG002 Pore-C rho = 0.38; CHM13 Hi-C rho = 0.72; HG002 Hi-C rho = 0.66; HG002 CiFi rho = 0.19; "all p < 10^-20" | The rho values are descriptive correlations over inter-chromosomal PHR sequence pairs. Pairs are not independent because the same PHRs, arms, haplotypes, and contact-map rows/columns recur across many dots. | effect size / remove | Keep the rho values and assay replication. Delete "all p < 10^-20". Suggested replacement: "with the same direction in a second genome and two further assays ...; Extended Data Fig. ED1." If an inferential statement is required, point to the cached arm-level row+column Mantel/permutation family in Methods, not to pointwise p-values. |
| `submission/paper.tex:247-248` Results, pedigree paragraph | 494/538 high-quality patches (92%) within sequence community vs marginal-aware permutation null mean 77.0%, p < 10^-4 | This p-value is not a pointwise Spearman p-value. It is a Monte Carlo null preserving marginal structure, so it is the right inferential class for the pedigree fraction. It is still bounded by the number of permutations and by patch-calling dependence within families/arms. | caveated | Keep as permutation evidence, but report with the permutation floor and CI. Suggested wording: "above a marginal-aware Monte Carlo null (null mean 77.0%; permutation p < 10^-4; Wilson 95% CI for observed fraction 89.2-93.9%)." Avoid calling it exact. |
| `submission/paper.tex:387-388` Fig. 4A caption | HG002 Pore-C pointwise Spearman rho = 0.381, n = 2,830, p = 1.2 x 10^-98 | Each dot is a PHR-pair measurement, but dots reuse PHRs and arms; the nominal p-value assumes independent pairs and is astronomically small mostly because n is large and dependent. | effect size / remove | Caption should keep "rho = 0.381, n = 2,830" and drop "p = 1.2 x 10^-98". Suggested replacement: "contact rises with similarity (descriptive pointwise Spearman rho = 0.381, n = 2,830; line, linear fit)." |
| `submission/scripts/figures/make_fig4a_human_scatter.R:75-79` active Fig. 4A label generator | Burns label: n, "pointwise Spearman rho", and "p = ..." for HG002 Pore-C | Same pointwise-pair dependence as Fig. 4A caption; the figure label gives the p-value visual authority. | effect size / remove | Remove the `sprintf("p = %s", fmt_p(it$p))` legend line. Keep n and rho. Optionally relabel rho as "descriptive pointwise Spearman rho". |
| `submission/paper.tex:391-394` Fig. 4B caption | HG002 Pore-C community matrix: B/W = 0.056, p = 3.9 x 10^-85 | The displayed statistic is a within-vs-between contact contrast over matrix cells/arm-haplotype pairs; cells share rows/columns and community labels. A raw/global p-value is not transparent unless described as a permutation/Mann-Whitney/bootstrap test with its dependency limits. The same cached global table also has Mantel rho = 0.486 with row+column permutations, but B0 notes p was reported as 0.0 rather than a floor. | Mantel p / caveated | In the caption, prefer the effect size and visual matrix: "within-community blocks dominate (B/W = 0.056)." If retaining inference, change to a permutation-floor statement from the global-test table, e.g. "row/column-permutation Mantel support in Methods; report p < 1/10,001 rather than p = 0.0." Do not headline 3.9 x 10^-85. |
| `submission/scripts/figures/make_fig4b_porec_community.R:500-522` active Fig. 4B statistic-check text | Records source B/W and p-value; displays p-value as `sprintf("%.1e", wb$p_value)` in the generated audit note | This is a provenance note, not a manuscript label, but it perpetuates the astronomical B/W p-value without the row/column-dependence caveat. | caveated | Keep source provenance if useful, but add a note that the p-value is a cached global-test value and should not be used as headline inference; prefer B/W plus row+column Mantel/permutation wording in manuscript. |
| `submission/paper.tex:399-404` Fig. 4C caption | Mouse zygotene pointwise Spearman rho = 0.614, n = 1,135, p = 1.2 x 10^-118; per-stage rhos 0.419/0.614/0.576/0.245; all p < 1 x 10^-16 | The per-PHR-pair mouse dots reuse PHRs/arms, so nominal pointwise p-values assume independence that is not present. The stage trajectory itself is descriptive and biologically interpretable. Cached arm-level Mantel with 10,000 row+column permutations exists for the same stages and is the better inferential family. | effect size + Mantel p | Drop the pointwise p-values from the caption. Keep rho and n for the scatter and keep the four-stage rho trajectory. Add, if space permits: "arm-level row+column Mantel also peaks at zygotene (Methods)." Do not report "all p < 1 x 10^-16" for the pointwise stage rhos. |
| `submission/scripts/figures/make_fig4c_mouse_zygotene.R:175-179` active Fig. 4C label generator | Burns label: n, "pointwise Spearman rho", and "p = ..." for mouse zygotene | Same pointwise-pair dependence as Fig. 4C caption. | effect size / remove | Remove the `sprintf("p = %s", fmt_p(ct$p.value))` legend line. Keep n and rho. |
| `submission/paper.tex:553-555` Methods, Hi-C/Pore-C/CiFi pipeline | "W/B ratio computed by bootstrap on 10,000 permutations; Mann-Whitney global test; Mantel Spearman with 10,000 row+column permutations" | This is the correct methods hook for inferential statistics, but the wording conflates W/B bootstrap, Mann-Whitney, and Mantel in one sentence. | Mantel p | Clarify that row+column Mantel/permutation tests, not pointwise Spearman p-values, are the inferential family for matrix-level sequence-contact association. Suggested addition: "Nominal pointwise Spearman p-values are treated as descriptive diagnostics because sequence-pair observations share arms/PHRs." |
| `submission/paper.tex:557-564` Methods, pointwise Spearman paragraph | Human pointwise rhos and n are listed, with no pointwise p-values in this paragraph | This is already mostly correct: it separates descriptive rho/n from p-values. It should explicitly state why p-values are omitted from Methods. | effect size | Keep as-is but append a caveat sentence: "These pointwise correlations are reported as descriptive effect sizes; matrix-level inference used row+column Mantel/permutation tests above." |
| `submission/paper.tex:595-597` Methods, PBMC Dip-C negative control | PBMC W/B = 0.983, p = 0.305 | Not a pointwise Spearman p-value and not astronomical. It is a negative-control p-value over cells/summary output. It carries limited inferential weight because hg19-projected PHR boundaries lack a pairwise Jaccard matrix and PBMC is a control. | caveated | Keep only if single-cell controls remain. Add "negative control" and source/caveat. No Mantel replacement needed unless B/F fan-in keeps Dip-C as a major claim. |
| `submission/paper.tex:607-614` Methods, mouse pipeline | Per-PHR-pair stage rhos 0.419/0.614/0.576/0.245; all p < 1 x 10^-16; arm-level Mantel trajectory 0.687/0.718/0.683/0.577; all p < 1 x 10^-4 from 10,000 row+column permutations | The first p-value family is pointwise and dependent. The second family uses arm-level matrices and row+column permutations and is the correct inference class; however, with 10,000 permutations the smallest defensible report is a permutation floor. | Mantel p / remove | Remove "all p < 1 x 10^-16" for the pointwise rhos. Keep the arm-level Mantel trajectory and report its p-values as permutation floors, e.g. "all row+column permutation p < 1/10,001" if zero exceedances, or exact cached permutation p-values if available. |
| `submission/paper.tex:659-663` Methods, limitations/statistical confidence sentence | Claims confidence intervals are reported for headline correlations; includes "observed p < 10^-4" for pedigree Monte Carlo null | This sentence overpromises CIs for headline correlations and still gives the pedigree p-value without floor language. It is not a pointwise Spearman issue, but it is part of the p-value audit. | caveated | Revise after B/F fan-in decides final stats. Suggested wording: "Where retained, headline correlations are reported as effect sizes with compatible matrix/permutation or interval summaries; the pedigree within-community fraction is compared with a marginal-aware Monte Carlo null (p bounded by the permutation count)." |
| `submission/paper.tex:685-688` Extended Data Fig. ED1 caption | CHM13 Hi-C rho = 0.716, n = 652, p = 1.2 x 10^-103; HG002 Hi-C rho = 0.662, n = 2,544, p < 10^-300; HG002 CiFi rho = 0.191, n = 2,757, p = 3.7 x 10^-24 | Same pointwise-pair dependence as Fig. 4A. The HG002 Hi-C p < 10^-300 is especially misleading as an independence-based astronomical value. | effect size / remove | Keep rho and n for each panel; remove all three p-values. Suggested replacement: "Contact rises with similarity in every panel: CHM13 Hi-C (descriptive rho = 0.716, n = 652), HG002 Hi-C (rho = 0.662, n = 2,544), and HG002 CiFi (rho = 0.191, n = 2,757; sparsest assay)." |
| `submission/scripts/figures/make_ed1_human_contacts.R:76-80` active ED1 label generator | Burns label: n, Spearman rho, and p for CHM13 Hi-C, HG002 Hi-C, HG002 CiFi | Same pointwise-pair dependence as ED1 caption; the labels make nominal p-values appear as panel-level evidence. | effect size / remove | Remove the `sprintf("p = %s", fmt_p(it$p))` legend line. Keep n and rho. Optionally label "descriptive Spearman rho" for consistency with Fig. 4A. |

## Legacy Paper-prep Assets to Guard Against Reuse

These are not the active `submission/paper.tex` captions, but they contain the
same reporting hazards and should be updated or ignored if figures are
regenerated from `paper_prep/figures/`.

| Location | Current statistic | Independence problem | Recommended reporting class | Exact integrator action |
|---|---|---|---|---|
| `paper_prep/figures/fig3/caption.md:9` and `paper_prep/figures/fig3/figure_fig3.R:77-80` | HG002 Pore-C B/W = 0.056, p = 3.9 x 10^-85 | Matrix cells/arm-haplotypes share rows, columns, and community labels. | Mantel p / caveated | If this older Fig. 3 returns, show B/W as effect size and move p-value language to Methods/SI as row+column Mantel/permutation support with a permutation floor. |
| `paper_prep/figures/fig3/caption.md:11-13` and `paper_prep/figures/fig3/figure_fig3.R:100-203` | Forest plot labels 14 tests with ratio, p, and convention across Hi-C, Pore-C, CiFi, Dip-C, sperm, mouse | Mixed tests have different units: matrix cells, cells, stages, and arm pairs. Listing p-values side-by-side implies comparability that is not justified. | caveated / remove | If retained, make the forest plot an effect-size/technology robustness plot. Remove p-values from row labels or put them in a separate supplement table with test-specific caveats. |
| `paper_prep/figures/fig3/caption.md:22` and `paper_prep/figures/fig3/figure_fig3.R:345-378` | Dip-C radial inset p = 1.6 x 10^-35 | Radial particles/cells are dependent and this is a defensive/contact-location control, not the lead sequence-contact claim. | caveated | Retain only as a control effect size ("radial 0.504 vs 0.556") and move/drop the astronomical p-value. |
| `paper_prep/figures/fig4/caption.md:47-51` and `paper_prep/figures/fig4/figure_fig4.R:249-276` | Older mouse zygotene per-PHR-pair Spearman rho = 0.715, p = 4.4 x 10^-55, n = 344 inter-chromosomal pairs | Older 50 kb/344-pair mouse statistic is superseded by active submission Fig. 4C 20 kb/1,135-pair analysis. It is still pointwise and dependent. | remove | Do not reuse this older p-value in the manuscript. If the older panel is shown historically, label rho/n as descriptive and point to current Methods for row+column Mantel inference. |
| `paper_prep/figures/ed5/caption.md:5` and `paper_prep/figures/ed5/figure_ed5.R:211-241` | Mantel rho before vs after acrocentric+sex exclusion; no p-values in caption | This is a row+column matrix statistic family and is exactly the kind of evidence that should carry inferential weight if p-values are reported. Caption currently uses effect sizes only. | Mantel p / effect size | Keep as effect-size robustness. If p-values are added, use cached row+column permutation Mantel p-values with floor language; do not add pointwise p-values. |
| `paper_prep/figures/ed5/caption.md:9` and `paper_prep/figures/ed5/figure_ed5.R:149-169,288-348` | Per-community reproducibility heatmap marks `* q < 0.05`, `** q < 0.001` from 10,000 random-label permutations | Community-level permutation/q-values are closer to the right inference class, but labels across 15 communities x 11 datasets still involve repeated, partly dependent tests. | caveated | Retain as robustness if needed; call these random-label permutation q-values, not independent pointwise p-values. |
| `paper_prep/figures/ed8/caption.md:3` and `paper_prep/figures/ed8/figure_ed8.R:131,191` | Feedback-loop panel cites CHM13 Hi-C rho = 0.674 and HG002 Pore-C rho = 0.485 as direct measurement, described as Mantel/community-free in legend | B0 identifies these values as per-arm-pair Spearman from `*_phr_pair_correlation.tsv`, not Mantel. Arm pairs share arms/rows/columns. | caveated / effect size | Correct label if reused: "per-arm-pair Spearman" rather than Mantel. Use row+column Mantel values from `*_global_test.tsv` if the point is matrix-level inference. |
| `paper_prep/figures/ed8/caption.md:5` | D4Z4 panel reports Mann-Whitney p = 5.3e-6 for DUX4L outliers | Not part of sequence-contact Mantel issue, but it is a p-value in a legacy caption. Observations likely share communities/arms. | caveated | Retain only if ED8 returns, and describe as a descriptive enrichment/outlier comparison; do not use as central inference. |
| `paper_prep/figures/ed8/caption.md:7` and `paper_prep/figures/ed8/figure_ed8.R:89-95,300-326` | Recombination vs cross-arm affinity: all arms rho = -0.35, p = 0.017; well-callable rho = -0.01, p = 0.97 | Per-arm points are fewer but still share ascertainment and callable-variant confounding; the caption itself identifies the confound. | caveated | Keep the "honest null" framing and emphasize the loss of effect after callable filtering. Do not treat the all-arm p = 0.017 as positive evidence. |
| `paper_prep/figures/fig2/caption.md:8-9` and `paper_prep/figures/fig2/figure_fig2.py:90-116` | Allele-vs-paralog Wilcoxon overall p < 1e-300; C7 p = 2.0e-7; per-community labels include p-values | Not a pointwise Spearman statistic; pairs are many and likely dependent within communities/arms. It is a legacy caption for a cut/reviewer-era figure. | remove / caveated | Do not reintroduce the astronomical p < 1e-300 into the manuscript. If the panel returns, lead with effect sizes/fractions and put test details in supplement with dependence caveat. |
| `paper_prep/figures/ed3/caption.md:7` and `paper_prep/figures/ed3/figure_ed3.R:143` | Telomere length by community: Kruskal-Wallis H = 100.89, p = 3.2e-15 across 13,668 sequences | Not a sequence-contact statistic, but an astronomical p-value over many sequences that share individuals/arms/communities. | caveated / remove | If ED3 returns, report community medians/ranges; remove or downweight the p-value unless a clustered model is supplied. |
| `paper_prep/figures/ed4/caption.md:3-5` and `paper_prep/figures/ed4/figure_ed4.R:36-80` | GO enrichment p_adj = 1.45e-3; copy_pvalue approximately 0; dedup p = 0.040 | Gene-copy entries are not independent and B0/task context says gene enrichment was cut from the active body. | remove / caveated | Do not reintroduce as headline evidence. If kept in supplement, use adjusted enrichment terms with copy/de-dup caveats. |

## Statistics That Should Carry Inferential Weight

Cached sources identified by B0 and manuscript Methods are enough for this
revision stage:

| Evidence family | Cached source | How to report |
|---|---|---|
| Human arm-level sequence-contact concordance | `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/*_global_test.tsv`; repo/audit summaries in `paper_prep/_brainstorming/stats_audit/mantel_multires_si_table.tsv`, `scripts/ci/mantel_fisher_z_ci.tsv`, `scripts/ci/mantel_bootstrap_ci.tsv` | Use Mantel rho as matrix-level effect size and row+column permutation p-values as inference. Avoid p = 0.0; report the permutation floor, e.g. p < 1/(permutations + 1), or exact cached nonzero p where present. |
| Human flanking unique-sequence control | `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/flanking/{5000,10000,20000,50000,100000}bp/*_global_test.tsv` and community-free flanking summaries | Use as defensive control against mapping artifact. Prefer flanking Mantel/W/B effect sizes and row+column/random-label permutation language; avoid pointwise flanking p-values as headline inference. |
| Mouse stage-resolved sequence-contact relationship | `submission/paper.tex:610-616`; `scripts/mouse/mantel_d_m5.{py,R}`; upstream mouse global-test tables | Use per-stage pointwise rho trajectory descriptively, and use arm-level Mantel with 10,000 row+column permutations for inference. Report permutation floors rather than nominal pointwise p-values. |
| Pedigree within-community exchange fraction | `submission/paper.tex:419-420`; `scripts/pedigree/monte_carlo_null_d_m4.py` | This is already the correct null class. Report observed 92%, Wilson CI, null mean/CI, and permutation p bounded by the number of Monte Carlo draws. |

## Recommended Integrator Patch Checklist

1. Remove pointwise Spearman p-values from active Fig. 4A, Fig. 4C, and ED1
   captions in `submission/paper.tex`.
2. Remove pointwise p-value legend lines from
   `submission/scripts/figures/make_fig4a_human_scatter.R`,
   `submission/scripts/figures/make_fig4c_mouse_zygotene.R`, and
   `submission/scripts/figures/make_ed1_human_contacts.R`; regenerate PDFs/PNGs
   only in the integrator/edit task.
3. Keep descriptive rho and n for pointwise scatters. Use the phrase
   "descriptive pointwise Spearman" once, then use "rho" consistently.
4. Clarify Methods that pointwise Spearman p-values are not used for inference
   because sequence pairs are dependent.
5. Use cached row+column Mantel/permutation outputs for matrix-level inference,
   with permutation-floor wording. Do not print p = 0.0 or astronomical
   scientific-notation p-values for matrix cells.
6. Keep the pedigree p-value only as a Monte Carlo permutation result with
   finite-permutation caveat.
7. Treat legacy `paper_prep/figures` p-values as non-authoritative unless the
   corresponding figure is resurrected; if resurrected, apply the same rules.

## Validation Checklist

- Artifact exists: this file.
- Lists all active pointwise Spearman p-value locations found:
  `submission/paper.tex:219-222`, `387-388`, `399-404`, `607-614`,
  `685-688`; active figure-label generators
  `submission/scripts/figures/make_fig4a_human_scatter.R:75-79`,
  `make_fig4c_mouse_zygotene.R:175-179`, and
  `make_ed1_human_contacts.R:76-80`.
- Separates descriptive rho from inferential p-values: see Reporting classes,
  Primary Audit Table, and Recommended Integrator Patch Checklist.
- Recommends exact manuscript changes for integrator: see the final column of
  the Primary Audit Table and checklist above.
- No new p-values computed; no heavy head-node work or Slurm run required.
