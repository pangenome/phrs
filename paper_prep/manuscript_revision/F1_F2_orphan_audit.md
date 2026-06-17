# F1-F2 Orphan Audit: F_ST and cM/Mb

Date: 2026-06-17
Task: `manuscript-revision-f1`

## Executive Recommendation

1. **F_ST: present, not absent; keep demoted/qualified, do not promote as a main causal result.**
   The analysis exists in upstream outputs, figure assets, local CI scripts, and a matched-control script. Its safe interpretation is that cross-arm/self-arm structural haplotypes preserve ordinary continental ancestry structure at approximately genome-wide background levels. It does **not** support a strong "subtelomere-specific population-genetic signature" claim, and should not be promoted as proof of ongoing exchange. If retained, use it as a short population-structure qualifier or Methods/Extended Data support for Fig. 2c/2d, with the matched-control caveat attached.

2. **cM/Mb anti-correlation: provenance resolved; clarify or cut from active manuscript.**
   The recombination-rate input is from the Lalli 2025 T2T-CHM13 short-read recombination map, but the join to PHR cross-arm affinity and the low-callability filter are our analysis. The apparent negative correlation supports the paper only as an **honest null/limitation**: it collapses after excluding low-callability acrocentric/PAR arms, so current short-read maps cannot test whether local recombination protects arm identity. The active `submission/paper.tex` currently includes only a compressed limitations clause; that is acceptable only if it is rephrased to make clear that the **Lalli input is cited work and the collapse-after-filter is our reanalysis**. Otherwise, cut the clause because, as written, it can read as an unexplained external claim.

## Scope Checked

- Active manuscript: `submission/paper.tex`
- End-to-end report sections: `end-to-end-report/report/04_heterogeneity.md`, `07_integrated.md`, `12_literature.md`
- Local scripts and outputs: `scripts/popgen/`, `scripts/ci/`
- Figure assets/manifests: `paper_prep/figures/fig2/`, `paper_prep/figures/ed8/`
- Synthesis decision records: `paper_prep/synthesis/ANALYSIS_D_M6.md`, `ANALYSIS_D_M12.md`, `NATURE_DRAFT_v6.md`, `OPEN_REVIEWER_CONCERNS.md`

## F_ST Audit

### Source Paths

Primary upstream analysis and outputs:

- `end-to-end-report/report/04_heterogeneity.md:101-113` defines the binary allele as `self-arm = 0` vs `cross-arm = 1`, reports Hudson Fst averaged across the 10 strongest arm/community pairs, and gives AFR/non-AFR Fst 0.10-0.15 with non-AFR pairs near zero.
- `end-to-end-report/report/04_heterogeneity.md:288-313` lists the upstream script `/moosefs/guarracino/HPRCv2/scripts/community/compute_fst_superpop.py` and output `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/fst_superpop_matrix.tsv`.
- `paper_prep/figures/fig2/figure_fig2.py:225-299` builds Fig. 2c from cross-arm superpopulation enrichment plus the Hudson Fst matrix.
- `paper_prep/figures/fig2/figure_fig2.py:299-422` builds the out-of-Africa UPGMA dendrogram from the Fst matrix.
- `paper_prep/figures/fig2/sources.tsv:8-10` names `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/fst_superpop_matrix.tsv` as the source for Fig. 2c/2d.

Local CI/reviewer-response scripts and outputs:

- `scripts/ci/bootstrap_ci_d_m12.py:29-34` states that D-M12 computes Hudson F_ST per superpopulation pair and block-jackknife CIs over 10 arm/community blocks.
- `scripts/ci/bootstrap_ci_d_m12.py:75-92` embeds the 10 cross-arm/self-arm count rows from `04_heterogeneity.md`.
- `scripts/ci/bootstrap_ci_d_m12.py:306-338` writes `scripts/ci/fst_block_jackknife.tsv` and `scripts/ci/fst_per_arm_per_pair.tsv`.
- `scripts/ci/fst_block_jackknife.tsv` exists and reports AFR-AMR 0.102 [0.022, 0.182], AFR-EAS 0.152 [0.065, 0.240], AFR-EUR 0.108 [-0.024, 0.240], AFR-SAS 0.103 [-0.001, 0.208], with all non-AFR pair CIs bracketing zero.
- `scripts/ci/fst_per_arm_per_pair.tsv` exists and gives the per-arm contribution for every population pair.

Matched-control analysis:

- `scripts/popgen/matched_fst_d_m6.py:1-45` states the M6 task: compare subtelomeric Hudson F_ST against matched non-subtelomeric autosomal F_ST in the same superpopulation frame.
- `scripts/popgen/matched_fst_d_m6.py:63-90` names the default moosefs input and the published 1000G/HGDP baseline F_ST values.
- `scripts/popgen/matched_fst_d_m6.py:97-109` implements the same Hudson-style F_ST formula.
- `scripts/popgen/matched_fst_d_m6.py:300-369` emits the matched-control comparison table and optional TSV.
- `paper_prep/synthesis/ANALYSIS_D_M6.md:9-23` records what was computed and why the published 1000G/HGDP baselines were used.
- `paper_prep/synthesis/ANALYSIS_D_M6.md:51-76` gives the key verdict: AFR/non-AFR subtelomeric F_ST is statistically indistinguishable from matched genome-wide background; no pair is elevated.
- `paper_prep/synthesis/ANALYSIS_D_M12.md:57-121` documents the block-jackknife CIs and explicitly says this CI does not by itself answer the matched-background question.

Active manuscript state:

- `submission/paper.tex:650-663` mentions only that headline CIs include an `$F_{\mathrm{ST}}$` block-jackknife CI. There is no active Results paragraph or Methods subsection in the current LaTeX that explains the F_ST analysis.
- `paper_prep/synthesis/NATURE_DRAFT_v6.md:52` contains the fuller intended Results sentence: cross-arm sequences carry population structure consistent with genome-wide patterns, not a subtelomere-specific signature.
- `paper_prep/synthesis/NATURE_DRAFT_v6.md:94-96` contains the fuller intended Methods-only F_ST explanation.

### Validation Run

I ran:

```bash
python3 scripts/popgen/matched_fst_d_m6.py --out-tsv /tmp/f1_fst_check.tsv
```

Result summary:

- Input read: `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_superpop_enrichment.tsv`
- Blocks: 10 arm/community pairs with `p_adjusted < 0.05`
- AFR/non-AFR pairs:
  - AFR-AMR: 0.1020, 95% jackknife CI [0.0203, 0.1837], matched baseline 0.0710, verdict equivalent
  - AFR-EAS: 0.1525, CI [0.0634, 0.2416], matched baseline 0.1440, verdict equivalent
  - AFR-EUR: 0.1080, CI [-0.0264, 0.2424], matched baseline 0.1500, verdict equivalent
  - AFR-SAS: 0.1034, CI [-0.0030, 0.2098], matched baseline 0.1100, verdict equivalent
- Non-AFR pairs are equivalent or depressed relative to the matched genome-wide baseline; none is elevated.

Note: running `python scripts/popgen/matched_fst_d_m6.py` fails in this environment because `python` is Python 2 and the script has non-ASCII text in the header. Use `python3`.

### Interpretation

The F_ST analysis is real and reproducible from the current worktree plus moosefs inputs. It is also already partly integrated into the manuscript-prep trail:

- It supports the modest claim that subtelomeric structural haplotypes preserve population structure.
- It does **not** support a stronger claim that subtelomeres are unusually differentiated above genome-wide background.
- It is not direct evidence for ongoing exchange; the direct exchange evidence remains the pedigree analysis.
- It should not be sold as a standalone manuscript result unless a main text paragraph also includes the matched-control caveat.

### Exact Recommendation for F_ST

Do **not** mark F_ST absent. Mark it as **present but demoted/qualified**.

Recommended downstream action:

- If Fig. 2 population panels remain in the manuscript, add a short Methods/legend sentence equivalent to:

> Hudson F_ST was computed on a binary cross-arm/self-arm structural-haplotype state across the 10 significant arm/community blocks; AFR/non-AFR values are 0.10-0.15 but are statistically indistinguishable from matched 1000G/HGDP genome-wide autosomal Hudson F_ST, so the result is an ancestry-preservation/background-population-structure control rather than a subtelomere-specific differentiation signal.

- If the manuscript is being compressed, move the F_ST result to Extended Data or Methods and keep only the Fig. 2c heatmap if needed.
- Do not promote the F_ST tree as a main biological endpoint. The out-of-Africa topology is expected background, not a new subtelomeric mechanism.

## cM/Mb Anti-Correlation Audit

### Source Paths

Report/survey substrate:

- `end-to-end-report/report/07_integrated.md:97-105` frames the recombination-rate correlation as a testable prediction, reports rho = -0.43, p = 0.006 across all 39 shared arms, then says the signal is driven by seven low-callability acrocentric/PAR arms and vanishes after exclusion (rho = 0.00, p = 0.98, N = 32).
- `end-to-end-report/report/12_literature.md:73-99` gives the more detailed table and states the same confound control.

Figure asset:

- `paper_prep/figures/ed8/sources.tsv:5-7` records the actual inputs:
  - `/moosefs/guarracino/HPRCv2/PHR_III/recombination_maps/subtelomeric_recomb_rates.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_affinity_sequences.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv`
- `paper_prep/figures/ed8/figure_ed8.R:295-353` plots recombination rate vs cross-arm affinity, highlights low-callability points, and renders all-arms vs well-callable-only panels.
- `paper_prep/figures/ed8/caption.md:7` states the ED8c result: Lalli 2025 cM/Mb input vs our per-arm cross-arm sharing fraction; all 46 arms rho = -0.35, p = 0.017; well-callable arms only n = 40 rho = -0.01, p = 0.97; survey-reported values rho = -0.43 (39 arms) to rho = 0.00 (N = 32).

Active manuscript state:

- `submission/paper.tex:650-658` includes only the limitations clause: "the short-read cM/Mb anti-correlation reported by Lalli 2025 collapses once seven low-callability arms are excluded."
- `paper_prep/synthesis/NATURE_DRAFT_v6.md:116` has the same limitations phrasing.

### Provenance Resolution

This line is **not purely cited work** and **not purely our de novo recombination map**.

Resolved provenance:

- Lalli 2025 supplies the T2T-CHM13 subtelomeric recombination-rate/cM/Mb map.
- Our analysis joins that map to PHR cross-arm affinity and sequence-community assignments.
- Our analysis performs the low-callability filter and concludes the anti-correlation collapses.

Therefore, "reported by Lalli 2025" is too compressed. It makes the collapse sound like a Lalli result, when the tested relationship is the PHR manuscript's derived correlation between Lalli's recombination-rate table and our cross-arm affinity metric.

### Interpretation

The cM/Mb anti-correlation does **not** undercut the paper's core claims, because the core claim is now supported by:

- sequence communities,
- 3D contact enrichment,
- mouse zygotene Hi-C,
- pedigree-resolved exchange.

It does undercut any claim that current short-read recombination maps can prove a recombination-rate-vs-exchange relationship at subtelomeres. In other words:

- As positive evidence: weak/unsafe; do not use.
- As an honest null: useful, if space allows.
- As a limitation: acceptable, but clarify provenance and interpretation.

### Exact Recommendation for cM/Mb

Preferred action: **clarify if retained; otherwise cut.**

If the limitations paragraph keeps it, replace the compressed clause with a provenance-clear version such as:

> A PHR-level comparison of Lalli 2025 T2T-CHM13 short-read cM/Mb estimates with our cross-arm affinity metric shows an apparent anti-correlation only before filtering; it collapses after seven low-callability acrocentric/PAR arms are excluded, so current short-read recombination maps cannot test the recombination-rate relationship at PHRs.

If word count is tight, cut the cM/Mb clause entirely. It is not needed for the manuscript's positive argument and currently creates more reader burden than evidentiary value.

Do **not** cite it as "Lalli reported the anti-correlation" unless the manuscript explicitly states that the anti-correlation is our analysis using Lalli's map.

## Orphan Status Summary

| Item | Exists? | Intended for manuscript? | Support level | Recommendation |
|---|---:|---|---|---|
| Hudson F_ST cross-arm/self-arm analysis | Yes | Yes, but only demoted/qualified | Population structure preserved; not subtelomere-specific elevation | Keep only with matched-control caveat; do not promote |
| F_ST block-jackknife CI | Yes | Yes, if F_ST retained | Reporting-standard support | Keep with F_ST Methods/legend text |
| F_ST matched-control analysis | Yes | Should control any F_ST claim | Shows no elevation over genome-wide background | Required caveat if F_ST appears |
| cM/Mb anti-correlation | Yes | Only as limitation/honest null | Positive signal collapses after callability filter | Clarify provenance or cut |

## Final Decision Record

F_ST is **not absent**. It is a real analysis with local scripts and outputs, but its manuscript role should be conservative: background-level ancestry structure, not a mechanistic or subtelomere-specific claim. The cM/Mb line is **resolved** as our reanalysis of a Lalli 2025 recombination-rate input; it should be clarified if retained and can be cut without weakening the paper.
