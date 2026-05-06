# Hi-C visual redesign candidates for review zoom slides 10a-12b

Task: `review-zoom-hic-figure-redesign`
Date: 2026-05-06 UTC
Scope: isolated candidate assets only. The deck source is not edited here.

## Outputs

Run:

```bash
Rscript slides/v2-review-zoom/_revision_assets/hic_visual_redesign/make_hic_visual_redesign.R
```

Generated candidate assets:

| File | Purpose |
|---|---|
| `slide_10a_square_matrix_candidate.png` | Square 1:1 HG002 Pore-C contact matrix, regenerated from the source matrix and community table with the labels, community boxes, B/W statistic, p-value, and log contact color scale preserved. |
| `slide_10b_mantel_exclusions_clarity.png` | Replacement Mantel exclusion plot with explicit full-arm-set x-axis, acrocentric+sex-excluded y-axis, and an on-plot explanation that points above `x=y` mean the signal increases after exclusions. |
| `slide_11_single_cell_purpose_candidate.png` | Simpler single-cell purpose slide: Fig. 3c content re-expressed as a readable per-cell test with the negative-control `S_all` result kept visible. |
| `slide_12_mouse_zygotene_trajectory_pairing.png` | Combined mouse zygotene scatter plus smaller stage-trajectory inset. Intended to replace the separate 12a/12b pair or serve as the source for a merged layout. |
| `make_hic_visual_redesign.R` | Reproducible generator for the four assets above. |

## Recommendation summary

1. Slide 10a: replace the stretched 4:3 panel crop with the square 1:1 matrix candidate.
2. Slide 10b: replace the full ED5 image with the clarified focused Mantel exclusion plot.
3. Slide 11: keep the current Fig. 3c concept, but retitle it around purpose. Do not split into a separate Mantel slide unless the deck gains another 3D-methods slot.
4. Slides 12a/12b: merge the mouse zygotene scatter and compact trajectory into one paired visual. If integration must retain two pages, use the same title/caption language below and shrink the trajectory typography substantially.

## Exact slide replacements

The paths below are relative to `slides/v2-review-zoom/_typst/`.

### Slide 10a

Current source lines: `slides/v2-review-zoom/_typst/zoom_review_deck.typ:251`.

Replace:

```typst
#figure-slide(
  "10a",
  "Fig 3 panel a - Hi-C/Pore-C contact matrix",
  "assets/s10_fig3_panel_a.png",
  source: "focused panel crop copied into this deck from local zoom assets; source figure s10_fig3.png",
)
```

With:

```typst
#figure-slide(
  "10a",
  "HG002 Pore-C contact matrix - square 1:1 community layout",
  "../_revision_assets/hic_visual_redesign/slide_10a_square_matrix_candidate.png",
  source: "candidate regenerated from hg002_porec_contact_matrix.tsv and hg002_porec_hic.arm-leiden.communities.tsv",
)
```

Title text:

`HG002 Pore-C contact matrix - square 1:1 community layout`

Caption / spoken line:

`77 arm-haplotypes are ordered by sequence community; diagonal community blocks are the physical-contact signal. B/W = 0.056 and p = 3.9e-85, with the original log contact scale preserved.`

Design note:

The current asset is `1600 x 1200`, so it can read as a stretched rectangular matrix when placed full-bleed on a 16:9 slide. The candidate is `1800 x 1800`; the matrix plotting region is square and the colorbar remains outside the matrix.

### Slide 10b

Current source lines: `slides/v2-review-zoom/_typst/zoom_review_deck.typ:260`.

Replace:

```typst
#figure-slide(
  "10b",
  "Mantel exclusions",
  "assets/s10b_ed5.png",
  source: "canonical review asset: s10b_ed5.png",
)
```

With:

```typst
#figure-slide(
  "10b",
  "Mantel exclusions: signal increases after acrocentric + sex-arm removal",
  "../_revision_assets/hic_visual_redesign/slide_10b_mantel_exclusions_clarity.png",
  source: "candidate regenerated from community_based/50000bp and no_acrocentric/50000bp global_test.tsv files",
)
```

Title text:

`Mantel exclusions: signal increases after acrocentric + sex-arm removal`

Caption / spoken line:

`Each point is one HPRC sample at 50 kb. X is the Mantel rho for the full arm set; Y is the same test after excluding acrocentric and sex arms. Points above x=y mean the sequence-similarity-by-contact correlation is stronger after the likely confounds are removed. HG002 CiFi is not in this exclusion run, so say n = 7, not 8/8.`

Design note:

The current slide title "Mantel exclusions" does not tell the audience how to read the diagonal. The replacement makes the interpretation visual: above `x=y` means the exclusion has increased rho, i.e. the evidence is not being driven by acrocentric/nucleolar or sex-arm contacts.

### Slide 11

Current source lines: `slides/v2-review-zoom/_typst/zoom_review_deck.typ:269`.

Replace:

```typst
#figure-slide(
  "11",
  "Fig 3 panel c - single-cell 3D",
  "assets/s10_fig3_panel_c.png",
  source: "focused panel crop copied into this deck from local zoom assets; source figure s10_fig3.png",
)
```

Preferred replacement:

```typst
#figure-slide(
  "11",
  "Single-cell 3D tests whether the bulk signal is per-cell",
  "../_revision_assets/hic_visual_redesign/slide_11_single_cell_purpose_candidate.png",
  source: "candidate regenerated from Dip-C and sperm per-cell community enrichment TSVs",
)
```

Title text:

`Single-cell 3D tests whether the bulk signal is per-cell`

Subtitle / caption:

`C-community arms are closer within individual GM12878 Dip-C and sperm nuclei; the zero-sharing S_all control goes the other way. This rules out the simplest bulk-average and chromosome-territory explanations.`

Recommendation:

Keep current Fig. 3c's conceptual content. The problem is not the underlying result; it is that the slide title does not state the purpose. A separate simple Mantel slide would add another statistic but would weaken the flow unless the deck has room for an additional 3D-methods explainer. The recommended fix is to keep the per-cell C-community versus `S_all` comparison and add a short purpose title/subtitle.

### Slides 12a and 12b

Current source lines: `slides/v2-review-zoom/_typst/zoom_review_deck.typ:278` and `slides/v2-review-zoom/_typst/zoom_review_deck.typ:287`.

Current 12a:

```typst
#figure-slide(
  "12a",
  "Fig 4 panel d - mouse meiotic Hi-C",
  "assets/s12_fig4_panel_d.png",
  source: "focused panel crop copied into this deck from local zoom assets; source figure s12_fig4.png",
)
```

Current 12b:

```typst
#figure-slide(
  "12b",
  "Meiotic stage trajectory inset",
  "assets/s12_trajectory.png",
  source: "canonical review asset: s12_trajectory.png",
)
```

Preferred replacement is one merged slide:

```typst
#figure-slide(
  "12",
  "Mouse zygotene: the bouquet-stage 3D signal",
  "../_revision_assets/hic_visual_redesign/slide_12_mouse_zygotene_trajectory_pairing.png",
  source: "candidate regenerated from zuo2021_zygotene_phr_pair_correlation.tsv plus the four-stage Mantel rho trajectory from slide_12_stage_trajectory.R",
)
```

Title text:

`Mouse zygotene: the bouquet-stage 3D signal`

Subtitle / caption:

`At zygotene, similar mouse subtelomeres contact more often (Spearman rho = 0.715, p = 4.4e-55, n = 344 inter-chromosomal pairs). The compact trajectory explains why zygotene is the focal stage: Mantel rho peaks at zygotene, the bouquet stage when telomeres cluster at the nuclear envelope.`

Fallback if two slides must remain:

- 12a title: `Mouse zygotene: similar subtelomeres contact more`
- 12a caption: `Show the zygotene scatter large; add one subtitle sentence that this is the bouquet stage.`
- 12b title: `Stage trajectory: zygotene is the peak`
- 12b caption: `Use the compact trajectory only, with font sizes matching an inset, not a standalone headline figure.`

Design note:

The current 12b trajectory asset is useful but typographically dominant when used as a full slide. The candidate places the trajectory as supporting evidence beside the scatter, with a smaller title and labels so it explains the zygotene choice without competing with panel d.

## Git and source provenance

The upstream provenance artifact was read first:

- `slides/v2-review-zoom/_revision_assets/git_provenance/README.md`
- Dependency log reports commit `a2ac7d1` for `review-zoom-git-provenance-audit`.

Additional git commands run for this task included:

```bash
git log --oneline --decorate --all -- paper_prep/figures/fig3 paper_prep/figures/ed5 paper_prep/figures/fig4 slides/v2/_typst/slide_12_stage_trajectory.R
git log --oneline --decorate --all -- slides/v2-review-zoom/_typst/assets/s10_fig3_panel_a.png slides/v2-review-zoom/_typst/assets/s10b_ed5.png slides/v2-review-zoom/_typst/assets/s10_fig3_panel_c.png slides/v2-review-zoom/_typst/assets/s12_fig4_panel_d.png slides/v2-review-zoom/_typst/assets/s12_trajectory.png
git show --stat --oneline --decorate 8a52549 -- paper_prep/figures/fig3
git show --stat --oneline --decorate 09f6e50 -- paper_prep/figures/ed5
git show --stat --oneline --decorate 4a1ee16 -- paper_prep/figures/fig4
git show --stat --oneline --decorate 28c2337 -- slides/v2/_typst/slide_12_stage_trajectory.R
git show --stat --oneline --decorate cb973ab -- slides/v2-zoom/_typst/assets/slide_10_fig3_panel_a.png slides/v2-zoom/_typst/assets/slide_11_fig3_panel_c.png slides/v2-zoom/_typst/assets/slide_12_fig4_panel_d.png
git show --stat --oneline --decorate 10bee88 -- slides/v2-review-zoom/_typst/assets/s10_fig3_panel_a.png slides/v2-review-zoom/_typst/assets/s10b_ed5.png slides/v2-review-zoom/_typst/assets/s10_fig3_panel_c.png slides/v2-review-zoom/_typst/assets/s12_fig4_panel_d.png slides/v2-review-zoom/_typst/assets/s12_trajectory.png
git blame --line-porcelain -- paper_prep/figures/fig3/figure_fig3.R
git blame --line-porcelain -- paper_prep/figures/ed5/figure_ed5.R
git blame --line-porcelain -- paper_prep/figures/fig4/figure_fig4.R
git blame --line-porcelain -- slides/v2/_typst/slide_12_stage_trajectory.R
```

### Commit lineage

| Item | Entered repo | Later changes relevant here |
|---|---|---|
| Fig. 3 source script and full PNG/PDF | `8a52549` (`feat: figure-3-3d-convergence`, agent-700) added `paper_prep/figures/fig3/figure_fig3.R`, caption, sources, PDF, and PNG. `git blame` attributes the inspected script lines to this commit. | Agent-878 `cb973ab` created zoom crops `slide_10_fig3_panel_a.png` and `slide_11_fig3_panel_c.png`; agent-951 `10bee88` copied them into the current review-zoom deck as `s10_fig3_panel_a.png` and `s10_fig3_panel_c.png`. Squash merges are `fd9a250` and `4862ec7`. |
| Extended Data Fig. 5 source script and full PNG/PDF | `09f6e50` (`feat: figure-ed5-hic-robustness`, agent-702) added `paper_prep/figures/ed5/figure_ed5.R`, caption, sources, PDF, and PNG. `git blame` attributes the inspected script lines to this commit. | Agent-951 `10bee88` copied `figure_ed5.png` into the current review-zoom deck as `s10b_ed5.png`; blob hash matches the source PNG. |
| Fig. 4 source script and full PNG/PDF | `4a1ee16` (`feat: figure-4-pedigree-mouse`, agent-693) added `paper_prep/figures/fig4/figure_fig4.R`, caption, sources, PDF, and PNG. `git blame` attributes the inspected script lines to this commit. | Agent-878 `cb973ab` created `slide_12_fig4_panel_d.png`; agent-951 `10bee88` copied it into the current review-zoom deck as `s12_fig4_panel_d.png`. |
| Stage trajectory inset | `28c2337` (`feat: build-bog-v2-2`, agent-813) added `slides/v2/_typst/slide_12_stage_trajectory.R` plus PNG/PDF outputs. `git blame` attributes the trajectory script to this commit. | Agent-951 `10bee88` copied the PNG into the current review-zoom deck as `s12_trajectory.png`. |
| v2 slide narratives | `0e76f59` / branch commit `bca2237` added `slide_10_hic_bulk_mantel_exclusions.md`; `b784e3e` / branch commit `1255c69` added `slide_11_single_cell_3d.md`; `0772fbe` / branch commit `964ce4c` added `slide_12_mouse_meiotic_zygotene_bouquet.md`. | Used as conceptual input for titles and captions. |

### Blob identity checks

Exact matches observed with `git hash-object`:

| Current review-zoom asset | Matching source | Blob hash |
|---|---|---|
| `slides/v2-review-zoom/_typst/assets/s10_fig3_panel_a.png` | `slides/v2-zoom/_typst/assets/slide_10_fig3_panel_a.png` | `adedf4bc3931134dc99c7e1164acd4599fae5720` |
| `slides/v2-review-zoom/_typst/assets/s10_fig3_panel_c.png` | `slides/v2-zoom/_typst/assets/slide_11_fig3_panel_c.png` | `b1d3ab4b6c9ace4a3f34ce5b1bd6758213f08d65` |
| `slides/v2-review-zoom/_typst/assets/s12_fig4_panel_d.png` | `slides/v2-zoom/_typst/assets/slide_12_fig4_panel_d.png` | `d1b70e2f4d51d24cef76c43a356230add55d9d4c` |
| `slides/v2-review-zoom/_typst/assets/s10b_ed5.png` | `paper_prep/figures/ed5/figure_ed5.png` | `b299dc2b6f149e8730baf7ebe3fc7fc776b7abe0` |
| `slides/v2-review-zoom/_typst/assets/s12_trajectory.png` | `slides/v2/slide_12_stage_trajectory.png` | `98dec76c3cd639d526fde61ad7ccd84350d417ac` |
| `slides/v2-review-zoom/_typst/assets/s10_fig3.png` | `paper_prep/figures/fig3/figure_fig3.png` | `d01991044447ac951cc68a199b31dc335128f0bd` |
| `slides/v2-review-zoom/_typst/assets/s12_fig4.png` | `paper_prep/figures/fig4/figure_fig4.png` | `799f9039f68ae20e4368aad96d2a0f6e75b9395b` |

## Source scripts inspected

- Fig. 3 panel a source: `paper_prep/figures/fig3/figure_fig3.R:32` starts the HG002 Pore-C contact matrix renderer; lines 34-40 read `hg002_porec_contact_matrix.tsv` and `hg002_porec_hic.arm-leiden.communities.tsv`; lines 63-80 draw community borders and the B/W/p-value label.
- Fig. 3 panel c source: `paper_prep/figures/fig3/figure_fig3.R:206` starts the `S_all` negative-control panel; lines 209-220 read GM12878 Dip-C and sperm per-cell/per-community TSVs; lines 261-274 add the below-unity counts and `S_all` interpretation.
- ED5 panel b source: `paper_prep/figures/ed5/figure_ed5.R:211` starts the Mantel exclusion panel; lines 214-218 merge full and no-acrocentric 50 kb Mantel rho values; lines 220-231 draw the identity diagonal and the full-vs-exclusion axes.
- Fig. 4 panel d source: `paper_prep/figures/fig4/figure_fig4.R:249` starts the mouse zygotene scatter; lines 251-260 filter inter-chromosomal pairs and compute Spearman rho/p-value; lines 262-278 plot the scatter and annotation.
- Stage trajectory source: `slides/v2/_typst/slide_12_stage_trajectory.R:1` describes the four-stage Mantel rho inset; lines 8-12 hard-code the stage values `0.687, 0.718, 0.683, 0.577`; lines 23-28 annotate zygotene as the bouquet.
- Current review-zoom deck references: `slides/v2-review-zoom/_typst/zoom_review_deck.typ:251`, `:260`, `:269`, `:278`, and `:287`.

## Methods glossary status

The requested Hi-C methods glossary task was checked via `wg list`; as of this work it was still open:

`review-zoom-hic-methods-glossary - Review zoom Hi-C methods: explain O/E, mcool, Mantel, Pore-C, Dip-C, mouse`

Because it was not available, this redesign used the upstream provenance audit, the v2 slide narrative files, and the source scripts listed above as conceptual input.

## Candidate generation provenance

The candidate generator reads these source data files directly:

- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_contact_matrix.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_hic.arm-leiden.communities.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/<sample>_global_test.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/no_acrocentric/50000bp/<sample>_global_test.tsv`
- `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_16cells_500kb_per_cell.tsv`
- `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_16cells_500kb_per_community_per_cell.tsv`
- `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_per_cell.tsv`
- `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_per_community_per_cell.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp/zuo2021_zygotene_phr_pair_correlation.tsv`

These are off-tree analysis paths, so final integration should either retain this generator as provenance or copy in a stable data snapshot/checksum if the deck must be fully reproducible without `/moosefs/guarracino`.

## Validation

- `Rscript slides/v2-review-zoom/_revision_assets/hic_visual_redesign/make_hic_visual_redesign.R` completed successfully.
- `file slides/v2-review-zoom/_revision_assets/hic_visual_redesign/*.png` confirmed:
  - `slide_10a_square_matrix_candidate.png`: `1800 x 1800`
  - `slide_10b_mantel_exclusions_clarity.png`: `1800 x 1200`
  - `slide_11_single_cell_purpose_candidate.png`: `1800 x 1200`
  - `slide_12_mouse_zygotene_trajectory_pairing.png`: `1800 x 1200`
- Visual checks were performed on all four generated PNGs.
