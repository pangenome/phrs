# Hi-C / Pore-C / Dip-C / Mouse Methods Glossary

Task: `review-zoom-hic-methods-glossary`

This note is for the review-zoom 3D evidence block: slides 10a, 10b,
11, 12a, and 12b. It does not edit deck source. It explains what the
current panels actually plot, which numbers belong to which panel, and
which labels need tightening before final slide assembly.

## One-screen takeaway

- Slide 10a is **HG002 Pore-C**, not generic Hi-C: a 77 x 77 arm-haplotype
  contact matrix at 50 kb, ordered by sequence community. Red diagonal
  blocks mean same-community arms contact more often.
- Slide 10b currently shows the **full Extended Data Fig 5 PNG**, even
  though the slide title says "Mantel exclusions". The Mantel exclusion
  plot is only panel b of that image: x = full-arm Mantel rho, y =
  no-acrocentric/sex Mantel rho, one point per sample/technology.
- Slide 11 is **single-cell 3D distance**, not contact frequency. The y-axis
  is W/B distance: within-community distance divided by between-community
  distance. Values below 1 mean same-community arms are closer.
- Slide 12a is **mouse zygotene per-PHR-pair Spearman**, not the stage
  trajectory. Each point is one inter-chromosomal mouse arm pair; x =
  subtelomeric sequence similarity, y = zygotene Hi-C contact.
- Slide 12b is the **meiotic-stage Mantel trajectory**. It is the evidence
  that zygotene is a peak, not just a hand-picked stage.

## Source hierarchy

Use this order when resolving conflicts:

1. Figure scripts and `sources.tsv`: exact plotted data, order, axes, and
   transformations.
2. The absolute source TSVs listed in `sources.tsv`: exact numeric values.
3. Figure captions: intended figure-level interpretation.
4. Survey/crosswalk notes: method background, caveats, and talk framing.
5. Slide notes and deck asset names: useful for narrative, but not always
   exact about what a specific plotted crop shows.

Authoritative local source files inspected:

- [Fig 3 script](../../../../paper_prep/figures/fig3/figure_fig3.R) and
  [caption](../../../../paper_prep/figures/fig3/caption.md), with
  [sources](../../../../paper_prep/figures/fig3/sources.tsv).
- [ED5 script](../../../../paper_prep/figures/ed5/figure_ed5.R) and
  [caption](../../../../paper_prep/figures/ed5/caption.md), with
  [sources](../../../../paper_prep/figures/ed5/sources.tsv).
- [Fig 4 script](../../../../paper_prep/figures/fig4/figure_fig4.R) and
  [caption](../../../../paper_prep/figures/fig4/caption.md), with
  [sources](../../../../paper_prep/figures/fig4/sources.tsv).
- [Survey 05](../../../../paper_prep/surveys/SURVEY_05_hic_validation.md)
  for Hi-C/Pore-C/CiFi methods, mcool resolutions, O/E, Mantel, and
  exclusions.
- [Survey 06](../../../../paper_prep/surveys/SURVEY_06_dipc_validation.md)
  for Dip-C, sperm scHi-C, S_all, and W/B distance.
- [Survey 07](../../../../paper_prep/surveys/SURVEY_07_integrated.md) for
  the convergent-evidence and meiotic-bouquet interpretation.
- [Survey 08](../../../../paper_prep/surveys/SURVEY_08_mouse.md) for mouse
  T2T, zygotene Hi-C, stage trajectory, windows, and caveats.
- [CROSSWALK](../../../../paper_prep/synthesis/CROSSWALK.md) for manuscript
  scope and the note that mouse is best framed as a talk/SI meiotic bridge,
  not necessarily a canonical main-text claim.
- [Git provenance note](../git_provenance/README.md) for asset lineage and
  the missing crop recipes for slides 10a, 11, and 12a.

## Glossary

### PHR

PHR means the projected homologous region detected in subtelomeric flanks.
For the 3D panels, Hi-C/Pore-C contact is extracted at PHR-specific
coordinates, not at arbitrary whole-chromosome windows, except where a
fixed fallback window is explicitly stated for Dip-C/S_all controls.

### Arm-haplotype

An arm-haplotype is a chromosome arm kept separately by haplotype, such as
`chr4_MATERNAL_q` or `chr4_PATERNAL_q`. Diploid Hi-C/Pore-C samples keep
maternal and paternal arms separate; CHM13 is haploid. Fig 3A uses 77
HG002 arm-haplotypes present in the Pore-C matrix.

### Sequence community

The C1...C15 communities are Leiden communities discovered from
subtelomeric sequence similarity. In the 3D block, the core question is
whether these sequence-defined groups also behave as physical groups in
Hi-C/Pore-C/Dip-C/mouse 3D data.

### B/W contact ratio

For bulk contact maps, the manuscript/talk usually reports:

`B/W = mean between-community contact / mean within-community contact`

B/W < 1 means within-community contacts are stronger. For example, Fig 3A
uses `hg002_porec_global_test.tsv`: within mean = 0.02746, between mean =
0.001525, so B/W = 0.001525 / 0.02746 = 0.056, with p = 3.9e-85.

### W/B distance ratio

For Dip-C and sperm single-cell 3D distance, the direction is inverted:

`W/B = mean within-community 3D distance / mean between-community 3D distance`

W/B < 1 means within-community arms are closer in Euclidean 3D space. Do
not call these "contacts"; they are distances derived from 3D coordinates.

### O/E contact

O/E means observed/expected contact for an arm pair. In this repo's ED5
implementation, expected contact is computed from inter-arm matrix
marginals:

`E_ij = row_sum_i * col_sum_j / total_inter`

The O/E entry is observed contact divided by that expected value. This
controls for arm-level contactability, chromosome-size effects, and
marginal mappability. ED5C then partitions O/E values into within-community
and between-community trans arm pairs, excluding cis pairs from the mean.

Important distinction: Fig 3A is a raw/log contact matrix, not an O/E
matrix. O/E appears in ED5C and in the independent Hi-C Leiden/ARI workflow.

### .mcool / "5 mcool"

`.mcool` is a multi-resolution Cooler file: an HDF5 container holding the
same contact map binned at multiple resolutions. In this project, "5 mcool"
means the five queried resolutions, not five separate biological datasets:

`5 kb, 10 kb, 20 kb, 50 kb, 100 kb`

The human Hi-C/Pore-C/CiFi survey says each sample's `.mcool` contact
matrix is loaded at all five resolutions, then PHR intervals are queried at
each resolution. Example `.mcool` inputs include:

- `/moosefs/guarracino/HPRCv2/PHR_III/HiC/HG002/hicpro_output/cool/HG002.mcool`
- `/moosefs/guarracino/HPRCv2/PHR_III/HiC/HG002/porec_output/porec.mcool`
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/HiC/zuo2021_zygotene.mcool`

### Mantel rho

In the bulk Hi-C/Pore-C panels, Mantel rho is the Spearman correlation
between two arm x arm matrices:

- sequence similarity or distance from subtelomeric Jaccard values
- 3D contact or distance from Hi-C/Pore-C/Dip-C

The global-test TSVs store Mantel rho in the `U_statistic` column on the
row where `test == "mantel"`. It is not a pointwise raw-contact scatter.
For the ED5B plot, each plotted point is one sample/technology, not one arm
pair.

### Pore-C

Pore-C is the ONT multi-way contact assay used here as an independent
platform from standard pairwise Hi-C. Fig 3A is HG002 Pore-C. The matrix is
rendered as a contact matrix after Pore-C processing, but the slide should
not let the audience think it is a conventional Hi-C-only panel.

### Dip-C

Dip-C reconstructs single-cell diploid 3D genome coordinates. Here, 17
GM12878 cells were remapped to T2T-CHM13v2.0; one duplicate was excluded,
leaving 16 cells. The analysis samples PHR or fallback terminal-window
coordinates and compares Euclidean 3D distances within and between
sequence communities.

### S_all

S_all is a negative-control pseudo-community made from seven arms with zero
inter-chromosomal subtelomeric sequence sharing:

`chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q`

In slide 11, C-community arms are closer in 3D, but S_all is farther:
GM12878 S_all W/B = 1.106 and sperm S_all W/B = 1.397. This is the best
quick rebuttal to "any arbitrary terminal arms cluster."

### Zygotene

Zygotene is the meiotic prophase stage when telomeres are clustered at the
nuclear envelope in the bouquet and homolog alignment is active. In the
mouse data, the sequence-similarity/contact relationship peaks at zygotene:
per-PHR-pair Spearman rho = 0.715, p = 4.4e-55, n = 344 inter-chromosomal
pairs. The companion trajectory uses Mantel rho values:

`leptotene 0.687, zygotene 0.718, pachytene 0.683, diplotene 0.577`

## Slide 10a: Fig 3 panel A

### What the panel shows

Rows and columns are 77 HG002 arm-haplotypes from the Pore-C contact
matrix:

`/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_contact_matrix.tsv`

The same arm order is used on both axes. The figure script matches matrix
rows to:

`hg002_porec_hic.arm-leiden.communities.tsv`

and then orders arms by `community`, then `arm`. The blue boxes are
community blocks. Each matrix cell is an arm-pair Pore-C contact value at
50 kb. The visible red diagonal blocks are not chromosome diagonals; they
are sequence-community blocks.

### Platform and resolution

- Platform: HG002 Pore-C, processed as an arm-level contact matrix.
- Resolution: 50 kb.
- Unit: contact frequency/contact value at PHR-derived regions.
- Population: HG002 diploid arm-haplotypes; maternal/paternal arms remain
  separate.

### Color scale

The script plots `log10(contact)`. Zero contacts are replaced with a small
floor, defined as the 5th percentile of positive contacts, so zeros can be
drawn rather than dropped. Redder means higher contact. The matrix is not
O/E-normalized.

### Number to say

"This is HG002 Pore-C at 50 kb. The sequence-community blocks are on both
axes. Same-community arms have much higher contact: B/W = 0.056, p =
3.9e-85."

If translating B/W verbally: B/W = 0.056 means within-community contact is
about 18x the between-community contact.

### Source references

- Script: `paper_prep/figures/fig3/figure_fig3.R`, panel A.
- Caption: `paper_prep/figures/fig3/caption.md`, Fig 3A.
- Sources: `paper_prep/figures/fig3/sources.tsv`, rows `3a`.

## Slide 10b: ED5 / Mantel exclusions

### What the current zoom asset contains

The current deck asset `slides/v2-review-zoom/_typst/assets/s10b_ed5.png`
is the full ED5 image, not a focused crop of panel b. Therefore the slide
contains four distinct panels:

- ED5A: within/between contact across five mcool resolutions.
- ED5B: Mantel rho full arm set vs no-acrocentric/sex arm set.
- ED5C: O/E within vs between sequence community.
- ED5D: per-community enrichment heatmap across datasets.

If the slide title stays "Mantel exclusions", crop to ED5B or retitle the
slide to reflect the full robustness figure.

### What ED5B shows

Each point is one HPRC sample/technology at 50 kb. HG002 CiFi is absent
from this panel because it was not run for this control. The axes are:

- x-axis: Mantel rho from the full arm set,
  `community_based/50000bp/<sample>_global_test.tsv`
- y-axis: Mantel rho after the no-acrocentric/sex exclusion,
  `no_acrocentric/50000bp/<sample>_global_test.tsv`

The plotted rho is read from the row `test == "mantel"`, column
`U_statistic`. The diagonal line is y = x. A point above the line means the
sequence-similarity/contact correlation is stronger after the exclusion.

Plotted ED5B values from the actual source TSVs:

| Sample | Full rho | ED5B no-acro/sex rho | Interpretation |
|---|---:|---:|---|
| CHM13 Hi-C | 0.656 | 0.799 | stronger |
| HG002 Hi-C | 0.657 | 0.742 | stronger |
| HG002 Pore-C | 0.486 | 0.554 | stronger |
| HG00658 Hi-C | 0.276 | 0.368 | stronger |
| HG02148 Hi-C | 0.152 | 0.209 | stronger, still weak/borderline |
| HG02559 Hi-C | 0.397 | 0.500 | stronger |
| NA19036 Hi-C | 0.266 | 0.486 | stronger |

### Speaker-ready explanation

"This is not a raw arm-pair scatter; each point is one dataset. The x value
is the Mantel correlation between subtelomeric sequence similarity and
bulk contact using all arms. The y value repeats the same test after
dropping the acrocentric/sex-chromosome confound set. Every point moves
above the identity line, so the signal is not created by acrocentric
nucleolar clustering or pseudoautosomal sex-chromosome contacts."

### Do not mix these two exclusion sources

Some slide notes and surveys quote a broader exclusion table with stronger
headline numbers, for example HG002 0.657 -> 0.790, CHM13 0.656 -> 0.850,
HG02148 0.152 -> 0.720. Those values come from the multi-exclusion summary
in `end-to-end-report/report/05_hic_validation.md`, especially the
`no acro pq + sex` / `no strong` columns. They are not the exact values in
the plotted ED5B panel, which reads the `analysis/human/no_acrocentric/`
TSVs.

If Erik points at the current ED5B plot, use the ED5B values above. If the
final slide wants the stronger 0.79/0.85/0.72 story, it needs a different
plot or a label that explicitly says it is using the broader exclusion
table, not ED5B as currently rendered.

### Source references

- Script: `paper_prep/figures/ed5/figure_ed5.R`, panel B.
- Caption: `paper_prep/figures/ed5/caption.md`, ED5B.
- Sources: `paper_prep/figures/ed5/sources.tsv`, row `ed5b`.
- Method survey: `paper_prep/surveys/SURVEY_05_hic_validation.md`.

## Slide 11: Fig 3 panel C single-cell 3D

### What the panel shows

The panel compares per-cell distance ratios for two cell systems and two
arm sets:

- GM12878 Dip-C, C-community arms.
- GM12878 Dip-C, S_all zero-sharing pseudo-community.
- Sperm scHi-C, C-community arms.
- Sperm scHi-C, S_all zero-sharing pseudo-community.

Each point is one cell. The y-axis is per-cell W/B distance:

`mean within-community 3D distance / mean between-community 3D distance`

The dashed horizontal line at 1 means no difference. Values below 1 mean
within-community arms are closer.

### Exact slide logic

- GM12878 C-community: 16/16 cells have W/B < 1.
- GM12878 S_all: 0/16 cells have W/B < 1.
- Sperm C-community: 20/20 cells have W/B < 1.
- Sperm S_all: 1/20 cells has W/B < 1.

The intended message is not "large single-cell rho". It is that the
within-community proximity signal survives the transition from bulk
Hi-C/Pore-C to single-cell 3D coordinates, and the zero-sharing control
reverses the direction.

### Speaker-ready explanation

"Bulk Hi-C averages many cells, so the next question is whether individual
3D genomes show the same effect. Here the y-axis is distance, not contact.
Below 1 means same-community arms are physically closer. In GM12878 Dip-C,
all 16 cells put C-community arms closer. In sperm, all 20 do. But the
seven zero-sharing arms pooled as S_all go the opposite direction: they are
11 percent farther in GM12878 and 40 percent farther in sperm. So the
single-cell signal follows sequence sharing, not just terminal-arm crowding."

### Source references

- Script: `paper_prep/figures/fig3/figure_fig3.R`, panel C.
- Caption: `paper_prep/figures/fig3/caption.md`, Fig 3C.
- Sources: `paper_prep/figures/fig3/sources.tsv`, rows `3c`.
- Method survey: `paper_prep/surveys/SURVEY_06_dipc_validation.md`.

## Slide 12a: Fig 4 panel D mouse meiotic Hi-C

### What the panel shows

Fig 4D is a mouse zygotene scatter:

- x-axis: mean Jaccard similarity for the PHR pair.
- y-axis: Hi-C contact in zygotene, 50 kb, plotted on a log scale.
- each point: one inter-chromosomal mouse arm pair after filtering to
  positive contact.
- source rows: `zuo2021_zygotene_phr_pair_correlation.tsv`.

The figure script recomputes Spearman correlation over the plotted points:

`rho = 0.715, p = 4.4e-55, n = 344 inter-chromosomal pairs`

This is per-PHR-pair Spearman, not the same statistic as the Mantel
trajectory on slide 12b.

### Platform, resolution, and genome

- Genome: mouse T2T B6 + CAST assemblies from Francis et al. 2025.
- Hi-C: Zuo et al. 2021 mouse meiotic Hi-C.
- Stage: zygotene.
- Window: 1 Mb mouse subtelomeric window for this panel.
- Resolution: 50 kb.
- Analysis mode: arm-level / per-PHR-pair sequence similarity vs Hi-C
  contact; inter-chromosomal pairs only.

### Why zygotene matters

Human bulk Hi-C is somatic; the recombination mechanism is meiotic. Mouse
zygotene matters because it is the bouquet stage: telomeres are clustered
at the nuclear envelope while homologs align. The zygotene peak provides a
stage-specific mechanistic bridge from "subtelomeric sequence similarity
predicts contact" to "the relevant meiotic architecture can put homologous
subtelomeres together."

Important caveat: this is mouse, not human meiotic Hi-C. It supports the
meiotic interpretation but should not be phrased as direct human meiotic
proof.

### Speaker-ready explanation

"This is the meiotic bridge. In zygotene mouse spermatocytes, each point is
one inter-chromosomal subtelomeric arm pair. More similar PHRs have higher
Hi-C contact, with Spearman rho 0.715 and p = 4.4e-55. Zygotene is the
bouquet stage, when telomeres are concentrated at the nuclear envelope, so
the same sequence-to-3D relationship peaks exactly where a meiotic
recombination model predicts it should."

### Source references

- Script: `paper_prep/figures/fig4/figure_fig4.R`, panel D.
- Caption: `paper_prep/figures/fig4/caption.md`, Fig 4D.
- Sources: `paper_prep/figures/fig4/sources.tsv`, rows `4d`.
- Method survey: `paper_prep/surveys/SURVEY_08_mouse.md`.

## Slide 12b: meiotic stage trajectory inset

### What the panel shows

The trajectory inset plots Mantel rho across four mouse meiotic prophase
stages at 50 kb, 1 Mb window:

| Stage | Mantel rho |
|---|---:|
| Leptotene | 0.687 |
| Zygotene | 0.718 |
| Pachytene | 0.683 |
| Diplotene | 0.577 |

The R script hard-codes these four values from the mouse report/survey and
highlights zygotene as the bouquet peak. The inset is load-bearing because
Fig 4D alone only shows the zygotene scatter; it does not prove zygotene is
the stage maximum.

### Speaker-ready explanation

"The inset is the reason zygotene is not cherry-picked. Run the same
Mantel comparison across prophase: leptotene 0.687, zygotene 0.718,
pachytene 0.683, diplotene 0.577. The contact-similarity relationship rises
into zygotene and falls as the bouquet resolves. That stage profile is the
biological point."

### Source references

- R script: `slides/v2/_typst/slide_12_stage_trajectory.R`.
- Slide source note: `slides/v2/slide_12_mouse_meiotic_zygotene_bouquet.md`.
- Mouse survey: `paper_prep/surveys/SURVEY_08_mouse.md`, sections 1.7 and
  7.
- CROSSWALK row: `paper_prep/synthesis/CROSSWALK.md`, chapter 08 row.

## Ambiguous labels to fix before final slides

1. **Slide 10a title should not say generic "Hi-C/Pore-C" if the only
   panel is Fig 3A.** The plotted panel is HG002 Pore-C. Suggested title:
   "HG002 Pore-C contacts, 50 kb, ordered by sequence community."

2. **Slide 10b title and asset mismatch.** The deck title says "Mantel
   exclusions", but the asset is the full ED5 PNG. Either crop panel b with
   a committed crop recipe, or retitle as "Hi-C robustness: resolution,
   exclusions, O/E, community reproducibility."

3. **Do not quote 0.79/0.85/0.72 while pointing at ED5B.** ED5B plotted
   values are HG002 0.657 -> 0.742, CHM13 0.656 -> 0.799, HG02148 0.152 ->
   0.209. The stronger values come from a broader exclusion table, not the
   current ED5B plot.

4. **Clarify the ED5B exclusion set.** The ED5 caption says
   acrocentric+sex exclusion and the `no_acrocentric` path is used. The arm
   count drop is consistent with acrocentric p arms plus sex arms, not the
   broader acro p+q+sex table. Final labels should say exactly which arms
   were excluded.

5. **ED5A ratio convention is confusing.** Fig 3/survey talk often says
   B/W < 1, but ED5A plots within/between contact > 1. Both mean the same
   biology, but the label must make the inversion explicit.

6. **Slide 11 filename starts with `s10_`, but it is slide 11 / Fig 3C.**
   This is only a lineage/naming issue, but notes should call it "Fig 3C
   single-cell 3D" to avoid confusing it with slide 10.

7. **Do not call slide 12a's rho "Mantel 0.715".** Fig 4D's 0.715 is
   per-PHR-pair Spearman on the plotted scatter. Slide 12b's zygotene
   Mantel value is 0.718.

8. **Mouse framing should be "meiotic bridge", not direct human proof.**
   The mouse panel is powerful because it is meiotic and stage-resolved,
   but it is cross-species. Say "supports the meiotic interpretation" or
   "mouse meiotic bridge", not "proves human meiotic proximity."

## Per-slide pocket notes

### 10a

"This is HG002 Pore-C at 50 kb. Rows and columns are 77 arm-haplotypes,
ordered by their sequence community. The red blocks are same-community
arms with elevated contact. B/W = 0.056, so within-community contact is
about 18x between-community, p = 3.9e-85."

### 10b

"Each point in panel b is a dataset, not an arm pair. x is the full-arm
Mantel rho; y repeats the test after the no-acrocentric/sex exclusion.
Every point is above the diagonal, so removing the obvious acrocentric and
sex-chromosome confounds strengthens rather than erases the correlation."

### 11

"This y-axis is distance. Below 1 means closer. C-community arms are closer
in every GM12878 Dip-C cell and every sperm cell, while the seven
zero-sharing S_all arms are farther. This is the single-cell and negative
control version of the bulk contact claim."

### 12a

"Each point is one mouse inter-chromosomal PHR arm pair at zygotene:
sequence similarity on x, Hi-C contact on y. The correlation is rho =
0.715 with p = 4.4e-55. Zygotene is the bouquet stage, so this is the
meiotic architecture version of the sequence-to-3D relationship."

### 12b

"The inset prevents cherry-picking. Mantel rho is 0.687, 0.718, 0.683,
0.577 from leptotene through diplotene. It peaks at zygotene, exactly when
telomeres are clustered at the nuclear envelope."
