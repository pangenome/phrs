# HPRCv2 Karyogram Candidate For Review Zoom Slide 04

Task: `review-zoom-04-hprcv2-karyogram`

## Candidate Asset

- Local candidate PNG: `p_interchrom_karyogram_count_rainbow_inset.100kb.png`
- Source request: <https://github.com/pangenome/HPRCv2/blob/main/p_interchrom_karyogram_count_rainbow_inset.100kb.png>
- Stable raw source used for this copy: <https://raw.githubusercontent.com/pangenome/HPRCv2/d14883c314e683063abe8b461992f12825ccd5ed/p_interchrom_karyogram_count_rainbow_inset.100kb.png>
- HPRCv2 remote: <https://github.com/pangenome/HPRCv2.git>
- Observed `main` commit: `d14883c314e683063abe8b461992f12825ccd5ed`
- Observed related branch: `inter-chr-plot` at `ab5c9bed552ac33dd63919a3b62ab624ceba2f3f`
- Local PNG dimensions: 4800 x 2100 RGB PNG
- Local PNG git blob hash: `f4f8f7c46f00679b742308af658707f7500d2898`
- Local PNG SHA-256: `858198c40119bef0b696a056fbea9f4bdc9af5649b51463b0903e3b7019f5689`

Upstream provenance from `scripts/plot-impg-coverage.inter-chr-map.R` at
`d14883c314e683063abe8b461992f12825ccd5ed`:

- The script reads `hprc25272-wf.CHM13.100kb-xm5-id098-l50000.tsv.gz` from the
  HPRCv2 data directory path in that repository script.
- The plotted data are 100 kb CHM13 windows for 466 HPRCv2 haplotypes.
- The main karyogram colors each CHM13 window by `num_other_chroms`, after
  zeroing out self-chromosome mappings.
- The rainbow version bins `num_other_chroms` as `0`, `1`, `2`, `3`, `4`, `5`,
  `6-10`, and `11+`, using a grey zero bin plus a turbo/viridis hue scale.
- The saved candidate is produced by the script as
  `p_interchrom_karyogram_count_rainbow_inset.100kb.png` with `ggsave(width =
  16, height = 7, dpi = 300, bg = "white")`.
- The insets are not generic legends. They zoom specific source-target count
  patterns at chr1 start, chr4q end, and chr13 p arm.

## Current Zoom Deck Comparison

Current zoom slide references are in
`slides/v2-review-zoom/_typst/zoom_review_deck.typ`; this task did not edit
that deck.

| Zoom label | Current asset | Current role | Relation to candidate |
| --- | --- | --- | --- |
| `04` | `assets/s04_fig1_panel_a.png` | Fig 1 panel a genome-wide identity heatmap, 1650 x 1600. Shows maximum alignment identity to another chromosome in 100 kb windows and includes the chr18q inset from the manuscript figure. | Candidate is not the same measurement. It shows count of other chromosomes with interchromosomal mappings, not maximum identity. It can replace the slide only if the slide title/caption/speaker framing also change from "identity heatmap" to "interchrom karyogram/count landscape". |
| `05` | `assets/s05_interchrom.png` | Interchromosomal similarities / number of chromosomes per region, 2400 x 1350. This is the closest current semantic match. | Candidate is a stronger HPRCv2 replacement or companion for this slide because both answer "how many chromosomes are mixing here?". Candidate adds a chromosome-by-chromosome karyogram layout and inset source-chromosome count heatmaps. |
| `07a` | `assets/s07_fig1_panel_c.png` | Fig 1 panel c arm-level all-vs-all Jaccard heatmap, 1650 x 1650. | Do not replace. It is a downstream arm-level clustering/distance panel, not a genome-wide 100 kb interchromosomal count/identity panel. It should only reference the candidate verbally as a transition from genome-wide signal to arm-level structure. |

No later zoom pages after `07a` repeat this genome-wide 100 kb interchromosomal
view. Slides `08a`, `08b`, `09a`, and `09c` are sequence-level MDS/PCA
projections; slides `10a` and later move to Hi-C/Pore-C, meiosis, pedigree, and
gene biology. They should not use this karyogram asset directly.

## Placement Recommendation

Recommended integration choice: add the HPRCv2 karyogram as a new zoom slide
between current zoom `04` and `05`, then either remove or demote the current
`05` repeated count plot in the final deck edit.

Exact recommendation:

1. Keep current `04` / `s04_fig1_panel_a.png` if the talk still needs the
   Fig 1a identity claim: "PAR2-scale high-identity interchromosomal homology".
   The requested HPRCv2 PNG is not a 1:1 replacement for that claim.
2. Add a new slide immediately after `04`, provisionally `04b` or `05a`, titled
   `HPRCv2 interchrom karyogram - # other chromosomes per 100 kb window`.
   Use this candidate PNG full-width.
3. Replace or skip current zoom `05` if deck length matters. If retained, make
   `05` a short reference/comparison slide rather than a second full count
   landscape. The HPRCv2 karyogram is the better count-landscape visual.
4. Do not replace `07a`. Its arm-level heatmap is a distinct downstream matrix
   used to introduce community/clade structure.

Fallback if Erik's request is interpreted as a hard replacement of slide `04`:

- Replace `04` with this PNG, but update the slide label/source text and speaker
  framing away from "Fig 1 panel a - genome-wide heatmap" and toward "HPRCv2
  interchrom karyogram - number of other chromosomes".
- Also remove or substantially compress current `05`, because the candidate
  already covers the genome-wide count view more directly.
- Do not present the replacement as Fig 1 panel a unless the final figure/caption
  source is also updated.

## Cropping And Scaling Guidance

- Use the PNG as a full-width image with `fit: "contain"` in the 16:9 zoom deck.
  It is already a 16 x 7 inch, 300 dpi export and is slide-ready.
- Aspect ratio is 4800:2100, or 2.286:1. The current zoom image frame is closer
  to 16:9 after the title/footer area, so full-width placement will leave some
  vertical whitespace. That is preferable to cropping because the right-side
  insets and legends carry essential interpretation.
- Do not crop the right edge. The chr1/chr4/chr13 insets are part of the
  requested asset and explain which chromosomes/haplotypes support selected
  high-signal regions.
- Do not crop the left or bottom axes. The chromosome labels and Mbp scale are
  needed for the karyogram to remain self-contained.
- If integration needs a closer conference-distance view, create a second
  derived zoom later rather than replacing this copy: one derived crop for the
  main karyogram body and one derived crop for the inset column. Keep this
  original PNG as the provenance-preserving candidate.

## Deck Edit Boundary

This task intentionally leaves `slides/v2-review-zoom/_typst/zoom_review_deck.typ`
unchanged. Final integration should copy or reference the candidate from this
revision asset directory and make any slide numbering, title, footer, and page
count changes in the fan-in deck edit.
