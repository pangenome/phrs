# Slide 06a violin/censor candidate

Task: `review-zoom-v3-slide06-violin-censor`

## What was inspected

- Current v2 rendered 06a: `slides/v2-review-zoom/_typst/page-09.png`, which uses `_revision_assets/06_length_redesign/candidate_06a_ranked_arm_summary.png`.
- Previous slide-06 fanout outputs in `slides/v2-review-zoom/_revision_assets/06_length_redesign/`, especially:
  - `candidate_06a_ranked_arm_summary.png`
  - `candidate_06a_focused_clade_facets.png`
  - `candidate_06b_clade_story_matrix.png`
  - `make_06_length_redesign.R`
  - `arm_length_summary.tsv`
  - `clade_length_summary.tsv`
- Provenance note: `slides/v2-review-zoom/_revision_assets/git_provenance/README.md`, which warns that the old 06a crop lineage came from copied/split PNGs and that slide 06 should be regenerated from source data.

## Data source

The candidate is regenerated from the canonical length TSV documented by the prior fanout and task metadata:

`/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv`

The generator reads this TSV directly. It does not copy or crop any existing PNG.

## Unit of observation

Each plotted observation is one non-empty inter-chromosomal PHR interval from the TSV. The plotted length is:

`(region_end - region_start) / 1000`

Rows where `region_start` and `region_end` are `.` are counted only for the no-signal arm audit and are not plotted in the violin distribution because they have no measured interval length.

## Grouping choice

The primary candidate groups intervals by named community/clade plus the remaining signaled background:

- `Other signaled arms`: all signaled arms not in the named clade callouts.
- `C2 10p-18p`: `10p`, `18p`.
- `C1 DUX4/D4Z4`: `4q`, `10q`.
- `C14 PAR2`: `Xq`, `Yq`.
- `C15 PAR1`: `Xp`, `Yp`.
- `C7 acrocentric p`: `13p`, `14p`, `15p`, `21p`, `22p`.

This grouping keeps the plot readable on one slide and connects 06a directly to the named community/clade vocabulary used by 06b and the later tree/MDS slides. The no-signal arms are documented in the caption and summary table rather than shown as an empty violin.

## 500 kb cap interpretation

The analysis/search window stops at `500 kb`. Values above `500 kb` are not measured in this TSV.

The candidate therefore treats `500 kb` as a right-censor/search cap, not as evidence that the biological distribution truly ends there. The figure marks this in three ways:

- The y-axis is capped just above `500 kb`.
- A dashed horizontal line marks `500 kb`.
- A top label states `>500 kb unobserved: analysis/search stops at 500 kb`.

Rows reported at exactly `500 kb` are shown with cap-hit markers and summarized as "at cap"; they should be interpreted as reaching the measurement/search limit, not as proven exact biological endpoint lengths.

## Generated assets

| Asset | Purpose |
| --- | --- |
| `make_06_violin_censor.R` | Reproducible generator using base R plus `ggplot2`; reads the canonical TSV directly. |
| `named_clade_violin_summary.tsv` | Per-group counts, quartiles, reported maximum, and number/percent reported at the 500 kb cap. |
| `candidate_06a_named_clade_violin_censor.png` | Talk-ready 16:9 PNG candidate for slide 06a. |
| `candidate_06a_named_clade_violin_censor.pdf` | PDF version of the same candidate. |

## Recommended deck use

Use `candidate_06a_named_clade_violin_censor.png` as the v3 06a replacement candidate if the final fan-in wants a distribution shape rather than a ranked-arm summary. It keeps the named clade story visible while making the 500 kb search cap explicit enough for a live audience.

Suggested speaker framing:

"This is a measured-length distribution, not a full biological endpoint distribution: the search stops at 500 kb, so clades piling up at 500 kb are right-censored at the measurement cap."

## Validation

- Ran `Rscript slides/v2-review-zoom/_revision_assets/v3/06_violin_censor/make_06_violin_censor.R`.
- Confirmed the candidate PNG and PDF were generated in this directory.
- Reviewed the candidate PNG for slide-scale legibility.
- Checked that no stale worktree paths are embedded in this directory, including rendered assets.
- Did not edit the Typst deck source; final fan-in owns deck edits.
