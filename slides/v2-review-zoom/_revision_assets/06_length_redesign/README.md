# Review zoom 06 length redesign candidates

Task: `review-zoom-06-length-redesign`

## Inputs Inspected

- `slides/v2/slide_06_length_distributions.md`
- `slides/v2/figure_manifest.md`
- Current review-zoom pages/assets:
  - `slides/v2-review-zoom/_typst/page-07.png` (`06a`, four-crop full faceted split view)
  - `slides/v2-review-zoom/_typst/page-08.png` (`06b`, clade callout table)
  - `slides/v2-review-zoom/_typst/assets/s06_length_dist*.png`
  - `slides/v2-review-zoom/_typst/assets/s06_clade_callouts.png`
- Source TSV: `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv`
- Cross-slide clade vocabulary from `paper_prep/synthesis/CROSSWALK.md`, `end-to-end-report/report/01_pipeline.md`, and the slide-09 review notes.

## Data Check

The generator reads the TSV directly and uses `region_end - region_start` for non-empty rows.

- Non-empty PHR intervals: `15,668`
- Arms with signal: `41/48`
- Global length scale: median `105 kb`, mean `143.8 kb`, range `5-500 kb`
- Current TSV zero-signal arms: `2p, 3p, 5p, 8q, 11q, 14q, 18q`

The zero-arm list needs final integration attention. The existing slide-06 and slide-09 notes repeatedly mention six introvert arms (`2p, 3p, 5p, 8q, 11q, 14q`) and slide-09 notes also mention `18q` as a C15 outlier. The current TSV used here has `18q` at `n=0`, and the current 06a rendered page also shows `chr18 qarm` as a pink `n=0` panel. I preserved the data-derived seven-arm list in the candidate assets rather than silently copying the older note.

## Generated Assets

All outputs are isolated in this directory. Final deck integration can pick any subset without changing canonical slide source here.

| Asset | Purpose |
| --- | --- |
| `make_06_length_redesign.R` | Reproducible generator using base R plus `ggplot2`. |
| `arm_length_summary.tsv` | Per-arm `n`, median, mean, p25, p75, p90, p95, max, and clade class. |
| `clade_length_summary.tsv` | Group-level summaries for C7, C14, C15, C1, C2, no-signal arms, and all other signaled arms. |
| `candidate_06a_ranked_arm_summary.png` / `.pdf` | Talk-friendly ranked arm summary. |
| `candidate_06a_focused_clade_facets.png` / `.pdf` | Focused histogram replacement for the dense 48-facet 06a. |
| `candidate_06b_clade_story_matrix.png` / `.pdf` | Stronger 06b callout matrix tying length evidence to the biological/community story. |

## Recommendation

Best two-slide talk sequence:

1. Use `candidate_06a_ranked_arm_summary.png` as the primary 06a replacement.
2. Use `candidate_06b_clade_story_matrix.png` as the 06b replacement.
3. Move the current full faceted split view to backup/appendix.

This gives the audience one slide that explains the arm-level length hierarchy and one slide that names the biological mechanism behind each highlighted row. It also avoids asking people to read 48 small histograms during a live talk.

If slide 06 must be only one spoken slide, use `candidate_06a_ranked_arm_summary.png` and fold the C7/C14/C15/C1/C2/introvert wording into speaker notes. Keep `candidate_06b_clade_story_matrix.png` as the appendix/Q&A map.

If the team wants to preserve histogram shapes in the spoken deck, use `candidate_06a_focused_clade_facets.png` instead of the ranked summary. It retains distributions but collapses the panel set to bulk background plus the named clades and zero-signal arms.

## Candidate Framing

| Candidate | Talk vs. appendix | One-sentence speaker framing |
| --- | --- | --- |
| `candidate_06a_ranked_arm_summary` | Primary talk recommendation for 06a. | "Ranked by arm, the length scale is not random: full-window C7/PAR1, PAR2 at the canonical 330 kb scale, a DUX4/D4Z4 tail, a short specific 10p-18p pair, and seven arms with no inter-chromosomal PHR signal." |
| `candidate_06a_focused_clade_facets` | Alternate talk 06a if histogram shape must remain visible; otherwise appendix. | "Instead of 48 tiny facets, this view compares the bulk background with the five named clade groups and the zero-signal arms, so the distribution shapes set up the community partition." |
| `candidate_06b_clade_story_matrix` | Primary talk recommendation for 06b if two slide-06 pages survive; appendix if only one slide survives. | "The long tails and blank panels are named biology: acrocentric rDNA-adjacent homogenization, PAR1/PAR2, DUX4/D4Z4 copy-number diversity, the Linardopoulou 10p-18p pair, and introvert arms absent from the arm matrix." |
| Current dense 06a split view | Appendix/reference only. | "This is the complete audit view: every arm-level histogram is preserved for per-arm checks, but it is too dense to carry the live story." |

## Strengthened 06b Story

Use the following language consistently with slide 09 and downstream community slides:

| Clade/group | Data-derived length signal | Biological/community point |
| --- | --- | --- |
| C7 acrocentric p (`13p, 14p, 15p, 21p, 22p`) | `n=763`, median/p90 `500/500 kb` | rDNA-adjacent acrocentric p-arms are the clearest full-window homogenization signal. |
| C14 PAR2 (`Xq, Yq`) | `n=431`, median/p90 `330/330 kb` | PAR2 is the familiar pseudoautosomal q-end scale and the anchor for the abstract's "comparable to PAR2" wording. |
| C15 PAR1 (`Xp, Yp`) | `n=419`, median/p90 `500/500 kb` | PAR1 makes the sex-chromosome p-ends a full-window sharing clade, visible before clustering. |
| C1 DUX4/D4Z4 (`4q, 10q`) | `n=714`, median/p90 `145/215 kb`, max `500 kb` | Moderate center with a real 500 kb tail, matching DUX4/D4Z4 copy-number diversity. |
| C2 10p-18p (`10p, 18p`) | `n=889`, median/p90 `70/105 kb` | A named community and historical exchange pair, not a long-tail outlier; do not oversell it as length-heavy. |
| No interchrom PHR (`2p, 3p, 5p, 8q, 11q, 14q, 18q`) | `n=0` in current TSV | The blank/pink panels are biological absence from the inter-chromosomal sharing landscape, not missing render data. |

## Appendix-Style Dense Spec

If final integration keeps a dense appendix slide, use the current full faceted split view as the base, but revise the callout/notes as follows:

- Treat the dense facet grid as a per-arm audit, not the spoken explanation.
- Keep p-arm blue, q-arm orange, and no-signal arms pink for continuity with v1.
- Add a compact overlay or caption that maps C7, C14, C15, C1, and C2 to their arms.
- Reconcile the zero-signal arm list against the exact TSV/render being used; current data says seven arms including `18q`.
- Do not describe C2 as a fat-right-tail group. It is a named specific exchange pair with short/tight lengths in this TSV.

## Validation

- Ran `Rscript slides/v2-review-zoom/_revision_assets/06_length_redesign/make_06_length_redesign.R`.
- Confirmed the script generated PNG and PDF variants for all three candidate figures.
- Confirmed `arm_length_summary.tsv` and `clade_length_summary.tsv` were generated from the source TSV.
- Reviewed rendered PNGs for legibility.
- Did not edit canonical deck source or current review-zoom Typst files.
