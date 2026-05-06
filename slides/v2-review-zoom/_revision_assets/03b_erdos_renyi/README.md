# Slide 03b Erdos-Renyi Replacement Recommendation

## Recommendation

**Yes, replace the current slide 03b ER callout.** The existing
`slides/v2-review-zoom/_typst/assets/s03_er_callout.png` is a local two-bar
callout whose rendered labels overlap on the zoom slide. A replacement should
use the `github.com:ekg/erdos_renyi` threshold logic, but adapt it to the actual
HPRC point (`n = 18,827`, `p ~ 0.12`) rather than dropping in the generic
simulation view unchanged.

Recommended candidate:

- `erdos_renyi_connectivity_candidate.png`
- `erdos_renyi_connectivity_candidate.pdf`
- generation script: `make_03b_erdos_renyi_plot.R`

## Source Inspection

External source inspected:

```bash
git ls-remote https://github.com/ekg/erdos_renyi.git HEAD refs/heads/*
git clone https://github.com/ekg/erdos_renyi.git /tmp/review-zoom-03b-erdos_renyi
git -C /tmp/review-zoom-03b-erdos_renyi log --oneline --decorate --max-count=1
```

Observed source state:

- Remote: `https://github.com/ekg/erdos_renyi.git`
- Commit: `d9ec48f1945d14f38d0131f56ed288cfc7883e73`
- Files: `erdos_renyi_viz.R`, `viz.html`

Best source element for slide 03b:

- Use `erdos_renyi_viz.R`, specifically the threshold definitions
  `giant_threshold <- 1/n` and `connected_threshold <- log(n)/n`.
- Do not use `viz.html` for the review deck: it is an interactive D3 circle
  graph demo, useful for intuition but not a static slide-ready explanation of
  the HPRC scale point.
- Do not use the repo's generic simulation-cloud plot unchanged: it samples
  `n = 10..500`, while slide 03b needs to justify `n = 18,827` flanks and
  `p ~ 0.12`.

## Regeneration Recipe

From the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/03b_erdos_renyi/make_03b_erdos_renyi_plot.R \
  slides/v2-review-zoom/_revision_assets/03b_erdos_renyi
```

This writes:

- `slides/v2-review-zoom/_revision_assets/03b_erdos_renyi/erdos_renyi_connectivity_candidate.png`
- `slides/v2-review-zoom/_revision_assets/03b_erdos_renyi/erdos_renyi_connectivity_candidate.pdf`

The script is self-contained except for `ggplot2`, which is already available in
this worktree's R installation.

## What The Plot Means

- The x-axis is graph size `n`; for slide 03b, `n = 18,827` subtelomeric
  flanks.
- The y-axis is edge probability `p` on a log scale. The black curve is the
  Erdos-Renyi connectivity threshold, `p* = log(n) / n`; at `n = 18,827`,
  `p* = 5.23e-4`.
- The blue point marks the deck's method claim: wfmash evaluates approximately
  `12%` of pair space, or `p ~ 0.12`.
- The vertical gap between `5.23e-4` and `0.12` is the important slide beat:
  `0.12 / p* = 229.5`, so the sampling rate is about `230x` above the random
  graph connectivity threshold.
- This belongs on slide 03b because it visually supports the methods argument
  that the implicit pangenome graph is densely connected enough for
  transitive-closure queries without chromosome partitioning.

## Integration Notes

- Use the candidate as a full-slide replacement or as the dominant panel on
  03b; it is 16:9 and should fit the existing `figure-slide` layout without
  label collisions.
- Keep the caveat in a caption or speaker note: wfmash pair selection is
  k-mer-selective, not literally random. The ER curve is a connectivity sanity
  check, not a biological null model.
- If the fan-in task wants a smaller callout rather than a full-slide plot, use
  the same script as the exact recipe and reduce label density after confirming
  the rendered slide remains legible.
