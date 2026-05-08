# Legibility Checklist

Task: `review-zoom-v8-typography-legend-cleanup`

Checked against the rule: every slide should be quickly legible in a talk, with
no dense tables or tiny unused legends.

## Checklist

| Slide | Check | Result |
| --- | --- | --- |
| 12 | Main plot title, axes, statistic text, callout, and stage trajectory labels are materially larger than v7. | Pass |
| 12 | Important figure content preserved: zygotene scatter plus four-stage trajectory. | Pass |
| 13b | Unused bottom legend is removed from the materialized crop. | Pass |
| 13b | Position axes and direct in-panel event labels remain visible in the materialized crop. | Pass |
| 14b | Long v7 row annotations removed; all eight signals have short labels and large bar-end values. | Pass |
| 14b | No legend is required; support labels are direct. | Pass |
| 14b | Statistical caveat remains visible: support counts are not q-values or BH-significant calls. | Pass |
| 14c | Dense 9 x 9 map condensed to seven communities x six signal rows. | Pass |
| 14c | Tile labels are direct and large enough to read at slide scale. | Pass |
| 14c | No color legend is required; fill intensity is secondary to direct labels. | Pass |
| 14c | Statistical-proof caveat remains visible. | Pass |

## Visual QA Notes

- Generated PNG assets are 3200 x 1800 and nonblank by visual inspection.
- Slide 14b and 14c replacements remove legends rather than shrinking legends.
- Slide 13b is delivered as a materialized deterministic crop because the
  existing pedigree image has no local regeneration recipe. Fan-in should verify
  the final rendered page after integration.
- Deck header/footer text remains controlled by the shared v7 Typst macros; this
  task only changes figure assets/snippets.
