# Evaluation: fig5-flex-merge-homologous-inheritance-blocks

Task: `fig5-flex-merge-homologous-inheritance-blocks`  
Evaluator: `agent-2974`  
Evaluation time: 2026-06-30T21:47Z-22:02Z UTC  
Overall grade: 0.71 / 1.00  
Confidence: 0.83  
Rubric underspecified: yes - no explicit scoring rubric or `## Validation` checklist was provided, so I graded against the concrete task description.

## Basis

The task asked the actor to investigate and implement a more flexible display
merge for Fig5 homologous father-child chains. The motivating defect was that
the exact end-to-end 2 kb merge fragments long inherited paternal-haplotype
blocks when windows are missing or slightly discontinuous. The requested work
included quantifying query/donor gap distributions, choosing a conservative
homolog-only merge tolerance, preserving haplotype switches and
interchromosomal winner overlays, regenerating PNG/PDF/SVG outputs, auditing
outputs, and committing/pushing.

I evaluated commit `9e4cd42` (`fix: merge fig5 homolog blocks for display
(agent-2974)`) and the contents it writes in
`paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft`.

## Evidence Checked

- The exact 2 kb end-to-end merge remains explicitly zero-gap:
  `CONTIGUOUS_MERGE_GAP_BP = 0` in
  `scripts/plot_whole_genome_ribbon_draft.py:79`, preserving the prior exact
  run layer.
- The actor added display-only homolog merge tolerances:
  `HOMOLOG_DISPLAY_MERGE_MAX_QUERY_GAP_BP = 5_000_000` and
  `HOMOLOG_DISPLAY_MERGE_MAX_DONOR_GAP_BP = 10_000_000` at
  `scripts/plot_whole_genome_ribbon_draft.py:80-81`.
- `merge_homolog_display_blocks()` buckets runs by
  `(query_seq, query_chrom, donor_haplotype, target_chrom, donor_seq)` and only
  coalesces within those buckets when query and donor gaps are within tolerance
  (`scripts/plot_whole_genome_ribbon_draft.py:602-632`). This is a genuine
  homolog-only display merge path layered on top of exact end-to-end runs.
- Interchromosomal grouping still uses the original exact
  `group_end_to_end_runs()` path and is not passed through the homolog display
  merger (`scripts/plot_whole_genome_ribbon_draft.py:527-532` and
  `main()` diff around `scripts/plot_whole_genome_ribbon_draft.py:1241-1248`).
- Span validation was adjusted so display blocks may have query spans larger
  than their covered bp while still rejecting first-window-only or shorter-than-
  bp spans (`scripts/plot_whole_genome_ribbon_draft.py:543-559`).
- The homologous summary records the selected display tolerances and resulting
  output: `homolog_display_merge_max_query_gap_bp = 5000000`,
  `homolog_display_merge_max_donor_gap_bp = 10000000`,
  `homolog_end_to_end_runs = 66592`, `drawn_homolog_runs = 94`,
  `drawn_homolog_bp = 347828487`, and `drawn_max_query_span_bp = 244276000`.
- The merge audit marks the display tolerances as `NA` for the interchromosomal
  layer and records the nonzero homolog display tolerances only for the homolog
  layer, while keeping interchromosomal `drawn_runs = 52` and
  `drawn_bp = 940000`.
- The figure script text now labels the light-gray layer as
  "display-merged same-chromosome father-child homologous blocks" and notes
  the `<=5 Mb` query-gap display merge in the rendered footnote
  (`scripts/plot_whole_genome_ribbon_draft.py:1098-1103` and `1170-1175` in
  the evaluated commit).
- Commit `9e4cd42` regenerated the homologous-context PDF/PNG/SVG, base PDF,
  homolog runs TSV, homolog summary TSV, merge audit TSV, and ribbon summary
  TSV. `conversion_status.txt` reports successful conversion with
  `rsvg-convert version 2.54.5`.
- The README in the evaluated tree was not updated for the new display-merge
  behavior and still describes the older exact/no-gap grouping and 50 kb drawn
  threshold. That creates stale package-level documentation despite correct
  script/output summaries.

## Dimension Scores

| Dimension | Score | Rationale |
|---|---:|---|
| Gap distribution quantification | 0.30 | The actor records exact-run counts, display-block counts, max spans, and tolerance values, but there is no actual query/donor gap distribution: no candidate gap table, quantiles, histogram, or before/after absorbed-gap summary explaining how 5 Mb / 10 Mb was selected. |
| Conservative tolerance choice | 0.45 | A concrete homolog display tolerance was chosen and exposed in summary/audit outputs. However, 5 Mb query and 10 Mb donor are large tolerances, and the commit does not provide evidence that these are conservative relative to observed gap distributions or switch boundaries. |
| Homolog-only flexible merge implementation | 0.82 | The new `merge_homolog_display_blocks()` path is scoped to homologous runs, uses same query/donor/haplotype keys, bridges nonzero query/donor gaps, and feeds the gray display layer. This substantially satisfies the core implementation request. |
| Preservation of switches and overlays | 0.78 | The merge key includes donor haplotype and donor sequence, which should preserve paternal haplotype switches, and the interchromosomal layer remains separately audited with 52 drawn calls. The gap tolerance is not explicitly validated against switch examples, so this is good but not fully proven. |
| Regenerated visual deliverables | 0.92 | The commit includes regenerated homolog-context PNG/PDF/SVG and related figure outputs. Conversion status reports successful PDF/PNG rendering. |
| Audit and documentation | 0.62 | Summary and audit TSVs expose the new tolerances, layer-specific counts, and span metrics. The README remains stale and there is no tolerance-selection note, so package documentation is only partially complete. |
| Validation evidence | 0.66 | The script validates exact runs and permits longer display spans only for display blocks, and the audit records that interchrom outputs remain unchanged. Missing pieces are a deterministic fixture, explicit switch-preservation checks, and direct gap-distribution validation. |
| Source/provenance discipline | 0.95 | The work stays within the Fig5 ribbon package and uses the established 10:10 IMPG class-winner source. I saw no unrelated source/data churn. |
| WG/process completion | 0.85 | The actor logged a concrete validation summary, committed `9e4cd42`, and pushed it to `main`. Task-specific artifacts are present via generated outputs, though the evaluator branch had to fetch the latest main to see the final commit. |
| Non-destructive behavior | 0.95 | Changes are scoped to the Fig5 ribbon draft package and generated outputs. No unrelated files were modified. |

## Overall Grade

`0.71 / 1.00`

This is a materially useful implementation of the requested idea. The actor
added a separate homolog-only display-block merge on top of exact 2 kb
end-to-end chains, kept interchromosomal winners out of that flexible merge,
regenerated the visual outputs, and exposed the new tolerances and resulting
counts in the summary/audit TSVs. The output reduction from 66,592 exact
homolog chains to 94 drawn homolog display blocks directly addresses the
fragmented gray inheritance-chain display.

The main reasons this is not a higher score are evidentiary and calibration
gaps. The task explicitly asked to quantify query/donor gap distributions and
choose a conservative tolerance; the commit chooses 5 Mb query and 10 Mb donor
tolerances but does not show the distribution or rationale behind those values.
The README also remains stale and still describes the previous exact/no-gap
behavior. Finally, while the merge keys are sensible for preserving haplotype
switches and interchromosomal overlays, there is no explicit validation against
switch cases or interchrom overlay regressions beyond unchanged layer counts.

## Calibration Notes

- A result with only display-width scaling and no flexible merge would score
  around 0.30-0.40.
- This result earns substantial credit because it implements a real homolog-only
  gap-tolerant display merge and regenerates/audits the outputs.
- A stronger result with gap quantiles, before/after absorbed-gap distributions,
  a written tolerance rationale, README updates, and explicit switch/overlay
  checks would score around 0.85-0.90.
- A near-perfect result would add a small deterministic regression fixture for
  discontinuous homolog windows, haplotype switches, and interchromosomal
  winners under the chosen tolerance.
