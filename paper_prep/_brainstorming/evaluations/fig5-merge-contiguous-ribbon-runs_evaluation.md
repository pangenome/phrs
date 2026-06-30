# Evaluation: fig5-merge-contiguous-ribbon-runs

Task: `fig5-merge-contiguous-ribbon-runs`  
Evaluator: `agent-2960`  
Evaluation time: 2026-06-30T19:10Z-19:18Z UTC  
Overall grade: 0.44 / 1.00  
Confidence: 0.86  
Rubric underspecified: yes - no explicit scoring rubric or `## Validation` checklist was provided, so I graded against the concrete task description.

## Basis

The task asked the actor to optimize the Fig5 whole-genome ribbon figures by
merging perfectly contiguous or overlapping run ranges between the same query
interval source and donor interval target before rendering. The task explicitly
called out both homologous donor-to-child chains and colored interchromosomal
runs where appropriate, required use of the existing SweepGA/F32 10:10 IMPG
class-winner source and plotting package only, and requested regenerated
PNG/PDF/SVG plus summary/docs updates if run counts changed.

I evaluated commit `4e01d35` (`feat: merge fig5 child-centered homolog ribbons
(agent-2958)`) and the current contents of
`paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft`.

## Evidence Checked

- `git show --stat 4e01d35` shows regenerated SVG/PDF/PNG outputs, updated
  `README.md`, updated `whole_genome_homologous_context_runs.tsv`, updated
  `whole_genome_homologous_context_summary.tsv`, and edits to
  `scripts/plot_whole_genome_ribbon_draft.py`.
- The plotting script still reads the existing 10:10 class-winner source:
  `PAN027pat_vs_PAN011_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz`,
  with a local path fallback and `/moosefs/erikg/phrs/...` source path. I found
  no evidence of a new alignment.
- The homologous-context rendered outputs exist:
  `fig5_homologous_recombination_context_ribbon_draft.svg`,
  `fig5_homologous_recombination_context_ribbon_draft.pdf`, and
  `fig5_homologous_recombination_context_ribbon_draft.png`.
- The base whole-genome ribbon outputs also exist and were regenerated:
  `fig5_homolog_vs_interchrom_whole_genome_ribbon_draft.svg`,
  `.pdf`, and `.png`.
- The homologous summary changed from `drawn_homolog_runs = 8183` and
  `homolog_min_bp = 10000` in the parent to `drawn_homolog_runs = 364` and
  `homolog_min_bp = 50000` after the commit.
- The colored interchromosomal summary is unchanged before and after the commit:
  `all_inter_beats_same_runs = 918`, `drawn_high_conf_runs = 59`,
  `drawn_high_conf_bp = 1118000`, and all category counts remain identical.
- The merge predicates in `group_runs()` and `group_homolog_runs()` were already
  present in the parent commit and were not changed by `4e01d35`.
- Those existing predicates are not the requested "perfectly contiguous or
  overlapping" semantics. They merge across tolerated gaps:
  `segment.query_start <= current.query_end + 2_000` and donor gaps up to
  `10_000` for interchromosomal runs, and `+10_000` query / `50_000` donor
  tolerance for homologous runs.
- The observed homologous run-count reduction appears to come primarily from
  increasing the display threshold from 10 kb to 50 kb, not from implementing a
  stricter contiguous/overlap merge.

## Dimension Scores

| Dimension | Score | Rationale |
|---|---:|---|
| Correct contiguous/overlap merge semantics | 0.15 | The task's core requirement was to merge perfectly contiguous/overlapping ranges with the same source/target. The actor did not change the merge functions and the retained predicates merge across sizable gaps, so the specified semantics are mostly unmet. Some pre-existing grouping by query/donor sequence remains relevant. |
| Homologous donor-to-child chain handling | 0.45 | The homologous-context figure became much less cluttered and now uses native interval endpoints, with drawn homologous runs dropping from 8183 to 364. However, this was achieved mainly by raising the minimum displayed run length to 50 kb rather than by implementing the requested merge behavior. |
| Colored interchromosomal run handling | 0.20 | Colored run counts and bp are unchanged before/after. Existing interchromosomal grouping remains, but the actor did not apply a new contiguous/overlap merge or document that no appropriate additional merges existed. |
| Source/provenance constraint | 0.95 | The script uses the existing 10:10 IMPG class-winner TSV and there is no evidence of a new alignment. The added local/moosefs fallback improves reproducibility. |
| Regenerated visual deliverables | 0.95 | SVG/PDF/PNG artifacts for both the base and homologous-context figures were regenerated and committed. |
| Summary/docs updates | 0.75 | The README and homologous summary were updated and reflect the changed displayed run counts and 50 kb threshold. The wording still overstates "full same-chromosome chains" relative to the tolerant grouping/filtering actually implemented. |
| Validation evidence | 0.35 | The committed outputs and summaries provide some self-checking evidence, but there is no explicit validation of contiguous/overlap merge correctness, no before/after audit of merge candidates, and no test or script proving that only same query-source/donor-target contiguous ranges are coalesced. |
| WG/process completion | 0.45 | The implementation commit is present on `main`, but the WG task record available to this evaluator has sparse actor logs/artifacts and the task itself was still shown as `in-progress` in this worktree. |
| Non-destructive behavior | 0.90 | The actor kept changes scoped to the Fig5 ribbon draft directory and did not damage unrelated project files. |

## Overall Grade

`0.44 / 1.00`

This is a useful but incomplete result. The actor used the right source data,
regenerated the required raster/vector outputs, and updated summary/docs. The
homologous-context figure is materially smaller and probably more readable.

The main problem is that the requested algorithmic change was not actually
implemented. The merge functions were already present in the parent commit, and
their predicates still allow non-contiguous gaps, contradicting the requirement
to merge only perfectly contiguous or overlapping ranges. The colored
interchromosomal run counts did not change and there is no audit showing that no
additional exact merges were appropriate. The homologous count reduction is
mostly a display-threshold change, not the requested contiguous-run merge.

## Calibration Notes

- A result that implemented exact `next.start <= current.end` coalescing with
  explicit source/target keys for both homologous and interchromosomal runs,
  regenerated outputs, and documented before/after counts would score around
  0.85.
- A strong result with a small validation fixture or audit proving no gap-based
  merges remain would score 0.90 or higher.
- This result earns substantial credit for regenerated deliverables and source
  discipline, but loses most of the core algorithm score.
