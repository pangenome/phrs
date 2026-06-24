# Evaluation: fig5-query-grid-panels-regenerate

Date: 2026-06-24
Evaluator: agent-2707
Rubric status: sufficiently specified by task acceptance criteria.
Underspecification flag: no

## Grade

Overall score: 0.00 / 1.00
Confidence: 0.99

## Dimension Scores

| Dimension | Score | Rationale |
|---|---:|---|
| New query-grid figure package | 0.00 | The required package `paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_query_grid_chop_filter_panels/` is absent, and no equivalently explicit query-grid panel directory was found under `paper_prep/_brainstorming/`. |
| Rendered PDF, PNG, and SVG outputs | 0.00 | Repository search found only older Fig5 raw FASTA panel PDFs. No query-grid Fig5 panel PDF, PNG, or SVG outputs were present. |
| Segment and summary TSVs with genomic query coordinates | 0.00 | No task-specific segment-level TSV or summary TSV exists in a query-grid panel package. The dependency rerun summaries exist, but they are not the requested panel package summaries and do not satisfy the downstream figure deliverable. |
| Manifest with sources, commands, and checksums | 0.00 | The dependency rerun manifest `query_grid_chop_filter_manifest.tsv` exists, but no figure-package manifest records the raw many:many source PAFs, query-grid chopped PAFs, filtered PAFs, SweepGA settings, pafchop settings, and sha256 checksums for the panel deliverable. |
| Plot content and labeling requirements | 0.00 | There is no new query-grid plot to inspect for PAR1, PAN027/PAN028 chr9q-to-chr3q windows, 10 kb/5 kb/2 kb/1 kb chop lengths, query-grid labels, SweepGA flag labels, separation of raw many:many vs post-filter support, or legend overprint fixes. |
| Validation evidence | 0.00 | No validation script output, rendered panel validation, or task-specific WG artifact was found. |
| README explanation | 0.00 | No README for a query-grid panel package exists, so there is no explanation of how query-grid chopping differs from older row-start chopped panels. |
| Commit and provenance | 0.00 | The task worktree has zero commits ahead of `main`; the most recent commit is the upstream dependency rerun `feat: fig5-f16-query-grid-chop-filter-rerun (agent-2700)`, not the required `feat: fig5-query-grid-panels-regenerate (agent-NNN)`. |

## Evidence Checked

- `wg show fig5-query-grid-panels-regenerate` reported the task assigned to this evaluator worktree with zero commits ahead and no registered deliverable artifacts.
- `git log --oneline --decorate -10` showed no task-specific implementation commit; `HEAD` is `97e3e96 feat: fig5-f16-query-grid-chop-filter-rerun (agent-2700)`.
- `find paper_prep/_brainstorming -maxdepth 2 -type d -name '*query*grid*panel*' -o -name '*fig5*query*grid*'` found no matching package directory.
- `find paper_prep/_brainstorming -maxdepth 2 -type d -iname '*fig5*raw*fasta*sweepga*f16*'` listed only prior non-query-grid packages: chop-filter sensitivity, chopped panels, multiway panels, no-chop merge panels, and scaffold-jump sensitivity.
- Targeted searches for query-grid Fig5 rendered PDF/PNG/SVG outputs, query-grid segment TSVs, query-grid summary TSVs, query-grid README files, and query-grid panel manifests found no task-specific deliverables.
- Dependency artifacts from `fig5-f16-query-grid-chop-filter-rerun` are present under `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/`, but they are upstream chopped/filter inputs rather than the requested regenerated Fig5 panel package.

## Rationale

The task rubric is concrete: create a new query-grid Fig5 panel package, render PDF/PNG/SVG panels, produce segment and summary TSVs with genomic query coordinates, include a complete manifest with command settings and checksums, document the query-grid-vs-row-start distinction, validate rendering, and commit with the specified provenance convention. None of the task-specific deliverables are present in the repository state available to this evaluator.

This is therefore a non-submission for the requested downstream figure-generation task. The upstream dependency rerun appears to have completed successfully, but the actor did not convert those artifacts into the required Fig5 query-grid panel package. Because all acceptance criteria depend on persisted outputs that are absent, the calibrated score is 0.00.

## Residual Risk

Very low. It is possible that an actor generated files outside the repository or failed before committing, but the task explicitly required committed outputs and a provenance-preserving commit. Uncommitted or external work would not satisfy the acceptance criteria.
