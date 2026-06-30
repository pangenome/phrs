# Evaluation: fig5-whole-genome-10to10-ribbon-overlay-draft

Task: `fig5-whole-genome-10to10-ribbon-overlay-draft`  
Evaluator: `agent-2943`  
Evaluation time: 2026-06-30T14:33Z-14:38Z UTC  
Overall grade: 0.02 / 1.00  
Confidence: 0.94  
Rubric underspecified: yes - no explicit scoring rubric or `## Validation` checklist was provided, but the task description gave concrete deliverables and constraints.

## Basis

The task asked for a draft whole-genome ribbon overlay from the existing
`PAN027pat_vs_PAN011_joint` SweepGA/F32 10:10 filtered IMPG class-winner data,
with no new alignments. Required outputs were:

- a deliverable directory at `paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft`;
- length-scaled genome tracks for PAN027 paternal child query and PAN011 father donor haplotypes;
- ribbons for high-confidence interchrom-over-same runs;
- SVG, PDF, PNG exports;
- a run table;
- use of the same class-winner source as the Fig5 zoom/ribbon and whole-genome homolog-vs-interchrom overview.

I evaluated the current task state, repository tree, and nearby Fig5 artifacts in
`/moosefs/erikg/phrs/.wg-worktrees/agent-2943`.

## Evidence Checked

- `wg show fig5-whole-genome-10to10-ribbon-overlay-draft` showed the task still `in-progress`, with 0 commits ahead and no uncommitted files in this evaluator worktree.
- `wg show .assign-fig5-whole-genome-10to10-ribbon-overlay-draft` showed only an assignment wrapper completed by `agent-2942`; it did not list artifacts or implementation logs.
- `paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft` does not exist.
- `find` found no whole-genome ribbon overlay outputs or run table matching the requested deliverable. The only ribbon table found was the pre-existing dependency `paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/donor_recipient_runs.tsv`.
- Existing dependency directories are present:
  `fig5_donor_recipient_ribbon_draft`, `fig5_homolog_vs_interchrom_whole_genome`, and `fig5_whole_genome_length_scaled_tracks`.
  These are inputs or related prior products, not the requested combined whole-genome 10:10 ribbon overlay deliverable.
- `rg` found no repository references to `fig5_homolog_vs_interchrom_whole_genome_ribbon_draft`, `interchrom-over-same`, or a new 10:10 whole-genome ribbon overlay implementation.
- Git history search around the assignment window found no commit for this task before the evaluation artifact.

## Dimension Scores

| Dimension | Score | Rationale |
|---|---:|---|
| Correct source data and no-new-alignment constraint | 0.05 | No evidence of new alignments, but also no evidence that the required 10:10 class-winner source was consumed for this task. |
| Whole-genome length-scaled PAN027pat/PAN011 tracks | 0.05 | A dependency directory contains length-scaled tracks, including a PAN027 hap2 paternal vs PAN011 product, but no new combined overlay deliverable was created. |
| High-confidence interchrom-over-same ribbon extraction | 0.00 | No run table or implementation identifying high-confidence interchrom-over-same runs was produced for the requested whole-genome overlay. |
| Visualization deliverables | 0.00 | The requested deliverable directory and required SVG/PDF/PNG exports are absent. Existing SVG/PDF/PNG files belong to prior dependency artifacts. |
| Run table deliverable | 0.00 | No new run table exists for the whole-genome 10:10 ribbon overlay; only the older donor-recipient run table exists. |
| Reproducibility/provenance | 0.00 | No script, README, manifest, validation note, or artifact registration was added for this task. |
| WG/process completion | 0.00 | No actor commit or recorded artifacts were found, and the implementation task remained `in-progress` at evaluation time. |
| Non-destructive behavior | 0.25 | There is no evidence the actor damaged existing artifacts or ran prohibited new alignments, but this is absence of harm rather than task progress. |

## Overall Grade

`0.02 / 1.00`

The actor did not produce the requested deliverable. The strongest available evidence is that prerequisite artifacts already existed, but the task specifically required a new whole-genome ribbon overlay draft combining the 10:10 class-winner source, length-scaled PAN027 paternal/PAN011 tracks, high-confidence interchrom-over-same ribbons, SVG/PDF/PNG exports, and a run table. None of those requested outputs exist in the named deliverable location or elsewhere under an identifiable task-specific name.

I assign a small nonzero score only because the repository already contains related dependency artifacts that could have supported the task, and there is no evidence of destructive behavior or prohibited alignment reruns. The performance is otherwise a near-total non-completion.

## Calibration Notes

- A minimally acceptable draft with a correct directory, generated SVG/PDF/PNG, a run table, and clear provenance from the existing 10:10 class-winner TSV would score around 0.70 even if the visual polish were rough.
- A strong result with validated run filtering, clear scripts/README, and visibly useful whole-genome ribbons would score 0.90 or higher.
- This task lacked a formal rubric, so `Rubric underspecified` is flagged. The low score is not due to ambiguity in the deliverables; it is due to the absence of the requested outputs.
