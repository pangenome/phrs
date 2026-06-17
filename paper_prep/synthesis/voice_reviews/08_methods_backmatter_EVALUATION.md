# Evaluation: voice-review-methods

Task evaluated: `voice-review-methods`

Evaluator context: manual evaluator pass. The task is currently assigned as an
evaluation-scheduled item, but the required actor deliverable is absent from the
main worktree at evaluation time.

## Rubric Status

The rubric is sufficiently specified. The task defines a concrete manuscript
scope (`submission/paper.tex:408` to `submission/paper.tex:667`), forbids edits
to `submission/paper.tex`, names the exact required deliverable
(`paper_prep/synthesis/voice_reviews/08_methods_backmatter.md`), requires
integration of memos `01` to `07`, and lists four required deliverable sections:

1. Keep/remove/relocate table with line references.
2. Required support map from each Results/Figure claim to Methods text.
3. Drift audit for FST, CEPH1463, RPE-1, single-cell controls, exclusion
   controls, gene enrichment, and other cut analyses.
4. Integration risks that must be resolved before implementation edits.

## Grade

Overall score: **0.08 / 1.00**

Confidence: **0.93**

This is a non-delivery grade. The required file
`paper_prep/synthesis/voice_reviews/08_methods_backmatter.md` does not exist in
the worktree, and `wg show voice-review-methods` lists no actor artifacts. As a
result, the evaluator cannot credit the requested keep/remove/relocate table,
required-support map, drift audit, or integration-risk synthesis. The only
positive credit is for apparent manuscript safety: `git status --short` shows no
tracked modification to `submission/paper.tex`.

## Dimension Scores

| Dimension | Score | Rationale |
| --- | ---: | --- |
| Deliverable existence | 0.00 | Required file `paper_prep/synthesis/voice_reviews/08_methods_backmatter.md` is absent. |
| Concrete line references | 0.00 | No deliverable exists, so no line-numbered audit can be assessed. |
| Manuscript safety | 0.75 | `submission/paper.tex` has no tracked local modification at evaluation time. This supports the no-edit constraint but cannot compensate for missing delivery. |
| Integration of memos 01-07 | 0.00 | No artifact integrates the upstream voice-review memos. |
| Required Methods support map | 0.00 | Missing. |
| Reviewer-era drift distinction | 0.00 | Missing; no audit of FST, CEPH1463, RPE-1, single-cell controls, exclusion controls, gene enrichment, or other cut analyses is present. |
| Back-matter and legend consistency audit | 0.00 | Missing. |
| Implementation-risk synthesis | 0.00 | Missing. |

## Validation Checklist

| Requirement | Result |
| --- | --- |
| Deliverable exists at `paper_prep/synthesis/voice_reviews/08_methods_backmatter.md`. | Fail |
| It references concrete line numbers. | Fail |
| It does not edit `paper.tex`. | Pass, based on tracked git status at evaluation time. |
| It integrates section-review memos 01-07. | Fail |
| It distinguishes required Methods support from reviewer-era drift. | Fail |

## Notes for Follow-Up

The next actor should perform the actual Methods/back-matter review rather than
editing this evaluation file. The main manuscript span to inspect is
`submission/paper.tex:408` to `submission/paper.tex:667`, with particular
attention to drift signals already visible in the scoped Methods text:

- FST appears in the limitations sentence at `submission/paper.tex:636` to
  `submission/paper.tex:640` even though reviewer-era population-genetics
  analyses are described in the project guide as cut from the body.
- Exclusion controls remain as a full Methods subsection at
  `submission/paper.tex:556` to `submission/paper.tex:564`.
- Single-cell 3D controls remain as a full Methods subsection at
  `submission/paper.tex:566` to `submission/paper.tex:575`.
- The task explicitly asks whether these analyses should remain cut unless
  needed for internal consistency; no delivered audit answers that question.
