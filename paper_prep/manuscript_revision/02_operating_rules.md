# Manuscript Revision 02 Operating Rules

Date: 2026-06-17

Scope: coordination rules for downstream `manuscript-revision-*` workers. This
artifact does not edit the manuscript, captions, figures, bibliography, or
analysis outputs. It defines how workers should use WG fanout and Slurm while
preserving the edit authority locked in
`paper_prep/manuscript_revision/01_fanout_graph.md`.

## Required Authority Model

All downstream workers must follow the fanout graph authority model:

- Mechanical audit tasks may inspect files, collect evidence, and propose exact
  patches. They must not apply manuscript edits unless their task description
  explicitly grants edit authority. Otherwise, integration happens in a later
  fan-in task.
- Mechanical fan-in tasks may apply only clearly supported mechanical fixes
  when their task description explicitly authorizes editing.
- Judgment tasks must not silently edit manuscript claims. They must write
  decision records or decision packages that identify the evidence, the
  recommended claim boundary, and author-decision points.
- Compute tasks may create scripts, Slurm files, logs, result tables, plots, and
  reports under `paper_prep/manuscript_revision/`, but scientific wording still
  flows through decision records, fan-ins, and the guarded manuscript patch
  task.
- `submission/paper.tex`, `submission/bibliography.bib`, figure captions, and
  figure assets are not open editing surfaces unless the current task explicitly
  says it can edit them.

## WG Fanout Rules

Workers are authorized to create their own WG subtasks when the assigned task
naturally decomposes into separable stages such as:

- data discovery or source inventory;
- script writing;
- Slurm submission setup;
- Slurm job polling;
- result validation;
- synthesis into the requested report, audit, or decision record.

For larger data processing, workers should prefer visible workgraph subtasks
over monolithic hidden work inside one agent turn. This is especially important
for C-0/C-D continuum tasks, F-3 mouse analysis, and any task that may need
cluster compute, repeated polling, or independent validation.

Subtasks must be structured so the graph remains recoverable:

- Use `wg add "<title>" --after <predecessor>` for every child task. Do not
  create unordered flat subtasks.
- Give each subtask a clear scope, expected artifact, and `## Validation`
  checklist.
- Put dependencies in WG edges, not only prose. A script-writing task should
  precede a Slurm-submit task; a submit task should precede polling; polling
  should precede validation; validation should precede synthesis.
- If multiple subtasks fan out, create a downstream integrator or synthesis task
  with `--after` edges from all parts.
- Do not create parallel subtasks that write the same file. If two steps need
  the same file, make them sequential.
- Register every durable file with `wg artifact <task-id> <path>`.

Subtasks are not a way to bypass manuscript-edit restrictions. A child task
inherits the parent task's edit authority unless the parent task explicitly
creates a narrower coordination/compute subtask.

## Head-Node Safety Rules

Workers may submit Slurm jobs as needed from compute/tux nodes, but must keep
heavy processing off the head node.

Do not run these operations on the head node:

- heavy matrix scans;
- full 15,668 x 15,668 sequence-level matrix traversal;
- large RDS loads or conversions;
- multi-core permutations, bootstrap scans, Leiden sweeps, or Mantel loops;
- long-running memory-intensive R or Python jobs;
- any exploratory command likely to monopolize CPU, memory, or I/O.

Head-node work should be limited to lightweight inspection, file-size checks,
small previews, manifest writing, script editing, `sbatch` submission, `squeue`
or `sacct` polling, and markdown synthesis. When in doubt, write a Slurm script
and submit it instead of running the workload interactively.

Small analyses may stay inline only when the inputs are clearly tiny, such as
the 41 x 41 arm-level matrix described for `manuscript-revision-c0a`.

## Slurm Logging Requirements

Every worker that submits or supervises Slurm work must leave enough information
for another agent to reproduce, poll, validate, or recover the run. Record the
following in the WG log and in the task artifact report:

- Slurm job ID.
- Submission command, including the exact `sbatch` command or wrapper command.
- Script path and git-relative path when the script is committed.
- Working directory.
- Resource request: partition/queue, nodes, tasks, CPUs, memory, time limit, and
  any array range.
- Environment or module setup, including Guix manifest path if used.
- Input paths, especially any `/moosefs/...` source-of-truth data.
- Stdout path.
- Stderr path.
- Expected output artifact paths.
- Completion status from `squeue`, `sacct`, or equivalent polling.
- Validation command or check used after completion.

Prefer deterministic output paths under
`paper_prep/manuscript_revision/<task-or-topic>/` for scripts, logs, and result
tables. If a Slurm script writes temporary files elsewhere, the final task
report must point to both the temporary location and the durable artifact copied
or summarized in the repo.

## Decision Records For Judgment Tasks

Judgment tasks must produce explicit decision records instead of silent edits.
A decision record should include:

- the manuscript claim or figure/caption claim under review;
- the evidence inspected, with local paths or upstream `/moosefs/...` paths;
- what is established by the evidence;
- what is not established or remains overclaimed;
- recommended wording direction or claim boundary;
- author-decision points, clearly labeled;
- downstream fan-in or patch task that should consume the decision.

Judgment records may include exact suggested language, but that language is a
proposal. It is not an applied manuscript change unless the task explicitly
grants editing authority or a later fan-in/guarded patch task applies it.

## Mechanical Patch Proposals

Mechanical tasks may propose exact patches when the correction is supported by
the audit, for example:

- sample-count reconciliation;
- duplicate bibliography cleanup;
- typo, acronym, or symbol-collision fixes;
- replacement of verified DOI metadata;
- removal of construction artifacts.

Unless the task explicitly says it can edit, the worker should write the patch
proposal into the requested artifact and leave integration to the relevant
fan-in task. Patch proposals should identify the target file, current text or
line context, replacement text, and validation needed after integration.

## Minimum Artifact Checklist

Every downstream worker should leave an artifact that answers:

- What was inspected or run?
- Which files, tables, scripts, or external data paths were used?
- What WG subtasks, if any, were created and why?
- For Slurm work, which jobs were submitted and where are the logs?
- What outputs were created?
- What validation passed?
- What author decision or fan-in action remains?

These rules are binding for all downstream manuscript-revision workers unless a
later human-authored WG task description explicitly overrides them.
