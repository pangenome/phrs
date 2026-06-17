# Manuscript Revision 01 Fanout Graph

Date: 2026-06-17

Scope: graph validation and dependency sanity only. I read
`wg_manuscript_revision_prompt.md`, `paper_prep/manuscript_revision/00_inventory.md`,
and the current `manuscript-revision-*` WG task descriptions. I did not rewrite
scientific text, run manuscript builds, run R/Python analyses, or submit Slurm
jobs.

## Executive Summary

The current fanout graph represents every prompt item A-1 through G-3. No new WG
tasks were required during this validation pass. The graph already has the
correct high-level shape:

1. `manuscript-revision-01-generate` locks this graph.
2. `manuscript-revision-02-operating` writes shared operating and Slurm rules.
3. Cluster workers write audits, compute reports, or decision records.
4. Cluster fan-ins synthesize the worker outputs.
5. `manuscript-revision-g0` drafts the abstract/introduction/title decision
   package only after A, B/F, C/D, and E evidence are available.
6. `manuscript-revision-final-fanin` creates the final edit plan and decision
   register.
7. `manuscript-revision-paper-patch` is the guarded manuscript-edit node.
8. `manuscript-revision-qa` compiles and audits the final manuscript-facing
   outputs.

No unsafe J-task direct edits were found in the current descriptions. Judgment
tasks are decision-record only unless the downstream guarded patch task has a
named author decision or a prior decision record selecting the wording. Heavy
compute is split into explicit Slurm-capable tasks and is not assigned to the
fan-in or QA nodes.

## Graph Edits Made

No `wg add`, dependency, or task-description edits were made by this task.

Reason: each prompt item is already represented, direct manuscript editing is
reserved for mechanical fan-in and the guarded patch node, and heavy analyses
are isolated in compute-aware tasks. The only correction needed is documentary:
downstream workers should treat `manuscript-revision-c0c` as the D-1 task even
though the task ID does not include `d1`.

## Global Authority Rules

- Mechanical tasks may audit, propose exact fixes, and in the explicit
  mechanical fan-in node apply non-judgment fixes.
- Judgment tasks must write decision records or decision packages. They must not
  edit `submission/paper.tex`, `submission/bibliography.bib`, figure captions,
  or figure assets unless a named author decision has already selected the
  change.
- Compute tasks may write scripts, Slurm submission files, logs, result tables,
  plots, and reports under `paper_prep/manuscript_revision/`.
- Heavy data processing must use Slurm/compute nodes. Do not stream the full
  sequence-level similarity matrix, load large RDS objects, or run multi-core
  permutations on the head node.
- Fan-in tasks synthesize existing worker artifacts. They may do lightweight
  table/markdown synthesis, but should not become hidden monolithic analysis
  jobs.
- Manuscript edits belong to `manuscript-revision-a-fanin` only for supported
  mechanical fixes, and to `manuscript-revision-paper-patch` for the final
  guarded patch.

## Dependency Graph By Cluster

### Cluster A: Mechanical Reconciliation

Prompt items represented:

| Prompt item | WG task | Class | Manuscript edit authority |
|---|---|---|---|
| A-1 sample counts | `manuscript-revision-a1` | Mechanical audit | No direct edits; writes `A1_sample_counts.md` |
| A-2 duplicate references | `manuscript-revision-a2` | Mechanical audit | No direct edits; writes `A2_bibliography_audit.md` |
| A-3 construction artifacts | `manuscript-revision-a3` | Mechanical audit | No direct edits; writes `A3_artifacts_audit.md` |
| A-4 DOI verification | `manuscript-revision-a4` | Mechanical verification | No direct edits; writes `A4_guarracino_doi.md` |
| A-5 XTR acronym | `manuscript-revision-a5` | Mechanical audit | No direct edits; writes `A5_A7_mechanical_audit.md` |
| A-6 empty sections | `manuscript-revision-a5` | Mechanical audit | No direct edits; included with A-5/A-7 |
| A-7 symbol collisions | `manuscript-revision-a5` | Mechanical audit | No direct edits; included with A-5/A-6 |

Fan-in:

- `manuscript-revision-a-fanin` depends on `a1`, `a2`, `a3`, `a4`, `a5`, and
  `02-operating`.
- This is the first node authorized to apply clearly supported mechanical
  corrections to `submission/paper.tex` or `submission/bibliography.bib`.
- It must not resolve J-marked framing questions.
- It feeds `manuscript-revision-g0` and `manuscript-revision-final-fanin`.

Status: complete representation. No deferred A item.

### Cluster B: 3D-Contact Control Apparatus

Prompt items represented:

| Prompt item | WG task | Class | Manuscript edit authority |
|---|---|---|---|
| B-1 lead with flanking unique-sequence evidence | `manuscript-revision-b1` | Judgment record | Decision-record only |
| B-2 reframe MAPQ disclosure | `manuscript-revision-b1` | Judgment record | Decision-record only |
| B-3 label PHR-internal scatters downstream | `manuscript-revision-b1` | Judgment record | Decision-record only |
| B-4 pointwise p-values vs Mantel | `manuscript-revision-b4` | Mechanical-to-Judgment audit | No direct edits; recommendations only |
| B-5 apparatus essentiality | `manuscript-revision-b5` | Judgment support | Decision-record only |

Supporting source inventory:

- `manuscript-revision-b0` inventories human 3D inputs, source tables, scripts,
  statistics, and whether each item is primary evidence or defensive apparatus.
- `b1`, `b4`, and `b5` all depend on `b0`, so their recommendations should be
  anchored in the same evidence inventory.

Fan-in:

- `manuscript-revision-bf-fanin` depends on `b1`, `b4`, `b5`, `f1`, `f3`, and
  `02-operating`.
- This fan-in synthesizes the 3D-contact, Mantel/p-value, mouse, and FST/cM-Mb
  outputs into `BF_3d_contact_synthesis.md`.
- It feeds `manuscript-revision-g0` and `manuscript-revision-final-fanin`.

Status: complete representation. No deferred B item. B tasks are safe because
they write inventories, audits, and decision records rather than applying
scientific rewrites directly.

### Cluster C: Communities, Continuum, NJ Tree, Bootstrap

Prompt items represented:

| Prompt item | WG task | Class | Manuscript edit authority |
|---|---|---|---|
| C-0 characterize the continuum | `manuscript-revision-c0a`, `manuscript-revision-c0b`, `manuscript-revision-c0c` | Compute + Judgment support | No direct edits; reports only |
| C-1 tree essentiality | `manuscript-revision-c1` | Judgment record | Decision-record only |
| C-2 0.1% bootstrap question | `manuscript-revision-c2` | Code/data audit, Judgment support | No direct edits |
| C-3 q-arm list as illustration | `manuscript-revision-c3` | Judgment-support text audit | Decision-record only |

C-0 split:

- `manuscript-revision-c0a` is the arm-level, compute-light characterization
  using the 41 x 41 arm matrix. It is allowed on the head node because the
  matrix is tiny.
- `manuscript-revision-c0b` is the full sequence-level continuum analysis on
  the 15,668 x 15,668 object. It must use Slurm if heavy processing is needed.
- `manuscript-revision-c0c` locates Leiden resolution scans and handles the
  sampling-stability/D-1 bridge.

Dependencies:

- `c1` depends on `c0a`, so the tree-essentiality record can incorporate the
  arm-level continuum evidence.
- `c3` depends on both `c0a` and `c0b`, so q-arm list language is conditional
  on the heatmap-density/continuum evidence.
- `c2` is independent of C-0 because it audits the existing bootstrap code and
  cache. It still feeds the C/D fan-in.

Fan-in:

- `manuscript-revision-cd-fanin` depends on `c0a`, `c0b`, `c0c`, `c1`, `c2`,
  `c3`, and `02-operating`.
- It writes `CD_continuum_community_synthesis.md` and feeds `g0` and
  `final-fanin`.

Status: complete representation. No deferred C item. C tasks are safe because
the potentially dangerous scientific changes are held for decision records and
later guarded patching.

### Cluster D: Sampling Justification

Prompt item represented:

| Prompt item | WG task | Class | Manuscript edit authority |
|---|---|---|---|
| D-1 state sampling guarantee once and correctly | `manuscript-revision-c0c` | Compute/planning + Judgment support | No direct edits |

`manuscript-revision-c0c` explicitly separates Leiden resolution dependence
from sampling stability and writes one honest D-1 manuscript statement. It is
then synthesized by `manuscript-revision-cd-fanin`.

Status: represented, not missing. The item is intentionally folded into C-0
because the same scans and stability checks support both the continuum claim
and the sampling caveat.

### Cluster E: Pedigree

Prompt items represented:

| Prompt item | WG task | Class | Manuscript edit authority |
|---|---|---|---|
| E-1 de-circularize | `manuscript-revision-e1` | Judgment support | Decision-record only |
| E-2 right-size claim/title | `manuscript-revision-e1` | Judgment support | Decision-record only |
| E-3 bound artifact | `manuscript-revision-e1` | Judgment support | Decision-record only |

Dependencies:

- `e1` depends on `01-generate` and `02-operating`.
- `e1` feeds `g0` and `final-fanin`.

Status: complete representation. No deferred E item. It is safe because it
writes `E_pedigree_audit.md` and does not edit the manuscript.

### Cluster F: Orphan Defenses And Vestigial Elements

Prompt items represented:

| Prompt item | WG task | Class | Manuscript edit authority |
|---|---|---|---|
| F-1 FST with no analysis | `manuscript-revision-f1` | Mechanical-to-Judgment audit | No direct edits |
| F-2 cM/Mb line | `manuscript-revision-f1` | Mechanical-to-Judgment audit | No direct edits |
| F-3 mouse Mantel double-reporting and curve shape | `manuscript-revision-f3` | Compute + Judgment support | No direct edits |

Fan-in:

- `f1` and `f3` feed `manuscript-revision-bf-fanin`.
- F is joined with B because the mouse curve, Mantel reporting, and orphan FST
  issues affect the same 3D/contact-evidence synthesis and abstract/title
  strength.

Status: complete representation. No deferred F item.

### Cluster G: Title, Abstract, Introduction Alignment

Prompt items represented:

| Prompt item | WG task | Class | Manuscript edit authority |
|---|---|---|---|
| G-0 abstract/introduction two-tier continuum | `manuscript-revision-g0` | Judgment package | Decision-package only |
| G-1 title vs evidence | `manuscript-revision-g0` | Judgment package | Decision-package only |
| G-2 "unorthodox recombination" framing | `manuscript-revision-g0` | Judgment package | Decision-package only |
| G-3 denominator clarity | `manuscript-revision-g0` | Judgment package | Decision-package only |

Dependencies:

- `g0` depends on `a-fanin`, `bf-fanin`, `cd-fanin`, `e1`, and
  `02-operating`.
- This correctly makes G last among scientific decision packages, because the
  abstract, introduction, title, and denominator sentence must reflect what
  survives in A-F.
- `g0` feeds `final-fanin`.

Status: complete representation. No deferred G item.

## Mechanical, Judgment, And Compute Classification

Mechanical-only or primarily mechanical:

- `manuscript-revision-a1`
- `manuscript-revision-a2`
- `manuscript-revision-a3`
- `manuscript-revision-a4`
- `manuscript-revision-a5`
- `manuscript-revision-a-fanin` for supported mechanical edits only
- `manuscript-revision-b0`
- `manuscript-revision-qa`

Mechanical-to-Judgment audits:

- `manuscript-revision-b4`
- `manuscript-revision-f1`
- `manuscript-revision-c2`

Judgment decision records or decision packages:

- `manuscript-revision-b1`
- `manuscript-revision-b5`
- `manuscript-revision-c1`
- `manuscript-revision-c3`
- `manuscript-revision-e1`
- `manuscript-revision-g0`
- `manuscript-revision-final-fanin`

Compute tasks:

- `manuscript-revision-c0a`: compute-light, head-node allowed.
- `manuscript-revision-c0b`: heavy sequence-level compute, Slurm required for
  matrix scans or large object loads.
- `manuscript-revision-c0c`: compute/planning; lightweight cached scans allowed
  locally, large sampling/stability work requires Slurm.
- `manuscript-revision-f3`: compute/statistics; cached small-table work may be
  local, nontrivial permutation/bootstrap work requires Slurm.

Coordination/operating tasks:

- `manuscript-revision-01-generate`
- `manuscript-revision-02-operating`
- `manuscript-revision-paper-patch`

## Manuscript Edit Authority

May edit manuscript-facing files:

- `manuscript-revision-a-fanin`: only non-judgment mechanical corrections that
  are clearly supported by A1-A5 artifacts.
- `manuscript-revision-paper-patch`: final guarded patch after
  `final_revision_package.md`; may apply J-marked changes only when a prior
  decision record contains a named author decision or otherwise explicitly
  selected wording.
- `manuscript-revision-qa`: should not make scientific edits; may only perform
  final QA/reporting unless a trivial build/audit fix is explicitly in scope.

Must not edit manuscript-facing files:

- All A worker audits: `a1`, `a2`, `a3`, `a4`, `a5`.
- All B workers: `b0`, `b1`, `b4`, `b5`.
- All C/D workers: `c0a`, `c0b`, `c0c`, `c1`, `c2`, `c3`, `cd-fanin`.
- All E/F workers: `e1`, `f1`, `f3`, `bf-fanin`.
- `g0` and `final-fanin` write decision packages and edit plans, not the
  manuscript itself.

## Slurm And Head-Node Safety

General rules:

- Head node is acceptable for `rg`, metadata inspection, small TSV summaries,
  markdown writing, and tiny matrices such as 41 x 41 arm-level Jaccard data.
- Head node is not acceptable for streaming the 10+ GB full sequence-level
  similarity TSV, loading large RDS objects, full-matrix permutations, or
  multi-core bootstrap/permutation jobs.
- Slurm tasks must record the exact command, job ID, partition/node class,
  resources, stdout path, stderr path, wall time, exit status, and output
  artifacts in their report.
- If Slurm is unavailable, the task should leave a ready `sbatch` script and a
  clear not-run reason rather than running heavy work locally.

Recommended resource shapes:

| Task | Safe local work | Slurm-required work | Recommended Slurm shape |
|---|---|---|---|
| `c0a` | 41 x 41 arm matrix summaries, small plots | None expected | No Slurm expected; if using R, 1 CPU, 4-8 GB, <30 min is enough |
| `c0b` | `ls -lh`, manifest inspection, script writing, cached small metadata | Any scan/load of the full 15,668 x 15,668 similarity object or 10+ GB TSV | 1 node, 4-8 CPUs, 64-128 GB RAM, 4-12 h, local scratch if available; start conservative and log actual peak memory |
| `c0c` / D-1 | Locate existing scan tables, summarize small cached matrices | New sequence-level subsampling, repeated Leiden scans, or large matrix stability checks | 1 node, 4-8 CPUs, 32-96 GB RAM, 2-8 h depending on cached input size |
| `f3` | Read cached mouse stage summary TSVs, write contrast script, run tiny summaries | Nontrivial permutation/bootstrap over PHR-pair contacts or large mouse contact tables | 1 node, 4-8 CPUs, 16-64 GB RAM, 1-6 h; increase only if cached table sizes require it |

Specific safety notes:

- `c0b` must not `zcat` or stream the full sequence-level similarity TSV on the
  head node. It should first inspect file sizes and existing caches, then submit
  Slurm if computation is needed.
- `c0c` must not let Leiden-vs-UPGMA agreement substitute for sampling
  stability. If true sampling stability requires large inputs, it should write
  Slurm-ready scripts and clearly separate not-run plans from completed
  evidence.
- `f3` must test the shape claim directly, especially zygotene versus flanking
  stages, rather than reporting only per-stage `rho != 0` p-values. Heavy
  permutation tests belong on Slurm.

## Explicit Fan-In And Synthesis Dependencies

Cluster-level fan-ins:

- `manuscript-revision-a-fanin`
  - After: `a1`, `a2`, `a3`, `a4`, `a5`, `02-operating`.
  - Output: `A_mechanical_fixes.md`.
  - Feeds: `g0`, `final-fanin`.

- `manuscript-revision-bf-fanin`
  - After: `b1`, `b4`, `b5`, `f1`, `f3`, `02-operating`.
  - Indirectly after `b0` through `b1`, `b4`, and `b5`.
  - Output: `BF_3d_contact_synthesis.md`.
  - Feeds: `g0`, `final-fanin`.

- `manuscript-revision-cd-fanin`
  - After: `c0a`, `c0b`, `c0c`, `c1`, `c2`, `c3`, `02-operating`.
  - Output: `CD_continuum_community_synthesis.md`.
  - Feeds: `g0`, `final-fanin`.

Cross-cluster decision package:

- `manuscript-revision-g0`
  - After: `a-fanin`, `bf-fanin`, `cd-fanin`, `e1`, `02-operating`.
  - Output: `G_abstract_intro_title_package.md`.
  - Feeds: `final-fanin`.

Final integration:

- `manuscript-revision-final-fanin`
  - After: `a-fanin`, `bf-fanin`, `cd-fanin`, `e1`, `g0`.
  - Output: `final_revision_package.md`.
  - Feeds: `paper-patch`.

- `manuscript-revision-paper-patch`
  - After: `final-fanin`.
  - Output: `manuscript_patch_report.md`.
  - Feeds: `qa`.

- `manuscript-revision-qa`
  - After: `paper-patch`.
  - Output: `qa_report.md`.

## Representation Checklist

All prompt items are represented:

- A-1: represented by `a1`.
- A-2: represented by `a2`.
- A-3: represented by `a3`.
- A-4: represented by `a4`.
- A-5: represented by `a5`.
- A-6: represented by `a5`.
- A-7: represented by `a5`.
- B-1: represented by `b1`.
- B-2: represented by `b1`.
- B-3: represented by `b1`.
- B-4: represented by `b4`.
- B-5: represented by `b5`.
- C-0: represented by `c0a`, `c0b`, and `c0c`.
- C-1: represented by `c1`.
- C-2: represented by `c2`.
- C-3: represented by `c3`.
- D-1: represented by `c0c` and synthesized in `cd-fanin`.
- E-1: represented by `e1`.
- E-2: represented by `e1`.
- E-3: represented by `e1`.
- F-1: represented by `f1`.
- F-2: represented by `f1`.
- F-3: represented by `f3`.
- G-0: represented by `g0`.
- G-1: represented by `g0`.
- G-2: represented by `g0`.
- G-3: represented by `g0`.

No item is intentionally deferred outside this revision graph. The only
conditional deferrals are inside worker tasks: if heavy compute cannot be run,
the worker must leave Slurm-ready scripts, exact commands, and a not-run reason;
if a J task lacks a named author decision, it must leave a decision point rather
than editing manuscript text.

## Residual Risks

- `wg_manuscript_revision_prompt.md` is untracked in the working tree. This is
  a provenance risk for the revision brief, but it is outside this task's edit
  scope. Downstream work should continue to cite this prompt path until the
  coordinator decides whether to commit it.
- D-1 is represented by `c0c`, not by a task ID containing `d1`. This is safe
  but easy to miss; fan-in documents should mention the mapping explicitly.
- `final-fanin` does not depend directly on `02-operating`, but all worker and
  synthesis paths feeding it do. This is acceptable.
- The graph relies on downstream workers respecting decision-record boundaries.
  `02-operating` should restate those boundaries before workers proceed.
