# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

<!-- workgraph-managed -->
# Workgraph

Use workgraph for task management.

**At the start of each session, run `wg quickstart` in your terminal to orient yourself.**
Use `wg service start` to dispatch work — do not manually claim tasks.

## For All Agents (Including the Orchestrating Agent)

CRITICAL: Do NOT use built-in TaskCreate/TaskUpdate/TaskList/TaskGet tools.
These are a separate system that does NOT interact with workgraph.
Always use `wg` CLI commands for all task management.

CRITICAL: Do NOT use the built-in **Task tool** (subagents). NEVER spawn Explore, Plan,
general-purpose, or any other subagent type. The Task tool creates processes outside
workgraph, which defeats the entire system. If you need research, exploration, or planning
done — create a `wg add` task and let the coordinator dispatch it.

ALL tasks — including research, exploration, and planning — should be workgraph tasks.

### Orchestrating agent role

The orchestrating agent (the one the user interacts with directly) does ONLY:
- **Conversation** with the user
- **Inspection** via `wg show`, `wg viz`, `wg list`, `wg status`, and reading files
- **Task creation** via `wg add` with descriptions, dependencies, and context
- **Monitoring** via `wg agents`, `wg service status`, `wg watch`

It NEVER writes code, implements features, or does research itself.
Everything gets dispatched through `wg add` and `wg service start`.

## Useful workgraph commands

- **`wg publish --wcc <task-id>`** — Publish every task in the weakly-connected component of TASK in one call (treats deps as undirected, unpauses the whole fan-out + synthesizer subgraph in topological order). **Use this for diamond-pattern dispatches**: don't loop `wg publish` over N draft tasks — each call writes to graph.jsonl and on moosefs the locking adds ~5s per call. One `--wcc` invocation publishes the entire connected component instantly. Seed task can be any paused node in the component.
- **`wg kill <agent-id>` + `wg abandon <task-id>`** — Stop a running task that's pursuing the wrong shape/strategy. Kill the agent process first (releases the worktree lock), then abandon the task with a `--reason`.
- **Untracked files: invisible at relative paths in worktrees, but readable via absolute paths.** Agent worktrees (`.wg-worktrees/agent-NNN/`) are fresh git checkouts — anything not in `git ls-files` does not exist at the worktree-relative path. BUT agents have `Read`/`Bash` access to the filesystem, so an absolute path like `/moosefs/erikg/phrs/slides/foo.pdf` resolves to the main checkout copy and they can read it that way. Practical implications:
  - Slide / figure / data tasks that consume the file via absolute path: **work** (often correctly) even if untracked.
  - Tasks that need the file to be *embedded* in their worktree's commit (rendering pipelines that include figures, integrator/synthesizer tasks that do `ls path/to/file`, anything that copies into the worktree): **fail or flag missing**.
  - Symptom of this failure mode: `figure_manifest.md`-style outputs marking inputs as "MISSING — not in worktree" while the file exists for you in the main checkout.
  - Fix: `git add` and commit critical inputs before dispatching downstream rendering / integration tasks. If you `Write` a file intending it to anchor downstream tasks, commit it in the same turn. Common culprits: PDFs, screenshots, slide decks, anchor docs (e.g. `ABSTRACT.md`), large data files.
- **`wg publish` returns the IDs of created tasks** in the form `(<id>)` at the end of the line — `sed -n 's/.*(\([a-z0-9-]*\))$/\1/p'` extracts cleanly. Use this to chain `--after` deps when scripting fan-outs.

<!-- WG-managed -->
# WG (project-specific guide)

This file is the **layer-2** project guide for agents working in this
WG project. It is NOT the universal chat-agent / worker-agent
contract — that is bundled inside the `wg` binary and emitted by:

```
wg agent-guide
```

Run `wg agent-guide` at session start (or read its output from a previous
session) to get the universal role contract: chat agent vs dispatcher vs worker
distinction, `## Validation` requirement, smoke-gate, cycle handling, git
hygiene, worktree isolation, "no built-in Task tool" rules, etc.

This file only covers things specific to this project. Add project-specific
build commands, test commands, architecture notes, and service recipes here.

**At the start of each session, run `wg quickstart` in your terminal to orient yourself.**
Use `wg service start` to dispatch work — do not manually claim tasks.

## Project: subtelomeric PHR analysis

Research repo, not a software product. Inter-chromosomal subtelomeric sequence
sharing across 233 HPRCv2 samples (465 haplotypes, 15,668 PHRs, 41 signal-bearing
arms, 15 arm-level / 50 sequence-level communities), plus 3D-genome validation
(Hi-C, Pore-C, CiFi, Dip-C, sperm sc, RPE-1, mouse meiosis). Outputs are a
14-section analysis report, paper figures, and slide decks.

GitHub: `https://github.com/ekg/phrs`. The same checkout exists at
`/moosefs/erikg/phrs/` and is the canonical mirror used by workgraph
(`/moosefs/erikg/phrs/.wg/`).

## Where data lives

Source-of-truth bioinformatics data is NOT in the repo. Scripts read absolute
paths:

- `/moosefs/guarracino/HPRCv2/PHR_III/` — similarity matrices, Leiden/UPGMA
  community assignments, 3D-validation tables, enrichment outputs.
- `/moosefs/erikg/phrs/` — mirror of this repo.

A worktree without moosefs access cannot rebuild figures or run figure scripts.
Pre-rendered upstream artifacts that figures vendor are committed at the repo
root (e.g. `p_genome_wide_identity_heatmap.pdf`, `chm13.phrs.bed`,
`identity_heatmap_chr*.pdf`).

## Layout

- `end-to-end-report/report/01_pipeline.md` … `14_pedigree_recombination.md` —
  14-section analysis report; `end-to-end-report/README.md` is the TOC.
- `paper_prep/figures/{fig1..fig4, ed1..ed5, ed8, nj_tree_arms}/` — one directory
  per paper figure. `ed6` and `ed7` slots are intentionally empty; do not
  create stubs. Each `fig*/`+`ed*/` has a single `figure_<id>.R` (`fig2/`
  uses `figure_fig2.py`), `sources.tsv` (input/output manifest), `caption.md`,
  and rendered PDF+PNG. `nj_tree_arms/` is an exception: script is
  `nj_tree.R`, outputs `nj_tree_annotated.{pdf,png}`, no `sources.tsv` /
  `caption.md`.
- `paper_prep/synthesis/` — paper-level artifacts (`ABSTRACT_BoG.md`,
  `ABSTRACT_nature.md`, `REFERENCES_v3.bib`, `CROSSWALK.md`).
- `paper_prep/lit_review/topic_NN_*.{md,bib}` — 14 per-topic literature
  reviews; `SYNTHESIS_v2.typ` renders the unified review.
- `paper_prep/surveys/SURVEY_NN_*.md` — section-by-section evidence
  inventories. Figure `sources.tsv` lines cite them as `SURVEY_NN §X.Y`.
- `paper_prep/_brainstorming/` — scratch (analysis experiments, validation
  reports); authoritative content is upstream in `paper_prep/figures/` or
  `synthesis/`.
- `slides/v2-review-zoom/` — active deck (BoG 2026). Other `slides/`
  subdirs are earlier rounds (`v2`, `v2-zoom`, `v2-review`) or background
  decks (`chm13-phr-ucsc-browser/`, `meiosis-prophase-background/`).
- `plot-impg-coverage.R` — standalone R script (run upstream) that emits the
  repo-root genome-wide heatmap PNGs/PDFs.
- `scripts/merge_archive_into_graph.py` — wg utility (not analysis code).

## Slide decks (Typst)

Active deck: `slides/v2-review-zoom/_typst/zoom_review_deck.typ` ->
`slides/v2-review-zoom/BoG_2026_review_zoom_v9.pdf`.

Build from `slides/v2-review-zoom/_typst/`:

```
typst compile --root .. zoom_review_deck.typ ../BoG_2026_review_zoom_v9.pdf
typst compile --root .. --ppi 144 zoom_review_deck.typ page-{0p}.png
```

`typst` is at `~/.cargo/bin/typst`. Each round (`v2` … `v9`) is a frozen
artifact: own `REVISION_NOTES_V<N>.md` (v2 is the unsuffixed
`REVISION_NOTES.md`; v3+ follow the `_V<N>` convention), plus per-slide
experiments under `_revision_assets/v<N>/<slide-id>/` (`README.md`,
`SLIDE_PATCH.md` handoff, `make_*.R` source, candidate PDF/PNG). Bump to
`v<N+1>` for a new revision round; do not edit a frozen `v<N>` in place.

## Figures (R / Python)

Each `paper_prep/figures/<id>/` is self-contained: one script reads the
`sources.tsv` inputs and writes `figure_<id>.{pdf,png}` next to it. R scripts
run under guix:

```
guix shell -m /tmp/manifest_fig<id>.scm -- Rscript paper_prep/figures/<id>/figure_<id>.R
```

`guix` lives on the moosefs deploy environment, not on every dev machine. The
exact manifest path is in each script's header comment.

`sources.tsv` schemas are not standardized across figures (columns differ:
`panel/status/role/path/note`, `panel/content/status/source_path/survey_section`,
or `panel/source_path/role/headline_metric`). Common statuses: `READY` (already
on disk, vendored), `GENERATE` (script will produce it). `/moosefs/...` paths
are external inputs.

## Naming and units

- Arms: `chr1_p`, …, `chr22_q`, `chrX_p`, `chrX_q`, `chrY_p`, `chrY_q`
  (48 total, 41 with signal). Some upstream tables use `_parm`/`_qarm`; figure
  scripts normalise before joining.
- Communities: arm-level Leiden k=15 (cross-checked by UPGMA k=14);
  sequence-level Leiden k=50.
- PHR (pseudohomologous region): telomere-anchored 500 kb flank from a contig
  >= 1 Mb that aligned at >= 95% identity to another arm. 15,668 PHRs across
  18,827 flanks.
- CHM13-coordinate PHR BEDs at repo root: `chm13.phrs.bed`,
  `chm13.phrs.no_acro.bed`, `CHM13-HG002.sub-telo-phrs.bed`.

## Commit convention

Commits originate from workgraph worker agents. Messages look like
`feat: <task-slug> (agent-NNN)`. The `(agent-NNN)` suffix is workgraph
provenance for the graph; keep it.

This guide is written to both `CLAUDE.md` and `AGENTS.md` and kept in
lock-step. The two files exist because Claude Code and Codex CLI look for
different filenames, but they should never drift in content. Any divergence is
a bug. Update both together.
