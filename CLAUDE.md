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
- **CRITICAL: Untracked files are invisible to agent worktrees.** Agent worktrees are fresh git checkouts — anything not in `git ls-files` is effectively missing from the agent's view, even if it exists on disk in the main checkout. Symptoms: agents flag inputs as "MISSING — not in worktree" when they're plainly there for you. Fix: `git add` and commit the file BEFORE dispatching tasks that depend on it. Common culprits: PDFs, screenshots, slide decks, anchor docs (e.g. `ABSTRACT.md`), large data files. If you write a file with the Write tool intending it to anchor downstream tasks, commit it in the same turn.
- **`wg publish` returns the IDs of created tasks** in the form `(<id>)` at the end of the line — `sed -n 's/.*(\([a-z0-9-]*\))$/\1/p'` extracts cleanly. Use this to chain `--after` deps when scripting fan-outs.
