# Workgraph Failure Modes: Copy-Number-Aware Enrichment Analysis

**Date:** 2026-04-01  
**Subgraph:** Copy-number-aware enrichment analysis  
**Timeline:** 14:46 - 19:01 (4h 15m duration)  
**Tasks:** Started ~85, exploded to 286+, 303 completed  
**Agent Churn:** 392 dead agents vs 2 alive (massive churn indicator)  
**Status:** Core work completed, significant system failures requiring manual intervention

## Executive Summary

This document provides feedback to workgraph system developers on failure modes encountered during a copy-number-aware enrichment analysis project. While the core research and methodology work completed successfully, we experienced five distinct system failures that required ~15 minutes of manual coordinator intervention and wasted significant compute resources. The most serious issue was an undetected crash loop that spawned 20+ agents over 30 minutes without circuit breaking.

## Failure Modes

### 1. Task Explosion / Over-Decomposition

**Description:** The `implement-copy-number` task was intended to implement 2-3 enrichment methods. Instead, the assigned agent decomposed it into 150+ subtasks covering theoretical statistics (mathematical formulations, type I error validation, null distribution validation, edge case analysis, parameter constraint validation, etc.). This transformed a practical implementation task into an academic research program.

**Impact:** 
- Graph size exploded from ~85 tasks to 286+ tasks
- Most subtasks were unnecessary for the actual goal
- Significant resource waste on over-engineering

**Root Cause Hypothesis:** 
The agent exhibited an 'autopoietic' decomposition pattern that recursively creates subtasks without a bound on depth or breadth. No guardrail prevented a single task from spawning dozens of children. This appears to be an emergent behavior where agents optimize for thoroughness over practicality.

**Suggested Fix:** 
- Consider max subtask limits per task (e.g., 10-20 direct children)
- Require coordinator approval for decompositions beyond N subtasks
- Add cost estimation for task decomposition patterns
- Implement depth limits on recursive task creation

**Timeline:** ~15:10-15:35

---

### 2. Claude CLI Crash Cascade (~16:02)

**Description:** Around 16:02, a wave of simultaneous Claude CLI failures hit multiple task types. This was likely a transient API rate limit or outage, but the cascade was disproportionate:
- All `.flip-*` tasks failed: `'FLIP inference LLM call failed — Claude CLI call failed (exit Some(1))'`
- All `.evaluate-*` tasks failed: `'Evaluation LLM call failed — Claude CLI call failed (exit Some(1))'`
- Multiple agent tasks failed with `'Agent exited with code 1'`

**Impact:** 
- ~25 tasks failed simultaneously
- Required manual retry of each failed task
- No automatic recovery despite likely transient cause

**Root Cause Hypothesis:** 
Likely a transient API rate limit or service outage. The eval/FLIP tasks lack built-in retry logic for transient API failures. A single API hiccup permanently fails every in-flight evaluation task.

**Suggested Fix:** 
- Add automatic retry with exponential backoff for Claude CLI failures in eval/FLIP tasks
- Implement retry logic: 3 attempts with increasing delays before marking as failed
- Distinguish between permanent failures (bad auth) and transient ones (rate limits)
- Consider circuit breaker patterns for API dependencies

**Timeline:** 16:02

---

### 3. Eval Tasks Racing Ahead of Parent Tasks

**Description:** After mass retry of failed tasks, several `.evaluate-*` tasks failed with: 'Task X has status Open — must be done or failed to evaluate'. The evaluation tasks were scheduled and executed before their parent work tasks had completed.

**Impact:** 
- Additional spurious failures requiring manual cleanup
- Evaluation tasks wasted compute on incomplete work
- Coordination between parent and evaluation tasks broken

**Root Cause Hypothesis:** 
When a parent task is retried, it resets to 'open' status, but already-queued evaluation tasks don't get re-blocked. They execute, find the parent isn't done, and immediately fail.

**Suggested Fix:** 
- Evaluation tasks should automatically re-block themselves if their parent task is not in a terminal state (done/failed)
- Implement dynamic dependency checking at task execution time
- Consider making evaluation tasks more resilient to parent state changes

**Timeline:** 16:02-16:03 (immediately after crash cascade)

---

### 4. Graph Context Scope Crash Loop (`integration-testing-with`)

**Description:** The task `integration-testing-with` had `context_scope: graph`. With 300+ tasks in the graph, loading the full context caused every spawned agent to crash within ~60 seconds. The coordinator kept respawning agents (20+ times over 30 minutes), each one dying immediately. This created a silent resource waste pattern - the system would continuously retry without detection or circuit breaking.

**Impact:** 
- Wasted compute: 20+ agent spawns over 30 minutes
- Blocked downstream tasks for 30+ minutes
- No circuit breaker intervention - silent resource waste
- Task remained stuck until manual abandonment

**Root Cause Hypothesis:** 
`context_scope: graph` loads ALL task descriptions/logs into the agent's context. At 300+ tasks, this exceeded context limits or token budgets, causing immediate OOM/crash conditions.

**Suggested Fix:** 
- **Critical:** Add circuit breaker - if an agent crashes N times (e.g., 3) on the same task, pause the task and alert
- Cap graph context loading: max 50 most relevant tasks, or use summarization instead of full dump
- Warn users when setting `context_scope: graph` on large graphs (>100 tasks)
- Implement smart context selection algorithms for large graphs

**Timeline:** 18:30-19:00 (30-minute undetected waste)

---

### 5. Verification System Command Execution Bug

**Description:** Task `type-i-error` was blocked by a verification bug where the verify criteria `'Type I error rates within 1% of nominal α levels; simulation results documented; false positive rate validation complete'` was being executed as a shell command instead of being evaluated as human-readable acceptance criteria.

**Impact:** 
- Task appeared failed despite all work being completed (commit 63d8c5b)
- False negative blocking legitimate task completion - task was actually complete
- Manual intervention required to identify and resolve the system misinterpretation

**Root Cause Hypothesis:** 
The verification system doesn't distinguish between machine-checkable commands and human-readable acceptance criteria. It attempts to execute all verification text as shell commands.

**Suggested Fix:** 
- Parse verification criteria intelligently:
  - If it looks like a shell command (starts with known binary, contains pipes, etc.), execute it
  - Otherwise, use LLM evaluation against the criteria text
- Provide clear syntax for distinguishing command vs. criteria verification
- Add validation of verification criteria at task creation time

**Timeline:** Task blocked throughout execution period

---

## Timeline Summary

| Time | Event | Type |
|------|-------|------|
| 14:46 | User requests copy-number-aware enrichment methods | Start |
| ~15:00 | Research task completes, implementation task begins | Normal |
| 15:10-15:35 | Implementation agent creates 150+ subtasks | **Failure #1** |
| 16:02 | Claude CLI crash cascade hits ~25 tasks | **Failure #2** |
| 16:02-16:03 | Eval tasks race ahead of incomplete parents | **Failure #3** |
| 18:10 | User notices failures, coordinator initiates retries | Recovery |
| 18:12 | Additional retries required | Recovery |
| 18:30-19:00 | `integration-testing-with` crash loops 20+ times | **Failure #4** |
| 19:01 | Coordinator abandons crash-looping task, manual cleanup | Manual intervention |

## Overall Assessment

### What Worked
- Core work completed successfully: methodology research, recommendations, and synthesis
- Graph demonstrated self-healing capabilities 
- Most failures were in validation/testing infrastructure, not primary work
- Workgraph's parallel execution model enabled progress despite failures

### What Failed
- **Systemic Risk:** Crash loop without circuit breaker (#4) - wastes resources silently
- **Reliability:** Cascade failures from external API issues (#2, #3)
- **Scalability:** Context loading breaks down at scale (#4)
- **Usability:** Verification system confusion (#5)
- **Resource Management:** No bounds checking on task decomposition (#1)

### Manual Intervention Required
- ~15 minutes of coordinator cleanup time
- Manual identification and abandonment of crash-looping task
- Manual retry of 25+ failed tasks
- Manual resolution of verification system bug

## Recommendations for Workgraph Development

### Critical (Immediate Implementation Needed)
1. **Circuit Breaker for Task Crash Loops** - If an agent crashes N times on the same task, pause and alert (most dangerous failure - wastes resources silently)
2. **Context Size Limits and Warnings** - Cap `context_scope: graph` in large graphs, warn when >100 tasks

### High (Next Release)  
3. **Task Decomposition Limits** - Max subtask limits per task (e.g., 10-15 direct children), require coordinator approval beyond N
4. **API Failure Retry Logic** - Exponential backoff for Claude CLI failures in eval/FLIP tasks (retry 3x before failing)
5. **Dynamic Dependency Checking** - Re-validate dependencies before task execution, eval tasks should re-block if parent not terminal

### Medium (Ongoing Improvement)
6. **Intelligent Verification Parsing** - Parse verify criteria: if looks like shell command (starts with binary, pipes, etc.) execute it, otherwise use LLM evaluation
7. **Agent Lifecycle Monitoring** - Better detection and handling of dead agents, 392 dead vs 2 alive indicates major churn issue
8. **Smart Context Selection** - For large graphs, summarize context instead of full dump, or select most relevant subset

## Pattern Analysis

The failures fall into three categories:

1. **Scalability Issues** (#1, #4): Systems that work fine at small scale break down with large graphs
2. **External Dependencies** (#2, #3): Cascade failures from external API issues  
3. **Configuration Bugs** (#5): System misinterpretation of user configuration

**Pattern Keywords Observed:**
- **Autopoietic decomposition** - Recursive task creation without bounds
- **Cascade failures** - Single points of failure affecting multiple systems
- **Context explosion** - Unbounded context growth causing crashes
- **Race conditions** - Tasks executing before dependencies satisfied
- **Circuit breaker gaps** - Missing protection against runaway processes

**Biggest Systemic Risk:** 
Failure #4 (crash loop without circuit breaker) represents the highest systemic risk - it wastes resources silently and would go unnoticed without human monitoring.

The workgraph system demonstrates strong self-organizing capabilities but needs better guardrails to prevent pathological behaviors like unbounded decomposition and silent resource waste.

## Appendix: Key Statistics

- **Initial Graph Size:** ~85 tasks
- **Final Graph Size:** 286+ tasks (final count: 303 completed)
- **Task Explosion Ratio:** 3.4x growth from single decomposition
- **Agent Churn:** 392 dead agents vs 2 alive (196:1 ratio - massive churn indicator)
- **Simultaneous Failures:** 25 tasks (16:02 cascade)
- **Crash Loop Duration:** 30 minutes undetected resource waste
- **Agent Respawns:** 20+ spawns for single failing task
- **Manual Intervention Time:** ~15 minutes coordinator cleanup
- **Service Uptime:** 22h 14m (system remained operational despite failures)
- **Core Work Completion:** Successful despite infrastructure failures