## Title
Implicit interval tree — the data structure under the implicit pangenome graph

## Bullets
- Each pairwise alignment becomes an interval `[start, end)` on its target sequence (CIGAR + target range)
- Li & Rong (2020) "implicit interval tree": store intervals in a sorted array; tree is implied by index, not pointers — O(log n) overlap queries with near-zero memory overhead
- Build one interval tree per target sequence; the union over all sequences is an **interval forest**
- That forest, queried by transitive closure, **is** our implicit pangenome graph — no graph is ever materialized explicitly
- Foundation for the rest of the talk: every plot you see comes from queries against this structure

## Primary figure
`slides/20260204_Subtelomics_overview_EG.pdf` (page 2 — kept verbatim from v1):
horizontal interval plot of 10 example alignment intervals; explicit-tree diagram on the left showing nodes labelled `[Start, End) Index, MaxEnd` across levels 0–3, with a highlighted node `[300, 320)` linked back to its alignment record; compact ordered-array form on the right showing the same intervals as the implicit representation (no pointers, just a sorted run with implied tree structure).

## Speaker notes
One slide of plumbing before the biology — the rest of the talk is just queries against this object.

An alignment is an interval on a target sequence: start, end, CIGAR. Li & Rong (2020) showed you can index a set of intervals as an *implicit* tree — a sorted array where tree topology is recovered from the index, not stored as pointers. O(log n) overlap queries, essentially zero memory beyond the intervals.

We build one tree per target sequence. The union is an interval forest. Walked by transitive closure, that forest **is** the implicit pangenome graph — the same object a `pggb` graph represents, never materialized. The implementation is **IMPG** (https://github.com/pangenome/impg, locally `~/impg`); `impg query -x` does the transitive lookup that powers everything downstream.

So when a later slide says "we queried the pangenome at chr18 q-arm," what's literally happening is `impg query` walking these trees. Keep this picture in your head — it's the substrate for every plot that follows.

## Time budget
50s

## Notes for synthesizer
- Slide 01 (title) sets up "we surveyed subtelomeres at HPRC v2 scale"; this slide answers *how is that even tractable* before slide 03 walks through the IMPG workflow figure. The transition from 01 → 02 → 03 is: motivation → data structure → pipeline.
- Slide 03 already shows the per-sequence interval-tree stack as part of the IMPG workflow; this slide deliberately *introduces* that visual vocabulary so 03 reads as "now zoom out to the full pipeline" rather than re-explaining trees.
- The "implicit pangenome graph" phrasing here matches **C2** of the canonical abstract (`paper_prep/synthesis/ABSTRACT.md`) — the speaker should land that exact phrase so it threads through the rest of the deck.
- Callback opportunity: when slide 04 (HPRC query) or slide 07 (all-vs-all heatmap) appears, the speaker can gesture back to "every cell is an `impg query`."
- Figure is unchanged from v1 — no new R/ggplot2 work needed. If a future revision wants a clearer ordered-array highlight, that is a minor v3 polish, not blocking for BoG.
- IMPG citation: Garrison et al., `pangenome/impg` (https://github.com/pangenome/impg). Underlying algorithm: Cordes, Li & Rong, "cgranges: a C/C++ library for fast interval overlap queries", 2020 (the implicit-interval-tree formulation).
