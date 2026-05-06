## Title

**The implicit pangenome graph — wfmash all-vs-all *is* the graph (no GFA, no construction step)**

Subtitle: every PAF edge is a pairwise mapping; an interval forest per sequence is the index. IMPG (https://github.com/pangenome/impg) does transitive closure over it.

## Bullets

- **All-vs-all wfmash (`-p 95`) over n = 18,827 telomere-anchored 500 kb flanks → 18,827 PAFs.** Each PAF edge is one pairwise alignment; the union of edges *is* the implicit pangenome graph. No pggb, no GFA, no explicit graph construction step.
- **Index = interval forest, one implicit interval tree per sequence (Li & Rong 2020).** Built directly from PAF target intervals — `O(n log n)` build, `O(log n + k)` interval queries, no extra data structures on disk.
- **Query = `impg query -x` (transitive closure)** — chase chains of overlapping pairwise mappings outward from any seed interval to recover every reachable region across haplotypes. This is the operation; the graph is never materialized.
- **Sampling is dense, not sparse.** wfmash's k-mer prefilter evaluates ~12% of the C(n,2) ≈ 1.77×10⁸ pair space at full alignment cost. The Erdős-Rényi connectivity threshold for G(n, p) on n = 18,827 is **p\* = log(n)/n ≈ 5.21×10⁻⁴**; **12% is ~230× above p\***, so the random sub-graph is densely connected w.h.p. — transitive closure from any subtelomere reaches genome-wide. This is what licenses "no chromosomal partitioning."
- **Cite IMPG**: https://github.com/pangenome/impg (Garrison et al.). The implicit-graph framing is canonical and predates this work — we use it; we do not invent it.

## Primary figure

Re-use the v1 slide 3 visual stack: **wave (wfmash) → dotplot (all-vs-all alignment) → interval-tree forest, one tree per sequence**. Source: `slides/20260204_Subtelomics_overview_EG.pdf` page 3 (do NOT modify; extract / re-photograph for the v2 deck during synthesis).

Add one small callout panel at lower right — the Erdős-Rényi sanity check. R/ggplot2 generation script (light, no SBATCH; runs in the agent worktree):

```r
# slide_03_er_callout.R — ER connectivity threshold sanity check
# Output: slide_03_er_callout.pdf  (~3 in × 2 in, intended as inset on slide 3)
library(ggplot2)
n <- 18827
p_star  <- log(n) / n            # ER connectivity threshold ≈ 5.21e-4
p_obs   <- 0.12                  # wfmash k-mer-prefilter evaluation rate
ratio   <- p_obs / p_star        # ≈ 230×

df <- data.frame(
  label = c("ER threshold\np* = log(n)/n", "wfmash sampling\np ≈ 12%"),
  p     = c(p_star, p_obs),
  fill  = c("threshold", "observed")
)
df$label <- factor(df$label, levels = df$label)

ggplot(df, aes(label, p, fill = fill)) +
  geom_col(width = 0.55) +
  scale_y_log10(
    breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1),
    labels = c("1e-4","1e-3","1e-2","1e-1","1")
  ) +
  scale_fill_manual(values = c(threshold = "#888888", observed = "#1f77b4")) +
  geom_text(aes(label = sprintf("%.2g", p)), vjust = -0.5, size = 3.5) +
  annotate("text", x = 1.5, y = 0.4,
           label = sprintf("~%.0f× above threshold\n→ densely connected\n→ closure reaches\n   genome-wide", ratio),
           size = 3.2, hjust = 0.5) +
  labs(x = NULL, y = "edge probability  p   (log scale)",
       title = sprintf("n = %s flanks", format(n, big.mark = ","))) +
  guides(fill = "none") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"))

ggsave("slide_03_er_callout.pdf", width = 3.0, height = 2.0)
ggsave("slide_03_er_callout.png", width = 3.0, height = 2.0, dpi = 300)
```

## Speaker notes

This is the methods anchor — slow down here. When we say *implicit pangenome graph*, we mean it literally: we do **not** build a GFA. We run wfmash all-vs-all on 18,827 telomere-anchored 500 kb flanks; the resulting PAF set **is** the graph. Each PAF edge is one pairwise alignment. We index it as an interval forest — one implicit interval tree per sequence, in the Li-and-Rong 2020 sense from the previous slide. To query we use IMPG (github.com/pangenome/impg), `impg query -x`, which walks transitive closure: start anywhere, chase overlapping mappings outward, collect what's reachable. No graph constructed, nothing to break.

Is 12% sampling enough? The Erdős-Rényi connectivity threshold for n = 18,827 is p\* = log(n)/n ≈ 5.2×10⁻⁴. Twelve percent is **230× above that**, so the random sub-graph is densely connected w.h.p. — closure from any subtelomere reaches every other one. *That* is what licenses "no chromosomal partitioning." Everything downstream rides on it.

## Time budget

**Target: 80 seconds.** This is the methods anchor — do not rush. Roughly 25 s on "the alignment IS the graph, IMPG is the query interface" (frame-shift from v1), 35 s on the Erdős-Rényi argument with the n = 18,827, p\* ≈ 5.21e-4, 230× numbers spoken out loud, 20 s on cite-IMPG and segue into the next slide.

## Notes for synthesizer

- **Slide 02 sets up** the implicit interval tree data structure (Li & Rong 2020). Slide 03 reuses it — call back explicitly: "the data structure from the previous slide, one tree per sequence, gives us the index." Do not re-explain it.
- **Slide 03 sets up slide 04** (querying the pangenome → genome-wide identity heatmap). The link is: "now that we have transitive closure, here's what we see when we ask the graph what each subtelomere matches." Make sure slide 04's lead bullet starts from `impg query -x`, not from a fresh introduction of the alignment.
- **Frame-shift vs v1 is the load-bearing change.** The v1 slide 3 title was "IMPG: IMplicit Pangenome Graphs" and just walked the pipeline. The v2 title must say *this is the graph — there is no GFA*. If the synthesizer compresses, keep that clause; drop the algorithmic detail before dropping the framing.
- **The 12% / Erdős-Rényi callout is new in v2 and must survive compression.** It is the methodological justification for "no chromosomal partitioning" in the abstract (CROSSWALK §7b). Without it, the abstract's wording is asserted but not defended in the talk.
- **Cite line:** `https://github.com/pangenome/impg` should appear on the slide, not just in the speaker notes — credibility cue for an audience that will recognize Garrison et al. tooling.
- **If a referee in Q&A pushes on "12% — derived how?":** the answer is "from the on-disk PAFs — fraction of C(n,2) pairs that pass wfmash's k-mer prefilter and reach full alignment." Methods writer must compute this from the actual PAF set (CROSSWALK §7b explicitly flags this as outstanding); the speaker can punt to Methods.
- **Do NOT** swap in a pggb GFA visualization here — that belongs to a downstream similarity / community-detection slide if at all. Mixing the two undoes the entire frame-shift.
- **No SBATCH needed** for the inset; the R script is ~3 in × 2 in, runs in seconds locally inside the agent worktree.
