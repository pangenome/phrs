# C/D Continuum and Community Synthesis

Date: 2026-06-17
Task: `manuscript-revision-cd-fanin`

## Scope

This fan-in record synthesizes the completed C/D artifacts:

- `paper_prep/manuscript_revision/C0_continuum/C0a_arm_level_report.md`
- `paper_prep/manuscript_revision/C0_continuum/C0b_sequence_level_report.md`
- `paper_prep/manuscript_revision/C0c_D1_resolution_sampling.md`
- `paper_prep/manuscript_revision/C1_tree_essentiality.md`
- `paper_prep/manuscript_revision/C2_bootstrap_audit.md`
- `paper_prep/manuscript_revision/C3_qarm_language.md`

No heavy analysis was run in this fan-in node. The quantitative claims below are
copied from the upstream reports and their generated TSV summaries.

## Decision Summary

The data support a two-tier continuum framing. At both arm and sequence levels,
there is a broad continuum of low-to-intermediate similarity, with locally dense
community neighborhoods that recover known subtelomeric systems. The manuscript
should therefore avoid "discrete blocks" or "closed groups" as the primary
conceptual language. The safer lead phrase is:

> locally dense sequence-similarity communities embedded in a broader
> inter-chromosomal subtelomeric similarity continuum

The q-arm sextet language should be retained only in softened form. The current
manuscript's "linked q-arm group", "previously uncharacterized q-arm grouping",
and "tight {22q, 21q, 19q, 1q, 13q, 17q} q-arm grouping" language should be
replaced with "enriched q-arm neighborhood" or "q-arm neighborhood exemplified
by 22q, 21q, 19q, 1q, 13q and 17q." It should not be presented as a closed
class, clade, or sequence-level stability result.

The community-method language is defensible if it is precise: Leiden community
detection on the fixed arm-level matrix gives 15 arm communities, the chosen
sequence-level Leiden partition is a constrained 50-community operating point,
and UPGMA is an algorithmic/display comparison rather than proof of sampling
stability. NJ and the character bootstrap should not be load-bearing in the
main text. If retained at all, they should be relabeled as exploratory or
surrogate analyses.

## Evidence for the Two-Tier Continuum

### Arm-Level Evidence

C0a used the 41 x 41 arm-level Jaccard distance matrix and the arm-level Leiden
15-community assignment. Similarity is `1 - Jaccard distance`, excluding the
matrix diagonal.

Key arm-level distribution:

| Category | Pairs | Mean similarity | Median similarity | q05 | q95 | Max |
|---|---:|---:|---:|---:|---:|---:|
| all off-diagonal arm pairs | 820 | 0.0685 | 0.0097 | 0.0000 | 0.3567 | 0.9847 |
| within Leiden community | 58 | 0.4506 | 0.4312 | 0.2572 | 0.7993 | 0.9847 |
| between Leiden communities | 762 | 0.0394 | 0.0078 | 0.0000 | 0.1921 | 0.4155 |

Interpretation: within-community arm pairs are much more similar than
between-community arm pairs, but the between-community and all-pair background
are not zero classes. This is the core arm-level support for "two-tier
continuum": dense communities on a lower, variable background.

Named-system examples from `community_similarity_summary.tsv` and
`named_system_peak_similarities.tsv` support the figure narrative directly:

| Community/system | Arms | Mean arm-level similarity | Peak pair |
|---|---|---:|---|
| C1 / D4Z4-DUX4 | 4q, 10q | 0.5830 | 10q-4q = 0.5830 |
| C2 / TUBB8B | 10p, 18p | 0.6202 | 10p-18p = 0.6202 |
| C7 / acrocentric p-arms | 13p, 14p, 15p, 21p, 22p | 0.4428 | 15p-22p = 0.4794 |
| C11 / OR4F core | 5q, 6q | 0.6027 | 5q-6q = 0.6027 |
| C14 / PAR2-like sex q pair | Xq, Yq | 0.9847 | Xq-Yq = 0.9847 |
| C15 / PAR1-like sex p pair | Xp, Yp | 0.7114 | Xp-Yp = 0.7114 |

These values justify saying that the partition recovers known systems and
highlights additional dense neighborhoods. They do not require tree or clade
language.

### Sequence-Level Evidence

C0b scanned the cached 15,668 x 15,668 sequence-level distance object under
Slurm, not on the head node.

Key sequence-level facts:

- Matrix size: 15,668 sequences; 122,735,278 upper-triangle non-self pairs.
- All-pair approximate median similarity: 0.0125.
- All-pair q90: 0.3075; q99: 0.9525.
- Pairs with similarity >= 0.50: 7,182,560 of 122,735,278, or 0.05852.
- Same sequence-community pairs have median similarity 0.9225 and density
  0.88141 at similarity >= 0.50.
- Different sequence-community pairs have median similarity 0.0125 and density
  0.03735 at similarity >= 0.50.

Interpretation: the full sequence object has a broad low-similarity background
and localized high-similarity peaks. It does not support a single clean
threshold that separates all pairs into "related" and "unrelated" classes.
This independently supports the continuum side of the framing while retaining
community-level structure.

## Figure Support

Use the following figures or figure elements to back the revised wording.

| Claim | Figure support | Why this is the right support |
|---|---|---|
| The graph is connected and too tangled to interpret node by node. | Fig. 2A PGGB graph/hairball. | Establishes why arm-level and community summaries are needed. |
| Arm-level similarities are structured but continuous. | Fig. 2B left heatmap, relabeled as UPGMA average-linkage order; C0a `arm_level_similarity_diagnostic.{pdf,png}` if a supplement/revision figure is allowed. | Shows non-uniform structure without requiring phylogeny. |
| Leiden 15-community ordering exposes locally dense neighborhoods. | Fig. 2B right community-ordered heatmap; `submission/fig/MainFigures/Fig2bc_jaccard.*`; `community_blocks.tsv`/colour bands. | This is the load-bearing community figure. |
| Known systems are recovered. | Fig. 2B right heatmap plus Fig. 3 browser examples for C1 4q/10q, C2 10p/18p, and C11 5q/6q. | Heatmap supports the partition; browser panels show sequence architecture. |
| q-arm pattern exists as enriched local density, not a closed class. | Fig. 2B right heatmap plus C0a/C0b density tables. | The visual block and density enrichment support "neighborhood"; they do not support "clade." |
| Sequence-level continuum has high-similarity peaks. | C0b distribution TSVs; optional histogram/table in response or supplement if space permits. | Sequence-level evidence backs the continuum framing more directly than arm-only heatmaps. |

Do not use the NJ tree as a primary figure for the continuum/community claim.
C1 found that the active Fig. 2 heatmap is generated from UPGMA average-linkage
ordering, not NJ, and that no current biological claim needs NJ or bootstrap
support.

## q-Arm Language Decision

### Evidence

C6/q-arm arms: 1q, 13q, 17q, 19q, 21q, 22q.

Arm-level C6 diagnostics:

- C6 within-sextet mean similarity: 0.5003.
- C6 within-sextet median similarity: 0.4338.
- C6-to-non-C6 median similarity: 0.0165.
- C6 within-sextet mean is 15.71x the C6-to-non-C6 mean.
- Exact six-arm set diagnostic over all 4,496,388 six-arm subsets gives
  observed mean 0.5003 versus null mean 0.0685 and null q95 0.1414.

Sequence-level C6 diagnostics:

- Within-C6 sequence-pair density at similarity >= 0.50: 0.5299.
- Outside-C6 different-arm background density at similarity >= 0.50: 0.03733.
- Enrichment over outside-C6 different-arm background at >= 0.50: 14.19x.

These are strong density diagnostics. They are not independent discovery
p-values because C6 was defined from the same similarity structure being
summarized, and they are not a sample-resampling stability test.

### Decision

Retain the q-arm list only as an illustrative enriched neighborhood. Soften in
the abstract, Results, Fig. 2 legend, and any Methods text that mentions the
bootstrap. Do not move the q-arm pattern entirely out of the paper, because both
arm-level and sequence-level summaries support a dense neighborhood. Do move or
delete "tight", "linked group", "previously uncharacterized grouping", "clade",
"phylogeny", and "closed class" language.

Recommended terminology:

- Use: "enriched q-arm neighborhood"
- Use: "exemplified by 22q, 21q, 19q, 1q, 13q and 17q"
- Use: "locally dense region within a broader similarity continuum"
- Avoid: "tight q-arm grouping"
- Avoid: "linked q-arm group"
- Avoid: "closed group", "bounded class", "clade", "monophyletic"
- Avoid: "the q-arm sextet" as a biological noun; reserve it for internal
  shorthand if needed.

## Defensible Community-Method Language

The following statements are defensible:

- "Leiden community detection on the fixed 41-arm Jaccard similarity matrix
  partitions the 41 signal-bearing arms into 15 arm-level communities."
- "The arm-level Leiden resolution scan selected resolution 1.16 by maximum
  mean silhouette; the 15-community state spans resolutions 1.13-1.18."
- "Average-linkage UPGMA on the same arm-level matrix gives 14 clusters and
  recovers the same major biological systems, agreeing exactly with 12 of the
  15 Leiden communities."
- "At sequence level, the reported 50-community partition is a constrained
  k-nearest-neighbor Leiden operating point, k_NN = 75 and resolution = 0.8,
  rather than the unconstrained modularity maximum."
- "UPGMA and resolution scans assess method and resolution dependence on fixed
  matrices; they do not measure sampling stability."

The following statements are not defensible without new Slurm-scale work or
author decision:

- "The Leiden communities are sample-stable."
- "The q-arm neighborhood is a closed sequence-level class."
- "The NJ tree supports a phylogeny of subtelomeric PHRs."
- "The character bootstrap uses the 15,668 x 15,668 sequence-level matrix."
- "The 0.1% q-arm support proves instability of the full sequence-level q-arm
  signal."

## Stability and Subsampling

### What Exists

Arm-level Leiden resolution scan:

- Resolution grid: 0.1-3.0 in 0.01 increments.
- Selected row: resolution 1.16, 15 communities, modularity 0.00714, mean
  silhouette 0.34734.
- The 15-community state appears from resolution 1.13 through 1.18.

Sequence-level Leiden scan:

- k_NN grid: 10, 25, 50, 75, 100, 125.
- Resolution grid: 0.1-3.0 in 0.1 increments.
- Unconstrained modularity maximum: k_NN = 10, resolution = 1.1, 148
  communities, modularity 0.98364.
- Selected constrained operating point: k_NN = 75, resolution = 0.8, 50
  communities, modularity 0.97053.

UPGMA comparisons:

- Arm-level UPGMA: 14 clusters, silhouette 0.342, exact agreement with 12 of
  15 Leiden communities.
- Sequence-level UPGMA: best k = 150 by silhouette 0.63469.

Character bootstrap:

- Existing D-M9 bootstrap resamples rows of a collapsed PHR involvement table
  and rebuilds a small arm-level surrogate distance.
- C2 confirmed it does not read the 15,668 x 15,668 matrix.
- The q-arm 0.1% support is instability under that surrogate, not evidence that
  the full sequence-level q-arm signal is unstable.

### What Does Not Exist

No completed artifact proves formal sampling stability of the Leiden
assignments. A valid D-1 sampling-stability analysis would resample PHRs,
samples, or haplotypes, rebuild the similarity matrix or k-NN graph, rerun the
Leiden selection rule, and report ARI/NMI or co-clustering probabilities against
the baseline partition. C0c created a Slurm-ready wrapper/specification for this
future work, but the implementation and outputs are not present.

Manuscript-facing conclusion:

> These analyses show that the community assignments are not tied to a single
> unexamined resolution or clustering family, but they do not establish formal
> sample-resampling stability.

## Exact Suggested Manuscript Language

### Abstract Replacement

Replace the current abstract sentence containing "identify a linked q-arm group"
with:

```latex
These pseudo-homolog regions define 15 sequence-similarity communities that
recover Xp/Yp, Xq/Yq, acrocentric, 10p--18p and 4q--10q DUX4 systems; one
enriched q-arm neighborhood includes 22q, 21q, 19q, 1q, 13q and 17q.
```

If the abstract needs to be shorter or more conservative, use:

```latex
These pseudo-homolog regions define 15 sequence-similarity communities that
recover Xp/Yp, Xq/Yq, acrocentric, 10p--18p and 4q--10q DUX4 systems, with
additional q-arm similarity visible as a locally enriched neighborhood.
```

### Results: Community Paragraph

Replace the current paragraph at `submission/paper.tex:151-169` with:

```latex
Read directly, the implicit graph of all 15,668 PHRs forms a single connected
component (Fig.~\ref{fig:fig2}A) that is too tangled to interpret node by node.
We therefore measured pangenome similarity between chromosome ends as the
Jaccard overlap of the graph nodes they traverse, and reduced the result to an
arm-level similarity matrix. The matrix is structured rather than uniform: it
contains locally dense blocks embedded in a broader continuum of lower
inter-community similarity (Fig.~\ref{fig:fig2}B, left). Leiden community
detection on this fixed arm-level matrix partitions the 41 signal-bearing arms
into 15 arm-level communities (Fig.~\ref{fig:fig2}B, right). The communities
recover every previously described case of inter-chromosomal subtelomere
homology: the two pseudoautosomal pairs Xp/Yp and Xq/Yq
\cite{sexchrompars_acquaviva2020}, the five acrocentric short arms
\cite{acrocentric_Altemose2022}, the 10p/18p pair first reported by
Linardopoulou and colleagues \cite{Linardopoulou2005}, and the 4q/10q D4Z4
pair \cite{dux4_d4z4_fshd_lemmers2010worldwide}. The same matrix highlights an
enriched q-arm neighborhood, illustrated by 22q, 21q, 19q, 1q, 13q and 17q,
within the broader continuum. Average-linkage UPGMA provides a secondary
display and method comparison on the same matrix, but we do not interpret the
ordering as a phylogeny.
```

Optional second sentence if the revision wants the method-agreement number in
Results:

```latex
At arm level, UPGMA gives 14 clusters and agrees exactly with 12 of the 15
Leiden communities, supporting the major named systems while leaving boundary
assignments as method-dependent.
```

### Results: Remove Clade Wording

Replace the current `submission/paper.tex:257` phrase:

```latex
population-scale extension of the cytogenetically defined clades that anchor it
```

with:

```latex
population-scale extension of the cytogenetically defined subtelomeric systems
that anchor it
```

### Fig. 2 Legend Replacement

Replace the current Fig. 2B legend text around `submission/paper.tex:327-335`
with:

```latex
\textbf{(B)} Arm-level Jaccard similarity heatmap ($41 \times 41$), shown
twice. \emph{Left}, ordered by an average-linkage UPGMA dendrogram used as a
display order for similar subtelomeric PHR profiles (p-arms red, q-arms blue),
showing locally dense blocks on a broader similarity continuum. \emph{Right},
the same matrix ordered by the 15-community Leiden arm-level partition, with
community colour bands. The blocks recover the known systems shown here: PAR1,
PAR2, the acrocentric p-arms, 10p/18p and the 4q/10q DUX4 pair, and show an
enriched q-arm neighborhood exemplified by 22q, 21q, 19q, 1q, 13q and 17q.
```

### Methods: Leiden Community Detection

Replace or extend the current `Community detection (Leiden)` Methods text with:

```latex
Arm-level communities were inferred from the 41-arm Jaccard distance matrix
after converting distances to graph weights and scanning Leiden resolution from
0.1 to 3.0. The selected arm-level operating point was the maximum-silhouette
15-community partition at resolution 1.16; the same community count was obtained
across resolutions 1.13--1.18. Sequence-level communities were inferred from a
k-nearest-neighbor graph on the cached 15,668-sequence distance matrix. The
reported 50-community sequence partition uses the constrained operating point
$k_{\mathrm{NN}} = 75$, resolution 0.8; the unconstrained modularity maximum in
the scan gives a finer 148-community partition. These scans evaluate resolution
dependence on fixed similarity matrices and are not a formal sample-resampling
stability analysis.
```

### Methods: UPGMA

Replace the current terse `UPGMA dendrogram` Methods text with:

```latex
For display and algorithmic comparison, we also ordered the arm-level matrix by
average-linkage hierarchical clustering using \texttt{hclust(..., method =
"average")}. Cutting the same arm-level distance matrix into 14 UPGMA clusters
gave mean silhouette 0.342 and exact agreement with 12 of the 15 Leiden
communities. We use this dendrogram as a heatmap ordering and method comparison,
not as a phylogenetic model.
```

### Methods: NJ and Character Bootstrap

Preferred action: delete the `Neighbour-joining tree and character-level
bootstrap` subsection from the active manuscript. It is not needed for any
current claim and risks confusing the continuum framing.

If author preference is to retain a short caveat, replace the subsection with:

```latex
We also audited an exploratory neighbour-joining and PHR-resampling analysis
used in earlier drafts. That analysis resamples rows of a collapsed per-PHR
cross-chromosome involvement table, rebuilds an arm-level surrogate distance
matrix, and recomputes NJ and UPGMA trees. It does not resample the full
15,668-sequence Jaccard matrix and is therefore not used as evidence for the
Leiden communities or for a closed q-arm class. Under this chromosome-granular
surrogate, the illustrative six-arm q-neighborhood has low support, so we treat
the q-arm pattern as an enriched neighborhood in the similarity matrix rather
than as a bootstrap-supported clade.
```

### Response-Letter Language

Use this if replying directly to reviewer concerns about discreteness,
community methods, or q-arm stability:

```latex
We revised the framing from discrete clades to locally dense neighborhoods on a
continuum. At arm level, within-Leiden pairs are much more similar than
between-community pairs (median 0.431 versus 0.0078), but the off-diagonal
background remains variable rather than zero. At sequence level, the full
15,668-sequence matrix shows the same pattern: a broad low-similarity
background with localized high-similarity peaks. The q-arm set
22q/21q/19q/1q/13q/17q is therefore described as an enriched q-arm
neighborhood, not as a closed clade. We also clarified that Leiden resolution
and UPGMA comparisons assess method and resolution dependence on fixed
similarity matrices; formal PHR/sample-resampling stability remains future
work and is not claimed in the manuscript.
```

## Downstream Edit Checklist

- Replace abstract q-arm language with "enriched q-arm neighborhood" wording.
- Replace Results "discrete blocks" with "locally dense blocks on a broader
  continuum."
- Replace "reveal a previously uncharacterized q-arm grouping" with the
  enriched-neighborhood sentence.
- Replace Fig. 2 "UPGMA tree" / "`phylogeny'" with "average-linkage UPGMA
  dendrogram/order."
- Replace Fig. 2 "tight {...} q-arm grouping" with "enriched q-arm
  neighborhood exemplified by..."
- Replace "clades" at `submission/paper.tex:257` with "subtelomeric systems."
- Keep UPGMA only as display/method comparison; do not use it as sampling
  stability.
- Delete or demote NJ/bootstrap Methods. If retained, call it a collapsed
  per-PHR involvement surrogate and state that it is not a full matrix
  bootstrap.
- Do not claim formal sample-resampling stability for the Leiden communities
  unless a later Slurm task produces ARI/NMI/co-clustering outputs.

## Validation

- Addresses whether the data support a two-tier continuum framing: yes, see
  "Evidence for the Two-Tier Continuum."
- Identifies figures that should back the framing: yes, see "Figure Support."
- Decides q-arm sextet language: retain only as softened enriched-neighborhood
  language, not a closed group.
- States defensible community-method language: yes, see "Defensible
  Community-Method Language."
- States stability/subsampling results and limits: yes, see "Stability and
  Subsampling."
- Provides exact suggested manuscript language: yes, see "Exact Suggested
  Manuscript Language."
- Avoids heavy analysis in this fan-in node: yes; only upstream markdown and
  lightweight TSV outputs were read.
