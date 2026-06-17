# C1 Tree Essentiality Decision Record

Date: 2026-06-17
Task: `manuscript-revision-c1`

## Question

Determine whether any current manuscript claim still depends on the NJ/UPGMA
tree or character bootstrap after the continuum reframing. This record audits
the active manuscript, active Fig. 2 assets, figure-package captions, and the
tree/bootstrap scripts. It does not edit the manuscript or delete any J-task
material.

## Files Audited

- Active manuscript: `submission/paper.tex`
- Active Fig. 2 manuscript assets and provenance:
  `submission/fig/MainFigures/Fig2bc_jaccard.{pdf,png}`,
  `submission/fig/MainFigures/arm_order_tree.tsv`,
  `submission/fig/MainFigures/arm_order_community.tsv`,
  `submission/fig/MainFigures/source_audit.tsv`, and
  `submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R`
- Figure-package captions and script:
  `paper_prep/figures/fig1/caption.md`,
  `paper_prep/figures/fig2/caption.md`,
  `paper_prep/figures/fig2/sources.tsv`,
  `paper_prep/figures/fig2/figure_fig2.py`
- NJ tree package: `paper_prep/figures/nj_tree_arms/README.md`,
  `paper_prep/figures/nj_tree_arms/nj_tree.R`,
  `paper_prep/figures/nj_tree_arms/nj_tree.newick`, and
  `paper_prep/figures/nj_tree_arms/nj_tree_annotated.{pdf,png}`
- Character-bootstrap audit and script:
  `scripts/cladistics/char_bootstrap_d_m9.R` and
  `paper_prep/manuscript_revision/C2_bootstrap_audit.md`
- Continuum dependency context:
  `paper_prep/manuscript_revision/C0_continuum/C0a_arm_level_report.md`

## Bottom Line

No current biological claim has to depend on the NJ tree or character bootstrap.
The essential claim can be supported by the arm-level Jaccard matrix, the Leiden
k = 15 partition, and direct heatmap/block evidence. UPGMA remains useful only
as a row/column ordering or a secondary agreement check. The NJ tree is not used
by the active Fig. 2 asset and is not needed for the manuscript's continuum
framing. The character bootstrap is not load-bearing and is potentially
confusing because C2 found it operates on a collapsed per-PHR involvement
surrogate, not the full 15,668 x 15,668 sequence-level Jaccard matrix.

Named author decision required: Erik Garrison must decide whether the revision
should (1) delete the NJ/bootstrap Methods subsection entirely, (2) demote it to
a clearly labeled exploratory/sensitivity note, or (3) keep it only if relabeled
as a non-phylogenetic ordering/surrogate analysis with the C2 limitation stated.

## Active Manuscript Mentions

| File:line | Mention | Claim supported now | Load-bearing for current claims? | Decision options |
|---|---|---|---|---|
| `submission/paper.tex:55` | Abstract says the data "support a model" of nuclear architecture and recombination opportunity. | Overall model claim from sequence sharing, 3D contact, mouse meiosis, and pedigree evidence. | No for tree/bootstrap. This is a generic support verb, not statistical branch support. | Keep as-is unless G/E tasks change model strength; not a C1 tree issue. |
| `submission/paper.tex:162` to `submission/paper.tex:165` | Matrix "resolves into discrete blocks" and Leiden partitions 41 arms into 15 communities. | Main community claim. | No for NJ. The claim is supported directly by the 41 x 41 matrix and Leiden partition. UPGMA ordering can help show blocks but is not the proof. | Keep community claim; optionally add continuum wording from C0. |
| `submission/paper.tex:166` to `submission/paper.tex:172` | Communities recover known cases and identify the q-arm grouping. | Named-system recovery and q-arm grouping claim. | No for NJ/bootstrap. The manuscript names these as Leiden communities and heatmap blocks, not NJ clades. The q-arm grouping has separate C0/C3 decision records. | Keep if author accepts q-arm wording; avoid using bootstrap as validation of the q-arm list. |
| `submission/paper.tex:173` to `submission/paper.tex:175` | "The partition is robust to method, with hierarchical clustering agreeing on 12 of the 15 communities, and we use the ordering only as a grouping device rather than as a phylogenetic claim." | Secondary robustness/method-agreement claim; disclaims phylogeny. | Weakly load-bearing only if the manuscript wants a method-robustness sentence. Not essential to the core result. It may conflict with Fig. 1 package caption claiming 14/15. | Delete, or demote to Methods/figure legend; if kept, relabel as "UPGMA average-linkage ordering agrees with 12/15 Leiden communities" and reconcile 12/15 vs 14/15. |
| `submission/paper.tex:335` to `submission/paper.tex:337` | Fig. 2 caption: left heatmap ordered by the "UPGMA tree" and called the "`phylogeny' of subtelomeric PHRs"; it resolves blocks rather than a uniform gradient. | Visual ordering of the left 41 x 41 heatmap. | No for biological inference. The ordering is a display choice. Calling it "phylogeny" makes it sound load-bearing despite the Methods disclaimer. | Preferred: relabel to "UPGMA average-linkage dendrogram/order" and remove "`phylogeny'"; or delete the tree label entirely and call it "hierarchical-clustering order." |
| `submission/paper.tex:338` to `submission/paper.tex:342` | Fig. 2 caption: community-ordered matrix with bands; blocks recover known cases and q-arm grouping. | Main figure support for communities. | No for NJ/bootstrap. This is the cleaner load-bearing Fig. 2 statement. | Keep, subject to C0/C3 q-arm wording. |
| `submission/paper.tex:508` to `submission/paper.tex:511` | Methods subsection "UPGMA dendrogram"; average linkage, k = 14, mean silhouette 0.342, agreement with Leiden 12/15. | Documents hierarchical clustering check and ordering source. | Not essential. Useful provenance if the left heatmap remains tree-ordered. | Demote/relabel as "Average-linkage heatmap ordering"; keep only if Fig. 2 retains UPGMA-ordered panel. Reconcile 12/15 vs other artifacts before retaining as a robustness metric. |
| `submission/paper.tex:513` to `submission/paper.tex:527` | Methods subsection "Neighbour-joining tree and character-level bootstrap"; `ape::nj()` on 41 x 41 matrix; 1,000-replicate character bootstrap; NJ + UPGMA support values; q-arm support 0.1%; NJ used as ordering, not phylogenetic. | Legacy tree/bootstrap support layer for named groupings. | No. The NJ tree is not the active Fig. 2 ordering, and the bootstrap does not support the continuum-framed community claims. C2 shows the character bootstrap uses a collapsed involvement surrogate, not the full sequence-level matrix. | Preferred: delete from active manuscript unless an author explicitly wants a caveated sensitivity note. If retained, relabel as "surrogate PHR-involvement bootstrap" and state it is not a full sequence-matrix bootstrap. |
| `submission/paper.tex:553` | "W/B ratio computed by bootstrap on 10,000 permutations" in 3D-contact Methods. | Contact-map within/between significance. | No for C1. This is a permutation/bootstrap phrase unrelated to NJ/UPGMA tree support. | Keep or handle under 3D/statistics tasks; do not delete as tree material. |

## Active Fig. 2 Asset Audit

The current manuscript Fig. 2 panel B asset is produced by
`submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R`, not by the
`paper_prep/figures/fig2/figure_fig2.py` package.

- `submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R:21` to
  `submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R:24` read the
  vendored arm-level distance matrix and Leiden assignments.
- `submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R:114` to
  `submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R:116` compute
  `hclust(as.dist(D), method = "average")` and use its labels as `tree_order`.
- `submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R:120` to
  `submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R:124` compute the
  community-ordered view by Leiden community, with UPGMA leaf position only used
  to order arms within each community.
- `submission/fig/MainFigures/source_audit.tsv:3` states explicitly that the
  tree-order source is UPGMA average linkage, "not an NJ tree."
- `submission/fig/MainFigures/source_audit.tsv:9` states the community order is
  C1-C15, then UPGMA position within community.

Implication: the active Fig. 2 heatmap does not require the NJ tree at all.
UPGMA is currently a layout choice for the left panel and a within-community
tie-breaker for the right panel. The manuscript should not describe the left
panel as a phylogeny unless an author wants to make a stronger, explicitly
defended phylogenetic claim.

## Figure-Package Caption and Script Mentions

These are not the active `submission/paper.tex` figure legends, but they are
relevant upstream material and may be copied into revision text if not audited.

| File:line | Mention | Claim supported | Load-bearing? | Decision options |
|---|---|---|---|---|
| `paper_prep/figures/fig1/caption.md:12` to `paper_prep/figures/fig1/caption.md:17` | 41 x 41 heatmap ordered by Leiden; UPGMA k = 14 dendrogram recovers 14/15 Leiden communities. | Older figure-package robustness statement. | Not load-bearing in active manuscript. Also conflicts with active text/Methods saying 12/15. | Do not copy blindly. If reused, reconcile the agreement statistic and label as secondary. |
| `paper_prep/figures/fig2/caption.md:52` to `paper_prep/figures/fig2/caption.md:54` | "Out-of-Africa tree" as UPGMA dendrogram from cross-arm Fst matrix. | Population-history panel in older Fig. 2 package. | Not load-bearing for C1. It is an Fst/population dendrogram, not the subtelomeric PHR NJ/UPGMA tree. | Leave for popgen/figure decisions; do not treat as PHR tree support. |
| `paper_prep/figures/fig2/sources.tsv:10` | Source manifest for "Out-of-Africa UPGMA dendrogram from Fst." | Same as above. | Not C1 load-bearing. | Leave unless Fig. 2 package is revived. |
| `paper_prep/figures/fig2/figure_fig2.py:3` to `paper_prep/figures/fig2/figure_fig2.py:9` and `:299` to `:412` | Implements UPGMA dendrogram for Fst matrix. | Older population-history panel. | Not C1 load-bearing. | Leave; no manuscript edit. |

## NJ Tree Package Mentions

| File:line | Mention | Claim supported | Load-bearing? | Decision options |
|---|---|---|---|---|
| `paper_prep/figures/nj_tree_arms/README.md:3` to `:12` | Builds an unrooted NJ tree from the 41 x 41 Jaccard matrix; named clades correspond to Leiden communities. | Legacy visual summary of named systems as NJ monophyletic clades. | Not in active manuscript figures. Not essential after continuum reframing. | Keep as archived supporting material, or demote to supplementary/exploratory if retained. |
| `paper_prep/figures/nj_tree_arms/README.md:13` to `:20` | Perturbation bootstrap gives 100% support for named clades; true character bootstrap not possible from derived distance summary. | Distance-perturbation stability of named NJ clades. | Not load-bearing and superseded/complicated by D-M9/C2 character-bootstrap audit. | Do not cite in active text without caveat; if retained, call it distance-perturbation sensitivity, not character support. |
| `paper_prep/figures/nj_tree_arms/nj_tree.R:66` to `:69` | `ape::nj()` on the arm-level distance matrix. | Script implementation of NJ tree. | Not load-bearing for active Fig. 2. | Archive or supplementary only. |
| `paper_prep/figures/nj_tree_arms/nj_tree.R:102` to `:167` | Perturbation bootstrap loop. | Sensitivity of NJ edges to distance noise. | Not load-bearing for manuscript claims. | Relabel if used; do not merge with character-bootstrap language. |
| `paper_prep/figures/nj_tree_arms/nj_tree.R:221` to `:228` | Plot title and subtitle print "NJ tree" and bootstrap support at nodes. | Figure labeling. | Not active manuscript material. | If figure is ever reused, change title/caption to avoid implying phylogeny. |
| `paper_prep/figures/nj_tree_arms/nj_tree.R:253` to `:285` | Sidecar summary and per-clade bootstrap support. | README support table. | Not load-bearing. | Keep as provenance; do not use as decisive evidence. |

## Character Bootstrap Mentions

`scripts/cladistics/char_bootstrap_d_m9.R` is relevant only if the active
manuscript keeps the Methods subsection at `submission/paper.tex:513` to
`submission/paper.tex:527`.

C2's audit is the controlling interpretation:

- `paper_prep/manuscript_revision/C2_bootstrap_audit.md:33` to `:38` confirms
  the script describes PHR-row resampling, recomputation of an arm-level Jaccard
  matrix, and NJ plus UPGMA trees.
- `paper_prep/manuscript_revision/C2_bootstrap_audit.md:184` to `:206` finds
  that the script never reads the 15,668 x 15,668 sequence-level matrix.
- `paper_prep/manuscript_revision/C2_bootstrap_audit.md:220` to `:249` states
  that the 0.1% q-arm support is instability of the audited chromosome-level
  involvement surrogate, not evidence that the six-q-arm grouping is unstable in
  the full sequence-level similarity matrix.

Implication: the active manuscript's "character-level bootstrap" wording is
too strong unless it is explicitly scoped to the collapsed involvement cache.
It should not be used to support or reject the q-arm community in the continuum
framing.

## Load-Bearing Assessment By Claim

| Claim | Current better support | Tree/bootstrap essential? | Notes |
|---|---|---|---|
| 41 signal-bearing arms form 15 sequence-similarity communities. | Leiden k = 15 on the arm-level Jaccard matrix; community-ordered heatmap. | No. | UPGMA can remain as a secondary display order. |
| Known systems are recovered: PAR1, PAR2, acrocentric p-arms, 10p/18p, 4q/10q. | Community membership and direct heatmap blocks. | No. | NJ monophyly is redundant. |
| Six q-arm grouping exists or can be named. | C0 arm-level density, C3 wording decision, Leiden community, heatmap block. | No. | Character bootstrap should not adjudicate this claim. |
| The matrix is not a uniform gradient but has stronger within-community structure on a continuum. | C0a arm-level continuum report and heatmap/distribution summaries. | No. | Tree language may obscure continuum framing. |
| Ordering in Fig. 2 left panel. | UPGMA average-linkage row/column order. | Yes only as layout provenance. | This is not biological support; call it an ordering. |
| Robustness to method. | UPGMA/Leiden agreement metric, if reconciled. | Weak and optional. | Current 12/15 vs 14/15 discrepancy needs resolution before using. |

## Decision Options

### Option 1: Delete NJ/bootstrap from the active manuscript

Recommended if the revision goal is a clean continuum/community story.

- Delete or do not carry forward the Methods subsection at
  `submission/paper.tex:513` to `submission/paper.tex:527`.
- Relabel Fig. 2 left panel from "UPGMA tree" / "`phylogeny'" to "UPGMA
  average-linkage order" or "hierarchical-clustering order."
- Keep the UPGMA dendrogram Methods subsection only if needed to document the
  heatmap order.

### Option 2: Demote tree/bootstrap material

Recommended if an author wants to preserve the work without making it central.

- Move NJ/bootstrap to a supplementary note or revision-support note.
- State that NJ is an exploratory view of the arm-level distance matrix.
- State that the D-M9 bootstrap is a surrogate PHR-involvement bootstrap and
  does not resample the full 15,668 x 15,668 matrix.

### Option 3: Relabel and keep in Methods

Use only if Erik Garrison decides the manuscript needs the tree material in the
main Methods.

- Replace "phylogeny" with "non-phylogenetic ordering" everywhere.
- Replace "character-level bootstrap" with "collapsed per-PHR involvement
  surrogate bootstrap" or an equivalent exact phrase.
- Keep the q-arm 0.1% result only with C2's limitation: it is low support under
  the surrogate, not a full-matrix stability test.
- Reconcile all agreement numbers before publication: active manuscript says
  12/15, while `paper_prep/figures/fig1/caption.md` says 14/15.

## Author Decision Point

Erik Garrison must choose one of the three options above before any downstream
fan-in edits `submission/paper.tex` or figure legends. My technical
recommendation is Option 1 for the active manuscript, with Option 2 as the
fallback if authors want to preserve the NJ/bootstrap analyses as provenance.
Option 3 is defensible only with explicit caveats and reconciled agreement
statistics.

## Downstream Consumers

- `manuscript-revision-cd-fanin` should consume this record together with C0,
  C2, and C3 before drafting continuum/community manuscript edits.
- No J material was deleted, renamed, or edited by this task.
