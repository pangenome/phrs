# A5-A7 Mechanical Audit

Date: 2026-06-17

Scope: `submission/paper.tex`, plus the Fig. 1 source path called out in
`paper_prep/manuscript_revision/00_inventory.md`. This is an audit artifact
only; `paper.tex` was not edited.

## Summary

| Item | Status | Minimal downstream action |
| --- | --- | --- |
| A-5 undefined `XTR` acronym | Absent from active manuscript text and Fig. 1 caption; still present in Fig. 1 plotting source and original revision prompt | No manuscript edit needed unless regenerated Fig. 1/legend exposes `XTR`; if exposed, expand as `X-transposed region (XTR)` on first use or remove from the visible legend |
| A-6 empty required sections | Present | Replace TODO-only Author Contributions and Competing Interests stubs with author-approved text |
| A-7 symbol collisions | Present in Methods | Rename nearest-neighbor graph parameter to `k_{\mathrm{NN}}`; describe community counts as "15 arm communities", "50 sequence communities", and "14 UPGMA clusters" instead of reusing bare `k` |

## Line-Numbered Findings

### A-5. Undefined `XTR` Acronym

1. `submission/paper.tex` contains no literal `XTR` occurrence. Command used:
   `rg -n "\bXTR\b" submission/paper.tex`.
2. The active Fig. 1 caption already avoids `XTR`: lines
   `submission/paper.tex:307-312` describe the region track as
   "centromere, PAR, PHR, terminal repeat". No acronym expansion is needed in
   the current caption.
3. The source-side risk is still real. The revision prompt records the original
   issue at `wg_manuscript_revision_prompt.md:127`: Fig. 1 listed `XTR`
   meaning "X-transposed region". The upstream plotting script still collapses
   `XTR1` and `XTR2` to a visible label `XTR` at
   `scripts/plot-impg-coverage.R:41`.

Minimal fix if the visible manuscript figure or legend is regenerated with
`XTR`: change the first visible occurrence to `X-transposed region (XTR)` or
change the visible legend label to `X-transposed region`. If the current
`submission/paper.tex:311` caption remains authoritative and the figure image
does not show `XTR`, record A-5 as already absent from the active manuscript.

### A-6. Empty Required Sections

1. `submission/paper.tex:282` starts `\subsection*{Author Contributions}`.
   Lines `submission/paper.tex:284-286` contain only a TODO comment:
   `TODO: confirm CRediT roles with all authors...`; no printable author
   contribution statement follows before the next subsection.
2. `submission/paper.tex:288` starts `\subsection*{Competing Interests}`.
   Line `submission/paper.tex:290` contains only `% TODO: fill in`; no
   printable competing-interest statement follows before `\clearpage`.
3. `submission/README.md:244-247` independently confirms that Acknowledgments
   are filled but Author Contributions and Competing Interests remain TODO
   stubs.

Exact minimal fixes:

```tex
\subsection*{Author Contributions}

Andrea Guarracino and Erik Garrison conceived the study. Andrea Guarracino,
Angela Gyamfi and Erik Garrison performed analyses, interpreted the results and
wrote the manuscript.

\subsection*{Competing Interests}

The authors declare no competing interests.
```

Author Contributions should be treated as author-confirmation-sensitive because
the manuscript source itself says to confirm CRediT roles. If the fan-in worker
is applying only non-judgment mechanical fixes, the safer exact replacement for
the first section is:

```tex
\subsection*{Author Contributions}

Author contributions will be finalized before submission.
```

That avoids a blank required section but is not submission-ready. The Competing
Interests sentence above is the standard minimal declaration if true.

### A-7. Reused or Conflicting Symbols in Methods

1. `submission/paper.tex:501-503` uses bare `$k = 15$` for the selected
   arm-level Leiden partition size.
2. `submission/paper.tex:504-506` reuses bare `k` twice with different meanings
   in one sentence: first as the sequence-level k-nearest-neighbor graph degree
   (`selected $k = 75$`), then as the resulting number of sequence-level
   communities (`$k = 50$`). This is the main collision.
3. `submission/paper.tex:510-511` uses bare `$k = 14$` for the UPGMA dendrogram
   cut, adjacent to the Leiden community-detection method. This is less severe
   because it is in a separate subsection, but it reinforces the same ambiguity.
4. `submission/paper.tex:338` uses `Leiden $k = 15$ community` in the Fig. 2
   legend; it is understandable but would be clearer as a community count.
5. Other symbols are not collisions needing mechanical edits:
   `submission/paper.tex:475-477` uses `p^{*}` and `n` for the Erdos-Renyi
   connectivity threshold; `submission/paper.tex:553-565` uses `B/W`, `W/B`,
   `n`, and `p` as conventional statistical quantities; and
   `submission/paper.tex:628-631` uses `$B = 10000$` for Monte Carlo
   permutations. These are context-specific and do not conflict with the
   Leiden `k` issue.

Exact minimal replacement for `submission/paper.tex:501-506`:

```tex
Arm-level: edge weight $w_{ij} = \exp(-d_{ij} / \text{median}(d))$; resolution
scanned 0.1--3.0 at 0.01 step; selected resolution 1.16 (mean silhouette
0.347), yielding 15 arm-level communities.
Sequence-level: k-NN graph with $k_{\mathrm{NN}} \in \{10, 25, 50, 75, 100, 125\}$ and
resolution 0.1--3.0; selected $k_{\mathrm{NN}} = 75$, resolution 0.8
(modularity 0.97, mean silhouette 0.602), yielding 50 sequence-level
communities.
```

Exact minimal replacement for `submission/paper.tex:510-511`:

```tex
\texttt{hclust(..., method = "average")}; tree cut into 14 clusters
(mean silhouette 0.342); agreement with Leiden 12 of 15.
```

Exact minimal replacement for `submission/paper.tex:338-339`:

```tex
\emph{Right}, the same matrix ordered by the 15-community Leiden arm-level
partition, with community colour bands.
```

These edits preserve all numerical results while assigning each number one
plain meaning: nearest-neighbor graph degree (`k_{\mathrm{NN}}`), arm-level
community count (15), sequence-level community count (50), and UPGMA cluster
count (14).

## Validation Checklist

- Artifact exists: `paper_prep/manuscript_revision/A5_A7_mechanical_audit.md`.
- A-5 records that `XTR` is absent from `submission/paper.tex` and identifies
  the remaining source-side occurrence at `scripts/plot-impg-coverage.R:41`.
- A-6 reports the Author Contributions and Competing Interests status with
  line numbers.
- A-7 lists symbol collisions and proposes exact disambiguations.
