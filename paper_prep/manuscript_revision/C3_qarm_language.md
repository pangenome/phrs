# C3 q-arm Language Audit

Date: 2026-06-17

## Scope

This audit covers the active manuscript at `submission/paper.tex`. It finds the
q-arm sextet/list and nearby language that could make the C6/q-arm pattern read
as a bounded discovery, closed clade, or isolated group. It does not edit the
manuscript.

The current C0 evidence supports a cautious framing:

- C0a arm-level report: the C6/q-arm block is denser than arm-level background,
  but the tests are descriptive because C6 was defined from the same matrix
  (`paper_prep/manuscript_revision/C0_continuum/C0a_arm_level_report.md:51-53`).
- C0b sequence-level report: the sequence-level distribution is a continuum with
  localized high-similarity peaks; C6/q-arm pairs are enriched at fixed
  thresholds but should be treated as an enriched neighborhood, not an isolated
  bounded class
  (`paper_prep/manuscript_revision/C0_continuum/C0b_sequence_level_report.md:55-57`).

## Direct q-arm Sextet/List Occurrences

| Location | Current language | Label | Risk | Proposed dense-neighborhood wording | Proposed artifact/insufficient wording |
|---|---|---|---|---|---|
| `submission/paper.tex:55` | Abstract: "identify a linked q-arm group involving 22q, 21q, 19q, 1q, 13q and 17q." | Overclaim | "identify" + "linked ... group" in the abstract reads as a new closed entity. The list is presented as a result on the same footing as known systems. | "These pseudo-homolog regions define 15 sequence-similarity communities that recover Xp/Yp, Xq/Yq, acrocentric, 10p--18p and 4q--10q DUX4 systems; one enriched q-arm neighborhood includes 22q, 21q, 19q, 1q, 13q and 17q." | Remove the sextet from the abstract: "These pseudo-homolog regions define 15 sequence-similarity communities that recover Xp/Yp, Xq/Yq, acrocentric, 10p--18p and 4q--10q DUX4 systems, while additional q-arm structure remains under evaluation." |
| `submission/paper.tex:171-172` | Results: "they also reveal a previously uncharacterized q-arm grouping of 22q, 21q, 19q, 1q, 13q and 17q." | Overclaim | "reveal" and "previously uncharacterized ... grouping" imply a discrete discovery. The line immediately follows known-system recovery, so it can be read as equal-certainty classification. | "The same ordering highlights an enriched q-arm neighborhood, illustrated by 22q, 21q, 19q, 1q, 13q and 17q, within the broader arm-level continuum." | "The same ordering shows additional q-arm similarity, but we do not name a bounded q-arm group from the present heatmap." |
| `submission/paper.tex:341-342` | Fig. 2 legend: "a tight \{22q, 21q, 19q, 1q, 13q, 17q\} q-arm grouping." | Overclaim | "tight" + braces strongly signals a closed set. In a stand-alone figure legend, readers may treat the brace list as the discovered class boundary. | "and an enriched q-arm neighborhood exemplified by 22q, 21q, 19q, 1q, 13q and 17q." | "and additional q-arm similarity visible in the community-ordered heatmap." |
| `submission/paper.tex:523-526` | Methods: "the 6-arm tight q-arm grouping is 0.1% because..." | Overclaim, but useful caveat source | The phrase still names a "tight" 6-arm grouping, even while explaining the surrogate-bootstrap limitation. This creates mixed messaging: a tight group with essentially no support under the described bootstrap. | "the illustrative six-arm q-neighborhood has 0.1% support under this chromosome-granular surrogate because..." | "q-arm support from this surrogate is not interpretable because the per-PHR involvement column is chromosome-granular..." |

## Bounded/Closed-Group Language Nearby

| Location | Current language | Label | Risk | Proposed dense-neighborhood wording | Proposed artifact/insufficient wording |
|---|---|---|---|---|---|
| `submission/paper.tex:162-164` | "The matrix is not a uniform gradient but resolves into discrete blocks... partitions the 41 signal-bearing arms into 15 arm-level communities." | Acceptable with revision | It correctly distinguishes matrix structure from a uniform gradient, but "resolves into discrete blocks" can overstate hard boundaries after C0's continuum finding. | "The matrix is structured rather than uniform, with locally dense blocks on a broader continuum... partitions..." | "The matrix shows non-uniform similarity structure... partitions..." |
| `submission/paper.tex:166-172` | "The communities recover every previously described case... they also reveal..." | Mixed: known systems acceptable; q-arm clause overclaim | "Recover" is appropriate for known systems. "Every" is strong but tied to named cases. The q-arm clause should not share the same closed-discovery cadence. | Keep known-system recovery; replace the q-arm clause with "and provide an illustrative q-arm neighborhood..." | Keep known-system recovery; remove q-arm clause or demote to "additional q-arm similarity is visible but not named here." |
| `submission/paper.tex:173-175` | "The partition is robust to method... ordering only as a grouping device rather than as a phylogenetic claim." | Acceptable caveat, but adjacent to overclaim | The caveat helps. It does not by itself solve the q-arm closed-list issue because the preceding sentence has already named the sextet as a grouping. | "The UPGMA ordering agrees with 12 of the 15 Leiden communities, and we use that ordering only as a display device, not as phylogeny or a boundary test." | Same as dense-neighborhood wording; if q-arm list is removed, this remains useful as a general caveat. |
| `submission/paper.tex:263` | "population-scale extension of the cytogenetically defined clades that anchor it" | Overclaim outside q-arm paragraph | "clades" reintroduces phylogenetic/closed-group language in the synthesis. It can retroactively strengthen the q-arm list as a clade-like discovery. | "population-scale extension of the cytogenetically defined subtelomeric systems that anchor it" | Same replacement. This is not conditional on the q-arm outcome. |
| `submission/paper.tex:335-337` | Fig. 2 legend: "UPGMA tree (the ``phylogeny'' of subtelomeric PHRs...)... resolves the matrix into discrete blocks..." | Overclaim | Scare quotes do not sufficiently neutralize "phylogeny" in a stand-alone legend. "Discrete blocks" is stronger than the C0 continuum result. | "ordered by a UPGMA dendrogram used only to group similar subtelomeric PHR profiles... shows locally dense blocks on a broader similarity continuum." | "ordered by a UPGMA dendrogram used only to display similarity structure... shows non-uniform arm-level similarity." |
| `submission/paper.tex:339-342` | Fig. 2 legend: "The blocks recover every known case... and a tight ... q-arm grouping." | Mixed: known systems acceptable; q-arm clause overclaim | "Blocks recover every known case" is acceptable if limited to known systems. The q-arm brace list is the risky part. | "The blocks recover the known systems shown here... and show an enriched q-arm neighborhood exemplified by..." | "The blocks recover the known systems shown here... with additional q-arm similarity left unnamed." |
| `submission/paper.tex:513-527` | Methods subsection "Neighbour-joining tree and character-level bootstrap"; NJ rooted at acrocentric grouping; q-arm support explanation; "NJ used as an ordering device, not a phylogenetic claim." | Mixed; mostly acceptable caveat, but q-arm "tight grouping" is overclaim | The final sentence is the right caveat. The subsection still gives tree/bootstrap machinery enough prominence that the q-arm list can seem like a clade whose only problem is a technical bootstrap artifact. | Rename prose around C6 to "illustrative q-neighborhood" and state the bootstrap is not an adjudication of a closed q-arm class. | If C0 decision rejects naming, keep only the general NJ/UPGMA ordering provenance and remove named q-arm support language. |

## Recommended House Wording

Use these terms consistently if the C0 heatmap-density decision accepts C6 as a
real local density:

- "enriched q-arm neighborhood"
- "illustrated by 22q, 21q, 19q, 1q, 13q and 17q"
- "locally dense region within a broader similarity continuum"
- "community-ordered heatmap highlights..." rather than "reveals..."
- "grouping" only when preceded by "illustrative" or "display/order-based"

Avoid these terms for C6/q-arm unless later evidence changes the decision:

- "closed group"
- "bounded class"
- "clade", "monophyletic", "cladistic recovery"
- "tight \{...\} q-arm grouping"
- "identify/reveal a linked q-arm group" in abstract-level prose
- "the q-arm sextet" as a biological noun rather than an analysis shorthand

## Conditional Replacement Blocks

### If C0 decision is dense-neighborhood

Abstract replacement:

> These pseudo-homolog regions define 15 sequence-similarity communities that
> recover Xp/Yp, Xq/Yq, acrocentric, 10p--18p and 4q--10q DUX4 systems; one
> enriched q-arm neighborhood includes 22q, 21q, 19q, 1q, 13q and 17q.

Results replacement:

> The communities recover every previously described case of inter-chromosomal
> subtelomere homology: the two pseudoautosomal pairs Xp/Yp and Xq/Yq, the five
> acrocentric short arms, the 10p/18p pair, and the 4q/10q D4Z4 pair. The same
> matrix highlights an enriched q-arm neighborhood, illustrated by 22q, 21q,
> 19q, 1q, 13q and 17q, within a broader continuum of lower inter-community
> similarity.

Fig. 2 legend replacement:

> Right, the same matrix ordered by Leiden k = 15 community, with community
> colour bands. The blocks recover the known systems shown here: PAR1, PAR2,
> the acrocentric p-arms, 10p/18p and the 4q/10q DUX4 pair, and show an
> enriched q-arm neighborhood exemplified by 22q, 21q, 19q, 1q, 13q and 17q.

Methods replacement:

> The illustrative six-arm q-neighborhood has 0.1% support under this
> chromosome-granular surrogate because the per-PHR involvement column
> conflates q-neighborhood signal with chr2_q and chr10_q; this bootstrap is
> therefore not a test of a closed q-arm class.

### If C0 decision is artifact or insufficient support

Abstract replacement:

> These pseudo-homolog regions define 15 sequence-similarity communities that
> recover Xp/Yp, Xq/Yq, acrocentric, 10p--18p and 4q--10q DUX4 systems.

Results replacement:

> The communities recover every previously described case of inter-chromosomal
> subtelomere homology: the two pseudoautosomal pairs Xp/Yp and Xq/Yq, the five
> acrocentric short arms, the 10p/18p pair, and the 4q/10q D4Z4 pair. Additional
> q-arm similarity is visible in the matrix but is not named as a bounded group
> here.

Fig. 2 legend replacement:

> Right, the same matrix ordered by Leiden k = 15 community, with community
> colour bands. The blocks recover the known systems shown here: PAR1, PAR2,
> the acrocentric p-arms, 10p/18p and the 4q/10q DUX4 pair.

Methods replacement:

> q-arm support from this surrogate is not interpretable because the per-PHR
> involvement column is chromosome-granular; a closed q-arm class is not tested
> by this bootstrap.

## Bottom Line

The manuscript should not present 22q/21q/19q/1q/13q/17q as a discovered
closed clade or bounded q-arm group. If C0 is accepted as a heatmap-density
decision, the list can remain as an illustration of an enriched C6
neighborhood. If C0 is judged artifact-prone or insufficient, the list should be
removed from the abstract and figure legend and left only, if needed, as
internal analysis shorthand in Methods or revision notes.
