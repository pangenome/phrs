# Voice review 05: communities and representative architecture

Scope reviewed: `submission/paper.tex:153` to `submission/paper.tex:192` and
figure legends `submission/paper.tex:311` to `submission/paper.tex:354`.
House-voice basis: `paper_prep/synthesis/ABSTRACT_TEXTURE_SYNTHESIS.md`,
`paper_prep/synthesis/ABSTRACT_nature.md`, and the graph/community talk segment
at `submission/notes/Session7-PopulationGenomics.en.srt:767` to
`submission/notes/Session7-PopulationGenomics.en.srt:1045`.

## 1. Community claim spine

The section mostly matches the house abstract voice. The claim spine is
descriptive, evidence-led, and close to the target abstract sentence that PHRs
"define 15 sequence-similarity communities" which recover known systems and
identify a linked q-arm group.

Recommended spine to preserve:

1. `submission/paper.tex:155` to `submission/paper.tex:159`: the all-PHR graph
   is connected and too tangled to interpret directly, so the analysis moves to
   arm-level Jaccard similarity. This follows the talk's "hairball" logic at
   `submission/notes/Session7-PopulationGenomics.en.srt:767` to
   `submission/notes/Session7-PopulationGenomics.en.srt:843`.
2. `submission/paper.tex:160` to `submission/paper.tex:163`: the similarity
   matrix resolves into blocks and Leiden partitions 41 signal-bearing arms into
   15 arm-level communities. This is the central descriptive community claim.
3. `submission/paper.tex:164` to `submission/paper.tex:170`: known systems
   validate the grouping, while the 22q/21q/19q/1q/13q/17q q-arm group carries
   the discovery. This exactly matches the house rule to pair known positives
   with new discoveries.
4. `submission/paper.tex:171` to `submission/paper.tex:173`: the robustness
   sentence correctly prevents over-reading the ordering as phylogeny. This
   qualification should remain adjacent to any mention of tree-like ordering.
5. `submission/paper.tex:177` to `submission/paper.tex:192`: representative
   architectures are framed as sequence-level examples of communities C1, C2 and
   C11, then generalized back to pseudogenes, non-coding RNAs, olfactory
   receptors, subtelomeric pseudogene duplicons and shared orientation. This
   keeps the loci connected to the global result rather than reading as an
   unrelated locus catalog.

Minor spine risk: `submission/paper.tex:177` says "each community carries a
characteristic duplicon" immediately before discussing only three communities.
As written, this sounds like a claim about all 15 communities, but the evidence
presented locally is representative. A safer house-voice version would be
"At the sequence level, representative communities are organized around
characteristic duplicons" or "The representative communities shown in Fig. 3
are organized around characteristic duplicons."

## 2. Term-standardization list

- **community / communities**: Use as the default term for Leiden/Jaccard
  sequence-similarity groupings. Preferred forms: "sequence-similarity
  communities", "arm-level communities", "Leiden communities" when the algorithm
  matters, and "within-community" for comparisons. The body at
  `submission/paper.tex:161` to `submission/paper.tex:163` is good, but Fig. 2
  legend `submission/paper.tex:326` to `submission/paper.tex:327` should clarify
  "Leiden sequence-similarity community" if edited later.
- **block / grouping**: Good descriptive alternatives for visual structure in
  heatmaps. `submission/paper.tex:160` to `submission/paper.tex:163` and
  `submission/paper.tex:327` to `submission/paper.tex:330` use these well. Keep
  "tight q-arm grouping" for 22q/21q/19q/1q/13q/17q.
- **tree**: Acceptable only for the UPGMA ordering object, not for biological
  descent. Use "UPGMA ordering tree", "UPGMA dendrogram" or "tree used for
  ordering" rather than bare "tree" when possible. The body avoids the term;
  Fig. 2 legend `submission/paper.tex:323` needs qualification.
- **phylogeny / phylogenetic**: Use only in negative or scare-quoted form if it
  remains at all. The body sentence `submission/paper.tex:171` to
  `submission/paper.tex:173` is appropriately qualified: "ordering only as a
  grouping device rather than as a phylogenetic claim." The Fig. 2 legend phrase
  "the ``phylogeny'' of subtelomeric PHRs" at `submission/paper.tex:323` to
  `submission/paper.tex:324` is the main unsupported implication because the
  qualifier arrives only by punctuation, not by explanatory prose.
- **clade / cladal**: Avoid in the scoped text. No "clade" term appears in
  `submission/paper.tex:153` to `submission/paper.tex:192` or Fig. 2/Fig. 3
  legends, which is good. The talk used "cladal structure" immediately after a
  phylogeny caveat (`submission/notes/Session7-PopulationGenomics.en.srt:849`
  to `submission/notes/Session7-PopulationGenomics.en.srt:875`), but manuscript
  prose should prefer "block structure" or "community structure."
- **recover / identify / reveal**: "Recover" is correct for known systems at
  `submission/paper.tex:164` and `submission/paper.tex:327`. For the novel q-arm
  pattern, prefer "identify" or "find" over "reveal" if edited for abstract
  consistency; the current body's "reveal" at `submission/paper.tex:169` is
  acceptable but slightly more promotional than the house voice.
- **characteristic duplicon**: Good as a linking term for Fig. 3, but avoid
  implying all 15 communities have been shown in the figure. Qualify with
  "representative" unless the sentence is supported by a full community-wide
  analysis.

## 3. Body-to-figure-legend consistency notes

- **Fig. 2A graph connectivity**: Body `submission/paper.tex:155` to
  `submission/paper.tex:156` says all 15,668 PHRs form a single connected
  component. Fig. 2 legend `submission/paper.tex:318` to
  `submission/paper.tex:321` says all 15,668 PHR flanks fall into a single
  component. These are directionally consistent. If edited, standardize whether
  the object is "PHRs", "PHR flanks" or "subtelomeric PHR flanks" so the number
  does not seem to count two different entities.
- **Jaccard matrix and blocks**: Body `submission/paper.tex:157` to
  `submission/paper.tex:163` and Fig. 2 legend `submission/paper.tex:322` to
  `submission/paper.tex:327` align well: the figure shows the same 41 x 41
  matrix ordered two ways, and both body and legend interpret the pattern as
  discrete blocks rather than a gradient.
- **Known-system validation plus novel q-arm grouping**: Body
  `submission/paper.tex:164` to `submission/paper.tex:170`, Fig. 2 legend
  `submission/paper.tex:327` to `submission/paper.tex:330`, and
  `ABSTRACT_nature.md:9` match the desired abstract pattern: recover Xp/Yp,
  Xq/Yq, acrocentric, 10p/18p and 4q/10q DUX4 systems; identify the linked
  22q/21q/19q/1q/13q/17q group. This is a strong consistency point.
- **PAR naming**: Body uses "the two pseudoautosomal pairs Xp/Yp and Xq/Yq" at
  `submission/paper.tex:164` to `submission/paper.tex:166`; Fig. 2 legend uses
  "PAR1, PAR2" at `submission/paper.tex:328`. Both are correct, but the body is
  more immediately readable. Consider pairing both once in a future edit:
  "PAR1/Xp-Yp and PAR2/Xq-Yq."
- **UPGMA phylogeny risk**: Body `submission/paper.tex:171` to
  `submission/paper.tex:173` has the right caveat. Fig. 2 legend
  `submission/paper.tex:323` to `submission/paper.tex:324` partially undermines
  it by calling the UPGMA ordering "the ``phylogeny'' of subtelomeric PHRs."
  This should be revised in any manuscript-editing pass because the legend is
  likely to be read independently from the body.
- **Fig. 3 representative architecture**: Body `submission/paper.tex:177` to
  `submission/paper.tex:187` and Fig. 3 legend `submission/paper.tex:342` to
  `submission/paper.tex:352` are consistent for C1/D4Z4-DUX4, C2/TUBB8B and
  C11/OR4F. The legend's opening at `submission/paper.tex:343` to
  `submission/paper.tex:345` correctly says "three communities," which is safer
  than the body's broader "each community" phrasing.
- **Global-return sentence after examples**: Body `submission/paper.tex:188` to
  `submission/paper.tex:192` is important because it turns the three examples
  back into a community-level conclusion. Fig. 3 legend has no equivalent
  global-return phrase, but that is acceptable for a figure legend. If space
  allows later, the Fig. 3 legend could add "These examples illustrate how
  community identity is carried by shared subtelomeric duplicon architecture"
  to prevent catalog-like reading.

## 4. Phrase-level risks and suggested alternatives

| Location | Risk | Suggested alternative |
| --- | --- | --- |
| `submission/paper.tex:161` to `submission/paper.tex:163`: "a community-detection algorithm partitions" | Accurate but algorithm-first. The house abstract voice makes the biological object the subject. | "These blocks define 15 arm-level sequence-similarity communities among the 41 signal-bearing arms..." |
| `submission/paper.tex:169` to `submission/paper.tex:170`: "they also reveal a previously uncharacterized q-arm grouping" | "Reveal" is acceptable but a little high-gloss; "previously uncharacterized" is defensible if no prior citation exists. | "they also identify a previously uncharacterized q-arm grouping of 22q, 21q, 19q, 1q, 13q and 17q." |
| `submission/paper.tex:171` to `submission/paper.tex:173`: "hierarchical clustering agreeing on 12 of the 15 communities" | Good caveat, but "hierarchical clustering" could be linked more explicitly to UPGMA for Fig. 2 consistency. | "An UPGMA ordering agrees on 12 of the 15 communities, and we use that ordering only as a grouping device rather than as a phylogenetic claim." |
| `submission/paper.tex:177` to `submission/paper.tex:178`: "each community carries a characteristic duplicon" | Overgeneralizes from three displayed communities and may sound like every community has one diagnostic locus. | "At the sequence level, representative communities are organized around characteristic duplicons." |
| `submission/paper.tex:179` to `submission/paper.tex:187`: C1/C2/C11 sentence sequence | Strong as examples, but could read as three disconnected locus entries if the opening sentence is not tightened. | Start the paragraph with "The representative communities in Fig. 3 show how the arm-level blocks correspond to shared subtelomeric architectures." |
| `submission/paper.tex:188` to `submission/paper.tex:191`: "gene content is dominated by..." | Good global-return sentence. Keep the contrast with "rather than community-specific protein-coding genes." | No required change; if shortened, preserve the contrast and the community tie-back. |
| `submission/paper.tex:323` to `submission/paper.tex:324`: "ordered by the UPGMA tree (the ``phylogeny'' of subtelomeric PHRs...)" | Highest-risk phrase in scope. Scare quotes are not enough; it invites unsupported phylogenetic interpretation in a stand-alone legend. | "ordered by a UPGMA dendrogram used only to group similar subtelomeric PHR profiles..." |
| `submission/paper.tex:324` to `submission/paper.tex:325`: "resolves the matrix into discrete blocks" | "Resolves" is acceptable here because it refers to visual block structure, not mechanism. | Keep, or use "shows" if reducing rhetorical force: "shows discrete blocks of co-similar arms rather than a uniform gradient." |
| `submission/paper.tex:327` to `submission/paper.tex:330`: "The blocks recover every known case..." | Strong and consistent with abstract. Add "shown here" only if a reviewer asks about possible omitted historical systems. | "These blocks recover the known systems shown here... and identify a tight..." |
| `submission/paper.tex:342` to `submission/paper.tex:345`: "Architecture of representative communities... so the duplicon content is legible" | Good legend voice; clearly marks examples as representative. | Keep. |

Bottom line: the body already keeps community language descriptive and includes
the necessary anti-phylogeny caveat. The only unsupported phylogenetic
implication that should be flagged for later editing is the Fig. 2 legend's
stand-alone phrase "the ``phylogeny'' of subtelomeric PHRs" at
`submission/paper.tex:323` to `submission/paper.tex:324`. Replace it with
"UPGMA dendrogram used only to group similar profiles" or an equivalent
descriptive phrase.
