# Abstract promise voice review

## Verdict

Keep.

The title/author/abstract/keywords region in `submission/paper.tex:42` to
`submission/paper.tex:56` now matches the house abstract almost exactly:
`submission/paper.tex:53` is the LaTeX version of `ABSTRACT_nature.md:9`, with
the expected LaTeX escapes for percent signs and en-dashes. It follows the house
sequence from `ABSTRACT_TEXTURE_SYNTHESIS.md:62` to
`ABSTRACT_TEXTURE_SYNTHESIS.md:73`: biological obstacle, graph-enabled survey,
extent, communities, 3D contact, pedigree caution and model-level close. The
main review action is not to rewrite the abstract, but to ensure the body
continues to satisfy the promises it creates.

## Sentence-by-sentence notes

- Title, `submission/paper.tex:42`: keep. "Concerted evolution" names the
  biological thesis, and "unorthodox recombination" is broad but not vague in
  this context because the abstract and body specify non-allelic,
  interchromosomal, subtelomeric exchange. If the title is revisited later, do
  so for journal taste rather than abstract-voice mismatch.

- Author and affiliation block, `submission/paper.tex:44` to
  `submission/paper.tex:50`: keep for voice purposes. This region creates no
  abstract promises. The only non-voice style point is that
  `submission/paper.tex:50` uses "Memphis, 38163, Tennessee, USA"; if copyediting
  later touches affiliations, the more conventional order would put Tennessee
  before the postal code.

- Abstract sentence 1, `submission/paper.tex:52` to
  `submission/paper.tex:53`; source text at `ABSTRACT_nature.md:9`: keep.
  "Human subtelomeres carry recurrent duplications, gene families and structural
  variation" has the Nature texture recommended at
  `ABSTRACT_TEXTURE_SYNTHESIS.md:23`: compact biology first, no decorative
  "landscape" or "chart" language. "Difficult to characterize" is acceptable
  here because it describes the prior limitation, matching the guidance at
  `ABSTRACT_TEXTURE_SYNTHESIS.md:35`.

- Abstract sentence 2, `submission/paper.tex:53`; source text at
  `ABSTRACT_nature.md:9`: keep. The 466/465+CHM13 framing is coherent:
  "466 near-complete assemblies" is immediately parenthesized as "465 HPRC v2
  haplotypes together with CHM13", matching the stable-number rule at
  `ABSTRACT_TEXTURE_SYNTHESIS.md:70` and the body support at
  `submission/paper.tex:111` to `submission/paper.tex:113`. "Implicit pangenome
  graph" and "approximately 12%" are technical, but the sentence keeps them
  subordinate to the survey action and does not over-explain the method.

- Abstract sentence 3, `submission/paper.tex:53`; source text at
  `ABSTRACT_nature.md:9`: keep. The 41/48 arms and 3.51 Mb claims are coherent
  with the Results and Methods support at `submission/paper.tex:132` to
  `submission/paper.tex:134`, `submission/paper.tex:148` to
  `submission/paper.tex:150`, and `submission/paper.tex:465`. The precision is
  also consistent with the review-plan preference for 3.51 Mb in Results-level
  claims (`VOICE_REVIEW_PLAN.md:98` to `VOICE_REVIEW_PLAN.md:103`), while the
  talk conclusion rounds the value to about 3.5 Mb
  (`submission/notes/Session7-PopulationGenomics.en.srt:1415` to
  `submission/notes/Session7-PopulationGenomics.en.srt:1421`).

- Abstract sentence 4, `submission/paper.tex:53`; source text at
  `ABSTRACT_nature.md:9`: keep. "Define 15 sequence-similarity communities"
  uses the recommended descriptive community language from
  `ABSTRACT_TEXTURE_SYNTHESIS.md:65` to `ABSTRACT_TEXTURE_SYNTHESIS.md:66` and
  avoids implying phylogeny. The known-positive-plus-discovery move is strong:
  Xp/Yp, Xq/Yq, acrocentrics, 10p--18p and 4q--10q DUX4 validate the analysis,
  while the linked 22q, 21q, 19q, 1q, 13q and 17q q-arm group carries the
  discovery promise. LaTeX dash style is correct in `paper.tex`:
  `10p--18p` and `4q--10q` on `submission/paper.tex:53` render as en-dashes.
  The Markdown source uses plain hyphens at `ABSTRACT_nature.md:9`, which is
  acceptable for the source memo but should not be copied verbatim into LaTeX.

- Abstract sentence 5, `submission/paper.tex:53`; source text at
  `ABSTRACT_nature.md:9`: keep. The evidence ladder is correct: human Pore-C
  and Hi-C "show" increased 3D contact, while mouse meiotic Hi-C "shows" the
  zygotene peak. This matches the house rule to keep human and mouse evidence
  distinct (`ABSTRACT_TEXTURE_SYNTHESIS.md:67`) and the talk's separation of
  mouse zygotene evidence from human Pore-C/Hi-C
  (`submission/notes/Session7-PopulationGenomics.en.srt:1185` to
  `submission/notes/Session7-PopulationGenomics.en.srt:1223`). No direct human
  meiotic measurement is claimed.

- Abstract sentence 6, `submission/paper.tex:53`; source text at
  `ABSTRACT_nature.md:9`: keep. "Provides evidence consistent with recent
  exchange preferentially within these sequence communities" is the right
  strength. It incorporates the caution recommended at
  `ABSTRACT_TEXTURE_SYNTHESIS.md:44` and avoids the older "directly observe"
  overclaim flagged in `VOICE_REVIEW_PLAN.md:150` to
  `VOICE_REVIEW_PLAN.md:157`. It also matches the transcript's "positive" but
  "not conclusive" pedigree cadence
  (`submission/notes/Session7-PopulationGenomics.en.srt:1335` to
  `submission/notes/Session7-PopulationGenomics.en.srt:1357`).

- Abstract sentence 7, `submission/paper.tex:53`; source text at
  `ABSTRACT_nature.md:9`: keep. "Near-ubiquitous feature" matches the transcript
  conclusion (`submission/notes/Session7-PopulationGenomics.en.srt:1373` to
  `submission/notes/Session7-PopulationGenomics.en.srt:1385`). "Support a model"
  is the right mechanism-level verb under `ABSTRACT_TEXTURE_SYNTHESIS.md:40` to
  `ABSTRACT_TEXTURE_SYNTHESIS.md:43`; it does not convert nuclear architecture
  and recombination opportunity into a demonstrated causal loop.

- Keywords, `submission/paper.tex:56`: keep, with one optional coordination
  note. The list covers the core search terms and uses "pseudo-homolog region"
  consistently with the definition at `submission/paper.tex:129`. The only
  possible omission is "Pore-C", because the abstract names Pore-C and the main
  3D figure relies on it (`submission/paper.tex:211` to
  `submission/paper.tex:215`). This is optional metadata cleanup, not a voice
  defect.

## Body-section obligations created by the abstract

- The implicit-graph Results must support the abstract's 466/465+CHM13 and
  12%-sampling promise without drifting into method-first prose. Current support
  is at `submission/paper.tex:110` to `submission/paper.tex:125`; Methods support
  starts at `submission/paper.tex:413`.

- The PHR extent section must preserve the exact scale anchors: 15,668 PHRs,
  41 of 48 arms, 3.51 Mb outside acrocentric short arms and PARs, and flank-size
  truncation at 500 kb. Current support is at `submission/paper.tex:132` to
  `submission/paper.tex:150`, `submission/paper.tex:306`, and
  `submission/paper.tex:465`.

- The community section must continue to pair known validations with the new
  linked q-arm group, and must not treat community ordering as phylogeny. Current
  support is at `submission/paper.tex:160` to `submission/paper.tex:173`, with
  figure-legend support at `submission/paper.tex:322` to
  `submission/paper.tex:330`.

- The representative-architecture section must make Fig. 3 feel like examples of
  the 15-community result, not a separate locus catalog. Current support is at
  `submission/paper.tex:175` to `submission/paper.tex:192`.

- The 3D section must maintain the abstract's split evidence ladder: human
  Pore-C/Hi-C show association in three dimensions; mouse meiotic Hi-C supplies
  zygotene bouquet timing; neither proves human meiotic recombination. Current
  support is at `submission/paper.tex:194` to `submission/paper.tex:225`, with
  figure support at `submission/paper.tex:361` to `submission/paper.tex:386` and
  Methods support at `submission/paper.tex:510` to `submission/paper.tex:590`.

- The pedigree section must satisfy the abstract's cautious promise by keeping
  "evidence consistent with" close to the assembly-artifact caveat. Current
  support is at `submission/paper.tex:229` to `submission/paper.tex:243`; Methods
  support is at `submission/paper.tex:595` to `submission/paper.tex:608`.

- The final synthesis must keep "support a model" weaker than "demonstrate a
  mechanism." Current lines `submission/paper.tex:245` to
  `submission/paper.tex:257` should therefore be reviewed carefully by the
  downstream pedigree/conclusion task, because "maintained by ongoing
  non-allelic homologous recombination" at `submission/paper.tex:252` is stronger
  than the abstract's "contribute to" formulation.

- The body and legends should not reintroduce abstract-level terms the house
  memo avoids: "landscape" at `submission/paper.tex:301`, "phylogeny" at
  `submission/paper.tex:323` to `submission/paper.tex:324` even in quotes, and
  "caught in the act" at `submission/paper.tex:394` are downstream review points
  because they are stronger or more atmospheric than the abstract.

## Exact replacement suggestions

These are suggestions only; this task does not edit manuscript files.

- Suggestion for `submission/paper.tex:56`, if the keyword list is allowed to
  expand:

  ```tex
  \keywords{pangenome, subtelomere, pseudo-homolog region, concerted evolution, recombination, Hi-C, Pore-C, meiotic bouquet, HPRC v2}
  ```

- Suggestion for `submission/paper.tex:301`, to align Fig. 1 legend language
  with the abstract's avoidance of "landscape":

  ```tex
  high-identity blocks mark interchromosomal subtelomeric sharing, and the
  acrocentrics and PARs appear as positive controls.
  ```

- Suggestion for `submission/paper.tex:394`, to make the Fig. 5 legend match the
  abstract's pedigree caution:

  ```tex
  Evidence consistent with recent interchromosomal exchange in the WashU
  three-generation T2T pedigree.
  ```

- Suggestion for `submission/paper.tex:245` to `submission/paper.tex:253`, for
  the downstream pedigree/conclusion review rather than the abstract itself:

  ```tex
  Together these observations support a model in which shared sequence tracks
  three-dimensional proximity, bouquet-stage proximity creates recombination
  opportunity between non-homologous chromosome ends, and occasional exchange
  can add new shared sequence to the population.
  ```
