# Voice review 03: implicit pangenome graph Results

Scope: `submission/paper.tex:108` to `submission/paper.tex:134`, with Methods
support from `submission/paper.tex:435` to `submission/paper.tex:465`.

## 1. Claim spine for the section

The section has the right job: introduce the enabling comparison that lets the
paper leave chromosome-by-chromosome subtelomere analysis and measure
interchromosomal sharing across all chromosome ends. Its strongest Results-level
spine is:

1. `submission/paper.tex:110`: "We treat every haplotype as its own reference."
   This is the cleanest Results entry point because it names the conceptual
   constraint without beginning in command-line detail.
2. `submission/paper.tex:111` to `submission/paper.tex:115`: the data scale and
   no-chromosomal-prior alignment. This should remain visible, but it can be
   more compressed: 464 HPRC v2 haplotypes plus CHM13, 18,827 500 kb
   telomere-anchored flanks, all-to-all alignment at 95% identity, no homologous
   chromosome restriction.
3. `submission/paper.tex:116` to `submission/paper.tex:120`: the implicit graph
   idea. Results only needs the reader to understand that alignments can be
   queried as graph reachability without building a single explicit graph. The
   equivalence proof and citations belong primarily in Methods.
4. `submission/paper.tex:121` to `submission/paper.tex:125`: the computational
   constraint-to-reach move. This is important for the talk logic, but the
   current sentence imports too much sampling theory. Results should say that
   sampled all-to-all alignment was sufficient to recover the subtelomeric graph;
   Methods can carry the `10^{24}` base-pair comparison estimate, 11.6% realised
   sampling rate, and `230x` Erdos-Renyi threshold.
5. `submission/paper.tex:129` to `submission/paper.tex:134`: the biological
   payoff. This should be the paragraph's landing point: PHRs are defined from
   transitive interchromosomal reachability and are found on 41 of 48 arms. The
   current final sentence is good Results material, though the seven-arm list may
   be optional if the next subsection immediately shows the genome-wide plot.

Recommended claim spine in prose form:

> We treat each haplotype assembly as its own reference and compare
> telomere-anchored flanks without chromosome labels. This produces an implicit
> pangenome graph: rather than build one explicit graph, we query transitive
> relationships among high-identity alignments. The approach overcomes the
> computational constraint that previously forced subtelomeric graph building
> into chromosome-partitioned components. Using these relationships, we define
> PHRs as telomere-proximal windows with high-identity reachability to multiple
> non-homologous chromosomes and find 15,668 PHRs on 41 of 48 chromosome arms.

This preserves the talk move from "we could not build a single graph reliably"
to "the implicit data structure lets us put all chromosome ends into one
relationship space" (`submission/notes/Session7-PopulationGenomics.en.srt:235`
to `submission/notes/Session7-PopulationGenomics.en.srt:290`) and then to
"understand first the scale" of the PHRs
(`submission/notes/Session7-PopulationGenomics.en.srt:671` to
`submission/notes/Session7-PopulationGenomics.en.srt:715`).

## 2. Results texture versus Methods detail

Essential to Results:

- `submission/paper.tex:110`: Keep. This is the section's best sentence and
  matches the house voice: compact, active, and conceptual.
- `submission/paper.tex:111` to `submission/paper.tex:115`: Keep, but compress.
  The assembly/flank scale and the absence of homologous-chromosome restriction
  are necessary to understand why the result is pangenomic rather than a
  chromosome-local comparison.
- `submission/paper.tex:116` to `submission/paper.tex:119`: Keep the idea, not
  the proof. A Results sentence should define "implicit pangenome graph" as a
  queryable alignment-reachability representation. The equivalence statement
  can be softer here and formalized in Methods.
- `submission/paper.tex:121` to `submission/paper.tex:123`: Keep only the
  computational constraint and the sampled-search solution. This is the hinge
  from method limitation to biological reach.
- `submission/paper.tex:129` to `submission/paper.tex:134`: Keep. The PHR
  definition and 41-of-48-arm result are the section's payoff. Consider moving
  the full threshold wording to Methods and giving Results a plainer definition.

Better left to Methods:

- `submission/paper.tex:116` to `submission/paper.tex:120`: "A pangenome graph
  and the set of alignments that define it are equivalent representations of
  the same object..." is too formal for the first Results section. Methods
  already has a concise version at `submission/paper.tex:450` to
  `submission/paper.tex:454`.
- `submission/paper.tex:121` to `submission/paper.tex:125`: "$10^{24}$",
  "11.6% of flank pairs", "`230x` the Erdos and Renyi threshold", and "single
  connected component" are Methods texture. They are important, but in Results
  they delay the biological result and make the section sound like a methods
  validation note. Methods already carries these facts at
  `submission/paper.tex:455` to `submission/paper.tex:458`.
- `submission/paper.tex:126` to `submission/paper.tex:128`: the genome-wide
  query QC is not part of the claim spine unless it supports a displayed figure
  or a reviewer-sensitive validation. It reads like a Methods/QC sentence and
  competes with the PHR result. If retained in Results, it should be reduced to
  one dependent clause or moved after the PHR definition as reassurance.
- `submission/paper.tex:129` to `submission/paper.tex:131`: the exact "at least
  five segments on at least two other chromosomes" threshold is Methods detail.
  Results needs the conceptual definition; Methods already gives the fuller rule
  at `submission/paper.tex:460` to `submission/paper.tex:464`, including the
  total aligned length threshold that the Results sentence omits.

## 3. Phrase-level risks and suggested alternatives

- Risk at `submission/paper.tex:111`: "From the 232 individuals we used 464 HPRC
  v2 haplotype assemblies together with CHM13v2.0" is correct but less aligned
  with the house abstract anchor than "465 near-complete assemblies (232 HPRC v2
  individuals and CHM13)." Suggested alternative: "Across 465
  near-complete assemblies -- 464 HPRC v2 haplotypes together with CHM13v2.0 --
  we extracted..." If manuscript style avoids dashes, use parentheses.
- Risk at `submission/paper.tex:113`: "aligned them all-against-all" is clear to
  specialists but can sound absolute when the next sentence says only 11.6% of
  flank pairs were evaluated. Suggested alternative: "queried all-to-all
  high-identity alignments with wfmash" in Results, while Methods explains the
  sampled realisation.
- Risk at `submission/paper.tex:116` to `submission/paper.tex:118`:
  "equivalent representations of the same object" is defensible but overly
  formal and may invite a reviewer to litigate graph theory in Results.
  Suggested alternative: "We used these alignments as an implicit pangenome
  graph, querying reachability among high-identity intervals rather than
  building a single explicit graph."
- Risk at `submission/paper.tex:118` to `submission/paper.tex:119`: "query it
  through their transitive closure" has a technical noun stack. Suggested
  alternative: "query which subtelomeric intervals are connected through chains
  of high-identity alignments." Keep "transitive closure" in Methods at
  `submission/paper.tex:450` to `submission/paper.tex:454`.
- Risk at `submission/paper.tex:121`: "Exhaustive all-against-all comparison is
  intractable" is a useful transition, but paired with `10^{24}` it becomes a
  proof. Suggested alternative: "Because exhaustive comparison across all
  flanks is computationally prohibitive, we used the sampled wfmash alignment
  graph to recover transitive relationships across chromosome ends."
- Risk at `submission/paper.tex:122`: "panmictic population" is unnecessary
  Results texture and may distract readers into population-genetic assumptions.
  Suggested alternative: move the population redundancy rationale to Methods,
  where `submission/paper.tex:446` to `submission/paper.tex:448` already houses
  it.
- Risk at `submission/paper.tex:124`: "Erdos and Renyi threshold" is not
  Results-level voice. Suggested alternative: omit from Results; Methods can
  retain the threshold and formula at `submission/paper.tex:455` to
  `submission/paper.tex:458`.
- Risk at `submission/paper.tex:126` to `submission/paper.tex:128`: the QC
  sentence introduces rDNA, centromeres, and sex chromosomes before the first
  PHR result. Suggested alternative if kept: "Genome-wide queries recovered the
  expected haplotype breadth, supporting use of the sampled graph for
  subtelomeric reachability." Otherwise move to Methods.
- Risk at `submission/paper.tex:129` to `submission/paper.tex:131`: the current
  PHR definition is threshold-first. Suggested alternative: "We define a
  pseudo-homolog region (PHR) as a telomere-proximal window connected by
  high-identity alignments to multiple non-homologous chromosomes." Then cite
  Methods for the exact segment and aligned-length thresholds.
- Risk at `submission/paper.tex:132` to `submission/paper.tex:134`: listing the
  seven no-signal arms is exact but slightly stalls the transition into the
  genome-wide figure. Suggested alternative: "This yields 15,668 PHRs on 41 of
  48 arms, leaving only seven arms without detectable interchromosomal
  homology." The explicit arm list can live in the Fig. 1 text, legend, or
  Methods if not needed at this moment.

## 4. Consistency obligations for Methods

- Assembly counts must stay synchronized across abstract, Results, and Methods.
  Use "465 near-complete assemblies" for the dataset-scale framing and define it
  as "464 HPRC v2 haplotypes together with CHM13" when enumerating inputs.
  `submission/paper.tex:111` to `submission/paper.tex:112` is internally
  consistent with the house abstract, but any rewrite should avoid switching to
  "465 HPRC haplotypes" or implying CHM13 is one of the 464 HPRC v2 haplotypes.
- If Results rounds to "approximately 12%" or says "sampled all-to-all
  comparison," Methods should retain the realised 11.6% value at
  `submission/paper.tex:455` and state how it relates to the all-vs-all
  alignment setup at `submission/paper.tex:442` to `submission/paper.tex:448`.
- If Results drops "Erdos and Renyi threshold," Methods should keep the
  connectivity justification at `submission/paper.tex:455` to
  `submission/paper.tex:458`. This lets Results stay biological while preserving
  the technical defense for reviewers.
- If Results uses a plain-language PHR definition, Methods must carry the exact
  operational thresholds. At present, Results (`submission/paper.tex:129` to
  `submission/paper.tex:131`) says five segments on at least two other
  chromosomes at 95% identity, while Methods (`submission/paper.tex:460` to
  `submission/paper.tex:464`) also adds total aligned length per window
  `>= 3 kb`. A later rewrite should either keep the Results definition
  conceptual or make clear that the complete threshold specification is in
  Methods.
- Methods should preserve the distinction between flank extraction and PHR
  detection. `submission/paper.tex:437` to `submission/paper.tex:439` defines
  the 18,827 telomere-bearing 500 kb flanks; `submission/paper.tex:462` to
  `submission/paper.tex:465` defines the scan that yields 15,668 PHRs. Results
  should not collapse these into "15,668 flanks" or imply all flanks are PHRs.
- The Methods phrase "absence of chromosomal partitioning" at
  `submission/paper.tex:457` to `submission/paper.tex:458` should support, not
  replace, the Results claim "without restricting alignment to homologous
  chromosomes" at `submission/paper.tex:114` to `submission/paper.tex:115`.
  Keep both concepts aligned: no chromosome labels as priors in Results; sampled
  connectivity justification in Methods.
- If the Results QC sentence at `submission/paper.tex:126` to
  `submission/paper.tex:128` is moved out of the section, Methods needs a
  destination for the genome-wide recovery check or it should be omitted
  consistently. Do not leave a Results claim unsupported by Methods or a Methods
  QC detail that appears to support no current claim.
