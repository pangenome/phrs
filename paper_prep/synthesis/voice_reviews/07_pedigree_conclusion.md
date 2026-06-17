# Voice review 07: pedigree exchange and final synthesis

## 1. Pedigree claim-strength ladder

The current pedigree subsection has the right core order, but it should keep a
clear ladder between observed objects and inferred biology.

1. **Direct observation: untangle patch pattern.** Lines
   `submission/paper.tex:229` to `submission/paper.tex:236` directly support
   the claim that the WashU T2T pedigree analysis reconstructs parental
   haplotype transitions and finds 538 high-quality interchromosomal patches,
   of which 494/538 (92%) fall within the same sequence community against a
   marginal-aware null. This is the strongest claim and can use "we observe" or
   "we find".
2. **Class-level interpretation: apparent exchange categories.** Lines
   `submission/paper.tex:237` to `submission/paper.tex:239` identify
   "apparent gene-conversion tracts" and "crossover-like reciprocal exchanges".
   The qualifiers are doing useful work. Keep "apparent" and "crossover-like"
   because the data classify patch geometry, not a fully validated individual
   germline event.
3. **Event-level caution: positive but not definitive.** Lines
   `submission/paper.tex:240` to `submission/paper.tex:243` are well aligned
   with the house voice: "positive but not definitive" keeps the talk endpoint
   active without overclaiming, and the assembly-artifact caveat is placed next
   to the claim it qualifies.
4. **Population-level inference: recent exchange preferentially within
   communities.** The abstract phrase at `submission/paper.tex:53`, "provides
   evidence consistent with recent exchange preferentially within these sequence
   communities," is the correct upper bound for broad claims. In the body, the
   safest take-home is that the patch enrichment is evidence consistent with
   recent interchromosomal exchange, not direct proof that each patch is a
   germline recombination product.

Recommended local calibration: preserve the evidence sentence at
`submission/paper.tex:233` to `submission/paper.tex:236`, keep the caveat at
`submission/paper.tex:240` to `submission/paper.tex:243`, and avoid upgrading
the result beyond "evidence consistent with recent exchange preferentially
within sequence communities."

## 2. Final synthesis mechanism audit

Lines `submission/paper.tex:245` to `submission/paper.tex:249` currently state
a self-reinforcing loop as if each link has been directly demonstrated: shared
sequence tracks 3D proximity, meiotic bouquet proximity creates opportunity,
and recombination generates new shared sequence that reinforces proximity at
the next meiosis. The first link is supported by the human 3D association; the
meiotic opportunity link is supported by the mouse bouquet timing and the model;
the final population feedback link is plausible but inferred.

The safest house-voice framing is to make the loop explicitly a supported
model before listing its parts. This also harmonizes with the abstract synthesis
language: "support a model in which nuclear architecture and recombination
opportunity contribute to their concerted evolution" (`submission/paper.tex:53`;
`paper_prep/synthesis/ABSTRACT_TEXTURE_SYNTHESIS.md`, "Recommended final house
abstract"). A revision target could be:

> Together these observations support a model in which shared sequence, nuclear
> proximity and recombination opportunity reinforce one another: sequence
> similarity tracks three-dimensional proximity; telomere clustering during the
> meiotic bouquet creates opportunities for ectopic exchange; and recent
> exchange can add shared sequence back into the population.

The sentence at `submission/paper.tex:250` to `submission/paper.tex:253`
currently says the 3.5 Mb of sequence is "maintained by ongoing non-allelic
homologous recombination." That is the highest-risk mechanism phrase in the
scope because "maintained by" reads as a demonstrated causal mechanism, while
`submission/paper.tex:254` to `submission/paper.tex:257` immediately reopens
causal direction and notes that human 3D data are somatic and the meiotic-stage
measurement is mouse. Prefer "consistent with contribution from ongoing
non-allelic homologous recombination" or "supporting a model in which ongoing
non-allelic homologous recombination contributes to their persistence."

The final sentence at `submission/paper.tex:258` to `submission/paper.tex:261`
has the right positive close, but "population-scale confirmation" is too
definitive after "which way the central link runs remains open" at
`submission/paper.tex:254` to `submission/paper.tex:255`. The confirmed object
is not the full mechanism; it is the population-scale presence and organization
of relationships first inferred from cytogenetics. Safer alternatives:

- "Even so, the pattern provides a population-scale genomic counterpart, in
  complete assemblies, to relationships first inferred from cytogenetics decades
  ago."
- "Even so, complete genomes now place those cytogenetic relationships on a
  population scale."
- "Even so, the data extend relationships first inferred from cytogenetics to a
  population-scale view in complete genomes."

## 3. Fig. 5/body/Methods consistency issues

Fig. 5 explicitly conflicts with the current body in two places.

1. **Legend overstates the event-level evidence.** The Fig. 5 caption says
   "Inter-chromosomal recombination caught in the act" at
   `submission/paper.tex:394` to `submission/paper.tex:395`. The body is more
   careful: it says the analysis looks "for the recombination itself" at
   `submission/paper.tex:229`, finds patch patterns at
   `submission/paper.tex:233` to `submission/paper.tex:239`, and then states
   "positive but not definitive" at `submission/paper.tex:240` to
   `submission/paper.tex:243`. The legend should match the body by describing
   "pedigree-resolved interchromosomal patch patterns" or "evidence consistent
   with recent subtelomeric exchange," not validated recombination caught
   directly in the act.
2. **Legend promises missing CEPH1463/RPE-1 text.** The caption says
   "Cross-assembler validation (CEPH1463) and the RPE-1 t(X;10) positive
   control are described in the main text" at `submission/paper.tex:403` to
   `submission/paper.tex:404`. A current-body search finds no other occurrence
   of "CEPH1463" or "RPE-1" outside this legend. The scoped body section
   `submission/paper.tex:227` to `submission/paper.tex:261` does not discuss
   either analysis, and the pedigree Methods support at
   `submission/paper.tex:595` to `submission/paper.tex:608` also does not
   mention CEPH1463 or RPE-1. This is a hard consistency issue. Either delete
   the promise from the legend or restore only the minimal body/Methods support
   if those analyses remain part of Fig. 5.

Methods support is otherwise aligned with the present body: lines
`submission/paper.tex:597` to `submission/paper.tex:601` give the untangle,
quality-filter and classification basis; `submission/paper.tex:605` to
`submission/paper.tex:608` support the permutation null and Wilson interval
reported in the body and legend. One wording risk remains in
`submission/paper.tex:599`: "within-Leiden filter as credibility constraint"
could imply circularity because the headline result is the within-community
fraction. If this means filtering candidates to assess credibility, clarify the
operation before manuscript editing; if it means post hoc annotation rather
than exclusion, use "within-Leiden status used as a credibility annotation" or
"within-Leiden enrichment tested as a credibility criterion."

## 4. Phrase-level risks and suggested alternatives

| Phrase | Current location | Risk | Suggested alternative |
| --- | --- | --- | --- |
| "directly" / "observed directly" | The question at `submission/paper.tex:105` to `submission/paper.tex:106` asks whether recombination can be "observed directly"; the Fig. 5 legend at `submission/paper.tex:394` to `submission/paper.tex:395` says it is "caught in the act". | Outruns the event-level caveat at `submission/paper.tex:240` to `submission/paper.tex:243`. The direct observation is an untangle patch pattern, not definitive proof of each germline recombination event. | "can recent exchange be detected in pedigrees"; "pedigree-resolved evidence consistent with recent exchange"; "interchromosomal patch patterns in a complete pedigree". |
| "ongoing exchange" | The evidence ladder supports the abstract-level phrase at `submission/paper.tex:53`; the body close implies ongoing NAHR at `submission/paper.tex:252` to `submission/paper.tex:253`. | Acceptable as inference when tied to "suggests" or "evidence consistent with"; too strong as a standalone maintained-by mechanism. | "evidence consistent with ongoing exchange"; "suggesting recent exchange"; "compatible with ongoing non-allelic homologous recombination". |
| "maintained by" | `submission/paper.tex:252` to `submission/paper.tex:253`. | Converts the model into a causal conclusion, then conflicts with unresolved directionality at `submission/paper.tex:254` to `submission/paper.tex:257`. | "supporting a model in which ongoing non-allelic homologous recombination contributes to their persistence"; "consistent with contribution from ongoing non-allelic homologous recombination". |
| "confirmation" | `submission/paper.tex:260` to `submission/paper.tex:261`. | Sounds definitive for mechanism, even though the manuscript still needs germline-stage contacts and trio recombination maps at `submission/paper.tex:258` to `submission/paper.tex:259`. | "population-scale genomic counterpart"; "population-scale extension"; "complete-genome support for relationships first inferred from cytogenetics". |
| "self-reinforcing loop" | `submission/paper.tex:245` to `submission/paper.tex:249`. | Strong if stated as an observed process; safer if introduced as a supported model. | "support a model in which sequence homology, nuclear proximity and recombination opportunity can reinforce one another". |
| "apparent gene-conversion tracts" and "crossover-like reciprocal exchanges" | `submission/paper.tex:237` to `submission/paper.tex:239`; `submission/paper.tex:401` to `submission/paper.tex:403`. | These are appropriately cautious in the body. In the legend, they sit under an overstrong "caught in the act" lead. | Keep the qualifiers, but pair them with a cautious lead: "patches include apparent gene-conversion tracts and crossover-like reciprocal patterns." |
