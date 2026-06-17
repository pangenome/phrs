# Abstract texture synthesis

## Recommendation

Use the Nature-compressed probe as the house abstract voice, with two targeted imports: the Talk probe's explicit caution around the pedigree evidence and the Biology probe's phrase "recombination opportunity" as a model-level term. The Nature probe has the best texture for a Nature submission because it moves in the target order already established by the committed abstract: biological obstacle, enabling graph method, genome-wide scale, community structure, 3D correspondence, pedigree evidence and model-level synthesis. It is compressed without sounding skeletal, keeps the quantitative anchors visible, and uses active evidentiary verbs rather than an explanatory talk script.

The Talk probe is the best cadence source. It has the cleanest evidence ladder for the risky endpoint: "positive but not definitive evidence" and "preferentially within the same sequence communities" match the transcript's caution that the pedigree is compatible with ongoing exchange but not definitive proof for an individual recombination event. Its weakness is length and spoken pacing. It explains each step rather than compressing them into publication form, so it should guide transitions and caution but not set the abstract's density.

The Biology probe has the strongest opening premise: subtelomeres create a setting where sequence homology, nuclear position and recombination opportunity can reinforce one another. That gives the paper a biological reason to exist before the method appears. Its weakness is that the method sentence becomes too technical for the abstract and the pedigree sentence says "direct evidence" before adding the caveat. That construction is defensible but riskier than the Talk probe's "evidence consistent with" formulation.

The current committed abstract in `submission/paper.tex` and `ABSTRACT_nature.md` is already close to the desired form, but it should inherit the synthesis' lower-risk endpoint. The phrase "further suggests recent interchromosomal exchange" is appropriately cautious; avoid returning to "directly observe interchromosomal exchange" in abstract-level prose unless the body has explicitly defined the observed object as an untangle patch rather than a validated germline event. The transcript should remain the cadence source: it ends with near ubiquity, population-scale observation, 3.5 Mb outside acrocentrics and PARs, sequence homology mirroring nuclear proximity, unresolved directionality and pedigree evidence that suggests ongoing exchange.

## Direct comparison of the three probes

| Probe | Best use | Main strength | Main risk | Verdict |
| --- | --- | --- | --- | --- |
| Nature-compressed | Primary house voice | Best Nature density; clean problem-to-result progression; strongest quantitative discipline | Slightly clipped; "resolve" and "support a model" need careful downstream use | Adopt as default texture |
| Talk-faithful | Cadence and claim-strength guardrail | Best spoken logic; best pedigree caution; makes the evidence ladder explicit | Too long and explanatory for the abstract; "ask how far" framing is less polished than the Nature version | Borrow transitions and caution |
| Biology-first | Opening biology and model vocabulary | Best biological premise; useful phrase "recombination opportunity" | Method sentence is over-technical; pedigree "direct evidence" phrasing invites overclaim | Borrow selectively |

## Phrases to borrow

- From Nature-compressed: "Human subtelomeres carry recurrent duplications, gene families and structural variation" is a strong, compact opener if the current "dynamic and structurally complex" phrase feels too generic. "without chromosomal or positional priors" is sharper than "without imposing chromosome or position labels." "define 15 sequence-similarity communities" is more economical than "form 15..." when the graph result is the grammatical subject.
- From Nature-compressed: "recover Xp/Yp, Xq/Yq, acrocentric, 10p-18p and 4q-10q DUX4 systems" is the right validation-plus-discovery move, as long as the manuscript keeps the LaTeX en-dash convention in `paper.tex`.
- From Talk-faithful: "positive but not definitive evidence" should be the house phrase for the pedigree endpoint outside very technical Methods prose. "fall preferentially within the same sequence communities" is useful because it reports the enrichment pattern before the inference.
- From Talk-faithful: "We then test whether this sequence structure is mirrored in three dimensions" is the best Results-section handoff into the 3D evidence. It is too procedural for the final abstract, but excellent for section review.
- From Biology-first: "sequence homology, nuclear position and recombination opportunity" is the best triad for the model, provided it is introduced as a model, not a demonstrated causal loop.
- From Biology-first: "contact one another more often in three dimensions" should remain the house wording for human Pore-C and Hi-C. It reports association and avoids claiming that human meiotic recombination was observed directly.
- From the current committed abstract: "near-ubiquitous feature of human chromosome ends" is the right closing scale phrase. Keep it.

## Phrase-level guidance

- "chart": Avoid. It sounds antiquated and can imply exploratory cartography rather than a specific computational measurement. Use "survey" for the dataset-wide measurement or "define" for PHR boundaries and communities.
- "survey": Keep for the abstract method sentence. It is accurate for population-scale sequence relationship measurement and has Nature-appropriate restraint. Do not overuse it in Results sections where a stronger measured verb is available.
- "characterize": Use sparingly. It is acceptable for the prior limitation, as in relationships being difficult to characterize, but it becomes vague as a result verb. Prefer "find", "define", "estimate", "recover" or "test" when the action is known.
- "resolve": Use only for structure that the analysis actually separates or orders, such as a linked q-arm group. It should not imply mechanistic resolution or proof of evolutionary direction.
- "map": Avoid in the abstract unless the manuscript is literally mapping coordinates. For relationships, "survey", "query", "define" or "compare" are more precise.
- "landscape": Avoid. It is a familiar but vague genomics abstraction and adds atmosphere without evidence. The house voice should be plainer.
- "reveal": Use rarely. It is acceptable for a genuinely new pattern, but it can sound ornate or promotional. Prefer "identify" or "find"; use "reveal" only when paired with a specific object such as the linked q-arm group.
- "demonstrate": Avoid for mechanism, human germline exchange and 3D-to-recombination interpretation. It is too strong unless the experiment directly proves the claim. Use "show" for measured associations and "support" or "suggest" for inference.
- "support": Keep for model-level synthesis. "Support a model in which..." is the correct strength for the nuclear architecture plus recombination interpretation. Avoid "support that" when it would imply proof.
- "suggest": Keep for the pedigree and unresolved directionality. It is the right caution level for human germline exchange. Avoid stacking too many "suggests" in one paragraph because it can make the paper sound weaker than the data are.
- "recombination opportunity": Keep as a model phrase. It captures telomere proximity plus homology without claiming that every observed shared tract is a validated recombination product.
- "direct evidence": Use only with an immediate qualifier, and preferably replace with "evidence consistent with" in the abstract. "Direct evidence consistent with" is less clean than "evidence consistent with" and may invite reviewer pressure.

## Overclaim audit

The Nature-compressed probe does not overclaim human germline exchange: "provides evidence consistent with recent exchange" is appropriately cautious. Its final sentence, "support a model in which nuclear architecture and recombination reinforce their persistence," is acceptable because "support a model" preserves uncertainty. Downstream reviewers should not strengthen this to "show that nuclear architecture maintains PHRs."

The Talk-faithful probe is the safest on claim strength. It explicitly says the pedigree provides "positive but not definitive evidence" and keeps mouse meiotic Hi-C separate from human Pore-C and Hi-C. Its conclusion says nuclear proximity and ongoing recombination "help maintain" concerted evolution; that is acceptable if preceded by "support a model," but too strong if used as a standalone claim.

The Biology-first probe is strongest rhetorically but has the highest overclaim risk. The opening "can reinforce one another" is fine as a framing hypothesis, but not as a settled mechanism. The phrase "direct evidence consistent with recent exchange" is slightly unstable because "direct evidence" leads before the qualifier. Prefer "evidence consistent with recent exchange preferentially within these communities." Its 3D wording is otherwise careful: human data show increased contacts, mouse data show a zygotene peak.

No candidate should claim direct human meiotic measurement. Human Pore-C and Hi-C support a contact association among sequence-similar subtelomeres. Mouse meiotic Hi-C supports stage timing around the bouquet. The pedigree supports recent exchange as an inference from high-quality interchromosomal patch patterns, not definitive proof of individual human germline recombination events.

## Recommended final house abstract

Human subtelomeres carry recurrent duplications, gene families and structural variation, yet their interchromosomal relationships have been difficult to characterize because chromosome ends were incompletely assembled and standard alignments partition homologous sequence by chromosome. Here we survey subtelomeric sequence sharing across 466 near-complete assemblies (465 HPRC v2 haplotypes together with CHM13) with an implicit pangenome graph that samples approximately 12% of haplotype-pair comparisons and queries transitive relationships without chromosomal or positional priors. We find extended interchromosomal homology on 41 of 48 chromosome arms, spanning tens to hundreds of kilobases and totaling 3.51 Mb outside acrocentric short arms and pseudoautosomal regions. These pseudo-homolog regions define 15 sequence-similarity communities that recover Xp/Yp, Xq/Yq, acrocentric, 10p-18p and 4q-10q DUX4 systems, and identify a linked q-arm group involving 22q, 21q, 19q, 1q, 13q and 17q. Human Pore-C and Hi-C show that sequence-similar subtelomeres contact one another more often in three dimensions, while mouse meiotic Hi-C shows that this relationship peaks at zygotene, when telomeres cluster in the bouquet. A telomere-to-telomere human pedigree provides evidence consistent with recent exchange preferentially within these sequence communities. Together, these data identify pseudo-homologous subtelomeres as a near-ubiquitous feature of human chromosome ends and support a model in which nuclear architecture and recombination opportunity contribute to their concerted evolution.

## House-style rules for section-review fan-out

1. Preserve the evidence order from the talk and current abstract: biological obstacle, graph-enabled measurement, PHR scale, community structure, 3D proximity, pedigree evidence and qualified model.
2. Start Results paragraphs from the biological question or figure-scale claim, then introduce the technical move only as far as needed to understand the evidence.
3. Use "we find" for sequence homology, PHR extent and community structure; "we observe" or "show" for measured contacts; "suggests" or "evidence consistent with" for human exchange; and "support a model" for mechanism.
4. Pair known systems with new discoveries. Xp/Yp, Xq/Yq, acrocentrics, 10p-18p and 4q-10q validate the approach; the linked q-arm group carries discovery value.
5. Keep community language descriptive. Do not imply phylogeny, ancestry or clade-level evolution unless the local sentence explicitly says the ordering is for grouping.
6. Keep human and mouse 3D evidence distinct. Human Pore-C and Hi-C show contact association; mouse meiotic Hi-C shows the zygotene bouquet timing.
7. Place caveats next to the claims they qualify, especially for the pedigree, then restate the defensible take-home rather than ending in apology.
8. Avoid ornate abstraction. Replace "chart", "landscape", broad "map" and routine "reveal" with measured verbs tied to the analysis.
9. Maintain the factual anchors unless a later data audit changes them: 466 near-complete assemblies, 465 HPRC v2 haplotypes plus CHM13, approximately 12% sampled comparisons, 41 of 48 arms, 3.51 Mb outside acrocentric short arms and PARs, and 15 sequence-similarity communities.
10. Use one spelling and one terminology layer per section: "interchromosomal" in prose, "pseudo-homolog region (PHR)" at first use where needed, and "pseudo-homologous subtelomeres" for the abstract-level biological object.
11. Do not promote the final model into a demonstrated mechanism. The house conclusion is that nuclear architecture and recombination opportunity are consistent with, and may contribute to, concerted evolution.
12. Keep the Nature target compressed. Borrow the talk's cadence for transitions and caution, but remove spoken scaffolding such as "we then ask", "finally" and "to close" unless it is doing real argumentative work.
