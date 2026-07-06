# Voice review 04: genome-wide PHR extent and Fig. 1

Scope reviewed: `submission/paper.tex:136` to `submission/paper.tex:151`
and `submission/paper.tex:291` to `submission/paper.tex:309`, against the
house abstract voice in `paper_prep/synthesis/ABSTRACT_TEXTURE_SYNTHESIS.md`
and `paper_prep/synthesis/ABSTRACT_nature.md`.

## 1. Fig. 1 argument spine

The body section is already close to the house voice. Its strongest structure
is the sequence in `submission/paper.tex:138` to `submission/paper.tex:151`:
survey the full set of chromosome ends, name the expected positive controls,
then make the broader claim that the signal recurs across almost every other
end and reaches megabase scale outside acrocentrics and PARs.

Recommended spine to preserve in any edit:

1. Fig. 1a establishes near-ubiquity: telomere-anchored high-identity blocks
   are present at 41 of 48 chromosome arms, not only at historically expected
   loci. The phrase "nearly every chromosome end" at `submission/paper.tex:138`
   to `submission/paper.tex:140` is acceptable because the preceding section
   has already listed the seven no-signal arms at `submission/paper.tex:132`
   to `submission/paper.tex:134`.
2. Known positives validate the measurement: acrocentric short arms and PARs
   are correctly framed as expected positive controls at
   `submission/paper.tex:141` to `submission/paper.tex:142`. This pairs known
   biology with the broader new result, matching the house rule to use known
   systems as validation rather than as the whole claim.
3. The broader result follows immediately: "the same signal recurs at almost
   every other chromosome end" at `submission/paper.tex:142` to
   `submission/paper.tex:143` is the key discovery sentence. Keep this paired
   with the positive controls rather than moving it later.
4. Fig. 1b quantifies extent: the median and mean PHR lengths at
   `submission/paper.tex:144` to `submission/paper.tex:146` make the scale
   legible without implying exact biological termini.
5. The final scale sentence at `submission/paper.tex:148` to
   `submission/paper.tex:151` is useful because it excludes acrocentric short
   arms and PARs before giving 3.51 Mb, then uses PAR2 as an intuitive
   comparison. This is consistent with the abstract's "3.51 Mb outside
   acrocentric short arms and pseudoautosomal regions."

Small risk: the body section depends on the previous Methods-enabling
paragraph for the exact "41 of 48 arms" and seven no-signal list. That is
probably fine in manuscript flow, but if this subsection is read alone,
`submission/paper.tex:138` to `submission/paper.tex:143` could be tightened by
adding "41 of 48 arms" in the first sentence or by preserving the nearby
sentence at `submission/paper.tex:132` to `submission/paper.tex:134`.

## 2. Number and precision style recommendations

Use `3.51 Mb` in Results prose and figure legends. The current body sentence
uses `3.51 Mb` at `submission/paper.tex:149`, while the Fig. 1 legend uses
`3.510 Mb` at `submission/paper.tex:306`. The extra trailing zero in the legend
creates a false precision cue and diverges from the abstract-level house anchor.
Recommended style: `3.51 Mb` everywhere in this section and the Fig. 1 legend,
with `3.5 Mb` reserved only for broad talk-style summary prose.

Keep the median and mean as `median 105 kb, mean 144 kb`. This appears in the
body at `submission/paper.tex:145` to `submission/paper.tex:146` and the legend
at `submission/paper.tex:306`; it is compact, readable and consistent. Avoid
adding decimal places unless the figure itself reports them.

Keep the PAR2 comparison as a relative comparison, not a second headline
quantity. `submission/paper.tex:149` to `submission/paper.tex:150` says the
total is comparable in scale to PAR2 and that the median PHR is 31% of PAR2's
length. This works as an intuition aid. Do not expand it into a claim that PHRs
are pseudoautosomal-like in function; the comparison is length scale only.

Keep the seven no-signal arms explicit somewhere near the first genome-wide
claim. The exact list at `submission/paper.tex:132` to `submission/paper.tex:134`
is valuable because it prevents "near-ubiquity" from sounding like universality.
If the line is moved or compressed, preserve both the count and the limitation:
7 of 48 arms carry no detectable interchromosomal homology under this assay.

Use 464 HPRC v2 haplotypes plus CHM13 (= 465 near-complete assemblies) in the Fig. 1 legend. The current legend
at `submission/paper.tex:296` to `submission/paper.tex:297` does this correctly.
Use 465 near-complete assemblies only where a more abstract dataset summary is
needed.

## 3. Caveat language for flank-size truncation and PHR extent

The flank-size caveat is necessary and currently well placed. The body says
"the upper bound set by the 500 kb flank window rather than by biology" at
`submission/paper.tex:146` to `submission/paper.tex:147`, and the legend says
"saturation at 500 kb reflects the input flank size, not biology" at
`submission/paper.tex:306` to `submission/paper.tex:307`. This is the right
claim strength.

Recommended caveat pattern:

- Use "span tens to hundreds of kilobases within the assayed terminal 500 kb"
  when summarizing PHR extent.
- Use "right-censoring at 500 kb reflects the flank size" if a more technical
  figure-legend phrase is acceptable.
- Avoid "walking inward ... until inter-chromosomal homology disappears" unless
  immediately qualified, because `submission/paper.tex:144` to
  `submission/paper.tex:147` can sound as if the true biological endpoint was
  always observed. The current second clause fixes this; keep the caveat in the
  same sentence or the next one.
- Avoid "complete extent", "full extent", "entire PHR", "maximum PHR length" or
  "terminates" in this section. The assay bounds the observation to terminal
  500 kb flanks.

Best compact version for the body if edited later:

> Within the assayed terminal 500 kb, PHRs span tens to hundreds of kilobases
> (median 105 kb, mean 144 kb; Fig.~\ref{fig:fig1}b), with values at 500 kb
> right-censored by the flank size rather than by biology.

Best compact version for the legend if edited later:

> Non-acrocentric, non-PAR PHRs total 3.51 Mb (median 105 kb, mean 144 kb);
> bins reaching 500 kb are right-censored by the input flank size.

## 4. Phrase-level risks and suggested alternatives

`submission/paper.tex:138`: "A genome-wide survey reveals..."

Risk: "reveals" is a little more promotional than the house default. It is not
wrong here, but the house voice prefers measured verbs.

Suggested alternatives:

- "A genome-wide survey identifies telomere-anchored, high-identity blocks..."
- "We find telomere-anchored, high-identity blocks..."
- "Across the genome, telomere-anchored high-identity blocks occur..."

`submission/paper.tex:139` to `submission/paper.tex:140`: "each window coloured
by the number of other chromosomes it aligns to elsewhere in the pangenome."

Risk: good figure-reading language, but it describes the visual before the
claim. If tightened, put the claim first and keep the visual encoding in the
legend.

Suggested alternative:

- "We find telomere-anchored, high-identity blocks at nearly every chromosome
  end (Fig.~\ref{fig:fig1}a), with colour indicating the number of partner
  chromosomes."

`submission/paper.tex:141` to `submission/paper.tex:143`: "positive controls"
plus "almost every other chromosome end."

Risk: no major risk. This is the section's best house-style move because it
pairs known positives with the new genome-wide result.

Suggested alternative only if compression is needed:

- "Acrocentric short arms and PARs serve as expected positive controls, but the
  signal extends far beyond them."

`submission/paper.tex:144` to `submission/paper.tex:147`: "Walking inward from
each telomere until inter-chromosomal homology disappears..."

Risk: without the caveat, "until ... disappears" implies complete biological
resolution. The current sentence includes the caveat, but the phrasing still
starts from a boundary-finding image that is stronger than the 500 kb assay.

Suggested alternatives:

- "Within the assayed terminal 500 kb, PHRs span tens to hundreds of
  kilobases..."
- "Estimating PHR length within each terminal 500 kb flank..."

`submission/paper.tex:148` to `submission/paper.tex:151`: "a non-trivial
fraction of the genome."

Risk: "non-trivial" is vague and less Nature-compressed than the surrounding
quantitative prose. It also weakens the clean 3.51 Mb and PAR2 comparison by
adding a generic magnitude judgment.

Suggested alternatives:

- End the sentence after the PAR2 comparison.
- "Excluding acrocentric short arms and PARs, these PHRs total 3.51 Mb of
  subtelomeric sequence, putting the recurrent homology outside known exchange
  systems on megabase scale."
- "Excluding acrocentric short arms and PARs, these PHRs total 3.51 Mb, a
  megabase-scale reservoir of subtelomeric homology outside known exchange
  systems."

`submission/paper.tex:296`: "Interchromosomal" versus body
`submission/paper.tex:144`: "inter-chromosomal."

Risk: local inconsistency. The house guidance prefers "interchromosomal" in
prose. This should be standardized in a later manuscript edit, but this review
does not edit `paper.tex`.

Suggested alternative:

- Use "interchromosomal" consistently in prose and legends unless journal style
  requires the hyphen.

`submission/paper.tex:301`: "inter-chromosomal exchange landscape."

Risk: "landscape" is explicitly disfavored by the house style because it is
vague. "Exchange" also adds mechanistic flavor in the Fig. 1 legend before the
paper has reached 3D proximity or pedigree evidence. Fig. 1 measures
high-identity sequence sharing, not exchange events directly.

Suggested alternatives:

- "Telomere-anchored high-identity blocks mark the chromosome ends involved in
  recurrent interchromosomal homology..."
- "Telomere-anchored high-identity blocks identify the arms contributing to
  subtelomeric sequence sharing..."
- "Telomere-anchored high-identity blocks mark the genome-wide distribution of
  PHRs..."

Best compact legend replacement:

> Telomere-anchored high-identity blocks mark the genome-wide distribution of
> PHRs, with acrocentrics and PARs appearing as positive controls.

`submission/paper.tex:303` to `submission/paper.tex:305`: "within-end share per
10 kb bin across the terminal 500 kb."

Risk: technically compact, but "within-end share" may be opaque to readers.
If the figure axis already makes this clear, it is acceptable. If not, the
legend could spell out whether the share is the fraction of PHR bases or the
fraction of PHRs occupying each bin.

Suggested alternative:

- "PHR occupancy per 10 kb bin across the terminal 500 kb..."

`submission/paper.tex:306`: "3.510 Mb."

Risk: false precision and mismatch with `3.51 Mb` in the body and abstract
voice.

Suggested alternative:

- "3.51 Mb"

Overall recommendation: keep the section's argument order and quantitative
anchors, but revise later manuscript prose toward "we find/identify", standard
`3.51 Mb`, explicit 500 kb right-censoring, and a Fig. 1 legend that describes
PHR distribution rather than an "exchange landscape."
