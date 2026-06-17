# G0-G3 Abstract, Introduction, Title and Denominator Decision Package

Date: 2026-06-17
Task: `manuscript-revision-g0`

This is a judgment package for downstream manuscript editing. It does not edit
`submission/paper.tex`. It synthesizes the A, B/F, C/D and E fan-in records for
the abstract, introduction, title, unorthodox-recombination framing, and
denominator wording.

## Source Records

- `paper_prep/manuscript_revision/A_mechanical_fixes.md`
- `paper_prep/manuscript_revision/BF_3d_contact_synthesis.md`
- `paper_prep/manuscript_revision/CD_continuum_community_synthesis.md`
- `paper_prep/manuscript_revision/C3_qarm_language.md`
- `paper_prep/manuscript_revision/E_pedigree_audit.md`
- `paper_prep/manuscript_revision/F1_F2_orphan_audit.md`
- Active text inspected in `submission/paper.tex`

## Accepted-Decision Placeholders

Use these placeholders in the final decision register. They are deliberately
phrased as editable author decisions rather than silent manuscript edits.

| ID | Decision placeholder | Recommended value for G package | Dependency |
|---|---|---|---|
| G-0 | `[ACCEPT_G0_TWO_TIER_CONTINUUM]` | Accept. Use "locally dense high-similarity peaks embedded in a broader subtelomeric similarity continuum" as the abstract/introduction concept. | C/D fan-in |
| G-0a | `[ACCEPT_G0_QARM_DENSE_NEIGHBORHOOD]` | Accept only if C0 heatmap-density result is accepted as descriptive evidence. The list 22q, 21q, 19q, 1q, 13q and 17q may remain only as an enriched q-arm neighborhood, not a closed group. | C0/C3 |
| G-0b | `[ACCEPT_G0_QARM_REMOVE_LIST]` | Alternative if C0 is judged insufficient/artifact-prone. Remove the q-arm list from abstract and figure legend; keep only unnamed q-arm similarity if needed. | C0/C3 |
| G-1 | `[ACCEPT_G1_TITLE_SOFTENING]` | Accept. Keep title-level mechanism broad; do not promote F_ST into the title or headline claim. | B/F and F1/F2 |
| G-2 | `[ACCEPT_G2_UNORTHODOX_RECOMBINATION_FRAME]` | Accept with hedging. "Unorthodox recombination" can remain as a model/framing term if the abstract says "evidence consistent with recent exchange" rather than "caught recombination" or "measured recombination rate." | E audit |
| G-3 | `[ACCEPT_G3_DENOMINATOR_SENTENCE]` | Accept. Use a single denominator sentence that distinguishes CHM13 span, HPRC-arm prevalence, and Leiden community count. | A and C/D fan-ins |

## Executive Recommendation

The safest front-matter framing is:

> Human subtelomeres show promiscuous inter-chromosomal sharing that forms a
> two-tier pattern: a broad continuum of low-to-intermediate sequence similarity
> with localized high-similarity peaks corresponding to known and newly
> highlighted sequence-similarity neighborhoods.

This language preserves the core result without implying that every community
is a closed clade. The words should be used consistently:

- **Continuum**: the global background of nonzero, variable inter-chromosomal
  subtelomeric similarity.
- **Peaks**: localized high-similarity neighborhoods on that continuum,
  including known systems and the conditional q-arm neighborhood.
- **Promiscuous**: the biological behavior of subtelomeric sequence sharing
  across non-homologous chromosome ends; not a synonym for "random" or
  "unstructured."

Avoid replacing this with "discrete blocks," "closed groups," "clades," or
"bounded classes" in front-matter prose. "Blocks" can remain only as visual
heatmap language if modified by "locally dense" and embedded in "continuum."

## Abstract Sentence Variants

Current abstract pressure points:

- It says the PHRs "define 15 sequence-similarity communities" without making
  clear that the 15 communities are arm-level Leiden communities among the 41
  signal-bearing arms.
- It says "identify a linked q-arm group," which overstates the q-arm evidence.
- It says mouse contact "peaks at zygotene," which B/F recommends softening to
  broad prophase-I support with the largest arm-collapsed point estimate at
  zygotene.
- It could leave the denominator relationship between 3.51 Mb, 41/48 arms and
  15 communities unclear.

### Variant A: C0 Dense-Neighborhood Accepted

Use if `[ACCEPT_G0_QARM_DENSE_NEIGHBORHOOD] = yes`.

```latex
We find extended inter-chromosomal homology on 41 of 48 chromosome arms,
spanning tens to hundreds of kilobases and totaling 3.51 Mb outside acrocentric
short arms and pseudoautosomal regions in CHM13. Across the 41 signal-bearing
arms, these pseudo-homolog regions form 15 arm-level sequence-similarity
communities: locally dense high-similarity peaks embedded in a broader
subtelomeric similarity continuum. The communities recover Xp/Yp, Xq/Yq,
acrocentric, 10p--18p and 4q--10q DUX4 systems; one enriched q-arm neighborhood
includes 22q, 21q, 19q, 1q, 13q and 17q.
```

Why this is acceptable: C0a/C0b support the q-arm list as an enriched
neighborhood. The sentence names the list only after the continuum/peaks
framing, avoiding "linked group," "tight grouping," and "closed class."

### Variant B: C0 Dense-Neighborhood Accepted, Short Abstract

Use if word count is tight but the q-arm list remains acceptable.

```latex
We find extended inter-chromosomal homology on 41 of 48 chromosome arms,
totaling 3.51 Mb outside acrocentric short arms and pseudoautosomal regions in
CHM13. These pseudo-homolog regions form 15 arm-level sequence-similarity
communities, with known systems and an enriched q-arm neighborhood appearing as
localized high-similarity peaks on a broader subtelomeric similarity continuum.
```

This keeps "peaks," "continuum," and "q-arm neighborhood" but removes the
specific q-arm list from the abstract. The list can still appear in Results or
Fig. 2 if authors want it.

### Variant C: C0 Insufficient or Artifact-Prone

Use if `[ACCEPT_G0_QARM_REMOVE_LIST] = yes`.

```latex
We find extended inter-chromosomal homology on 41 of 48 chromosome arms,
spanning tens to hundreds of kilobases and totaling 3.51 Mb outside acrocentric
short arms and pseudoautosomal regions in CHM13. Across the 41 signal-bearing
arms, these pseudo-homolog regions form 15 arm-level sequence-similarity
communities, with Xp/Yp, Xq/Yq, acrocentric, 10p--18p and 4q--10q DUX4 systems
appearing as localized high-similarity peaks on a broader subtelomeric
similarity continuum.
```

This removes the q-arm list entirely from abstract-level claims while preserving
the two-tier continuum result.

### Abstract 3D/Recombination Sentence

Replace "mouse meiotic Hi-C shows that this relationship peaks at zygotene" with
one of these, depending on final figure emphasis:

Preferred:

```latex
Human Pore-C and Hi-C support increased three-dimensional contact among
sequence-similar subtelomeres, and mouse meiotic Hi-C supports the same
sequence-to-contact relationship across prophase I, with the largest
arm-collapsed point estimate at zygotene.
```

Shorter:

```latex
Human and mouse contact maps support increased three-dimensional proximity among
sequence-similar subtelomeres during interphase and meiotic prophase I.
```

If the abstract needs one mechanism sentence:

```latex
Together with pedigree-resolved patch calls consistent with recent exchange,
these data support a model in which nuclear architecture creates opportunities
for unorthodox recombination among promiscuously shared subtelomeric sequences.
```

Use "support a model" and "consistent with recent exchange." Avoid "prove,"
"catch recombination in the act," "recombination rate," or "zygote-specific
peak" language.

## Introduction Framing Outline

The introduction should follow a four-step funnel. This keeps promiscuous
subtelomeric sharing, continuum/peaks, 3D proximity and exchange in the same
logical order as the evidence.

1. **Known biology: subtelomeres are promiscuous sequence mosaics.**

   State that human subtelomeres contain duplicated sequence shared among
   non-homologous chromosome ends. Use "promiscuous inter-chromosomal sharing"
   for the phenomenon and "pseudo-homolog regions" for the operational units.
   Keep examples: PAR1/PAR2, acrocentric short arms, D4Z4/DUX4, 10p/18p TUBB8B,
   OR4F. Do not call these "clades."

2. **Knowledge gap: old methods saw examples, not the genome-wide continuum.**

   The gap is not merely that the groups were unknown; it is that fragmentary
   cytogenetic/BAC/reference-era data could not measure how much sharing exists
   across all arms, whether it is a hard set of discrete classes or a continuum
   with peaks, and whether sequence similarity predicts nuclear proximity.

3. **New capability: population-scale T2T/pangenome view.**

   Introduce HPRC v2 plus CHM13, telomere-anchored flanks, implicit pangenome
   graph/transitive closure, and community summaries as a way to organize a
   connected graph. Use "summarize" or "organize" for community detection; avoid
   "discover closed groups."

4. **Hypothesis: nuclear architecture and unorthodox recombination opportunity.**

   Present the mechanistic model as an opportunity model: sequence-similar,
   promiscuously shared subtelomeric ends are more likely to contact one another,
   especially under telomere-clustering contexts, creating opportunities for
   ectopic exchange. The pedigree evidence is then framed as patch-level
   consistency with recent exchange, not as a rate estimate or independent
   meiotic replicate count.

Suggested last paragraph of the introduction:

```latex
We therefore asked whether promiscuous subtelomeric sharing is an exceptional
property of a few named chromosome ends or a near-genome-wide feature organized
as local high-similarity peaks on a broader continuum; whether sequence
similarity predicts three-dimensional proximity; and whether pedigree-resolved
assemblies contain patch-level signatures consistent with recent exchange.
```

## Title Options

### Recommended Title

```text
Concerted evolution and unorthodox recombination of human subtelomeres
```

Rationale: keep the current title if the abstract/intro soften the evidentiary
claims. "Concerted evolution" is supported by pervasive high-identity sharing
and community structure; "unorthodox recombination" is acceptable as a broad
mechanistic frame when the body says "consistent with recent exchange" and does
not overclaim rates or independent meioses.

### Softer Mechanism Title

```text
Concerted evolution and sequence exchange among human subtelomeres
```

Use if authors want to reduce the risk around "unorthodox recombination." This
title is safer but less distinctive. It matches the pedigree evidence without
requiring the title to carry a recombination-mechanism claim.

### Continuum-Forward Title

```text
A continuum of promiscuous sequence sharing among human subtelomeres
```

Use if the revision wants G-0 to dominate and reduce the mechanistic headline.
It is accurate but underplays the 3D and pedigree evidence.

### Avoid: F_ST-Promoted Title

Do not promote F_ST into the title, subtitle, or abstract lead. The F1/F2 audit
shows F_ST is present but demoted: AFR/non-AFR subtelomeric values are
approximately genome-wide background and support ancestry-preservation context,
not subtelomere-specific differentiation or ongoing exchange. A title such as
"population structure of human subtelomeric recombination" would point readers
toward the wrong evidence hierarchy.

### Avoid: Over-Discrete Title

Avoid titles using "clades," "closed groups," or "discrete communities." The
C/D fan-in supports a two-tier continuum with peaks, not hard partitions.

## Title Softening vs F_ST Promotion Decision

The decision is not "soften title or promote F_ST." F_ST should not be promoted
either way. The choice is:

- Keep the current title if authors accept "unorthodox recombination" as a
  model-level phrase backed by 3D proximity and pedigree patch evidence.
- Use the softer "sequence exchange" title if authors want the title to reflect
  only directly observed patch-level evidence.
- Keep F_ST as Methods/Extended Data context only if retained, with the matched
  genome-wide-background caveat attached.

Recommended accepted-decision text:

> We keep the title focused on concerted evolution and unorthodox recombination,
> but the abstract and introduction define the evidence level: sequence sharing
> is pervasive, 3D contact is concordant with sequence similarity, and pedigree
> patch calls are consistent with recent exchange. We do not elevate F_ST into
> the title because the matched-control audit supports background ancestry
> structure rather than subtelomere-specific differentiation.

## Denominator Reconciliation

Use one sentence early in the Results or abstract-side package:

```latex
In CHM13 coordinates, non-acrocentric, non-pseudoautosomal PHRs total 3.51 Mb;
across the population graph, inter-chromosomal PHR signal is detected on 41 of
48 chromosome arms, and Leiden community detection summarizes those 41
signal-bearing arms into 15 arm-level sequence-similarity communities.
```

This reconciles three denominators:

- **3.51 Mb**: a CHM13-coordinate span after excluding acrocentric short arms
  and pseudoautosomal regions.
- **41/48 arms**: prevalence of detected inter-chromosomal PHR signal across all
  chromosome arms in the population graph.
- **15 communities**: arm-level Leiden summary over the 41 signal-bearing arms,
  not a count over all 48 arms and not a sequence-level community count.

Optional longer denominator sentence if sequence-level count is also nearby:

```latex
The 3.51 Mb value is the CHM13-coordinate span of non-acrocentric,
non-pseudoautosomal PHRs; the 41/48 value is the number of chromosome arms with
detectable inter-chromosomal PHR signal across the population graph; and the 15
communities are the arm-level Leiden partition of those 41 signal-bearing arms,
distinct from the constrained 50-community sequence-level operating point.
```

## q-Arm Conditional Decision

If C0 heatmap-density and sequence-threshold results are accepted, use:

```latex
The same matrix highlights an enriched q-arm neighborhood, illustrated by 22q,
21q, 19q, 1q, 13q and 17q, within the broader continuum.
```

If C0 is not accepted, use:

```latex
Additional q-arm similarity is visible in the matrix, but we do not name a
bounded q-arm group from the present heatmap.
```

Do not use:

- "linked q-arm group"
- "tight q-arm grouping"
- "q-arm clade"
- "closed q-arm class"
- "the q-arm sextet" as a biological noun

## Downstream Patch Notes

These are not edits performed by this task; they are the exact choices a final
integrator should make after author/J decisions.

- Abstract: choose Variant A, B or C above.
- Introduction: add the three-question sentence so the reader expects
  continuum/peaks, 3D proximity and pedigree exchange evidence.
- Results/Fig. 2: replace "discrete blocks" with "locally dense blocks on a
  broader similarity continuum" and make the q-arm list conditional.
- 3D/mouse: replace "peaks at zygotene" with "broad prophase-I association"
  plus "largest arm-collapsed point estimate at zygotene" if the point estimate
  is retained.
- Pedigree: use "patch calls consistent with recent exchange" and state that
  the denominator is patch calls nested within three WashU child-haplotype
  transmissions if needed outside the abstract.
- F_ST: retain only as demoted context with matched-background caveat; do not
  use as title or abstract support for ongoing exchange.
- Denominators: insert the reconciliation sentence wherever 3.51 Mb, 41/48
  arms, and 15 communities appear close together.

## Validation

- Artifact exists: this file.
- Uses "peaks," "continuum," and "promiscuous" consistently as defined terms.
- Makes the q-arm list conditional on C0 heatmap-density/sequence-threshold
  acceptance.
- Addresses title softening versus F_ST promotion and recommends against F_ST
  promotion.
- Provides a denominator reconciliation sentence for 3.51 Mb, 41/48 arms, and
  15 arm-level communities.
