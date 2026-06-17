# B1-B3 3D Evidence Decision Record

Date: 2026-06-17

Scope: judgment support for manuscript-revision B-1, B-2, and B-3. This record
does not edit `submission/paper.tex`, captions, figures, bibliography, or figure
assets. It proposes decision paths and manuscript-ready wording blocks for an
author or later guarded patch task to accept, modify, or reject.

## Questions Under Decision

- **B-1:** Should the Results lead with flanking unique-sequence evidence, or
  keep PHR-internal pointwise sequence/contact scatters as the lead evidence?
- **B-2:** How should MAPQ filtering and retained multi-mappers be framed:
  apology, limitation, or preemptive design choice?
- **B-3:** How should PHR-internal scatters be labeled if the text foregrounds
  clean flanking controls?

## Evidence Inspected

Primary upstream decision context:

- `paper_prep/manuscript_revision/B0_3d_inventory.md`
- `paper_prep/manuscript_revision/01_fanout_graph.md`
- `paper_prep/manuscript_revision/02_operating_rules.md`
- Current manuscript context in `submission/paper.tex` around the human 3D
  Results, Fig. 4 caption, and Hi-C/Pore-C/CiFi Methods.

Relevant B0 evidence summary:

- PHR-internal pointwise scatters are the cleanest direct human measurement of
  the sequence-to-contact relationship: within one genome, per inter-chromosomal
  PHR sequence pair, community-free, with both axes from the same sample.
- HG002 Pore-C pointwise Spearman: `rho = 0.381`, `n = 2,830`, `p = 1.2e-98`.
  Replicates: CHM13 Hi-C `rho = 0.716`, HG002 Hi-C `rho = 0.662`, HG002 CiFi
  `rho = 0.191`.
- PHR-internal contacts remain mapping-sensitive because highly similar
  subtelomeric PHR tracts cannot be uniquely assigned by any MAPQ threshold.
- Flanking unique-sequence analyses are the cleanest anti-multi-mapping control.
  Community-free flanking Spearman is weaker but positive in key datasets:
  CHM13 Hi-C `rho = 0.136`, HG002 Hi-C `rho = 0.131`, HG002 Pore-C
  `rho = 0.038`, with CiFi not significant and fragmented assemblies showing
  sparsity/NaN behavior.
- Flanking arm-level Mantel remains positive in multiple datasets, including
  CHM13 `rho = 0.522`, HG002 Hi-C `rho = 0.520`, HG002 Pore-C `rho = 0.314`,
  with at least one weak/non-significant sample.
- Strict-MAPQ comparison is currently source-incomplete in the repo audit:
  `scripts/hic/mapq_strict_d_peerq1.py` exists and defines expected outputs, but
  B0 did not find committed result tables for the quantified strict-MAPQ claim.
- Current authority rules require this task to leave a decision record only.
  It must not silently resolve J tasks as manuscript edits.

## Recommended Decision

**Recommended path: keep the PHR-internal pointwise scatters as the lead human
3D evidence, but immediately frame MAPQ retention as an explicit preemptive
design choice and present flanking unique-sequence evidence as the primary
control that protects the interpretation. Label PHR-internal scatters as the
finer, downstream view of contact at the homologous tracts themselves, not as
clean read-origin proof.**

This is the strongest balance of biological relevance and defensibility:

- PHR-internal scatters directly test the paper's biological object of interest:
  whether more similar PHRs contact more often.
- Flanking unique-sequence signal is essential for credibility, but it is a
  control on nearby subtelomeric domains, not the same measurement as contact
  inside PHR sequence.
- Leading only with flanks would reduce the multi-mapping attack surface but
  would also make the central sequence-to-contact result look indirect and
  weaker.
- MAPQ0/random primary placement should be disclosed early as a necessary
  strategy for paralogous subtelomeres, paired with the flanking result as a
  pre-specified guardrail against a broad mapping artifact.
- The strict-MAPQ result should not be used as a quantified pillar unless B4/B5
  or a later fan-in locates or regenerates the missing result tables.

**Author decision required:** accept, modify, or reject this recommended
ordering. The downstream fan-in should not treat this recommendation as an
applied manuscript edit until an author or guarded patch task selects a path.

## Alternatives

### Alternative A: Lead With Flanking Unique-Sequence Evidence

Decision: make the flanking unique-sequence result the first human 3D evidence,
then present PHR-internal scatters as a downstream finer view.

Use this path if the authors want maximum reviewer preemption against
multi-mapping criticism, even at the cost of a less direct opening measurement.

Advantages:

- Starts from uniquely mappable sequence and is harder to dismiss as a MAPQ0
  artifact.
- Makes the MAPQ caveat feel designed into the analysis instead of appended
  after the result.
- Cleanly supports the claim that sequence-similar subtelomeric domains tend to
  occupy similar 3D neighborhoods.

Costs:

- Flanking effect sizes are weaker than PHR-internal sequence-level effects.
- Flanks are adjacent controls, not the homologous PHR tracts themselves.
- CiFi and fragmented-assembly behavior are less clean than the main
  PHR-internal pointwise scatter.
- The opening result may no longer map as directly to Fig. 4A unless the figure
  order or caption emphasis is changed later.

Author decision required:

- Decide whether the revision should prioritize defensive interpretability over
  direct biological measurement in the Results opening.

### Alternative B: Keep Current PHR-Internal Lead With Minimal Reframing

Decision: leave the PHR-internal pointwise scatter as the unqualified lead and
keep flanking/MAPQ material mainly in Methods.

Use this path only if the authors judge that the current Fig. 4A framing is
already sufficiently careful.

Advantages:

- Preserves the direct, intuitive Fig. 4 story.
- Minimizes manuscript restructuring.
- Keeps Results concise.

Costs:

- Reviewers can read MAPQ retention as a post hoc caveat rather than a
  preemptive design choice.
- Flanking controls may look buried despite being the best defense against a
  broad multi-mapping artifact.
- The strict-MAPQ sentence remains vulnerable unless source tables are verified.

Author decision required:

- Decide whether reviewer-risk reduction is worth adding visible Results
  language around the control apparatus.

### Alternative C: Two-Tier Lead

Decision: open the human 3D paragraph with one sentence on the flanking
unique-sequence control, then immediately transition to PHR-internal scatters as
the finer view at the homologous tracts.

This is a practical compromise if the authors like the recommended path but want
the first sentence to be maximally preemptive.

Advantages:

- Puts the clean-mapping control before the mapping-sensitive scatter.
- Still lets Fig. 4A remain the main visual result.
- Gives the reader the intended hierarchy: domain-level clean control first,
  PHR-internal high-resolution measurement second.

Costs:

- Adds a small amount of conceptual complexity to the Results opening.
- Requires careful wording so the flank result is not misread as the primary
  biological endpoint.

Author decision required:

- Decide whether the human 3D paragraph should start with the direct
  PHR-internal result or with the flanking control.

## B-1 Decision: Flanking Evidence As Lead Or Control

Recommended decision:

- Do **not** make flanking unique-sequence evidence the sole lead.
- Do make it a first-class control in the Results, close to the first human 3D
  claim.
- Use language that says flanks test whether the signal extends into adjacent
  uniquely mappable subtelomeric sequence and argues against broad multi-mapping
  inflation.

Claim boundary:

- Supported: sequence-similar subtelomeric regions show increased 3D contact,
  and this relationship is visible both at PHR tracts and in adjacent uniquely
  mappable flanks.
- Not supported as a lead claim: flanking data alone prove the exact
  PHR-internal contact rate or fully bound uniform MAPQ0 placement bias inside
  paralogous PHR tracts.

Author decision required:

- Select one of these evidence orders:
  1. PHR-internal scatter first, flanking control immediately after
     (recommended).
  2. Flanking control first, PHR-internal scatter as finer downstream view.
  3. PHR-internal scatter first, flanking control left mainly in Methods.

## B-2 Decision: MAPQ Filtering As Preemption

Recommended decision:

- Reframe MAPQ retention as a necessary preemptive analysis design rather than
  a defect discovered after the fact.
- State plainly that default MAPQ-filtered deposited maps cannot represent the
  relevant paralogous PHR contacts because the high-identity tracts are the
  signal-bearing sequence.
- Pair the disclosure with two constraints:
  1. PHR-internal pairwise contact is mapping-sensitive and should be interpreted
     as aggregate concordance, not as unique molecule origin.
  2. Flanking unique-sequence controls argue against a broad artifact from
     retaining multi-mappers.

Strict-MAPQ handling:

- Do not rely on the strict-MAPQ claim as quantified evidence unless downstream
  tasks locate or regenerate the promised `comparison_summary.tsv` and related
  output tables.
- If retained before verification, phrase it as a reproducibility expectation or
  analysis plan, not as a settled numeric result.

Author decision required:

- Decide whether the manuscript should keep a strict-MAPQ numerical claim. B0
  found the script/provenance, but not committed result tables.

## B-3 Decision: Label PHR-Internal Scatters As Finer Downstream Views

Recommended decision:

- Label the PHR-internal scatter as the finer, coordinate-level measurement at
  the homologous tracts themselves.
- Avoid wording that implies each PHR-internal contact read pair is uniquely
  assigned to its chromosome of origin.
- Use "aggregate", "coordinate-level", "PHR-window", or "PHR-internal" language
  rather than "cleanly mapped" language.

Claim boundary:

- Supported: at PHR coordinates, length-normalized aggregate contact increases
  with sequence similarity across inter-chromosomal PHR pairs.
- Not supported: MAPQ0/random primary alignment gives a read-level truth set for
  every paralogous PHR molecule.

Author decision required:

- Decide whether Fig. 4A caption and Results should explicitly call the scatter
  a "finer downstream view" after flanking control, or whether that label should
  stay in prose only.

## Manuscript-Ready Wording Blocks

These blocks are proposals. They are not applied edits.

### Accepted Path 1: Recommended Ordering

Use when the authors keep Fig. 4A / PHR-internal pointwise scatter as the lead
human 3D result, with visible flanking and MAPQ preemption.

Proposed Results block:

> In human contact-map data, contact also increases with PHR sequence
> similarity. We first measure this at PHR coordinates, per inter-chromosomal
> PHR pair within one genome and without community assignment or cross-sample
> averaging: HG002 Pore-C gives a pointwise Spearman correlation of
> $\rho = 0.38$, with the same relationship in CHM13 Hi-C, HG002 Hi-C, and
> HG002 CiFi. Because the most similar PHR tracts are also the least uniquely
> mappable, we treat these PHR-window measurements as aggregate
> coordinate-level evidence rather than read-origin assignments. As a
> preemptive control, we repeated the analysis in adjacent centromere-ward
> unique-sequence flanks; the weaker but concordant flanking signal argues
> against a broad multi-mapping artifact and supports the interpretation that
> sequence-similar subtelomeric domains contact one another more often.

Proposed Methods block:

> We disabled MAPQ filtering for the primary PHR-window contact analyses because
> default MAPQ-filtered processed maps remove the paralogous sequence that is
> central to the test. For each retained multi-mapping read, one primary
> placement was used, so PHR-internal contacts should be interpreted as
> aggregate contact at PHR coordinates rather than unique read-origin calls. To
> preempt broad inflation from multi-mapper retention, we repeated the
> community-free and arm-level analyses in adjacent centromere-ward
> unique-sequence flanks, where mapping is substantially less ambiguous.

Proposed Fig. 4A caption phrase:

> Each point is an aggregate PHR-window measurement for one inter-chromosomal
> PHR pair; this finer coordinate-level view is interpreted together with
> adjacent unique-sequence flank controls because highly similar PHR tracts are
> not uniquely mappable by MAPQ thresholding.

### Accepted Path 2: Flanking-First Ordering

Use when the authors decide B-1 should lead with clean flanking evidence.

Proposed Results block:

> We first asked whether the sequence/contact relationship is visible in
> uniquely mappable sequence adjacent to PHRs. In centromere-ward flanks, where
> the strongest paralogous PHR mapping ambiguity is avoided, sequence similarity
> remains positively associated with 3D contact in key human Hi-C and Pore-C
> datasets, although with smaller effect sizes than inside the PHR windows.
> This flanking signal is the primary control against a broad multi-mapping
> artifact. We then used the PHR windows themselves as a finer downstream view
> of the homologous tracts: in HG002 Pore-C, per-PHR-pair contact increases with
> PHR sequence similarity, and the relationship recurs in CHM13 Hi-C, HG002
> Hi-C, and HG002 CiFi.

Proposed Methods block:

> The human 3D analysis uses two tiers. Adjacent centromere-ward flanks provide
> the clean-mapping control, while PHR-window measurements provide the finer
> coordinate-level view of the homologous tracts. PHR-window analyses retain
> multi-mapping reads with one primary placement because a MAPQ threshold would
> remove much of the signal-bearing paralogous sequence. We therefore interpret
> PHR-window contacts as aggregate concordance and use the flanking analyses to
> test whether the result can be explained by broad multi-mapping inflation.

Proposed Fig. 4A caption phrase:

> PHR-window scatter showing the finer downstream view at the homologous tracts;
> adjacent unique-sequence flanks provide the clean-mapping control described in
> the text.

### Accepted Path 3: Minimal Reframing

Use only if the authors decide to preserve the current manuscript structure and
make the smallest defensibility edit.

Proposed Results sentence:

> Because the PHR tracts that define the strongest sequence relationships are
> not uniquely mappable by MAPQ thresholding, we interpret the PHR-internal
> scatters as aggregate coordinate-level concordance and use adjacent
> unique-sequence flanks as the control against broad multi-mapping-driven
> inflation.

Proposed Methods sentence:

> MAPQ-filtered deposited maps are inappropriate for the primary PHR-window
> analysis because they discard the paralogous sequence under test; the retained
> multi-mapper analysis is therefore paired with adjacent unique-sequence flank
> controls and should not be read as assigning every PHR-internal read to a
> unique chromosome end.

### Strict-MAPQ Conditional Wording

Use only if downstream B4/B5/BF fan-in verifies the strict-MAPQ result tables:

> A strict-MAPQ re-binning provides the expected negative control: PHR-internal
> contact collapses when paralogous placements are removed, whereas the
> adjacent-flank estimates are preserved within sampling noise. We therefore use
> strict MAPQ as a reproducibility boundary, not as the primary estimator of
> PHR-internal contact.

Use if strict-MAPQ tables remain unverified:

> Because strict MAPQ filtering removes the paralogous sequence that defines the
> PHR-window test, we do not use default MAPQ-filtered maps as the primary
> estimator. Instead, we pair retained-multi-mapper PHR-window measurements with
> adjacent unique-sequence flank controls that are less sensitive to PHR
> multi-mapping.

## Downstream Fan-In Guidance

For `manuscript-revision-bf-fanin`:

- Treat this as a decision record, not as an applied author decision.
- Preserve the recommended path unless a later author note selects an
  alternative.
- Do not promote strict-MAPQ numerical claims unless B4/B5 or another verified
  artifact locates the missing committed output tables.
- If drafting an integration recommendation, distinguish:
  - lead biological measurement: PHR-internal pointwise sequence/contact
    scatters;
  - primary mapping-artifact defense: adjacent flanking unique-sequence analyses;
  - support/robustness apparatus: community matrices, W/B, Mantel, O/E,
    multi-resolution, exclusion walks, and single-cell corroboration.

## Validation Notes

- Artifact exists: `paper_prep/manuscript_revision/B1_B3_3d_decision_record.md`.
- No manuscript, caption, figure, bibliography, or analysis-output file was
  edited.
- The record explicitly marks B-1, B-2, and B-3 author decisions required.
- It contains a concise recommended decision plus alternatives.
- It provides manuscript-ready wording blocks for the recommended path,
  flanking-first path, minimal-reframing path, and strict-MAPQ conditional paths.
