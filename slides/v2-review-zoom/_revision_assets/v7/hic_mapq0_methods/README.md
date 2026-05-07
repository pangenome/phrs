# Hi-C MAPQ0 Methods Slide

Task: `review-zoom-v7-hic-mapq0-method-slide`

Deck placement: insert before the first Hi-C/3D validation result slide, after
the existing method transition slide `10m` and before slide `10a`.

## Purpose

This asset gives fan-in a concise methods/background slide explaining why the
3D validation had to be reanalyzed for subtelomeric repeats. The slide should
make three points:

- The original standard filtering was biased by ignoring MAPQ0 reads.
- MAPQ=0 is needed because much of the contact signal lies in repetitive
  subtelomeric PHR sequence.
- Multimappers were retained carefully: one primary/randomly chosen placement
  per read or segment, not duplicate copies across every possible placement.

## Slide-Ready Assets

- `hic_mapq0_methods_flow.svg`: 16:9-compatible visual flow for a
  `figure-slide` body.
- `slide_text.md`: editable slide title, flow text, caveat, speaker notes, and
  footer source text.
- `SLIDE_PATCH.md`: exact fan-in insertion instructions and Typst snippet.

## Report Anchors

- `end-to-end-report/report/05_hic_validation.md:5-11` states that Hi-C,
  Pore-C, and CiFi disable MAPQ/multimapper filters, keep exactly one
  randomly chosen alignment for each multimapped read, do not duplicate across
  all placements, and interpret this as symmetric noise with unreliable
  pair-level contacts in repetitive regions.
- `end-to-end-report/report/05_hic_validation.md:7-9` gives the technology
  details: Hi-C uses `MIN_MAPQ = 0`, `RM_MULTI = 0`, Bowtie2 default mode with
  no `-k`/`--all`; Pore-C and CiFi use minimap2 default one-primary alignment
  behavior with no `--min-mapq` filter.
- `end-to-end-report/report/05_hic_validation.md:44-70` documents the
  flanking unique-sequence control and says flanking regions have no
  multi-mapping risk.
- `end-to-end-report/report/05_hic_validation.md:105-111` supports the
  interpretation that flanking Mantel signal reflects broader domain proximity
  rather than multimapping artifacts, and that independent Hi-C communities use
  O/E-normalized inter-chromosomal contact matrices.
- `end-to-end-report/report/06_dipc_validation.md:5-7` documents Dip-C MAPQ=0
  throughout (`sam2seg -q 0`, `hickit --min-mapq=0`), says the default
  MAPQ >=20 filter would discard 60-99% of subtelomeric reads, and states that
  BWA-MEM2 reports one primary alignment per read.
- `end-to-end-report/report/07_integrated.md:7-32` summarizes convergent
  evidence across Hi-C, Pore-C, CiFi, Dip-C, sperm 3D, RPE-1, and mouse.
- `end-to-end-report/report/07_integrated.md:34-42` explains the flanking
  paradox: PHR intervals create multimapping ambiguity, `RM_MULTI=1` discards
  these reads, 100 kb flanking regions are unique sequence, and the flanking
  signal is stronger than the duplicated PHR signal.
- `end-to-end-report/report/10_limitations.md:35-37` gives the correct caveat:
  all technologies retain multimappers at MAPQ0 with one randomly chosen
  alignment; this creates symmetric noise, preserves aggregate enrichments, and
  leaves individual pair-level contacts in repetitive regions unreliable.
- `subtelomeric_analysis_report.md:962-993` gives the main Hi-C/Pore-C/CiFi
  validation setup and result context.
- `subtelomeric_analysis_report.md:997-1023` documents the flanking
  unique-sequence control in the long report.
- `subtelomeric_analysis_report.md:1025-1058` documents the Mantel/contact
  normalization logic and comparable flanking signal.
- `subtelomeric_analysis_report.md:1516-1520` documents T2T Dip-C remapping and
  MAPQ=0 for maximum subtelomeric coverage.
- `subtelomeric_analysis_report.md:1618-1653` summarizes convergent 3D evidence
  and the flanking paradox.
- `subtelomeric_analysis_report.md:2081-2095` gives long-report limitations for
  3D validation, including multimapping and confound controls.

## Wording Constraints

Use these exact ideas on the slide:

- Required phrase: "the original standard filtering was biased by ignoring
  MAPQ0 reads"
- Standard filtering failure: "MAPQ >=20 / removing multimappers deletes most
  subtelomeric signal"
- MAPQ0 handling: "reanalyze with MAPQ=0"
- Multimapper handling: "one primary/randomly chosen placement per
  read/segment, not all placements"
- Interpretation: "aggregate community signal, not precise pair-level claims
  in PHR repeats"
- Control/normalization: "flanking unique-sequence controls and O/E/contact
  normalization"
- Caveat: "MAPQ0 adds symmetric noise"

Copy-safe one-line phrases for fan-in:

- "the original standard filtering was biased by ignoring MAPQ0 reads"
- "MAPQ >=20 / removing multimappers deletes most subtelomeric signal"
- "reanalyze with MAPQ=0"
- "one primary/randomly chosen placement per read/segment, not all placements"
- "aggregate community signal, not precise pair-level claims in PHR repeats"
- "flanking unique-sequence controls and O/E/contact normalization"
- "MAPQ0 adds symmetric noise"

Avoid these misleading phrasings:

- Do not say multimappers were copied to every valid genomic position.
- Do not imply MAPQ0 makes individual PHR pair contacts precise.
- Do not present flanking controls as a separate discovery method; they are the
  clean-mapping control for the PHR contact interpretation.
- Do not describe standard unique-map filtering as neutral for subtelomeric
  PHRs; the point is that it is biased here.

## Recommended One-Slide Claim

Title: "Making Hi-C work at subtelomeric repeats"

Headline: "3D validation required reanalysis because the original standard
filtering was biased by ignoring MAPQ0 reads."

Flow:

1. Standard unique-map filtering fails here: MAPQ >=20 / removing multimappers
   deletes most subtelomeric signal.
2. Reanalyze with MAPQ=0 and retain multimappers as one primary/randomly chosen
   placement per read/segment, not all placements.
3. Interpret aggregate community signal, not individual repetitive pair
   contacts; validate with flanking unique-sequence controls and O/E/contact
   normalization.

Caveat: "MAPQ0 adds symmetric noise; it does not support precise pair-level
claims in PHR repeats."

## Validation Notes

- Line anchors above were checked from the current worktree with `nl -ba` and
  `rg`.
- The SVG contains the required MAPQ0, one-placement, aggregate-signal, and
  pair-level caveat language.
- The patch file does not modify the final deck; it only instructs
  `review-zoom-v7-fanin-render` where to insert the slide.
