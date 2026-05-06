# Slidelet 14a spec: DUX4 / D4Z4

Status: candidate slidelet spec, not a rendered deck asset.

## Review Goal

Replace the current weak DUX4 crop with a content-first slidelet that answers:
What is DUX4, what is D4Z4, and why does this locus matter for PHRs?

## One-Line Takeaway

DUX4 is the gene embedded in D4Z4 repeat units; in this pangenome analysis, the
disease-relevant D4Z4/DUX4L architecture is concentrated in the C1 PHR
community at chr4q and chr10q.

## On-Slide Copy

Title:

DUX4 marks the D4Z4 PHR community

Left column, short definitions:

- DUX4: double-homeobox transcription factor; DUX4L copies sit inside D4Z4 repeat units.
- D4Z4: subtelomeric macrosatellite repeat array at chr4q and chr10q.
- FSHD relevance: pathogenicity depends on contracted D4Z4 on a permissive 4qA haplotype, not on every scattered DUX4-like annotation.

Right column, PHR evidence:

- C1 = chr4q + chr10q, the D4Z4 community.
- C1 haplotypes carry median 22 DUX4L copies.
- Non-C1 outliers carry only 0-2 copies on their own arm.
- Inter-chromosomal signal peaks at 0-15 kb from the telomere, where D4Z4 sits.

Footer:

Why this matters: PHRs recover a clinically named subtelomeric exchange system,
but the slide should use FSHD only as the biological anchor for chr4q/chr10q
sharing.

## Candidate Layout

Use a clean two-panel composition instead of the existing decorative crop:

1. Left: simple linear subtelomere schematic for chr4q and chr10q.
   Draw telomere at the left, D4Z4 repeat units immediately internal to it,
   then a fading PHR region extending centromere-ward. Label D4Z4 units as
   "DUX4L copies".
2. Right: compact box/point summary of DUX4L copy count by group:
   "C1 chr4q/chr10q" vs "non-C1 q-arm outliers". Put the median 22 callout
   next to C1 and "0-2 copies" next to outliers.
3. Bottom strip: "C1: chr4q/chr10q PHR sharing" plus the FSHD note in small
   text. Keep disease language secondary to PHR architecture.

Do not reuse the current single cropped panel unless its text is rebuilt. The
existing crop has provenance risk and does not explain enough on its own.

## Speaker Notes

DUX4 is the disease-named way into the chr4q/chr10q subtelomeric community.
The point is not to teach FSHD genetics in detail. The point is that the same
community the PHR graph calls C1 is the D4Z4/DUX4L system: chr4q and chr10q
carry the full macrosatellite array, with a median of 22 DUX4L copies per
haplotype, while other annotated q-arm hits are sparse. That makes DUX4 a useful
biological proof point for PHRs: a clinically visible subtelomeric exchange
locus drops out of the sequence-sharing analysis.

Suggested spoken version:

"DUX4 is the gene embedded in the D4Z4 macrosatellite repeat. We see DUX4-like
annotation beyond chr4q and chr10q, but the full D4Z4 array is the C1 community:
chr4q plus chr10q. That is why this belongs in the PHR story. FSHD tells us this
subtelomeric repeat system is biologically real; the pangenome shows its
population-scale architecture."

## Source Anchors

- `slides/v2/slide_14_gene_biology.md:7` gives the current slide claim:
  DUX4 annotated on 18 q-arms, C1 carries the full D4Z4 array, median 22
  DUX4L copies, other arms 0-2.
- `slides/v2/slide_14_gene_biology.R:5` through `slides/v2/slide_14_gene_biology.R:8`
  define the original panel intent.
- `subtelomeric_analysis_report.md:155` identifies C1 as chr4q/chr10q D4Z4
  macrosatellite sharing.
- `subtelomeric_analysis_report.md:531` maps Ambrosini block 7 D4Z4 to C1 and
  the DUX4L gene-level signal.
- `subtelomeric_analysis_report.md:573` explains that DUX4L pseudogenes are
  copies of DUX4 within D4Z4 repeat units and gives the FSHD/CTCF context.
- `subtelomeric_analysis_report.md:1671` through
  `subtelomeric_analysis_report.md:1680` gives the strongest C1 mechanistic
  support: 0-15 kb D4Z4 peak, median 22 DUX4L, 0-2 outliers, peripheral
  positioning, and FSHD-modifying chr4q/chr10q exchanges.
- `paper_prep/figures/ed8/caption.md:5` has a compact D4Z4-CTCF-lamin caption
  that can inform the schematic, but should not replace this slidelet.
- `subtelomeric_analysis_report.md:2384` through
  `subtelomeric_analysis_report.md:2387` list the local literature anchors:
  Lemmers 2010, Ottaviani 2009, and Masny 2004.

## Provenance Notes

The current visible `s14_dux4.png` is an exact-copy crop from the agent-878
zoom deck, not a reproducibly generated asset. The full source panel was
generated earlier from `slides/v2/slide_14_gene_biology.R`. See
`README.md` in this directory for the full git timeline and blob hashes.
