# Slidelet 14c spec: TAR1

Status: candidate slidelet spec, not a rendered deck asset.

## Review Goal

Improve TAR1 readability and background. Avoid hard-to-read text-on-image
overlays. Explain TAR1 in plain terms and connect it to PHRs without claiming a
mechanism that the data cannot prove.

## One-Line Takeaway

TAR1 is Telomere-Associated Repeat 1: a subtelomeric satellite repeat that is
nearly universal in non-PAR PHRs and nearly absent from PAR1, making it a useful
marker of subtelomeric repeat architecture.

## On-Slide Copy

Title:

TAR1: a telomere-adjacent repeat marker

Definition block:

TAR1 = Telomere-Associated Repeat 1, a subtelomeric repeat element first
described near human telomeres.

Evidence block:

- 21,544 TAR1 entries.
- 14,816 of 15,668 PHR sequences carry TAR1 (94.6%).
- Present on all 41 arms with PHR signal.
- PAR1 is the exception: chrXp 0.3%, chrYp 1.1%, community C15 0.5%.
- 66.9% of TAR1 entries sit within 10 kb of the telomere; 70.3% within 25 kb.

Interpretation block:

TAR1 marks the telomere-proximal edge of the subtelomeric duplicated zone. It is
consistent with a shared repeat architecture, but the present data cannot prove
that TAR1 itself causes exchange.

## Candidate Layout

Use a readable chart-plus-text layout:

1. Left 65%: sorted bar plot of TAR1 prevalence by arm from ED3 panel a or a
   simplified rerender. Keep axis labels horizontal or use abbreviated arm
   labels with enough spacing.
2. Highlight only PAR1 arms in blue and acrocentric p-arms in orange; all other
   arms neutral gray.
3. Right 35%: three short metric cards:
   "94.6% PHR sequences", "PAR1 0.5%", "66.9% within 10 kb".
4. Optional bottom micro-strip: "TAR1 near telomere; PHR extends inward".
5. No paragraph overlays on the plot or background image.

If the deck needs one visual instead of a new chart, use ED3 panel a as the
visual source and ED3 panel d only as a speaker-note backup. Panel d is useful
for Q&A, but too detailed for the main slidelet.

## Speaker Notes

TAR1 is the repeat part of the gene-background slide. It is a
telomere-associated repeat, not a protein-coding gene. In these PHR calls, it is
almost everywhere: 21,544 entries, present in 94.6% of retained sequences across
all 41 arms. The exception is PAR1, where Xp and Yp recombine by the standard
pseudoautosomal mechanism and are essentially TAR1-free. TAR1 also sits very
near the telomere: two-thirds of entries are within 10 kb. So the clean claim is
that TAR1 marks the telomere-proximal repeat architecture of PHRs. Do not say it
drives exchange; the report explicitly treats that as hard to establish because
TAR1 is so widespread outside PAR.

Suggested spoken version:

"TAR1 is Telomere-Associated Repeat 1. It is not a gene; it is a subtelomeric
repeat marker. In our PHRs it is almost universal outside PAR1, and it sits very
close to the telomere. That makes it useful for explaining the architecture:
the terminal telomere, then TAR1-rich subtelomeric repeat, then the broader PHR
duplicated zone extending inward."

## Source Anchors

- `slides/v2/slide_14_gene_biology.md:9` gives the current TAR1 slide claim:
  21,544 entries, 94.6% of 15,668 PHR sequences, all 41 arms, PAR1 0.5%.
- `slides/v2/slide_14_gene_biology.R:100` through
  `slides/v2/slide_14_gene_biology.R:130` implements the current TAR1 bar plot.
- `paper_prep/figures/ed3/caption.md:3` gives ED3 panel a: 14,816 of 15,668
  sequences with TAR1, 94.6%, PAR1 low, acrocentric p-arms 73-79%, autosomal
  arms mostly over 99%.
- `paper_prep/figures/ed3/caption.md:9` gives ED3 panel d: most arms place
  TAR1 within 1 kb, with deeper positions on acrocentric p-arms and chr9q.
- `paper_prep/figures/ed3/sources.tsv:2` and
  `paper_prep/figures/ed3/sources.tsv:6` list the TAR1 prevalence and
  positional source tables.
- `subtelomeric_analysis_report.md:359` through
  `subtelomeric_analysis_report.md:367` explains how gene and TAR1 annotations
  were extracted and gives the 21,544 / 94.6% result.
- `subtelomeric_analysis_report.md:371` through
  `subtelomeric_analysis_report.md:381` defines TAR1, gives per-arm prevalence,
  and states the PAR1 interpretation and mechanism caveat.
- `subtelomeric_analysis_report.md:383` through
  `subtelomeric_analysis_report.md:391` gives the positional result:
  66.9% within 10 kb, 70.3% within 25 kb, telomere-proximal marker.
- `subtelomeric_analysis_report.md:397` through
  `subtelomeric_analysis_report.md:409` gives the per-community summary,
  including C15/PAR1 at 0.5%.
- `paper_prep/surveys/SURVEY_10_11_12_limits_summary_lit.md:368` through
  `paper_prep/surveys/SURVEY_10_11_12_limits_summary_lit.md:371` gives the
  local synthesis: PAR1 near-absence and passenger-status caveat.
- `subtelomeric_analysis_report.md:2379` lists Brown et al. 1990 as the local
  literature anchor for TAR1.

## Provenance Notes

The current visible `s14_tar1.png` is an exact-copy crop from agent-878. It is
not wrong, but it is too easy to make unreadable when combined with explanatory
text. Any replacement should rerender or clearly cite ED3-style TAR1 charts and
keep text separate from the plot. See `README.md` in this directory for the full
git timeline and blob hashes.
