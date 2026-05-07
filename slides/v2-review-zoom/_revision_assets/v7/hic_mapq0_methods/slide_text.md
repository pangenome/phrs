# Making Hi-C work at subtelomeric repeats

## Header Claim

3D validation required reanalysis: the original standard filtering was biased
by ignoring MAPQ0 reads.

## Three-Step Flow

1. **Standard unique-map filtering fails here.** MAPQ >=20 / removing
   multimappers deletes most subtelomeric signal; Dip-C defaults would discard
   60-99% of reads at subtelomeric tips.
2. **Reanalyze at MAPQ=0.** Hi-C uses `MIN_MAPQ = 0`, `RM_MULTI = 0`;
   Pore-C/CiFi use no min-MAPQ filter; Dip-C uses `sam2seg -q 0` and
   `hickit --min-mapq=0`.
3. **Keep one placement, read aggregates.** Each multimapped read/segment keeps
   one primary/randomly chosen placement, not all placements. Interpret
   aggregate community signal using flanking unique-sequence controls and
   O/E/contact normalization, not precise pair-level claims in PHR repeats.

## Plain-Language Caveat

MAPQ0 adds symmetric noise. It supports aggregate community enrichment; it does
not support precise pair-level claims inside repetitive PHRs.

## Speaker Notes

This is a methods honesty slide. Standard unique-map filtering is usually a
defensible default, but it is biased for subtelomeric PHRs because the signal is
inside repetitive DNA. The reanalysis keeps the needed MAPQ0 reads while
avoiding systematic double counting: aligners emit one primary/randomly chosen
placement per read or segment, rather than expanding reads across every valid
locus. That choice makes pair-level PHR contacts noisy, so the deck should point
to aggregate community enrichments, O/E-normalized contact analyses, Mantel
tests, and the flanking unique-sequence controls.

## Footer Source Text

v7/hic_mapq0_methods/README.md; reports 05:5-11, 06:5-7, 07:34-42,
10:35-37
