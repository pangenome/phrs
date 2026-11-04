# chr1p / chr5q: what it is, and why it is not the 9q/3q case

Sources: `response/01_evidence_log_pedigree.txt` rev 10-12; `response/notes/9q3q_locus_full_report.md` §6-8;
cohort numbers from `hprc_heterolog_divergence_pairs.tsv`.

## What was found

A discrete chr1p-derived segment sits at the chr5q subtelomere:

- G2 candidate (PAN027 paternal hap2 vs PAN011): chr5q **182.052-182.080 Mb (28 kb)**.
- G3 (PAN028 maternal vs PAN027): a **~54-56 kb** region at chr5q 182.984-183.040 Mb that maps to
  PAN027's **paternal** chr5q in **four split blocks**, and to the maternal chr5q copy **not at all**.
- PAN011's chr5q carries the same arrangement, i.e. it is **standing**, not new in SELF.
- Population scale: chr1p-chr5q has **52.6% of aligned bases at >= 99.5% identity** - one of only three
  heterologous pairs above 50% - with a **49.5 kb** block at 99.529%, weighted d = 0.568% (T ~ 4.7 Myr).

## Why the child shows it without a new event

PAN028's maternal chr5q subtelomere descends from **PAN027's paternal chr5q** via an ordinary meiotic
crossover near the chr5q telomere. The crossover carried the standing chr1p-derived tract along with it.
So the child documents **transmission of a pre-existing arrangement**, not a second transfer.

## The decisive difference from 9q/3q: no boundary

Multi-alignment dice along the chr5q window:

| region | chr5 (homolog) | chr1p (paralog) |
|---|---|---|
| outside the candidate tract | 0.9960 | 0.9970 |
| inside the candidate tract | 0.9967 | 0.9990 |

The paralog leads **both inside and outside**, with no sign change and no step - a uniform offset, not a
mosaic boundary. At 9q/3q the homolog goes from 1.0000 outside to 0.9933 inside while chr3 stays flat at
0.9994, with the flip at the tract edge.

## Conclusion

**chr1p/chr5q is best described as an ancient small transfer/homogenisation tract that has become a
standing, widely shared arrangement** - same physical unit as 9q/3q (~50 kb subtelomeric segment from a
non-homolog), **different age regime**: old enough to have spread and lost its boundary, rather than a
discrete event datable to a meiosis. It is therefore *not* a second pedigree event, and calling it
"ancestral homology" (rev 12) understates it: homology of that specificity at 52.6% of bases above 99.5%
identity is itself evidence of past transfer, just not of a recent one.

## Aligner caveat (recorded in rev 10)

The chr5q/chr1p call is a **sweepga 0.1.1 (FastGA)** result. wfmash with default 50 kb chunking loses it
(the chr5q chunks map 100% to the father's chr5 hap2, so the embedded chr1p homology is diluted), while a
literal `-p95` targeted at the window recovers chr1 hap2 at 99.7-1.0000. The call is therefore
parameter-dependent and reproducible only with settings designed for short embedded tracts - which is
independent support for using FastGA rather than pangenome-default wfmash for cross-species work.

## Child phasing: why this cannot be pushed further

PAN028's **paternal** haplotype belongs to PARTNER, who was never sequenced, so nothing on that side is
interpretable - the reason the earlier child-side "events" (228 of 310 patch calls, 74%, on the wrong
copy) were withdrawn. In SELF both parents are sequenced, which is why the 9q/3q claim is anchored there:
it is the only chromosome copy in the pedigree where every required chromosome is sequenced.
