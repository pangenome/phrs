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

## The symmetric test, run on the parental basis: PAN011 is NOT heterozygous at chr5q

The earlier demotion of 1p/5q was justified by reasoning that did not apply to the parental comparison
(no boundary-locked switch; population sharing; grandparents). The child-side demotion is separate and
rests on phasing, because PARTNER is unsampled. So the 9q/3q standard was applied here directly: in a
single individual, does the **donor-derived block sit on one haplotype only**?

Direct alignment of PAN011's four terminal 200 kb windows (FastGA -f512; best single alignment per pair):

| comparison | aligned bp | identity |
|---|---|---|
| chr5 **hap2** vs chr1 hap1 | 30,974 | 99.364% |
| chr5 **hap2** vs chr1 hap2 | 42,503 | 99.475% |
| chr5 **hap1** vs chr1 hap1 | 35,489 | 99.265% |
| chr5 **hap1** vs chr1 hap2 | 42,282 | 99.468% |
| **chr5 hap2 vs chr5 hap1 (allelic)** | 56,463 | **99.463%** |
| chr1 hap2 vs chr1 hap1 (allelic) | 42,040 | 99.113% |

**Both chr5 haplotypes carry the chr1p-like sequence**, at comparable span and identity, and their
paralogous identity to chr1p (99.27-99.48%) is **indistinguishable from the chr5 allelic identity
(99.46%)**. There is no haplotype asymmetry.

Contrast 9q/3q in the same individual: chr9 hap2 at 100% coverage and 99.79% identity to chr3 hap2, chr9
hap1 at **7.5% coverage** (fragments only).

**Verdict, on a stated and symmetric standard:** the discriminator is **haplotype asymmetry within one
individual**. 9q/3q has it; 5q/1p does not. 1p/5q is therefore a shared family present on both homologues
- allelic-like in structure and in divergence - not a donor swap onto one haplotype.

**Correction to earlier reasoning in this note:** the demotion does NOT rest on the grandparents being
unsampled, and does not rest on PARTNER. Those are irrelevant to a within-individual comparison. The only
relevance of the grandparents is to *when* a module entered a lineage, which is a different question and
should not be attached to the donor claim.

## Read support was also examined (previously never done at these loci)

Same ONT analysis as for 9q/3q, PAN011, four terminal windows:

| contig | reads in terminal 200 kb | reads in last 5 kb | alignments spanning >20 kb |
|---|---|---|---|
| chr5 hap1 | 298 | 62 | 167 |
| chr5 hap2 | 332 | 86 | 194 |
| chr1 hap1 | 318 | 85 | 160 |
| chr1 hap2 | 33,322 | 48 | 165 |

Longest single molecules: 243,232 bp (chr5 hap2, MAPQ 60), 233,832 bp (chr1 hap1, MAPQ 60), 226,835 bp
(chr1 hap2, MAPQ 60), 192,760 bp (chr5 hap1, MAPQ 60). So all four ends are well supported to the
telomere by single molecules - the assemblies are not unanchored there.

**Caveat:** chr1 hap2 accumulates 33,322 reads in its first 200 kb, ~100x the other three ends. That end
is a multi-mapping/repeat-rich region, so any divergence statistic computed from it should be treated
with suspicion. The identity calls above use best-block alignments, which are less affected, but this is
recorded.
