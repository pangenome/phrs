# Reconciling Andrea Guarracino's DAD chr9q patch with the 9q/3q call

Branch merged: `pedigree-fig5-validation-and-patching` (merge `46d584d`). Her README:
`scripts/pedigree/fig5_v131_check/README.md`.

> **Correction (2026-11-18, A. Guarracino; checked on the raw v1.1 files as shipped).** The k-mer containment in
> §2-3 counted 21-mers on the forward strand only. The chr9 h2 block is a reverse-complement copy of chr1 hap1,
> so forward-only counting cannot see it. With canonical 21-mers (both strands), the block
> `PAN011.chr9.haplotype2:133,682,643-133,744,561` is **100.0%** contained in chr1 hap1 first 200 kb (forward-only
> 0.6%, which reproduces §2), 1.2% in chr1 hap2 (forward-only 0.5%) and 89.4% in chr3 hap2 terminal 579 kb (same
> both ways). The 20 kb tract `133,702,042-133,722,041` and the reverse complement of `PAN011.chr1.haplotype1:
> 22,528-42,527` (`samtools faidx -i`) have the same md5 (`e075205bcb5e`): 0 of 20,000 differences; forward-strand
> comparison gives 14,902. So the chr9 h2 block is a chr1p copy **and** chr3-like (chr1p and chr3q ends share
> sequence). The §2 verdict "not chr1-like" and the matching rows in §7 do not hold. Later work superseded §7:
> DAD chr9 corrected with a DAD-only hifiasm contig gives same-chromosome identity 1.0000 and the 9q/3q call is
> retracted (see `fig5_retraction_and_concerted_evolution_story.md`).

## §0 — Verified restatement (supersedes §2, §3 and §7 below on the chr1 question)

I verified A. Guarracino's correction independently. **She is right and my §2-3 were wrong**: I counted
21-mers on the forward strand only, and the inserted block is a **reverse-complement** copy, so a
forward-only count cannot see it. With canonical (both-strand) 21-mers:

| reference | canonical 21-mer containment of the chr9 h2 block (133,682,643-133,744,561) |
|---|---|
| **chr1 hap1 first 200 kb** | **100.0%** |
| chr3 hap2 terminal 579 kb | 89.4% |
| chr1 hap2 first 200 kb | 1.2% |

(positional check, independent of the containment count: block positions 0, 1,000, 20,000 and 60,000 map
to chr1 hap1 reverse-strand positions 138,032; 139,074; 158,074 and 197,948 - a constant offset,
collinear over 60 kb.)

So the released chr9q h2 terminus contains an **rc copy of a chr1p segment**, and it is simultaneously
89.4% chr3-like, because the chr1p and chr3q terminal blocks themselves share sequence. "chr1-like" and
"chr3-like" are therefore not alternatives: the block is a chr1p copy, and the family it belongs to is
shared with chr3q. Her identification was correct; my refutation was a strand error.

**Correction to her coordinates.** She places the copy at chr1 h1 8-61,926 and the 20 kb tract at the
reverse complement of chr1 h1 22,528-42,527. I do not reproduce that window: the collinear mapping I
measure is chr1 h1 ~138,032-197,948 (reverse strand), and the 20 kb window from 22,529 differs from the
tract at 14,064 of 20,000 positions. Her md5 `e075205bcb5e` does match the tract sequence itself, so the
difference is in the chr1 coordinates she paired it with - most likely a different release or contig. This
does not affect any conclusion; it only means the published description should either use the coordinates
verified here or cite her exact command.

**What this changes elsewhere.** Nothing in the retraction. The corrected DAD assembly (hifiasm contig from
DAD HiFi reads only) still matches the daughter's haplotype at graph dice 1.0000, chr3q is still second at
0.9990, and the 9q/3q call is still retracted - because the *released* segment was an inserted chr1p copy,
not the chromosome's own sequence.

## Two separate loci - and only one of her conclusions survives

| locus | her claim | verdict after independent testing |
|---|---|---|
| **chr5q (the 1p call)** | not an assembly error; the chr1p winner is a **FastGA `-f32` seed-loss artefact** | **CONFIRMED** (see §1) |
| **chr9q (the 3q call)** | the terminal 62 kb is a **copy of chr1 hap1**, an assembly error, and the 9q/3q signal is artefactual | **NOT SUPPORTED** (see §2) |

## 1. The chr5q call is an alignment artefact - she is right

Her `fastga_chr5_f.summary.tsv`, FastGA alone (no sweepga) on `PAN027#2#chr5` against the v1.1 joint
target:

| version | -f | chain end in query | target interval reached |
|---|---|---|---|
| v11 | **32** | **181,298,272** | 31,823,271-181,086,256 |
| v11 | 48 | **182,144,241** | 31,823,271-**181,932,123** |
| v11 | 64..255 | 182,144,241 | same |
| v131 | 32 | 182,144,260 | 31,823,271-181,932,123 |

At `-f32` - the value Fig. 5 uses - the chain stops ~846 kb short of the query end, so the terminal
chr5 window never receives a chr5:h2 candidate and chr1p wins by default. At `-f48` and above the chain
reaches 181,932,123, **the last base of the chr5 hap2 contig**, at identity 1.0000. Her diagnostic-site
check agrees: the SELF.p chr5 tract matches chr5 hap2 at 48 sites against 1 for chr1p.

**Consequence:** the 1p/5q "donor" was already withdrawn on pedigree grounds; this supplies the
mechanism, and it is an alignment-threshold artefact, not biology. The paper must not present it as a
candidate, and Methods must state the `-f` sensitivity, because the same artefact class could affect
other terminal calls.

## 2. The chr9q block is chr3-like, not chr1-like - her identification does not reproduce

21-mer containment (k=21, decisive, independent of any aligner parameter). Query: the chr9 hap2
terminal block v1.1 133,682,643-133,744,561 (61,919 bp, the region she says is a chr1 copy).

| reference | shared 21-mers |
|---|---|
| **chr3 hap2 terminal 579 kb** | **89.4%** |
| chr1 hap1 first 200 kb | 0.6% |
| chr1 hap2 first 200 kb | 0.5% |

Her claim that the block is "a copy of DAD chr1 haplotype 1 positions 8-61,926" and that the 20 kb
tract is "identical letter by letter to the reverse complement of DAD chr1 h1 22,528-42,527" does not
hold on the v1.1 assemblies as shipped. A direct FastGA search of the chr9 tract against chr1p at
`-f8`, `-f32` and `-f512` returns **zero alignments**, and naive positional comparison of the 20 kb
tract against the reverse complement gives 14,064/20,000 differences (~random). So this is not a
threshold effect like §1; the sequence simply is not chr1.

> **Correction (2026-11-18):** not reproduced with canonical (both-strand) 21-mers or with the reverse complement:
> 100.0% of the block's 21-mers are in chr1 hap1 first 200 kb, and the 20 kb tract equals the reverse complement of
> chr1 hap1 22,528-42,527 base for base (md5 `e075205bcb5e` on both). See the note at the top.

Her chr1 coordinates may belong to a different contig or release, and that should be resolved with her
directly before the claim is used anywhere. It is not reproducible from the shipped assemblies.

## 3. Her own patch proves the reads are chr3-like

`patch_v11/PAN011_h2_chr9_v11_patch.fa` is the 69,687 bp read consensus she substituted. 21-mer
containment of **the patch itself**:

| reference | shared 21-mers |
|---|---|
| **chr3 hap2 terminal 70 kb** | **96.4%** |
| chr3 hap1 terminal 70 kb | 91.2% |
| the OLD chr9 h2 70 kb it replaced | 90.7% |
| chr9 h1 terminal 70 kb | 0.5% |
| chr1 hap1 first 70 kb | 0.3% |

The consensus built from DAD ONT reads over that region is **more** chr3-like than the assembly it
replaced. The reads therefore support a chr3-like terminus; they do not support a chr1-like one.

## 4. Why her patched rerun loses the signal - the patch is circular for this test

Her README's method, step 2-3: reads are collected over **SELF.p chr9 136,139,447-136,209,149**
(`PAN027#2#chr9`) and the consensus is made with "SELF.p is only the frame". But a consensus framed on
SELF.p's coordinates, restricted to reads that map into that frame, can only reproduce SELF.p's
sequence. So the patched DAD chr9q end is guaranteed to agree with SELF.p, and her own result shows it:

> "In the 10 tract windows the same-chromosome match to DAD chr9 h2 now has identity
> 0.999250-1.000000 (Run B: 0.978153-0.998250)."

The interchromosomal win disappears (10/10 -> 2/10) because the same-chromosome match was driven to
~1.0 **by construction**. The patch therefore cannot discriminate between "the assembly was wrong" and
"the assembly was right and we overwrote it with the transmitted sequence". Her 2 remaining
interchromosomal windows are at margins -0.000500 and -0.000250 with identity 1.0000 to chr3, i.e. ties.

## 5. What remains genuinely unexplained, and belongs in the record

Her independent observations that our analysis must still account for:

1. **nucflag MISJOIN** at chr9 h2 133,673,490-133,686,589 - the start of the block. This is a
   release-QC call and is not explained by anything above.
2. **Depth arithmetic**: 16 reads on the chr9 block vs 73 at 100-200 kb before it; 60 in the
   corresponding zone elsewhere; 16+60 = 76 ~ one genomic copy. The standard explanation is that reads
   multi-map between the chr3 and chr9 copies and are lost from both - which is also what the HiFi
   `Dup` call reports - but the deficit is real and should be stated.
3. **106 read-assembly disagreements** in the last 330 kb, all inside the block (her criterion: >=3
   MAPQ>=20 reads, >=80% agreeing, base != assembly). Consistent with either polymorphism or
   mis-mapping between the two copies.
4. **MOM also carries the module** on one chr9q (matches SELF.p at 0.9988), with three MAPQ 60 MOM reads
   spanning 270-297 kb that carry chr9-type bases at 75-78 tract sites. This is evidence **for** a real
   inherited haplotype, not against it.

## 6. Correction to my own earlier test (her criticism is valid)

My read re-assignment (59 of 90 reads preferring "hap2 as assembled") did not include chr1p hap1 as a
candidate reference, so it could not distinguish "hap2 correct" from "hap2 equals a chr1p copy". That
was a real design flaw and the test should not be cited. The 21-mer containment in §2-3 replaces it and
is not subject to that objection.

## 7. Where this leaves the claim

- **chr5q / chr1p**: dead as a candidate; cause identified (`-f32` seed loss). Already withdrawn.
- **chr9q / chr3q**: the block is chr3-like (89.4% k-mer, 99.79% identity over 100% coverage), the reads
  are chr3-like (her own patch: 96.4%), and it is chr3-like and not chr1-like by three orders of
  magnitude. Her artefact explanation rests on a chr1p identity that does not reproduce and on a patch
  that is circular for this comparison. **The retained call stands**, but the manuscript must now state
  three things: the nucflag MISJOIN, the depth deficit, and that the terminal few kb are capped by a
  low-confidence segment.
- **Both** candidates now have a documented artefact risk in the *alignment* layer, which is a stronger
  reason than before to report the `-f` sensitivity in Methods.
