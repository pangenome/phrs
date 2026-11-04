# What the CpG island at the chromosome end actually is: the TERRA promoter

Our measurement: the first ~2 kb inside the telomere is a CpG island (observed/expected 1.02, GC 62.3%;
genome-wide CpG O/E is ~0.2-0.3), present at 88-95% of ends, independent of whether a protein-coding gene
is annotated within 5 kb. **This is a known, named element: the CpG-island promoter that drives TERRA
transcription.**

## What the literature says

**The element.** Human subtelomeres of XqYq and 10q (originally TelBam3.4 and TelSau2.0) share a conserved
repetitive region that runs centromere-to-telomere as a **61 bp repeat, a 29 bp repeat and a 37 bp repeat**
- the **61-29-37 repeats**. The 29 bp and 37 bp elements together form a DNA island **rich in CpG
dinucleotides**, a feature of many mammalian RNAPII promoters. At least **20 human chromosome ends** share
these repeats, in different tandem copy numbers (refs 14,16 of the Transcription review below).

**It is a promoter, functionally.** GFP reporter assays with progressive 5' deletions of a 61-29-37 tract
from XqYq show the **29 bp and 37 bp elements alone are sufficient to drive expression**, and the repeats
are bound in vivo by **phosphorylated (active) RNAPII**. On XqYq most TERRA transcription start sites lie
**27 nt downstream of the last 37 bp repeat and about 250 bp upstream of the telomeric tract**.

**It is what makes TERRA.** Nergadze et al. (RNA 2009) showed CpG-dinucleotide-rich islands "**shared among
multiple human chromosome ends**" promote transcription of TERRA, the telomeric repeat-containing RNA
transcribed from subtelomeres into the (TTAGGG)n tract.

**Prevalence matches our number.** Roughly **75% of analysed chromosome ends carry a CpG island** within the
first 15 kb, and TERRA promoters come in two classes - with and without a CpG island - of which only the CpG
class is DNA-methylation regulated. We measure 88-95% in the first 2 kb, i.e. slightly higher because we
sample the extreme tip rather than the first 15 kb.

**Chromatin.** Most human subtelomeres carry a **CTCF- and cohesin-binding site within ~1-2 kb of the TTAGGG
tract, adjacent to a CpG island** implicated in TERRA control; RNAPII enrichment extends from the CTCF site
toward the telomere. In CHM13, telomeres harbouring **TAR1** repeats are "almost uniformly punctuated" by
CTCF-bound elements and CpG-island promoters. An older element, the 3.6 kb **HST repeat** immediately
adjacent to (TTAGGG)n, likewise contains one conserved CpG-rich region.

**GC content is known to be extreme here.** Subtelomeric one-copy regions average 47.9% GC against 41.6%
genome-wide, but vary widely - the **1p one-copy region reaches 62.5%** GC while 3p falls to 37.5%. Our
62.3-65.5% at the tip sits at that upper extreme.

## Consequences for our analysis

1. **It is not a protein-coding gene, and it is not "telomere".** It is the TERRA promoter - a lncRNA
   promoter - built from a low-copy dispersed repeat family. The DDX11L genes we annotated at 1.7-2.5 kb are
   separate, further in.
2. **The element at the exchange hotspot is itself a shared repeat family.** 61-29-37 repeats (and TAR1) occur
   at 20 or more ends in different copy numbers. So the sequence immediately inside the telomere is a
   dispersed repeat shared between non-homologous ends - exactly the kind of substrate that licenses ectopic
   exchange. That is a mechanistic link between the terminal element and the inter-arm homogenisation we
   measure, and it was not visible from our data alone.
3. **It explains the allelic alignment collapse.** Different individuals carry **different tandem copy numbers
   of the 29/37 bp repeats** at a given end, so same-arm alignments in the first few kb are aligning repeat
   arrays of differing length and phase. That is why we saw two orders of magnitude fewer individuals and
   blocks ~18x shorter there, and why the ratio is undefined - not because the sequence is unreadable.
4. **It gives the empty zone a positive description.** Instead of "unmeasurable", the terminal few kb is
   "the TERRA promoter array, a copy-variable shared repeat family, whose variation in copy number and
   position between haplotypes is itself the phenomenon."

## Open

- **Whether Tan et al. discuss this is unverified.** Their full text was retrieved earlier but the cached copy
  has expired, and a fresh retrieval returned only the abstract. Their introduction does discuss TERRA and
  haplotype-specific (TTAGGG)n tract lengths, so it is plausible they treat it; it must not be assumed.
- Whether the CpG island we detect on the ~10% of ends without one is a variant promoter class, as the
  75%-with-islands / two-class literature would predict.
- Whether repeat copy number of 29/37 bp elements correlates with the paralog/allelic ratio of the same arm -
  the direct test of point 2.

## References

- Nergadze SG, Farnung BO, Wischnewski H, et al. CpG-island promoters drive transcription of human
  telomeres. RNA 2009;15:2186-2194. doi:10.1261/rna.1748309
- "Promoting transcription of chromosome ends." Transcription 2010;1(3). doi:10.4161/trns.1.3.13191
  (describes the 61-29-37 architecture, the 29/37 CpG island, >=20 ends, TSS positions)
- Deng Z, Wang Z, Stong N, et al. A role for CTCF and cohesin in subtelomere chromatin organization, TERRA
  transcription, and telomere end protection. EMBO J 2012.
- Expression and differential regulation of human TERRA at several chromosome ends. RNA 2019;25:1470
  (~75% of ends carry a CpG island within 15 kb; two promoter classes; Koskas et al. 2017 for the prevalence)
- Riethman lab, Mapping and initial analysis of human subtelomeric sequence assemblies (subtelomeric
  one-copy GC 47.9% vs 41.6% genome-wide; 1p 62.5%, 3p 37.5%)
- Single-molecule architecture and heterogeneity of human telomeric DNA and chromatin. bioRxiv 2022
  (TAR1 at CHM13 telomeres punctuated by CTCF elements and CpG-island promoters; 46 telomeres)
- Structural and transcriptional analysis of a human subtelomeric repeat. NAR 1991;19:149 (HST repeat,
  3.6 kb immediately adjacent to TTAGGG, one conserved CpG-rich region)
