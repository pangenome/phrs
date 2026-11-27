# Response-letter material on the pedigree (NOT manuscript text)

Everything here belongs in the response to reviewers. The manuscript carries only the corrected analysis and claims.
Naming: SELF = PAN027, DAD = PAN011, MOM = PAN010, KID = PAN028 (her father, PARTNER, is not sequenced).

## Point: the pedigree exchange candidates

**What we did.** We redid the whole pedigree analysis on assemblies that we corrected with each person's own reads.

1. **Assembly correction.** We corrected the version 1.1 assemblies of all four people with each person's own reads only.
   - We flagged bins where the HiFi depth falls or the HiFi reads disagree with the assembly.
   - We reassembled each candidate with hifiasm from the same person's reads.
   - We replaced a region only when more of the person's HiFi reads and more of their ONT reads support the new sequence.
   - We kept a replacement only when it did not raise the number of assembly k-mers absent from the person's Illumina reads (Merqury).
   - We polished the replaced pieces with the person's HiFi reads.
   - We fixed two chromosome ends: the father's chr9q haplotype 2 and a 10,434 bp non-telomeric tail on the daughter's maternal chr22q.

   We made 318 changes in total. The Merqury quality value rose in all eight haplotypes, and the error k-mers fell from 773,971 to 712,344. One script rebuilds the four corrected assemblies from the released ones, base for base. We applied the same procedure to v1.3.1 (311 changes) and share both corrected releases, with the list of changes, with the data producers.
2. **Transmission scan.** We repeated the whole-genome comparison of each transmitted haplotype with the transmitting parent's two haplotypes on the corrected assemblies.
   - We cut the query into 2 Mb pieces, because FastGA stopped short of some chromosome ends with whole-chromosome queries.
   - We applied the 18:18 mapping filter and classified 2 kb windows.
   - We drew runs of at least 10 kb at identity at least 0.95.
3. **Result.** No run joins two non-homologous chromosome ends outside PAR1 and the acrocentric short arms.
   - **Father to daughter:** PAR1 (138 kb), 47 runs between acrocentric short arms (830 kb), and two thin runs from interior sequence to acrocentric short arms. These two lie in the chr3 centromere region and in the chr1q pericentromeric satellite at 131 Mb, not at chromosome ends.
   - **Mother to daughter:** 9 runs, all between acrocentric short arms (182 kb).
   - **Daughter to granddaughter (her maternal copy, resolved per chromosome):** 22 runs, all between acrocentric short arms (396 kb).
   - **Leiden communities:** 80 of the 82 runs join arms of the same community, against 0.929 expected under a permutation null that keeps the per-arm run counts (95% interval 0.927-0.951, p = 0.0003). The two exceptions are the interior runs.

## Point: the chr9q/chr3q candidate

**It was an assembly error, and we retract it as an exchange.**
- The last ~60 kb of the father's released chr9q haplotype 2 are a reverse-complement copy of a segment of his own chr1p. The release QC marks a misjoin there, and his own long reads disagree with the assembled sequence.
- With his chr9q end rebuilt from his own HiFi reads alone, without any data from the daughter, her paternal chr9q matches his chr9q at graph dice 1.0000 across the called tract.
- His chr3q is second at 0.9985.
- His other chr9q haplotype does not carry this sequence at all: the 20 kb tract has no alignment within the last 3 Mb of that haplotype.

We now present this locus for what it shows: a non-homologous end (chr3q) carries the block that one of the father's homologous ends lacks, which is the concerted-evolution signal of Fig. 5d, not a transmission event.

## Point: the chr1p/chr5q candidate

**It was an alignment artefact.** At the seed-frequency cutoff of the original scan, the father's chr5q candidate never entered the alignment, so chr1p won by default. At higher sensitivity the chr5q homologue reaches the last base of its contig at identity 1.0000. On the corrected assemblies, with the query in 2 Mb pieces, no chr1p/chr5q run remains.

## Point: child-side calls and haplotype labels

The released PAN028 files name the two copies hap1/hap2 per chromosome. The maternal copy is hap2 on chr1, 3, 4, 5, 7, 9, 12, 14, 15, 18, 21, 22 and X, and hap1 on the rest. Our earlier child-side scan read hap1 as the maternal copy, which is correct on only 10 of 23 chromosomes. We now resolve the maternal copy per chromosome. We report only transmissions whose donor genome is sequenced, so PAN028's paternal copy is excluded.

## Point: assembly quality of the called intervals

The earlier answer intersected the 538 odgi-untangle patches with the release QC calls. That patch set is withdrawn together with the untangle analysis. We now answer in two ways:
- **Assembly correction.** We corrected the assemblies themselves with each person's reads, as described above.
- **Merqury error k-mers in the runs.** We counted assembly 21-mers absent from the person's Illumina reads inside the 82 final runs, on both the query and the donor side. We compared the counts with the genome-wide rate of each haplotype (numbers in Methods).

## Point: what this does not change

- The genome-wide result, inter-chromosomal subtelomeric sharing organised into communities across 232 individuals, does not depend on the pedigree.
- The pedigree now gives a clean negative result outside the obligate compartments, with PAR1 recovered as the positive control.

## Tone rules for the letter

- State the errors plainly; the reviewers identified real defects and we thank them for it.
- Keep the naming explicit (SELF = PAN027, DAD = PAN011, MOM = PAN010, KID = PAN028, PARTNER unsampled).
