# Response-letter material (NOT manuscript text)

Everything here belongs in the response to reviewers. The manuscript carries only the
corrected analysis and claims; it never comments on what was wrong before.

## Point: the pedigree exchange candidates

**What we found on re-examination.** Two defects in our own analysis, not in the data:

1. **Haplotype labels were read as parent-of-origin.** The released PAN028 assemblies
   name the two copies `hap1`/`hap2` per chromosome; the v1.0 parent-of-origin assignment
   differs between chromosomes (maternal copy = hap2 on chr1, 3, 4, 5, 7, 9, 12, 14, 15,
   18, 21, 22, X; hap1 on the remainder). Our child-side scans queried hap1 and treated it
   as the maternally inherited copy, which is correct on only 10 of 23 chromosomes:
   **228 of the 310 child-side patch calls (74%) are on chromosomes where that copy is
   paternal**, i.e. from the child's father, whose genome is not in the pedigree.
2. **Two events depended on those calls.** The "largest crossover-like event", a 27,969 bp
   chr3q->chr7p run, and a 2,089 bp chr3q->chr9q run are both on that paternal copy. On the
   daughter's actual maternal chr3 the same interval is 99.94% identical to her mother's
   paternal chr3q - ordinary transmission - with chr7p losing by 3.3%. Both are withdrawn.

**So the reviewers were right** that the array of pedigree "events" was an artefact. It was
an artefact of our labelling, and only the father-to-daughter comparisons (where both of
the daughter's parents are sequenced and her own labels are verified by parental k-mers)
support a claim.

**What survives, and why it is stronger than before.** One ectopic exchange: the daughter's
paternal chr9q terminus derives from her father's chr3q (99.79-99.84% identity, 100% coverage
of the terminal 60 kb), while his chr9q hap1 has no match there and his chr9q hap2 carries
only a drifted descendant (99.45-99.77%). Because he is heterozygous for the module, the
exchange predates his soma and arose in the paternal lineage's germline. No flagger call
contradicts it: HiFi flagger QC is clean at every interval, and the ONT `Col` calls covering
the terminal tens of kb appear at every chromosome end in every sample.

**chr1p/chr5q is not an event.** It is ancestral homology: the chr5q homologue explains the
sequence cleanly (99.58% over the full terminal 60 kb versus 97.9% for the chr1p copies) and
the chr1p-derived module is present in both grandparents' own chr5q.

**Re-counted statistics.** Applying the same filters to the correct copies gives ~316
patches with ~309 within a Leiden community (**98%**, previously 92% of 538). The
wrong-copy calls were disproportionately cross-community because they described an
individual outside the pedigree, so removing them sharpens the claim. The manuscript
Fig. 4 caption and the pedigree section were corrected accordingly; the child-side panels
were retired.

**Residual limitation we accept.** 54 of the 228 intervals cannot be tested at all without
the paternal genome (they fall outside the maternal flank because the haplotypes differ in
arm length). We state this in Methods and report only the transmissions with complete donor
genomes.

## Point: assembly-quality QC of the called patches

We re-ran the QC of all 538 WashU patch intervals against the assembly-QC calls for these
assemblies:

| query class | patches | overlap with a HiFi flagger call | rate |
|---|---|---|---|
| acrocentric | 459 | 1 | 0.2% |
| non-acrocentric | 79 | 0 | 0.0% |
| all | 538 | 1 | 0.2% |

The previously quoted "66 of 538" figure came from the ONT-based calls, which mark the
terminal tens of kilobases of every chromosome end in every individual, including arms
where no patch is called; that overlap is a property of ONT assemblies, not of the calls.
At HiFi resolution the called intervals are essentially untouched. We now state this in
Methods with the HiFi/ONT distinction explicit.

## Point: the chr1p/chr5p candidate "recurs, so it must be an artefact"

The recurrence is inheritance of a standing arrangement, not two independent events, and the
call is not a margin-supported exchange at all:

- On the maternal side the chr5q homologue explains the sequence cleanly: 99.58% identity
  over the full terminal 60 kb versus 97.9% for the chr1p copies.
- The chr1p-derived module is present in **both grandparents' own chr5q**, so it predates the
  family; it is ancestral homology, and we now describe it that way rather than as an
  exchange candidate.
- Our earlier claim that a different aligner "did not reproduce" the call was itself an
  artefact of that aligner's default chunking filter; with permissive parameters the chr1p
  homology is recovered at 99-100%. We have removed that argument from the manuscript, since
  the correct reason to downgrade the candidate is the standing-homology comparison, not
  aligner disagreement.

Net effect on the paper: the chr5q/chr1p candidate is no longer presented as an exchange, so
the "it happens twice" objection no longer applies to anything we claim.

## Point: what this does not change

- The genome-wide result - inter-chromosomal subtelomeric sharing organised into
  communities across 232 individuals - does not depend on the pedigree.
- The one retained exchange is supported by the donor's own competing chromosomes, not by
  a narrow identity margin alone.

## Tone rules for the letter

- State the error plainly and without hedging; the reviewers identified a real defect and
  we thank them for it.
- Do not characterise the surviving claim as "caveated" - it is the better-supported claim.
- Keep the naming explicit (SELF = PAN027, DAD = PAN011, KID = PAN028, PARTNER unsampled)
  so the label issue is unambiguous.
