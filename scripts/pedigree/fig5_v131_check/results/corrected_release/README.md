# Corrected WashU pedigree assemblies (v1.1 and v1.3.1)

We corrected the released v1.1 and v1.3.1 assemblies of the four sequenced members of the WashU pedigree with each person's own reads only.
We never used another family member's data to change a person's assembly.

| sample | role in our paper |
|---|---|
| PAN010 | MOM (grandmother) |
| PAN011 | DAD (grandfather) |
| PAN027 | SELF (their daughter) |
| PAN028 | KID (PAN027's daughter, father unsampled) |

Haplotype 1 and 2 keep the names and contig order of the release (`SAMPLE#1`, `SAMPLE#2`).

## What we changed

**1. Local reassembly of flagged regions (G1).**
- We flagged 10 kb bins where the person's HiFi depth fell to at most a quarter of the contig median, or where the HiFi reads carried at least three single-nucleotide differences at depth at least 3.
- We merged flagged bins within 20 kb and skipped groups lying mostly in centromeric satellite arrays.
- For each candidate we reassembled the candidate +- 400 kb with hifiasm 0.25.0 from the person's HiFi reads, ONT reads or both. We took as replacement a contig that matched the assembly exactly over at least 5 kb on each side.
- We applied a replacement only when more of the person's HiFi reads and more of their ONT reads aligned better to it than to the assembled sequence.
- We kept a replacement only when it did not raise the number of assembly 21-mers absent from the person's Illumina reads (Merqury error k-mers) within 1 kb.
- We polished the replaced pieces with the person's HiFi reads (minimap2 map-hifi, bcftools haploid calls with quality at least 20, bcftools consensus).
- We did not apply any replacement whose window held at least 100 kb of N (rDNA and satellite placeholders). This skips no locus in v1.1 and 5 loci in v1.3.1.

**2. Two chromosome-end fixes.**
- **PAN011 haplotype 2, chr9q.** We replaced the terminal 150,287 bp with a contig that hifiasm assembled from PAN011's HiFi reads alone, joined at the end of a 120,784 bp exact match to the assembly. We then polished the inserted piece with his HiFi reads, which corrected 1 indel at depth 8. The released terminus is a reverse-complement copy of a segment of his own chr1p, and his long reads disagree with it. The same contig is used in v1.1 and v1.3.1.
- **PAN027 haplotype 1, chr22q.** We removed a 10,434 bp tail without telomeric repeats that lay beyond the q-arm telomere.

## Numbers

| version | changes | G1 pieces | end fixes |
|---|---|---|---|
| v1.1 | 318 | 316 | 2 |
| v1.3.1 | 311 | 309 | 2 |

Merqury QV (k = 21, the person's own Illumina reads) rose in all 16 haplotypes.

| version | sample | hap | QV before | QV after | error k-mers before | error k-mers after |
|---|---|---|---|---|---|---|
| v1.1 | PAN010 | 1 | 56.8302 | 57.1754 | 131,561 | 121,507 |
| v1.1 | PAN010 | 2 | 57.5283 | 57.5626 | 112,530 | 111,643 |
| v1.1 | PAN011 | 1 | 59.4089 | 59.6775 | 72,788 | 68,423 |
| v1.1 | PAN011 | 2 | 56.8758 | 57.2459 | 125,825 | 115,547 |
| v1.1 | PAN027 | 1 | 59.6227 | 60.0815 | 69,118 | 62,189 |
| v1.1 | PAN027 | 2 | 60.7802 | 60.8453 | 53,042 | 52,253 |
| v1.1 | PAN028 | 1 | 58.5256 | 59.9373 | 89,426 | 64,610 |
| v1.1 | PAN028 | 2 | 57.2455 | 57.3747 | 119,681 | 116,172 |
| v1.3.1 | PAN010 | 1 | 56.8222 | 57.1664 | 131,818 | 121,771 |
| v1.3.1 | PAN010 | 2 | 57.5283 | 57.5626 | 112,530 | 111,643 |
| v1.3.1 | PAN011 | 1 | 59.4083 | 59.6756 | 72,810 | 68,463 |
| v1.3.1 | PAN011 | 2 | 57.6663 | 58.1139 | 104,882 | 94,610 |
| v1.3.1 | PAN027 | 1 | 59.6687 | 60.1248 | 68,406 | 61,586 |
| v1.3.1 | PAN027 | 2 | 60.8181 | 60.8949 | 52,626 | 51,704 |
| v1.3.1 | PAN028 | 1 | 58.5194 | 59.9345 | 89,281 | 64,454 |
| v1.3.1 | PAN028 | 2 | 57.2298 | 57.3791 | 120,177 | 116,115 |

In v1.1, G1 filled three small gaps (PAN011 chr9 11,426 N, PAN028 chr19 5,000 N, PAN028 chr14 11,794 N) and removed one (PAN010 chr14 5,000 N). In v1.3.1 the number of N is unchanged in all 8 haplotypes.

## Files

- `v1.1/v11.SAMPLE.changes.tsv`, `v1.3.1/v131.SAMPLE.changes.tsv`: one row per change. Columns:
  - `step`: `G1`, `end_replace` or `tip_trim`.
  - `change_id`.
  - `orig_start`, `orig_end`, `orig_bp`, `orig_N`: interval in the released assembly.
  - `final_start`, `final_end`, `final_bp`, `final_N`: the same change in the corrected assembly.
  - `mode`: hifiasm input (hifi, ont or both), NA for the end fixes.
  - `hifi_votes_old_new` and `ont_votes_old_new`: reads closer to the assembly / reads closer to the replacement (for example 3/29).
  - `merqury_kmers_orig` and `merqury_kmers_final`: error k-mers within 1 kb.
  - `reason`.

  Coordinates are 0-based, end excluded. An insertion has `orig_start = orig_end`. `tip_trim` has `final_start = final_end` (the contig end).
- `v1.1/v11.SAMPLE.summary.tsv`, `v1.3.1/v131.SAMPLE.summary.tsv`: one row per haplotype. Columns:
  - contigs;
  - bp, N, N %, N runs (before and after);
  - non-ACGTN letters after;
  - Merqury error k-mers and QV (before and after);
  - change counts;
  - changes without final coordinates (0 everywhere).
- `original_contig_lengths.tsv`: contig lengths of the released assemblies.
- `manifest.tsv`: one row per corrected FASTA. Columns: version, sample, haplotype, path on our cluster, bytes, contigs, bp, md5 of the bgzipped file.
- `figures/`: 16:9 summary images (method, QV, gaps, where the changes are, the two end fixes).

## Reproducibility

One script (`scripts/pedigree/fig5_v131_check/g1_rebuild.sh`, one job per person, `g1_rebuild_person.sbatch`) rebuilds each corrected assembly from the released one. For all 8 assemblies the rebuild is identical to the files listed in `manifest.tsv`, compared contig by contig by md5 (PAN011 rebuilt again on 2026-11-26 after adding the HiFi polish of the chr9q piece).

## Notes for the release

- **PAN010 haplotype 2, chrXp, 0-4 kb.** The (CCCCAA)n tract is absent from her Illumina reads. We left it as is.
- **D4Z4 on chr10q.** The PAN011 haplotype-2 array has two lengths in his ONT reads: 8 of 16 spanning reads match the assembly (about 8 units) and 8 carry a 32.6-39.5 kb longer array. PAN027's paternal copy has the long array in all 45 spanning reads. This agrees with the 10q D4Z4 change you described (PAN011 42 and 9 units, PAN027 21). We did not change these arrays.
