# Sequence similarity vs 3D proximity of PHRs: checks for the revision (R1-M4, R1-M5, R3-M3)

Status 2026-11-27 09:30 MST. All numbers come from the paper's scripts and pipelines, rerun on octopus. Tables and scripts are in `paper_prep/_brainstorming/sim_hic_null/`, `paper_prep/_brainstorming/matched_sd_control/` and `scripts/hic/`.

## 1. The question

The paper reports that chromosome ends whose PHRs are more similar in sequence are closer in the nucleus. The test is a Mantel test between arm similarity (Jaccard) and arm contacts. At 50 kb the Mantel is 0.657 for HG002 Hi-C, 0.656 for CHM13 Hi-C, 0.486 for HG002 Pore-C and 0.308 for HG002 CiFi.

The referees ask two things.

- R1-M4 and R3-M3: the correlation could be a mapping artefact. A read from one PHR can be placed on the similar PHR of another chromosome, and such errors create false contacts exactly between the most similar ends. They ask for analyses with uniquely mapped reads and for simulations.
- R1-M5: is the effect specific to PHRs, or does it hold for any pair of similar segmental duplications?

## 2. What the paper's pipelines do

- Hi-C (HiC-Pro 3.0.0, `MIN_MAPQ = 0`, `RM_MULTI = 0`, Bowtie2 in local mode) keeps multi-mapping reads. Bowtie2 gives them MAPQ 0 or 1 and places each at one of its equally good positions at random.
- Pore-C and CiFi keep only uniquely mapped reads. `pairtools parse2` 1.1.3 has `--min-mapq 1` by default, and our commands do not change it.
- The Methods (lines 1062-1069) say that we keep multi-mapping reads for all three assays. This is true only for Hi-C. A correction is proposed in `clean28_paper_edits.md`.
- The paper's Python environment (porec/1.0) had lost compiled modules. We rebuilt it from the exact package list, and the rebuilt environment reproduces the paper's tables byte for byte.

## 3. Hi-C with a mapping-quality filter

We rebuilt the Hi-C matrices keeping only pairs with MAPQ >= 1 or >= 2 on both reads. With Bowtie2 here, MAPQ >= 2 means uniquely mapped.

| | Pairs kept, MAPQ >= 1 | Pairs kept, MAPQ >= 2 | Mantel, all reads | Mantel, MAPQ >= 1 | Mantel, MAPQ >= 2 |
|---|---|---|---|---|---|
| HG002, 50 kb | 54.7% | 7.48% | 0.657 | 0.630 | 0.487 |
| HG002, 20 kb | | | 0.648 | 0.615 | 0.460 |
| CHM13, 50 kb | 97.1% | 87.1% | 0.656 | 0.636 | -0.023 (p 0.34) |
| CHM13, 20 kb | | | 0.721 | 0.717 | 0.021 (p 0.36) |

- HG002 keeps only 7.48% of pairs at MAPQ >= 2 genome-wide. The reference is diploid and the maternal and paternal copies are almost identical, so almost every read maps equally well to both.
- The CHM13 MAPQ >= 2 result is a loss of data, not a loss of signal. CHM13 PHRs are so similar across chromosomes that almost no read inside them is unique. For example, raw cis contacts per million drop from 57.7 to 4.1 in the chr22 p PHR and from 62.4 to 9.0 in the chr15 p PHR. Cooler balancing then masks these bins as too sparse. 24 of 42 CHM13 PHR loci lose every balanced bin, and the paper's statistic sees only 18 loci (152 of 835 pairs).
- The masked loci are not selected by similarity. Of the 84 most similar pairs (top 10% of Jaccard), 72 are lost (86%), against 683 of 835 pairs overall (82%). The mean Jaccard is 0.061 for the pairs lost and 0.050 for the pairs kept. (Jaccard values are explained in section 4.)
- On raw contact densities, similarity still follows contact at CHM13 MAPQ >= 2 (arm-level Spearman rho 0.438, from 3,327 contacts between PHRs).
- In HG002 the MAPQ >= 2 filter lowers coverage everywhere, so PHR bins are not outliers, keep their weights, and the Mantel stays at 0.487.

Balancing itself is not the problem for the paper's analysis, which keeps all reads: there PHR bins have normal coverage and balancing corrects coverage biases as intended. The problem appears only when a filter removes most reads inside PHRs. Then cooler's default bin filters treat PHR bins as missing data, and the test can no longer see the regions it is about. Any unique-read analysis of PHRs needs either raw densities or a balancing that keeps PHR bins.

## 4. Matched control: non-subtelomeric segmental duplications (R1-M5)

How we chose the controls.

1. We started from the SEDEF segmental-duplication tables of CHM13 and HG002 v1.1. Each row is a pair of copies, with identity and aligned length. We kept the pairs whose two copies lie on different chromosomes, as for PHRs.
2. We merged the copies on each chromosome into loci and described each pair of loci by three numbers: identity (weighted by aligned length), total aligned length, and copies (the mean number of other chromosomes that the two loci share duplications with).
3. A PHR pair has both loci overlapping a PHR. CHM13 has 347 PHR pairs and HG002 has 734.
4. A control candidate has both loci at least 1 Mb from either chromosome end.
5. For each PHR pair we took up to 3 candidates closest in identity, log length and log copies, each used once, within 0.5 standard deviations (caliper). 158 of the 347 CHM13 PHR pairs found a match, giving 426 controls (HG002 1,069 controls). The 189 unmatched CHM13 PHR pairs are the most identical and longest (mean identity 0.972), and no duplication away from the ends is as similar and as long.
6. The matched PHR pairs and the controls are close. In CHM13 the mean identity is 0.954 for matched PHR pairs and 0.951 for controls, and the mean log10 length is 3.99 and 3.94.

How we measured similarity and contact.

- We extracted the sequences of all PHR loci and control loci (CHM13 240 loci, HG002 537 loci), built one new pggb graph per genome with the paper's call (`pggb -p 95`) and computed the Jaccard of every pair of loci with `odgi similarity`.
- For contact we used the paper's per-pair function (balanced contact over the two regions divided by the number of valid bin pairs). For each set we computed Spearman rho between Jaccard and contact over all pairs of loci on different chromosomes, with 1,000 bootstrap resamples of loci for the 95% intervals.

A note on the Jaccard values. They are similarities, not distances (`odgi similarity`, column `jaccard.similarity`): the base pairs two loci share in the graph divided by the base pairs of their union. A PHR locus shares a lot with a few partner loci and little or nothing with the rest. For example, the chr1 p PHR (210 kb) has Jaccard 0.519 with the chr6 q PHR and 0.267 with the chr9 q PHR. Across all 835 pairs of CHM13 PHR loci on different chromosomes, 330 share nothing, so the median Jaccard is low (0.0060 among the 683 pairs lost at MAPQ >= 2, 0.0063 among the 152 kept) even though partner pairs are highly similar. Over all 57,360 ordered pairs of loci in the CHM13 graph (PHR and control loci), the 90th percentile is 0 and the 99th is 0.213.

Results at 50 kb unless stated.

| Matrix | PHR rho | Control rho | PHR minus control (95% CI) |
|---|---|---|---|
| Hi-C HG002, all reads | 0.488 | 0.493 | -0.004 (-0.119 to 0.118) |
| Hi-C HG002, MAPQ >= 1 | 0.485 | 0.417 | 0.068 (-0.048 to 0.188) |
| Hi-C HG002, MAPQ >= 2 | 0.359 | 0.250 | 0.110 (-0.007 to 0.224) |
| Hi-C HG002, MAPQ >= 2, 20 kb | 0.362 | 0.232 | 0.130 (0.003 to 0.248) |
| Hi-C CHM13, all reads | 0.469 | 0.539 | -0.070 (-0.247 to 0.120) |
| Hi-C CHM13, MAPQ >= 1 | 0.446 | 0.514 | -0.068 (-0.245 to 0.124) |
| Hi-C CHM13, MAPQ >= 2 | 0.138 (152 pairs) | 0.416 | -0.279 (-0.625 to 0.018) |
| Hi-C CHM13, MAPQ >= 2, 20 kb | -0.009 (478 pairs) | 0.393 | -0.402 (-0.578 to -0.214) |
| Pore-C HG002 (unique reads) | 0.420 | 0.202 | 0.219 (0.119 to 0.311) |
| CiFi HG002 (unique reads) | 0.191 | 0.071 | 0.120 (0.032 to 0.213) |
| CiFi HG002, multi-mapping reads kept | 0.335 | 0.180 | 0.155 (0.038 to 0.273) |

Many pairs of loci share no sequence in the graph (Jaccard 0). For example, 2,120 of 3,514 PHR pairs and 14,499 of 43,334 control pairs share sequence in HG002 Pore-C. The overall rho mixes two questions, whether two loci share sequence and how much they share. Restricting to pairs that share sequence (Jaccard > 0) asks only the second one, whether more similar loci are closer.

| Matrix, Jaccard > 0 only | PHR rho | Control rho |
|---|---|---|
| Pore-C HG002 | 0.365 | -0.005 |
| CiFi HG002 | 0.168 | -0.056 |
| CiFi HG002, multi-mapping reads kept | 0.436 | -0.053 |
| Hi-C HG002 | 0.622 | 0.486 |

- In Hi-C with all reads or MAPQ >= 1, PHRs behave like any segmental duplication (differences within the interval).
- In HG002 Hi-C with unique reads (MAPQ >= 2), PHRs lie above the controls (significant at 20 kb). In CHM13 with unique reads, PHRs lie below the controls, but most CHM13 PHR bins are masked there (section 3) and only 152 or 478 pairs remain.
- In Pore-C and CiFi, PHRs lie above the controls. Among pairs that share sequence, more similar PHRs are closer, and more similar control duplications are not.

## 5. Confounders (R3-M3)

We repeated the Mantel test controlling for arm type, acrocentric p arms, X/Y arms, PHR length, chromosome size and contact propensity (partial Mantel). It stays at 0.624 for HG002 Hi-C, 0.676 for CHM13 Hi-C, 0.393 for Pore-C and 0.281 for CiFi, all with p 0.0001. It stays between 0.26 and 0.81 within same-type pairs, within p-q pairs and without acrocentric p and X/Y arms (`scripts/hic/contact_controls_r3m3.*`).

## 6. Null simulations

We simulated reads that all come from one chromosome end, so every contact between two different ends is a mapping error, and we ran them through the paper's pipelines.

| Assay | Reads | Errors | Depth |
|---|---|---|---|
| Hi-C | pairs of 250 bp (HG002) or 150 bp (CHM13), contact distance P(s) ~ 1/s from 1 to 500 kb | 0.2% substitutions | 40 pairs per kb |
| Pore-C | real NlaIII fragments, 4 per read on average | 1% (0.5% substitutions, 0.25% insertions, 0.25% deletions) | 20 reads per kb |
| CiFi | real DpnII fragments, 4 per read on average | 0.1% | 20 reads per kb |

The read lengths match the real data. HG002 Hi-C reads are 248-249 bp on average (maximum 251) in all 8 input files, and CHM13 Hi-C reads are 150 bp in all 4 files.

Fraction of contact sides placed in the PHR of another arm:

| Simulation | All reads | Unique reads only |
|---|---|---|
| Hi-C HG002 | 19.4% | 0.094% (MAPQ >= 2) |
| Hi-C CHM13 | 14.7% | 0.012% (MAPQ >= 2) |
| Pore-C | 15.9% (`--min-mapq 0`) | 0.24% (`--min-mapq 1`, as in the paper) |
| CiFi | 11.4% (`--min-mapq 0`) | 0.17% (`--min-mapq 1`, as in the paper) |

In HG002 Hi-C, the error rate of each pair of ends follows their similarity (Spearman rho 0.774 over 2,830 pairs).

The paper's Mantel on the simulated matrices (50 kb):

| Assay | Null simulation | Real data |
|---|---|---|
| Hi-C HG002 | 0.755 | 0.657 |
| Hi-C CHM13 | 0.781 | 0.656 |
| Pore-C, unique reads | 0.580 | 0.486 |
| CiFi, unique reads | 0.334 | 0.308 |

The null Mantel is as high as the real one even with 0.24% errors. In the simulation those errors are the only contacts between ends, and they all follow similarity. The Mantel test uses only the ranking, not the number of contacts, so on its own it cannot separate errors from biology. We therefore counted.

## 7. Calibration on depth

For each end we measured how much deeper the real data are than the simulation (real over simulated raw contacts inside its PHR window). We scaled the simulated false contacts to that depth and obtained the false contacts expected in the real data, with a lower and an upper bound (the error can come from either end of a pair). We compared them with the contacts observed between the same PHRs, on different chromosomes. The last two columns give the arm-level Spearman rho between similarity and contact density, before and after removing the upper bound of expected false contacts.

| Assay | Observed contacts between PHRs | Explained by errors | rho before | rho after |
|---|---|---|---|---|
| Hi-C HG002 | 339,856 | 30.0-39.4% | 0.758 | 0.721 |
| Hi-C HG002, MAPQ >= 1 | 166,352 | 35.1-45.0% | 0.732 | 0.697 |
| Hi-C HG002, MAPQ >= 2 | 6,426 | 0.65-1.04% | 0.325 | 0.323 |
| Hi-C CHM13 | 209,864 | 49.8-75.3% | 0.780 | 0.584 |
| Hi-C CHM13, MAPQ >= 1 | 160,340 | 56.2-84.4% | 0.765 | 0.504 |
| Hi-C CHM13, MAPQ >= 2 | 3,327 | 0.86-1.42% | 0.438 | 0.433 |
| Pore-C (paper) | 8,039 | 3.2-8.1% | 0.450 | 0.440 |
| CiFi (paper) | 942 | 2.8-7.9% | 0.308 | 0.300 |
| CiFi, multi-mapping reads kept | 86,218 | 34.8-45.0% | 0.545 | 0.517 |

## 8. The paper's Mantel on matrices without the expected false contacts

We subtracted the upper bound of expected false contacts from each real matrix, pixel by pixel, at the resolution of each statistic. We then balanced it as in the paper and ran the paper's script. Two controls passed.

- With the simulation in place of the real data, the subtraction removes all 488,882 false contacts and leaves none between different ends.
- With nothing subtracted, the pipeline reproduces the paper's CiFi tables byte for byte at 50 and 20 kb.

| Assay | Mantel 50 kb, original | Mantel 50 kb, cleaned | Mantel 20 kb, original | Mantel 20 kb, cleaned | Expected false contacts removed (50 kb) |
|---|---|---|---|---|---|
| Hi-C HG002 | 0.657 | 0.632 | 0.648 | 0.587 | 91.0% of 290,009 |
| Hi-C HG002, MAPQ >= 1 | 0.630 | 0.575 | 0.615 | 0.560 | 89.7% of 163,591 |
| Hi-C CHM13 | 0.656 | 0.563 | 0.721 | 0.684 | 85.7% of 237,749 |
| Hi-C CHM13, MAPQ >= 1 | 0.636 | 0.430 | 0.717 | 0.663 | 84.0% of 203,436 |
| Pore-C | 0.486 | 0.482 | 0.476 | 0.475 | 52.8% of 2,218 |
| CiFi | 0.308 | 0.308 | 0.321 | 0.321 | 30.6% of 125 |
| CiFi, multi-mapping reads kept | 0.440 | 0.427 | 0.512 | 0.502 | 89.7% of 116,735 |

All Mantel tests keep p < 0.0001 (10,000 permutations). The part of the expected false contacts that cannot be removed falls on pixels where the real matrix has fewer contacts than expected.

## 9. CiFi rerun with multi-mapping reads

We reprocessed CiFi from raw reads with two `parse2` branches on the same alignments.

- The paper's branch (`--min-mapq 1`) reproduces the published matrix exactly (18,028,298 of 18,028,298 pixels, 35,560,326 contacts) and the published statistics byte for byte.
- The `--min-mapq 0` branch keeps 240,300,440 pairs (6.76 times more). The Mantel rises from 0.308 to 0.440 at 50 kb and stays at 0.427 after removing the expected false contacts.

The Pore-C rerun is running (expected 2026-11-29 around 00:00 MST).

## 10. What the checks say together

1. Mapping errors are real and large in Hi-C. In the calibration they explain up to 39% of the contacts between PHRs in HG002 and up to 75% in CHM13 (150 bp reads, against 250 bp in HG002). This is why, in Hi-C, PHRs behave like the control duplications.
2. After removing the expected errors, the correlation remains. The cleaned Mantel is 0.43-0.68 in Hi-C. In Pore-C and CiFi the errors are at most 8% and the Mantel barely changes.
3. With unique reads only (Pore-C, CiFi), PHRs lie above matched non-subtelomeric duplications. Among pairs that share sequence, more similar PHRs are closer, and more similar control duplications are not.
4. A high null Mantel only shows that the Mantel alone is not evidence. The evidence is the counting (calibration) and the cleaned matrices.
5. Strict filters remove almost every read inside PHRs (MAPQ >= 2, or unique reads in CiFi), so the tests lose power. CHM13 at MAPQ >= 2 shows this. Keeping multi-mapping reads and removing the expected errors gives more contacts and a stronger test (CiFi 0.427 against 0.308).

## Open issue: real Hi-C reads have more errors than the simulated ones

We mapped 200,000 real R1 reads and the simulated reads with the paper's Bowtie2 options and index.

| | Mismatches per aligned base | Soft-clipped read bases | Reads with 0 mismatches | Error implied by base qualities |
|---|---|---|---|---|
| HG002 real | 2.42% | 30.3% | 36.4% | 1.18% |
| HG002 simulated | 0.20% | 0.01% | 60.8% | |
| CHM13 real | 1.49% | 11.5% | 46.5% | 0.81% |
| CHM13 simulated | 0.19% | 0.02% | 75.1% | |

Real reads carry 7.8 (CHM13) and 12.1 (HG002) times more mismatches, and many cross a ligation junction, which local alignment soft-clips (30% of read bases in HG002). Both effects shorten the unique part of a read and make wrong placements more likely, so the Hi-C simulations probably underestimate the errors, and sections 6-8 for Hi-C may need a rerun with an error profile taken from the real reads (base qualities and junction clipping). We are running the same check of the Pore-C and CiFi error rates (1% and 0.1%) against real reads.
