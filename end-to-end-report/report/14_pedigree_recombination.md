# Pedigree subtelomeric inter-chromosomal exchanges

Analysis of inter-chromosomal patches in odgi untangle results (`nth.best=1`) for two pedigrees, validated against the Leiden arm-level community structure from the HPRCv2 analysis.

## Pedigree quality summary

| Pedigree | Samples | Assembly | Within-community patches | Reliability |
|---|---|---|---|---|
| **WashU** | 4 (PAN010 grandmother, PAN011 grandfather, PAN027 mother, PAN028 granddaughter) | **T2T (Cechova et al. 2025)** | **494 / 538 (92%)** | **High — primary evidence** |
| CEPH1463 hifiasm | 28 (Porubsky et al. 2025) | Fragmented contigs (~780/sample) | 324 / 2,775 (12%) | Low — supplementary |
| CEPH1463 verkko | 14 (NA12877-NA12887, NA12889-NA12892) | Fragmented contigs (~144/sample) | 359 / 2,671 (13%) | Low — supplementary |

**Key principle**: Only inter-chromosomal patches whose source and query arms are in the SAME Leiden community (from `hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`) are considered biologically credible. The Leiden community structure was derived from 233 HPRCv2 samples and represents the reference truth for which subtelomeric arm pairs share sequence.

**WashU is the primary dataset** because of its T2T quality. CEPH1463 results are reported only where BOTH hifiasm AND verkko detect the same parent feature (cross-assembler validation), to control for assembly fragmentation artifacts.

## Method

- **Patch**: contiguous run of untangle segments where the child's flank maps to the same parent chromosome + haplotype.
- **HQ filter**: `is_interchr=True` + `min_score >= 0.8` + `500 bp <= size <= 100 kb`.
- **Quality** = `min_score` (lowest alignment identity within the patch; 1.0 = perfect).
- **Pattern classification** (from immediate left/right neighbor):

| Pattern | Definition | Interpretation |
|---|---|---|
| `gene_conversion_like` | Sandwich: `chrN:hX → chrM:hY → chrN:hX` (same chr + same hap on both sides, patch from different chr + different hap) | Ectopic gene conversion tract |
| `crossover_like` | `chrN:h1 → chrM:hZ → chrN:h2` (same chr but different hap left vs right) | Meiotic crossover with inter-chr sequence at the breakpoint |
| `acros_like` | Patch in a flank with ≥5 inter-chr patches from ≥3 different source chromosomes | Extensive NAHR signature (acrocentric-like) |

- **Leiden community validation**: each patch is cross-referenced against the 15 arm-level Leiden communities from HPRCv2. Only **within-community** patches are reported below.

Full per-patch TSV: `PHR_III/pedigrees/all_pedigrees_patches.tsv` (all 5,984 HQ patches with community columns; cross-community and unknown patches are filtered out from this report but available in the TSV).

# Part 1: WashU pedigree (T2T quality)

WashU has 4 T2T-resolved samples (PAN010, PAN011, PAN027, PAN028) sequenced and assembled to telomere-to-telomere quality. **The WashU results are the primary evidence in this report.**

WashU HQ inter-chr patches: 538 total, **494 (92%) within-community** (validated by Leiden community structure from HPRCv2 analysis).

Pattern breakdown of WashU within-community patches:
- `acros_like`: 229
- `gene_conversion_like`: 133
- `sandwich_same_hap`: 115
- `crossover_like`: 16
- `complex`: 1

## WashU `gene_conversion_like` (all 133 within-community patches)


| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | PAN027 maternal (hap1) vs PAN010 | chr13p | 418,382-427,257 | 8,875 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 2 | PAN027 maternal (hap1) vs PAN010 | chr21p | 392,614-401,099 | 8,485 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 3 | PAN027 maternal (hap1) vs PAN010 | chr21p | 236,804-243,127 | 6,323 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 4 | PAN027 maternal (hap1) vs PAN010 | chr13p | 449,858-454,600 | 4,742 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 5 | PAN027 maternal (hap1) vs PAN010 | chr13p | 226,088-229,801 | 3,713 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 6 | PAN027 maternal (hap1) vs PAN010 | chr13p | 85,116-88,723 | 3,607 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 7 | PAN027 maternal (hap1) vs PAN010 | chr13p | 98,877-102,414 | 3,537 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 8 | PAN027 maternal (hap1) vs PAN010 | chr13p | 77,223-79,926 | 2,703 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 9 | PAN027 maternal (hap1) vs PAN010 | chr13p | 379,044-381,738 | 2,694 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 10 | PAN027 maternal (hap1) vs PAN010 | chr13p | 66,165-68,754 | 2,589 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 11 | PAN028 maternal (hap1) vs PAN027 | chr13p | 430,970-433,553 | 2,583 | chr22p:h1 | chr13p:h2 | chr13p:h2 | 1.000/1.000 | gene_conv | in | C7 |
| 12 | PAN027 maternal (hap1) vs PAN010 | chr21p | 253,177-255,718 | 2,541 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 13 | PAN027 maternal (hap1) vs PAN010 | chr13p | 105,322-107,786 | 2,464 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 14 | PAN027 maternal (hap1) vs PAN010 | chr13p | 125,962-128,121 | 2,159 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 15 | PAN027 maternal (hap1) vs PAN010 | chr13p | 144,698-146,857 | 2,159 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 16 | PAN027 maternal (hap1) vs PAN010 | chr13p | 207,214-209,373 | 2,159 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 17 | PAN027 maternal (hap1) vs PAN010 | chr13p | 257,558-259,717 | 2,159 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 18 | PAN027 maternal (hap1) vs PAN010 | chr13p | 395,776-397,935 | 2,159 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 19 | PAN027 maternal (hap1) vs PAN010 | chr13p | 414,665-416,824 | 2,159 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 20 | PAN027 maternal (hap1) vs PAN010 | chr13p | 446,145-448,304 | 2,159 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 21 | PAN027 maternal (hap1) vs PAN010 | chr13p | 477,625-479,784 | 2,159 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 22 | PAN027 maternal (hap1) vs PAN010 | chr21p | 182,975-185,134 | 2,159 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 23 | PAN027 maternal (hap1) vs PAN010 | chr21p | 201,864-204,023 | 2,159 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 24 | PAN027 paternal (hap2) vs PAN011 | chr19p | 74,936-76,974 | 2,038 | chr7p:h2 | chr19p:h1 | chr19p:h1 | 1.000/1.000 | gene_conv | in | C3 |
| 25 | PAN027 maternal (hap1) vs PAN010 | chr3q | 329,643-331,678 | 2,035 | chr9q:h2 | chr3q:h1 | chr3q:h1 | 1.000/1.000 | gene_conv | in | C3 |
| 26 | PAN027 paternal (hap2) vs PAN011 | chrXp | 329,647-331,663 | 2,016 | chrYp:h2 | chrXp:h1 | chrXp:h1 | 1.000/1.000 | gene_conv | out | C15 |
| 27 | PAN027 maternal (hap1) vs PAN010 | chr13p | 45,125-46,983 | 1,858 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 28 | PAN027 maternal (hap1) vs PAN010 | chr21p | 45,361-47,219 | 1,858 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 29 | PAN027 maternal (hap1) vs PAN010 | chr22p | 102,417-104,189 | 1,772 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 1.000/1.000 | gene_conv | in | C7 |
| 30 | PAN027 maternal (hap1) vs PAN010 | chr22p | 108,858-110,626 | 1,768 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 1.000/1.000 | gene_conv | in | C7 |
| 31 | PAN027 maternal (hap1) vs PAN010 | chr21p | 328,865-330,630 | 1,765 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 32 | PAN027 maternal (hap1) vs PAN010 | chr21p | 292,132-293,862 | 1,730 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 33 | PAN027 maternal (hap1) vs PAN010 | chr22p | 222,553-224,111 | 1,558 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 1.000/1.000 | gene_conv | in | C7 |
| 34 | PAN027 maternal (hap1) vs PAN010 | chr13p | 190,495-192,049 | 1,554 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 35 | PAN027 maternal (hap1) vs PAN010 | chr22p | 241,446-243,000 | 1,554 | chr15p:h1 | chr22p:h2 | chr22p:h2 | 1.000/1.000 | gene_conv | in | C7 |
| 36 | PAN027 paternal (hap2) vs PAN011 | chr19p | 102,352-103,840 | 1,488 | chr7p:h2 | chr19p:h1 | chr19p:h1 | 1.000/1.000 | gene_conv | in | C3 |
| 37 | PAN027 maternal (hap1) vs PAN010 | chr22p | 193,760-195,207 | 1,447 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 1.000/1.000 | gene_conv | in | C7 |
| 38 | PAN027 maternal (hap1) vs PAN010 | chr22p | 225,247-226,690 | 1,443 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 1.000/1.000 | gene_conv | in | C7 |
| 39 | PAN027 maternal (hap1) vs PAN010 | chr3q | 366,707-368,133 | 1,426 | chr9q:h2 | chr3q:h1 | chr3q:h1 | 1.000/1.000 | gene_conv | in | C3 |
| 40 | PAN027 maternal (hap1) vs PAN010 | chr21p | 65,054-66,397 | 1,343 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 41 | PAN027 paternal (hap2) vs PAN011 | chr3q | 422,686-423,908 | 1,222 | chr9q:h2 | chr3q:h1 | chr3q:h1 | 1.000/1.000 | gene_conv | in | C3 |
| 42 | PAN027 maternal (hap1) vs PAN010 | chr13p | 82,399-83,545 | 1,146 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 43 | PAN027 maternal (hap1) vs PAN010 | chr13p | 321,225-322,371 | 1,146 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 44 | PAN027 maternal (hap1) vs PAN010 | chr13p | 340,097-341,243 | 1,146 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 45 | PAN027 maternal (hap1) vs PAN010 | chr13p | 440,862-442,008 | 1,146 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 46 | PAN027 maternal (hap1) vs PAN010 | chr21p | 76,323-77,469 | 1,146 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 47 | PAN027 maternal (hap1) vs PAN010 | chr22p | 158,440-159,586 | 1,146 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 1.000/1.000 | gene_conv | in | C7 |
| 48 | PAN027 paternal (hap2) vs PAN011 | chr21p | 234,468-235,614 | 1,146 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 49 | PAN027 paternal (hap2) vs PAN011 | chr21p | 259,602-260,748 | 1,146 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 50 | PAN028 maternal (hap1) vs PAN027 | chr13p | 195,646-196,792 | 1,146 | chr15p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 51 | PAN028 maternal (hap1) vs PAN027 | chr13p | 428,270-429,416 | 1,146 | chr22p:h1 | chr13p:h2 | chr13p:h2 | 1.000/1.000 | gene_conv | in | C7 |
| 52 | PAN027 maternal (hap1) vs PAN010 | chr13p | 91,440-92,576 | 1,136 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 53 | PAN027 maternal (hap1) vs PAN010 | chr13p | 223,505-224,641 | 1,136 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 54 | PAN027 maternal (hap1) vs PAN010 | chr13p | 236,091-237,227 | 1,136 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 55 | PAN027 maternal (hap1) vs PAN010 | chr13p | 267,563-268,699 | 1,136 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 56 | PAN027 maternal (hap1) vs PAN010 | chr13p | 443,562-444,698 | 1,136 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 57 | PAN027 maternal (hap1) vs PAN010 | chr13p | 468,750-469,886 | 1,136 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 58 | PAN027 maternal (hap1) vs PAN010 | chr21p | 110,884-112,020 | 1,136 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 59 | PAN027 maternal (hap1) vs PAN010 | chr21p | 186,688-187,824 | 1,136 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 60 | PAN027 maternal (hap1) vs PAN010 | chr21p | 192,985-194,121 | 1,136 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 61 | PAN027 maternal (hap1) vs PAN010 | chr21p | 199,281-200,417 | 1,136 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 62 | PAN027 maternal (hap1) vs PAN010 | chr21p | 205,578-206,714 | 1,136 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 63 | PAN027 maternal (hap1) vs PAN010 | chr22p | 148,546-149,682 | 1,136 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 1.000/1.000 | gene_conv | in | C7 |
| 64 | PAN027 paternal (hap2) vs PAN011 | chr21p | 243,431-244,567 | 1,136 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 65 | PAN027 maternal (hap1) vs PAN010 | chr21p | 363,661-364,684 | 1,023 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 66 | PAN027 maternal (hap1) vs PAN010 | chr13p | 75,064-76,077 | 1,013 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 67 | PAN027 maternal (hap1) vs PAN010 | chr13p | 94,044-95,057 | 1,013 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 68 | PAN027 maternal (hap1) vs PAN010 | chr21p | 87,962-88,975 | 1,013 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 69 | PAN027 maternal (hap1) vs PAN010 | chr21p | 113,467-114,480 | 1,013 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 70 | PAN027 maternal (hap1) vs PAN010 | chr21p | 195,568-196,581 | 1,013 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 71 | PAN027 maternal (hap1) vs PAN010 | chr21p | 333,024-334,037 | 1,013 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 72 | PAN028 maternal (hap1) vs PAN027 | chr13p | 75,065-76,078 | 1,013 | chr15p:h2 | chr13p:h1 | chr13p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 73 | PAN027 maternal (hap1) vs PAN010 | chr21p | 100,586-101,597 | 1,011 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 74 | PAN027 maternal (hap1) vs PAN010 | chr21p | 107,024-108,035 | 1,011 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 1.000/1.000 | gene_conv | in | C7 |
| 75 | PAN027 paternal (hap2) vs PAN011 | chrXp | 189,599-194,850 | 5,251 | chrYp:h2 | chrXp:h1 | chrXp:h1 | 1.000/1.000 | gene_conv | out | C15 |
| 76 | PAN027 maternal (hap1) vs PAN010 | chr13p | 332,792-337,633 | 4,841 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/0.999 | gene_conv | in | C7 |
| 77 | PAN027 maternal (hap1) vs PAN010 | chr13p | 347,551-350,245 | 2,694 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.999/0.999 | gene_conv | in | C7 |
| 78 | PAN027 maternal (hap1) vs PAN010 | chr22p | 153,286-154,844 | 1,558 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.999/0.999 | gene_conv | in | C7 |
| 79 | PAN027 maternal (hap1) vs PAN010 | chr13p | 365,303-371,603 | 6,300 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 1.000/0.999 | gene_conv | in | C7 |
| 80 | PAN027 maternal (hap1) vs PAN010 | chr22p | 111,762-113,222 | 1,460 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.999/0.999 | gene_conv | in | C7 |
| 81 | PAN027 maternal (hap1) vs PAN010 | chr13p | 327,517-331,345 | 3,828 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.999/0.998 | gene_conv | in | C7 |
| 82 | PAN027 maternal (hap1) vs PAN010 | chr21p | 189,271-191,430 | 2,159 | chr22p:h2 | chr21p:h1 | chr21p:h1 | 0.999/0.998 | gene_conv | in | C7 |
| 83 | PAN027 maternal (hap1) vs PAN010 | chr13p | 261,271-262,407 | 1,136 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.998/0.998 | gene_conv | in | C7 |
| 84 | PAN027 maternal (hap1) vs PAN010 | chr13p | 188,336-189,349 | 1,013 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.998/0.998 | gene_conv | in | C7 |
| 85 | PAN027 maternal (hap1) vs PAN010 | chr13p | 194,632-195,645 | 1,013 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.998/0.998 | gene_conv | in | C7 |
| 86 | PAN027 maternal (hap1) vs PAN010 | chr13p | 213,506-214,519 | 1,013 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.998/0.998 | gene_conv | in | C7 |
| 87 | PAN028 maternal (hap1) vs PAN027 | chr9q | 370,238-379,107 | 8,869 | chr3q:h1 | chr9q:h2 | chr9q:h2 | 0.999/0.998 | gene_conv | in | C3 |
| 88 | PAN027 maternal (hap1) vs PAN010 | chr22p | 92,583-94,047 | 1,464 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.997/0.997 | gene_conv | in | C7 |
| 89 | PAN027 maternal (hap1) vs PAN010 | chr13p | 361,707-364,290 | 2,583 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.999/0.997 | gene_conv | in | C7 |
| 90 | PAN027 maternal (hap1) vs PAN010 | chr13p | 343,945-346,405 | 2,460 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.999/0.997 | gene_conv | in | C7 |
| 91 | PAN027 maternal (hap1) vs PAN010 | chr22p | 137,104-138,549 | 1,445 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.997/0.997 | gene_conv | in | C7 |
| 92 | PAN028 maternal (hap1) vs PAN027 | chr22p | 24,901-40,349 | 15,448 | chr14p:h1 | chr22p:h2 | chr22p:h2 | 0.997/0.997 | gene_conv | in | C7 |
| 93 | PAN028 maternal (hap1) vs PAN027 | chr13p | 452,445-454,604 | 2,159 | chr15p:h2 | chr13p:h1 | chr13p:h1 | 0.998/0.997 | gene_conv | in | C7 |
| 94 | PAN028 maternal (hap1) vs PAN027 | chr15p | 278,997-280,555 | 1,558 | chr22p:h1 | chr15p:h2 | chr15p:h2 | 0.996/0.996 | gene_conv | in | C7 |
| 95 | PAN027 maternal (hap1) vs PAN010 | chr13p | 351,692-360,149 | 8,457 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.999/0.996 | gene_conv | in | C7 |
| 96 | PAN027 maternal (hap1) vs PAN010 | chr22p | 120,680-123,372 | 2,692 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.998/0.996 | gene_conv | in | C7 |
| 97 | PAN027 maternal (hap1) vs PAN010 | chr22p | 207,803-208,816 | 1,013 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.996/0.996 | gene_conv | in | C7 |
| 98 | PAN027 maternal (hap1) vs PAN010 | chr22p | 174,860-178,466 | 3,606 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.998/0.996 | gene_conv | in | C7 |
| 99 | PAN027 paternal (hap2) vs PAN011 | chr13p | 212,436-213,572 | 1,136 | chr21p:h1 | chr13p:h2 | chr13p:h2 | 0.995/0.995 | gene_conv | in | C7 |
| 100 | PAN027 maternal (hap1) vs PAN010 | chr22p | 124,508-125,955 | 1,447 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.995/0.995 | gene_conv | in | C7 |
| 101 | PAN027 maternal (hap1) vs PAN010 | chr22p | 244,136-245,579 | 1,443 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.993/0.993 | gene_conv | in | C7 |
| 102 | PAN027 paternal (hap2) vs PAN011 | chr13p | 130,139-132,728 | 2,589 | chr15p:h2 | chr13p:h1 | chr13p:h1 | 0.995/0.993 | gene_conv | in | C7 |
| 103 | PAN027 maternal (hap1) vs PAN010 | chr13p | 129,667-130,803 | 1,136 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.993/0.993 | gene_conv | in | C7 |
| 104 | PAN027 maternal (hap1) vs PAN010 | chr22p | 130,808-132,255 | 1,447 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.990/0.990 | gene_conv | in | C7 |
| 105 | PAN027 maternal (hap1) vs PAN010 | chr22p | 155,980-157,427 | 1,447 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.990/0.990 | gene_conv | in | C7 |
| 106 | PAN028 maternal (hap1) vs PAN027 | chr22q | 490,400-491,408 | 1,008 | chr21q:h2 | chr22q:h1 | chr22q:h1 | 0.990/0.990 | gene_conv | in | C6 |
| 107 | PAN027 maternal (hap1) vs PAN010 | chr22p | 126,968-129,672 | 2,704 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.993/0.990 | gene_conv | in | C7 |
| 108 | PAN027 maternal (hap1) vs PAN010 | chr13p | 252,271-253,417 | 1,146 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.990/0.990 | gene_conv | in | C7 |
| 109 | PAN027 maternal (hap1) vs PAN010 | chr13p | 292,737-293,873 | 1,136 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.990/0.990 | gene_conv | in | C7 |
| 110 | PAN027 maternal (hap1) vs PAN010 | chr13p | 374,303-375,439 | 1,136 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.990/0.990 | gene_conv | in | C7 |
| 111 | PAN027 maternal (hap1) vs PAN010 | chr13p | 323,921-326,504 | 2,583 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.994/0.988 | gene_conv | in | C7 |
| 112 | PAN027 maternal (hap1) vs PAN010 | chr13p | 286,437-287,573 | 1,136 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.988/0.988 | gene_conv | in | C7 |
| 113 | PAN027 maternal (hap1) vs PAN010 | chr22p | 218,951-220,394 | 1,443 | chr21p:h1 | chr22p:h2 | chr22p:h2 | 0.988/0.988 | gene_conv | in | C7 |
| 114 | PAN027 maternal (hap1) vs PAN010 | chr13p | 184,195-185,753 | 1,558 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.985/0.985 | gene_conv | in | C7 |
| 115 | PAN027 maternal (hap1) vs PAN010 | chr13p | 238,674-248,675 | 10,001 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.998/0.984 | gene_conv | in | C7 |
| 116 | PAN027 maternal (hap1) vs PAN010 | chr22p | 200,058-201,505 | 1,447 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.981/0.981 | gene_conv | in | C7 |
| 117 | PAN028 maternal (hap1) vs PAN027 | chr14p | 153,300-154,316 | 1,016 | chr15p:h1 | chr14p:h2 | chr14p:h2 | 0.977/0.977 | gene_conv | in | C7 |
| 118 | PAN027 maternal (hap1) vs PAN010 | chr22p | 341,147-344,371 | 3,224 | chr21p:h1 | chr22p:h2 | chr22p:h2 | 0.975/0.975 | gene_conv | in | C7 |
| 119 | PAN027 maternal (hap1) vs PAN010 | chr22p | 236,703-237,838 | 1,135 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.971/0.971 | gene_conv | in | C7 |
| 120 | PAN027 paternal (hap2) vs PAN011 | chr15p | 232,996-234,545 | 1,549 | chr22p:h1 | chr15p:h2 | chr15p:h2 | 0.966/0.966 | gene_conv | in | C7 |
| 121 | PAN028 maternal (hap1) vs PAN027 | chr10q | 402,748-403,967 | 1,219 | chr4q:h2 | chr10q:h1 | chr10q:h1 | 0.957/0.957 | gene_conv | out | C1 |
| 122 | PAN027 maternal (hap1) vs PAN010 | chr13p | 230,937-234,541 | 3,604 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.985/0.955 | gene_conv | in | C7 |
| 123 | PAN027 maternal (hap1) vs PAN010 | chr22p | 162,272-163,721 | 1,449 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.954/0.954 | gene_conv | in | C7 |
| 124 | PAN028 maternal (hap1) vs PAN027 | chr5q | 394,731-395,732 | 1,001 | chr6q:h1 | chr5q:h2 | chr5q:h2 | 0.951/0.951 | gene_conv | in | C11 |
| 125 | PAN027 maternal (hap1) vs PAN010 | chr22p | 181,148-182,599 | 1,451 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.950/0.950 | gene_conv | in | C7 |
| 126 | PAN027 maternal (hap1) vs PAN010 | chr13p | 218,351-221,955 | 3,604 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.983/0.947 | gene_conv | in | C7 |
| 127 | PAN027 maternal (hap1) vs PAN010 | chr13p | 161,987-163,433 | 1,446 | chr22p:h2 | chr13p:h1 | chr13p:h1 | 0.940/0.940 | gene_conv | in | C7 |
| 128 | PAN028 maternal (hap1) vs PAN027 | chr5q | 400,738-402,742 | 2,004 | chr6q:h1 | chr5q:h2 | chr5q:h2 | 0.940/0.928 | gene_conv | in | C11 |
| 129 | PAN027 maternal (hap1) vs PAN010 | chr22p | 231,543-232,990 | 1,447 | chr13p:h1 | chr22p:h2 | chr22p:h2 | 0.904/0.904 | gene_conv | in | C7 |
| 130 | PAN028 maternal (hap1) vs PAN027 | chr1p | 35,250-36,605 | 1,355 | chr8p:h1 | chr1p:h2 | chr1p:h2 | 0.863/0.863 | gene_conv | in | C11 |
| 131 | PAN027 maternal (hap1) vs PAN010 | chr18p | 1,008-5,255 | 4,247 | chr10p:h2 | chr18p:h1 | chr18p:h1 | 0.843/0.843 | gene_conv | in | C2 |
| 132 | PAN028 maternal (hap1) vs PAN027 | chr22p | 394,906-395,939 | 1,033 | chr21p:h1 | chr22p:h2 | chr22p:h2 | 0.817/0.817 | gene_conv | in | C7 |
| 133 | PAN028 maternal (hap1) vs PAN027 | chr22p | 475,642-476,652 | 1,010 | chr21p:h1 | chr22p:h2 | chr22p:h2 | 0.814/0.814 | gene_conv | out | C7 |

## WashU `crossover_like` (all 16 within-community patches)


| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | PAN027 paternal (hap2) vs PAN011 | chr3q | 345,337-348,751 | 3,414 | chr9q:h2 | chr3q:h2 | chr3q:h1 | 1.000/1.000 | crossover | in | C3 |
| 2 | PAN028 maternal (hap1) vs PAN027 | chr13p | 126,976-128,122 | 1,146 | chr15p:h2 | chr13p:h2 | chr13p:h1 | 1.000/1.000 | crossover | in | C7 |
| 3 | PAN028 maternal (hap1) vs PAN027 | chr13p | 145,712-146,858 | 1,146 | chr15p:h2 | chr13p:h2 | chr13p:h1 | 1.000/1.000 | crossover | in | C7 |
| 4 | PAN028 maternal (hap1) vs PAN027 | chr13p | 296,334-297,480 | 1,146 | chr15p:h2 | chr13p:h2 | chr13p:h1 | 1.000/1.000 | crossover | in | C7 |
| 5 | PAN028 maternal (hap1) vs PAN027 | chr13p | 421,974-423,120 | 1,146 | chr15p:h2 | chr13p:h1 | chr13p:h2 | 1.000/1.000 | crossover | in | C7 |
| 6 | PAN028 maternal (hap1) vs PAN027 | chr13p | 459,750-460,896 | 1,146 | chr22p:h1 | chr13p:h2 | chr13p:h1 | 1.000/1.000 | crossover | in | C7 |
| 7 | PAN027 paternal (hap2) vs PAN011 | chr13p | 126,430-128,589 | 2,159 | chr15p:h2 | chr13p:h2 | chr13p:h1 | 0.998/0.998 | crossover | in | C7 |
| 8 | PAN028 maternal (hap1) vs PAN027 | chr22p | 47,449-65,278 | 17,829 | chr14p:h1 | chr22p:h2 | chr22p:h1 | 0.998/0.998 | crossover | in | C7 |
| 9 | PAN028 maternal (hap1) vs PAN027 | chr1p | 5,588-20,462 | 14,874 | chr8p:h1 | chr1p:h1 | chr1p:h2 | 0.999/0.997 | crossover | in | C11 |
| 10 | PAN028 maternal (hap1) vs PAN027 | chr13p | 284,892-287,574 | 2,682 | chr21p:h1 | chr13p:h2 | chr13p:h1 | 0.998/0.997 | crossover | in | C7 |
| 11 | PAN028 maternal (hap1) vs PAN027 | chr13p | 293,874-295,321 | 1,447 | chr22p:h1 | chr13p:h1 | chr13p:h2 | 0.995/0.995 | crossover | in | C7 |
| 12 | PAN028 maternal (hap1) vs PAN027 | chr13p | 143,252-144,699 | 1,447 | chr21p:h1 | chr13p:h1 | chr13p:h2 | 0.977/0.977 | crossover | in | C7 |
| 13 | PAN028 maternal (hap1) vs PAN027 | chr3q | 366,221-368,310 | 2,089 | chr9q:h2 | chr3q:h2 | chr3q:h1 | 0.970/0.970 | crossover | in | C3 |
| 14 | PAN028 maternal (hap1) vs PAN027 | chr3q | 262,953-290,922 | 27,969 | chr7p:h2 | chr3q:h2 | chr3q:h1 | 0.961/0.943 | crossover | in | C3 |
| 15 | PAN028 maternal (hap1) vs PAN027 | chr5q | 392,728-393,730 | 1,002 | chr8p:h2 | chr5q:h1 | chr5q:h2 | 0.943/0.943 | crossover | in | C11 |
| 16 | PAN027 maternal (hap1) vs PAN010 | chr13p | 176,753-179,454 | 2,701 | chr22p:h2 | chr13p:h2 | chr13p:h1 | 0.960/0.920 | crossover | in | C7 |

## WashU `acros_like` (229 within-community patches, top 30 by quality)


| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | PAN028 maternal (hap1) vs PAN027 | chr13p | 0-23,406 | 23,406 | chr22p:h1 | edge | chr15p:h2 | 1.000/1.000 | acros_like | in | C7 |
| 2 | PAN028 maternal (hap1) vs PAN027 | chr13p | 24,431-39,881 | 15,450 | chr22p:h1 | chr15p:h2 | chr13p:h1 | 1.000/1.000 | acros_like | in | C7 |
| 3 | PAN028 maternal (hap1) vs PAN027 | chr13p | 481,083-492,096 | 11,013 | chr22p:h1 | chr13p:h1 | chr21p:h1 | 1.000/1.000 | acros_like | in | C7 |
| 4 | PAN028 maternal (hap1) vs PAN027 | chr15p | 16,187-22,351 | 6,164 | chr13p:h2 | chr22p:h1 | chr15p:h2 | 1.000/1.000 | acros_like | in | C7 |
| 5 | PAN028 maternal (hap1) vs PAN027 | chr13p | 492,096-497,456 | 5,360 | chr21p:h1 | chr22p:h1 | chr22p:h1 | 1.000/1.000 | acros_like | in | C7 |
| 6 | PAN027 maternal (hap1) vs PAN010 | chr13p | 428,270-433,549 | 5,279 | chr22p:h2 | chr13p:h1 | chr15p:h1 | 1.000/1.000 | acros_like | in | C7 |
| 7 | PAN028 maternal (hap1) vs PAN027 | chr13p | 409,382-414,661 | 5,279 | chr22p:h1 | chr13p:h1 | chr15p:h2 | 1.000/1.000 | acros_like | in | C7 |
| 8 | PAN028 maternal (hap1) vs PAN027 | chr13p | 404,236-408,369 | 4,133 | chr22p:h1 | chr15p:h2 | chr13p:h1 | 1.000/1.000 | acros_like | in | C7 |
| 9 | PAN028 maternal (hap1) vs PAN027 | chr13p | 101,346-104,187 | 2,841 | chr22p:h1 | chr15p:h2 | chr13p:h1 | 1.000/1.000 | acros_like | in | C7 |
| 10 | PAN028 maternal (hap1) vs PAN027 | chr13p | 497,456-499,997 | 2,541 | chr22p:h1 | chr21p:h1 | chr21p:h2 | 1.000/1.000 | acros_like | in | C7 |
| 11 | PAN028 maternal (hap1) vs PAN027 | chr13p | 98,878-101,346 | 2,468 | chr15p:h2 | chr13p:h1 | chr22p:h1 | 1.000/1.000 | acros_like | in | C7 |
| 12 | PAN028 maternal (hap1) vs PAN027 | chr15p | 212,433-214,893 | 2,460 | chr13p:h2 | chr22p:h1 | chr15p:h2 | 1.000/1.000 | acros_like | in | C7 |
| 13 | PAN027 maternal (hap1) vs PAN010 | chr13p | 389,484-391,643 | 2,159 | chr22p:h2 | chr13p:h1 | chr15p:h1 | 1.000/1.000 | acros_like | in | C7 |
| 14 | PAN027 maternal (hap1) vs PAN010 | chr21p | 164,091-166,250 | 2,159 | chr22p:h2 | chr21p:h1 | chr13p:h1 | 1.000/1.000 | acros_like | in | C7 |
| 15 | PAN027 maternal (hap1) vs PAN010 | chr21p | 214,459-216,618 | 2,159 | chr22p:h2 | chr13p:h1 | chr21p:h1 | 1.000/1.000 | acros_like | in | C7 |
| 16 | PAN028 maternal (hap1) vs PAN027 | chr9q | 329,410-331,445 | 2,035 | chr3q:h1 | chr16q:h2 | chr9q:h2 | 1.000/1.000 | acros_like | in | C3 |
| 17 | PAN028 maternal (hap1) vs PAN027 | chr13p | 397,936-399,494 | 1,558 | chr22p:h1 | chr15p:h2 | chr13p:h1 | 1.000/1.000 | acros_like | in | C7 |
| 18 | PAN027 maternal (hap1) vs PAN010 | chr13p | 391,643-393,197 | 1,554 | chr15p:h1 | chr22p:h2 | chr22p:h2 | 1.000/1.000 | acros_like | in | C7 |
| 19 | PAN027 maternal (hap1) vs PAN010 | chr13p | 410,532-412,086 | 1,554 | chr15p:h1 | chr22p:h2 | chr22p:h2 | 1.000/1.000 | acros_like | in | C7 |
| 20 | PAN028 maternal (hap1) vs PAN027 | chr15p | 216,039-217,593 | 1,554 | chr22p:h1 | chr15p:h2 | chr13p:h1 | 1.000/1.000 | acros_like | in | C7 |
| 21 | PAN028 maternal (hap1) vs PAN027 | chr15p | 234,927-236,477 | 1,550 | chr22p:h1 | chr13p:h2 | chr13p:h2 | 1.000/1.000 | acros_like | in | C7 |
| 22 | PAN028 maternal (hap1) vs PAN027 | chr13p | 334,952-336,498 | 1,546 | chr22p:h1 | chr21p:h1 | chr13p:h1 | 1.000/1.000 | acros_like | in | C7 |
| 23 | PAN028 maternal (hap1) vs PAN027 | chr13p | 105,323-106,776 | 1,453 | chr22p:h1 | chr13p:h1 | chr15p:h2 | 1.000/1.000 | acros_like | in | C7 |
| 24 | PAN027 maternal (hap1) vs PAN010 | chr13p | 457,294-458,741 | 1,447 | chr15p:h1 | chr22p:h2 | chr22p:h2 | 1.000/1.000 | acros_like | in | C7 |
| 25 | PAN027 maternal (hap1) vs PAN010 | chr13p | 463,590-465,037 | 1,447 | chr15p:h1 | chr22p:h2 | chr22p:h2 | 1.000/1.000 | acros_like | in | C7 |
| 26 | PAN028 maternal (hap1) vs PAN027 | chr13p | 438,406-439,853 | 1,447 | chr15p:h2 | chr22p:h1 | chr13p:h2 | 1.000/1.000 | acros_like | in | C7 |
| 27 | PAN028 maternal (hap1) vs PAN027 | chr13p | 444,702-446,149 | 1,447 | chr15p:h2 | chr22p:h1 | chr13p:h2 | 1.000/1.000 | acros_like | in | C7 |
| 28 | PAN028 maternal (hap1) vs PAN027 | chr9q | 366,750-368,176 | 1,426 | chr3q:h2 | chr3q:h1 | chr9q:h2 | 1.000/1.000 | acros_like | in | C3 |
| 29 | PAN028 maternal (hap1) vs PAN027 | chr9q | 328,203-329,410 | 1,207 | chr16q:h2 | chr9q:h1 | chr3q:h1 | 1.000/1.000 | acros_like | in | C3 |
| 30 | PAN027 maternal (hap1) vs PAN010 | chr13p | 271,159-272,305 | 1,146 | chr21p:h1 | chr22p:h2 | chr15p:h1 | 1.000/1.000 | acros_like | in | C7 |


# Part 2: CEPH1463 pedigree (supplementary, cross-assembler validated only)

CEPH1463 has 28 (hifiasm) / 14 (verkko) samples assembled from fragmented contigs. **The assemblies are NOT T2T**, so untangle results are noisier:
- CEPH1463 hifiasm: 2,775 HQ patches, only 12% within-community
- CEPH1463 verkko: 2,671 HQ patches, only 13% within-community

To control for assembly artifacts, only inter-chromosomal exchanges detected by **BOTH** hifiasm AND verkko assemblies are reported here. These are robust because they replicate across two independent assembly methods.

**Cross-assembler validation criterion**: same parent + same chromosome pair detected by both assemblers in at least one child each (within Leiden community).

## CEPH1463 cross-assembler validated parent features (11 total)

| Parent | Chr pair | Leiden community | Hifi children | Verk children | Hifi best score | Verk best score |
|---|---|---|---|---|---|---|
| NA12877 | chr1/chr19 | C6 | NA12879,NA12883 | NA12884 | 0.982 | 0.860 |
| NA12877 | chr10/chr18 | C2 | NA12883,NA12884 | NA12884,NA12885 | 0.998 | 0.997 |
| NA12877 | chr17/chr19 | C6 | NA12883 | NA12884 | 0.984 | 0.868 |
| NA12877 | chr6/chr9 | C5 | NA12886 | NA12881 | 0.975 | 0.989 |
| NA12878 | chr10/chr18 | C2 | NA12884,NA12885,NA12887 | NA12882,NA12886 | 0.972 | 0.998 |
| NA12878 | chr19/chr22 | C6 | NA12881,NA12882,NA12886 | NA12879,NA12881,NA12883,NA12887 | 0.897 | 0.978 |
| NA12878 | chr21/chr22 | C7 | NA12879 | NA12881,NA12882 | 0.992 | 0.996 |
| NA12878 | chr6/chr9 | C5 | NA12881,NA12883,NA12884,NA12885 | NA12882 | 0.947 | 0.959 |
| NA12889 | chr12/chr9 | C5 | NA12877 | NA12877 | 0.957 | 0.978 |
| NA12890 | chr12/chr9 | C5 | NA12877 | NA12877 | 0.971 | 0.976 |
| NA12892 | chr21/chr22 | C6 | NA12878 | NA12878 | 0.994 | 0.993 |

### Detailed patches per validated CEPH1463 parent feature


#### NA12877: chr1/chr19 (C6)

**Hifiasm:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12879 paternal (hap1) vs NA12877 | chr1q | 495,288-496,344 | 1,056 | chr19q:h2 | chr19q:h1 | chr8q:h1 | 0.982/0.982 | acros_like | in | C6 |
| 2 | NA12883 paternal (hap1) vs NA12877 | chr1q | 495,058-496,114 | 1,056 | chr19q:h2 | chr19q:h1 | chr8q:h1 | 0.982/0.982 | acros_like | in | C6 |
| 3 | NA12879 paternal (hap1) vs NA12877 | chr1q | 481,449-483,542 | 2,093 | chr19q:h2 | chr19q:h1 | chr4p:h2 | 0.872/0.872 | acros_like | in | C6 |
| 4 | NA12883 paternal (hap1) vs NA12877 | chr1q | 481,218-483,311 | 2,093 | chr19q:h2 | chr19q:h1 | chr4p:h2 | 0.872/0.872 | acros_like | in | C6 |
| 5 | NA12883 paternal (hap1) vs NA12877 | chr1q | 489,225-491,274 | 2,049 | chr19q:h2 | chr13p:h2 | chr19q:h1 | 0.872/0.872 | acros_like | in | C6 |
| 6 | NA12879 paternal (hap1) vs NA12877 | chr1q | 489,456-491,504 | 2,048 | chr19q:h2 | chr13p:h2 | chr19q:h1 | 0.871/0.871 | acros_like | in | C6 |
| 7 | NA12879 paternal (hap1) vs NA12877 | chr1q | 491,504-495,288 | 3,784 | chr19q:h1 | chr19q:h2 | chr19q:h2 | 0.916/0.853 | acros_like | in | C6 |
| 8 | NA12883 paternal (hap1) vs NA12877 | chr1q | 491,274-495,058 | 3,784 | chr19q:h1 | chr19q:h2 | chr19q:h2 | 0.916/0.853 | acros_like | in | C6 |

**Verkko:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12884 paternal (hap2) vs NA12877 | chr19q | 251,649-252,781 | 1,132 | chr1q:h2 | chr17q:h2 | chr19p:h2 | 0.860/0.860 | acros_like | in | C6 |

#### NA12877: chr10/chr18 (C2)

**Hifiasm:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12883 paternal (hap1) vs NA12877 | chr18p | 28,047-29,172 | 1,125 | chr10p:h2 | chr18p:h2 | chr18p:h2 | 0.998/0.998 | sandwich_same_hap | in | C2 |
| 2 | NA12884 paternal (hap1) vs NA12877 | chr18p | 28,047-29,172 | 1,125 | chr10p:h2 | chr18p:h2 | chr18p:h2 | 0.998/0.998 | sandwich_same_hap | in | C2 |
| 3 | NA12883 paternal (hap1) vs NA12877 | chr18p | 0-1,003 | 1,003 | chr10p:h2 | edge | chr18p:h2 | 0.871/0.871 | complex | in | C2 |
| 4 | NA12884 paternal (hap1) vs NA12877 | chr18p | 0-1,003 | 1,003 | chr10p:h2 | edge | chr18p:h2 | 0.871/0.871 | complex | in | C2 |

**Verkko:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12884 paternal (hap2) vs NA12877 | chr10p | 24,510-25,548 | 1,038 | chr18p:h2 | chr10p:h1 | chr10p:h1 | 0.997/0.997 | gene_conv | in | C2 |
| 2 | NA12885 paternal (hap2) vs NA12877 | chr10p | 24,537-25,575 | 1,038 | chr18p:h2 | chr10p:h1 | chr10p:h1 | 0.997/0.997 | gene_conv | in | C2 |

#### NA12877: chr17/chr19 (C6)

**Hifiasm:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12883 paternal (hap1) vs NA12877 | chr17q | 488,911-490,712 | 1,801 | chr19q:h2 | chr19q:h1 | chr6q:h2 | 0.984/0.984 | acros_like | in | C6 |
| 2 | NA12883 paternal (hap1) vs NA12877 | chr17q | 495,321-496,378 | 1,057 | chr19q:h2 | chr6q:h2 | chr8q:h1 | 0.978/0.978 | acros_like | in | C6 |
| 3 | NA12883 paternal (hap1) vs NA12877 | chr17q | 486,255-488,911 | 2,656 | chr19q:h1 | chr6q:h2 | chr19q:h2 | 0.968/0.953 | acros_like | in | C6 |

**Verkko:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12884 paternal (hap2) vs NA12877 | chr19q | 248,747-251,649 | 2,902 | chr17q:h2 | chr9q:h2 | chr1q:h2 | 0.868/0.868 | acros_like | in | C6 |

#### NA12877: chr6/chr9 (C5)

**Hifiasm:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12886 paternal (hap1) vs NA12877 | chr9p | 4,577-9,146 | 4,569 | chr6p:h1 | chr12p:h1 | chr3p:h2 | 0.975/0.975 | acros_like | in | C5 |

**Verkko:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12881 paternal (hap2) vs NA12877 | chr9p | 2,197-4,462 | 2,265 | chr6p:h1 | chr7p:h1 | chr16p:h1 | 0.989/0.989 | acros_like | in | C5 |

#### NA12878: chr10/chr18 (C2)

**Hifiasm:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12884 maternal (hap2) vs NA12878 | chr10p | 15,934-38,247 | 22,313 | chr18p:h1 | chr10q:h2 | chr10q:h2 | 0.988/0.972 | gene_conv | in | C2 |
| 2 | NA12885 maternal (hap2) vs NA12878 | chr10p | 16,679-38,992 | 22,313 | chr18p:h1 | chr10q:h2 | chr10q:h2 | 0.988/0.972 | gene_conv | in | C2 |
| 3 | NA12887 maternal (hap2) vs NA12878 | chr10p | 66,491-67,506 | 1,015 | chr18p:h1 | chr10q:h1 | chr10q:h1 | 0.912/0.912 | sandwich_same_hap | in | C2 |
| 4 | NA12884 maternal (hap2) vs NA12878 | chr10p | 65,990-67,003 | 1,013 | chr18p:h1 | chr10q:h2 | chr10q:h2 | 0.912/0.912 | gene_conv | in | C2 |
| 5 | NA12885 maternal (hap2) vs NA12878 | chr10p | 66,735-67,748 | 1,013 | chr18p:h1 | chr10q:h2 | chr10q:h2 | 0.912/0.912 | gene_conv | in | C2 |
| 6 | NA12887 maternal (hap2) vs NA12878 | chr10p | 0-38,723 | 38,723 | chr18p:h1 | edge | chr10q:h1 | 0.976/0.853 | complex | in | C2 |

**Verkko:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12882 maternal (hap1) vs NA12878 | chr10p | 35,121-36,127 | 1,006 | chr18p:h2 | chr10p:h1 | chr10p:h1 | 0.998/0.998 | gene_conv | in | C2 |
| 2 | NA12886 maternal (hap1) vs NA12878 | chr10p | 35,121-36,127 | 1,006 | chr18p:h2 | chr10p:h1 | chr10p:h1 | 0.998/0.998 | gene_conv | in | C2 |
| 3 | NA12882 maternal (hap1) vs NA12878 | chr10p | 0-17,855 | 17,855 | chr18p:h2 | edge | chr10p:h1 | 0.921/0.856 | complex | in | C2 |
| 4 | NA12886 maternal (hap1) vs NA12878 | chr10p | 0-17,855 | 17,855 | chr18p:h2 | edge | chr10p:h1 | 0.921/0.856 | complex | in | C2 |

#### NA12878: chr19/chr22 (C6)

**Hifiasm:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12882 maternal (hap2) vs NA12878 | chr19q | 478,035-499,219 | 21,184 | chr22q:h2 | chr19p:h2 | chr19p:h2 | 0.961/0.897 | sandwich_same_hap | in | C6 |
| 2 | NA12881 maternal (hap2) vs NA12878 | chr19q | 478,064-499,219 | 21,155 | chr22q:h2 | chr19p:h2 | chr19p:h2 | 0.959/0.897 | sandwich_same_hap | in | C6 |
| 3 | NA12886 maternal (hap2) vs NA12878 | chr19q | 479,029-493,319 | 14,290 | chr22q:h2 | chr19p:h2 | chr20q:h1 | 0.950/0.850 | complex | in | C6 |

**Verkko:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12879 maternal (hap1) vs NA12878 | chr19q | 486,790-487,852 | 1,062 | chr22q:h1 | chr21q:h2 | chr21q:h2 | 0.978/0.978 | acros_like | in | C6 |
| 2 | NA12883 maternal (hap1) vs NA12878 | chr19q | 486,810-487,872 | 1,062 | chr22q:h1 | chr21q:h2 | chr21q:h2 | 0.978/0.978 | acros_like | in | C6 |
| 3 | NA12881 maternal (hap1) vs NA12878 | chr19q | 479,330-482,450 | 3,120 | chr22q:h1 | chr8q:h2 | chr21q:h2 | 0.969/0.958 | acros_like | in | C6 |
| 4 | NA12887 maternal (hap1) vs NA12878 | chr19q | 479,299-482,419 | 3,120 | chr22q:h1 | chr8q:h2 | chr21q:h2 | 0.969/0.958 | acros_like | in | C6 |
| 5 | NA12879 maternal (hap1) vs NA12878 | chr19q | 480,301-483,422 | 3,121 | chr22q:h1 | chr8q:h2 | chr10q:h1 | 0.969/0.957 | acros_like | in | C6 |
| 6 | NA12883 maternal (hap1) vs NA12878 | chr19q | 480,321-483,442 | 3,121 | chr22q:h1 | chr8q:h2 | chr10q:h1 | 0.969/0.957 | acros_like | in | C6 |

#### NA12878: chr21/chr22 (C7)

**Hifiasm:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12879 maternal (hap2) vs NA12878 | chr21p | 4,042-16,908 | 12,866 | chr22p:h2 | chr4p:h2 | chr4p:h2 | 0.994/0.992 | complex | in | C7 |

**Verkko:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12882 maternal (hap1) vs NA12878 | chr22p | 20,546-26,353 | 5,807 | chr21p:h2 | chr22p:h2 | chr22p:h1 | 0.996/0.996 | crossover | in | C7 |
| 2 | NA12882 maternal (hap1) vs NA12878 | chr22p | 245,634-316,725 | 71,091 | chr21p:h2 | chr22p:h1 | chr22p:h2 | 0.996/0.995 | crossover | in | C7 |
| 3 | NA12881 maternal (hap1) vs NA12878 | chr21p | 127,336-129,919 | 2,583 | chr22p:h1 | chr21p:h2 | chr21p:h2 | 0.995/0.995 | gene_conv | in | C7 |
| 4 | NA12882 maternal (hap1) vs NA12878 | chr22p | 3,347-17,185 | 13,838 | chr21p:h2 | chr22p:h2 | chr22p:h2 | 0.995/0.994 | sandwich_same_hap | in | C7 |
| 5 | NA12881 maternal (hap1) vs NA12878 | chr21p | 111,594-114,177 | 2,583 | chr22p:h1 | chr21p:h2 | chr21p:h2 | 0.994/0.994 | gene_conv | in | C7 |
| 6 | NA12881 maternal (hap1) vs NA12878 | chr21p | 119,465-122,048 | 2,583 | chr22p:h1 | chr21p:h2 | chr21p:h2 | 0.994/0.994 | gene_conv | in | C7 |
| 7 | NA12881 maternal (hap1) vs NA12878 | chr21p | 135,207-137,790 | 2,583 | chr22p:h1 | chr21p:h2 | chr21p:h2 | 0.994/0.994 | gene_conv | in | C7 |
| 8 | NA12881 maternal (hap1) vs NA12878 | chr21p | 150,957-153,540 | 2,583 | chr22p:h1 | chr21p:h2 | chr21p:h2 | 0.994/0.994 | gene_conv | in | C7 |
| 9 | NA12882 maternal (hap1) vs NA12878 | chr22p | 329,959-341,682 | 11,723 | chr21p:h2 | chr22p:h2 | chr22p:h2 | 0.996/0.993 | sandwich_same_hap | in | C7 |
| 10 | NA12881 maternal (hap1) vs NA12878 | chr21p | 4,103-7,851 | 3,748 | chr22p:h2 | chr4p:h1 | chr21p:h2 | 0.993/0.993 | acros_like | in | C7 |
| 11 | NA12881 maternal (hap1) vs NA12878 | chr21p | 44,057-48,091 | 4,034 | chr22p:h2 | chr22p:h1 | chr21p:h2 | 0.995/0.992 | acros_like | in | C7 |
| 12 | NA12881 maternal (hap1) vs NA12878 | chr21p | 30,860-44,057 | 13,197 | chr22p:h1 | chr21p:h2 | chr22p:h2 | 0.991/0.991 | acros_like | in | C7 |
| 13 | NA12881 maternal (hap1) vs NA12878 | chr21p | 143,078-145,669 | 2,591 | chr22p:h1 | chr21p:h2 | chr21p:h2 | 0.991/0.991 | gene_conv | in | C7 |
| 14 | NA12881 maternal (hap1) vs NA12878 | chr21p | 85,063-86,212 | 1,149 | chr22p:h2 | chr21p:h2 | chr21p:h2 | 0.991/0.991 | sandwich_same_hap | in | C7 |
| 15 | NA12882 maternal (hap1) vs NA12878 | chr22p | 225,293-228,742 | 3,449 | chr21p:h2 | chr22p:h2 | chr22p:h1 | 0.993/0.989 | crossover | in | C7 |
| 16 | NA12881 maternal (hap1) vs NA12878 | chr21p | 63,771-67,247 | 3,476 | chr22p:h2 | chr21p:h2 | chr21p:h2 | 0.989/0.988 | sandwich_same_hap | in | C7 |
| 17 | NA12881 maternal (hap1) vs NA12878 | chr21p | 232,441-235,020 | 2,579 | chr22p:h1 | chr21p:h2 | chr21p:h2 | 0.988/0.988 | gene_conv | in | C7 |
| 18 | NA12881 maternal (hap1) vs NA12878 | chr21p | 352,711-359,037 | 6,326 | chr22p:h1 | chr21p:h2 | chr14q:h2 | 0.988/0.988 | acros_like | in | C7 |
| 19 | NA12881 maternal (hap1) vs NA12878 | chr21p | 259,192-272,142 | 12,950 | chr22p:h2 | chr21p:h2 | chr14q:h2 | 0.991/0.987 | acros_like | in | C7 |
| 20 | NA12881 maternal (hap1) vs NA12878 | chr21p | 440,089-453,365 | 13,276 | chr22p:h2 | chr14q:h2 | chr21p:h2 | 0.987/0.987 | acros_like | in | C7 |
| 21 | NA12882 maternal (hap1) vs NA12878 | chr22p | 355,745-364,631 | 8,886 | chr21p:h2 | chr22p:h2 | chr22p:h2 | 0.988/0.986 | sandwich_same_hap | in | C7 |
| 22 | NA12881 maternal (hap1) vs NA12878 | chr21p | 225,585-228,164 | 2,579 | chr22p:h1 | chr21p:h2 | chr21p:h2 | 0.985/0.985 | gene_conv | in | C7 |
| 23 | NA12882 maternal (hap1) vs NA12878 | chr22p | 345,639-351,609 | 5,970 | chr21p:h2 | chr22p:h2 | chr22p:h2 | 0.993/0.982 | sandwich_same_hap | in | C7 |
| 24 | NA12881 maternal (hap1) vs NA12878 | chr21p | 70,390-73,563 | 3,173 | chr22p:h2 | chr21p:h2 | chr21p:h2 | 0.981/0.973 | sandwich_same_hap | in | C7 |
| 25 | NA12881 maternal (hap1) vs NA12878 | chr21p | 76,706-79,879 | 3,173 | chr22p:h2 | chr21p:h2 | chr21p:h2 | 0.981/0.973 | sandwich_same_hap | in | C7 |
| 26 | NA12882 maternal (hap1) vs NA12878 | chr22p | 166,610-170,594 | 3,984 | chr21p:h2 | chr22p:h2 | chr22p:h2 | 0.973/0.973 | sandwich_same_hap | in | C7 |
| 27 | NA12881 maternal (hap1) vs NA12878 | chr21p | 91,395-92,549 | 1,154 | chr22p:h1 | chr21p:h2 | chr21p:h2 | 0.972/0.972 | gene_conv | in | C7 |
| 28 | NA12881 maternal (hap1) vs NA12878 | chr21p | 466,263-467,881 | 1,618 | chr22p:h2 | chr22p:h1 | chr14q:h2 | 0.968/0.968 | acros_like | in | C7 |
| 29 | NA12882 maternal (hap1) vs NA12878 | chr22p | 63,755-65,897 | 2,142 | chr21p:h2 | chr22p:h2 | chr22p:h2 | 0.975/0.965 | sandwich_same_hap | in | C7 |
| 30 | NA12882 maternal (hap1) vs NA12878 | chr22p | 398,864-399,867 | 1,003 | chr21p:h2 | chr22p:h1 | chr22p:h2 | 0.965/0.965 | crossover | in | C7 |
| 31 | NA12881 maternal (hap1) vs NA12878 | chr21p | 484,251-488,362 | 4,111 | chr22p:h1 | chr21p:h2 | chr14q:h2 | 0.962/0.962 | acros_like | in | C7 |
| 32 | NA12881 maternal (hap1) vs NA12878 | chr21p | 491,927-499,863 | 7,936 | chr22p:h2 | chr14q:h2 | chr21p:h2 | 0.959/0.959 | acros_like | in | C7 |
| 33 | NA12882 maternal (hap1) vs NA12878 | chr22p | 418,944-419,971 | 1,027 | chr21p:h2 | chr22p:h1 | chr22p:h1 | 0.957/0.957 | gene_conv | in | C7 |
| 34 | NA12881 maternal (hap1) vs NA12878 | chr21p | 297,254-301,276 | 4,022 | chr22p:h2 | chr21p:h2 | chr14q:h2 | 0.954/0.954 | acros_like | in | C7 |
| 35 | NA12882 maternal (hap1) vs NA12878 | chr22p | 388,779-389,780 | 1,001 | chr21p:h2 | chr22p:h1 | chr22p:h1 | 0.944/0.944 | gene_conv | in | C7 |
| 36 | NA12882 maternal (hap1) vs NA12878 | chr22p | 441,160-442,161 | 1,001 | chr21p:h2 | chr22p:h1 | chr22p:h1 | 0.942/0.942 | gene_conv | in | C7 |
| 37 | NA12882 maternal (hap1) vs NA12878 | chr22p | 443,163-445,192 | 2,029 | chr21p:h2 | chr22p:h1 | chr22p:h1 | 0.938/0.933 | gene_conv | in | C7 |
| 38 | NA12882 maternal (hap1) vs NA12878 | chr22p | 432,040-433,041 | 1,001 | chr21p:h2 | chr22p:h1 | chr22p:h1 | 0.932/0.932 | gene_conv | in | C7 |
| 39 | NA12882 maternal (hap1) vs NA12878 | chr22p | 390,781-392,811 | 2,030 | chr21p:h2 | chr22p:h1 | chr22p:h1 | 0.933/0.932 | gene_conv | in | C7 |
| 40 | NA12882 maternal (hap1) vs NA12878 | chr22p | 437,136-440,158 | 3,022 | chr21p:h2 | chr22p:h1 | chr22p:h1 | 0.946/0.923 | gene_conv | in | C7 |
| 41 | NA12882 maternal (hap1) vs NA12878 | chr22p | 383,745-387,777 | 4,032 | chr21p:h2 | chr22p:h2 | chr22p:h1 | 0.916/0.898 | crossover | in | C7 |
| 42 | NA12882 maternal (hap1) vs NA12878 | chr22p | 411,910-413,922 | 2,012 | chr21p:h2 | chr22p:h1 | chr22p:h1 | 0.901/0.895 | gene_conv | in | C7 |
| 43 | NA12881 maternal (hap1) vs NA12878 | chr21p | 462,097-466,263 | 4,166 | chr22p:h1 | chr21p:h2 | chr22p:h2 | 0.932/0.895 | acros_like | in | C7 |
| 44 | NA12882 maternal (hap1) vs NA12878 | chr22p | 427,014-429,026 | 2,012 | chr21p:h2 | chr22p:h1 | chr22p:h1 | 0.890/0.889 | gene_conv | in | C7 |
| 45 | NA12882 maternal (hap1) vs NA12878 | chr22p | 378,718-379,719 | 1,001 | chr21p:h2 | chr22p:h1 | chr22p:h2 | 0.883/0.883 | crossover | in | C7 |
| 46 | NA12882 maternal (hap1) vs NA12878 | chr22p | 380,736-381,743 | 1,007 | chr21p:h2 | chr22p:h2 | chr22p:h2 | 0.880/0.880 | sandwich_same_hap | in | C7 |
| 47 | NA12881 maternal (hap1) vs NA12878 | chr21p | 186,420-190,446 | 4,026 | chr22p:h2 | chr14q:h2 | chr21p:h2 | 0.939/0.879 | acros_like | in | C7 |
| 48 | NA12881 maternal (hap1) vs NA12878 | chr21p | 193,568-197,594 | 4,026 | chr22p:h2 | chr14q:h2 | chr21p:h2 | 0.939/0.879 | acros_like | in | C7 |
| 49 | NA12881 maternal (hap1) vs NA12878 | chr21p | 215,281-219,291 | 4,010 | chr22p:h2 | chr14q:h2 | chr21p:h2 | 0.933/0.875 | acros_like | in | C7 |
| 50 | NA12882 maternal (hap1) vs NA12878 | chr22p | 404,876-407,884 | 3,008 | chr21p:h2 | chr22p:h2 | chr22p:h1 | 0.882/0.875 | crossover | in | C7 |
| 51 | NA12882 maternal (hap1) vs NA12878 | chr22p | 480,544-484,566 | 4,022 | chr21p:h2 | chr22p:h1 | chr22p:h2 | 0.899/0.862 | crossover | in | C7 |
| 52 | NA12882 maternal (hap1) vs NA12878 | chr22p | 80,826-88,868 | 8,042 | chr21p:h2 | chr14q:h2 | chr22p:h1 | 0.860/0.860 | complex | in | C7 |
| 53 | NA12882 maternal (hap1) vs NA12878 | chr22p | 369,662-376,703 | 7,041 | chr21p:h2 | chr22p:h2 | chr22p:h2 | 0.941/0.859 | sandwich_same_hap | in | C7 |
| 54 | NA12882 maternal (hap1) vs NA12878 | chr22p | 446,210-451,273 | 5,063 | chr21p:h2 | chr22p:h1 | chr22p:h2 | 0.921/0.859 | crossover | in | C7 |
| 55 | NA12881 maternal (hap1) vs NA12878 | chr21p | 208,152-212,161 | 4,009 | chr22p:h2 | chr14q:h2 | chr21p:h2 | 0.927/0.857 | acros_like | in | C7 |
| 56 | NA12882 maternal (hap1) vs NA12878 | chr22p | 393,829-397,845 | 4,016 | chr21p:h2 | chr22p:h1 | chr22p:h1 | 0.882/0.823 | gene_conv | in | C7 |
| 57 | NA12882 maternal (hap1) vs NA12878 | chr22p | 74,284-76,575 | 2,291 | chr21p:h2 | chr22p:h2 | chr22p:h1 | 0.841/0.817 | crossover | in | C7 |
| 58 | NA12882 maternal (hap1) vs NA12878 | chr22p | 461,315-462,316 | 1,001 | chr21p:h2 | chr22p:h2 | chr22p:h1 | 0.815/0.815 | crossover | in | C7 |

#### NA12878: chr6/chr9 (C5)

**Hifiasm:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12885 maternal (hap2) vs NA12878 | chr6p | 141,162-142,753 | 1,591 | chr9p:h1 | chr8q:h1 | chr20p:h1 | 0.947/0.947 | acros_like | in | C5 |
| 2 | NA12881 maternal (hap2) vs NA12878 | chr6p | 138,182-139,794 | 1,612 | chr9p:h1 | chr8q:h1 | chr20p:h1 | 0.945/0.945 | acros_like | in | C5 |
| 3 | NA12884 maternal (hap2) vs NA12878 | chr6p | 138,185-139,797 | 1,612 | chr9p:h1 | chr8q:h1 | chr20p:h1 | 0.945/0.945 | acros_like | in | C5 |
| 4 | NA12883 maternal (hap2) vs NA12878 | chr6p | 132,474-134,058 | 1,584 | chr9p:h1 | chr8q:h1 | chr8q:h1 | 0.942/0.942 | acros_like | in | C5 |

**Verkko:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12882 maternal (hap1) vs NA12878 | chr9p | 8,582-11,352 | 2,770 | chr6p:h2 | chr6p:h1 | chr17p:h2 | 0.959/0.959 | acros_like | in | C5 |
| 2 | NA12882 maternal (hap1) vs NA12878 | chr9p | 7,576-8,582 | 1,006 | chr6p:h1 | chr10p:h2 | chr6p:h2 | 0.950/0.950 | acros_like | in | C5 |

#### NA12889: chr12/chr9 (C5)

**Hifiasm:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12877 paternal (hap1) vs NA12889 | chr9p | 29,934-32,892 | 2,958 | chr12p:h2 | chr9q:h2 | chr6q:h2 | 0.969/0.957 | acros_like | in | C5 |

**Verkko:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12877 paternal (hap2) vs NA12889 | chr12p | 31,715-33,943 | 2,228 | chr9p:h2 | chr5p:h2 | chr20p:h2 | 0.978/0.978 | acros_like | in | C5 |
| 2 | NA12877 paternal (hap2) vs NA12889 | chr12p | 20,975-24,181 | 3,206 | chr9p:h2 | chr20p:h2 | chr20q:h1 | 0.966/0.966 | acros_like | in | C5 |

#### NA12890: chr12/chr9 (C5)

**Hifiasm:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12877 maternal (hap2) vs NA12890 | chr9p | 4,022-8,589 | 4,567 | chr12p:h2 | chr20p:h1 | chr20p:h1 | 0.971/0.971 | acros_like | in | C5 |

**Verkko:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12877 maternal (hap1) vs NA12890 | chr12p | 5,749-7,088 | 1,339 | chr9p:h2 | chr15p:h1 | chr15p:h2 | 0.976/0.976 | acros_like | in | C5 |
| 2 | NA12877 maternal (hap1) vs NA12890 | chr12p | 22,296-25,500 | 3,204 | chr9p:h2 | chr6p:h1 | chr12q:h1 | 0.961/0.961 | acros_like | in | C5 |

#### NA12892: chr21/chr22 (C6)

**Hifiasm:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12878 maternal (hap2) vs NA12892 | chr22q | 496,517-497,711 | 1,194 | chr21q:h2 | chr9q:h2 | chr1q:h2 | 0.994/0.994 | acros_like | in | C6 |
| 2 | NA12878 maternal (hap2) vs NA12892 | chr22q | 489,219-491,028 | 1,809 | chr21q:h2 | chr19q:h1 | chr16q:h2 | 0.986/0.986 | acros_like | in | C6 |
| 3 | NA12878 maternal (hap2) vs NA12892 | chr22p | 1,205-15,340 | 14,135 | chr21p:h2 | chr4p:h1 | chr4p:h1 | 0.992/0.985 | complex | in | C7 |

**Verkko:**

| # | Pairing | Flank | Position | Size | Source | L bg | R bg | Score | Pattern | PHR | Leiden |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NA12878 maternal (hap1) vs NA12892 | chr22p | 22,233-41,245 | 19,012 | chr21p:h2 | chr22p:h2 | chr14p:h1 | 0.994/0.993 | acros_like | in | C7 |
| 2 | NA12878 maternal (hap1) vs NA12892 | chr22p | 1,295-11,328 | 10,033 | chr21p:h2 | chr4p:h2 | chr22p:h2 | 0.994/0.990 | acros_like | in | C7 |
| 3 | NA12878 maternal (hap1) vs NA12892 | chr22p | 63,269-64,404 | 1,135 | chr21p:h2 | chr22p:h2 | chr22p:h2 | 0.990/0.990 | sandwich_same_hap | in | C7 |
| 4 | NA12878 maternal (hap1) vs NA12892 | chr22p | 380,262-381,400 | 1,138 | chr21p:h2 | chr14p:h1 | chr14p:h1 | 0.989/0.989 | acros_like | in | C7 |
| 5 | NA12878 maternal (hap1) vs NA12892 | chr22p | 174,628-179,614 | 4,986 | chr21p:h2 | chr22q:h1 | chr22q:h1 | 0.968/0.963 | gene_conv | in | C7 |
| 6 | NA12878 maternal (hap1) vs NA12892 | chr22p | 65,405-115,121 | 49,716 | chr21p:h2 | chr22p:h2 | chr22q:h1 | 0.971/0.894 | crossover | in | C7 |
| 7 | NA12878 maternal (hap1) vs NA12892 | chr22p | 412,233-414,266 | 2,033 | chr21p:h2 | chr14p:h1 | chr14p:h1 | 0.904/0.879 | acros_like | in | C7 |
| 8 | NA12878 maternal (hap1) vs NA12892 | chr22p | 469,807-470,808 | 1,001 | chr21p:h2 | chr14p:h1 | chr14p:h1 | 0.871/0.871 | acros_like | in | C7 |
| 9 | NA12878 maternal (hap1) vs NA12892 | chr22p | 440,479-441,480 | 1,001 | chr21p:h2 | chr14p:h1 | chr14p:h1 | 0.862/0.862 | acros_like | in | C7 |
| 10 | NA12878 maternal (hap1) vs NA12892 | chr22p | 437,443-438,446 | 1,003 | chr21p:h2 | chr14p:h1 | chr14p:h1 | 0.830/0.830 | acros_like | in | C7 |
| 11 | NA12878 maternal (hap1) vs NA12892 | chr22p | 465,786-468,806 | 3,020 | chr21p:h2 | chr14p:h1 | chr14p:h1 | 0.852/0.802 | acros_like | in | C7 |

# Conclusions

1. **WashU pedigree (T2T) provides strong evidence for ectopic gene conversion in subtelomeric PHRs**: 133 `gene_conversion_like` patches at perfect alignment scores (1.000/1.000), 92% within Leiden communities, predominantly involving:
   - Acrocentric p-arms (community C7: chr13p/14p/15p/21p/22p) — ~120 patches, expected from extensive rDNA/satellite NAHR
   - **Non-acrocentric exchanges**: chr3q/chr9q (C3, f7501 cluster), chr18p/chr10p (C2, Linardopoulou pair), chr8p/chr11p (C3 region), chr19p/chr7p (C3), chrXp/chrYp (C15, PAR1)

2. **Inheritance across generations is captured**: PAN027 maternal patches (from PAN010) and PAN027 paternal patches (from PAN011) are independently detected. PAN028 maternal patches show transmission from PAN027 with patches both retained and expanded between generations.

3. **CEPH1463 cross-assembler validation is sparse but informative**: only 11 unique parent features are detected by both hifiasm AND verkko within Leiden communities. The most robust are:
   - **chr10/chr18 (C2)** in NA12877 paternal AND NA12878 maternal — independent observations of the same Linardopoulou pair, reinforcing the WashU finding
   - **chr19/chr22 (C6)** in NA12878 maternal — transmitted to multiple children in both assemblers
   - **chr12/chr9 (C5)** in BOTH NA12889 (paternal grandfather) AND NA12890 (paternal grandmother) of NA12877 — community C5 (RPL23A/WASH/DDX11L) signature visible in both G1 individuals
   - **chr6/chr9 (C5)** in NA12877 paternal AND NA12878 maternal — same community, independent

4. **The CEPH1463 single-assembler results are dominated by graph topology noise** (12-13% within-community vs 92% in WashU) and should not be used as primary evidence. Use only the cross-assembler-validated parent features.

5. **All findings reinforce the HPRCv2 Leiden community structure**: every cross-assembler validated CEPH1463 finding maps to a known community (C2, C3, C5, C6, C7), and the WashU patches predominantly fall within communities. The pedigree-level analysis directly observes the inter-chromosomal exchange events that produce the population-level community structure.

## Methods

- **Untangle**: `odgi untangle -e 50000 -m 1000 -j 0 -n 100`, `nth.best=1` only.
- **Score**: column 7 of untangle BED (alignment identity, 0-1).
- **Pattern classification**: based on immediate left/right neighbor patches within the same flank.
- **Leiden community validation**: cross-referenced against `PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv` (15 communities, 41 arms).
- **Cross-assembler validation (CEPH1463)**: same parent + same chromosome pair detected by both hifiasm and verkko in at least one child each.
- **Scripts**: `scripts/pedigree/analyze-pedigree-recombination.py`, `scripts/pedigree/plot-pedigree-untangle.R`.
