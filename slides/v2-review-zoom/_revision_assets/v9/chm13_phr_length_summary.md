# CHM13 PHR Length Summary

Question: total CHM13 PHR length excluding acrocentric short arms and X/Y PAR-terminal PHRs.

Primary source:

- `chm13.phrs.bed`
- 37 raw CHM13 PHR intervals
- No chrY intervals are present in this file.
- No chr13p interval is present in this file.

Exclusion rule:

- Exclude acrocentric p arms: chr13p, chr14p, chr15p, chr21p, chr22p.
- Exclude sex-chromosome PAR-terminal arms: chrXp, chrXq, chrYp, chrYq.

Result from raw PHR calls:

| set | intervals | bp | kbp | Mbp |
| --- | ---: | ---: | ---: | ---: |
| all raw CHM13 PHRs | 37 | 6,014,981 | 6,014.981 | 6.014981 |
| excluded acro p + X/Y PAR-terminal | 6 | 2,589,999 | 2,589.999 | 2.589999 |
| retained non-acro, non-PAR PHRs | 31 | 3,424,982 | 3,424.982 | 3.424982 |

Companion figure-analysis footprint:

- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/chm13_subtelomeric_regions.bed`
- 38 intervals, rounded/expanded to 50 kb analysis bins.
- This file includes chr13p and rounds several raw intervals to analysis-bin boundaries.

| set | intervals | bp | kbp | Mbp |
| --- | ---: | ---: | ---: | ---: |
| all 50 kb analysis intervals | 38 | 6,840,000 | 6,840.000 | 6.840000 |
| excluded acro p + X/Y PAR-terminal | 7 | 3,330,000 | 3,330.000 | 3.330000 |
| retained non-acro, non-PAR analysis footprint | 31 | 3,510,000 | 3,510.000 | 3.510000 |

One-line talk number:

The raw CHM13 PHR calls sum to **3.425 Mbp** after excluding acrocentric short arms and X/Y PAR-terminal PHRs. The corresponding 50 kb binned analysis footprint is **3.510 Mbp**.
