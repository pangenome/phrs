# Figure 3 — Three-dimensional nuclear organisation mirrors sequence communities

**Sequence-defined communities are physical.** Conventions: **B/W** =
between/within contacts (Hi-C); **W/B** = within/between 3D distance
(Dip-C / sperm scHi-C); both < 1 when within-community arms cluster
(`SURVEY_07 §5 #11`).
(**a**) HG002 Pore-C inter-arm contact matrix, 50 kb; 77 arm-haplotypes
ordered by sequence community; diagonal blocks enriched.
**B/W = 0.056, p = 3.9 × 10⁻⁸⁵**
(`analysis/human/community_based/50000bp/hg002_porec_global_test.tsv`).
(**b**) Forest plot — 14 tests across Hi-C, Pore-C, CiFi, Dip-C, sperm scHi-C
and mouse meiotic Hi-C; every effect lies left of unity. Range 0.020–0.93
(rows label ratio, p, convention; sources in `sources.tsv`).
(**c**) Negative control — per-cell C-community W/B vs **S_all** (7
non-sharing arms): 16/16 GM12878 + 20/20 sperm C-cells have W/B < 1; 0/16 +
1/20 S_all cells do — S_all is 11 % / 40 % *farther* in 3D
(`community_enrichment_16cells_500kb_per_{cell,community_per_cell}.tsv`).
(**d**) Flanking paradox — paired PHR-50 kb vs flanking-100 kb B/W; flanking
is stronger though it has no duplicated sequence (HG002 Hi-C: PHR 0.027 →
flanking 0.0031 ≈ 9×; `SURVEY_07 §1.2` reports 0.002 / 13× by alternate
aggregation). Inset: GM12878 Dip-C flanking particles are more interior than
non-flanking terminal (radial 0.504 vs 0.556; p = 1.6 × 10⁻³⁵;
`PHR_III/dipc_flanking_radial.tsv`). HG002 CiFi flanking: NA.
