# Fig5 whole-genome existing-PAF IMPG-like scan

## Scope

This scan reuses existing whole-genome WFMASH p95 and SweepGA/FastGA f32 PAFs.
It does not build a seqwish/ODGI graph. Query alignments are tiled into 2 kb
query-space bins, summarized by target chromosome/arm, and then aggregated into
10 kb and 50 kb plotting bins.

## Execution

- Slurm job id: `1706571`
- Host: `octopus07`
- `SLURM_CPUS_PER_TASK`: `48`
- Process workers: `12`
- `pigz` threads per worker: `4`
- Accounted helper threads: `48`
- Wall seconds: ``

The worker count and decompression threads are derived from `SLURM_CPUS_PER_TASK`.
On the octopus run this accounts for the full 48-CPU allocation across worker
processes and `pigz` helper threads.

## Key focal comparisons

| region | method | evidence | comparison | query | target | bins | aligned bp | weighted identity | match distance |
|---|---|---|---|---|---|---:|---:|---:|---:|
| PAR_XY | sweepga_fastga_f32 | filtered_one_to_one | PAN027pat_vs_PAN011_joint | chrX | chrY | 84 | 164,035 | 0.999672 | 0.000328 |
| PAR_XY | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chrX | chrY | 4,984 | 10,898,982 | 0.936414 | 0.063586 |
| PAR_XY | wfmash_p95_updated_bin | filtered_one_to_one | PAN027pat_vs_PAN011_joint | chrX | chrY | 75 | 147,001 | 0.999823 | 0.000177 |
| PAR_XY | wfmash_p95_updated_bin | raw_many_many | PAN027pat_vs_PAN011_joint | chrX | chrY | 106 | 225,476 | 0.991728 | 0.008272 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN027mat_vs_PAN010_joint | chr22 | chr15 | 1 | 2 | 1.000000 | 0.000000 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN027mat_vs_PAN010_joint | chr13 | chr15 | 1 | 2,000 | 0.999500 | 0.000500 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN027mat_vs_PAN010_joint | chr21 | chr13 | 1 | 2,000 | 1.000000 | 0.000000 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN027mat_vs_PAN010_joint | chr14 | chr22 | 2 | 4,000 | 1.000000 | 0.000000 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN027pat_vs_PAN011_joint | chr22 | chr21 | 2 | 4,000 | 1.000000 | 0.000000 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr22 | chr15 | 6 | 12,000 | 0.999583 | 0.000417 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr21 | chr15 | 8 | 10,161 | 0.999311 | 0.000689 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr14 | chr15 | 1 | 1,035 | 1.000000 | 0.000000 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr22 | chr14 | 14 | 26,192 | 0.997557 | 0.002443 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr21 | chr14 | 18 | 34,001 | 0.998735 | 0.001265 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr22 | chr13 | 1 | 2,000 | 1.000000 | 0.000000 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr21 | chr13 | 1 | 1,284 | 1.000000 | 0.000000 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr21 | chr22 | 5 | 6,052 | 0.999670 | 0.000330 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr15 | chr13 | 55 | 110,000 | 0.999355 | 0.000645 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr14 | chr22 | 1 | 2,000 | 0.984508 | 0.015492 |
| acrocentric_p_cross | sweepga_fastga_f32 | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr15 | chr21 | 5 | 8,014 | 0.999626 | 0.000374 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr22 | chr21 | 193 | 2,512,596 | 0.986994 | 0.013006 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr15 | chr21 | 168 | 814,041 | 0.951076 | 0.048924 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr14 | chr21 | 227 | 1,632,123 | 0.959530 | 0.040470 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr13 | chr21 | 241 | 4,499,224 | 0.986234 | 0.013766 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr22 | chr15 | 214 | 2,245,078 | 0.981754 | 0.018246 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr21 | chr15 | 243 | 1,888,989 | 0.980176 | 0.019824 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr14 | chr15 | 164 | 533,619 | 0.975261 | 0.024739 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr13 | chr15 | 234 | 2,988,456 | 0.984202 | 0.015798 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr22 | chr14 | 129 | 505,615 | 0.982881 | 0.017119 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr21 | chr14 | 191 | 649,085 | 0.963584 | 0.036416 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr15 | chr14 | 172 | 1,286,058 | 0.957551 | 0.042449 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr13 | chr14 | 51 | 234,123 | 0.986156 | 0.013844 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr22 | chr13 | 179 | 3,466,365 | 0.982077 | 0.017923 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr21 | chr22 | 216 | 3,020,111 | 0.983091 | 0.016909 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr21 | chr13 | 213 | 4,840,543 | 0.981719 | 0.018281 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr15 | chr22 | 167 | 889,312 | 0.954574 | 0.045426 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr15 | chr13 | 168 | 609,567 | 0.956835 | 0.043165 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr14 | chr22 | 234 | 2,745,665 | 0.972180 | 0.027820 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr14 | chr13 | 138 | 1,159,189 | 0.960463 | 0.039537 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr13 | chr22 | 238 | 3,434,388 | 0.985197 | 0.014803 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr22 | chr21 | 177 | 1,493,239 | 0.978433 | 0.021567 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr15 | chr21 | 137 | 668,236 | 0.971510 | 0.028490 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr14 | chr21 | 130 | 1,164,940 | 0.965025 | 0.034975 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr13 | chr21 | 117 | 996,489 | 0.971724 | 0.028276 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr22 | chr15 | 105 | 555,973 | 0.972774 | 0.027226 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr21 | chr15 | 139 | 914,138 | 0.952294 | 0.047706 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr14 | chr15 | 113 | 907,728 | 0.969998 | 0.030002 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr13 | chr15 | 244 | 4,031,996 | 0.986079 | 0.013921 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr22 | chr14 | 96 | 274,107 | 0.981225 | 0.018775 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr21 | chr14 | 90 | 236,429 | 0.969609 | 0.030391 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr15 | chr14 | 95 | 833,970 | 0.973199 | 0.026801 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr13 | chr14 | 83 | 250,544 | 0.983160 | 0.016840 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr22 | chr13 | 136 | 601,830 | 0.982056 | 0.017944 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr21 | chr22 | 188 | 1,460,419 | 0.978911 | 0.021089 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr21 | chr13 | 177 | 1,190,028 | 0.974899 | 0.025101 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr15 | chr22 | 240 | 1,940,725 | 0.982520 | 0.017480 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr15 | chr13 | 240 | 6,192,734 | 0.985055 | 0.014945 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr14 | chr22 | 90 | 386,974 | 0.976345 | 0.023655 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr14 | chr13 | 80 | 273,909 | 0.976863 | 0.023137 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr13 | chr22 | 244 | 2,707,120 | 0.984526 | 0.015474 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr22 | chr21 | 180 | 949,690 | 0.946073 | 0.053927 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr15 | chr21 | 244 | 1,799,590 | 0.981478 | 0.018522 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr14 | chr21 | 165 | 624,284 | 0.959069 | 0.040931 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr13 | chr21 | 222 | 1,619,896 | 0.989192 | 0.010808 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr22 | chr15 | 180 | 646,939 | 0.949513 | 0.050487 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr21 | chr15 | 205 | 7,409,906 | 0.964391 | 0.035609 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr14 | chr15 | 157 | 490,315 | 0.948087 | 0.051913 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr13 | chr15 | 123 | 726,260 | 0.986363 | 0.013637 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr22 | chr14 | 130 | 3,722,783 | 0.949351 | 0.050649 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr21 | chr14 | 207 | 11,213,280 | 0.959970 | 0.040030 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr15 | chr14 | 130 | 467,922 | 0.978981 | 0.021019 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr13 | chr14 | 76 | 361,581 | 0.979346 | 0.020654 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr22 | chr13 | 80 | 339,128 | 0.983113 | 0.016887 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr21 | chr22 | 201 | 957,322 | 0.966931 | 0.033069 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr21 | chr13 | 109 | 570,633 | 0.975158 | 0.024842 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr15 | chr22 | 228 | 1,978,468 | 0.984703 | 0.015297 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr15 | chr13 | 185 | 5,412,584 | 0.988998 | 0.011002 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr14 | chr22 | 117 | 437,815 | 0.977535 | 0.022465 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr14 | chr13 | 116 | 431,131 | 0.975217 | 0.024783 |
| acrocentric_p_cross | sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr13 | chr22 | 240 | 1,867,327 | 0.989584 | 0.010416 |
| acrocentric_p_cross | wfmash_p95_updated_bin | filtered_one_to_one | PAN027mat_vs_PAN010_joint | chr14 | chr22 | 1 | 2,000 | 0.999000 | 0.001000 |
| acrocentric_p_cross | wfmash_p95_updated_bin | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr22 | chr15 | 8 | 16,000 | 0.999375 | 0.000625 |
| acrocentric_p_cross | wfmash_p95_updated_bin | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr21 | chr15 | 8 | 16,000 | 0.987250 | 0.012750 |
| acrocentric_p_cross | wfmash_p95_updated_bin | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr14 | chr15 | 12 | 20,838 | 0.997121 | 0.002879 |
| acrocentric_p_cross | wfmash_p95_updated_bin | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr22 | chr14 | 32 | 60,456 | 0.981449 | 0.018551 |
| acrocentric_p_cross | wfmash_p95_updated_bin | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr21 | chr14 | 61 | 117,403 | 0.979514 | 0.020486 |
| acrocentric_p_cross | wfmash_p95_updated_bin | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr22 | chr13 | 2 | 4,000 | 0.996003 | 0.003997 |
| acrocentric_p_cross | wfmash_p95_updated_bin | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr15 | chr13 | 107 | 200,572 | 0.997907 | 0.002093 |
| acrocentric_p_cross | wfmash_p95_updated_bin | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr15 | chr22 | 9 | 17,993 | 0.996724 | 0.003276 |
| acrocentric_p_cross | wfmash_p95_updated_bin | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr14 | chr22 | 3 | 4,368 | 0.996795 | 0.003205 |
| acrocentric_p_cross | wfmash_p95_updated_bin | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr13 | chr22 | 26 | 50,976 | 0.998444 | 0.001556 |
| acrocentric_p_cross | wfmash_p95_updated_bin | filtered_one_to_one | PAN028mat_vs_PAN027_joint | chr15 | chr21 | 13 | 26,000 | 0.997388 | 0.002612 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN027mat_vs_PAN010_joint | chr22 | chr13 | 3 | 4,986 | 0.997193 | 0.002807 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN027mat_vs_PAN010_joint | chr21 | chr22 | 25 | 49,000 | 0.999694 | 0.000306 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN027mat_vs_PAN010_joint | chr14 | chr22 | 3 | 5,311 | 0.955426 | 0.044574 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN027mat_vs_PAN010_joint | chr13 | chr22 | 25 | 48,000 | 0.999563 | 0.000437 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN027mat_vs_PAN010_joint | chr22 | chr21 | 25 | 48,983 | 0.945088 | 0.054912 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN027mat_vs_PAN010_joint | chr13 | chr21 | 26 | 49,018 | 0.999327 | 0.000673 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr22 | chr15 | 4 | 6,385 | 0.488990 | 0.511010 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr21 | chr15 | 50 | 93,778 | 0.973811 | 0.026189 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr14 | chr15 | 14 | 23,838 | 0.962581 | 0.037419 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr22 | chr14 | 102 | 274,438 | 0.976465 | 0.023535 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr21 | chr14 | 172 | 480,387 | 0.954178 | 0.045822 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr22 | chr13 | 32 | 101,983 | 0.989784 | 0.010216 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr21 | chr22 | 3 | 5,000 | 0.997200 | 0.002800 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr21 | chr13 | 3 | 5,000 | 0.997200 | 0.002800 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr15 | chr13 | 165 | 442,499 | 0.994321 | 0.005679 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr15 | chr22 | 25 | 47,989 | 0.996270 | 0.003730 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr14 | chr13 | 44 | 85,355 | 0.996112 | 0.003888 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr14 | chr22 | 12 | 22,353 | 0.994636 | 0.005364 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr13 | chr22 | 26 | 50,976 | 0.998570 | 0.001430 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr22 | chr21 | 20 | 37,993 | 0.957373 | 0.042627 |
| acrocentric_p_cross | wfmash_p95_updated_bin | raw_many_many | PAN028mat_vs_PAN027_joint | chr15 | chr21 | 60 | 116,913 | 0.996865 | 0.003135 |

## Whole-genome target-support patterns

Top target-support totals by aligned query bp:

| method | evidence | comparison | query arm | target arm | class | bins | aligned bp | weighted identity |
|---|---|---|---|---|---|---:|---:|---:|
| sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr1internal | chr1_internal | same_chromosome | 123,609 | 28,292,200,000 | 0.770571 |
| sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr1internal | chr1_internal | same_chromosome | 123,087 | 14,277,000,000 | 0.788601 |
| sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr1internal | chr1_internal | same_chromosome | 123,883 | 13,974,100,000 | 0.780564 |
| sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr14internal | chr15_internal | interchromosomal | 8,923 | 9,833,460,000 | 0.826383 |
| sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr9internal | chr9_internal | same_chromosome | 67,184 | 8,321,260,000 | 0.792950 |
| sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr16internal | chr16_internal | same_chromosome | 45,788 | 6,093,670,000 | 0.760470 |
| sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr9internal | chr9_internal | same_chromosome | 68,102 | 5,263,970,000 | 0.797013 |
| sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr3internal | chr3_internal | same_chromosome | 100,707 | 5,227,560,000 | 0.841423 |
| sweepga_fastga_f32 | raw_many_many | PAN027pat_vs_PAN011_joint | chr16internal | chr16_internal | same_chromosome | 45,252 | 4,604,050,000 | 0.773273 |
| sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr14internal | chr22_internal | interchromosomal | 10,371 | 4,349,140,000 | 0.825052 |
| sweepga_fastga_f32 | raw_many_many | PAN027mat_vs_PAN010_joint | chr4internal | chr4_internal | same_chromosome | 96,197 | 3,418,650,000 | 0.871360 |
| sweepga_fastga_f32 | raw_many_many | PAN028mat_vs_PAN027_joint | chr5internal | chr5_internal | same_chromosome | 90,520 | 3,416,630,000 | 0.904168 |

## Outputs

- `summaries/bin_target_support_manifest.tsv`: manifest of the 12 full per-bin target-support shard TSVs in `summaries/tmp_worker_bin_support/`.
- `summaries/target_support_totals.tsv`: compact whole-genome support totals.
- `summaries/focal_region_summary.tsv`: Fig5 chr9q->chr3q, PAR, and acrocentric controls.
- `summaries/resource_usage.tsv`: Slurm allocation and helper-thread accounting.
- `scripts/summarize_existing_paf_impg_like_scan.py`: full monolithic combine/best-target implementation; the cluster cancelled the multi-GB materialization step, so the worker shards are retained as the full bin-level output.

## Parquet status

| TSV | Parquet sidecar | status |
|---|---|---|
| target_support_totals.tsv | target_support_totals.parquet | SKIP: pyarrow unavailable in default Slurm Python |
| focal_region_summary.tsv | focal_region_summary.parquet | SKIP: pyarrow unavailable in default Slurm Python |
| bin_target_support_manifest.tsv | bin_target_support_manifest.parquet | SKIP: pyarrow unavailable in default Slurm Python |

## Arm annotation rule

No centromere table is required for this PAF-only scan. Query and target arms are
called from CHM13 chromosome sizes using a 500 kb subtelomeric rule: alignments in
the first 500 kb are p-arm, alignments in the final 500 kb are q-arm, and all other
alignments are marked internal. This is appropriate for the Fig5 subtelomeric
candidate/control readout and keeps the full-genome summary auditable.

## Interpretation

The focal summary separates the PAN027/PAN028 chr9q-to-chr3q signal from PAR and
acrocentric-p positive controls while retaining the genome-wide target landscape.
The raw many:many layer captures broad direct-similarity support; the filtered
one-to-one layer records the corresponding best-chain support after SweepGA filtering.
