# D1 — Descriptive summary of PHR length and community membership by superpopulation

> Provenance for this note: see `response/PROVENANCE.md` (tool versions: `response/notes/tool_versions_and_pins.md`).

Referee #1, minor point 11. Verbatim (`PEER_REVIEW.md`:97):

> 11. The 232 individuals span five HPRC superpopulations, yet no population-stratified
> analysis was performed. Given the "population-scale" framing, and that subtelomeric
> duplicon copy number varies between populations, even a descriptive summary of PHR length
> and community membership by superpopulation would strengthen the claim.

Authority documents read: `PEER_REVIEW.md` (line 97, the verbatim comment above; lines 1–80,
the editor's letter). `paper_prep/synthesis/SIMULATED_PEER_REVIEW_v1.md` was **not** opened.
This note closes the internal tracking row **R1-m11** / item 13 in
`response/POINT_BY_POINT.md` (table at :48), whose state was "L1 (F_ST done)". The `POINT_BY_POINT.md`
file itself is NOT edited here (outside my write scope).

**What this note is.** A write-up of artefacts that already exist. No PHR was re-called for
this note and no statistic in Tables 1–3 was recomputed by hand: every cell is copied from
`scripts/cladistics/phr_by_superpopulation.tsv`, which is the output of
`scripts/cladistics/phr_by_superpopulation.py`. The only arithmetic I did was the
internal-consistency arithmetic in §5.

**This is descriptive, not a test of differentiation.** Nothing here tests whether
populations differ. The differentiation analysis is the Hudson $F_{ST}$ comparison produced by
`scripts/popgen/matched_fst_from_blocks.py` (writing `scripts/popgen/matched_fst_d_m6.tsv`,
reusing the block jackknife and the published continental baselines from
`scripts/popgen/matched_fst_d_m6.py`), and it is reported in the manuscript section
"Communities carry ordinary continental population structure" (`submission/paper.tex`:256–280).

---

## 1. Provenance of the superpopulation mapping (asked for explicitly)

| item | answer |
| --- | --- |
| mapping actually used by the script | **repo file** — `slides/v2-review-zoom/_revision_assets/v9/labels_superpop_stats_polish/nearest_same_superpop_mds_distances.tsv` (4,180,581 B, 15,669 lines). Columns include `Name`, `Sample`, `Haplotype`, `Chromosome`, `Arm`, `Superpopulation`. |
| where that repo file came from | `source_manifest.tsv` in the same directory records `seq_assignment_tsv = /moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv`, sha256 `17bb574dfd4ab661437d7de02c200d2b1abf7e745e9628ea273b77393d8a7f82`; the repo copy is the `Superpopulation` column taken through `make_labels_superpop_stats_polish.R` in the same directory. |
| did I have to read moosefs for the labels? | **No.** The label column is repo-local, so this note is reproducible without the `/moosefs` mount. |
| optional cross-checks on moosefs | `--crosscheck /moosefs/guarracino/HPRCv2/data/hprc-sequence-production.tsv` (15,435 B, header `ChildID N Sex Subpopulation Superpopulation Production Year`) and `--crosscheck-assign /moosefs/.../hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv` (1,410,827 B). Both existed and were readable when I checked them on 2026-10-23 (UTC); the script treats a missing/unreadable cross-check as a reported-and-skipped condition and a label **conflict** as fatal. |
| why the cross-check cannot be the primary source | `hprc-sequence-production.tsv` lacks ~29 of the 232 individuals (see `paper_prep/manuscript_revision/A1_sample_counts.md`), so it confirms a subset only. |

Other inputs of the script, for completeness:

| input | path | status |
| --- | --- | --- |
| PHR call table | `data/all-vs-all.1Mb.p95.id95.len.tsv` | in-repo |
| arm-level Leiden communities (15 communities, 41 arms) | `data/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv` | in-repo |

**CHM13.** CHM13v2.0 is one of the 465 assemblies but is not an HPRC individual; the upstream
heterogeneity tables carry it as `EUR`. The script reports it as its own `reference` section
(38 PHRs on 1 haplotype, 6,840,000 bp) and **excludes it from every superpopulation and cohort
aggregate** below. This is visible as an arithmetic check: the repo label table's own summary
gives EUR as 33 samples / 2,206 MDS points
(`slides/.../nearest_same_superpop_mds_summary.tsv`, quoted in that directory's `README.md`),
while Table 1 gives EUR as 32 individuals / 2,168 PHRs; the difference is exactly CHM13
(1 sample, 38 PHRs).

---

## 2. Table 1 — PHRs and PHR length by superpopulation

Copied from `scripts/cladistics/phr_by_superpopulation.tsv`, `section == superpopulation`
(and the `cohort` row for reference). Length IQR is the `q1_len_bp`–`q3_len_bp` interval;
quantiles are type-7 linear interpolation.

| superpopulation | individuals | haplotypes | arms with ≥1 PHR | PHRs | PHR bp (Mb) | median PHR length (IQR) | total PHR bp | PHR share of cohort | PHR-bp share of cohort | PHR share ÷ haplotype share | PHRs per haplotype, median (IQR) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AFR | 67 | 134 | 41 | 4,549 | 656.47 | 105.0 kb (40.0–205.0 kb) | 656,470,000 | 0.2910 | 0.2923 | 1.0078 | 34 (33–35) |
| AMR | 44 | 88 | 41 | 2,974 | 428.01 | 105.0 kb (40.0–200.0 kb) | 428,005,000 | 0.1903 | 0.1906 | 1.0033 | 34 (33–35) |
| EAS | 52 | 104 | 41 | 3,438 | 480.29 | 105.0 kb (40.0–185.0 kb) | 480,285,000 | 0.2200 | 0.2139 | 0.9814 | 33.5 (31–35) |
| EUR | 32 | 64 | 41 | 2,168 | 319.46 | 105.0 kb (40.0–205.0 kb) | 319,460,000 | 0.1387 | 0.1423 | 1.0056 | 34 (33–35) |
| SAS | 37 | 74 | 41 | 2,501 | 361.45 | 105.0 kb (40.0–200.0 kb) | 361,450,000 | 0.1600 | 0.1610 | 1.0033 | 34 (32.25–35) |
| **ALL_HPRC** | **232** | **464** | **41** | **15,630** | **2,245.67** | **105.0 kb (40.0–200.0 kb)** | **2,245,670,000** | **1** | **1** | **1** | **34 (33–35)** |
| CHM13 (reference, *excluded* above) | 1 | 1 | 38 | 38 | 6.84 | 132.5 kb (45.0–250.0 kb) | 6,840,000 | — | — | — | 38 |

Read plainly: **all five superpopulations carry the same 41 arms**, the same modal PHR length
(105 kb median in each), and a PHR count within 2 % of what their haplotype count predicts
(`PHR share ÷ haplotype share` = 0.98–1.01). The largest deviation is EAS, 1.8 % below
proportional.

## 3. Table 2 — PHR length tails and flank censoring by superpopulation

Same TSV rows; `p10_len_bp` / `p90_len_bp` are deciles, `pct_len_500kb` is the share of that
superpopulation's PHRs whose length equals the whole 500 kb telomere-anchored flank (the
caller's right-censoring point).

| superpopulation | p10 length | median | p90 length | PHRs at full 500 kb flank | share |
| --- | --- | --- | --- | --- | --- |
| AFR | 20.0 kb | 105.0 kb | 320.0 kb | 296 | 6.51 % |
| AMR | 20.0 kb | 105.0 kb | 330.0 kb | 200 | 6.72 % |
| EAS | 20.0 kb | 105.0 kb | 330.0 kb | 242 | 7.04 % |
| EUR | 20.0 kb | 105.0 kb | 330.0 kb | 166 | 7.66 % |
| SAS | 20.0 kb | 105.0 kb | 330.0 kb | 177 | 7.08 % |
| **ALL_HPRC** | **20.0 kb** | **105.0 kb** | **330.0 kb** | **1,081** | **6.92 %** |

The right-censored fraction is the only length statistic that varies appreciably between
superpopulations, and it varies over 6.5–7.7 %, i.e. by ~1.2 percentage points.

## 4. Table 3 — share of PHRs per community, within each superpopulation

Share of *that superpopulation's* PHRs that fall in each arm-level Leiden community
(`section == community`, `frac_phrs`, rounded to 0.1 %; columns sum to 100 % within rounding).
`ALL_HPRC` is the cohort column, for comparison.

| community | ALL_HPRC | AFR | AMR | EAS | EUR | SAS |
|---|---|---|---|---|---|---|
| C1 | 4.6% | 4.3% | 4.7% | 4.7% | 4.2% | 4.8% |
| C2 | 5.7% | 5.7% | 5.7% | 5.7% | 5.6% | 5.7% |
| C3 | 17.0% | 17.3% | 17.2% | 16.3% | 16.8% | 17.1% |
| C4 | 5.7% | 5.8% | 5.8% | 5.7% | 5.7% | 5.6% |
| C5 | 11.0% | 11.1% | 11.0% | 11.2% | 10.9% | 10.9% |
| C6 | 17.2% | 17.4% | 17.1% | 17.3% | 17.3% | 17.2% |
| C7 | 4.9% | 4.7% | 4.4% | 4.9% | 5.8% | 4.8% |
| C8 | 2.9% | 2.9% | 2.8% | 2.9% | 2.9% | 2.9% |
| C9 | 2.8% | 2.8% | 2.9% | 2.8% | 2.8% | 2.8% |
| C10 | 2.9% | 2.8% | 2.9% | 3.0% | 2.8% | 2.9% |
| C11 | 11.4% | 11.3% | 11.4% | 11.6% | 11.5% | 11.2% |
| C12 | 5.8% | 5.8% | 5.7% | 5.9% | 5.8% | 5.7% |
| C13 | 2.9% | 2.9% | 2.9% | 2.8% | 2.8% | 2.8% |
| C14 | 2.8% | 2.8% | 2.8% | 2.7% | 2.6% | 2.8% |
| C15 | 2.7% | 2.7% | 2.8% | 2.6% | 2.5% | 2.7% |

Community membership is essentially flat across superpopulations. The largest departures
from the cohort column are C7 in EUR (5.8 % vs 4.9 %) and C7 in AMR (4.4 % vs 4.9 %); C3
and C6, the two largest communities, track the cohort value within 0.7 percentage points in
every superpopulation. Per superpopulation the per-community PHR-bp shares are in the same
TSV (`frac_bp` column), and there C7 dominates in every superpopulation (15.2 % in AMR to
19.6 % in EUR of that superpopulation's PHR bp) because C7's members are the
nearly-flank-length PHRs (median 500 kb in every superpopulation; 91.6 % in AMR to 95.8 % in
AFR of them hit the 500 kb flank boundary). That bp concentration is a length effect, not a membership
effect, which is why Table 3 is given in PHR counts and §3 reports the censoring share.

## 5. Consistency gates I ran on the existing TSV

All four are `awk` over `scripts/cladistics/phr_by_superpopulation.tsv`; they test the file's
own arithmetic, not the pipeline.

```
$ awk -F'\t' 'NR>1 && $1=="superpopulation"{ni+=$4;nh+=$5;np+=$7;bp+=$8} END{printf "indiv=%d hap=%d phrs=%d bp=%d\n",ni,nh,np,bp}'
indiv=232 hap=464 phrs=15630 bp=2245670000        # == the cohort row exactly
$ awk -F'\t' 'NR>1 && $1=="community"{s[$2]+=$9} END{for(k in s) printf "%s %.6f\n",k,s[k]}'
AFR 1.000000
ALL_HPRC 0.999900
AMR 1.000000
EAS 1.000100
EUR 1.000100
SAS 1.000200
$ awk -F'\t' 'NR>1 && $1=="community"{s[$2]+=$10} END{for(k in s) printf "%s %.6f\n",k,s[k]}'
AFR 1.000100
ALL_HPRC 1.000000
AMR 1.000100
EAS 1.000000
EUR 0.999900
SAS 0.999900
$ awk -F'\t' 'NR>1{c[$1]++} END{for(k in c) print k,c[k]}'
superpopulation 5
reference 16
community 90
cohort 1
```

- The five superpopulation rows add to the cohort row exactly (232 / 464 / 15,630 /
  2,245,670,000), so no PHR is double-counted or dropped between the two sections.
- Per-population community `frac_phrs` and `frac_bp` each sum to 1.0000 ± 0.0002, the
  stated rounding.
- Row counts: 5 superpopulation + 90 community (15 communities × 6 populations:
  ALL_HPRC and the five superpopulations) + 1 cohort + 16 reference (CHM13 + 15
  communities) = 112 data rows, matching `wc -l` = 113 with the header.

## 6. Reproducing the table — and the reproduction I ran

```
python3 scripts/cladistics/phr_by_superpopulation.py            # writes the TSV
python3 scripts/cladistics/phr_by_superpopulation.py --nperm 20000 --seed 20260201
```

I reran the script (default seed and `nperm`, output redirected to a scratch path so the
in-repo TSV was not touched) and compared:

```
$ md5sum scripts/cladistics/phr_by_superpopulation.tsv /tmp/erikg-laneD/phr_by_superpopulation.rerun.tsv
ac0d6116abda54ebcaf159014416ba10  scripts/cladistics/phr_by_superpopulation.tsv
ac0d6116abda54ebcaf159014416ba10  /tmp/erikg-laneD/phr_by_superpopulation.rerun.tsv
$ diff scripts/cladistics/phr_by_superpopulation.tsv /tmp/erikg-laneD/phr_by_superpopulation.rerun.tsv
IDENTICAL
```

byte-identical, so the TSV quoted in Tables 1–3 is what the checked-in script produces from
the checked-in inputs. Wall time: the rerun wrote its output inside one minute of launch
(file mtimes 18:59 and 19:00 UTC, 2026-10-23); not instrumented, and other tenants on the box
were not recorded.

Stdout of that rerun (verbatim, this is the provenance and gate trail for Tables 1–3):

```
NOTE: flank table all-vs-all.1Mb.p95.id95.len.tsv: 18826 flanks, 15668 called PHRs, flank lengths observed = [500000]
NOTE: label cross-check vs hprc-sequence-production.tsv (sample-level metadata): 203 samples labelled by both sources, 0 conflicts, 29 of 232 HPRC individuals absent from that file
NOTE: label cross-check vs hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv (per-PHR assignment): 232 samples labelled by both sources, 0 conflicts, 0 of 232 HPRC individuals absent from that file
NOTE: per-sample PHR counts identical between that file and the call table across all 233 samples
GATE total called PHRs = 15668 (HPRC 15630 + CHM13 38); HPRC individuals with >=1 PHR = 232, haplotypes = 464
GATE haplotypes carrying >=1 called PHR (incl. CHM13) = 465
GATE per-haplotype PHR counts identical between all-vs-all.1Mb.p95.id95.len.tsv and nearest_same_superpop_mds_distances.tsv across all 465 haplotypes
WROTE /tmp/erikg-laneD/phr_by_superpopulation.rerun.tsv (112 rows)

== permutation checks (labels permuted across haplotypes) ==
n_haplotypes_tested=464 nperm=20000 seed=20260201
community x superpopulation (PHR bp, {'EUR': 64, 'EAS': 104, 'AMR': 88, 'AFR': 134, 'SAS': 74}): Cramers V=0.0407 chi2=14886869.9 df=56 p_perm=0.0000
per-haplotype median PHR length: H=16.245 df=4 p_perm=0.0024
Both are DESCRIPTIVE heterogeneity measures on haplotype-level permutations; haplotypes of the same individual are not independent and no population-genetic inference should be drawn from them.
```

So: no label conflict against either cross-check, the per-haplotype PHR counts agree between
the call table and the label table across all 465 haplotypes, and the two label-permutation
heterogeneity measures are Cramér's V = 0.0407 (community × superpopulation PHR-bp,
permutation p < 1/(20000+1) = 5e-5) and Kruskal–Wallis H = 16.245 on 4 df (permutation
p = 0.0024). Both are the script's own descriptive checks on haplotype-level permutations;
they are not a differentiation test and are not used for any claim in this note.

## 7. What this note does NOT establish

1. **It does not test population differentiation.** No test, no interval, no p-value for any
   population comparison appears above. The claim that subtelomeric sharing carries ordinary
   continental structure rests on the matched Hudson $F_{ST}$ comparison
   (`scripts/popgen/matched_fst_from_blocks.py` → `scripts/popgen/matched_fst_d_m6.tsv`;
   `submission/paper.tex`:256–280), which is a different analysis over different units
   (arm × community blocks against published genome-wide baselines).
2. **The two permutation statistics in §6 are descriptive heterogeneity measures, not tests of
   differentiation.** They permute superpopulation labels across haplotypes while preserving
   within-haplotype clustering, treat the 464 haplotypes as exchangeable (they are not: 232
   individuals × 2), and have no bearing on the manuscript's differentiation claim. They are
   not written to the TSV; I captured them by rerunning the script and quoting its stdout in
   §6.
3. **These are PHR-level and haplotype-level descriptives, not individual-level genotypes.**
   The units are telomere-anchored 500 kb flanks with a called PHR, so an arm with no call
   contributes nothing and the per-individual counts are floored at zero, not modelled.
4. **Length is right-censored at the 500 kb flank.** Table 2 reports the censoring share
   because every called region contains the telomeric boundary of its flank by construction
   (the caller walks inward and stops at the first failing window). Median and IQR are
   therefore of *called region length within a 500 kb flank*, not of PHR length in the genome.
5. **Copy number is not measured here.** Referee #1 motivates the request with
   "subtelomeric duplicon copy number varies between populations". What Table 1–3 show is PHR
   *count and length* by superpopulation; copy-number variation is addressed separately by the
   copy-number-aware enrichment analysis (`submission/paper.tex`:670–720), not by this table.
6. **Sample sizes differ between superpopulations** (32 EUR individuals to 67 AFR), and no
   correction or weighting for that is applied — the normalisation offered is the
   `PHR share ÷ haplotype share` ratio in Table 1, which is a proportional-share comparison,
   not a model.
