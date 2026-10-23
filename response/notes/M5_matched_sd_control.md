# M5 — specificity control for the sequence→contact association

Referee #1 major point 5, Referee #3 major point 3. Script: `scripts/hic/matched_sd_control.py`.
Results: `scripts/hic/matched_sd_control.tsv` (199 rows), `scripts/hic/matched_sd_control.json`.

**STATUS: answered as a WINDOW-CLASS control. The referee's literal "matched pericentromeric +
interstitial segmental duplication" census is NOT delivered — it cannot be built from what is on
disk (details in §1). The window-class comparison gives a clear, conservative verdict (§5).**

---

## 1. CAP STATEMENT — the literal request, and why the two named tables cannot carry it

The instruction was to build Axis I ("interstitial, ≥1 Mb from either telomere") from
`data/all-vs-all.1Mb.p95.id95.len.tsv` and `data/all-vs-all.1Mb.p95.id95.len.per_window.tsv`
(the per-window table is not vendored; it was read read-only from
`/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.per_window.tsv`, 521,457 lines,
93 MB). **Both tables are PHR-internal, so they cannot represent interstitial windows.**

Evidence, computed in the same pass (`axis1_containment` in the JSON):

- 521,456 data rows, 18,828 flank sequences, 48 arm labels, 5 kb windows,
  `dist_from_telomere` ∈ [0, 495,000] (the analysed flank is the ~500 kb left after telomere
  trimming, not 1 Mb).
- **450,502 / 508,724 rows (88.56 %) lie inside the sequence's own called PHR interval**
  (checked row by row against `region_start`/`region_end` of the PHR-level table); 12,732 rows
  belong to sequences with no called PHR at all.
- Per-arm window coverage therefore stops where the PHR stops: `chr1_q` has windows at
  0–35 kb only (445 haplotypes × 7 windows + 6 outliers); **14 of 48 arms have no window beyond
  100 kb**, median per-arm maximum distance from the telomere is 257.5 kb.

Consequence: "≥1 Mb from either telomere" cannot be built from either named table. What is
missing, precisely:

1. a **position-resolved, genome-wide, arm-labelled** neighbour table (the ones on disk are either
   PHR-internal, as above, or chromosome-resolved without partner-arm labels);
2. a **segmental-duplication annotation for CHM13** with per-copy identity and copy number — the
   resource the referee names (Jeong et al. / UCSC `hprc2annot` SD track) is **not on disk**:
   `find /moosefs/guarracino/HPRCv2 -iname "*segdup*" -o -iname "*sd.bed*" -o -iname
   "*chm13.dv*"` returns **zero hits**, and `/moosefs/.../PHR_III/annotations/` contains only
   `subtelomeric_annotations.1Mb.rds`;
3. a **per-pair sequence-similarity axis over non-subtelomeric SD copies**. The on-disk
   genome-wide alignments are CHM13-vs-haplotype PAFs whose target contigs are accessions with no
   contig→chromosome map for the 465 haplotypes; the CHM13 self-alignment file
   `aln-ba8599f.../chm13v2.0_maskedY_rCRS.fa.PanSN-vs-chm13v2.0_maskedY_rCRS.fa.PanSN.p95.Pinf.aln.paf.gz`
   is **empty (28 bytes)**. Producing that axis means running a new all-vs-all of extracted SD
   copies (a wfmash/pggb job of my own);
4. a matching table on (identity, length, copy number) between PHR pairs and SD pairs, which needs
   (2) and (3).

So the literal control needs (a) an SD annotation download plus (b) one new alignment run plus
(c) the Hi-C extraction. Everything else (contact machinery, statistics) is in place and gated.

## 2. What was delivered instead — two window classes, same arms, same contacts

Same principle in both substrates: **only the sequence windows change**; the contact matrix and
the pipeline are held fixed. Between-node similarity for a band is

```
share(A,B) = 0.5 * [ #{A-windows whose partner set contains B} / #{A-windows}
                   + #{B-windows whose partner set contains A} / #{B-windows} ]
```

computed by one code path for every band, plus a stricter `share_both_directions = min(f_A, f_B)`.

**AXIS 1 (secondary, descriptive) — arm level, PHR-internal depth bands.**
Substrate: the per-window table (5 kb windows, partner lists carry **arm** labels). Bands are
distance from the telomere *within the called PHR*: 0–100, 100–200, 200–300, 300–400,
400–500 kb. Because the windows are PHR-internal, this is a *depth gradient inside subtelomeric
duplications*, not an interstitial control, and the deep bands are contributed by the few arms
with long PHRs.

**AXIS 2 (primary) — chromosome level, genome-wide, literal ≥1 Mb geometry.**
Substrate: `data/hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz` (31,202 CHM13 100 kb windows
aligned to 25,272 haplotype contigs at ≥98 %; 625 chrY rows dropped because the Hi-C coordinate
space is the masked-Y CHM13 fasta). Bands: telomeric ≤100 kb and ≤500 kb from a telomere;
interstitial ≥1 Mb and ≥2 Mb from **both** telomeres; plus **duplicated-only** variants of the
telomeric-100 kb and interstitial-1 Mb bands (windows that actually carry ≥1 cross-chromosome
duplicate with ≥2 haplotype copies — the closest available analogue of the referee's SD set).
Arm labels are unavailable for partners in this table, so nodes are chromosomes (23), not arms.

**CONTACT (identical for every window class).** Public arm-level mean balanced contact from
`data/human_{CHM13_hic,HG002_hic,HG002_porec,HG002_cifi}_50000bp_seqlevel.tsv`, aggregated to
chromosome pairs for Axis 2 (`contact_to_chrom`, mean over arm pairs). Three mappers:
HiC-Pro+bowtie2 (CHM13, HG002), minimap2+pairtools (HG002 Pore-C, HG002 CiFi).

**Statistics.** Mantel Spearman on node-pair values with **average ranks** (`scipy.rankdata`,
`method="average"`) — explicitly not argsort-of-argsort, which is what the 129/344 zero-tie audit
in `scripts/hic/mantel_mouse_d_m5.py` hinged on; p from 10,000 node-label (row/column)
permutations; 95 % CI from 2,000 node-block bootstrap replicates (resampling nodes, which keeps
the dependence induced by one arm appearing in many pairs); paired bootstrap CI for
rho(T)−rho(I) on identical resampled node sets. Tie counts are written to every result row
(`n_distinct_sim`, `n_tied_sim`, `n_zero_sim`, same for contact).

## 3. Gates (verbatim from the final run)

```
GATE G1: published CHM13 pair-level rho reproduced = 0.71625 (published 0.716) on 652 pairs -> PASS
GATE G: telomeric-band similarity vs published arm matrix: rho=+0.8555 on 820 arm pairs -> PASS
GATE G: axis2 telomeric sharing vs axis1 telomeric sharing aggregated to chromosomes: rho=+0.6902 on 253 chromosome pairs -> PASS
```

G1 pins the contact contract (the published per-PHR-pair Spearman). The second gate shows the
per-window aggregation reproduces the published arm-level Jaccard ranking (1 − dist,
`data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`). The third shows the two substrates measure the
same thing where they overlap, which is what licenses using Axis 2 for the interstitial geometry.

## 4. Results

### 4.1 Axis 1 — PHR-depth gradient (arm level, `share_symmetric`)

| PHR band from telomere | arms | arm pairs | CHM13 Hi-C | HG002 Hi-C | HG002 Pore-C | HG002 CiFi |
|---|---|---|---|---|---|---|
| 0–100 kb | 48 | 630 | **+0.684** (p 1e-4) | **+0.681** (p 1e-4) | **+0.440** (p 1e-4) | +0.306 (p 1e-4) |
| 100–200 kb | 31 | 336 | +0.777 | +0.768 | +0.632 | +0.406 |
| 200–300 kb | 25 | 218 | +0.698 | +0.676 | +0.538 | +0.505 |
| 300–400 kb | 18 | 95 | +0.474 | +0.515 | +0.473 | +0.563 |
| 400–500 kb | 13 | 42 | +0.430 (CI −0.19…+0.88) | +0.470 | +0.397 | +0.676 |

Paired delta, outermost vs innermost band: **CHM13 Hi-C +0.214 [−0.216, +0.791] p=0.40**;
HG002 Hi-C +0.109 p=0.62; Pore-C +0.022 p=0.90; CiFi **−0.178 p=0.16**.

Reading: the association does **not** live only at the outermost window; it is present across the
whole PHR (strongest at 100–200 kb) and the outer-vs-inner contrast is not significant, with the
sign flipping for CiFi. This band is descriptive only — the inner band rests on 11–13 arms, of
which 5 are acrocentric p arms plus chrX/chrY, so it is composition-confounded (see §6).

**Axis-1 coverage table** (windows per arm; the collapse is the PHR extent, §1):

| arm | 0–100 kb | 100–200 | 200–300 | 300–400 | 400–500 |
|---|---|---|---|---|---|
| chr1_p | 8764 | 8758 | 7943 | 1747 | 20 |
| chr1_q | 3301 | 0 | 0 | 0 | 0 |
| chr10_p | 7968 | 12 | 0 | 0 | 0 |
| chr10_q | 7104 | 4020 | 457 | 41 | 2 |
| chr11_p | 8169 | 7481 | 2836 | 107 | 14 |
| chr11_q | 1804 | 0 | 0 | 0 | 0 |
| chr12_p | 8159 | 0 | 0 | 0 | 0 |
| chr12_q | 3973 | 0 | 0 | 0 | 0 |
| chr13_p | 1520 | 1520 | 1520 | 1520 | 1515 |
| chr13_q | 2732 | 10 | 0 | 0 | 0 |
| chr14_p | 4582 | 4560 | 4560 | 4530 | 4448 |
| chr14_q | 1560 | 0 | 0 | 0 | 0 |
| chr15_p | 2380 | 2380 | 2365 | 2360 | 2356 |
| chr15_q | 8979 | 3755 | 4 | 0 | 0 |
| chr16_p | 5369 | 2413 | 966 | 39 | 6 |
| chr16_q | 8360 | 5937 | 745 | 0 | 0 |
| chr17_p | 5676 | 196 | 0 | 0 | 0 |
| chr17_q | 5871 | 0 | 0 | 0 | 0 |
| chr18_p | 8878 | 1994 | 23 | 0 | 0 |
| chr18_q | 1844 | 20 | 20 | 20 | 20 |
| chr19_p | 8662 | 8422 | 1576 | 0 | 0 |
| chr19_q | 3733 | 0 | 0 | 0 | 0 |
| chr20_p | 9010 | 2629 | 15 | 0 | 0 |
| chr20_q | 8298 | 2103 | 1 | 0 | 0 |
| chr21_p | 2318 | 2300 | 2300 | 2263 | 2236 |
| chr21_q | 3598 | 0 | 0 | 0 | 0 |
| chr22_p | 4458 | 4440 | 4440 | 4418 | 4376 |
| chr22_q | 3814 | 0 | 0 | 0 | 0 |
| chr2_p | 1752 | 0 | 0 | 0 | 0 |
| chr2_q | 8987 | 3689 | 204 | 83 | 40 |
| chr3_p | 1812 | 0 | 0 | 0 | 0 |
| chr3_q | 8699 | 8639 | 6132 | 62 | 52 |
| chr4_p | 7236 | 7 | 0 | 0 | 0 |
| chr4_q | 7106 | 5021 | 785 | 108 | 27 |
| chr5_p | 1804 | 0 | 0 | 0 | 0 |
| chr5_q | 8934 | 8812 | 1261 | 0 | 0 |
| chr6_p | 8458 | 3404 | 2 | 0 | 0 |
| chr6_q | 8702 | 5333 | 186 | 0 | 0 |
| chr7_p | 8749 | 8365 | 774 | 4 | 0 |
| chr7_q | 4585 | 0 | 0 | 0 | 0 |
| chr8_p | 7722 | 5755 | 102 | 0 | 0 |
| chr8_q | 1808 | 0 | 0 | 0 | 0 |
| chr9_p | 5526 | 121 | 0 | 0 | 0 |
| chr9_q | 9000 | 8948 | 6298 | 1003 | 0 |
| chrX_p | 6304 | 6080 | 6080 | 6071 | 6060 |
| chrX_q | 6600 | 6600 | 6600 | 3341 | 40 |
| chrY_p | 1806 | 1780 | 1780 | 1758 | 1740 |
| chrY_q | 2020 | 2020 | 2020 | 1016 | 0 |

Arms dropped by the ≥20-window rule: 100–200 kb → chr10_p, chr13_q, chr4_p; 200–300 kb →
chr15_q, chr20_p, chr20_q, chr6_p; 300–400 kb → chr7_p; 400–500 kb → chr10_q, chr11_p,
chr16_p. The 13 arms retained in 400–500 kb are chr13_p, chr14_p, chr15_p, chr21_p, chr22_p,
chrX_p, chrX_q, chrY_p (1,515/4,448/2,356/2,236/4,376/6,060/40/1,740 windows) plus chr18_q,
chr1_p, chr2_q, chr3_q, chr4_q (20/20/40/52/27 windows).

### 4.2 Axis 2 — PRIMARY: telomeric vs genuinely interstitial windows (chromosome level)

| band | windows | nodes / pairs | CHM13 Hi-C | HG002 Hi-C | HG002 Pore-C | HG002 CiFi |
|---|---|---|---|---|---|---|
| telomeric ≤100 kb | 109 | 23 / 242 | **+0.540** (p 1e-4, CI +0.35…+0.66) | **+0.491** | **+0.506** | +0.325 (p 2e-4) |
| telomeric ≤500 kb | 293 | 23 / 242 | +0.548 | +0.536 | +0.512 | +0.349 |
| telomeric ≤100 kb, duplicated only | 87 | 23 / 242 | **+0.581** | **+0.513** | **+0.530** | +0.309 |
| interstitial ≥1 Mb | 30,077 | 23 / 242 | +0.104 (p 0.145, CI −0.17…+0.38) | +0.105 (p 0.13) | +0.154 (p 0.046) | +0.068 (p 0.19) |
| interstitial ≥2 Mb | 29,617 | 23 / 242 | +0.109 (p 0.14) | +0.114 (p 0.11) | +0.165 (p 0.033) | +0.075 (p 0.17) |
| interstitial ≥1 Mb, duplicated only | 1,280 | 22 / 220 | +0.147 (p 0.051, CI −0.13…+0.42) | +0.111 (p 0.098) | +0.156 (p 0.032) | +0.044 (p 0.28) |

Paired deltas rho(T) − rho(I):

| contrast | CHM13 Hi-C | HG002 Hi-C | HG002 Pore-C | HG002 CiFi |
|---|---|---|---|---|
| T(≤100 kb) − I(≥1 Mb) | **+0.436 [+0.14, +0.71] p=0.001** | **+0.386 [+0.10, +0.67] p=0.007** | **+0.352 [+0.09, +0.60] p=0.010** | **+0.257 [+0.008, +0.50] p=0.042** |
| T(≤100 kb, dup) − I(≥1 Mb, dup) | **+0.447 [+0.15, +0.72] p=0.003** | **+0.419 [+0.12, +0.71] p=0.010** | **+0.393 [+0.11, +0.68] p=0.006** | +0.238 [−0.012, +0.50] p=0.062 |

Acrocentric- and sex-excluded sensitivity (chr13, 14, 15, 21, 22, X, Y dropped; 124 pairs):

| band | CHM13 Hi-C | HG002 Hi-C | HG002 Pore-C | HG002 CiFi |
|---|---|---|---|---|
| telomeric ≤500 kb, dup only | +0.689 (p 1e-4) | +0.548 (p 2e-4) | +0.504 (p 1e-4) | +0.374 (p 2e-4) |
| interstitial ≥1 Mb, dup only | +0.115 (p 0.18) | +0.011 (p 0.47) | +0.129 (p 0.16) | −0.013 (p 0.55) |

Non-zero-similarity-only sensitivity (drops tied-at-zero pairs): telomeric +0.54/+0.49/+0.50/+0.32;
interstitial +0.014/+0.094/+0.140/+0.085 with bootstrap CIs spanning zero in all four regimes.

## 5. VERDICT

**Specificity survives, in the window-class sense, but it is a weaker specificity than the
manuscript's wording implies, and it is not the referee's decisive test.**

1. Between-chromosome sequence sharing at the **telomeric** end predicts between-chromosome
   contact: rho = +0.31…+0.58 across four contact datasets from two mapper families, all p ≤ 2e-4,
   bootstrap CIs excluding zero.
2. The **same quantity computed from genuinely interstitial windows (≥1 Mb from both telomeres)
   does not**: rho = +0.04…+0.17, p = 0.03…0.28, bootstrap CIs spanning zero, and it stays at that
   level when the comparison is restricted to interstitial windows that *do* carry a
   cross-chromosome duplicate (1,280 windows).
3. The telomeric-minus-interstitial contrast is **+0.24…+0.45 with 3 of 4 regimes significant at
   p ≤ 0.01** (CiFi marginal, p = 0.062), and it is unchanged, or slightly larger, when both sides
   are restricted to duplicated windows and when acrocentric + sex chromosomes are removed. So the
   result is not a nucleolar-association or arm-type artefact.
4. **But** within the subtelomeric flank the association is *not* confined to the outermost window
   (Axis 1: present at 0–500 kb from the telomere; the outer-vs-inner contrast is not significant
   and reverses sign for CiFi). The manuscript should therefore not claim that the effect is
   specific to the terminal ~100 kb; the defensible claim is "specific to the telomeric compartment
   as a class, relative to interstitial genome".
5. Consequences for the paper's wording: "subtelomere-specific" is supported at the class level
   (telomeric vs interstitial), **not** at the fine-scale level (outermost window), and **not**
   established against the referee's actual comparator (identity/length/copy-matched non-telomeric
   SDs run through the graph-Jaccard pipeline). Because the paper's reported rho = 0.66 is an
   arm-level number on the graph-Jaccard axis and Axis 2 is a chromosome-level number on the
   alignment-free partner-share axis, the absolute values here are **not** comparable to 0.66; only
   the within-substrate T-vs-I contrast is interpretable, and that is what §5.1-5.3 report.

### Scoped statement that can go in the response letter (no new number required)

> "To ask whether the association is specific to the telomeric compartment we replaced the sequence
> axis with one computed from interstitial windows of the same genome, holding the arms/chromosomes
> and the Hi-C contact matrix fixed, and running the identical Mantel machinery with average-rank
> tie handling, 10,000 label permutations and 2,000 node bootstrap replicates. Telomeric sequence
> sharing predicts contact (rho = 0.31–0.58, p ≤ 2 × 10⁻⁴, four contact datasets, two mapper
> families); interstitial sharing does not (rho = 0.04–0.17, p = 0.03–0.28, bootstrap CIs spanning
> zero), and this holds when the interstitial set is restricted to windows that carry
> cross-chromosome duplicates and when acrocentric and sex chromosomes are excluded. The
> telomeric-minus-interstitial contrast is +0.24 to +0.45 (3 of 4 regimes p ≤ 0.01). We therefore
> support a compartment-level specificity claim. We do not claim that the effect is confined to the
> terminal 100 kb: within the subtelomeric flank the association is present from 0 to 500 kb from
> the telomere. A literal identity/length/copy-matched pericentromeric + interstitial segmental
> duplication control through the graph-Jaccard pipeline remains outstanding ([§1] lists what it
> requires); we have not substituted a surrogate for it in the text."

## 6. Residual risks / what could overturn this

1. **Not the referee's test.** Axis 2 compares window *classes* on a partner-share axis; it does not
   compare matched SD *pairs* on the graph-Jaccard axis. A matched-SD version could in principle
   behave differently from a window-class version, in either direction.
2. **Axis 2 is chromosome-level, not arm-level** (partner arms are not resolved in the genome-wide
   table). Arm-level structure is exactly what the paper's communities are about, so Axis 2 is a
   coarser read-out; the third gate (rho = +0.69 against the arm-level matrix aggregated to
   chromosomes) bounds but does not remove this concern.
3. **Sparse x axis interstitially.** Median partners/window for interstitial windows is 0 (mean
   0.15) versus 7 for telomeric windows, i.e. the interstitial class is dominated by unique
   sequence. A null rho there is partly "no variation in x". The duplicated-only band (1,280
   windows) and the non-zero-pairs sensitivity are the safeguards, and both give the same answer,
   but the interstitial duplicated set is still far smaller than the telomeric one, so the CI on
   rho(I) is wide (−0.2…+0.4) and cannot exclude a weak positive association.
4. **Axis-1 composition confound.** The deep PHR bands are carried by 11–13 arms (5 acrocentric
   p arms + chrX/chrY + 5 others); Axis-1 conclusions must not be read as genome-wide.
5. **Axis-2 telomeric bands are coarse**: 109 windows (≤100 kb) and 293 (≤500 kb), i.e. ~4–12 per
   chromosome, so the share values are quantised (median 0.50 at ≤100 kb). Coarse but non-degenerate
   (255/276 pairs non-zero).
6. **Contact-side caveats carry over unchanged** from the published analysis: MAPQ-0 random
   placement of multi-mapping reads (R1-M4/R3-M3), and the arm-level contact matrix used here is the
   published one, whose arm-level Mantel value (0.656, n=38 arms) is *not* reproducible from the
   exported pair table (documented in the previous pass; my aggregation gives 0.725–0.750). The
   within-study contrast uses one contact matrix throughout, so this does not bias T vs I.

## 7. Journal of measurements

- Final run (deterministic, seeded 20260518):
  `python3 scripts/hic/matched_sd_control.py --perms 10000 --boot 2000`
  → **wall 4:59.10, max RSS 100 MB, single process, 199 result rows**
  (`/usr/bin/time -v`). Box: 256 cores, 1 TB RAM; **nothing else of mine was running** during the
  timed run (the two earlier 96/128-thread wfmash jobs had been killed; `ps aux | grep wfmash` = 0).
- Inputs read: the 521,456-data-row per-window table (93 MB, read read-only from moosefs,
  ~20 s), the 31,202-row genome-wide window table, four published pair tables
  (652 / 2,545 / 2,831 / 2,758 lines), and the PHR-level table (18,827 lines) for the containment
  check. No centromere or annotation file is needed for either axis.
- Deliverables: `matched_sd_control.tsv` = 200 lines (header + 199 result rows);
  `matched_sd_control.json` holds coverage, the containment evidence, both sets of similarity
  matrices, the three gates and provenance.
- Reproduce: the script rebuilds everything from the named tables; it writes only
  `scripts/hic/matched_sd_control.{tsv,json}`.

## 8. Aborted earlier approach (kept for the record)

The first attempt at this item derived non-subtelomeric SD windows de novo and planned to run a
`wfmash` all-vs-all over 2,326 extracted CHM13 100 kb windows to obtain a shared-sequence Jaccard
for every window pair. Reconnaissance established, and these findings stand:

- genome-wide CHM13-vs-haplotype PAFs exist (`/moosefs/.../aln-ba8599f.../`, 466 files named
  `*vs-chm13v2.0_maskedY_rCRS*`) but their target contigs carry no contig→chromosome map for the
  haplotypes;
- the CHM13 self-alignment PAF is empty (28 bytes), so pairwise paralog similarity cannot be read
  off the disk;
- no SD annotation anywhere under `/moosefs/guarracino/HPRCv2` (zero hits for `*segdup*`,
  `*sd.bed*`, `*chm13.dv*`);
- `wfmash` silently drops paralogous mappings when window names are PanSN-style (it groups by the
  last `#` and then skips self-group mappings: measured 0 mappings for 5 sequences including exact
  duplicates), keeps only ~1 mapping per segment unless `-f -M` is passed (pilot: 36 vs 351
  records), and writes temp files into the CWD unless `-B` is given;
- the 100 kb window-grid scatter for the SD classes was strongly copy-number imbalanced
  (median distinct partner chromosomes 14 for PHR windows vs 2 pericentromeric / 5 interstitial),
  confirming that the referee's matching requirement is non-trivial rather than cosmetic.

That approach was stopped on instruction as out of scope; the window-class comparison above is
what replaced it. Nothing from it is used in the numbers reported here.
