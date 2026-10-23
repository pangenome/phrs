# IMPG transitivity copy-enumeration of core sequences

Status: COMPLETE (analysis + numbers below are measured; every number is traceable to a
command in §9 or to `scripts/cladistics/impg_copy_enumeration.tsv`).
Author: sxgc.deepseek-flash-writer (child agent).
Authority: `PEER_REVIEW.md` (the real Nature referee report). The simulated review was not read.

## 0. Changelog
- t0 skeleton.
- t1 resources probed. The supplied 37 GB IMPG index is **orphaned from its source
  alignments and cannot answer any CIGAR-dependent query** (failure F1 in §8). Pivoted to a
  self-built, core-scale index; the pivot and its limits are stated in §2 and §8.
- t2 core table verified against `data/chm13.phrs.bed`, the community assignments, the
  Liftoff GFF and the centromere annotation (§3).
- t3 hop-1 enumeration (FastGA `-pafx -f8192 -c85 -s500 -T32`, cores → CHM13) + impg
  index/query; two filter settings (§4, §5).
- t4 hop-2 (transitive closure) measured (§4.4).
- t5 independent cross-check with wfmash for one core (§6).
- t6 TSV written, notes finalised (§7 establishes / does not).
- t7 parent-supplied IMPG-surface facts incorporated and the remaining gap closed by
  testing it: no index-consuming subcommand accepts `--index` alone (§8 F8), the analysis
  is aligner-based enumeration (§8 F9), and the identity-floor sensitivity was measured
  from the on-disk PAF (§4.6).

## 1. Question

The operator's idea: rather than compare *distances between whole duplicated windows*
(an indirect similarity view), take a CORE sequence of length L out of the reference,
follow it through an all-vs-all pangenome alignment representation (transitivity), and
enumerate every other reference position it connects to. That is a direct
copy-number / copy-enumeration view of a subtelomeric duplicon.

This analysis does that for 10 explicitly chosen cores (8 subtelomeric, 2 interstitial
controls), bounded: it never scans the index genome-wide and it reads the reference only
through windowed random access.

## 2. Method actually run (and why it is not the shipped index)

### 2.1 The shipped index is unusable (F1, §8)

`hprc25272.bidirectional.impg` (37,097,351,967 B, mtime 2026-01-19T22:28:06Z) stores, per
PAF record, a BGZF virtual offset into the *original* alignment file. Those files are gone
(`/scratch/hprc25272.aln.paf.gz`, and the whole `/scratch/all-vs-all.*` trees), and the
surviving copies in `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95/` (18,226 PAFs,
2026-01-21) and `.../all-vs-all.1Mb.p95/` (18,827 PAFs, 2026-03-20) are **newer than the
index**, so every stored offset is stale. Any query that needs a CIGAR (bed, bedpe, paf,
and every identity filter) dies; §8 F1 has the verbatim panics. Nothing was estimated in
its place: the alignment representation was rebuilt at core scale.

### 2.2 Reference

`/moosefs/guarracino/HPRCv2/PHR_III/trimmed-assemblies/chm13v2.0_maskedY_rCRS.fa.PanSN.trimmed.fa.gz`
(966,416,611 B; `.fai` + `.gzi` present) — 25 sequences whose names carry their own
chromosome coordinates, e.g. `CHM13#0#chr4:3087-193572740`. Its `.fai` chromosome lengths
are exactly `data/chm13.chrom.sizes` (chr4 = 193,574,945), i.e. the same coordinate system
as `data/chm13.phrs.bed` and the Liftoff GFF.

Staged once to a local uncompressed FASTA (`zcat` 966 MB → 3.10 GB, **24.3 s wall**,
single-threaded, moosefs read, measured 20:23Z) because (a) samtools cannot fetch a name
that itself contains `:`, and (b) FastGA/wfmash must write index sidecars next to the
target. Sequence names were **not** rewritten, so impg's
`--original-sequence-coordinates` converts results back to chromosome coordinates.
After staging, all reference access in `impg_copy_enumeration.py` is windowed
(`FaidxReader`, `.fai` byte arithmetic + `pread` per core/hop-2 span; no whole-genome read,
no text written except the ~150 kb core FASTA and the 5.06 Mb hop-2 query FASTA that the
aligner requires as files).

### 2.3 Alignment layers

1. **hop 1** — FastGA (`/home/erikg/.cargo/bin/FastGA`) `-pafx -f8192 -c85 -s500 -T32`,
   10 cores (140,321 bp total; one is 5 kb) as *query*, the 25 CHM13 sequences (3.10 Gb) as
   *target*. `-pafx` (PAF **with** `cg:Z`) is required — plain `-paf` has no CIGAR and impg
   then cannot compute identity at all. A second run at `-f1000` is reported as a
   seed-frequency sensitivity check (§4.5).
2. **impg layer** — `impg index -f` (forced rebuild) then `impg query -r <core> -d 100
   --min-result-identity 0.95 -l 500 -o bedpe --original-sequence-coordinates --self`.
   impg does the projection, merging, identity computation (gap-compressed `gi` and block
   `bi`) and coordinate remapping; this script only aggregates.
3. **hop 2 (the transitivity step)** — every hop-1 reference placement of a core (≥50 % core
   coverage) is cut out of the reference and aligned back to the reference with the same
   FastGA settings. Positions that are ≥95 % identical to a *copy* but no longer ≥95 %
   identical to the *core* are then reachable. That is the transitivity the operator asked
   for, expressed as an explicit second hop, because the 37 GB representation that would do
   it implicitly is unusable.

### 2.4 Definitions used (all reported, none hidden)

- **placement / `n_positions`**: one distinct reference start coordinate per
  (chromosome, strand) reached by the core, after FastGA fragmentation is merged with
  `--merge-gap 1000`. In a tandem array this counts *unit steps*, which is why it is
  reported separately from:
- **`n_loci`**: merged intervals with a 1 kb gap tolerance (isolated-copy view);
- **`ref_span_bp_all` / `copies_core_units`**: union of reference bp the core connects to,
  and that bp expressed in core lengths (the right copy-number proxy for arrays);
- **`n_arms` / `n_chroms`**: distinct arms (centromere partition:
  `data/chm13-annotations.bed` `Centromere` rows) and distinct chromosomes;
- **`n_on_community_arms` / `n_on_noncommunity_arms`**: placements on arms belonging to the
  core's arm community (`data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`) vs not;
- **`n_interstitial`**: placements that fall in no PHR interval of `data/chm13.phrs.bed`.
- Filters: gap-compressed identity ≥ 0.95, impg `-l 500`, and **core coverage ≥ 0.50**
  (default `--min-core-cov 0.5`); the number of placements thrown away by the coverage
  filter is reported as `n_dropped_short`, because in repeat-dense cores FastGA emits
  thousands of 500–700 bp shared sub-repeat hits (TUBB8B_18p: 4,316 dropped of 4,322).

## 3. The cores (10; 8 subtelomeric + 2 interstitial controls)

Verified programmatically by `cores` (fails loud: exactly one PHR overlap required for
`kind=phr`; zero PHR overlap and ≥5 Mb from both chromosome ends required for `kind=inter`;
arm must exist in the assignments file; the sequence must exist in the reference `.fai`).

| core | locus | len | kind | in PHR | arm | community | genes in core (Liftoff) |
|---|---|---|---|---|---|---|---|
| D4Z4_4q | chr4:193,469,000-193,484,000 | 15,000 | phr | chr4:193,392,741-193,572,740 | 4q | C1 (4q,10q) | 10 (DUX4L6 chain) |
| DUX4L_10q | chr10:134,620,000-134,635,000 | 15,000 | phr | chr10:134,574,995-134,754,994 | 10q | C1 (4q,10q) | 6 (DUX4L10/11/14/22/23) |
| TUBB8B_18p | chr18:195,000-210,000 | 15,000 | phr | chr18:2,017-267,017 | 18p | C2 (10p,18p) | 2 (TUBB8B) |
| TUBB8_10p | chr10:35,000-50,000 | 15,000 | phr | chr10:2,489-72,489 | 10p | C2 (10p,18p) | 2 (TUBB8) |
| OR4F_5q | chr5:181,900,000-181,915,000 | 15,000 | phr | chr5:181,848,925-182,043,924 | 5q | C11 (1p,5q,6q,8p) | 1 (OR4F3) |
| OR4F_6q | chr6:172,005,000-172,020,000 | 15,000 | phr | chr6:171,978,883-172,123,882 | 6q | C11 (1p,5q,6q,8p) | 2 (OR4F7P) |
| SEPTIN14_1p | chr1:142,000-147,000 | 5,000 | phr | chr1:2,706-297,706 | 1p | C11 (1p,5q,6q,8p) | 2 (SEPTIN14P22_5, SEPTIN14P13) |
| rDNA_13p | chr13:5,772,000-5,787,000 | 15,000 | acro | none (see below) | 13p | C7 (13p,14p,15p,21p,22p) | 4 (RNA45S1 unit + 18S/5.8S/28S) |
| INT_REXO1L_chr8 | chr8:86,087,000-86,102,000 | 15,000 | inter | none | 8q, interstitial | — | 2 (REXO1L1P/3P) |
| INT_16p11_2 | chr16:29,500,000-29,515,000 | 15,000 | inter | none | 16p, interstitial | — | 2 |

Notes on the two awkward cores:
- **rDNA_13p**: the 45S array on chr13p (chr13:5,774,303-9,295,790, 76 `RNA45S` genes) is
  **not inside any PHR interval of `data/chm13.phrs.bed`** — that file's chr13 entry
  (chr13:113,553,247-113,563,246) is the q-arm PHR. The core is therefore labelled
  `kind=acro`, `in_phr=False`, arm `13p` (from the centromere partition), and it is the C7
  positive control.
- **SEPTIN14_1p** is 5 kb, not 15 kb. Measured in an earlier run of the identical pipeline
  with a 15 kb core at chr1:138,000-153,000: 18 alignments, 17 positions, 12 arms — but
  **not one placement on 1p itself** (its own 2.6 kb SEPTIN14 paralogs inside the same PHR
  never reach 50 % core coverage). The core is therefore sized to the duplicon unit (5 kb);
  it then enumerates 18 positions of which 1 is on 1p. This is a real property of the
  family, not a threshold artefact, and it is why core length must match duplicon unit
  length or the copy number is silently under-counted at the locus of interest.

## 4. Results

Full table: `scripts/cladistics/impg_copy_enumeration.tsv` (37 columns, one row per core).
Headline per core (`-f8192` primary run):

| core | positions | loci | chroms | arms | ident. median (min–max) | core cov. median | ref span | copies (core units) | on community arms / not | hop-2 new positions / arms |
|---|---|---|---|---|---|---|---|---|---|---|
| D4Z4_4q | 61 | 1 | 2 | 2 | 0.9937 (0.9900–0.9985) | 1.000 | 220,536 | 14.70 | **61 / 0** | 0 / 0 |
| DUX4L_10q | 61 | 1 | 2 | 2 | 0.9940 (0.9926–0.9997) | 1.000 | 220,530 | 14.70 | **61 / 0** | 0 / 0 |
| TUBB8B_18p | 5 | 5 | 3 | 3 | 0.9990 (0.9674–0.9991) | 1.000 | 83,221 | 5.55 | 4 / 1 | 14 / 5 |
| TUBB8_10p | 6 | 6 | 3 | 3 | 0.9705 (0.9554–0.9709) | 1.000 | 91,658 | 6.11 | 4 / 2 | 18 / 6 |
| OR4F_5q | 16 | 16 | 12 | 14 | 0.9949 (0.9722–0.9971) | 0.872 | 213,504 | 14.23 | 2 / 14 | 86 / 14 |
| OR4F_6q | 6 | 6 | 6 | 6 | 0.9961 (0.9883–0.9965) | 1.000 | 100,714 | 6.71 | 2 / 4 | 9 / 5 |
| SEPTIN14_1p | 18 | 18 | 13 | 13 | 0.9826 (0.9706–0.9976) | 1.000 | 90,890 | 18.18 | 1 / 17 | 21 / 8 |
| rDNA_13p | 222 | 222 | 5 | 5 | 0.9973 (0.9873–0.9998) | 1.000 | 3,317,612 | 221.17 | **222 / 0** | 7 / 2 |
| INT_REXO1L_chr8 | 68 | 11 | 1 | 1 | 0.9940 (0.9926–0.9969) | 0.799 | 789,925 | 52.66 | — | 235 / **1** |
| INT_16p11_2 | 1 | 1 | 1 | 1 | 0.9713 | 0.580 | 23,722 | 1.58 | — | 0 / 0 |

### 4.1 The method recovers the known biology without being told it

- **D4Z4_4q and DUX4L_10q are mirror images**: 61 positions each, 30 on 4q + 31 on 10q
  (D4Z4) and 31 on 4q + 30 on 10q (DUX4L), identity 0.990–0.999, every placement covering
  the full core. These are the same 4q↔10q duplicon pair seen from both ends, i.e. the C1
  community recovered from sequence alone. The placements are unit-offset steps of the
  tandem array (one 3.3 kb D4Z4 unit apart, see §8 F4), which is why `n_loci`=1 while
  `n_positions`=61 and `ref_span_bp` ≈ 220 kb.
- **rDNA_13p** hits exactly the five acrocentric p arms — 13p:76, 15p:51, 21p:56, 22p:22,
  14p:17 — i.e. 222/222 placements on C7 arms, no other chromosome. This is the positive
  control and it comes out clean.
- **TUBB8/TUBB8B** cross-hit: the 10p core finds 18p:4, 16q:1, 9q:1; the 18p core finds
  10p:1, 16q:1, 18p:3. Copy number 6 (10p core) and 5 (18p core) — the same family seen
  from either end, plus 16q and 9q partners.

### 4.2 Per-arm copy-number distribution

`per_arm` in the TSV (format `arm:count;...`). Examples:
- D4Z4_4q `10q:31;4q:30`; DUX4L_10q `10q:30;4q:31`.
- rDNA_13p `13p:76;14p:17;15p:51;21p:56;22p:22`.
- OR4F_5q `10p:1;11p:1;16q:1;19p:1;1p:1;1q:2;20p:1;20q:1;2q:1;3q:1;4q:2;6q:1;7p:1;9q:1`.
- SEPTIN14_1p `10p:1;11p:1;16q:1;19p:1;1q:3;20q:1;2q:1;3q:1;4q:1;5q:1;7p:3;9q:1;Y:2`
  (Y has no centromere row in `data/chm13-annotations.bed`, so it is labelled
  `Y_unknown`; reported, not silently mapped).
- INT_REXO1L_chr8 `8q:68`; INT_16p11_2 `16q:1`.

### 4.3 Community members vs non-members

Per core, placements on arms of the core's own Leiden community vs all other arms
(`n_on_community_arms` / `n_on_noncommunity_arms`):

| core | community | on community arms | on other arms |
|---|---|---|---|
| D4Z4_4q | C1 | 61 | 0 |
| DUX4L_10q | C1 | 61 | 0 |
| TUBB8B_18p | C2 | 4 | 1 |
| TUBB8_10p | C2 | 4 | 2 |
| OR4F_5q | C11 | **2** | **14** |
| OR4F_6q | C11 | 2 | 4 |
| SEPTIN14_1p | C11 | **1** | **17** |
| rDNA_13p | C7 | 222 | 0 |
| interstitial controls | — (no arm community) | — | 69/69 on interstitial, non-PHR positions |

So the arm-community partition is strongly reproduced by the two "self-defining" families
(D4Z4/DUX4L, rDNA) and by TUBB8/TUBB8B, but **individual OR4F and SEPTIN14 cores reach
mostly arms outside their core's community**. These are not contradictory results: the
communities summarise whole-PHR similarity, and single families cross boundaries. It is a
caution against reading a community as a family grouping.

### 4.4 hop-2 (transitive closure) adds real positions

290 hop-1 copies were cut out and re-aligned (5,056,198 bp of query, 606,735 PAF records,
**226.1 s wall**). New positions reached only through a copy (≥95 % to a hop-1 copy,
≥1 kb block, ≥1 kb away from any hop-1 start):

| core | new positions | new arms | where |
|---|---|---|---|
| D4Z4_4q | 0 | 0 | — (array already saturated at hop 1) |
| DUX4L_10q | 0 | 0 | — |
| TUBB8B_18p | 14 | 5 | 18p:4, 9q:4, 3q:4, 16q:1, 4q:1 |
| TUBB8_10p | 18 | 6 | 18p:8, 9q:3, 3q:3, 4q:2, 10p:1, 16q:1 |
| OR4F_5q | 86 | 14 | expands the whole OR4F spread; adds 9p (hop-1 had none) and takes 5q from 1 to 7 |
| OR4F_6q | 9 | 5 | 15q:2, 3q:2, 5q:2, 6q:2, 19p:1 |
| SEPTIN14_1p | 21 | 8 | 1q:4, Y:4, 1p:3, 11p:2, 19p:2, 20p:2, 4q:2, 7p:2 |
| rDNA_13p | 7 | 2 | 13p:3, 21p:4 (array-internal extension) |
| INT_REXO1L_chr8 | 235 | **1** | 8q only |
| INT_16p11_2 | 0 | 0 | — |

Hop 2 adds 390 positions in total (464 hop-1 → 854), and it adds **arms** for 6 of the 8
subtelomeric cores. For the chromosome-restricted interstitial array (REXO1L) hop 2 adds
235 more positions but still **zero new arms** — the transitivity step amplifies copy
number without changing the distribution, which is exactly the distinction this analysis
needs.

### 4.5 Seed-frequency sensitivity (`-f1000` vs `-f8192`)

`alt_*` columns in the TSV. `n_positions` and `n_arms` are **identical for all 10 cores**
between `-f8192` and `-f1000`. The only difference is in the raw row count of
TUBB8B_18p (4,322 rows at `-f8192` vs 71 at `-f1000`) and hence in what the ≥50 %-coverage
filter discards (`n_dropped_short` 4,316 vs 65): the extra rows at the permissive setting
are short shared sub-repeat hits that the coverage filter removes, and both settings then
report the same 5 positions. Conclusion: at this threshold the copy-number answer is not
driven by the FastGA seed-frequency setting.

### 4.6 What the 95 % identity floor hides (measured, not estimated)

The FastGA PAF was produced with the aligner default `-i0.7`, so alignments down to 70 %
identity are present on disk and can be re-filtered. Same block (≥1 kb) and core-coverage
(≥50 %) filters, only the gap-compressed identity floor changed (`n_positions` / `n_arms`):

| core | ≥0.95 | ≥0.90 | ≥0.85 |
|---|---|---|---|
| D4Z4_4q | 61 / 2 | 61 / 2 | 61 / 2 |
| DUX4L_10q | 61 / 2 | 61 / 2 | 61 / 2 |
| TUBB8B_18p | 5 / 3 | 6 / 4 | 8 / 5 |
| TUBB8_10p | 6 / 3 | 7 / 4 | 9 / 5 |
| OR4F_5q | 16 / 14 | 16 / 14 | 16 / 14 |
| OR4F_6q | 6 / 6 | 6 / 6 | 6 / 6 |
| SEPTIN14_1p | 18 / 13 | 21 / 14 | 21 / 14 |
| rDNA_13p | 222 / 5 | 222 / 5 | 222 / 5 |
| INT_REXO1L_chr8 | 68 / 1 | 68 / 1 | 68 / 1 |
| INT_16p11_2 | 1 / 1 | 1 / 1 | 1 / 1 |

Six of ten cores are unchanged down to 85 %; the three that move (TUBB8B, TUBB8, SEPTIN14)
gain 1–3 positions and at most 2 arms. **So the 95 % floor is not the binding constraint in
this analysis** — the binding constraints are FastGA's seeding (`-f8192 -s500`) and the
≥50 % core-coverage filter (§3, SEPTIN14 at 15 kb; F7 for the TAF11L experiment). This also
means the comparison against the index's documented 95 % threshold is like-for-like for
these cores.

## 5. Timing and cost (honest)

Host: 256 cores, 1 TB RAM, shared — other jobs were running throughout, so every wall
below is an upper bound for this configuration. FastGA was given `-T32` (its own hard cap
is 32; `-T64` is rejected). `%CPU` reported by `/usr/bin/time -v` was 519–576 %, i.e. the
job was itself memory/IO-bound, not CPU-saturated.

| step | command (abridged) | wall | notes |
|---|---|---|---|
| reference staging | `zcat …PanSN.trimmed.fa.gz > chm13.ref.fa; samtools faidx` | 24.3 s | 966 MB → 3.10 GB, single thread, moosefs |
| core extraction + verification | `impg_copy_enumeration.py cores` | 56.8 s and 65.2 s (two runs) | dominated by streaming the 56 MB Liftoff GFF once per core |
| FastGA hop 1 (`-f8192`) | `-pafx -f8192 -c85 -s500 -T32` | 115.1, 115.7, 120.0, 120.8 s (four runs) | 5,008 PAF records for the final 10-core set (32,787–32,920 when the earlier TAF11L core was in it); maxRSS 2.37 GB; 24 MB PAF |
| `impg index` | `-f -t 32` | 0.0–0.1 s | 2.3 MB index (f8192), 360 KB (f1000) |
| `impg query` ×10 cores | as §2.3 | **0.01–0.31 s each** | 4,916 rows total on the primary run (4,322 of them, 88 %, are TUBB8B_18p's short sub-repeat hits, which the coverage filter then drops); max 0.31 s at 156 rows |
| FastGA hop 2 | same settings, 290 spans / 5.06 Mb query | 226.1 s | 606,735 PAF records, 695 MB PAF |
| wfmash cross-check | `-p 95 -w 10000 -s 5000 -c 500 -t 32` | 46.2 s | maxRSS 6.18 GB |
| aggregate | pure Python parse + join (no impg call) | 70.9 s, maxRSS 1.25 GB | dominated by parsing the 695 MB hop-2 PAF (606,735 rows) |

Whole pipeline for the 10 cores, warm reference staging: **≈6.5 min wall on 32 threads**.
The impg layer itself is essentially free (< 1 s for all 10 cores).

**Genome-wide feasibility (extrapolation, labelled as such).** I have two measured points on
the same target with the same settings: 0.15 Mb of query → 120 s, and 5.06 Mb → 226 s. The
marginal cost is ≈21.6 s per Mb of query at 32 threads, plus a ~100 s fixed cost (target
index build). A genome-wide set of, say, 3,000 cores × 15 kb = 45 Mb of query would
therefore be ≈16 min of alignment and produce PAF volume at the measured rate of ≈140 MB
per Mb of query, i.e. ≈6 GB of PAF and a proportional impg index — all of which is
tractable. **What is not tractable without new data is the full all-vs-all representation**:
rebuilding the 37 GB index requires the original 18k–37k window PAFs or a fresh all-vs-all
of 25,272 haplotypes, which was not attempted (and was explicitly out of scope: no
genome-wide scans). So a genome-wide *core* version is feasible; a genome-wide
*all-vs-all transitivity* version needs the alignment representation to be rebuilt first.

## 6. Cross-check of one core against independent methods

`D4Z4_4q` (chr4:193,469,000-193,484,000), same filter
(gap-compressed identity ≥ 0.95, block ≥ 1 kb, core coverage ≥ 0.50):

| source | alignments | positions | arms | per arm |
|---|---|---|---|---|
| impg bedpe on FastGA PAF (primary) | 62 | 61 | 2 | 10q:31; 4q:30 |
| raw FastGA PAF, re-filtered in this script (bypasses impg) | 61 | 61 | 2 | 10q:31; 4q:30 |
| **wfmash** `-p 95 -w 10000 -s 5000 -c 500` (independent aligner) | 22 | 22 | 2 | 10q:12; 4q:10 |

- impg vs raw PAF: identical positions and identical per-arm counts (1 alignment differs
  because impg's `-l 500` merge split one record). This validates the impg query layer.
- wfmash agrees on the **structure** (two arms, same two arrays, same identity regime,
  gap-compressed identity 0.990–1.000 after CIGAR re-scoring) but reports **22 instead of 61**
  positions. Reading the wfmash PAF shows why: wfmash chains across adjacent D4Z4 units
  into single blocks with large internal deletions (e.g. `matches=14,983`,
  `block=18,301`), so it enumerates roughly one placement per two–three 3.3 kb units where
  FastGA emits one per unit. The copy-number *ratio* between the two tools is therefore an
  algorithmic/definitional property of tandem arrays (§7, §8 F4), not a disagreement about
  where the copies are. **The brief's suggested cross-check tool was FastGA; FastGA is
  already the primary engine here, so an independent aligner (wfmash) was used instead** —
  this is a deliberate deviation, documented rather than hidden.

## 7. What this establishes, and what it does NOT

**Establishes (for these 10 cores, at ≥95 % gap-compressed identity, in CHM13):**
1. A direct copy-enumeration per core: number of distinct reference positions, arms,
   chromosomes, identity distribution (median and range) and core coverage — TSV columns
   `n_positions`, `n_arms`, `n_chroms`, `identity_median/min/max`, `core_cov_median`,
   plus the per-arm distribution.
2. That the two "self-defining" duplicon systems come out of sequence alone: D4Z4/DUX4L as
   a 4q↔10q pair (61 positions each, mirror-symmetric per-arm counts) and rDNA as exactly
   the five acrocentric p arms (222 positions, 222/222 on C7 arms).
3. That a transitive second hop is informative at this scale: +390 positions overall and
   +arms for 6 of 8 subtelomeric cores, while adding zero new arms for the
   chromosome-restricted interstitial array.
4. A bounded answer to the referee's "is this just a generic property of high-identity
   segmental duplications?" for these 10 cores: **copy number alone does not separate the
   classes** (the interstitial REXO1L array has 68 hop-1 copies and 52.7 core-units of
   reference span, more than six of the eight subtelomeric cores), but **chromosomal
   distribution does**: all 8 subtelomeric cores connect to ≥2 chromosomes (median 4,
   max 13) and both interstitial controls connect to exactly **1** chromosome and 1 arm
   (medians: 1 arm / 1 chromosome, versus 4 arms / 4 chromosomes for subtelomeric cores).
   The interstitial placements are also shorter on the reference side (median locus length
   44.8 kb for REXO1L, 8.7 kb for 16p11.2) and cover the core less completely (core-coverage
   median 0.690 across the two controls versus 1.000 for the eight subtelomeric cores).
   With n=2 interstitial controls and no identity/length/copy-number matching, this is a
   suggestive contrast, not a matched-control result.

**Does NOT establish:**
1. **Not an alignment-free truth.** Every number is conditional on: FastGA finding the
   alignment, the ≥95 % gap-compressed identity threshold, the ≥50 % core-coverage filter,
   the 1 kb merge gap, and CHM13 being the only reference. Lower-identity copies
   (90–95 %) are invisible here.
2. **Index/representation completeness bounds it, and this is aligner-based, not an index
   query** (F8, F9). The shipped all-vs-all index is unusable, so the representation was
   rebuilt from 10 cores only. This is a *query-anchored* enumeration, not an all-vs-all
   graph, and cross-arm relations that would need a haplotype bridge (a copy that is
   <95 % to the core but ≥95 % to an intermediate haplotype) are only partly recovered by
   hop 2, which uses the reference as its own bridge. The identity floor itself is not the
   binding constraint (§4.6: unchanged down to 85 % for 6/10 cores); the binding
   constraints are FastGA seeding and core-coverage.
3. **Nothing about which copies are older, and nothing about concerted evolution.** High
   identity between copies is equally explained by recent duplication and by
   homogenisation; this analysis has no outgroup, no haplotype comparison, no polarity. It
   must not be used as evidence for concerted evolution (agrees with referee #1 point 1).
4. **Nothing about copy-number variation between individuals.** Only CHM13 is enumerated.
   The referee's question "how variable are D4Z4/DUX4 arrays across the HPRC" (referee #2,
   L152-154) is not answered here; the hop-2 machinery would be the natural way to answer
   it, but the per-haplotype sequences would have to be added as a second target.
5. **Not repeat-masked.** RepeatMasker BEDs could not be fetched (F5), so repeat content
   per core is not reported; PHR/community survival under repeat masking (referee #1
   point 8) is untouched.

## 8. Failures and anomalies (running log)

### F1 — `hprc25272.bidirectional.impg` cannot resolve CIGAR data (orphaned index)
Evidence:
- `impg query -i …hprc25272.bidirectional.impg -a <placeholder> -r CHM13#0#chr10:38000-53000
  -d 1k -o bed` → `thread 'main' panicked at src/impg.rs:1217: index out of bounds:
  the len is 1 but the index is 3048` (the index references ≥3,049 alignment files; `-i`
  still requires an `--alignment-*` argument and `alignment_file_index` indexes that list).
- With the shipped list `all-vs-all.1Mb.p95.paf.list` (18,827 entries, all under
  `/scratch/...`, deleted): `Failed to open compressed file
  '/scratch/all-vs-all.1Mb.p95/hprcv2.1Mb…paf.gz'` (impg.rs:623).
- With a rewritten list pointing at the surviving copies in
  `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95/` (18,827 `*.paf.gz` present) and
  at `all-vs-all.p95/` (18,226 `*.paf.gz`), the query dies in BGZF seek:
  `slice index starts at 22159 but ends at 0` (noodles-bgzf `io/block/data.rs:45`).
- Timestamps: index 2026-01-19T22:28:06Z; `all-vs-all.p95/*.paf.gz` 2026-01-21;
  `all-vs-all.1Mb.p95/*.paf.gz` 2026-03-20 — both surviving PAF sets are **newer** than the
  index, so the stored BGZF virtual offsets cannot address them.
Conclusion: unusable for anything needing a CIGAR (bed/bedpe/paf, identity filtering).
Reported, not worked around silently.

### F2 — `impg index` silently reuses a stale index
After the core table changed, `impg index -i <existing>` reloaded the old index (only a
warning: "Alignment file … has been modified since impg index creation") and every query
then panicked with `Invalid CIGAR operation: `. Fixed by always passing `-f` (force
rebuild). Recorded because it produces a *plausible-looking* wrong answer path.

### F3 — impg's default filters are wrong for copy enumeration
`impg query` applies a plane sweep / scaffolding stack (`--overlap 0.95`,
`--scaffold-jump 50k`, `--scaffold-mass 10k`). Measured: it returned 3 loci for D4Z4_4q
where 62 alignments pass the same identity threshold. Raising `-l` to ≥1 kb and keeping
the default filter recovered all 62 records for this core, so no `--no-filter` was needed
here — but any copy-enumeration use of impg must check `-l`/merge behaviour explicitly.

### F4 — "copy number" for tandem arrays is definitional
In a tandem array every copy maps to a start offset one unit away from its neighbour, so
`n_positions` (FastGA, one per unit) and wfmash's chained blocks (one per 2–3 units) differ
by a factor ~2.8 for D4Z4. Both are reported (`n_positions`, `n_loci`, `ref_span_bp_all`,
`copies_core_units`); the arm/chromosome distribution is invariant to the choice.

### F5 — RepeatMasker index not fetchable
`https://human-pangenomics.s3.amazonaws.com/HPRCv2/PHR_III/repeat_masker_bed_hprc_r2_v1.0.index.csv`
returns `NoSuchKey`. Repeat content per core is therefore not reported (nothing estimated).

### F6 — FastGA argument/environment quirks (cost time, recorded for reproducibility)
- `-T64` rejected: "GIXmake: # of threads can be at most 32".
- `-P <dir>` with a space is rejected; must be `-P<dir>`.
- FastGA writes `.gdb/.1ano/.post` stubs beside the *query*, so the query must live in a
  writable private directory (an earlier run died with
  `FastGA: Cannot open post stub file …/cores.post`).
- `wfmash -s` must be ≤ `-w`; `-s 5000` with the default window (1000) is an immediate
  error, and wfmash needs a `.fai` for the query.

### F7 — a 600 bp-unit gene duplicon is invisible to a 15 kb core
Before INT_16p11_2 was chosen, `INT_TAF11L_chr5` (chr5:17,430,000-17,445,000, the TAF11L
array, 600 bp genes at ~3.4 kb spacing) was run through the identical pipeline: 268 PAF
records, but the largest alignment anywhere was 8,705 bp at 93.5 % identity to
chr5:17,607,620-17,616,272 (i.e. a non-array paralogue), and **nothing at the array locus
itself** passed the filters. Recorded because it is a property of the method, not a bug:
when the duplicon unit is much shorter than the core and the spacers are not conserved, a
long core cannot enumerate copies. Core length must be matched to unit length (compare the
SEPTIN14 15 kb vs 5 kb case in §3).

### F8 — IMPG surface probing: no subcommand accepts `--index` alone
Independently probed (parent-supplied facts plus my own tests), in this order:
- installed binary `impg 0.5.0`, source tree at `v0.5.0-552-gbe61a4e` (a dev build, not the
  tag); index magic bytes `IMPGIDX2` (format v2), confirmed with `xxd -l 16`.
- a pinned `v0.4.1` build is not available: `/home/erikg/impg-v0.4.1` fails at its vendored
  wfmash-rs `build.rs` cmake step (`gmake: *** [Makefile:136: all] Error 2`), so the older
  reader cannot be tried.
- `impg query` requires, in order: an alignment source (`Either --alignment-files or
  --alignment-list must be provided`), then `-d/--merge-distance is required`, then a target
  range/bed.
- `impg stats  -i <index>` → `Either --alignment-files or --alignment-list must be provided`.
- `impg similarity -i <index> -r …` → `-d/--merge-distance is required`; with `-d 0` →
  `The 'seqwish' and 'pggb' engines are not yet implemented for the similarity command`;
  with `--fast` → `Either --alignment-files or --alignment-list must be provided`.
- `impg partition -i <index> --window-size …` takes neither `-r` nor `-b`
  (`unexpected argument`), and its only alignment input is `-a`/`--alignment-list`.

Conclusion: **every index-consuming subcommand hard-requires the source alignment files.**
That is path (b) in the parent's guidance: confirmed, and it is exactly why the shipped
index cannot be used. The reason is not only "the PAFs are gone" but that the surviving PAF
sets are *newer* than the index, so even supplying them cannot work (F1).

### F9 — this analysis is aligner-based enumeration, not an index query
Stated plainly to avoid any misreading of §4: the numbers come from
`FastGA → (my own small index) → impg query`, i.e. **placements are found by FastGA**, and
are bounded by FastGA's settings (`-f8192 -c85 -s500 -T32`, plus a ≥95 % gap-compressed
identity post-filter and a ≥50 % core-coverage filter), not by the shipped index's stated
95 % threshold and not by its (unusable) all-vs-all graph. §4.6 measures how much the
identity floor matters: little. Flags differ slightly from the parent's suggestion
(`-T4 -f8192 -s500 -c70 -i0.6`): `-T32` is FastGA's own hard cap here and was used to keep
the run short; `-c85` was used with a 95 % post-filter instead of `-c70 -i0.6`, and since
FastGA's `-i` defaults to 0.7 the more permissive alignments are present on disk anyway,
which is what makes §4.6 possible. A cross-check against the shipped index was impossible
for the reasons above; the substitutes are the raw-PAF recount and the independent wfmash
run in §6.

## 8b. Verification gates (verbatim)

Cheap gates, in the order they were used:

1. **Windowed-read geometry** (`selftest`) — the `.fai` byte arithmetic must land exactly on
   a known gene start:
   ```
   seq CHM13#0#chr10:2489-134754994 chrom chr10 beg 2489
   chr10:39416+60 = cccaggccATCGTGATGATGGAGTGCCTGGGGTCCCACTTGCCCACCCTTCATTGCTGGT
   ```
   (upper-case, i.e. non-soft-masked, and an ORF-start-containing stretch: TUBB8 is at
   chr10:39,416-69,557 on the minus strand per the Liftoff GFF.)
2. **Self-placement invariant** — for every core exactly one placement must coincide with
   the core's own interval (±1 bp); `core_stats` aborts otherwise. Verbatim example row:
   ```
   CHM13#0#chr10  35001  50001  TUBB8_10p  0  15000  TUBB8_10p:0-15000  0  +  +  gi:f:1  bi:f:1
   ```
   (also fixes the coordinate convention: impg `--original-sequence-coordinates` returns
   1-based inclusive, i.e. the core written as 0-based half-open 35000-50000 comes back as
   35001-50001.)
3. **Known-truth gates** — rDNA_13p lands on 13p/14p/15p/21p/22p and nothing else
   (222/222 placements on C7 arms); D4Z4_4q and DUX4L_10q are per-arm mirror images
   (`10q:31;4q:30` vs `10q:30;4q:31`). Both were recovered without being told the answer.
4. **impg layer vs raw alignments** (cross-check, §6): `fastga-raw-paf` 61 positions
   `10q:31;4q:30` vs `impg-bedpe` 61 positions `10q:31;4q:30`.
5. **Independent aligner** (cross-check, §6): wfmash 22 positions, 2 arms.
6. **Whole-pipeline reproducibility** — a full re-run into a *fresh* work directory
   (`cores → fastga -f8192 → index → query → hop2 → aggregate`, 442 s wall) produced
   **byte-identical** hop-1 and hop-2 columns to the reported TSV:
   ```
   AGREE: two independent FastGA/impg runs produce identical hop1+hop2 columns
   ```

## 9. Reproduction

```bash
cd /home/erikg/phrs
S=scripts/cladistics/impg_copy_enumeration.py
python3 $S selftest                       # faidx window geometry sanity check
python3 $S cores                          # verify + write cores.tsv / cores.fa
python3 $S fastga  --run f8192 --freq 8192    # hop 1  (~2 min, 32 threads)
python3 $S index   --run f8192
python3 $S query   --run f8192                # impg bedpe per core (<1 s total)
python3 $S fastga  --run f1000 --freq 1000    # sensitivity run
python3 $S index   --run f1000
python3 $S query   --run f1000
python3 $S hop2    --run f8192                # transitive hop 2 (~4 min)
python3 $S crosscheck --core D4Z4_4q --wfmash <wfmash.paf>
python3 $S aggregate                      # -> scripts/cladistics/impg_copy_enumeration.tsv
```
Scratch goes to `$TMPDIR/impg_copy_enum/work` (~750 MB, dominated by the 695 MB hop-2 PAF);
the only repo artifacts are the script, the TSV and this note. Reference staging
(`chm13.ref.fa`, 3.10 GB) is external to the script and documented in §2.2.

## 10. Relation to the referee report (`PEER_REVIEW.md`)

- Referee #1 point 5 asks for non-subtelomeric duplications, matched on identity, length and
  copy number, run through the identical pipeline. §4.4/§7 give the *unmatched* version of
  exactly that comparison for two interstitial controls, with the honest finding that
  copy number does not separate the classes while chromosome distribution does. It is the
  right axis to match on next; it is not yet a matched control.
- Referee #1 point 1 (concerted evolution needs more than identity) and point 8 (repeat
  masking) are the two places where this analysis must not be over-read; §7 says so.
- Referee #2 (L152-154, copy-number variability of D4Z4/DUX4 across HPRC) is the natural
  next step and is explicitly out of scope here.
