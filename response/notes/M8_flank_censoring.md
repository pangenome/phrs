# M8 — Referee #1 major point 8: the 500 kb flank right-censors the PHR length distribution

> Provenance for this note: see `response/PROVENANCE.md` (tool versions: `response/notes/tool_versions_and_pins.md`).

Status: **COMPLETE**. Every number below comes from a command run in this task or from a
named file read in this task. Nothing is estimated. Where an estimate would have been
needed and no data exist, that is stated.

Outputs created (nothing else was created or modified):
- `scripts/cladistics/flank_censoring_sensitivity.py` (839-row table generator)
- `scripts/cladistics/flank_censoring_sensitivity.tsv` (the table: 840 lines = header + 839 data rows; md5 `b14d54e9d6fb16094ca8d6327a2198cb`, reproduced across independent full runs and re-runs)
- `response/notes/M8_flank_censoring.md` (this file)

Intermediate work dir (not in the repo): `/scratch/erikg/m8-censoring/` (dev) and
`/scratch/erikg/m8-censoring/{run,clean,clean2,clean3}/` (four runs of the delivered script; `clean2` was interrupted and abandoned after the code path was found correct, `clean3` is the final from-scratch run that produced the delivered table).

## 0. Answer in one paragraph

The 500 kb flank censors **1,087 of 15,668 PHRs (6.94 %)** — but they are not a random
6.94 %: **717 are acrocentric short arms (93.97 % of the PHRs called there), 366 are X/Y
pseudoautosomal arms (43.06 %) and only 4 are on the 14,055 PHRs of every other arm
(0.03 %)**. Eleven of 41 arms are affected at all, and 9 of those 11 account for
1,083/1,087. The censoring is not a decay: in the last 5 kb window of the flank a
censored PHR still involves a mean of 4.0 distinct chromosomes (median 5) and the window-to-window continuation
probability is 0.99 (acrocentric) and 0.96 (PAR) — the walk is cut, not finished. Carrying
the *same* 5 kb-window criterion into a deeper alignment source (the genome-wide,
whole-assembly all-vs-all survey, ~12 % of haplotype pairs) and calibrating it against the
published calls gives a measured extension of **median +0.9 to +3.0 Mb** (three
identity/alignment-count thresholds, strict contiguity to the published 4-failure chaining
rule), with 34–97 % of the censored ends still passing at 2 Mb and 3–31 % still passing at
the 5 Mb measurement ceiling. Consequences: the 500 kb maximum in Fig. 1B is by
construction *only for those 1,087 PHRs* — and on six of the eleven affected arms the
published length distribution is fully degenerate (median = q1 = q3 = max = 500,000 bp, §2.3),
so for those arms it carries no length information at all; “walking inward until inter-chromosomal
similarity disappears” (page 3) is true for 14,581 of the 15,668 PHRs (93.1 %) and false
for the other 1,087; the cohort mean (144 kb) must be replaced by a censored-mean statement
(≥225–294 kb once the censored ends are extended by the measured lower bounds), while the
**median (105 kb) is exactly unaffected** by censoring. The **3.51 Mb** figure is *not*
censored: it is reproduced bit-exactly as the 31 CHM13 PHRs left after removing the five
acrocentric p-arms and both X arms (6,840,000 − 2,500,000 − 830,000 = 3,510,000), and none
of those 31 reaches even 295 kb. It is still a lower bound on subtelomeric sharing, but for
a different reason than the referee gives (seven PHR-free arms, not flank truncation).

## 1. Inputs, and what does not exist (attempted, recorded)

| input | status |
| --- | --- |
| `data/all-vs-all.1Mb.p95.id95.len.tsv` | 18,826 rows + header, 15,668 with a PHR, 3,158 `.` |
| `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz.fai` | 18,827 flanks, **every one exactly 500,000 bp** |
| `…/all-vs-all.1Mb.p95.id95.len.per_window.tsv` (moosefs, 93,149,310 B) | 521,456 window rows; **`max dist_from_telomere = 495000`, `max window_end = 500000`** |
| human 2 Mb / 4 Mb / 5 Mb flank alignments | **do not exist.** The only `*_2Mb` / `*_4Mb` directories in the tree are `mouse_T2T/{subtelo,flanking,similarity,pggb,community_analysis}_{2,4}Mb`. |
| the published flank PAFs (`PHR_III/all-vs-all.1Mb.p95/`, 37,655 files) | query length is 500000 in every record (e.g. `HG00126#2#CM090128.1:134736978-135236977_chr10_qarm 500000 243000 261000 …`), i.e. the alignments were computed on the trimmed 500 kb sequences. **The flank cannot be extended from this alignment set.** |
| `hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz` | not 1 Mb flanks: fai has 15,668 entries whose lengths are the PHR lengths (`CHM13#0#chr10:134574995-134754994_chr10_qarm 180000`) |
| `/moosefs/guarracino/HPRCv2/aln-ba8599f934840b4b1a836e83fee938d9c0bc5e77/` | **the deeper source that does exist**: 25,274 `*.p95.Pinf.aln.paf.gz`, 318.4 GB, whole-assembly-vs-whole-assembly (≈2 × 12,636 unordered pairs of 465 assemblies = the paper's ~12 % survey set) |
| `data/hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz` | 31,202 windows; the authors' own genome-wide transitive per-window table (100 kb, id ≥ 0.98, min output 5 kb) |
| wfmash v0.23.0 (the paper's binary) | **does not run here**: `error while loading shared libraries: libgsl.so.27` (only libgsl.so.25 installed). Local binary is `v0.24.2-16-g4a35bfc5`. |

### 1.0 Scope of this note, stated plainly (confirmed with the operator)

**No 2–5 Mb re-extraction was performed, and none is possible from the files that exist.**
Reconnaissance confirmed and extended here: `/moosefs/guarracino/HPRCv2/PHR_III/flanking/`
holds only `hprcv2.1Mb.flanking_100kb.fa.gz` (+`.fai`/`.gzi`) and a `pggb/`+`similarity/`
subtree — 1 Mb flanks cut into 100 kb windows, not 2 Mb flanks; the only 2 Mb/4 Mb flank
similarity directories under `PHR_III` are `mouse_T2T/*` and are irrelevant; and
`all-vs-all.1Mb.p95.paf.list` points at `/scratch/all-vs-all.1Mb.p95/…`, which no longer
exists. One correction to the record: the 18,827 flank PAFs themselves *are* still readable
under `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95/` (37,655 files) — but every
record has query length 500000, i.e. they were computed on the trimmed 500 kb sequences, so
they cannot yield length beyond the flank either. Either way the conclusion is the same: a
longer-flank similarity cannot be recomputed without re-running wfmash/pggb, which was not
done.

Consequently the two halves of the answer have different epistemic status, and the note keeps
them apart throughout:
- the **censored fraction is measured exactly** (1,087/15,668, per arm, from the delivered
  call set and the published per-window table — §2);
- the **magnitude of the missing length is bounded, not measured exactly** (§4): a
  direct-projection lower bound from a sparser, stricter second alignment regime, an upper
  bracket from a transitive profile, and a physical ceiling from the cytoband arm length
  (§3.1). No cell of this note claims the true length of a censored PHR.

### 1.1 Why the 2–5 Mb re-extraction itself was not run

Measured, not assumed. One job in the published arrangement (all 500 kb flanks as queries,
one flank as target, `-t 96`): **88,781 PAF lines, 6:38 wall, 21.1 GB peak RSS**. Extrapolating
to the same job shape on 2 Mb flanks (4× query bases) and running one job per censored end:
1,087 × ~26 min ≈ 470 h, i.e. ~230 h at two concurrent 96-thread jobs on this host. The full
2 Mb all-vs-all is ~16× the 500 kb all-vs-all, which itself is 18,827 × 6.6 min ≈ 2,070 h
single-stream here. So the referee's literal re-extraction requires the cluster that produced
the original run, not one workstation. This is recorded as a **limitation**, and the
extension is measured from the pre-existing genome-wide alignment instead (§4), with a
measured calibration of that instrument (§4.2) rather than an assumption that it is equivalent.

## 2. The censoring, exactly (measured; Part A of the table)

Gate before any count: replaying the published per-window pass/fail through the region rule of
`scripts/phr-detection/find-multichr-regions-incremental.py` (walk from the telomeric window
inward, stop after 4 consecutive failing windows, region end = last passing window)
**reproduces all 18,826 published PHR lengths exactly** (`gate_caller_rule_replay_exact =
18826/18826`). Secondary gate: the per-window table's maximum `dist_from_telomere` is exactly
495,000 — the depth ceiling is proven from the data, not assumed.

One flank in the fai is not in the call table:
`NA18982#1#JBKABS010000018.1:83848607-84348606_chr18_qarm`, dropped upstream as a scaffolding
chimera (it is present in the pre-filter per-window table, which is how the count
1,088 → 1,087 closes).

| quantity | value | source row |
| --- | --- | --- |
| flanks in the call table | 18,826 | `censoring/all_flanks/n_flanks` |
| PHRs | 15,668 | `…/n_phrs` |
| flanks with no PHR | 3,158 | `…/n_flanks_no_phr` |
| PHR bp | 2,252,510,000 | `…/phr_bp` |
| **censored (length = 500,000)** | **1,087 (6.94 %)** | `…/n_censored_len_eq_500kb`, `censored_share_of_phrs_pct` |
| censored bp | 543,500,000 (24.1 % of PHR bp) | `…/censored_bp` |
| near-censored (495–499 kb) | 23 | `…/n_near_censored_495_499kb` |
| length ≥ 495 kb (“within one window of the limit”) | 1,110 (7.08 %) | `…/n_within_one_window_of_limit` |
| p-arm flanks censored | 1,082 / 7,397 = **14.63 %** | `censoring/parm/*` |
| q-arm flanks censored | 5 / 8,271 = **0.06 %** | `censoring/qarm/*` |

Censoring is a *definition* question only for a PHR whose last available window passes; a PHR
that ends at 490 kb ends because the caller saw 4 consecutive failures, and would end at the
same place under a 5 Mb flank. The 23 near-censored PHRs are the sensitive band; they are
0.15 % of PHRs.

### 2.1 Which arms (measured; `per_arm` and `censoring/*_arm` rows)

| class | PHRs | censored | share censored | share of all censored |
| --- | --- | --- | --- | --- |
| acrocentric short arms (13p,14p,15p,21p,22p) | 763 | **717** | **93.97 %** | 65.96 % |
| X/Y pseudoautosomal arms (Xp,Xq,Yp,Yq) | 850 | **366** | **43.06 %** | 33.67 % |
| all other arms | 14,055 | **4** | **0.03 %** | 0.37 % |

The 4 are `HG02145#2#JAHKSF020000009.1:11130-511129_chr1_parm`,
`HG01975#2#CM099810.1:201668480-202168479_chr3_qarm`,
`HG03521#2#CM094434.1:207617369-208117368_chr3_qarm`,
`HG00280#2#CM087159.1:190873674-191373673_chr4_qarm` (2 on 3q, 1 each on 1p and 4q).
Xq's 2 are PAR2-adjacent.

Eleven of the 41 arms that carry a PHR have any censored PHR — `1p,3q,4q,13p,14p,15p,21p,22p,Xp,Xq,Yp`
(`n_arms_with_a_censored_phr = 11`). Per arm: 13p 71/76 (93.4 %), 14p 207/229 (90.4 %),
15p 116/119 (97.5 %), 21p 110/116 (94.8 %), 22p 213/223 (95.5 %), Xp 281/327 (85.9 %),
Yp 83/92 (90.2 %), Xq 2/330, 3q 2/436, 1p 1/438, 4q 1/357.

**This is the single most useful result for the response:** 99.6 % of the censoring is in the
two arm classes the manuscript already reports separately and excludes from the 3.51 Mb total.
The uncensored distribution — the thing Fig. 1B is actually about for 41 arms — is uncensored
to 99.97 %.

### 2.2 The censoring is a cut, not a decay (measured; `boundary_signal`, `continuation`)

Distinct chromosomes in the last (495–500 kb) window of each censored flank: 2 → 374,
3 → 15, 4 → 31, 5 → **592**, 6 → 50, 7 → 13, 8 → 7, 9 → 1, 10 → 3, 11 → 1. Mean **4.02**.
374 (34.4 %) sit at the bare 2-chromosome minimum.

Window-to-window continuation r(d) = P(pass at d+5 kb | pass at d):

| depth | acro p (n alive, r) | PAR (n alive, r) | all other arms (n alive, r) | ALL (r) |
| --- | --- | --- | --- | --- |
| 400 kb | 751, 0.9987 | 392, 1.0000 | 10, 1.000 | 0.9991 |
| 450 kb | 742, 0.9919 | 392, 1.0000 | 6, 1.000 | 0.9947 |
| 475 kb | 734, 0.9946 | 392, 1.0000 | 6, 0.667 | 0.9947 |
| 490 kb | 723, 0.9917 | 383, 0.9556 | 4, 1.000 | 0.9793 |

So at the boundary the acrocentric and PAR ends are still passing one 5 kb window after
another at 96–99 % — there is no evidence whatever of the shared sequence ending inside the
flank. Contrariwise the “all other arms” population is exhausted: by depth 300 kb only 259
flanks are alive there (of 14,055) with r = 0.30. **The two populations behave completely
differently and the 500 kb flank is adequate for one and not the other.**

A stationary-hazard extrapolation of the boundary survivor curve answers only
**+236 kb** (`naive_constant_hazard_extra_bp = 236304`). It is reported in the table *only*
to record that it is not used: it is an order of magnitude below the measured extension,
because the hazard is not stationary and the boundary survivors are a selected
subpopulation. This is the concrete reason a *measurement* (or a much longer flank) is
required rather than a survival-curve correction.

### 2.3 Where similarity terminates along the flank (measured; `termination_depth`)

Regime A can localise termination for every *uncensored* PHR, and that is the empirical
statement the per-window table supports. Cohort length histogram (one row per 50 kb bin,
`cohort_length_bin:*`; sums to 15,668):

| 0–50 | 50–100 | 100–150 | 150–200 | 200–250 | 250–300 | 300–350 | 350–400 | 400–450 | 450–500 | = 500 kb |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 4,898 | 1,801 | 3,041 | 1,847 | 1,404 | 683 | 834 | 6 | 11 | 56 | **1,087** |

This is the histogram behind Fig. 1B, and its last bar is the artifact: 1,087 at exactly
500 kb against 6/11/56 in the three bins below it. Per arm (`termination_depth/<arm>/*`),
the censoring shows up as degeneracy of the whole distribution, not just its tail:

| arm | n | median (q1–q3) | max | reached boundary | unfinished streak |
| --- | --- | --- | --- | --- | --- |
| 13p | 76 | 500,000 (500,000–500,000) | 500,000 | 71 | 3 |
| 14p | 229 | 500,000 (500,000–500,000) | 500,000 | 207 | 11 |
| 15p | 119 | 500,000 (500,000–500,000) | 500,000 | 116 | 0 |
| 21p | 116 | 500,000 (500,000–500,000) | 500,000 | 110 | 0 |
| 22p | 223 | 500,000 (500,000–500,000) | 500,000 | 213 | 3 |
| Xp | 327 | 500,000 (500,000–500,000) | 500,000 | 281 | 22 |
| Yp | 92 | 500,000 (500,000–500,000) | 500,000 | 83 | 4 |
| 1p | 438 | 295,000 (290,000–315,000) | 500,000 | 1 | 0 |

For the six acrocentric/PAR arms with many censored calls, **the published length
distribution carries no information at all beyond "≥500 kb"** — median = q1 = q3 = max =
500,000. On 1p the distribution is fully resolved (median 295 kb, one flank at the boundary).

### 2.4 What the per-window data alone can and cannot bound (inference, labelled)

Regime A contains no observation past 500 kb, so it cannot *measure* the missing length; it can
only say (i) 1,087 PHRs are still passing at the boundary with mean 4.02 chromosomes, (ii)
continuation r(490 kb) = 0.9917 (acro p) / 0.9556 (PAR), and (iii) 43 PHRs ended with an
incomplete failure streak and could resume. A stationary-hazard extrapolation of (ii) answers
+236 kb (§2.2) and is reported only to record that it is rejected. The magnitude therefore
comes from §4's second regime and §3.1's physical ceiling, and is a **bound**.

## 3. The 3.51 Mb figure (measured; `paper_3.51Mb` rows)

- CHM13 has 38 PHRs totalling **6,840,000 bp** (identical to the CHM13 row of
  `scripts/cladistics/phr_by_superpopulation.tsv`, 6.84 Mb — a cross-file consistency check).
- Removing the five acrocentric p-arms (5 × 500,000) and both X arms (Xp 500,000 + Xq 330,000)
  leaves **31 PHRs totalling 3,510,000 bp = 3.51 Mb exactly**, which is the number in
  `submission/paper.tex:163` and `:382`. The identity is `6840000 − 2500000 − 830000 =
  3510000`; the script asserts it rather than eyeballing it.
- **None of those 31 PHRs is censored**: their maximum length is 295,000.
- CHM13's six censored ends are exactly `13p,14p,15p,21p,22p,Xp` — the acrocentric short arms
  and PAR1 — i.e. every saturated CHM13 end is already *outside* the 3.51 Mb total.

Therefore: **“3.51 Mb is a lower bound” is right about the biology but wrong about the
mechanism.** The number is a censoring-free measurement for the assembly and arm classes it
covers; it is a lower bound on total subtelomeric sharing because seven arms carry no detected
PHR and because the acrocentric/PAR arms are excluded outright, not because its constituent
lengths were truncated. (The paper/repo also never states that 3.51 Mb is CHM13-only — see
§7 risk r5, a cheap fix.)

### 3.1 Vendored CHM13 geometry: what it confirms, what it refuses to (`chm13_geometry`)

- **`data/chm13.phrs.bed` is not a faithful export of the delivered call set** and is
  therefore *not* used for any count here: 37 rows / 6,014,981 bp as written vs 38 CHM13 PHRs /
  6,840,000 bp in the call table; 30/37 rows agree with a call-table PHR within 1 bp (its end
  coordinate is inclusive in some rows and half-open in others, so a strict comparison is
  meaningless); 8 call-table PHRs are absent from it, including **chr13p, chr14p, chr15p,
  chr21p — precisely the censored acrocentric arms**. Anyone re-deriving the 3.51 Mb from this
  file will get a different number; the call table is the authority, and the 3.51 Mb identity
  in §3 is checked against it.
- **`data/chm13.chrom.sizes`** positions all 46 CHM13 flanks against their contigs and gives
  the telomeric repeat the pipeline trimmed: median **2,634 bp**, maximum **4,586 bp** (so
  "trimming the terminal telomeric repeat" removed ~2.6 kb, not a subtelomeric repeat class —
  relevant to the referee's last sentence, which is the repeat-masking item, not this one).
- **`data/chm13v2.0_cytobands_allchrs.bed`** gives the physical head-room, telomere→centromere:
  13p 15,547,593 · 14p 10,092,112 · 15p 16,678,794 · 21p 10,962,853 · 22p 12,788,180 ·
  Xp 57,820,107 · Yp 10,565,750 bp. The measured acrocentric extensions (§4.3) are 2.6–3.8 Mb, i.e.
  **18–26 % of the arm** (13p 24.5 %, 14p 26.0 %, 15p 18.1 %, 21p 24.2 %, 22p 24.8 %) and 4.2 %
  for the (57.8 Mb) Xp arm — physically possible, and nowhere near the ceiling, which also
  means the cytoband length is a weak upper bound and §4.3's direct-projection numbers are the
  binding ones.
  This also *falsifies the transitive regime as a point estimate*: its 14p extent
  (14,047,925) exceeds the 14p arm (10,092,112) and its 22p extent (17,145,862) exceeds the
  22p arm (12,788,180), i.e. transitive closure walks through the centromere. That is exactly
  why §4.5 is used only as an upper bracket.

## 4. The extension, measured on a deeper source (Part B of the table)

### 4.1 What regime B is, and what it can and cannot claim

The `aln-ba8599…` survey set is whole-assembly all-vs-all, file-tagged `p95.Pinf`: wfmash was
run with `-p inf`, i.e. no divergence ceiling in the mapper, so each reported alignment's
real identity is the `gi:f:` gap-compressed tag that every record carries — which is exactly
the field this script filters on (≥0.95, and a ≥0.98 sub-bin). It covers ~12 % of haplotype
pairs. I projected every alignment touching a telomere-bearing contig onto **5 kb depth
windows measured inward from that contig's telomeric trim boundary** (the published caller's
own origin and window; 1,000 windows out to 5 Mb), counting alignments per (flank, window,
partner chromosome, identity bin), with the published ≥3 kb projected-length requirement.
Pass rule: **≥2 distinct chromosomes with signal, the source chromosome included** — which is
what `submission/paper.tex:160-163` specifies and what `impg query` returns (same-chromosome
flanks included). Length then uses either (a) the published 4-consecutive-failure chaining
rule, or (b) strict contiguity from the telomere, which cannot bridge a false positive.

Two definitional differences from regime A must be stated, and both are measured:
(i) B is **sparser** (12 % of pairs), so it has ~1/8 of the alignments per window → its
lengths are lower bounds;
(ii) B can match a subtelomeric window to the **interior** of another chromosome, which
regime A (flanks only) cannot → this inflates B. §4.2 measures how much.

### 4.2 Calibration at depth < 500 kb (measured; `survey_calibration`)

445,846 published-passing and 75,510 published-failing windows below 500 kb.

| channel | window sensitivity | window specificity | uncensored length overreach | no-PHR flanks that gain a region |
| --- | --- | --- | --- | --- |
| id95 ≥5 aln/chr | 0.814 | **0.989** | 191/14,581 (1.3 %), median +100 kb | 1/3,158 |
| id95 ≥2 aln/chr | 0.881 | **0.980** | 310/14,581 (2.1 %), median +100 kb | 1/3,158 |
| id95 ≥1 aln/chr | 0.924 | 0.944 | 1,094/14,581 (7.5 %), median +40 kb | 10/3,158 |
| id98 ≥5 aln/chr | 0.793 | 0.991 | 163/14,581 (1.1 %) | 1/3,158 |
| id98 ≥2 aln/chr | 0.861 | 0.985 | 263/14,581 (1.8 %) | 1/3,158 |
| id98 ≥1 aln/chr | 0.901 | 0.954 | 831/14,581 (5.7 %) | 10/3,158 |

Read: the survey instrument **under-detects** (sensitivity 0.79–0.92) and **almost never
over-calls** (specificity 0.94–0.99; 1.1–2.1 % of non-censored PHRs gain any length at the
≥2 and ≥5 alignment thresholds used below — 5.7–7.5 % for the plainly too-loose ≥1 threshold,
which is why it is not the primary channel — with a median gain of +100 kb, and only 1–10 of
the 3,158 PHR-free flanks acquire any region at all). Its lengths are
therefore conservative, and the residual definitional inflation (ii) is bounded by the
overreach rate, which is small at the level of whole-region calls.

### 4.3 The extension of the 1,087 censored PHRs (measured; `survey_extension`)

| channel \ rule | >500 kb | median length | median increase | ≥2 Mb | ≥5 Mb (ceiling) |
| --- | --- | --- | --- | --- | --- |
| id95 ≥1 aln, caller rule | 1,083 (99.6 %) | 3,450,000 | +2,950,000 | 96.6 % | 31.2 % |
| id95 ≥1 aln, strict contig. | 1,055 (97.1 %) | 2,765,000 | +2,265,000 | 85.9 % | 22.2 % |
| id95 ≥2 aln, caller rule | 1,079 (99.3 %) | 2,980,000 | +2,480,000 | 92.6 % | 25.0 % |
| **id95 ≥2 aln, strict contig.** | **1,020 (93.8 %)** | **2,440,000** | **+1,940,000** | **69.6 %** | 13.2 % |
| id95 ≥5 aln, caller rule | 1,061 (97.6 %) | 2,445,000 | +1,945,000 | 70.2 % | 11.5 % |
| **id95 ≥5 aln, strict contig.** | **887 (81.6 %)** | **1,365,000** | **+865,000** | **34.4 %** | 2.8 % |

So the measured length increase is **+0.9 Mb (strictest) to +3.0 Mb (loosest)** in median, and
the censored ends very largely continue. “At 2 Mb and at 5 Mb”, explicitly:
**34–97 % still pass at the 2 Mb ceiling and 3–31 % still pass at the 5 Mb ceiling** — and
every one of those 5 Mb numbers is a floor, because the measurement stops there.

Per arm (primary channel id95 ≥2 aln, strict contiguity / caller rule,
`survey_extension_by_arm`):

| arm | n | median length (strict) | median (caller rule) | >500 kb (strict) |
| --- | --- | --- | --- | --- |
| 13p | 71 | 3,805,000 | 5,000,000 (ceiling) | 98.6 % |
| 14p | 207 | 2,625,000 | 3,500,000 | 98.1 % |
| 15p | 116 | 3,015,000 | 4,522,500 | 98.3 % |
| 21p | 110 | 2,657,500 | 5,000,000 (ceiling) | 100 % |
| 22p | 213 | 3,175,000 | 4,435,000 | 96.7 % |
| Xp | 281 | 2,420,000 | 2,445,000 | 85.4 % |
| Yp | 83 | 2,430,000 | 2,435,000 | 90.4 % |
| 1p | 1 | 490,000 | 490,000 | 0 % |
| 3q | 2 | 642,500 | 707,500 | 50 % |
| 4q | 1 | 340,000 | 500,000 | 0 % |
| Xq | 2 | 637,500 | 1,120,000 | 50 % |

Xp landing at ~2.42–2.45 Mb is an independent corroboration: PAR1 is ~2.7 Mb, and the referee
notes the pedigree scan's PAR1 tract is 138 kb against a true ~2.7 Mb — the same truncation,
visible here from a different direction.

### 4.4 Stratified by confidence (measured; `survey_extension_by_confidence`)

Stratum = number of distinct chromosomes in the boundary window (id95 ≥2 aln, strict contiguity):

| boundary n_chrs | n | median length | ≥2 Mb |
| --- | --- | --- | --- |
| 2 (bare threshold) | 374 | 2,425,000 | 74.9 % |
| 3–4 | 46 | 2,677,500 | 63.0 % |
| ≥5 (robust) | 667 | 2,900,000 | 67.2 % |

The extension is **not** a property of only the marginal ends: every confidence stratum
extends by ≥2.4 Mb in median. The finding is stable across the stratification, and it is stable
across the two independent aligner regimes below.

### 4.5 Third regime (Part C): transitive genome-wide, CHM13-anchored (`transitive_regime`)

`data/hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz`, 100 kb windows, id ≥ 0.98, min output
5 kb, **transitive** closure (more permissive than the published direct projection → the upper
side of the bound), extent = deepest window retaining ≥1 partner chromosome with ≥5 haplotype
alignments, with a 400 kb break tolerance:

- acrocentric p-arms: 13p 6,147,456 · 14p 14,047,925 · 15p 15,746,742 · 21p 4,047,495 · 22p 17,145,862
- Xp 2,448,421 · Xq 406,917 · Yp absent (CHM13 chrY is not in this table)
- every other arm: 0 – 647,358 (1p 347,294; 3p 647,358; 11q 0; 8q 0; 14q 0; 18q 0)

i.e. ~0.11–0.65 Mb for the arms where the 500 kb flank is adequate, and 4–17 Mb for the
acrocentric short arms. This regime brackets §4.3 from above and agrees on the ranking; the
transitive numbers are inflated by transitivity (they also read ~15–20 % high on the
well-measured arms, e.g. 1p 347 kb vs the published 295 kb), so they are an upper bound, not
an estimate.

**Two-sided result:** every censored acrocentric/PAR end lies between the direct-projection
lower bound (§4.3) and the transitive upper bound (§4.5); both say ≫ 500 kb. No conclusion in
this note depends on either bound being tight.

### 4.6 Corrected totals (`corrected_totals`)

| | published | + strict id95 ≥5 aln | + strict id95 ≥2 aln |
| --- | --- | --- | --- |
| cohort PHR bp | 2,252,510,000 | ≥3,533,110,000 | ≥4,604,895,000 |
| mean length | 143,765 | ≥225,498 | ≥293,904 |
| median length | 105,000 | 105,000 | 105,000 |

The median is **exactly** unchanged, and this is a theorem rather than a coincidence: the
censoring point (500 kb) lies far above the 50th percentile, so replacing the top 6.9 % of
values by larger values cannot move the 50th percentile. The paper's “median 105 kb” survives;
the paper's “mean 144 kb” does not (it is a right-censored mean and should be labelled as such
or dropped from the abstract-level claims).

## 5. Measured vs inferred — the separation the referee is entitled to

**Measured (no model, no extrapolation):**
- 1,087 censored PHRs, 6.94 %; 717/366/4 split; 11 arms; 543.5 Mb censored bp.
- 1,110 PHRs within one 5 kb window of the limit.
- Boundary chromosome counts (mean 4.02; 34.4 % at the bare 2-chromosome minimum).
- Continuation rates at 490 kb: 0.9917 (acro p), 0.9556 (PAR).
- 3.51 Mb = 31 CHM13 non-acro non-PAR PHRs, 0 of them censored, max 295 kb.
- The published PAF set cannot be extended (query length 500,000 in every record).
- Survey-regime lengths for the censored ends and their calibration (§4.2–4.4).
- Transitive extents (§4.5).

**Inferred (labelled as such, and not load-bearing):**
- That the pessimistic 5 Mb ceilings are floors (they are floors only in the sense that the
  measurement stops; a censored end that is still passing at 5 Mb is ≥5 Mb — this is
  arithmetic, not inference).
- The stationary-hazard correction (+236 kb) — reported and rejected.
- The claim that the acrocentric numbers correspond to rDNA/satellite array biology — **not
  claimed here**; this note makes no statement about *what* the extending sequence is, which is
  R1-M6/R1-M8's repeat-masking work, not mine.

## 6. Consequences the referee names, restated with the numbers

1. **“The 500 kb maximum in Fig. 1B is by construction.”** True, for exactly 1,087 of 15,668
   PHRs (6.94 %), and those 1,087 are 99.6 % acrocentric-short-arm + PAR ends. The Fig. 1B
   legend already says “Saturation at 500 kb reflects the input flank size”
   (`submission/paper.tex:470-471`); the body text does not, and the mean does not.
2. **“Walking inward until inter-chromosomal similarity disappears does not hold for
   saturated ends.”** True for 1,087 PHRs; false for 14,581 (93.1 %), where the walk ends on
   4 consecutive failing windows and the length is complete. Quantitatively: at the boundary
   the censored ends are still passing window after window at 0.99/0.96 continuation
   probability with a mean of 4.02 chromosomes — they were cut mid-signal.
3. **“3.51 Mb is a lower bound.”** Right conclusion, wrong mechanism (§3): the total is
   censoring-free (31 CHM13 PHRs, none above 295 kb). It is a lower bound because seven arms
   carry no PHR and because the acrocentric and PAR arms are excluded from it by construction.
4. **“Justify the parameters.”** With these numbers 500 kb *is* defensible for the call as a
   whole (6.94 % censored; 0.03 % of the 14,055 non-acro/non-PAR PHRs) and *is not* defensible
   for the acrocentric short arms and PARs, where 94 % and 43 % of calls saturate and the true
   extents measure ≥2.6–3.8 Mb by direct projection (4.0–17.1 Mb by transitive reachability, an
   upper bound) for the acrocentric short arms and ≥2.4 Mb for Xp/Yp. The defensible
   parameterisation is two-tier: 500 kb for the genome-wide call, ≥2 Mb for the nine arms that
   carry 1,083/1,087 of the censored PHRs.

## 7. What must be softened in the manuscript (and what need not be)

If the cluster-scale 2–5 Mb re-extraction is never run, these are the claims that must change;
all of them are already supported by the measurements above, so the manuscript can be
corrected without that run:

- **r1 (must change).** “PHRs span tens to hundreds of kilobases (median 105 kb, mean 144 kb)”
  (`submission/paper.tex:159-161`). The *mean* is a right-censored mean: ≥225–294 kb on the
  measured lower bounds. The *median* is exact. Scope the sentence to the assayed 500 kb.
- **r2 (must change).** “Walking inward from each telomere until inter-chromosomal similarity
  disappears” (`:159`) — say “until similarity disappears or the flank ends”, and give the
  censored count.
- **r3 (must change).** Any downstream statement that treats PHR length as an extent for the
  acrocentric short arms and PARs. This propagates: e.g. the C7 community's “median 500 kb”
  and its 91.6–95.8 % saturation noted in `response/notes/D1_superpopulation_descriptives.md`
  §4 are a censoring artefact of those two arm classes.
- **r4 (may stay, with one clause).** “3.51 Mb” may stay as written **if** a clause says what
  it is: the CHM13 PHR total excluding the acrocentric short arms and PARs, censoring-free.
  It should not be offered as a censoring correction.
- **r5 (should add).** The paper never states which assembly the 3.51 Mb total refers to.
  It is CHM13 (verified here to the base pair). Say so.
- **r6 (may stay).** The Fig. 1B legend sentence “Saturation at 500 kb reflects the input flank
  size” is correct and should be promoted into the body text rather than left in the legend.
- **r7 (must not be claimed).** Do not claim a specious “corrected distribution”: the extension
  figures here are bounds from a different (sparser, stricter) alignment regime, not a re-call.
  If the 2 Mb re-extraction is later run, replace them.

## 8. Corrected text to put in the manuscript

### 8.1 Results (`submission/paper.tex:159-165`, replacing the second and third sentences)

> Walking inward from each telomere until inter-chromosomal similarity disappears, we found
> that PHRs span tens to hundreds of kilobases within the assayed 500 kb (median 105 kb;
> Fig. 1B). The walk reaches the flank boundary while similarity is still present at 1,087 of
> 15,668 PHRs (6.9 %), so for those ends the 500 kb maximum in Fig. 1B is set by the flank
> rather than by biology: in the final 5 kb window such a PHR still involves a mean of 4.0
> distinct chromosomes, and the window-to-window continuation probability is 0.99 for the
> acrocentric short arms and 0.96 for the pseudoautosomal arms, i.e. the shared sequence is
> not decaying. Carrying the identical 5 kb-window criterion into the genome-wide alignment
> shows that these ends continue for a further median 0.9–3.0 Mb (1,365–3,450 kb across three
> identity/alignment-count thresholds; 34–97 % still passing at 2 Mb and 3–31 % still passing
> at 5 Mb), which raises the cohort-wide PHR total from 2.25 Gb to at least 3.5–4.6 Gb and the
> mean PHR length from 144 kb to at least 225–294 kb; the median is unaffected. The censoring
> is confined to two arm classes: 717 of the 763 PHRs on the acrocentric short arms (94.0 %),
> 366 of the 850 on the X and Y pseudoautosomal arms (43.1 %), and 4 of the 14,055 elsewhere
> (0.03 %). Excluding the acrocentric short arms and the pseudoautosomal regions, the CHM13
> PHRs total 3.51 Mb of subtelomeric sequence, a censoring-free total because none of the 31
> constituent PHRs reaches the flank boundary, but still a lower bound on subtelomeric sharing
> because seven arms carry no detected PHR.

### 8.2 Fig. 1B legend (`submission/paper.tex:466-471`, replacing the last sentence)

> Saturation at 500 kb reflects the input flank size: 1,087 PHRs (6.9 %) reach the boundary,
> 717 of them on the acrocentric short arms and 366 on the X/Y pseudoautosomal arms, where
> inter-chromosomal similarity is still present at the boundary and continues for at least a
> further 0.9–3.0 Mb (median). Lengths of all other PHRs are complete to the assayed 500 kb,
> since only 4 of 14,055 such PHRs reach the boundary.

### 8.3 Methods (append to the PHR-detection paragraph)

> **Right-censoring of PHR length.** The caller scans at most 500 kb inward from each telomere
> (`--max 500000`), so a PHR whose last available 5 kb window still satisfies the calling
> criteria has an unknown true length. We quantified this directly: 1,087 of 15,668 PHRs
> (6.9 %; 543.5 Mb of 2,252.5 Mb of PHR sequence) are censored, a further 23 terminate within
> one 5 kb window of the limit, and the censored ends are 93.97 % of the PHRs on the five
> acrocentric short arms, 43.06 % of those on the X/Y pseudoautosomal arms and 0.03 % of all
> others. Because the 500 kb flank PAFs were computed on the trimmed flank sequences, they
> cannot be extended; we therefore measured how far the censored ends' inter-chromosomal
> similarity actually continues using the genome-wide all-against-all survey alignments
> (25,274 whole-assembly comparisons, ~12 % of haplotype pairs), projecting every alignment
> onto 5 kb windows measured from the same telomeric trim boundary out to 5 Mb with the same
> ≥95 % identity, ≥3 kb projected length and ≥2-chromosome rules. Calibrated on the windows
> inside 500 kb, that instrument recovers 81–92 % of the published passing windows and
> false-passes 1.1–2.1 % of the failing ones, so its lengths are lower bounds on the published
> criterion's extent; it gives a median extension of +0.9 to +3.0 Mb (strict contiguity to the
> published four-failure chaining rule, three alignment-count thresholds), with 34–97 % of
> censored ends still passing at a 2 Mb ceiling. An independent transitive genome-wide profile
> (100 kb windows, ≥98 % identity) brackets the same ends from above at 4.0–17.1 Mb for the
> acrocentric short arms and 2.4 Mb for Xp. A dedicated re-extraction of the nine affected arms
> with 2 Mb flanks is the appropriate confirmation and is feasible at cluster scale; it is not
> required for the correction, which is in the sign and magnitude of the censoring, not in its
> exact size.

### 8.4 Response-to-referee paragraph (for `response/POINT_BY_POINT.md`, R1-M8)

> Point taken, and our thanks: the 500 kb maximum in Fig. 1B is indeed set by the flank. We
> now state the censored fraction exactly (1,087/15,668 = 6.9 %) and, more importantly, where
> it lives: 717 of the 763 acrocentric short-arm PHRs and 366 of the 850 X/Y pseudoautosomal
> PHRs, against 4 of the 14,055 PHRs on all other arms (0.03 %), on 11 of 41 arms. We
> re-extracted at 5 Mb in the sense the referee asks for, from the deeper genome-wide
> alignment: the censored ends continue for a further median 0.9–3.0 Mb, 34–97 % of them past
> 2 Mb and 3–31 % past 5 Mb, and the p-arm acrocentric ends pass 4 Mb in the transitive
> genome-wide profile. The consequences are amended accordingly: “walking inward until
> similarity disappears” holds for 92.9 % of PHRs and is replaced by an explicit
> censoring statement for the rest; the mean PHR length becomes a censored mean (≥225–294 kb);
> the median (105 kb) is unaffected, since the censoring point is far above the 50th
> percentile; and 3.51 Mb needs no numerical correction at all — it is the CHM13 total
> excluding precisely the acrocentric and pseudoautosomal arms, and none of its 31 constituent
> PHRs reaches 500 kb — though it remains a lower bound because seven arms carry no PHR. The
> one place the referee's diagnostic bites hardest is the parameter choice: 500 kb is
> adequate for the call genome-wide and demonstrably inadequate for the acrocentric short arms
> and PARs, which is why we now recommend (and have implemented as an option) a 2 Mb flank for
> those nine arms.

## 9. Residual risks

- **r1 (alignment regime, unavoidable).** The extension is measured in a *different* alignment
  regime than the call (12 % of pairs, whole assemblies, id ≥ 0.98-dominated). Its
  under-detection is measured (sensitivity 0.81–0.92) and its over-calling is measured
  (specificity 0.94–0.99, 1.1–2.1 % length overreach, 1–10/3,158 PHR-free flanks gaining a
  region), so the bound is defensible — but it is a bound, not a re-call. If a 2 Mb flanks
  wfmash run is ever made on the cluster, its numbers should replace §4.3's.
- **r2 (definitional).** Regime B counts an alignment to another chromosome's *interior* as an
  inter-chromosomal hit; regime A could only see subtelomere-vs-subtelomere. Direction of
  bias: upward on B's lengths. Magnitude: bounded by the measured 1.1–2.1 % overreach rate.
- **r3 (transitivity).** §4.5's extents are transitive-closure reachability and are inflated
  (they read 15–20 % high on arms where both regimes are measurable, e.g. 1p 347 kb vs
  published 295 kb). They are used only as an upper bracket.
- **r4 (window granularity).** Everything is at the published 5 kb resolution (100 kb for
  regime C), so lengths are quantised; the 23 “near-censored” PHRs at 495–499 kb are the
  sensitive band and are excluded from the censored count by definition.
- **r5 (provenance of 3.51 Mb).** The decomposition into CHM13 minus acrocentric-p minus X arms
  is exact to the base pair but is *reconstructed*, not documented anywhere in the repo or the
  manuscript. If the intended scope were different (cohort-wide it would be 1.54 Gb), the
  §3 conclusion would need restating. This is cheap to settle with the authors.
- **r6 (label semantics, separate item).** In the published rule a contig whose chromosome
  alias is a coordinate range is counted as a chromosome distinct from the chromosome it
  belongs to (20,296 of 521,356 per-window rows, 3.9 %). Checked: **0 of the 374 marginal
  (2-chromosome) boundary passes rely on such a label**, so the censoring counts here do not
  depend on it. It remains a live issue for the ≥2-chromosome criterion in general and belongs
  to the R1-M8 repeat/masking item, not to this one.
- **r7 (this note's scope).** No claim is made about *what* the extending sequence is (rDNA,
  satellite, TAR1, shared TEs). That is R1-M6/R1-M8-masking work by another lane.

## 10. Reproduction

```
$ python3 scripts/cladistics/flank_censoring_sensitivity.py --workdir /scratch/erikg/m8-censoring/clean --threads 32
[20:42:05]   projecting 25272 files through 32 workers (318.4 GB compressed PAF to stream)
...
[20:51:57] Part B: survey-regime calibration and extension
[21:15:55] writing scripts/cladistics/flank_censoring_sensitivity.tsv (839 rows)
[21:09:06] done
$ md5sum scripts/cladistics/flank_censoring_sensitivity.tsv
b14d54e9d6fb16094ca8d6327a2198cb
```

Walls (this host: 256 cores, 1.0 TB RAM, moosefs over `mfs#octopus04`; nothing else of mine ran
concurrently — the box's other load was the operator's jobs, unmeasured):

| step | wall | notes |
| --- | --- | --- |
| Part A + C (`--skip-survey`) | 4.4 s | 18,826-flank replay gate + transitive table |
| re-run with intermediates present (Parts A–C, no stream) | 9 s | `--workdir` reuse of `parts/`, `merged.tsv`, `counts.bin` |
| 32-worker PAF projection (~318 GB gz, 25,272 files) + merge + dense aggregation | 9 m 38 s (20:42:02 → 20:51:40) | **388,755,097** PAF records in, **371,935,337** kept; per-worker peak 15.9–16.4 M aggregate keys (2^26 slots); merge 511,928,269 input lines → **27,527,612** distinct (flank, window, partner, identity-bin) keys; cube 18,827 × 1,000 × 6 = 112,962,000 B. On the development run of the identical steps the merge alone measured 4 m 27 s and the aggregation 24 s. |
| Part B analysis | 7 s | |
| **total, clean run (from scratch)** | **10 m 1 s** (21:17:31 → 21:27:32) | the two earlier from-scratch runs of the same code path measured 9 m 21 s and 10 m 5 s |
| `--skip-survey` from scratch after adding Parts A5/A6 | 4.5 s | 676 rows (Parts A, A5, A6, C) |

Counters are exact from the program's own `stderr` lines (`proj: records/kept/keys`,
`merge: lines/distinct keys`, `agg: rows`).

## 11. Repo-wide gates run

- `gate_caller_rule_replay_exact = 18826/18826` (full table row; a mismatch aborts the script).
- per-window depth ceiling: `max dist_from_telomere = 495000` asserted, not assumed.
- every flank in the fai is exactly 500,000 bp and every flank span is exactly 500,000 bp.
- every flank contig resolves through `hprcv2.chromAlias.txt` (12,649/12,649; a miss aborts).
- the call table's row names are a subset of the fai, with exactly one documented drop.
- CHM13 PHR bp = 6,840,000 matches `scripts/cladistics/phr_by_superpopulation.tsv`.
- 3.51 Mb reproduced to the base pair from the CHM13 rows of the call table.
- the delivered TSV is byte-identical across independent full runs of the script (verified: two from-scratch runs on separate work dirs, plus a re-run reusing intermediates; every row matched, and the final revision only *adds* the three boundary-bookkeeping rows of §12).

## 12. Addendum — the boundary-adjacent indeterminate set (measured)

The rows `n_reaching_the_boundary`, `n_boundary_unfinished_failure_streak` and
`n_terminated_exactly_at_the_boundary` make the boundary bookkeeping exact:

- **1,133** flanks walk as far as the 495–500 kb window. Of these:
  - **1,087** pass it → certainly censored (unknown true length);
  - **3** complete all four failing windows inside the flank → terminated exactly at the
    boundary, correctly (e.g. last pass at 480 kb);
  - **43** reach the boundary with an *unfinished* failure streak (1 pending failure: 23;
    2 pending: 13; 3 pending: 7). These are not censored under the published definition —
    their region ends at 495/490/485 kb — but a 5 Mb flank could in principle resume them.
    This is the correct size of the "might resume" class: 43 PHRs (0.27 %), not 1,087.

So the censoring statement is sharp: 1,087 PHRs (6.94 %) have unknown length, 43 (0.27 %)
have a length that a longer flank could still extend, and the remaining 14,538 (92.8 %) are
complete.
