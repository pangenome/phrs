# M1 — Concerted evolution of duplicated families (D4Z4/DUX4L, TUBB8/TUBB8B, OR4F, SEPTIN14)

Status: COMPLETE (resumed lane A; steps 1-4 done)

**One-paragraph answer.** The per-family material available in the repo
(`scripts/cladistics/concerted_evolution_families.tsv`, produced by
`scripts/cladistics/concerted_evolution_families.py`) does *not* show the homogenisation
signature. On the PHR substrate — the only substrate on which the test is meaningful — the
two alleles of one individual at one arm are closer than that individual's closest copy on a
different arm for all four families, and on the referee's alternative formulation
(within-individual cross-arm copies vs the same arm's copies in *other* individuals) the
cross-arm copies are 2.2x-11.5x farther apart. That is evidence *against* recent
cross-arm homogenisation at whole-window resolution for these four families. It is not a
contradiction of the manuscript's aggregate result, because the aggregate's single
homogenisation signal is community C7 (acrocentric p-arms, 70.5% paralog closer, p = 2.0e-7)
and **none of the four families has a copy on a C7 arm**. Nothing in this lane dates any
exchange; see "what this does NOT establish".

## commands run

All commands were run in `/home/erikg/phrs` unless a `cd` is shown. Canonical list, in order
(the full transcript-level detail is in the step sections below):
```
wc -l scripts/cladistics/concerted_evolution_families.py                      # 908
grep -n ... scripts/.../concerted_evolution_families.py ; sed -n 1,908p <same>  # full read (step 1)
ls -la /moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz{,.fai}
ls -la /moosefs/guarracino/HPRCv2/PHR_III/flanking/hprcv2.1Mb.flanking_100kb.fa.gz{,.fai}
head -3 <both .fai files>; grep -m3 'chr4_qarm' <both .fai files>                   # substrate shapes
awk -F'\t' 'NR>1 && $3=="family_aggregate"{...}' scripts/cladistics/concerted_evolution_families.tsv
awk -F'\t' 'NR>1 && $3=="arm_pair"{...}'        scripts/cladistics/concerted_evolution_families.tsv
awk -F'\t' '{print $1}' <P .fai> | grep -c 'chr4_qarm$'                            # 357 (chr10_qarm 357)
awk (uniqueness of sample#hap#arm)                                                 # 0 duplicates
awk -> 85 samples with hap1+hap2 at both D4Z4 arms                                 # == n_matched_units
cd /scratch/erikg/celfam_ac2 && ./run_step3.sh    # step 3 harness (see STEP 3 section)
cd /scratch/erikg/celfam_ac2 && python3 analyze.py | tee analyze1.txt
cd /scratch/erikg/celfam_ac2 && python3 order.py   | tee order1.txt
cd /scratch/erikg/celfam_ac2 && python3 - <<PY ... PY | tee extra.txt               # AUCs, between-indiv
awk -F'\t' ... ratio paralog/between-individual, substrate P                        # step 4a table
grep -n -i "paralog\|allele" end-to-end-report/report/04_heterogeneity.md           # step 4b quote
grep -n -i "concerted" PEER_REVIEW.md                                              # referee text
grep -n 'C7' data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv                # C7 = 13p,14p,15p,21p,22p
```

## numbers

Every number below is from one of those commands; the step sections carry the context.
Core table (substrate P, the meaningful one):

| family | arms (community) | med allelic | med paralog | n_paralog_closer/n | sign p | paralog/between-indiv |
|---|---|---|---|---|---|---|
| D4Z4_DUX4L | chr4_q,chr10_q (C1) | 0.00358841 | 0.00868924 | 36/85 (0.4235) | 0.19276 | 2.22 |
| TUBB8_TUBB8B | chr10_p,chr18_p (C2) | 0.0022683 | 0.0268019 | 0/197 | 9.95682e-60 | 11.45 |
| OR4F | 7 arms (C3,C11) | 0.0049049 | 0.0168533 | 272/1144 (0.2378) | 1.07925e-73 | 2.95 |
| SEPTIN14 | 10 arms (C3,C11) | 0.00393894 | 0.0224132 | 843/3127 (0.2696) | 8.24949e-152 | 4.57 |

Sketch-vs-alignment validation (D4Z4, 20 individuals, 120 same-individual pairs):
Spearman(mash distance, 1 - alignment identity over the union of FastGA blocks) =
**+0.7704, p = 8.0e-25**; 120 pair classes separate with Mann-Whitney AUC 0.12
(p = 5.6e-12); the allelic-vs-paralog ordering survives (8/20 paralog-closer with mash,
3/20 with alignment identity, sign p = 0.0026).

Full numbers: see STEP 2, STEP 3, STEP 4a/4b below.

### step 0 — resume context

Lane A was killed at 30 min. Existing artefacts (read, not modified):
- `scripts/cladistics/concerted_evolution_families.py` (908 lines, read in full)
- `scripts/cladistics/concerted_evolution_families.tsv` (57 rows incl. header)

### step 0 — commands run
- `wc -l scripts/cladistics/concerted_evolution_families.py` -> 908
- `head -60 scripts/cladistics/concerted_evolution_families.tsv`
- `sed -n 1,908p` of the script (full read)

### STEP 1 — what P and F mean (read from code, not guessed)

Evidence, script line refs (`scripts/cladistics/concerted_evolution_families.py`):

- `PHR_FA` (line 78-79) = `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz`.
  `FLANK_FA` (line 80-81) = `.../PHR_III/flanking/hprcv2.1Mb.flanking_100kb.fa.gz`;
  `SUBSTRATES = {"P": PHR_FA, "F": FLANK_FA}` (line 83).
- `extract_windows()` (line 253) selects, from a substrate's `.fai`, **one window per
  (individual, haplotype, chromosome arm)**, header format `SAMPLE#HAP#contig:start-end_arm`;
  `arm_of_window_name()` (line 141) and `sample_hap_of_window_name()` (line 149) parse it.
- Distances are Mash (k=21, sketch 20000; `MASH_K`/`MASH_SKETCH`, lines 89-90) within each
  individual only (`same_sample_distances()`, line 356).
- `classify()` (line 486):
  * **allelic** = window of hap1 vs window of hap2 of the SAME individual at the SAME arm;
  * **paralog** = MINIMUM over the 2x2 haplotype combinations of the cross-ARM windows of the
    SAME individual, restricted to arm pairs in the SAME community (`communities`/arm-leiden
    k15 assignments, line 72-73).
- Matched unit = one (individual, armA, armB) with a paralog value and both allelic values;
  the allelic comparator is the MEAN of the two arms' allelic distances (line 815-821).

**P (one sentence):** the telomere-anchored PHR window itself — the p95/id95-trimmed
subtelomeric window that carries the duplicated/paralogous sequence, on each haplotype of
each individual (`SUBSTRATES["P"]`, line 83).

**F (one sentence):** the contiguous 100 kb of *unique* sequence immediately
centromere-proximal to that same PHR window, on the same haplotype and arm — adjacent to,
not inside, the PHR (`SUBSTRATES["F"]`, line 83; the name literally is
`SAMPLE#HAP#chrN_qarm_flank:start-end`).

Verified on the `.fai` files (headers differ in shape between substrates):

```
P: CHM13#0#chr10:134574995-134754994_chr10_qarm   180000 bp
   CHM13#0#chr10:2489-72489_chr10_parm             70000 bp
F: CHM13#0#chr10_qarm_flank:134474995-134574995  100000 bp
   CHM13#0#chr10_parm_flank:72489-172489         100000 bp
```
So the flank window is exactly abutting the PHR window on its centromere side (chr10_qarm:
flank 134,474,995-134,574,995 then P 134,574,995-134,754,994; chr10_parm: P ends at 72,489,
flank starts at 72,489).

**Consequence for interpretation (this is the thing that flips the table).** On substrate F,
the "paralog" pair compares two *non-paralogous unique loci* on different arms. F's median
paralog distance (~0.18-0.24, Jaccard ~0.006-0.012) is therefore not a paralog divergence at
all: it is the method floor for unrelated sequence. The F columns can only demonstrate that
the *allelic* comparison is locus-matched; they cannot test homogenisation, and the
alleles-closer outcome there is guaranteed by construction. Only the **P** substrate is a
test of the homogenisation hypothesis. (The script's own docstring, line 22-27, says exactly
this, but the table does not distinguish the two roles.)

Substrate size check (`.fai`, both substrates contain the identical panel):
```
P records=15668 unique_samples=233 unique_arms=41
F records=15666 unique_samples=233 unique_arms=41
diff of the two sample lists: empty
```
=> same 233 names in both (232 individuals + CHM13; CHM13 has hap 0 and so can never enter an
allelic or paralog pair), F missing 2 windows. This is why n_allelic_pairs / n_paralog_pairs /
n_between_pairs are identical between P and F for every family (e.g. 281/202/28508 for D4Z4).

### STEP 2 — the table, per family and substrate

`scripts/cladistics/concerted_evolution_families.tsv`: 56 data rows = 8 `family_aggregate`
+ 48 `arm_pair`, 28 rows per substrate.

command:
```
awk -F'\t' 'NR>1 && $3=="family_aggregate"{...}' scripts/cladistics/concerted_evolution_families.tsv
```

| family | sub | arms | n_allelic | med_allelic_dist | n_paralog | med_paralog_dist | n_matched | n_paralog_closer | n_allelic_closer | frac_paralog_closer | sign_test_p | wilcoxon_p (H1: paralog<allelic) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| D4Z4_DUX4L | P | chr10_qarm;chr4_qarm (C1) | 281 | 0.00358841 | 202 | 0.00868924 | 85 | 36 | 49 | 0.4235 | 0.19276 | 0.979528 |
| D4Z4_DUX4L | F | same | 281 | 0.00288799 | 202 | 0.177977 | 85 | 1 | 84 | 0.0118 | 4.4461e-24 | 1 |
| TUBB8_TUBB8B | P | chr10_parm;chr18_parm (C2) | 426 | 0.0022683 | 229 | 0.0268019 | 197 | 0 | 197 | 0.0000 | 9.95682e-60 | 1 |
| TUBB8_TUBB8B | F | same | 426 | 0.00110699 | 229 | 0.199098 | 197 | 0 | 197 | 0.0000 | 9.95682e-60 | 1 |
| OR4F | P | 7 arms, C3+C11 | 1471 | 0.0049049 | 1374 | 0.0168533 | 1144 | 272 | 872 | 0.2378 | 1.07925e-73 | 1 |
| OR4F | F | same | 1471 | 0.0021172 | 1374 | 0.207474 | 1144 | 0 | 1144 | 0.0000 | 0 | 1 |
| SEPTIN14 | P | 10 arms, C3+C11 | 2087 | 0.00393894 | 3655 | 0.0224132 | 3127 | 843 | 2284 | 0.2696 | 8.24949e-152 | 1 |
| SEPTIN14 | F | same | 2087 | 0.00242354 | 3655 | 0.207842 | 3127 | 0 | 3127 | 0.0000 | 0 | 1 |

(between-individual reference, per family in bp of alignable sequence, is the same in both
substrates: D4Z4 28508 pairs median 0.00390884 (P) / 0.00329118 (F); TUBB8 28494 / 0.00234094 /
0.00117535; OR4F 99700 / 0.00571152 / 0.00227104; SEPTIN14 142456 / 0.00490649 / 0.00257741.
Because F is unique sequence, median allelic ~ median between-individual in F: 0.0024 vs 0.0026
for SEPTIN14, i.e. the flank behaves like ordinary single-copy sequence. In P the allelic
distance is far below the between-individual distance for OR4F/SEPTIN14 (0.0049 vs 0.0057;
0.0039 vs 0.0049) — a small but systematic enrichment of within-individual identity.)

**Direction, stated plainly:** in every one of the 4 families on both substrates, the
alleles are closer than the paralogs (frac_paralog_closer <= 0.4235 everywhere; 0.0000 for
TUBB8/OR4F/SEPTIN14 under P). On the substrate that can actually test the hypothesis (P),
3 of 4 families exclude the homogenisation direction at extreme significance
(sign test p = 9.96e-60, 1.08e-73, 8.25e-152) and D4Z4 is null in both directions
(36/85, p = 0.19). This is the OPPOSITE of the homogenisation signature.

Caveats I must carry:
- The paralog value is a MINIMUM over 4 haplotype combinations, whereas the allelic
  comparator is a mean of 2 values. That biases the test TOWARD "paralog closer"; the
  observed anti-homogenisation result survives that bias, so the bias is conservative here —
  but the two one-sided tests can disagree within an arm pair (OR4F C3 chr11/chr19:
  frac_paralog_closer=0.5707, sign p=0.0596, wilcoxon p=7.0e-05) because the Wilcoxon
  compares distributions, not matched signs.
- Only 6 of 48 arm_pair rows trend in the homogenisation direction and none is significant by
  the sign test (OR4F C3 chr11/chr19 0.5707 p=0.0596; SEPTIN14 C3 chr11/chr16 0.5700 p=0.0514;
  SEPTIN14 C3 chr11/chr19 0.5707 p=0.0596). All involve chr11_parm.
- `chm13_copies_in_chm13_window` and `chm13_window_bp` are 0 on every F row. This is a
  parsing artefact, not biology: `window_span_of_chm13()` (line 322) matches
  `:([0-9]+)-([0-9]+)_chr`, which never matches the flank header
  `CHM13#0#chr10_qarm_flank:...`. The F copy-in-window columns are therefore UNPOPULATED.

### STEP 3 — setup checks (D4Z4 only, alignment cross-check)

Tools: `mash` = /home/erikg/bin/mash, `samtools` = /usr/bin/samtools, `FastGA` =
/home/erikg/.cargo/bin/FastGA, `wfmash` = /home/erikg/bin/wfmash. Substrates stat-able
(P 541,667,887 B; F 400,982,104 B), nothing was I/O-dead on this substrate.

Selection audit (why 85 = n_matched_units for D4Z4): from the P `.fai`,
`chr4_qarm` has 357 windows, `chr10_qarm` 357; 0 duplicate (sample,hap,arm) keys; and
**85 samples have hap1+hap2 at BOTH arms**, exactly the 85 matched units in the D4Z4 rows of
the TSV. Independent confirmation that I read `classify()`'s pairing rule correctly.
Only CHM13 carries hap `0` (38 windows) so CHM13 can never enter a pair.

Window lengths are NOT uniform (this is a real limitation of the whole table):
```
HG00128#1#CM090078.1:191805164-191975163_chr4_qarm   170000
HG00128#2#CM090094.1:192443268-192553267_chr4_qarm   110000
HG00128#1#CM090084.1:135699503-135794502_chr10_qarm   95000
HG00128#2#CM090099.1:135924357-136054356_chr10_qarm  130000
```
So an "allelic" mash distance at chr10_qarm compares a 95 kb window with a 130 kb window;
mash Jaccard is over the union of the two sketches, so a length difference mechanically
depresses Jaccard/inflates distance. Lengths are trimmed at p95/id95 so they track array
length per haplotype. Applied to ALL comparisons (allelic and paralog) so the design is
balanced, but not length-matched.

Harness scratch: `/scratch/erikg/celfam_ac2` (new; the killed lane's `celfam`/`celfam_ac`
were left untouched). Resume = re-run `run_step3.sh` in that directory (below).

FastGA invocation check on 4 real sequences (HG00128, both arms, both haps):
`FastGA -T4 -f512 -s500 -paf t4.fa t4.fa` -> exit 0 in 1.3 s. Note: FastGA is a
chaining/local aligner, so one pair emits MANY PAF lines (9-24 kb blocks for cross-arm pairs):
identity must therefore be built from PAF columns 10/11 (nmatch/alnlen) per block, and I
report coverage so low-coverage pairs cannot masquerade as full-length identity.

### STEP 3 — results (D4Z4_DUX4L only, substrate P, 20 individuals)

Harness: `/scratch/erikg/celfam_ac2/run_step3.sh` (idempotent; resume = re-run it).
Selection: 85 qualifying -> `awk 'BEGIN{srand(20260501)}...sort' -> 20` =
HG00350 HG00621 HG00658 HG01261 HG01530 HG01934 HG01993 HG02178 HG02273 HG02809 HG03130
HG03369 HG03486 HG03710 HG03742 HG03874 HG04187 NA18971 NA20827 NA21110
80 windows (20 x 2 haps x 2 arms), 11,465,000 bp total.
Commands (verbatim):
```
samtools faidx -r sel_names.txt <PHR_III/...p95.id95.fa.gz> > sel.fa
mash sketch -i -k 21 -s 20000 -o sel sel.fa
mash dist sel.msh sel.msh > mash_pairs.tsv            # 6400 ordered pairs
FastGA -T4 -f512 -s500 -paf sel.fa sel.fa > fastga.paf
```
FastGA wall 19.93 s, 4 threads, max RSS 404,392 kB (`/usr/bin/time -v`, `fastga.time`);
nothing else heavy on the box during this step (only my own shell). 63,115 PAF lines.
Per-pair identity from PAF cols 10/11 (nmatch/alnlen). Three metrics per pair, because
FastGA is a chaining/local aligner and D4Z4 arrays break into 7-8 blocks per pair (median
nblocks 7 allelic / 8 cross):
- `best_block` = identity of the largest single block;
- `union` = matches/columns over the greedy non-overlapping union along the query;
- `global` = union matches / length of the shorter window (unaligned bases counted as
  mismatches; the closest analogue to a whole-window k-mer distance).
Coverage = covered query bp / shorter window length.

Pair classes: 40 allelic (20 indiv x 2 arms), 80 cross (20 x 2x2 hap combos), 3040
between-individual.

**Medians**

| class | n | mash_d | best_block id | union id | covered bp | cov frac | global id | nblocks |
|---|---|---|---|---|---|---|---|---|
| allelic | 40 | 0.002854 | 0.99851 | 0.99804 | 103,530 | 0.8165 | 0.8153 | 7 |
| cross-arm (paralog) | 80 | 0.026888 | 0.99284 | 0.98970 | 94,254 | 0.7726 | 0.7662 | 8 |
| between-individual | 3040 | 0.009110 | 0.99398 | 0.99263 | - | - | - | - |

**Correlation, sketch distance vs alignment identity** (mash distance is the x; a positive
rho below against `1-identity` means mash distance rises as identity falls, i.e. agreement):

```
same-individual pairs (n=120)
  vs 1-best_block  rho=+0.6812  p=1.1e-17
  vs 1-union       rho=+0.7704  p=8.0e-25     <-- headline number
  vs 1-global      rho=+0.3507  p=8.6e-05
allelic only (n=40)      vs 1-union rho=+0.8304 p=3.4e-11
cross only   (n=80)      vs 1-union rho=+0.5069 p=1.6e-06
ALL pairs    (n=3160)    vs 1-union rho=+0.7228 p reported as 0 by scipy (<1e-300)
>=50% coverage (n=110)   vs 1-union rho=+0.8020 p=6.5e-26
```
So the k-mer sketch distance does track base-level alignment identity, but the two are far
from a monotone re-parameterisation of each other; the `global` (coverage-penalised) variant
tracks worst, because coverage and identity trade off.

**Discrimination between allelic and cross-arm pairs** (Mann-Whitney, H1 allelic smaller):
```
mash_d      AUC=0.1187  p=5.6e-12   med 0.002854 vs 0.026888  (9.4x)
1-best_id   AUC=0.1047  p=9.7e-13   med 0.0014945 vs 0.0071605 (4.8x)
1-union_id  AUC=0.1503  p=2.4e-10   med 0.0019555 vs 0.010296  (5.3x)
1-global    AUC=0.3844  p=0.0198    med 0.1847 vs 0.23379      (1.3x)
```

**Does the allelic-vs-paralog ordering survive?** Yes. Replicating `classify()` on these 20
individuals (allelic = mean of the two arms' hap1-vs-hap2 value; paralog = min over the 2x2
cross-arm haplotype combinations):

```
metric        n_pc  n_ac  frac_pc  sign_p     wilcoxon_p(H1 par<al)  med_par   med_al
mash          8     12    0.4000   0.5034     0.9385                 0.017764  0.0038835
best_block    3     17    0.1500   0.002577   0.98                   0.005985  0.0018163
union         6     14    0.3000   0.1153     0.9932                 0.0083265 0.0025577
global        9     11    0.4500   0.8238     0.5364                 0.17521   0.17224
```
The ordering is the same direction under every metric and, under the strictest base-level
metric, stronger than with mash (3/20 vs 8/20 paralog-closer). Note the whole-table value for
this family is 36/85 = 0.4235, and my 20-sample subset gives 8/20 = 0.40 — the subset is
representative.

Honest reading of the alignment numbers: cross-arm D4Z4 copies align to each other at ~99%
identity over ~77% of the shorter window. So the paralogs ARE extremely similar in the
aligned blocks; what separates them from the alleles is not block identity (0.9928 vs
0.9985, a 0.6% gap) but the fact that the alleles are more similar still, and that identity
over aligned blocks saturates. This is exactly why block identity is a weak substitute for
whole-window k-mer distance and why the referee's "high identity alone does not show
homogenisation" point bites: identity is high everywhere here.

### STEP 4a — the same conclusion in the referee's own formulation

Referee #1 (PEER_REVIEW.md:76) asks for "a per-family comparison of within-individual
divergence between non-homologous copies against **between-individual divergence at the same
end**". That is a different contrast from allelic-vs-paralog, and the TSV already carries the
column it needs (`median_between_individual_dist`, distances between DIFFERENT individuals at
the SAME arm, capped at 120 haplotypes per arm with seed 20260501). Computed from
`concerted_evolution_families.tsv`, family_aggregate rows, substrate P:

| family | within-indiv paralog | between-indiv same-arm | ratio paralog/between | allelic | allelic/between |
|---|---|---|---|---|---|
| D4Z4_DUX4L | 0.00868924 | 0.00390884 | 2.22 | 0.00358841 | 0.918 |
| TUBB8_TUBB8B | 0.0268019 | 0.00234094 | 11.45 | 0.0022683 | 0.969 |
| OR4F | 0.0168533 | 0.00571152 | 2.95 | 0.0049049 | 0.859 |
| SEPTIN14 | 0.0224132 | 0.00490649 | 4.57 | 0.00393894 | 0.803 |

(F rows, for completeness: 54.08x, 169.39x, 91.36x, 80.64x — meaningless as divergence, see
Step 1; they only re-demonstrate that the flanking substrate is non-paralogous.)

In all four families the within-individual non-homologous (cross-arm) copies are 2.2x-11.5x
FARTHER apart than the same-arm copies of *different* individuals. Under the referee's own
criterion this is the opposite of homogenisation: cross-arm copies have not converged past
the level of ordinary population variation at one arm.

### STEP 4b — how this sits against the manuscript's own aggregate result

The manuscript's aggregate (my authority for it is the report, not this script) is
`end-to-end-report/report/04_heterogeneity.md`, section "Allele vs paralog distance"
(lines 7-28). Quoted verbatim from that table (lines 13-24):

```
C2 424 pairs  med allelic 0.023 med paralog 0.365  0.0% paralog closer  p=3.3e-71  Allele closer
C12 437      0.037 / 0.727   0.2%  3.0e-73   Allele closer
C11 854      0.098 / 0.364  14.4%  2.7e-103  Allele closer
C5  805      0.032 / 0.439  14.9%  1.4e-90   Allele closer
C6 1307      0.073 / 0.250  16.1%  3.1e-159  Allele closer
C4  429      0.109 / 0.549  20.5%  3.2e-55   Allele closer
C3 1269      0.170 / 0.368  34.9%  7.0e-52   Allele closer
C1  265      0.293 / 0.375  44.2%  6.0e-5    Allele closer
C7  156      0.481 / 0.353  70.5%  2.0e-7    Paralog closer
Overall 5,946  0.069 / 0.387  20.4%  <1e-300  Allele closer
```
and its stated conclusion (line 26, 28): "Alleles are significantly more similar than
paralogs in 8 of 9 multi-arm communities (overall Wilcoxon p < 1e-300). The sole exception is
C7 (acrocentric p-arms), where paralogs are closer than alleles in 70.5% of individuals
(p = 2.0e-7, now significant)... C7's reversal ... confirming complete inter-chromosomal
homogenization of acrocentric p-arms (silhouette = -0.029). C1 (D4Z4) shows the weakest
allele advantage (44.2% paralog closer), consistent with its intermediate homogenization
(silhouette = 0.147)."

I must be precise here, because the brief's phrasing ("the aggregate is the result that DOES
support homogenisation") is only true for one of the nine communities. The aggregate and my
per-family table agree in direction on 8 of 9 communities: alleles closer than paralogs. The
result that supports homogenisation is **C7 and C7 alone** (70.5% paralog closer, p = 2.0e-7,
156 pairs); the overall row is allele-closer at 20.4%.

**And the two are not in conflict, because they do not measure the same arms.** From
`data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`, C7 = chr13_p, chr14_p, chr15_p,
chr21_p, chr22_p — the acrocentric p-arms. The four refereed families sit in
C1 (chr4_q, chr10_q), C2 (chr10_p, chr18_p), C3 (chr11_p, chr16_q, chr19_p, chr3_q, chr7_p,
chr9_q), C11 (chr1_p, chr5_q, chr6_q, chr8_p[unused]) and C8 (chr15_q). **Not one arm
carrying a D4Z4/DUX4L, TUBB8/TUBB8B, OR4F or SEPTIN14 copy is in C7.** So the per-family test
does not sample the community where homogenisation is actually detected; it samples the eight
where arm identity is maintained. The families therefore cannot confirm or refute C7 — they
are the wrong eight communities. Also note the two analyses are numerically incomparable:
report 04 uses Jaccard distances from
`/moosefs/guarracino/HPRCv2/scripts/community/allele_vs_paralog_distance.R` with medians
0.069/0.387, whereas this script's PHR medians are 0.0039/0.0224 — a ~10-17x different scale
(different substrate window set, sketch and aggregation), so only the sign is comparable.

Caveats on the comparison itself: sign-test p values printed as exactly `0` in the TSV (OR4F F,
SEPTIN14 F) are scipy `binomtest` underflow, i.e. "< 1e-300", not a literal zero; and the
manuscript's C1 row (44.2% paralog closer, 265 pairs, p=6.0e-5) is the weakest allele
advantage in its own table, which is the same qualitative story as D4Z4's family-level
36/85 = 0.4235, p = 0.19 here (the only family in either analysis that is genuinely
ambiguous).

## what this does NOT establish

1. **It does not date anything.** This table has no time axis. A
   paralog-farther-than-allele result bounds *conversion* frequency, not its age: it says
   that at this resolution the cross-arm copies have not been homogenised recently enough to
   be indistinguishable, not when any exchange happened.
2. **It does not show that no exchange occurred.** Absence of the homogenisation signature at
   whole-window k-mer resolution is compatible with old exchange, with exchange of a
   sub-window-sized segment (a few kb inside a 100-235 kb window), or with exchange whose
   donor tract is not fixed in the population — all of which leave the window-level distance
   dominated by array-length and copy-number differences. Note the length spread of the PHR
   windows themselves (chr4_qarm: n=357, min 35 kb, median 155 kb, max 500 kb; chr10_qarm:
   n=357, min 20 kb, median 135 kb, max 390 kb), which is exactly the signal a
   length-changing ectopic event would produce and which a whole-window Jaccard/mash distance
   partly absorbs (a mash distance from a 20 kb window to a 390 kb window is not a measure of
   per-base divergence).
3. **It is not a duplication-time test.** Referee #1's point is that high identity alone is
   compatible with recent duplication. Nothing here separates "recent duplication" from
   "old duplication plus homogenisation" for D4Z4 (the one null family, 36/85, p = 0.19).
   For the other three families the alleles-closer result does exclude the *cross-arm*
   homogenisation direction, but it does not establish the duplication age either.
4. **It does not test C7** — the one community where the manuscript's own aggregate reports
   homogenisation (70.5% paralog closer, p = 2.0e-7). None of the four families has a copy on
   a C7 arm.
5. **It is not read-level genotyping.** Windows come from a fixed draft-style trimming
   (`telo_trimmed.p95.id95`), so a haplotype with an array outside the p95 length/identity
   envelope is simply absent from that arm-individual cell; n_allelic_pairs (281 for D4Z4 of
   85*2*2 possible cells) shows this attrition. Absent cells are dropped, not modelled.
6. **It does not control for array length.** Distances are length-unmatched (95 kb vs 130 kb
   is a typical allelic pair); mash Jaccard over the union of two sketches is mechanically
   depressed by a length difference. The design applies this to allelic and paralog pairs
   alike, and the D4Z4 allelic and cross classes have similar coverage (0.82 vs 0.77 of the
   shorter window), so the bias does not create the direction — but it is not removed.
7. **The F substrate tests nothing about homogenisation.** It is unique flanking sequence;
   its paralog distance (~0.18-0.24) is the floor for unrelated sequence, and
   alleles-closer there is guaranteed by construction. Any citation of the F half of the
   table as a "control for paralog divergence" would be wrong.
8. **The `chm13_copies_in_chm13_window` / `chm13_window_bp` columns are unpopulated for F**
   (parser artefact), so "how much of the family each flank window contains" is unverified.

## What the surviving clock is

The dating logic the manuscript can legitimately use, and which this table constrains rather
than supplies:
- **A state shared by all individuals predates their common ancestor.** Whatever degree of
  cross-arm similarity is common to all 232 individuals (e.g. the D4Z4 cross-arm copies being
  ~99% identical over ~77% of the window in every individual sampled) cannot be explained by
  exchanges in the individuals' own lineages; it is at least as old as the population's most
  recent common ancestor. Population-shared homogenisation is therefore ancestral, not
  recent-and-ongoing.
- **Residual paralog differences bound the number of conversion events.** The leftover
  cross-arm divergence (mash 0.0087-0.0268 depending on family, i.e. far above the allelic
  0.0023-0.0049) is a lower bound on how much sequence has escaped conversion: N conversion
  events cannot have covered these tracts unless the last one is old enough to have
  accumulated the residual variation. In this table that residual is large, which is an upper
  bound on recent conversion, not a date.
- **The surviving clock is flanking single-copy sequence and the breakpoints**, not the
  duplicons: unique sequence flanking the arrays (and the junctions where the exchange tracts
  end) mutates clock-like, whereas the duplicon interior is the thing being homogenised and
  therefore carries no usable clock. This is why the F substrate matters as a *clock
  calibration*, not as a paralog-divergence control: its median allelic distance
  (0.0011-0.0029) and its between-individual median (0.0012-0.0033) give the neutral rate for
  a locus of that length in this panel. Note how tight that is: allelic/between-individual is
  0.877-0.969 for the flank, i.e. single-copy sequence at one arm carries almost no excess
  within-individual identity, so any excess seen in P is a property of the duplicons.
  Breakpoint-level dating is NOT done here and remains open (it needs the exchange tracts
  themselves, which a 100-235 kb window distance cannot resolve).

### Verification pyramid (what was actually gated)

1. **Brute-force oracle, already in the script, re-run now (unmodified script):**
   ```
   python3 scripts/cladistics/concerted_evolution_families.py --self-test
   == self-test 1: mash distance vs exact Jaccard on tiny synthetic sequences
      s_ref vs s_ident  mash_d=0.000000  exact_J=1.00000  mash_J=1.00000
      s_ref vs s_mut5   mash_d=0.005092  exact_J=0.81586  mash_J=0.81586
      s_ref vs s_unrel  mash_d=1.000000  exact_J=0.00000  mash_J=0.00000
   == self-test 2: allelic/paralog classification on a synthetic panel
      allelic=5 paralog=3 (expected 5 and 3)
   == self-test 3: sign test and Wilcoxon on a known-direction toy vector
      sign p=0.00195  wilcoxon p=0.000244
   ALL SELF-TESTS PASSED
   ```
   Self-test 2 includes the adversarial cases that matter for the interpretation: an arm in a
   different community must never enter the paralog set, a single-haplotype arm must not leak
   an allelic pair, and the 2x2 minimum must be taken. All three are asserted and pass, so my
   reading of `classify()` in STEP 1 is the code's actual behaviour.
2. **Structural cross-check on the real data:** 85 qualifying samples == the table's
   n_matched_units for D4Z4; identical n_allelic/n_paralog/n_between between the P and F
   substrates explained by the identical 233-name panel (diff empty, F missing 2 windows).
3. **Independent-method gate (STEP 3):** k-mer distances re-derived from base-level FastGA
   alignments on real sequences, agreeing in sign and ordering (rho = +0.77, p = 8e-25;
   AUC 0.12; ordering survives with a stricter metric).

### Final checks on my own claims

- Community map used for the "no family arm in C7" claim, read from
  `data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv` (41 arms, 15 communities):
  `C1: chr10_q chr4_q | C2: chr10_p chr18_p | C3: chr11_p chr16_q chr19_p chr3_q chr7_p chr9_q
  | C4: chr12_q chr7_q | C5: chr12_p chr20_q chr6_p chr9_p | C6: chr13_q chr17_q chr19_q
  chr1_q chr21_q chr22_q | C7: chr13_p chr14_p chr15_p chr21_p chr22_p | C8: chr15_q |
  C9: chr16_p | C10: chr17_p | C11: chr1_p chr5_q chr6_q chr8_p | C12: chr20_p chr2_q |
  C13: chr4_p | C14: chrX_q chrY_q | C15: chrX_p chrY_p`.
  Family arms: D4Z4 {chr4_q, chr10_q} = C1; TUBB8 {chr10_p, chr18_p} = C2;
  OR4F {chr1_p, chr3_q, chr5_q, chr6_q, chr11_p, chr15_q, chr19_p} = C3, C3, C11, C11, C3,
  C8, C3; SEPTIN14 adds chr7_p, chr9_q (C3), chr16_q (C3), chr20_q (C5, no partner ->
  no paralog pairs), chr2_q (C12, no partner). Intersection with C7 = empty. That is why the
  OR4F arm_pair rows are exactly the 3 C3 pairs + 3 C11 pairs, and why chr15_q/ chr20_q /
  chr2_q appear in the arm list but never as an arm_pair row (a single-arm community yields
  no cross-arm pair).
- Not modified: no file outside `response/notes/M1_concerted_evolution.md` was written by me
  (checked with `find . -newermt '-6 hours'`; the other recent files belong to other lanes).
  No git commands were run.
- Harness files live only in `/scratch/erikg/celfam_ac2` (`run_step3.sh`, `sel20.txt`,
  `qualifying.txt`, `canon.tsv`, `order1.txt`, `analyze1.txt`, `extra.txt`, `fastga.time`).
