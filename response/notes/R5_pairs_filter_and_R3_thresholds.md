# R5 L193–195 (the 2,830 PHR pairs) and R3 major 6 (PHR-call thresholds)

Items: **R5-L193** ("why only 2,830 pairs of 15,668 PHRs") and **R3-M6** ("41 of 48
arms may depend on the selected thresholds") from `response/POINT_BY_POINT.md`.

Authority documents read: `PEER_REVIEW.md` (lines 133 and 160 hold the two verbatim
comments), `response/POINT_BY_POINT.md` (lines 41 and 175 hold the two TODOs),
`submission/paper.tex` (Methods, lines 620–665), `scripts/phr-detection/PIPELINE.md`
(command-level record of the PHR call), `scripts/phr-detection/
find-multichr-regions-incremental.py` (the calling code).
`paper_prep/synthesis/SIMULATED_PEER_REVIEW_v1.md` was **not** opened.

Artefacts produced by this item:

| path | what it is |
| --- | --- |
| `scripts/cladistics/phr_call_threshold_sensitivity.py` | re-runs the PHR window-walking logic over the published per-window IMPG table and sweeps every threshold that table lets us vary (with a hard replication gate) |
| `scripts/cladistics/phr_call_threshold_sensitivity.tsv` | 21-row arm-count table across settings |
| `response/notes/R5_pairs_filter_and_R3_thresholds.md` | this note |

---

## 0. Inputs, and which numbers are borrowed vs re-derived

Throughout, **[quoted]** = copied from a named file/line; **[re-derived]** =
produced by a command I ran here, in this worktree.

| input | provenance | md5 / size |
| --- | --- | --- |
| `data/all-vs-all.1Mb.p95.id95.len.tsv` | per-flank PHR call of the published run, 18,826 data rows | `47cbaca27807a61a12f64606a9209864` |
| `data/all-vs-all.1Mb.p95.id95.len.per_window.tsv` | **not vendored**; read from `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.per_window.tsv`, 18,827 flanks / 521,456 windows / 93 MB | `5260c41c47659250168e4ba69fed365c` |
| `/moosefs/guarracino/HPRCv2/scripts/community/sequence_hic_correlation.py` | the code that builds every Fig. 4 / ED1 PHR pair | `88bada2a798298e43d996ece5cd3a1ea` |
| `data/human_seqlevel_sweep/human_HG002_porec_{5,10,20,50,100}000bp_seqlevel.tsv` | per-PHR-pair tables behind Fig. 4A | in-repo |
| `data/human_seqlevel_sweep/human_HG002_hic_*`, `human_HG002_cifi_*`, `human_CHM13_hic_*` | per-PHR-pair tables behind ED1 | in-repo |
| `data/chm13.phrs.bed`, `data/chm13.phrs.no_acro.bed` | CHM13 subset BEDs | `29809854…`, `631c795e…` |

The repo does **not** contain `data/all-vs-all.1Mb.p95.id95.len.per_window.tsv`
(the task listed it as a repo-local input; it is not there). It is readable on
`/moosefs` and I streamed it read-only; the sweep script takes `--per-window
PATH` and aborts with a fetch instruction if the file is absent, so it cannot
silently fall back to something else.

---

# Part 1 — Referee #5, L193–195

## 1.1 Verbatim (`PEER_REVIEW.md:160`)

> L193-195: The authors measure Jaccard similarity and contact strength at 2,830
> pairs of inter-chromosomal PHRs. If the authors identified over 15,000
> candidate PHRs, shouldn't there be far more than 2,830 possible pairs? Or was
> this collection of 2,830 pairs filtered in some way (that I may have missed)?

## 1.2 Short answer

Nothing about similarity or contact was filtered out. The filter is **which two
objects are being paired**: the 15,668 PHRs are the *called regions across 465
haplotypes* (one row per telomere-anchored flank), whereas a Fig. 4A point is a
pair of *PHR sequences belonging to one single sample*, because a Hi-C/Pore-C
bin must contain both members of the pair. HG002 contributes 77 PHR sequences
(41 arms × ~2 haplotypes, minus 5 arms with a single haplotype); C(77,2) = 2,926
unordered pairs, of which **96 are on the same chromosome and are skipped**,
leaving **2,830**. That is the whole of the arithmetic.

## 1.3 Where the pairs are built

`sequence_hic_correlation.py` (moosefs, md5 above) is the only code that writes
these tables (`end-to-end-report/report/05_hic_validation.md:153` and
`data/human_seqlevel_sweep/README.md` both name it). Relevant lines:

- `L272–279` sequence universe filtered by `--sample` and by `--min-region-bp` (default 5000)
- `L327–331` coordinates clamped to the cooler's chromosome sizes; sequence dropped if `end-start < min_region_bp`
- `L346` similarity matrix loaded, restricted to the surviving sequences
- `L360–366` **pair enumeration**: `for j in range(i+1, n)` over the sorted surviving sequence list; `if chr_a == chr_b: continue` → same-*base*-chromosome pairs skipped
- `L376` `if n_bins == 0: continue` → pairs whose every cooler bin is ICE-masked are dropped
- `L381–384` contact normalised by the number of valid bins; Jaccard taken from the similarity matrix with **default 0.0** when the pair is absent (i.e. no similarity filter)
- `L387–402` one output row per surviving pair

## 1.4 The filter chain, in order

1. **Sample filter.** Only HG002's own PHR sequences enter (`--sample HG002`).
   The 77 unique sequences in the Pore-C table are all `HG002#1#…` or
   `HG002#2#…` **[re-derived]**, so this flag was necessarily used; without it
   the universe is 15,668 sequences and the table would have ~10^8 rows.
2. **Region-size floor.** ≥ 5 kb per sequence (code default). Every queried
   region in the published table is ≥ 9,999 bp **[re-derived]**, consistent with
   the default; the run's own log is not retained, so this flag value is
   *inferred from the code default*, not read from a log.
3. **Cooler membership.** Sequences must map to a chromosome the cooler actually
   has, after clamping to chromosome length.
4. **Inter-chromosomal only.** Pairs on the same *base* chromosome are skipped —
   including the paternal-vs-maternal pair of the same chromosome, because the
   comparison is on base chromosome, not haplotype-labelled chromosome.
5. **At least one valid (non-ICE-masked) bin-pair.** Dropped pairs are "missing
   data, not no contact" (comment at `L375`).
6. **No other filters.** Zero Jaccard is kept (626 of the 2,830 Pore-C 50 kb rows
   have Jaccard exactly 0) **[re-derived]**; zero contact is kept (1,585 of 2,830)
   **[re-derived]**. There is no p-value, no top-N, no uniqueness/dedup step.

## 1.5 The arithmetic

**[re-derived]**, two independent ways, both giving exactly 2,830:

(a) from the published variant table `data/all-vs-all.1Mb.p95.id95.len.tsv`
(PHR-bearing flanks only, grouped by sample, base chromosome parsed from the
`_chrN_parm|_qarm` suffix):

```
HG002: nseq=77 all_pairs=2926 same_base_chr=96 inter=2830
```

(b) from the Fig. 4A input table itself
`data/human_seqlevel_sweep/human_HG002_porec_50000bp_seqlevel.tsv`: 2,830 data
rows, 77 unique sequences (39 paternal + 38 maternal), 41 distinct arms
(36 arms with both haplotypes, 5 with one) — so 2926 − 96 = 2830, and nothing
was dropped at step 5 for Pore-C.

Cross-checks that this is the same 2,830 as in the manuscript:
`submission/paper.tex:285` reports "pointwise Spearman ρ = 0.376, n = 2,830";
re-running Spearman on the vendored 20 kb table gives **n = 2,830, ρ = 0.3760**
and on the 50 kb table **n = 2,830, ρ = 0.3814** (the latter is the
`end-to-end-report/report/05_hic_validation.md:165` value) **[re-derived]**.

Genome-wide context **[re-derived]**: over all 233 samples in the published PHR
call there are **503,302** inter-chromosomal PHR pairs in total (520,936 less
17,634 same-base-chromosome pairs); HG002's 2,830 is 0.56 % of that. Pairing all
15,668 PHR sequences with each other, as the referee's phrasing invites, would
give C(15,668,2) = **122,735,278** pairs — but those mix sequences from 465
different assemblies, and no contact map contains two sequences from different
assemblies, so that set is not physical.

## 1.6 Why n differs across the Fig. 4A / ED1 panels

**[re-derived]**, by enumerating each panel's own sequence universe and its
same-base-chromosome pairs:

| panel dataset | file rows | unique seqs | all pairs | same-chr | inter-chr | rows in file | dropped (no valid bin) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| HG002 Pore-C 20 kb (Fig. 4A) | 2,830 | 77 | 2,926 | 96 | 2,830 | 2,830 | 0 |
| HG002 Pore-C 50 kb | 2,830 | 77 | 2,926 | 96 | 2,830 | 2,830 | 0 |
| HG002 Hi-C 20 kb | 2,413 | 77 | 2,926 | 96 | 2,830 | 2,413 | **417** |
| HG002 CiFi 50 kb | 2,757 | 76 | 2,850 | 93 | 2,757 | 2,757 | 0* |
| CHM13 Hi-C 20 kb | 688 | 38 | 703 | 15 | 688 | 688 | 0 |
| CHM13 Hi-C 50 kb | 652 | 37 | 666 | 14 | 652 | 652 | 0* |

\* one whole PHR sequence is absent from these tables (CiFi: `HG002#1#chr13:
16777-516777_chr13_parm`; CHM13 50 kb: one of the 38), so the file's own
universe is already one sequence smaller and its pair count is the complete
enumeration of that smaller universe.

So the four n in the referee's Fig. 4A comment (2,830; 2,413; 2,757; 688) are
three different things: the same 2,830-pair universe thinned by masked bins in
Hi-C (2,413), a 76-sequence universe (2,757), and a different sample altogether
(CHM13, 38 PHR sequences → 688). They are not evidence of additional filtering.

## 1.7 What the manuscript already says, and the one clause to add

`submission/paper.tex:284–286` already says "per inter-chromosomal PHR pair
**within a single genome** in HG002 Pore-C (Fig. 4A; pointwise Spearman ρ = 0.376,
n = 2,830)". "Within a single genome" is the answer the referee missed; what is
missing is the count. **Suggested addition** (proposal only, no new number
beyond those re-derived above):

> Each Fig. 4 and ED1 point is one unordered pair of that sample's own PHR
> sequences (HG002: 77 PHR sequences on 41 arms, i.e. two per arm except five
> arms with a single haplotype; CHM13: 38). Pairs of PHRs on the same chromosome
> are excluded (96 of the 2,926 HG002 pairs), and a pair is dropped only if every
> cooler bin spanning it is masked; no similarity or contact filter is applied,
> so pairs with zero shared graph nodes (626/2,830 in Pore-C) and with zero
> contact (1,585/2,830) remain. This is why n differs between panels: 2,830
> (Pore-C, complete enumeration), 2,413 (HG002 Hi-C, 417 pairs lost to masked
> bins), 2,757 (CiFi, one PHR sequence absent from that map) and 688 (CHM13,
> 38 PHR sequences).

## 1.8 Adjacent facts from the same referee comment (zero/floor points)

The same referee's Fig. 4A comment (`PEER_REVIEW.md:95`) asks about the flat line
of points on the vertical axis and about differing n. **[re-derived]** zero
counts at the canonical 50 kb setting:

| dataset | pairs | Jaccard = 0 | normalised contact = 0 |
| --- | --- | --- | --- |
| HG002 Pore-C 50 kb | 2,830 | 626 | 1,585 |
| HG002 CiFi 50 kb | 2,757 | 620 | 2,468 |
| HG002 Hi-C 20 kb | 2,413 | 519 | 286 |
| CHM13 Hi-C 20 kb | 688 | 149 | 85 |

The vertical-axis pile is the 1,585 Pore-C pairs with zero length-normalised
contact (drawn at the plot floor, `make_fig4a_human_scatter.R:51–52`,
`y_plot <- pmax(y, floor_y)`). Note the CiFi tie problem: 2,468/2,757 (89.5 %)
of CiFi pairs are at the contact floor, which is why its ρ = 0.191 is weaker than
Pore-C's. (Full treatment of that comment is a different item; recorded here
because it is the same table and the same arithmetic.)

---

# Part 2 — Referee #3, major point 6

## 2.1 Verbatim (`PEER_REVIEW.md:133`)

> 6. Page 3, second paragraph; Methods pages 14–15. The conclusion that PHRs
> occur on 41 of 48 arms may depend on the selected thresholds: 95% identity, at
> least five alignments from at least two chromosomes, 5-kb windows, a 500-kb
> flank, and four consecutive failed windows.

## 2.2 The thresholds as actually written in Methods

**[quoted]** `submission/paper.tex`:

- `L624–630`: telomere detection with `seqtk telo`; contigs ≥ 1 Mb classified by
  arm; "At every telomere-bearing end we trimmed the telomeric repeat and
  retained the 500 kb immediately internal to it as the flank, so a
  centromere-spanning contig yields both a p-arm and a q-arm flank, producing
  18,827 candidate flanks across the 48 arms."
- `L630–634`: one flank excluded (NA18982 hap1 chr18q scaffolding chimera) →
  "leaving 18,826 flanks for the all-against-all alignment and PHR calling."
- `L636–640`: "We aligned all 18,826 flank sequences all-against-all with wfmash
  v0.23.0 at 95% minimum identity … without restricting comparisons to
  homologous chromosomes."
- `L643–648`: "Starting at each telomere, we queried 5 kb windows in 5 kb steps,
  stopping after four consecutive windows failed the detection criteria. A
  window passed when the direct IMPG query returned at least five alignments
  from each of at least two chromosomes, including its own (so pseudoautosomal
  Xp/Yp and Xq/Yq, with a single partner, are also detected), at a minimum
  identity of 95% and minimum output length of 3 kb."
- `L649`: "This defined 15,668 PHRs from the 18,826 retained flanks (83.2%)."
  (15668/18826 = 83.24 % **[re-derived]**.)
- `L151–153`: headline: "15,668 PHRs on 41 of the 48 chromosome arms; the
  remaining seven arms (chr2p, chr3p, chr5p, chr8q, chr11q, chr14q, chr18q)
  contained no PHR."

**[quoted]** `scripts/phr-detection/PIPELINE.md` Stage 4 gives the flag values of
the production run, which is the authoritative parameter list:

```
--window 5000 --step 5000 --max 500000 \
--min-identity 0.95 --min-output-length 3000 \
--min-diff-chrs 2 --min-count 5 --min-consecutive-stop 4
```

Two thresholds the referee did **not** list are in the code and matter:
`--min-output-length 3000` (each projected alignment must be ≥ 3 kb) and the
`--min-count 5` being a *per-chromosome* alignment count evaluated independently
for each chromosome (the referee states this correctly as "at least five
alignments from at least two chromosomes", but it is worth separating: ≥5 is per
chromosome, ≥2 is the number of such chromosomes).

## 2.3 What can and cannot be re-derived from local inputs

The published PHR call is a walk over per-window IMPG query results. The
per-window table survives (521,456 rows; columns `n_chrs`, `n_arms`,
`chrs_involved`, `arms_involved`, `passes_threshold`, `dist_from_telomere`), so
the *walk* can be replayed exactly and three of the six thresholds can be varied
post hoc:

| threshold | variable locally? | why |
| --- | --- | --- |
| ≥ 2 chromosomes per window (`--min-diff-chrs`) | **yes** — `n_chrs` column | post-hoc re-threshold of a stored count |
| 4 consecutive failed windows (`--min-consecutive-stop`) | **yes, stricter only** (1, 2, 3) | the production scan stopped at the 4th failure, so no window beyond it was ever queried; loosening to >4 would need re-querying |
| 500 kb flank / inward scan limit (`--max`) | **yes, shorter only** (25 kb … 400 kb) | `dist_from_telomere` is stored; windows beyond 500 kb do not exist |
| 95 % identity (`--min-identity`, and `wfmash -p 95`) | **no** | below 95 % the alignments were never produced by wfmash; a different identity needs a fresh all-vs-all. Rebuilding the IMPG index is the blocker: the 18,827 all-vs-all PAFs are 88 GB on `/moosefs` and the `paf.list` points at `/scratch` paths that no longer exist. A subset index does **not** work either: `query_impg_window()` counts alignments **per target chromosome**, so with only some targets indexed every chromosome's count falls below 5 and all windows fail — the count is a whole-pangenome quantity. |
| ≥ 5 alignments per chromosome (`--min-count`) | **no** | the table stores *chromosomes that already passed the ≥5 rule* (`n_chrs`, `chrs_involved`), not raw per-chromosome alignment counts; the 5 can only be moved by re-querying IMPG (same 88 GB blocker) |
| 5 kb window / 5 kb step (`--window/--step`), 3 kb minimum output length (`--min-output-length`) | **no** | window geometry and the projection-length floor are baked into the query that wrote the table |

`wfmash` and `impg` binaries *are* installed here (`~/.cargo/bin/impg`,
`~/bin/wfmash`), so the limitation is the missing 88 GB intermediate + index
build, not the tools. This is stated as prescribed: the identity, window-size and
alignment-count axes are **out of scope for this item** and any claim about them
must be scoped accordingly.

## 2.4 Arm-count table

`scripts/cladistics/phr_call_threshold_sensitivity.tsv` (21 rows). Every row
carries the fixed values it did not vary (`min_count_alignments_per_chr = 5`,
`min_identity = 0.95`, `window_bp = 5000`, `min_output_len_bp = 3000`) so the
baked-in parameters are visible in the data itself. Arm count is "arms bearing at
least one PHR" over the 48-arm universe; `non_acro_of_38` removes the ten arms
chr13p/q, chr14p/q, chr15p/q, chr21p/q, chr22p/q, matching what
`data/chm13.phrs.no_acro.bed` actually drops relative to `data/chm13.phrs.bed`
**[re-derived]**.

| variant | value | PHRs | % of 18,826 flanks | arms (of 48) | non-acro (of 38) | median PHR | total Mb | arms gained / lost |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PUBLISHED_TABLE | — | **15,668** | 83.23 | **41** | 32 | 105,000 | 2,252.5 | — |
| REPLICATION_FROM_PER_WINDOW | — | 15,668 | 83.23 | 41 | 32 | 105,000 | 2,252.5 | — (gate) |
| min-diff-chrs | 1 | 18,825 | 99.99 | **48** | 38 | 90,000 | 2,606.0 | +chr2p,chr3p,chr5p,chr8q,chr11q,chr14q,chr18q |
| min-diff-chrs | 3 | 14,369 | 76.33 | **40** | 31 | 85,000 | 1,716.6 | −chr12q |
| min-diff-chrs | 4 | 14,357 | 76.26 | 40 | 31 | 70,000 | 1,643.5 | −chr12q |
| min-diff-chrs | 5 | 14,344 | 76.19 | 40 | 31 | 45,000 | 1,437.0 | −chr12q |
| min-diff-chrs | 6 | 13,422 | 71.30 | 38 | 29 | 60,000 | 1,170.4 | −chr10p,chr12q,chr18p |
| min-diff-chrs | 8 | 11,523 | 61.21 | 38 | 29 | 60,000 | 985.5 | −chr10p,chr12q,chr18p |
| min-diff-chrs | 10 | 11,516 | 61.17 | 38 | 29 | 60,000 | 977.1 | −chr10p,chr12q,chr18p |
| min-diff-chrs | 12 | 11,511 | 61.14 | 38 | 29 | 30,000 | 388.2 | −chr10p,chr12q,chr18p |
| min-diff-chrs | 15 | 10,030 | 53.28 | 37 | 29 | 15,000 | 132.6 | −chr10p,chr12q,chr15p,chr18p |
| consecutive-stop | 1 | 15,624 | 82.99 | **41** | 32 | 105,000 | 2,139.2 | — |
| consecutive-stop | 2 | 15,664 | 83.20 | **41** | 32 | 105,000 | 2,206.0 | — |
| consecutive-stop | 3 | 15,668 | 83.23 | **41** | 32 | 105,000 | 2,222.5 | — |
| flank limit | 25 kb | 15,668 | 83.23 | **41** | 32 | 25,000 | 368.6 | — |
| flank limit | 50 kb | 15,668 | 83.23 | **41** | 32 | 50,000 | 667.0 | — |
| flank limit | 100 kb | 15,668 | 83.23 | **41** | 32 | 100,000 | 1,147.4 | — |
| flank limit | 150 kb | 15,668 | 83.23 | **41** | 32 | 105,000 | 1,494.3 | — |
| flank limit | 250 kb | 15,668 | 83.23 | **41** | 32 | 105,000 | 1,894.8 | — |
| flank limit | 350 kb | 15,668 | 83.23 | **41** | 32 | 105,000 | 2,080.9 | — |
| flank limit | 400 kb | 15,668 | 83.23 | **41** | 32 | 105,000 | 2,138.7 | — |

ARMS count includes chrYp/chrYq PHR-bearing flanks from HPRC assemblies (CHM13 Y
is excluded upstream, `paper.tex:622–623`), so the 48-arm universe is 24
chromosomes × 2 arms regardless of which assembly supplies the flank.

## 2.5 Reading of the table

1. **The consecutive-failed-window rule is not load-bearing.** 1, 2, 3 and 4
   consecutive failures give 15,624 / 15,664 / 15,668 / 15,668 PHRs on **41 arms
   in all four cases**; only the PHR *extent* changes (total length 2,139 →
   2,253 Mb).
2. **The 500-kb flank is not load-bearing for the arm count at all.** Truncating
   the inward scan at 25 kb leaves the PHR count and the arm count *identical*
   (15,668 on 41 arms) and only shortens regions (median 105 kb → 25 kb; total
   2,253 → 369 Mb). Every PHR in this dataset is detected within the first 25 kb
   of the flank.
3. **The two-chromosome requirement is the load-bearing one, and it is a
   qualitative filter, not a tuning knob.** Dropping it (`min-diff-chrs = 1`,
   i.e. accepting a window with any signal at all) makes all **48** arms qualify
   — and the seven arms gained are exactly the seven the paper names as
   PHR-free. The claim "PHRs on 41 of 48 arms" is therefore a statement about
   *inter-chromosomal* sharing, and hinges on this criterion rather than on any
   continuous parameter.
4. **The arm count is flat from 2 to 5 chromosomes and degrades slowly above
   that**: 41 (k=2), 40 (k=3,4,5), 38 (k=6–12), 37 (k=15). Doubling the
   stringency from 2 to 4 costs exactly one arm (chr12q).
5. **The PHR count, by contrast, does move** — 15,668 → 14,344 at k=5 (−8.4 %) —
   so a reply should distinguish the stable arm-level conclusion from the
   threshold-sensitive region count.

Robustness statement that follows from 1–4: the arm-level conclusion ("nearly all
chromosome arms, and specifically 41 of 48") is invariant to the two thresholds
that shape region *extent* (consecutive-failure count, flank length) and is
stable within ±1 arm for a 2.5-fold change in the chromosome-count requirement;
what it is not invariant to is the requirement that sharing be inter-chromosomal
at all.

## 2.6 The identity axis — no number, and one rejected surrogate

No arm count as a function of identity can be produced here (see 2.3: 88 GB
all-vs-all + IMPG index rebuild).

I considered and **rejected** a surrogate: the independent CHM13 whole-genome
scan at 98 % identity (`data/hprc25272-wf.CHM13.100kb-xm5-id098-l50000.tsv.gz`,
100 kb windows) might have been read as "the signal survives at 98 %". It does
not measure the PHR rule: its terminal 100-kb bin is largely telomeric repeat —
**[re-derived]** 13 of the 48 terminal bins (the q-end bins of chr3, chr5, chr6,
chr9, chr11, chr12, chr16, chr18, chr19, chr20, chr22 and both ends of chrY) list
no chromosome at all in `chroms-num_haplotypes`, and where the PHR call does find
a region the 100-kb bin does not coincide with it. Using it would have
produced a number about the telomeric tip, not about the PHR call, so no
identity number is reported.

## 2.7 Suggested wording for the reply (proposal)

> We re-derived the arm count directly from the per-window IMPG table of the
> published call, replaying the detection walk exactly (it reproduces all 18,826
> called intervals and the 15,668 / 41-arm totals). Varying the two thresholds
> that govern region *extent* leaves the arm count unchanged: 1, 2, 3 or 4
> consecutive failing windows all give 41 arms, and truncating the flank from
> 500 kb to 25 kb also gives 41 arms (only the region lengths shrink, from a
> median 105 kb to 25 kb). Raising the number of chromosomes required per window
> from 2 to 3, 4 and 5 costs a single arm (chr12q, 41 → 40); 6–12 chromosomes
> give 38 arms and 15 gives 37. Removing the inter-chromosomal requirement
> altogether makes all 48 arms qualify. The count is therefore not sensitive to
> the window-support, flank or stopping parameters, and loses at most one arm
> under a 2.5-fold increase in the chromosome requirement; the identity and
> per-chromosome alignment-count thresholds enter before the per-window table is
> written and were not varied.

---

# 3. Gates

## Gate A — headline numbers reproduced from the vendored per-flank table

```
$ awk -F'\t' 'NR>1 && $4!="."{c++; n=split($1,p,"_"); arm=p[n-1] (p[n]=="parm"?"p":"q"); arms[arm]=1} END{print "regions="c", arms="length(arms)}' data/all-vs-all.1Mb.p95.id95.len.tsv
regions=15668, arms=41
```
and the seven PHR-free arms are exactly `chr2p chr3p chr5p chr8q chr11q chr14q
chr18q` **[re-derived]**, matching `paper.tex:151–153` **[quoted]**.

Independent length cross-check, same table, same walk **[re-derived]**: the
15,668 called intervals have median 105,000 bp and mean 143,765 bp, against
`paper.tex:159–160` "median 105 kb, mean 144 kb" **[quoted]**. This is a second,
independent quantity agreeing with the published text (the first being the
15,668 / 41 totals), which is why I regard the replication gate as tight rather
than merely count-matched.

## Gate B — full replay of the detection walk (in the script)

```
$ python3 scripts/cladistics/phr_call_threshold_sensitivity.py
reading per-window table /moosefs/.../all-vs-all.1Mb.p95.id95.len.per_window.tsv
  18827 flanks, 521456 windows
reading published per-flank call data/all-vs-all.1Mb.p95.id95.len.tsv: 18826 flanks
replication gate PASSED: 15668 PHRs on 41 arms, all 18826 region_start/region_end/arm sets identical to the published table
...
wrote scripts/cladistics/phr_call_threshold_sensitivity.tsv (21 rows)
```

The gate is exact, not statistical: for all 18,826 flanks, the replayed
`region_start`, `region_end` and the set of `arms_involved` equal the published
values, and the chimera flank named in `PIPELINE.md` Stage 4 is the only flank in
the per-window table absent from the published summary.

**Adversarial checks of the gate** (the gate must be able to fail):

```
# (i) perturb one published region_end by +5000 bp
$ python3 scripts/cladistics/phr_call_threshold_sensitivity.py --flank-summary /tmp/gt/tampered.tsv
MISMATCH ('CHM13#0#chr10:134254995-134754994_chr10_qarm', 'coords', (320000, 500000), (320000, 505000))
replication FAILED: 1 of 18826 flanks disagree with /tmp/gt/tampered.tsv at the published settings
PYTHON_EXIT=1

# (ii) point --per-window at a file that does not exist
$ python3 scripts/cladistics/phr_call_threshold_sensitivity.py --per-window /tmp/does_not_exist.tsv
per-window table not found: /tmp/does_not_exist.tsv
It is not vendored in the repo; re-fetch it with ... and pass --per-window.
```

The oracle here is not a small brute-force case but something stronger: an
independent stored product of the upstream job (the per-window table) replayed
into the upstream job's other stored product (the per-flank table). Agreement on
all 18,826 intervals, plus the independently agreeing length distribution
(median 105 kb / mean 143.8 kb) and PHR/arm totals, is what licenses the sweeps
below.

## Gate C — Fig. 4A pair count from an independent table

```
$ awk -F'\t' '
  NR>1 && $4!="." {
    k=split($1,p,"_"); chr=p[k-1]; nh=split($1,h,"#"); sm=(nh>=3?h[1]:$1)
    if(!(sm in n)){n[sm]=0; order[++ns]=sm}
    n[sm]++; ch[sm,n[sm]]=chr
  }
  END{
    for(i=1;i<=ns;i++){ s=order[i]; m=n[s]; a=0; b=0
      for(x=1;x<=m;x++) for(y=x+1;y<=m;y++){ a++; if(ch[s,x]==ch[s,y]) b++ }
      tot+=a-b; intra+=b
      if(s=="HG002") printf "HG002: nseq=%d all_pairs=%d same_base_chr=%d inter=%d\n",m,a,b,a-b
    }
    printf "samples=%d total_all_pairs=%d same_base_chr_pairs=%d inter_chromosomal_pairs=%d\n",ns,tot+intra,intra,tot
  }' data/all-vs-all.1Mb.p95.id95.len.tsv
HG002: nseq=77 all_pairs=2926 same_base_chr=96 inter=2830
samples=233 total_all_pairs=520936 same_base_chr_pairs=17634 inter_chromosomal_pairs=503302
```
agreeing with the row count of `human_HG002_porec_50000bp_seqlevel.tsv` (2,830)
and with the manuscript's n = 2,830 (`paper.tex:285`) **[quoted]**.

---

# 4. What this does NOT establish

1. **Nothing about the 95 % identity threshold.** No arm count vs identity is
   reported. The claim "the arm count is robust" is limited to the
   chromosome-count, consecutive-failure and flank-length axes.
2. **Nothing about the ≥5-alignments-per-chromosome threshold.** It is constant
   at 5 in every row; the `min_count_alignments_per_chr` column says so
   explicitly. A referee demanding "what if 3 alignments?" cannot be answered
   from local data.
3. **Nothing about window size (5 kb / 100 kb) or the 3 kb minimum projection
   length.** Both are baked in; this is Referee #5's L86 item, not solved here.
4. **No joint variation.** Each row moves one threshold from the published
   baseline (2 chromosomes, 4 failures, 500 kb). Interactions are unmeasured;
   e.g. 3 chromosomes *and* 250 kb was not computed.
5. **The re-thresholding is post hoc, not re-detection.** Raising `min-diff-chrs`
   filters stored window-level chromosome counts that were themselves computed
   under `--min-count 5 --min-identity 0.95 --min-output-length 3000`. It answers
   "what would the arm count have been if the caller had required more
   chromosomes", not "what does a fresh run at that setting give".
6. **Arm count is a union over 465 haplotypes**, not a per-assembly property.
   "41 of 48 arms" means "at least one flank somewhere bears a PHR on that arm".
   The table's `n_phrs` is likewise a flank count, not a locus count.
7. **The 2,830 is HG002's pair count at one contact-map resolution.** It is not a
   count of all PHR pairs in the study (503,302 genome-wide, re-derived), nor a
   count of pairs in ED1 (which are subsets/different samples).
8. **I did not re-run wfmash, pggb or the IMPG queries.** All PHR numbers are
   derived from the stored per-window and per-flank tables of the one published
   call.
9. **The `--sample`/`--min-region-bp` flags of the Fig. 4A run were not read from
   a retained log.** `--sample` is forced by the data; `--min-region-bp` is
   assumed at its code default because no log survives. If the production run
   overrode the default, the sequence *universe* would change and the 2,830 would
   not be reproduced — but it is reproduced exactly, which is consistent with the
   default.

# 5. Residual risks and provenance findings

1. **`data/chm13.phrs.bed` is not from the published PHR call.** Its 37
   intervals match the *earlier* run `/moosefs/.../all-vs-all.p95.id95.len.tsv`
   (18,226 flanks → 15,088 regions), not the published
   `data/all-vs-all.1Mb.p95.id95.len.tsv` (18,826 flanks → 15,668 regions).
   **[re-derived]** CHM13 has 38 PHR intervals in the published table, and six of
   the 37 vendored intervals differ from it, e.g. chr14 ends at 372,075 in the
   BED vs 502,075 in the published call; chr8 at 11,881 vs 31,881. Scripts that
   consume `data/chm13.phrs.bed` (Fig. 2b/2c heatmaps, Fig. 5 supporting material,
   `scripts/pedigree/crossover_rate_f34.py`, the permutation/shuffle controls)
   therefore use a slightly different PHR set from the headline 15,668. Worth
   flagging in the rebuttal's data-provenance answer (Referee #5's L152/L193
   questions are provenance questions) — and it may explain any small numeric
   drift between figures and text.
2. **`data/chm13.phrs.no_acro.bed` drops every interval on the acrocentric
   chromosomes (both arms of 13/14/15/21/22), not "the acrocentric short arms and
   the pseudoautosomal regions"** as `paper.tex:162–164` words it. The wording and
   the artefact agree in spirit but not literally; my `non_acro_of_38` column
   follows the artefact.
3. **`PIPELINE.md` says "15,668 PHRs (83.2% of the 18,827 flanks)"** while the
   paper says "from the 18,826 retained flanks". 15,668/18,826 = 83.24 % is the
   consistent arithmetic; `PIPELINE.md`'s 18,827 is the pre-chimera-removal count.
   Documentation nit only.
4. **The `paf.list` of the PHR call references `/scratch` paths** that no longer
   exist, so the identity-axis rerun would also need path rewriting on top of the
   88 GB index rebuild.
5. **`/moosefs` was read repeatedly without I/O errors in this session** (the
   per-window table, the calling script, the PHR directory listings, sizes). Two
   files the task expected in the repo (`data/all-vs-all.1Mb.p95.id95.len.
   per_window.tsv`) are absent and were read from `/moosefs` instead.
