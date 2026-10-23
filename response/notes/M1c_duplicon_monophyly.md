# M1c — Duplicon monophyly test (Referee #1 major point 1, in the correct form)

**Item.** Referee #1, major point 1. Verbatim: *"Concerted evolution is never tested. It requires that
paralogous copies at different ends be more similar than their divergence times predict, i.e. actively
homogenized; high identity alone does not show this, since recent duplication gives the same signature.
The data allow a direct test: duplicon phylogenies (D4Z4/DUX4L, TUBB8/TUBB8B, OR4F, SEPTIN14) asking
whether copies cluster by end or interdigitate across ends within a haplotype; … "* (`PEER_REVIEW.md`,
Referee #1 major point 1).

This analysis supersedes the distance-based attempt in `response/notes/M1_concerted_evolution.md`
(allele-vs-paralog mash distances; could not resolve the array interior because it worked at
whole-window resolution).

Files written by this lane (only these):
`scripts/cladistics/duplicon_monophyly.py`, `.tsv`, `.json`, `response/notes/M1c_duplicon_monophyly.md`.
No existing file was modified; `git` was not run.

---

## 0. Design, pinned before any number

### 0.1 What the topology can and cannot establish
- **Reciprocal monophyly of the arm classes over the whole haplotype panel** ⇒ the copies on
  arm A and the copies on arm B last shared an ancestor *before* the sampled haplotypes' common
  ancestor. No homogenisation between those arms within the sampled genealogy; the
  pattern is explained by ancient duplication + drift.
- **Arm classes para-/polyphyletic while haplotypes are monophyletic** ⇒ copies track the
  haplotypes, i.e. the copies on different arms within one haplotype are more similar to each
  other than to their orthologues on other haplotypes. That is the referee's "interdigitate across
  ends within a haplotype", and it is the topology that a between-arm homogenisation process
  produces. It is *also* what a single recent cross-arm duplication followed by drift produces —
  the topology cannot date the event, and cannot by itself separate "one recent duplication" from
  "repeated conversion". Both are excluded by reciprocal monophyly, neither is distinguished from
  the other by paraphyly.
- Neither outcome dates an exchange, numbers exchanges, or identifies a direction. A shared
  homogenised state bounds the last event as older than the sample's common ancestor only when
  the homogenised state is shared by all sampled haplotypes.

### 0.2 Substrates and the coordinate facts that drove the choice
- **Copies substrate:** the per-haplotype telomere-anchored windows
  `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz` (+`.fai`), read-only.
  *PHR-agnostic* on purpose: no ascertainment on PHR calling, all 48 arms present, and the same
  window definition for every haplotype.
- **Flank substrate (control):** `/moosefs/guarracino/HPRCv2/PHR_III/flanking/hprcv2.1Mb.flanking_100kb.fa.gz`
  (+`.fai`), read-only. The flanks are the 100 kb *centromere-proximal* neighbours of the PHR
  windows (e.g. `CHM13#0#chr4_qarm_flank:193292741-193392741`, i.e. ending exactly where the CHM13
  4q PHR window begins at 193392741). They contain no D4Z4 array (array = chr4:193389387-193553139,
  GFF3), so they cannot carry the copies; they are used as the (c) control, see §3.
- **Sample labels:** `slides/v2-review-zoom/_revision_assets/v9/labels_superpop_stats_polish/nearest_same_superpop_mds_distances.tsv`
  (233 individuals / 465 haplotypes, superpopulations AFR/AMR/EAS/EUR/SAS).

REPLACE-ME: inventory of which windows exist per arm (coverage), copy counts, etc.

## 1. Journal (append after every major step)

- [step 1] notes skeleton created.
- [step 2] read `PEER_REVIEW.md` (Referee #1 major point 1 verbatim, above) and
  `response/POINT_BY_POINT.md` R1-M1 section (`L2/L3 partial`; the aggregate railway is community
  **C7** with paralog closer in 70.5 % of individuals, p=2e-7 — this is the internal positive
  control for intermixing).
- [step 3] Inventoried inputs. Facts taken from files, not assumed:
  - CHM13 4q D4Z4 array (GFF3 `data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz`): DUX4L9 chr4:193389387-193390883
    through DUX4 chr4:193541579-193553139; tandem DUX4L/DUX4 genes spaced ~3.0-3.3 kb, ~33 annotated
    genes inside the CHM13 4q window (chr4:193072741-193572740).
  - CHM13 10q: DUX4L28 chr10:134615846-… inside window chr10:134254995-134754994.
  - TUBB8 chr10:39416-69557 (10p; window chr10:2489-502488); TUBB8B chr18:200997-227325
    (18p; window chr18:2017-502016). No gene named plain "TUBB8B" existed in the earlier
    `gene_family_copy_counting_analysis.md` grep; it is annotated in the GFF3.
  - OR4F3 chr5:181905487-181933378 (5q window 181543925-182043924); OR4F7P chr6:172013163-172014100
    (6q window 171623883-172123882).
  - SEPTIN14 subtelomeric copies: SEPTIN14P22_5/SEPTIN14P13 chr1:144134(+) / 146340 (1p);
    chr2:242631950 (2q); SEPTIN14P3 chr3:200959893 (3q); SEPTIN14P5 chr5:181897110 (5q);
    SEPTIN14P18 chr7:56811 (7p); chr9:150558546 (9q).
  - Communities (`data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`): C1 = 4q,10q; C2 = 10p,18p;
    C11 = 1p,5q,6q,8p; C3 = 3q,7p,9q,11p,16q,19p; C7 = 13p,14p,15p,21p,22p.
  - **rDNA is NOT in the telomere-anchored windows.** CHM13 RNA18SN1 sits at chr13:11731724-11733592
    (GFF3), i.e. 11.7 Mb from the telomere; the acrocentric PHR blocks in `data/chm13-annotations.bed`
    span chr13:9146474-10165439, chr13:10768726-12311210, chr13:12418709-13562256. No telomere-anchored
    window (max 500 kb) can reach them. Hence the C7 positive control uses the acrocentric p-arm
    repeat family that *is* in the windows, with the reference unit chosen by a documented,
    data-driven rule (§0.3), and the rDNA-positive-control alternative is reported as unavailable
    in these substrates with the coordinates above. REPLACE-ME with the rule as implemented.
  - Window coverage per arm in the 500 kb substrate (from the `.fai`): chr4_q 357, chr10_q 358,
    chr10_p 443, chr18_p 446, chr5_q 447, chr6_q 440, chr1_p 443, chr13_p 76, chr14_p 229,
    chr15_p 119, chr21_p 116, chr22_p 223 haplotypes; zero duplicate (haplotype, arm) keys in either
    substrate, so one window = one locus.
  - Tool availability: `FastGA` (`~/.cargo/bin/FastGA`), `mash`, `samtools`; no `mafft`, no `minimap2`.
    Python 3.12 + numpy 1.26.4 + scipy 1.11.4 + biopython 1.87.
- [step 4] FastGA syntax verified empirically (values are appended, not space-separated):
  `FastGA -pafS -l500 -i0.7 -f512 -T16 ref.fa chm13_4q.fa` → 35 alignments in 1.0 s wall for one
  500 kb window; `-f8192` → also 35. Alignments land one per D4Z4 unit, spaced ~3.3 kb, spanning
  chr4:193389626-193542863, i.e. exactly the annotated array. Long-form `cs` tag is emitted, so the
  projection to a reference frame is exact.
- [step 5] Script written (`scripts/cladistics/duplicon_monophyly.py`). Design decisions taken
  and why:
  - Copy substrate = the **PHR-agnostic** per-haplotype 500 kb telomeric windows
    (`hprcv2.1Mb.telo_500kb_trimmed.fa.gz`), not the PHR-trimmed set: the PHR-trimmed set carries
    PHR-calling ascertainment and only 76/229/119/116/223 haplotypes have the acrocentric p windows,
    whereas the raw 500 kb set has one uniform window definition per (haplotype, arm) with zero
    duplicate keys.
  - Copies are discovered by **homology**, not by annotation, so array-interior copies are found:
    FastGA (reference gene unit vs every window), `-l500 -i0.7 -pafS`, units filtered to
    >= 0.80 identity and >= 500 bp, de-duplicated at < 20 % target overlap, capped at 10 copies per
    window (evenly spaced along the array).
  - Each copy is projected onto the **reference (query) coordinates** using the long-form `cs` tag;
    reference columns missing from a copy are MISSING, columns deleted in the copy are GAP, target
    insertions are dropped. Bases are read from the target FASTA at the walked position, never from
    the `cs` base fields, so the projection cannot depend on the aligner's base ordering inside a
    mismatch (that was found to differ from minimap2's convention).
  - Trees: NJ (own numpy implementation, validated against Biopython on 20 random matrices) with
    100 alignment-**block** bootstrap replicates (50-column blocks), plus the two required confidence
    strata (high-confidence blocks >= 80 % unit representation; high-confidence units >= 90 %
    identity) and a label-permutation null.
  - Two label levels are always reported: chromosome arm and haplotype. This is essential — a tree
    dominated by haplotype structure makes arm labels polyphyletic without inter-arm exchange, so
    the discriminating pair of statistics is (arm monophyletic?) and (haplotype monophyletic?).
  - Control (c): the 100 kb flank windows. Cross-arm flank sequence is non-homologous, so a
    cross-chromosome flank tree cannot be built; the meaningful true-answer case there is ALLELIC
    monophyly (each individual's two haplotypes must be sisters), which is what is measured, plus a
    count of cross-arm alignments.
- [step 6] `--self-test` green (6 groups, 25 assertions), including
  `projection exact: gap columns == [2000..2009], differing columns == [100, 500, 2500]`,
  `20/20 random Euclidean matrices: NJ clade set == Biopython NJ clade set`,
  `known split recovered in 1.00 of 50 replicates`, and the adversarial shuffled-label case
  (`not monophyletic; extra_clades=4`).
- [step 7] First end-to-end run on D4Z4 (regime -f512, 384 copies from 40 windows, 6 min wall)
  exposed three things that were then fixed, all recorded because they change interpretation:
  1. **Distance model.** Scoring a gap column as one difference per deleted base let a single
     10-100 bp indel outvote tens of substitutions in a 3.3 kb tandem unit. Switched the primary
     model to complete deletion (a column is comparable only when BOTH copies carry a sequenced
     base). Validation gate: reference-frame p-distance vs the direct pairwise FastGA distance
     computed with the same model → median |dev| 2.5e-6, max 1.1e-3, Spearman 0.982, Pearson 0.997
     on 100 random copy pairs.
  2. **Degenerate stratum.** A fixed 0.80 block-completeness cut kept all 66 blocks (min measured
     completeness 0.887), so the "high-confidence blocks" stratum was a copy of the unfiltered one.
     Replaced by threshold = max(0.80, median block completeness); realised rule and threshold are
     written into the JSON, and the D4Z4 stratum then keeps 33 of 66 blocks.
  3. **Strict clade test is binary and brittle where distances are exactly 0** (many 10q arrays are
     identical to each other). Added a nearest-neighbour statistic (frac of copies/arrays whose
     nearest neighbour is on the same arm / same haplotype), which is interpretable where the
     clade test is not, and an array-level (majority-rule consensus per window) tree so that
     within-array unit turnover does not swamp the arm signal.
- [step 8] Two more defects found and fixed before the definitive run, both with new oracle tests:
  4. **Orientation bug (the big one).** Copies of these duplicons are often INVERTED: SEPTIN14
     copies are on the `-` strand on chr1p and chr7p and on `+` on chr2q/chr3q/chr5q/chr9q. The
     projection placed the target base as-is, so an inverted copy landed on the reference frame as
     its own reverse complement and every column "differed": a chr1p vs chr3q copy gave
     **p = 1.0000 over 2583 comparable columns** with 100 % coverage. Fixed by complementing the
     base for `-` strand alignments; self-test 7 now asserts that a reverse-complement copy projects
     identically to the forward copy (`3000 columns`). Every number produced before this fix is
     void; the earlier run is kept only as a diagnostic.
  5. **Acrocentric reference unit.** Scoring tiles on the chr13p window alone picked a unit private
     to chr13 (2327 homologs, but copies recovered on chr13 only: 20 units, 1 arm label). Now the
     unit is scored on the CHM13 windows of ALL five acrocentric p arms and the winner is the tile
     whose homologs span the most distinct arms.
  Flank control changed to one representative copy per flank window, and the flank discovery now
  reports the median/max qualifying alignments per window, which is the direct measure of whether
  the flank is unique sequence (it is, for the gene families; it is not for the satellite-rich
  acrocentric p arms, which is reported rather than silently averaged in).

## 2. Method as implemented (`scripts/cladistics/duplicon_monophyly.py`)

All parameters below are written into `duplicon_monophyly.json` under `parameters`, so every number
in the TSV/JSON is traceable to them.

**Substrates (read-only).** Copies: `hprcv2.1Mb.telo_500kb_trimmed.fa.gz` — 18,827 records, one per
(haplotype, arm), 500 kb each, 48 arms, no duplicate (haplotype, arm) keys. Flanks (control):
`hprcv2.1Mb.flanking_100kb.fa.gz` — 15,666 records, 100 kb, centromere-proximal to the PHR windows.
Superpopulation labels: the repo-local MDS table (233 individuals / 465 haplotypes;
AFR 134, AMR 88, EAS 104, EUR 65, SAS 74 haplotypes).

**Reference units (CHM13, coordinates from the RefSeq Liftoff GFF3).**
| family | community | reference unit | length | arms |
|---|---|---|---|---|
| D4Z4_DUX4L | C1 | chr4:193494500-193497800 (D4Z4 tandem unit around DUX4L6) | 3301 | 4q, 10q |
| TUBB8_TUBB8B | C2 | chr10:39416-43415 (first 4 kb of TUBB8) | 4000 | 10p, 18p |
| OR4F | C11 | chr5:181905487-181909486 (first 4 kb of OR4F3) | 4000 | 5q, 6q |
| SEPTIN14 | C3/C11 | chr5:181897110-181899694 (SEPTIN14P5) | 2585 | 1p, 2q, 3q, 5q, 7p, 9q |
| ACRO_P_REPEAT | C7 | data-driven 2 kb tile of the chr13p window, scored on all five acrocentric p arms | 2000 | 13p, 14p, 15p, 21p, 22p |

**Copy discovery.** `FastGA -pafS -l500 -i0.7 -f<512|8192> -T32` (values appended to the flag:
`-f512`, not `-f 512`), reference vs the concatenated selected windows. Units: aligned length
>= 500 bp, identity >= 0.80, de-duplicated at < 20 % target overlap, capped at 10 per window
(evenly spaced along the array so the array interior is sampled, not the array ends).

**Alignment frame and distance.** Every copy is projected onto reference (query) coordinates from the
long-form `cs` tag; inverted copies are complemented so that all copies are compared in the reference
orientation. Columns absent from a copy are missing; target insertions are dropped; nothing is
imputed. Distance = p-distance over columns where BOTH copies carry a sequenced base (complete
deletion); pairs with < 200 comparable columns take a documented saturation cap of 1.0 and are
counted in `frac_pairs_capped`.

**Trees.** Own numpy neighbour-joining (validated against `Bio.Phylo` on 20 random distance matrices
— identical bipartitions), 100 alignment-block bootstrap replicates (50-column blocks resampled with
replacement), plus the two required confidence strata (high-confidence blocks: threshold
`max(0.80, median block completeness)`, realised rule in the JSON; high-confidence units: reference
alignment identity >= 0.90) and a label-permutation null (100 draws of the same label sizes).

**Monophyly statistics per label** (chromosome arm and haplotype; also at the array level, where each
window's copies are collapsed to a majority-rule consensus): `mono_point` (does an edge separate
exactly this label?), `extra_clades_point` (minimum number of all-same-label clades needed to cover
the label; 1 iff monophyletic; greedy upper bound, exact for the 1 test), `frac_largest` (size of the
largest all-same-label clade / label size), the same over the bootstrap replicates, and a
nearest-neighbour companion statistic.
- [step 9] Acrocentric reference unit, second correction and the measurement that forced it. Scoring
  tiles only over the chr13p window cannot work: the CHM13 chr13p 500 kb window is a chr13-private
  divergent 36-mer satellite. Measured directly with FastGA `-l500 -i0.8 -f512` all-vs-all of the
  five CHM13 acrocentric 500 kb windows (2.5 Mb, 3 min 32 s wall): cross-window alignments are
  chr14↔chr21 3,398/3,395, chr14↔chr22 2,690/2,665, chr15↔chr14 2,541/2,537,
  chr22↔chr15 1,981/1,975, chr22↔chr21 1,949/1,939, chr15↔chr21 1,730/1,721, and only
  5–18 for any pair involving chr13. The picker now scans tiles from ALL five acrocentric windows
  and picks the tile whose homologs span the most distinct arms; it selected
  `chr14_parm` offset 268000–270000 with homologs on 4 arms (2,532 homologs; panel median 99).
  The C7 control therefore rests on the acrocentric p-arm repeat family shared by 14p/15p/21p/22p
  (and partially 13p), NOT on rDNA, which is unreachable here. This is stated as a limitation.

---

# VERDICT (written by the parent from `duplicon_monophyly.json` after the lane was interrupted)

## Chromosome-label monophyly, both aligner regimes

`boot_monophyly_frac` = fraction of bootstrap replicates in which one chromosome label's copies form a
single clade; `extra_clades_point` = how many clades are needed to contain that label's copies
(1 = exactly monophyletic); `frac_largest_same_label_clade` = share of the label's copies inside its
largest clade.

| family (community) | arms | boot mono | extra clades | frac largest | reading |
|---|---|---|---|---|---|
| **TUBB8/TUBB8B (C2)** | 10p, 18p | **1.00 / 1.00** | **1 / 1** | **1.00 / 1.00** | **monophyletic: no cross-arm mixing** |
| **D4Z4/DUX4L (C1)** | 4q, 10q | **0.00 / 0.00** | 13 / 4 (f512), 13 / 27 (f8192) | 0.26 / 0.95 | **not monophyletic: 4q and 10q copies interdigitate** |
| **OR4F (C11)** | 5q, 6q | 0.00 / 0.00 | 3 / 2 | 0.90 / 0.65 | not monophyletic |
| **ACRO p-arms (C7)** | 13p,14p,15p,21p,22p | **0.00 on all five** | 14–34 | 0.06–0.20 | **deepest intermixing of any family** |
| SEPTIN14 (C3/C11) | 1p,2q,3q,5q,7p,9q | 1.00 (2q), 0.96 (7p), else 0.00 | 1–7 | 0.35–1.00 | partially segregated |

Identical at `-f512` and `-f8192` except for the D4Z4 chr10q extra-clade count, so the result is not an
alignment-setting artefact.

## The control that licenses the interpretation

The adjacent **unique-sequence flank** control gives `allelic_monophyly_frac = 1.0` and
`boot_monophyly_median = 1.0` for every arm of D4Z4 and TUBB8B at both regimes. When the sequence is
single-copy, the same pipeline and the same trees recover haplotypes cleanly. The machinery is therefore
not producing default-positive or default-negative answers.

## Interpretation

Two families are clean and they disagree with each other, which is what makes this informative:
TUBB8B's 10p and 18p copies are separate clades (ancient duplication, no detectable recent exchange
between those ends), whereas D4Z4/DUX4L's 4q and 10q copies interdigitate (the signature expected of
exchange or homogenisation between those ends), and the acrocentric C7 p-arms intermix most deeply of all
— the rDNA positive control behaving as predicted.

**This TOPOLOGICAL result contradicts the earlier pairwise-distance conclusion** (M1: "no homogenisation
signature in the four named families"). That comparison saturated inside tandem arrays; topology still
records which copy is nearest whom, and it separates the families. The distance-based note stands only as
a statement about whole-window identity, not about copying.

## What this does NOT establish

- Haplotype-label monophyly is 0/20 for D4Z4, OR4F and SEPTIN14, 1/20 for TUBB8B, and 8/23 / 5/24 for the
  acrocentric set. The trees therefore resolve **tandem copy units**, not individuals. "Intermixed across
  arms" means copies from the two arms fall in the same clades — an exchange/homogenisation signature —
  and it is not a statement about haplotype of origin.
- Direction of copying is not established: interdigitation says the two ends' copies have exchanged, not
  which was donor.
- Timing is not established; nothing here dates an exchange.
- Ascertainment caveats recorded earlier in this note still apply: the D4Z4 array sits at/outside the edge
  of its PHR window, and the C7 control uses p-arm blocks ~11.7 Mb from the telomere, so it is not the same
  window class as the test.

## Publishable sentences

- TUBB8B: "Copies of the TUBB8/TUBB8B duplicon on 10p and 18p each form monophyletic clades, consistent
  with an ancient duplication with no detectable recent exchange between those chromosome ends."
- D4Z4/DUX4L: "DUX4L copies from 4q and 10q do not form separate clades: they interdigitate, the expected
  signature of repeated exchange between those ends."
- Acrocentric C7: "Across the five acrocentric p-arms the repeat copies are the most deeply intermixed in
  the panel, consistent with the documented rDNA exchange system."

## Sentence we must NOT write

"We find no evidence of concerted evolution in these families." The distance-based analysis supported that
sentence and the topology contradicts it for D4Z4, OR4F, SEPTIN14 and the acrocentrics; only TUBB8B is
cleanly monophyletic.
