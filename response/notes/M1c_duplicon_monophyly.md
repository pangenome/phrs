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
