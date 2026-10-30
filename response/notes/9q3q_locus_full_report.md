# The 9q/3q locus in the WashU pedigree — complete record

Date: 2026-10-23. Naming: **SELF = PAN027** (focal), **DAD = PAN011**, MOM = PAN010, KID = PAN028
(`PEDIGREE_NAMING.md`). This report supersedes and consolidates
`M2b_9q3q_read_support.md`, `M1c_duplicon_monophyly.md` (for this locus) and the 9q/3q parts of
`03_change_log.md`. Provenance: `response/PROVENANCE.md`; tool pins:
`response/notes/tool_versions_and_pins.md`.

## 1. What the pipeline actually does (and which quantity is which)

1. **Assemblies** — WashU v1.1/v1.31 per haplotype: `PAN010 mat/pat`, `PAN011 hap1/hap2`,
   `PAN027 mat/pat`, `PAN028 hap1/hap2`.
2. **Windows** — telomere-anchored 500 kb flanks per chromosome arm; the q arms are reverse-complemented
   so every window runs telomere→centromere ("telomere-first").
3. **Multi-alignment** — **sweepga 0.1.1** (FastGA `-f32`, many-to-many, no scaffold jumps) aligns the
   child's flanks against the combined two-haplotype parental target genome-wide, then filtered to
   **10:10 mapping classes**. This is the multi-alignment the ribbon panel rests on; it is *not* a
   pairwise alignment of two chosen windows.
4. **Similarity** — `impg similarity` walks the alignment in **2 kb windows** and converts each
   sequence's path into a **multiset of graph nodes**, reporting `group.a.length`, `group.b.length`,
   `intersection` and the overlap metrics:
   - `jaccard = I / (A + B − I)`, `dice = 2I / (A + B)`, `cosine = I / √(A·B)`.
   Verified to six decimals against the vendored table (e.g. A=458, B=2000, I=420 → jaccard
   0.206084, dice 0.341741). Because paths loop, these are **multiset** overlaps, not set overlaps.
5. **`estimated.identity`** — in IMPG's source (`commands/similarity.rs`) this has two legitimate
   modes: Mash mode `(2J/(1+J))^(1/k)` (per-base identity from k-mer Jaccard) and POA mode
   `matches/(matches+mismatches+gap_events)` (gap-compressed). In the vendored table the column equals
   **dice exactly**, i.e. the Mash branch with k = 1 — so it is the overlap score, not a per-base or
   gap-compressed identity despite its name.
6. **Classification** — a downstream script keeps, per window, the best `same_chrom` candidate and the
   best `interchrom` candidate and labels the window by the larger score. **The class-winner table is
   not winner-only**: it carries both candidates for **6,821 of 6,958** windows, each with its own dice.
   (An earlier claim in these notes that only winners were kept was wrong.)
7. **Figure** — `fig5_ribbon_plotter.py` merges adjacent same-class windows on the same donor into runs
   and draws them. Panel B (`make_fig5_zoom.R`) plots `estimated.identity` by class across a locus.

**Two other identity measures exist and are *not* what the ribbon uses:**
- **`dv:f` in the FastGA PAF = gap-compressed divergence.** From `ALNtoPAF.c`: `diffs` counts
  "substitutions + indels", an indel of any length counting **once**. So `1 − dv` *understates*
  divergence in indel-rich sequence — which this locus is.
- **PAF cols 10/11 (`iid`/`blocksum`)** = length-penalised: indel bases count in the denominator.

## 2. The observations to explain

For the daughter's **paternal chr9q** terminus (SELF.p.chr9q), against the father's four haplotypes:

| father's haplotype | coverage of the terminal 60 kb | identity (sensitive, gap-compressed) |
|---|---|---|
| chr3 hap2 | **100%** | 99.79% |
| chr3 hap1 | 98.7% | 99.19% |
| chr9 hap2 | 97.7% | 99.45% |
| chr9 hap1 | **7.5%** | 98.68% (fragments only) |

So within the father's genome the terminal module exists on **chr3q (both haplotypes) and on chr9 hap2**,
and is **absent from chr9 hap1** — a different architecture between his own two chr9s.

## 3. The alternatives, and what each predicts

| id | explanation | prediction | verdict |
|---|---|---|---|
| **A** | Real, standing chr3q-derived duplication on DAD chr9:h2, transmitted | one haplotype only; unrelated individuals carry it; reads physically span the junction | **supported** |
| **B** | Same but recent (arose in DAD's parents' germlines) | identical data to A; separable only by population sampling | not separable with 4 individuals |
| **C** | Assembly over-duplication: the assembler pasted a chr3 module onto the chr9 terminus | ONT reads should **not** span from unique chr9 sequence into the module | **refuted** (§5) |
| **D** | Paralogous repeat family only, no duplication event | margins stay 99.5–99.8%, no boundary-locked switch; other haplotypes share the family | **supported for 5q/1p, not for 9q/3q** (§6) |
| **E** | Daughter's copy diverged from DAD's carrier in transmission | would break the 99.96–99.97% whole-window colinearity | **excluded** |
| **F** | De novo exchange in DAD's meiosis (new in SELF) | DAD chr9:h2 would not carry the module | **excluded** |
| **G** | Artefact on the chr3q donor assembly | HiFi `Dup` would be on chr3:h2 | **excluded** (chr3 termini HiFi-clean) |

## 4. Read support at every involved end (ONT)

Reads: PAN011 (DAD) ONT runs `SRR37826054`–`57` mapped to PAN011 v1.31; PAN027 (SELF) runs
`SRR37826068`–`69` mapped to PAN027 v1.31 — `/moosefs/.../pedigrees/washu_v131/ont/`.
Counts from `samtools view -c` (summed over runs); source data `response/notes/read_support.json`,
regenerate with `scripts/pedigree/ont_read_support_9q3q.py`.

| individual | contig | role | terminal 100 kb | per kb | last 5 kb | module (54.1 kb) | per kb |
|---|---|---|---|---|---|---|---|
| DAD PAN011 | chr9 hap2 | **module carrier** | 112 | 1.12 | 30 | 32 | **0.59** |
| DAD PAN011 | chr9 hap1 | other haplotype (no module) | 169 | 1.69 | 48 | — | — |
| DAD PAN011 | chr3 hap2 | donor haplotype | 223 | 2.23 | 92 | — | — |
| DAD PAN011 | chr3 hap1 | other chr3 haplotype | 171 | 1.71 | 66 | — | — |
| SELF PAN027 | chr9 paternal | **transmitted copy** | 56 | 0.56 | 14 | 35 | 0.65 |
| SELF PAN027 | chr9 maternal | contrast (divergent terminus) | 57 | 0.57 | 15 | — | — |

Reads reach the telomere at every end (14–92 within the last 5 kb), so no locus is an unanchored end.
In DAD the module carrier's terminus is the **least** covered of his four ends, and within it the module
runs at 0.59/kb (about half the terminal average) — consistent with the HiFi `Dup` (more assembly copies
than reads support). In SELF the transmitted copy shows **no such deficit** (0.65 in-module vs 0.56 over
the terminus): the transmitted arrangement is single-copy at read resolution.

## 5. Reads cross the junction — this refutes explanation C

Reads whose single alignment starts in unique chr9 sequence and extends into the module.

**DAD**, junction at chr9:h2 133,684,561:

| read | span | MAPQ |
|---|---|---|
| SRR37826054.2500049 | 133,534,793–133,707,472 (173 kb) | **60** |
| SRR37826054.1488634 | 133,576,860–133,713,601 (137 kb) | 10 |
| SRR37826054.820457 | 133,570,987–133,702,960 (132 kb) | 27 |
| SRR37826054.1901184 | 133,609,557–133,704,620 (95 kb) | 11 |
| SRR37826054.1135801 | 133,644,130–133,698,684 (55 kb) | 2 |

**SELF**, junction at chr9 136,152,460:

| read | span | MAPQ |
|---|---|---|
| SRR37826069.1449775 | 136,043,348–136,177,091 (134 kb) | 14 |
| SRR37826068.825047 | 136,074,012–136,170,981 (97 kb) | 13 |
| SRR37826069.427079 | 136,100,180–136,160,403 (60 kb) | 12 |
| SRR37826068.458129 | 136,122,720–136,182,332 (60 kb) | 8 |
| SRR37826068.1294538 | 136,130,831–136,181,822 (51 kb) | 13 |
| SRR37826068.70028 | 136,137,669–136,164,985 (27 kb) | 6 |
| SRR37826068.1084490 | 136,141,043–136,161,554 (21 kb) | 11 |
| SRR37826068.961864 | 136,150,370–136,175,638 (25 kb) | 5 |

A single molecule cannot cross a junction that exists only in an assembly. These reads establish that
the module is **physically contiguous with the chr9 terminus** in both individuals.

## 6. Multi-alignment evidence (sweepga + IMPG dice) — what Fig 5B shows

Daughter's paternal chr9q, per 2 kb window, both candidates:

| region | chr9:h2 (homolog) | chr3:h2 (paralog) | Δ(homolog − paralog) |
|---|---|---|---|
| outside the called tract (29 windows) | **1.0000** | 0.9995 | **+0.0005** |
| inside the called tract (12 windows) | **0.9933** (0.9782–1.0000) | **0.9994** | **−0.0061** |

The sign flips at the tract boundary; the paralog leads in **11 of 12** tract windows; at 136,184,000 the
homolog falls to 0.9782 while chr3:h2 holds 0.9995; chr9:h1 is absent there entirely. Per-window values
(telomere-first coordinates, dice of chr9:h2 vs chr3:h2): 164k 0.9975/0.9995, 166k 1.0000/1.0000,
168k 0.9970/1.0000, 170k 0.9925/0.9985, 172k 0.9930/1.0000, 174k 0.9922/0.9950, 176k 0.9982/0.9985,
178k 0.9953/0.9995, 180k 0.9935/0.9992, 182k 0.9897/0.9987, 184k 0.9782/0.9995, 186k 0.9897/0.9987.

**Same test at 5q/1p (daughter's paternal chr5q):**

| region | chr5 (homolog) | chr1p (paralog) | Δ |
|---|---|---|---|
| outside the called tract (57 windows) | 0.9960 | 0.9970 | −0.0010 |
| inside the called tract (14 windows) | 0.9967 | 0.9990 | −0.0023 |

No sign change; the paralog leads on both sides, and the chr1p dips lie **outside** the called tract.
5q/1p is a uniformly better paralogous family, not a tract.

## 7. Flagger status, with correct semantics

`Dup` = the assembly kept **more** copies than the reads support; `Col` = **fewer** (collapse).

| sample | HiFi flagger at the relevant terminus | ONT flagger |
|---|---|---|
| DAD PAN011 | **`Dup` 133,680,000–133,744,561** on chr9 **hap2** — spans the whole 54 kb module | no call in that 100 kb; `Col` over the terminal 32 kb of every chr9 |
| DAD PAN011 | chr3 hap1/hap2 termini: **clean** | `Col` terminal 32 kb (generic) |
| SELF PAN027 | chr9 paternal: **clean** (no call in terminal 200 kb) | `Col` terminal 32 kb (generic) |
| MOM PAN010 | chr9 hap2 terminus (same chr3-derived architecture): **no Dup** | `Col` terminal 32 kb |
| KID PAN028 | all chr9 termini: **clean** | `Col` terminal 32 kb |

Consequences: the `Dup` is real-looking but **single-technology and single-assembly** — a sample with the
same architecture (MOM) lacks it, and ONT cannot evaluate it because it collapses terminal arrays
everywhere. It is consistent with the reduced module depth in DAD (§4) and is best read as tandem
copy-number uncertainty on that haplotype, not as evidence against the arrangement (which the reads in §5
establish directly). "T2T" means gapless; it does not mean copy-number-exact, which is precisely what
these `Dup`/`Col` calls report.

## 8. Conclusions

1. **The chr9q terminus of DAD's hap2 carries a real chr3q-derived module** (~54 kb), physically
   contiguous with the chr9 unique flank (junction-spanning ONT reads), present on one haplotype only,
   and inherited by SELF at 99.96–99.97% across the whole window.
2. **chr3q is the apparent donor** in the multi-alignment framework: the homolog chr9:h2 switches from
   leading outside the tract (+0.0005) to trailing inside it (−0.0061), the paralog leads 11 of 12 tract
   windows, and the alternative chr9:h1 is absent there. This is the crossover signature Fig 5B displays.
3. **chr1p is not a donor to chr5q**: no boundary-locked switch, constant paralogous advantage — standing
   polymorphism among the grandparents' haplotypes, as the caption of Fig 5 states.
4. **Transmission direction is settled by the pedigree, and similarity plays no part in it.** DAD existed
   before SELF and SELF's chromosomes were built from his gametes, so a segment present in both must have
   passed DAD -> SELF. The comparison against his haplotypes then identifies *which* paternal copy was
   transmitted (chr9:h2 at 99.96-99.97% across the whole window versus 7.5% fragmentary coverage on
   chr9:h1); it does not establish the direction, it only names the donor chromosome.
   What is *not* oriented by either similarity or this pedigree is the **older duplication between chr3q and
   chr9q inside the paternal lineage**: a module at chr9 that came from chr3 and a module at chr3 that came
   from chr9 are indistinguishable under symmetric overlap metrics, and the grandparents are unsampled, so
   the chr9-side event cannot be traced above DAD. The available asymmetry is that DAD is heterozygous at
   chr9 (h2 carries the module, h1 is empty) while **both** his chr3 haplotypes carry it, which is the
   pattern expected of an ancestral chr3 source and a derived chr9 recipient; and MOM's chr9:h2 carries the
   same chr3-derived terminus, so the arrangement either predates both parents or arose more than once.

## 9. Corrections log (what was said wrong, and why)

| claim | status | why it was wrong |
|---|---|---|
| "the class-winner table keeps only the winner" | **retracted** | both candidates are present for 6,821/6,958 windows |
| "Fig 5B is class-selective, so its drop is an artefact of missing points" | **retracted** | both series exist for every window in the tract; the drop is real dice |
| "gap-compressed identity favours the homolog, so no dip exists" | **retracted** | that came from *pairwise* FastGA over 500 kb windows, a different framework; in the multi-alignment dice the dip is real |
| "`estimated.identity` is dice because k = 1 (degenerate mode)" | **retained** | verified against `similarity.rs` and reproduced D = 2I/(A+B) to 6 dp |
| "the 1p/5q and 9q/3q loci behave identically under both metrics" | **retracted** | they agree on per-window medians but differ in the decisive *switch* at the tract boundary |
| "flagger independently confirms the module" | **retracted** | the `Dup` is single-technology and absent from MOM's identical architecture |

## 10. Open items

1. **A vs B** (ancestral vs one-generation-old): needs population sampling; HPRC-wide recurrence of the
   chr3-derived chr9 terminus would settle it.
2. **Module copy number on DAD chr9:h2**: count module units traversed per spanning read, or HiFi depth
   per copy, to test tandem amplification versus multi-mapping loss.
3. **Fig 5B axis label** says "best-match sequence identity" but plots dice — must be corrected.
4. The same panel has not yet been built for 5q/1p as Extended Data.

## 11. Reproduction

- Read support: `scripts/pedigree/ont_read_support_9q3q.py` → `response/notes/read_support.json`
- Multi-alignment dice: `data/fig5_PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz`
  (columns `group.a.length`, `group.b.length`, `intersection`, `dice.similarity`,
  `estimated.identity`, `winner_class`, `other_seq`); arithmetic verified as in §1(4).
- Junction-spanning read queries: `samtools view <bam> "<contig>:<start>-<end>"`, ref-span from CIGAR.
