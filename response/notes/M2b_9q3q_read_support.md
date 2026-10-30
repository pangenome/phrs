# R1-M2b — What the ONT reads say at the 9q/3q locus: alternatives and read support

> Provenance for this note: see `response/PROVENANCE.md` (tool versions: `response/notes/tool_versions_and_pins.md`).

Reads: PAN011 (DAD) ONT, 4 runs `SRR37826054`-`57`, mapped to the PAN011 v1.31 assembly; PAN027
(SELF) ONT, 2 runs `SRR37826068`-`69`, mapped to PAN027 v1.31. Source dir
`/moosefs/.../pedigrees/washu_v131/ont/`. Alignments queried directly with `samtools view`; no
re-analysis. Companion data: `read_support.json` (this directory's working copy is
`/tmp/flagger_v11/read_support.json`; regenerate with the script below).

## The alternatives, and what each predicts

Setting: DAD's chr9 **hap2** carries a ~54 kb chr3q-like terminal module; SELF inherited that copy
(SELF.p.chr9q vs DAD.chr9:h2 = 99.96-99.97% across the whole 500 kb window). Observations to explain:
(1) DAD chr9:h2 terminus matches his chr3:h2 at 99.79% (100% coverage), chr3:h1 at 99.19%, while his
chr9:h1 is absent there (7.5% coverage); (2) the module is present on ONE of his chr9 haplotypes;
(3) node-overlap (dice) prefers chr3q at the terminus, gap-compressed identity prefers chr9q;
(4) HiFi flagger marks a `Dup` spanning exactly the module on DAD chr9:h2 only, while ONT flagger
marks `Col` at the terminal ~32 kb of every chr9 in every sample; PAN010's equally chr3-derived
chr9 terminus has no HiFi call.

| id | explanation | prediction | verdict |
|---|---|---|---|
| **A** | Real, standing chr3q-derived duplication on DAD chr9:h2, transmitted | one haplotype only; unrelated individuals carry it; reads physically span the junction | **supported** |
| **B** | Same but recent (arose in DAD's parents' germlines) | identical data to A; separable only by population sampling | not separable with 4 individuals |
| **C** | Assembly over-duplication: the assembler pasted a chr3 module onto the chr9 terminus | ONT reads should NOT span from unique chr9 sequence into the module; HiFi depth over the module ~0.5x | **refuted** (see below) |
| **D** | Paralogous repeat family only, no duplication event | margins stay 99.5-99.8%, no sharp donor; dice/identity disagree; other haplotypes share the family | **partly supported**; A and D differ only in genealogical depth |
| **E** | Daughter's copy diverged from DAD's carrier in transmission | would break the 99.96% whole-window colinearity | **excluded** |
| **F** | De novo exchange in DAD's meiosis (i.e. new in SELF) | DAD chr9:h2 would NOT carry the module | **excluded** |
| **G** | Artefact on the chr3q donor assembly | HiFi `Dup` should be on chr3:h2 | **excluded** (chr3 termini HiFi-clean) |

## Read support at every involved end

`samtools view -c` over the terminal 100 kb, and over the last 5 kb (reads that reach the telomere),
summed over all runs per individual.

| individual | contig | role | terminal 100 kb | per kb | last 5 kb |
|---|---|---|---|---|---|
| DAD PAN011 | chr9 hap2 | **module carrier** | 112 | 1.12 | 30 |
| DAD PAN011 | chr9 hap1 | other haplotype (no module) | 169 | 1.69 | 48 |
| DAD PAN011 | chr3 hap2 | donor haplotype | 223 | 2.23 | 92 |
| DAD PAN011 | chr3 hap1 | other chr3 haplotype | 171 | 1.71 | 66 |
| SELF PAN027 | chr9 paternal | **transmitted copy** | 56 | 0.56 | 14 |
| SELF PAN027 | chr9 maternal | contrast (divergent terminus) | 57 | 0.57 | 15 |

Module-restricted depth (54.1 kb): DAD chr9:h2 **32 reads, 0.59/kb**; SELF paternal chr9
**35 reads, 0.65/kb**.

## The decisive observation: reads cross the junction

Reads whose single alignment begins in unique chr9 sequence and extends into the module:

DAD, junction at chr9:h2 133,684,561 (module boundary), reads starting 24-150 kb proximal:

| read | span | MAPQ |
|---|---|---|
| SRR37826054.2500049 | 133,534,793-133,707,472 (173 kb) | 60 |
| SRR37826054.1488634 | 133,576,860-133,713,601 (137 kb) | 10 |
| SRR37826054.820457 | 133,570,987-133,702,960 (132 kb) | 27 |
| SRR37826054.1901184 | 133,609,557-133,704,620 (95 kb) | 11 |
| SRR37826054.1135801 | 133,644,130-133,698,684 (55 kb) | 2 |

SELF, junction at chr9 136,152,460, reads starting 39-166 kb proximal and extending into the module:

| read | span | MAPQ |
|---|---|---|
| SRR37826069.1449775 | 136,043,348-136,177,091 (134 kb) | 14 |
| SRR37826068.825047 | 136,074,012-136,170,981 (97 kb) | 13 |
| SRR37826069.427079 | 136,100,180-136,160,403 (60 kb) | 12 |
| SRR37826068.458129 | 136,122,720-136,182,332 (60 kb) | 8 |
| SRR37826068.1294538 | 136,130,831-136,181,822 (51 kb) | 13 |
| SRR37826068.70028 | 136,137,669-136,164,985 (27 kb) | 6 |
| SRR37826068.1084490 | 136,141,043-136,161,554 (21 kb) | 11 |
| SRR37826068.961864 | 136,150,370-136,175,638 (25 kb) | 5 |

**C is refuted.** A single DNA molecule cannot cross a junction that exists only in the assembly.
Reads at MAPQ 60 (DAD) and 14 (SELF) span 100-170 kb from unique sequence into the module, so the
module is physically contiguous with the chr9 terminus in both individuals. Reads also reach the
telomere itself in every case (14-92 reads within the last 5 kb), so the loci are not unanchored ends.

**What the depths add.** In DAD, the module carrier's terminus is the *least* covered of his four
ends (1.12/kb vs 1.69-2.23), and within it the module runs at 0.59/kb — about half the terminal
average. That is exactly what the HiFi `Dup` reports: the assembly carries more module copies than
the reads justify, either because the module is tandemly amplified on this haplotype (A) or because
reads multi-map between the chr9 and chr3 copies and are lost from this locus (D). In SELF the
transmitted copy shows no such deficit (0.65/kb in the module vs 0.56/kb over the terminus), i.e.
the transmitted arrangement is a single-copy locus at read resolution.

## What remains open

1. **A vs D** — tandem amplification on DAD chr9:h2 versus a shared repeat family without extra
   copies. Separable by counting module units traversed per spanning read, or by depth per copy.
2. **Copy number of the module** — the `Dup` call says the assembly has more copies than reads;
   neither the exact count nor whether it is real amplification is established here.
3. **Nothing here dates the arrangement.** Both A and B remain consistent; population sampling across
   HPRC would separate them.

## Reproduction

`scripts/pedigree/ont_read_support_9q3q.py` (this note's companion script) recomputes the table
above and the junction-spanning read lists from the two BAM sets.

---

## Multi-alignment evidence (sweepga + IMPG) — the framework that produced Fig 5B

Correction to earlier text in this note and in the analysis notes: the class-winner table is **not**
winner-only. It carries one `same_chrom` and one `interchrom` candidate for **6,821 of 6,958**
windows, each with its own dice, so it *is* the multi-alignment comparison (sweepga many:many →
IMPG node-multiset dice), not a single alignment of two windows.

**Per-window dice, daughter's paternal chr9q (− PAN027#2#chr9 vs father PAN011):**

| region | chr9:h2 (homolog) | chr3:h2 (paralog) | Δ(homolog − paralog) |
|---|---|---|---|
| outside the called tract (29 windows) | **1.0000** | 0.9995 | **+0.0005** |
| inside the called tract (12 windows) | **0.9933** (range 0.9782-1.0000) | **0.9994** | **-0.0061** |

The sign flips at the tract boundary and the paralog leads in 11 of 12 tract windows; at
136,184,000 the homolog falls to 0.9782 while chr3:h2 holds 0.9995. chr9:h1 is absent there
entirely. **This is the crossover signature, and it is what Fig 5B shows.**

**Same test at 5q/1p (daughter's paternal chr5q − PAN027#2#chr5 vs father PAN011):**

| region | chr5 (homolog) | chr1p (paralog) | Δ |
|---|---|---|---|
| outside the called tract (57 windows) | 0.9960 | 0.9970 | −0.0010 |
| inside the called tract (14 windows) | 0.9967 | 0.9990 | −0.0023 |

No sign change, and the paralog leads both inside and outside; the dips in the chr1p series lie
*outside* the called tract. So 5q/1p is a uniformly better paralogous family, not a tract.

**Conclusions now supported by the multi-alignment evidence:**

1. chr3q is the apparent donor to the daughter's paternal chr9q terminus (dice; homolog switches
   from leading outside to trailing inside; the alternative chr9:h1 is absent). The ONT reads show
   the arrangement is real (junction-spanning reads at MAPQ 60 in DAD and 14 in SELF; telomere-reaching
   reads in both), so this is a sequence-level finding, not an assembly artefact.
2. chr1p is **not** a donor to chr5q: there is no boundary-locked switch, only a constant paralogous
   advantage, consistent with standing polymorphism among the grandparents' haplotypes.
3. **Transmission direction is settled by the pedigree, not by similarity.** DAD existed first and SELF's
   chromosomes derive from his gametes, so a shared segment must have gone DAD -> SELF; the comparison only
   names which paternal copy (chr9:h2) carried it. Undetermined is the *older* within-genome duplication
   between chr3q and chr9q: symmetric metrics cannot orient it, and the grandparents are unsampled. The
   informative asymmetry is DAD's heterozygosity at chr9 (h2 carries, h1 empty) against both chr3 haplotypes
   carrying the module, plus the same chr3-derived terminus in MOM's chr9:h2.

**Fig 5B defect that remains:** the axis is dice (node overlap) but the caption calls it
"best-match sequence identity". That label must change; the metric itself is legitimate and is the
correct one here, because it penalises incomplete coverage and is sensitive to indels and SVs.
