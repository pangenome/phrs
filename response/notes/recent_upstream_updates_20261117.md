# Upstream updates, 15-17 November: both Fig. 5 candidates look like artefacts

Checked after `dc4fbbc`. Two contributors pushed work after our last exchange.

## 1. Laura Pignata (`FarLa` branch; merged to main then deleted)

Uploaded four items to main (`eaa487d`) and then deleted all four (`0bac6eb`, `ad49d63`, `ec442b2`,
`dc4fbbc`), leaving main net-zero; the files live on `origin/FarLa`.

Recovered from history: `response/REPORT_PAN027_chr9_linear.md`, a 648-line work report on the PAN027
chr9 analysis. Its most useful content is a **sample-identity mapping** that our notes never had:

| sample | ONT accessions | label in `washu_v131/ont/` |
|---|---|---|
| HG06803 | SRR37826058-60, 66, 77-78 | PAN010 |
| HG06804 | SRR37826054-57, 74-76 | PAN011 |
| **HG06807** | SRR37826068-73 | **PAN027** |
| HG06808 | SRR37826061-65, 67 | PAN028 |

She also notes that only `SRR37826068` was basecalled for her work and the other five accessions remain
`.sra`; we used SRR37826054-57 (PAN011) and SRR37826068-69 (PAN027), which exist as BAMs, so her
statement describes her own scope rather than the limit of what is available.

## 2. Andrea Guarracino: a non-circular correction, and the calls do not survive it

Her earlier patch was a consensus framed on SELF.p, which is circular for the SELF-vs-DAD comparison
(see `andrea_patch_reconciliation.md` §4). She has now replaced that method. From
`results/corrections_dad_v11.tsv`:

> `PAN011#2#chr9  133,594,275-133,744,561  150,287 -> 147,583 bp  corrected  hifiasm contig from DAD
> HiFi only (no SELF data), polished with DAD HiFi`

**DAD HiFi only, explicitly no SELF data.** That removes my circularity objection and is a fair test.
The same table also flags other known assembly errors (chr3 98,480,001-98,580,000 ONT 1,197 SNV, HiFi
depth 0; chr10 D4Z4 134,420,001-134,592,605; chr7 1,076,489-1,120,704 corrected).

### Both candidates are unstable across the FastGA seed-frequency cutoff

From `fig5_cnd_runs_compare.txt`, comparing drawn interchromosomal runs at different `-f`:

| comparison | present at lower `-f`, ABSENT at higher `-f` |
|---|---|
| f32_v11 vs **f48_v11** | chr5 182,010,000-182,022,000 -> chr1 h2 0.9987; **chr5 182,052,000-182,080,000 -> chr1 h2 0.9983** (the chr5q/1p candidate) |
| f48_v11 vs **f96_v11** | **chr9 136,168,000-136,188,000 -> chr3 h2 0.9988** (the chr9q/3q candidate) |

So the chr5q/1p candidate exists only at `-f32` and vanishes at `-f48` (and returns at `-f96`), while the
chr9q/3q candidate survives `-f48` but vanishes at `-f96`. **Neither is stable across the parameter that
Methods has to state.**

### On the corrected assembly the run set collapses

`fig5_cnd_runs_compare.txt` for the corrected-DAD runs: `all_runs=122 drawn_runs=2` (against 51 drawn
runs on the original). Her G2 commit messages say the same: *"G2 on DAD corrected with own reads only
(PAR1 only)"*.

## 3. Synthesis: identity cannot assign a donor at the telomere

The number that made the 9q/3q call look strong - SELF.p.chr9q matching DAD chr3 hap2 at **99.79% over
100% coverage** - is not discriminating here. Our own cohort data put the median best non-homologous
partner at the tip at **0.52% divergence (3.5x allelic), with 11 of 45 ends inside 2x allelic**. A 99.79%
match is **0.21% divergence**, i.e. squarely inside the normal range for tip-family similarity between
unrelated ends. A near-identical match to a plausible donor is therefore the *default* at subtelomeric
tips and carries almost no information about direction or about whether a transfer occurred.

Combine that with (i) an assembly error at exactly this locus in the two releases that were both
assembled from the same data, (ii) run instability across `-f`, and (iii) the signal disappearing on a
correction built without SELF data, and the honest verdict is that **the Fig. 5 pedigree exchange is not
established**. It should be withdrawn or restated as an artefact candidate, and the same applies to the
chr5q/1p call.

## 4. What this means for the manuscript (not yet applied)

- Results `Ectopic subtelomeric recombination within a pedigree`: the retained single-event claim cannot
  stand as written.
- Fig. 5B/5C: the drawn run `chr9_run0054` (SELF.p 136,168,000-136,188,000 -> DAD h2 chr3) is
  parameter-dependent and assembly-dependent.
- Methods must state the `-f` sensitivity and the assembly-error scan regardless of the outcome, because
  this is a general property of terminal windows, not a peculiarity of one locus.
- What survives intact: the community/tier structure, the allelic-vs-paralogous framework, the exchange
  rate, the coalescence bound, the TERRA-promoter analysis, and ED Figs 5-6 - none of which depend on the
  pedigree call.

## 5. Outstanding

1. Confirm which 2 runs remain on the corrected assembly, and at which `-f`, from her run tables.
2. Establish whether MOM's chr9q module (her README point 6, three MAPQ 60 reads, chr9-type bases at
   75-78 tract sites) is real; if MOM carries it and SELF carries it, the arrangement is inherited rather
   than assembled de novo - but that is now a statement about a haplotype in the population, not about an
   exchange observable in this pedigree.
3. Decide the manuscript action: withdraw, or reframe the pedigree section as a worked example of how
   such calls fail.
