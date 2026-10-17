# Pedigree naming convention (WashU 3-generation T2T pedigree)

Canonical naming for all analyses, figures, and file names in this project.
**Focal individual = PAN027.** All nicknames are relative to her.

Source: Cechova et al. 2025, "Complete genomes of a multi-generational pedigree…"
(bioRxiv 10.64898/2025.12.14.693655, PMC12746033). v1.1 assemblies at
`/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/assemblies_v11/`.

## 1. The pedigree

```
   MOM (F)  PAN010 / HG06803  ─┐
                               ├──►  SELF (F)  PAN027 / HG06807  ──┬──►  KID (F)  PAN028 / HG06808
   DAD (M)  PAN011 / HG06804  ─┘                                   │
                                        PARTNER (M) —— UNSAMPLED ——┘
```

- SELF is the daughter of MOM x DAD.
- KID is the daughter of SELF x PARTNER.
- The paper's labels (grandmother/grandfather/mother/granddaughter) are relative
  to KID. We use SELF-relative names so the focal person is unambiguous.

## 2. Names

| nickname | PAN id | HG id | sex | relationship to SELF | v1.1 haplotype labels in files |
|---|---|---|---|---|---|
| **MOM** | PAN010 | HG06803 | F | mother | `PAN010_hap1`, `PAN010_hap2` (no parent-of-origin; her parents unsampled) |
| **DAD** | PAN011 | HG06804 | M | father | `PAN011_hap1`, `PAN011_hap2` (XY; his parents unsampled) |
| **SELF** | PAN027 | HG06807 | F | — (focal) | `PAN027_mat` (= from MOM), `PAN027_pat` (= from DAD) — parental k-mer assignment, both parents sequenced |
| **KID** | PAN028 | HG06808 | F | daughter | `PAN028_hap1`, `PAN028_hap2` — **plain labels; parent-of-origin is PER CHROMOSOME and not in the file names** |
| **PARTNER** | — | — | M | father of KID | **not sequenced** — no assembly exists |

## 3. Chromosome-copy notation (use this, never bare hap1/hap2)

Each individual has two chromosome copies. Write:

- `MOM.h1`, `MOM.h2` — MOM's copies (arbitrary internal numbering; fine because
  MOM is a founder in this pedigree)
- `DAD.h1`, `DAD.h2` — DAD's copies
- `SELF.m` (= `PAN027_mat`) and `SELF.p` (= `PAN027_pat`) — verified: SELF.p chr3q
  matches DAD.h1 chr3q at 99.98%; SELF.p chr9q matches DAD.h2 chr9q at 99.96%;
  SELF.m chr5q matches MOM.h1 chr5q at 99.58%.
- `KID.m` (from SELF) and `KID.p` (from PARTNER) — **per chromosome**. Resolve
  from v1.1 sequence (the copy matching SELF at transmission grade, >=99.8%)
  before every analysis; do not infer from file names and do not import the
  v1.0 per-chromosome table blindly.

### KID per-chromosome mapping (v1.1, verified by transmission-grade identity)

| chromosome | `KID.m` is | `KID.p` is |
|---|---|---|
| chr1 | hap2 | hap1 |
| chr3 | hap2 | hap1 |
| chr5 | hap2 | hap1 |
| chr9 | hap2 | hap1 |

Evidence (terminal 60 kb, per-position best identity to SELF):
chr9q hap2 99.97% (other copy 9.9%); chr5q hap2 99.81%; chr1p hap2 99.99%;
chr3q hap2 99.70%.

The remaining chromosomes are **not yet verified in v1.1**. The v1.0 analysis
(paper Fig 1B, and `maternal.PAN028.blocks.bed` / `paternal.PAN028.blocks.bed`
in github.com/biomonika/washu-pedigree) reports maternal = hap2 for
chr1,3,4,5,7,9,12,14,15,18,21,22,X and hap1 for chr2,6,8,10,11,13,16,17,19,20 —
treat as *provisional* until checked on v1.1 the same way.

## 4. Locus / window notation

Arms: `chr<NN><p|q>`. Windows are **telomere-first**, 500 kb by default, named

```
<individual>.<copy>.<arm>            e.g.  SELF.p.chr9q
```

Write the flank convention explicitly when it matters:
`chr9q:last500kb_revcomp_telifirst`, `chr1p:1-500kb_telifirst` (chr1p is already
telomere-first in forward orientation), and give native coordinates plus the
assembly when reporting a tract.

Haplotype-resolved tract notation:
`SELF.p.chr9q:13-57kb` = SELF's paternal chr9q copy, telomere-first 13-57 kb.

## 5. Figure / file naming

```
<arm>_<contrast>_<copies>[_<method>].<ext>
```

Examples:
- `chr9q_tip_SELF.p_vs_DAD.h2_chr3q_sweepga10to10.pdf`
- `chr9q_tip_KID.m_vs_SELF.p_pafplot_raw.pdf`
- `chr5q_tip_SELF.m_vs_MOM.h1_sweepga10to10.pdf`

Always state in the caption whether the PAF is **filtered** (sweepga
`--num-mappings 10:10`) or **raw/unfiltered** (FastGA `-paf` straight through),
because the filter removes multi-mapping copies.

## 6. Standing vs transmitted language

- **standing architecture**: an arrangement already present in an ancestor or in
  the population (e.g. the chr3q-derived chr9q terminus in DAD.h2 and MOM.h2).
- **transmitted**: a segment inherited within the pedigree, identity should be
  near-identical (>=99.8% on v1.1).
- **ectopic-derived**: sequence at an arm end that derives from a non-homologous
  arm (established by structure — e.g. the homologous copy is absent — not by a
  wide identity margin, since paralogous subtelomeric copies are ~99.5-99.8%
  similar).

## 7. Results by copy (status after the rev 10 correction)

| analysis | copies compared | status |
|---|---|---|
| chr9q terminus, SELF's paternal copy | `SELF.p.chr9q` vs `DAD.h1/h2.chr3q`, `DAD.h1/h2.chr9q` | **stands** — SELF.p chr9q terminus is chr3q-derived; DAD.h2.chr3q is the best explanation (99.79-99.84%, 100% covered); DAD.h1.chr9q excluded (7.5%) |
| chr9q terminus, KID's maternal copy | `KID.m.chr9q` vs `SELF.p.chr9q` | **stands (corrected)** — 99.97% identity, 100% covered: clean telomere-proximal crossover in SELF |
| chr5q terminus, KID's maternal copy | `KID.m.chr5q` vs `SELF.m.chr5q` | **stands (corrected)** — 99.81%: plain transmission; no chr1p event |
| chr1p terminus, KID's maternal copy | `KID.m.chr1p` vs `SELF.p.chr1p` | **stands (corrected)** — 99.99%: plain transmission |
| chr3q terminus, KID's maternal copy | `KID.m.chr3q` vs `SELF.p.chr3q` | **stands (corrected)** — 99.70% transmitted; chr9q only 99.44% (family similarity) |
| child-side analyses run on `KID.p` (= hap1) in rev 4-9 | wrong copy (PARTNER's) | **retracted**: chr3q G3 "ectopic patch"/"near-identical core" boxes, and the "KID chr5q tip is chr1p-like" observation |

## 8. Hard rules

1. Never write bare `hap1`/`hap2` for KID without stating the chromosome.
2. Never assume `PAN028_hap1` is maternal (it is paternal for chr1,3,5,9).
3. PARTNER is unsampled: any statement about KID's paternal copy is a statement
   about an unsequenced individual and cannot support a pedigree-event claim.
4. Filter depth changes what you see: report filtered and unfiltered separately.
5. flagger: HiFi-based flagger QC is the relevant one; ONT-based `Col` calls
   cover the terminal tens of kb of every chromosome end in every sample and
   carry no call-specific information.
