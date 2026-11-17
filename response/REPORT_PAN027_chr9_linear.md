# PAN027 / chr9 subtelomeric analysis — work report

Working directory: `/moosefs/farla` (world-writable, shared)

**Project:** Subtelomeric sequence sharing and candidate chr9q/chr3q recombination  
**Working directory:** `/moosefs/farla`  
**Assembly set:** WashU pedigree v1.1 (`assemblies_v11`)  
**Status:** SweepGA alignment complete; IMPG similarity and read-based validation pending

---

## 1. Objective and data identity

### 1.1 Sample identity

**Established early, and it reframes everything:**

| sample ID | ONT accessions | label used in `…/washu_v131/ont/` |
|---|---|---|
| HG06803 | SRR37826058-60, 66, 77-78 | PAN010 |
| HG06804 | SRR37826054-57, 74-76 | PAN011 |
| **HG06807** | **SRR37826068-73** | **PAN027** |
| HG06808 | SRR37826061-65, 67 | PAN028 |

So **PAN027 = HG06807**. All alignments described here are therefore *self-alignments*: the sample's own ONT reads against its own (or an earlier version of its own) assembly — not a child-vs-parent comparison. The directory is named `HG06807` while the pipeline labels the sample `PAN027`; the SRA accessions prove they are the same data.

Only **one accession (SRR37826068)** was ever basecalled for this work; the other five remain `.sra`. That was a deliberate choice — the user confirmed the other 5 were not needed.

### 1.2 Objective of the assembly comparison

Investigate sequence sharing between chromosomes 3 and 9 across eight haplotype-resolved WashU pedigree assemblies. In particular, evaluate the candidate **20 kb chr9q/chr3q exchange** described in Figure 5 of *Concerted evolution and unorthodox recombination of human subtelomeres*: a region of PAN027 paternal chr9 reportedly has a better sequence-identity match to the father PAN011's chr3 than to his chr9.

This stage generates an all-vs-all set of chromosome alignments; it does **not** by itself demonstrate recombination.

---

## 2. Environment and Slurm jobs

### 2.1 Environment

| item | value |
|---|---|
| host (login) | `octopus01` (48 cores, 27 users) |
| Slurm | `/usr/local/bin/sbatch`, partition `workers` (octopus05-11, 48 cpus, 257-515 GB) |
| minimap2 | 2.30-r1287 — `/export/local/home/lpignata/miniconda3/envs/genomics/bin/minimap2` |
| samtools | 1.22.1 (htslib 1.22.1) — same env |
| seqkit | same env |
| SRA toolkit | fasterq-dump 3.4.1 — `/moosefs/salehi/washu_pedigree/tools/sra_tools/bin` |
| bgzip | `/home/lpignata/miniconda3/bin/bgzip` |
| impg | 0.4.1 — `~/.cargo/bin/impg` |
| awk on compute nodes | **mawk 1.3.3** (login node has gawk 4.2.1 — this matters, see §6) |

### 2.2 Slurm jobs

| job | what | result | wall |
|---|---|---|---|
| 2867388 | align 1 accession → PAN027 **v1.3.1** hap1/hap2 (separate) | COMPLETED 0:0 | 6:34:51 |
| 2867806 | chr3/chr9 extraction | COMPLETED 0:0 | 13:39 |
| 2868363 | v1.1 combined ref align, PAF + BAM | COMPLETED 0:0 | 13:02:41 |

---

## 3. Reads: alignment to the PAN027 assemblies

### 3.1 SRA → FASTQ

Input: `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/ont_sra/HG06807/SRR37826068.sra` (52,349,088,077 B)

```bash
fasterq-dump -Z --concatenate-reads -e 16 \
    --temp /moosefs/farla/tmp/SRR37826068 \
    /moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/ont_sra/HG06807/SRR37826068.sra \
  2> fastq/SRR37826068.fasterq.log | pigz -p 8 > fastq/SRR37826068.fastq.gz
```

- `--concatenate-reads` is **mandatory**: fasterq-dump 3.4.1 defaults to `--split-3`, and `-Z` then aborts with *"requested mode ( FASTQ split 3 ) would produce multiple files"*, exit 3.
- `-Z` **does** stream, but only after a size-check pass (~2 min of 0 bytes, easy to mistake for a hang).

Output: `/moosefs/farla/fastq/SRR37826068.fastq.gz` — 52,609,630,881 B (49 GiB), **1,565,811 reads**, `pigz -t` clean.

### 3.2 Alignment vs PAN027 **v1.3.1** (per haplotype)

Inputs:
```
/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/assemblies/PAN027_hap1.v1.3.1.raw.fa.gz   (mat)
/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/assemblies/PAN027_hap2.v1.3.1.raw.fa.gz   (pat)
```

Index (110 s, 22.1 GB peak RSS):
```bash
minimap2 -x map-ont -t 32 -d index/PAN027_hap2.v1.3.1.raw.fa.mmi  <fasta>
```

Align, per haplotype:
```bash
minimap2 --eqx -cx map-ont -t 32 --secondary=yes -N 20 --cs=long \
    index/PAN027_hap2.v1.3.1.raw.fa.mmi fastq/SRR37826068.fastq.gz \
    > HG06807_vs_PAN027_hap2.paf
```

Outputs (`/moosefs/farla/`):

| file | size | records | peak RSS |
|---|---|---|---|
| `HG06807_vs_PAN027_hap1.paf` | 108,858,885,598 B | 2,520,643 | 20.7 GB |
| `HG06807_vs_PAN027_hap2.paf` | 107,084,383,307 B | 2,497,248 | 20.5 GB |
| `index/PAN027_hap1…mmi` | 7,347,007,040 B | — | — |
| `index/PAN027_hap2…mmi` | 7,362,269,636 B | — | — |

Verification: `samtools view mapped-counter` sums = 1,565,811 per hap; both PAFs end on read `SRR37826068.1565811` with a complete 24-field record; `wc -l` matched the job-logged counts exactly.

### 3.3 chr3 / chr9 extraction (v1.3.1)

Job 2867806, single streaming pass over both PAFs:

```bash
awk -F'\t' -v o3=… -v o9=… '
  $6 ~ /\.chr3(\.|$)/ { print >> o3; … }
  $6 ~ /\.chr9(\.|$)/ { print >> o9; … }
' HG06807_vs_PAN027_hap1.paf HG06807_vs_PAN027_hap2.paf
```

| file | records | reads | size |
|---|---|---|---|
| `PAN027_chr3.paf` | 305,779 | 113,196 | 13 GiB |
| `PAN027_chr9.paf` | 247,072 | 71,881 | 9.3 GiB |

`wrong_target=0` on both; mat/pat split verified (chr3 152,086/153,693; chr9 134,608/112,464).

A chimeric cross-chromosome subset was produced and then **deleted at the user's request** ("we want all alignments on chr3 and chr9, not reads aligning to both").

### 3.4 chr9 ROI work at v1.3.1

ROI: `PAN027.chr9.paternal:136,130,000–136,190,000` (v1.3.1 contig = 136,216,880).

```bash
awk -F'\t' '$6=="PAN027.chr9.paternal" && $8<136190000 && $9>136130000' PAN027_chr9.paf
seqkit grep -f roi_reads.txt fastq/SRR37826068.fastq.gz > roi_reads.fastq
```

- 206 alignments / **205 reads**, 67 fully inside; `roi_reads.fastq` 18 MB, IDs identical to the read list.
- Read lengths: min 449, median 37,365, max 202,852, mean 45,900; 9,409,576 bp total (~157× over 60 kb).

#### 3.4.1 Telomere investigation

Assemblies carry real telomeric repeat (v1.3.1):

```
PAN027.chr9.paternal  (TTAGGG)552  at 136,213,282-136,216,594  (3,312 bp)
PAN027.chr9.maternal  (TTAGGG)673  at 129,067,399-129,071,437  (4,038 bp)
```

**Raw signal is unavailable** (user confirmed no pod5/fast5), so Uncalled4 signal alignment — the only way to prove array length — is off the table.

Soft-clip test instead (PAF discards clipping, so re-align to SAM):
```bash
samtools faidx PAN027_hap2…fa.gz PAN027.chr9.paternal:136016880-136216880 > chr9_distal200kb.fa
minimap2 -ax map-ont -t 8 --secondary=no chr9_distal200kb.fa roi_reads.fastq | samtools sort -o roi_vs_distal.bam
```

Result: 63 reads with trailing soft clip ≥ 1 kb, max 22,366 bp. Telomeric repeat confirmed *in the clipped tails*:

```
SRR37826068.1488241  clip 16,077 bp   1,557 TTAGGG units  (9,342 bp tandem, 58.1%)
SRR37826068.658420   clip 14,368 bp   1,088 units         (6,528 bp, 45.4%)
SRR37826068.353759   clip 13,544 bp     989 units         (5,934 bp, 43.8%)
SRR37826068.38279    clip 12,548 bp     782 units         (4,692 bp, 37.4%)
```

Phase-shifted motifs (`TGGGTTTGGG`, `TTGGGGTTGGGG`) indicate ONT basecaller error in the repeat.

**Correction:** an earlier PAF-only count gave "max 86 tandem units in any read" → concluded the assembly array was unsupported. That was wrong: the repeat is in soft-clipped tails, invisible to PAF. The assembly array is real and *understated*.

**Critical caveat discovered later:** the reads reaching the contig end in that 200 kb reference were partly *forced* there — with only chr9 in the reference every read must land on chr9. See §3.5.

### 3.5 PAN027 **v1.1** combined reference

#### 3.5.1 Combined reference

```bash
pigz -dc PAN027_mat_HiFi_element_final_mat.v1.1.raw.fa.gz \
         PAN027_pat_HiFi_element_final_pat.v1.1.raw.fa.gz \
  | bgzip -@ 8 > /moosefs/farla/PAN027_v11/PAN027_v11_hap1hap2.fa.gz
samtools faidx PAN027_v11_hap1hap2.fa.gz
```

Inputs:
```
/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/assemblies_v11/PAN027_mat_HiFi_element_final_mat.v1.1.raw.fa.gz
/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/assemblies_v11/PAN027_pat_HiFi_element_final_pat.v1.1.raw.fa.gz
```

Output: 46 contigs (23 `.maternal` + 23 `.paternal`), 0 name collisions, 6,044,255,741 bp, 1.5 GB bgzf + `.fai` + `.gzi`.
v1.1 chr lengths differ from v1.3.1 (chr3 mat 201,412,273 vs 201,412,277; chr9 pat **136,209,149** vs 136,216,880) → **the v1.3.1 `.mmi` could not be reused**.

#### 3.5.2 Index + two-pass alignment (job 2868363)

```bash
minimap2 -x map-ont -t 32 -d index/PAN027_v11_hap1hap2.fa.mmi PAN027_v11_hap1hap2.fa.gz

# pass 1 - PAF
minimap2 --eqx -cx map-ont --secondary=yes -N 50 --cs=long -t 32 \
    index/PAN027_v11_hap1hap2.fa.mmi fastq/SRR37826068.fastq.gz \
    > PAN027_vs_PAN027v11_bothhap.paf

# pass 2 - BAM
minimap2 -ax map-ont --eqx --cs=long --MD --secondary=yes -N 50 -t 32 \
    index/PAN027_v11_hap1hap2.fa.mmi fastq/SRR37826068.fastq.gz \
  | samtools sort -@ 4 -o PAN027_vs_PAN027v11_bothhap.bam -
samtools index PAN027_vs_PAN027v11_bothhap.bam
```

Outputs (`/moosefs/farla/PAN027_v11/`):

| file | size | records |
|---|---|---|
| `PAN027_vs_PAN027v11_bothhap.paf` | 180 G | 4,341,697 (mat 2,177,270 / pat 2,164,427) |
| `PAN027_vs_PAN027v11_bothhap.bam` (+`.bai`) | 65 G | 4,343,884 (1,565,811 primary, 2,736,784 secondary, 41,289 supplementary) |
| `index/PAN027_v11_hap1hap2.fa.mmi` | 13,724,713,696 B | 46 targets, 22.1 GB RSS |

Peak RSS 29.5 GB; mapped 99.95 % (99.86 % of primaries); BAM `SO:coordinate`, 46 `@SQ`.
`-N 50` chosen so the *losing* haplotype's placement survives for the AS comparison.

#### 3.5.3 ROI best-match selection (v1.1 ROI = `chr9.paternal:136,130,000–136,190,000`)

```bash
samtools view BAM PAN027.chr9.paternal:136130000-136190000 | cut -f1 | sort -u   # 210 names
samtools view -h -N roi_reads.all.txt BAM > tmp/roi_all.sam                      # 3,448 records
python3 rank_roi_reads.py tmp/roi_all.sam PAN027.chr9.paternal 136130000 136190000 out 0
```

Ranking: `dAS = AS(ROI) − AS(best overall)`, `MARGIN = 0`. Secondary ≠ supplementary.

| class | n | % | meaning |
|---|---|---|---|
| keep | 15 | 7.1 | ROI alignment is primary **and** best-scoring |
| drop | 149 | 71.0 | a different placement scores higher |
| ambiguous | 1 | 0.5 | ties with a rival at identical AS |
| split | 45 | 21.4 | read carries a supplementary block |

`dAS` percentiles: p0 −43,419 · p25 −3,336 · p50 −1,308 · p75 −159 · p95 0.

Drop winners (telomeric repeat multi-mapping): chr16.mat 32, chr15.pat 21, chr3.mat 14, chr15.mat 13, chr3.pat 13, chr19.mat 12, chr11.pat 12.

MAPQ reality at the ROI: **189/210 = MAPQ 0**, 6 in 1-9, 15 ≥ 10 (of which 2 at 60).

#### 3.5.4 Paternal / maternal split

The ROI coordinates exist **only** on the paternal contig (maternal chr9 = 129,071,589 bp). The maternal ortholog was located from the data, not liftover: maternal q-arm alignments of these reads span **128,856,537–128,995,706** (24 reads, 21 of 24 at MAPQ 0).

```bash
samtools view -bh roi_all_reads.bam PAN027.chr9.paternal:136130000-136190000 > roi_paternal.bam
samtools view -bh roi_all_reads.bam PAN027.chr9.maternal:128850000-129000000  > roi_maternal.bam
```

**The split is nested, not disjoint**: paternal-only 186, maternal-only 0, both 24.

#### 3.5.5 Chimeric-read audit

`split` was only a flag. Query-space geometry shows:

```
disjoint blocks covering >=50% of read : 21
blocks overlapping in query space     : 24   (repeat artifacts)
```

and the 321 supplementary blocks are overwhelmingly one locus:

```
305 / 321  PAN027.chr22.maternal:48,453,934-48,470,802   (typically 48,470,775)
```

Supplementary AS 92-14,129 vs primary AS 40,000-235,000. **No chr3↔chr9 pair exists** — chr9-ROI reads partner with chr22 or chrX. So this ROI gives **no 3q/9q translocation evidence**.

#### 3.5.6 PAF subsets

```bash
awk -F'\t' '
FILENAME ~ /roi_reads.all.txt$/ { keep[$1]=1; next }
($1 in keep) {
    print >> "roi_all_reads.paf"
    if ($6=="PAN027.chr9.paternal" && $8<136190000 && $9>136130000) print >> "roi_paternal.paf"
    else if ($6=="PAN027.chr9.maternal" && $8<129000000 && $9>128850000) print >> "roi_maternal.paf"
}' roi_reads.all.txt PAN027_vs_PAN027v11_bothhap.paf

# primary only (tp:A:P)
awk -F'\t' '{p=0; for(i=13;i<=NF;i++) if($i=="tp:A:P") p=1; if(p) print}' roi_paternal.paf > roi_paternal.primary.paf
```

#### 3.5.7 FASTAs for impg

```bash
samtools fasta -F 0x900 roi_all_reads.bam > roi_reads.all.fa          # full-length SEQ from primaries
seqkit grep -f roi_reads.keep.txt roi_reads.all.fa > roi_reads.bestmatch.fa
samtools faidx PAN027_v11_hap1hap2.fa.gz PAN027.chr9.paternal:136130000-136190000 > chr9p_ROI_136.13-136.19.fa
samtools faidx PAN027_v11_hap1hap2.fa.gz PAN027.chr9.maternal:128850000-129000000  > chr9m_ortholog_128.85-129.00.fa
```

---

## 4. Assemblies: chr3 and chr9 all-vs-all alignment (SweepGA)

### 4.1 Input data

**Source directory:**

```text
/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/assemblies_v11
```

| Individual | Input assembly identifier | Chromosomes extracted |
|---|---|---|
| PAN010 | `PAN010_hap1` | chr3, chr9 |
| PAN010 | `PAN010_hap2` | chr3, chr9 |
| PAN011 | `PAN011_hap1` | chr3, chr9 |
| PAN011 | `PAN011_hap2` | chr3, chr9 |
| PAN027 | `PAN027_mat` | chr3, chr9 |
| PAN027 | `PAN027_pat` | chr3, chr9 |
| PAN028 | `PAN028_hap1` | chr3, chr9 |
| PAN028 | `PAN028_hap2` | chr3, chr9 |

A total of **16 full-length chromosome sequences** (eight chr3 and eight chr9) were extracted with `samtools faidx`, using contig names from each assembly's `.fai` index. Examples include `PAN010.chr3.haplotype1`, `PAN010.chr9.haplotype1`, `PAN027.chr3.paternal`, and `PAN027.chr9.paternal`.

The initial attempt used the v1.3.1 assemblies; the **final successful run used v1.1**. Results from differing assembly versions should not be combined without explicit reconciliation.

### 4.2 Alignment strategy and parameters

The 16 extracted sequences were concatenated into a single FASTA and aligned against themselves using **SweepGA 0.1.1 with the FastGA backend**. This is a many-to-many, combined self-alignment covering chr3↔chr3, chr9↔chr9, and chr3↔chr9 relationships across all haplotypes.

The alignment settings match those in the supplied paper's `fig5_align.sbatch` example:

```bash
sweepga \
  --fastga \
  --fastga-frequency 32 \
  --num-mappings many:many \
  --scaffold-jump 0 \
  --threads 32 \
  --temp-dir "$WORK/tmp" \
  --output-file "$RAW" \
  "$COMBINED"
```

**Important distinction:** The paper's Figure 5 pedigree procedure aligns an inherited child haplotype (query) against the transmitting parent's combined two-haplotype assembly (target), then applies mapping-class, window, and similarity filters. This project instead used a combined self-alignment of 16 sequences as an exploratory input to downstream analyses; it is **not yet an exact reproduction of Figure 5**.

### 4.3 Software environment and troubleshooting

| Component | Version / condition |
|---|---|
| SweepGA | 0.1.1 (updated from 0.1.0) |
| Rust/Cargo for build | 1.99.1 stable (global default remained 1.70.0) |
| CMake for build | 3.25.2 on a Slurm worker |
| IMPG installed during this work | 0.2.4 |
| FastGA helper binaries | Installed under `/home/salehi/.cargo/libexec/sweepga/` |

Resolved setup issues:

1. The initial SweepGA build encountered CMake 3.13.4, below the dependency requirement of 3.16. A newer CMake was used on a worker node.
2. FastGA initially failed with `FAtoGDB binary not found`; the built helpers were installed to SweepGA's helper directory.
3. FastGA helpers could resolve `libz.so.1` on worker nodes, but batch and interactive environments differed.
4. A globally forced system `LD_LIBRARY_PATH` caused GLIBC version errors for Guix-linked `samtools`. The working script left chromosome extraction in its compatible environment and applied the library override only when invoking SweepGA.
5. The initial v1.1 extraction failed to detect chromosome names. Matching contig headers of the form `PAN010.chr3.haplotype1` and `PAN010.chr9.haplotype1` resolved that problem.

The FastGA build reported AGC support disabled, but this analysis used FASTA input and did not require AGC.

The original SweepGA 0.1.0 executable was preserved as:

```text
/home/salehi/.cargo/bin/sweepga-0.1.0
```

### 4.4 Final successful run

**Slurm job:** `2868382`  
**State:** `COMPLETED`  
**Exit code:** `0:0`  
**Elapsed time:** `00:59:24`  
**Input:** 16 chromosome sequences from v1.1 assemblies

**Output directory:**

```text
/moosefs/farla/sweepga_chr3_chr9_v11/
```

| Output | Observed size | Contents |
|---|---:|---|
| `raw_paf/chr3_chr9_8haplotypes.frequency32.paf.gz` | 5.6 GB | Full SweepGA-reported all-vs-all alignment PAF |
| `raw_paf/chr3_vs_chr9_8haplotypes.frequency32.paf.gz` | 33 MB | Subset with chr3↔chr9 alignments |
| `input_contigs.tsv` | — | Mapping from sample identifier to extracted chr3 and chr9 names |
| `chr3_chr9_8haplotypes.fa` | — | Combined chromosome FASTA used for the run |
| `chr3_chr9_8haplotypes.fa.fai` | — | FASTA index |

The full compressed PAF contains **1,103,883 alignment records**, confirmed with:

```bash
zcat /moosefs/farla/sweepga_chr3_chr9_v11/raw_paf/chr3_chr9_8haplotypes.frequency32.paf.gz | wc -l
# 1103883
```

These records are not equivalent to independent sequence-sharing regions or recombination events.

### 4.5 Preliminary PAN027 paternal chr9 comparisons

In the full PAF, records involving `PAN027.chr9.paternal` and either haplotype of PAN011 were counted in both PAF orientations:

| Comparison | PAF alignment records |
|---|---:|
| PAN027 paternal chr9 ↔ PAN011 chr3 (hap1 or hap2) | **8,225** |
| PAN027 paternal chr9 ↔ PAN011 chr9 (hap1 or hap2) | **28,814** |

This verifies that alignments were detected to both homologous chr9 and non-homologous chr3. **Counts alone do not determine relative sequence identity**, and the results do not yet demonstrate the paper's reported 20 kb identity crossover.

### 4.6 Interpretation and limitations

- The all-vs-all chr3/chr9 alignment completed successfully and provides the PAF data needed for downstream analysis.
- Non-homologous chr3↔chr9 matches are present, consistent with sequence sharing, but not proof of a recent inter-chromosomal exchange.
- The chr9q/chr3q candidate must be evaluated using matched coordinates and sequence identity across consecutive windows.
- Reproducing the paper's pedigree analysis additionally requires its **10:10 mapping-class filter**, exclusion of centromeric/unsupported or excessively high-depth windows, **2 kb windows**, best homologous versus non-homologous comparisons, and conservative joining of directly adjacent windows into runs.
- The installed IMPG 0.2.4 provides `partition`, `query`, and `stats`. A version supporting the paper's similarity workflow should be verified before analysis; the paper used IMPG 0.4.1.
- Read-level or orthogonal validation remains necessary to distinguish possible biological exchange from assembly or alignment artifacts.

### 4.7 Conclusion

The **SweepGA 0.1.1 all-vs-all alignment of chr3 and chr9 from eight WashU pedigree v1.1 haplotypes is complete** (Slurm job `2868382`, exit code `0:0`, runtime `00:59:24`). It produced **1,103,883 PAF alignment records**. Of these, **8,225** linked PAN027 paternal chr9 to PAN011 chr3 and **28,814** linked PAN027 paternal chr9 to PAN011 chr9.

The alignment dataset is ready for downstream IMPG-based similarity analysis. The **20 kb chr9q/chr3q candidate recombination event is not yet confirmed** by this work.

---

## 5. Similarity of assemblies and reads to the parental haplotypes (impg)

Pipeline for computing per-window similarity of the PAN027 haplotypes and of ONT reads to the chr9 and chr3 haplotypes of the parents (PAN010, PAN011), using the **v11 assemblies** throughout.

### 5.0 Inputs

| What | Path |
|---|---|
| v11 assemblies (bgzipped, `.fai` + `.gzi`) | `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/assemblies_v11/` |
| Assembly all-vs-all PAF (sweepga, chr3 + chr9, 8 haplotypes incl. PAN028) | `/moosefs/farla/sweepga_chr3_chr9_v11/raw_paf/chr3_chr9_8haplotypes.frequency32.paf` |
| Read alignments (22 primary minimap2 alignments to v11 `PAN027.chr9.paternal`) | `/moosefs/farla/PAN027_v11/chr9_roi/roi_paternal.primary.paf` |
| Reads FASTA (22 reads, SRR37826068) | `/moosefs/farla/PAN027_v11/chr9_roi/roi_reads.primary.fa` |
| Region | `PAN027.chr9.paternal:136129999-136190000` |

```bash
DIR=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/assemblies_v11
RD=/moosefs/farla/chr9_all_vs_all
C=/moosefs/farla/sweepga_chr3_chr9_v11/raw_paf/chr3_chr9_8haplotypes.frequency32.paf
R=/moosefs/farla/PAN027_v11/chr9_roi/roi_paternal.primary.paf
READS=/moosefs/farla/PAN027_v11/chr9_roi/roi_reads.primary.fa
```

### 5.1 Windows

5 kb windows over the region, dropping the 1 bp leftover window at the end (the region is 60,001 bp).

```bash
printf "PAN027.chr9.paternal\t136129999\t136190000\n" \
  | bedtools makewindows -b - -w 5000 \
  | awk '$3-$2 > 1' \
  > $RD/chr9_recreg.5kb.clean.bed

wc -l $RD/chr9_recreg.5kb.clean.bed   # expect 12
```

### 5.2 Combined PAF (assemblies without PAN028 + reads)

Original v11 sequence names are kept. Do **not** rename sequences in the PAF: names must match the FASTAs the coordinates were computed on.

```bash
NEW=$RD/ASS_pat_primary_reads.v11.paf

grep -v 'PAN028' $C > $NEW
cat $R >> $NEW

wc -l $NEW   # expect 594059 (594037 assembly + 22 read alignments)
```

### 5.3 Reads FASTA in the same format as the assemblies

impg refuses to mix plain and compressed FASTAs, so the reads are bgzipped to match the `.fa.gz` assemblies.

```bash
bgzip -c $READS > $RD/roi_reads.primary.fa.gz
samtools faidx $RD/roi_reads.primary.fa.gz

wc -l $RD/roi_reads.primary.fa.gz.fai   # expect 22
```

### 5.4 Check that names and lengths in the PAF match the FASTAs

Every sequence in the PAF must exist in the FASTAs with the same length. Empty output means everything matches.

```bash
# Names and lengths in the PAF
awk '{print $1"\t"$2; print $6"\t"$7}' $NEW | sort -u > /tmp/paf_len.txt

# Names and lengths in the FASTAs passed to impg
cat $DIR/PAN010_*.fa.gz.fai $DIR/PAN011_*.fa.gz.fai $DIR/PAN027_*.fa.gz.fai \
    $RD/roi_reads.primary.fa.gz.fai \
  | cut -f1,2 | sort -u > /tmp/fa_len.txt

comm -23 /tmp/paf_len.txt /tmp/fa_len.txt   # expect no output
```

### 5.5 Similarity

```bash
impg similarity \
  -a $NEW \
  -b $RD/chr9_recreg.5kb.clean.bed \
  --sequence-files $DIR/PAN010_*.fa.gz $DIR/PAN011_*.fa.gz $DIR/PAN027_*.fa.gz \
                   $RD/roi_reads.primary.fa.gz \
  --gfa-engine poa \
  --no-merge \
  > /moosefs/farla/ASS_pat_primary_readssimilarity.v11.tsv
```

### 5.6 Extract similarity to the 8 parental haplotypes

References are given without coordinates; the script strips `:start-end` from group names before matching.

```bash
SIM=/moosefs/farla/ASS_pat_primary_readssimilarity.v11.tsv

python extract_similarity.py \
  -i $SIM \
  -r PAN010.chr9.haplotype1 PAN010.chr9.haplotype2 PAN011.chr9.haplotype1 PAN011.chr9.haplotype2 \
     PAN010.chr3.haplotype1 PAN010.chr3.haplotype2 PAN011.chr3.haplotype1 PAN011.chr3.haplotype2 \
  -o /moosefs/farla/ASS_pat_primary_reads_vs_parents.v11.tsv
```

Keep PAN027 and the reads:

```bash
awk 'NR==1 || $1 ~ /^PAN027\./ || $1 ~ /^SRR/' \
  /moosefs/farla/ASS_pat_primary_reads_vs_parents.v11.tsv \
  > /moosefs/farla/ASS_pat_primary_reads_vs_parents.v11.filtered.tsv
```

### 5.7 Plot

```bash
python viz_optimized.py \
  /moosefs/farla/ASS_pat_primary_reads_vs_parents.v11.filtered.tsv \
  /moosefs/farla/ASS_pat_primary_reads_vs_parents.v11 \
  --threshold 0.001 --log-transform \
  --max-seqs-per-page 10
```

Main output: `ASS_pat_primary_reads_vs_parents.v11_merged_position_similarity_log.pdf`. Run once more without `--log-transform` for absolute values (`_linear` suffix): with the log transform, colours are normalised within each window.

---

## 6. Bugs found and fixed (reproducibility notes)

| # | bug | symptom | fix |
|---|---|---|---|
| 1 | fasterq-dump default `--split-3` vs `-Z` | exit 3, "would produce multiple files" | add `--concatenate-reads` |
| 2 | `length(array)` is gawk-only | job 2867801 died instantly on mawk nodes | explicit counters |
| 3 | `/tmp` is per-node | cross-node test found no inputs | use `/moosefs/farla` |
| 4 | classifier treated **secondary** as **supplementary** | "89.5 % split" — meaningless | split = any `0x800` |
| 5 | presence test `(q SUBSEP "c") in best` when only composite keys stored | 0 rows | separate `has[]` array |
| 6 | summary read `$18` (cover_frac) instead of `$19` (class) | class table showed numbers | index fix |
| 7 | `chmod 666` applied to a **directory** | `index/` and later `tmp/` lost `+x`, everything "missing" | `find … -type f -exec chmod` |
| 8 | `samtools faidx -l` invented | usage error | read `.fai` field 1 |
| 9 | `grep -c '\.maternal$'` on whole `.fai` lines | reported `mat=0 pat=0` on a good file | `cut -f1` first |
| 10 | hand-rolled CIGAR span parser required CIGAR to start with `[0-9]+[MDN=X]` | soft-clipped-leading alignments dropped; 142 instead of 210 | `samtools view bam region` |
| 11 | edited a **running** bash script | bogus "syntax error at `esac`" (byte-offset) | never edit while running |
| 12 | `pkill -f roi_bestmatch.sh` matched own shell | killed self, exit 143 | kill by PID |
| 13 | awk `-v` args arrive as strings | `TypeError: int + str` | cast |

---

## 7. Deliverable inventory

### `/moosefs/farla/` (roots)
```
fastq/SRR37826068.fastq.gz                       49 GiB   1,565,811 reads
HG06807_vs_PAN027_hap1.paf                      102 GiB   2,520,643 records   (v1.3.1)
HG06807_vs_PAN027_hap2.paf                      100 GiB   2,497,248 records   (v1.3.1)
PAN027_chr3.paf                                  13 GiB   305,779 records
PAN027_chr9.paf                                 9.3 GiB   247,072 records
index/PAN027_hap{1,2}.v1.3.1.raw.fa.mmi          6.9 GB each
run_alignment.sh, map_HG06807.sbatch, align_HG06807.sbatch, extract_chr3_chr9.{sh,sbatch}
```

### `/moosefs/farla/PAN027_v11/` (current)
```
PAN027_v11_hap1hap2.fa.gz  +.fai +.gzi           1.5 GB   46 contigs, 6.04 Gb
index/PAN027_v11_hap1hap2.fa.mmi                12.8 GiB
PAN027_vs_PAN027v11_bothhap.paf                  180 G    4,341,697 records
PAN027_vs_PAN027v11_bothhap.bam (+.bai)           65 G    4,343,884 records
build_combined_ref.sh, run_alignment_pan027_v11.sh, align_PAN027_v11.sbatch
rank_roi_reads.py, roi_bestmatch.sh, describe_chimeric.py
```

### `/moosefs/farla/PAN027_v11/chr9_roi/`
```
BAM   roi_all_reads.bam (3,448 rec / 210 reads)   roi_region_only.bam (210)
      roi_paternal.bam (210)  roi_maternal.bam (24)
      roi_bestmatch.bam (162 / 15 reads)          roi_bestmatch_region.bam (15)
      roi_chimeric.bam (1,624 / 45 reads)         (+ .bai for each)
PAF   roi_all_reads.paf (3,448)   roi_paternal.paf (210)   roi_maternal.paf (24)
      roi_paternal.primary.paf (22)             roi_maternal.primary.paf (1)
      roi_bestmatch.paf (162)
FASTA roi_reads.all.fa (210 / 9,254,451 bp)     roi_reads.primary.fa (22 / 877,799 bp)
      roi_reads.bestmatch.fa (15 / 655,551 bp)
      chr9p_ROI_136.13-136.19.fa (60,001 bp)    chr9m_ortholog_128.85-129.00.fa (150,001 bp)
TSV   roi_read_decisions.tsv   roi_ranking.summary.txt
TXT   roi_reads.{all,keep,drop,ambiguous,split}.txt
```

### `/moosefs/farla/igv/`
```
chr9_paternal.fa (+.fai)      133 MB   single contig, exact v1.3.1 coords
HG06807_roi205_on_chr9paternal.bam    205 reads, chr9-only reference (forced placement)
HG06807_roi205_vs_PAN027_hap{1,2}.bam full-genome v1.3.1 BAMs
PAN027_hap{1,2}…fa.gz symlinks (+.fai/.gzi)
chr9_telomere_features.bed
```

### `/moosefs/farla/telomere_proof/`
```
chr9_distal200kb.fa           200,001 bp extract, contig renamed chr9_distal
roi_vs_distal.bam (+.bai)     soft-clip analysis substrate
clip_analysis.tsv             361 alignments, sorted by distal clip length
run_clip_test.sh, analyse_clips.py, minimap2.log
```

---

## 8. Scientific conclusions

1. **PAN027 = HG06807**; all work here is self-alignment. Any "child vs parent" framing is wrong.
2. **chr9 q-arm ROI (136.13-136.19 Mb) is a telomeric repeat sink.** 189/210 ROI alignments are MAPQ 0; 71 % of ROI reads map *better* elsewhere (chr16.mat 32, chr15 34, chr3 27, chr19 22, chr11 24).
3. **Only 15/210 reads have their ROI placement as the best alignment**, and only 2 reach MAPQ 60.
4. **Telomeric repeat is real in the assembly and in the reads**, but the evidence lives in soft-clipped tails (up to 1,557 tandem units / 9.3 kb), not in aligned sequence.
5. **No raw signal exists** (no pod5/fast5), so array length cannot be proven; Uncalled4 is not an option.
6. **No 3q/9q recombination evidence at this ROI.** "Chimeric" reads partner chr9 with chr22 (305/321 supplementary blocks, one 17 kb locus at chr22.mat:48.45-48.47 Mb) or chrX — repeat homology, not translocation.
7. **The paternal/maternal split is nested**, not disjoint: 24 maternal-window reads are a subset of the 210 paternal-ROI reads; 21/24 of those maternal placements are MAPQ 0.

## 9. Open items

- Whether the 22-read `primary` PAF or the 15-read `bestmatch` set is the intended impg input; 5 of 15 bestmatch reads are MAPQ < 10.
- Whether the maternal ortholog window should be data-derived (128.85-129.00 Mb, used here) or the terminal 60 kb of maternal chr9 (129,011,589-129,071,589).
- v1.3.1 vs v1.1 coordinate shift: the telomeric array position and ROI must be re-derived per assembly version — v1.1 chr9.pat is 7,731 bp shorter.
- Only 1 of 6 accessions used; the other five remain `.sra` by decision.
