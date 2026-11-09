# fig5_v131_check: the 9q/3q call in Fig. 5 and the DAD chr9q assembly

Names follow `PEDIGREE_NAMING.md`: DAD = PAN011, MOM = PAN010, SELF = PAN027 (SELF.p = `PAN027#2`, the copy from DAD), KID = PAN028.
Coordinates are 1-based. "v1.1" and "v1.3.1" are the two WashU assembly releases; the Fig. 5 calls use v1.1.
Data, BAMs and outputs are on octopus in `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131` (below: `$D`).

## Result

The 9q/3q signal in Fig. 5B comes from an error in the DAD assembly. The error is the same in v1.1 and v1.3.1.

1. The last 61,919 bp of DAD chr9 haplotype 2 (`PAN011#2#chr9`, v1.1 133,682,643-133,744,561; v1.3.1 133,683,117-133,745,035) is a copy of DAD chr1 haplotype 1 positions 8-61,926 (minus strand, 12 differences, minimap2 asm5).
2. The 20 kb Fig. 5 tract on DAD chr9 h2 (v1.1 133,702,042-133,722,041; v1.3.1 133,702,516-133,722,515) is identical letter by letter to the reverse complement of DAD chr1 h1 22,528-42,527 (md5 prefix `e075205bcb5e`, both versions).
3. DAD ONT reads (7 runs) disagree with the DAD assembly at 106 positions in the last 330 kb of chr9q h2 (v1.3.1: 97), all inside the copied block. SELF reads on the DAD assembly disagree at 185 positions (v1.3.1: 190) from 60 kb before the tract to the chromosome end. On SELF.p chr9q, DAD reads and SELF reads disagree at 0 positions. A position counts when at least 3 MAPQ >= 20 reads show a base and at least 80% of them share one base that is not the assembly base.
4. Read depth (all DAD reads) is 16 on the copied block of chr9q h2 and 73 at 100-200 kb before it. On chr1p h1 the copied zone has 60 against 85 at 162-262 kb. 16 + 60 = 76: the genome has one copy of this sequence, on chr1p.
5. v1.1 release QC: HiFi flagger calls Dup on chr9 h2 133,680,000-133,744,561 and on chr1 h1 0-48,000; nucflag calls MISJOIN on chr9 h2 133,673,490-133,686,589 (the start of the block). ONT flagger has no call on the chr9 h2 end.
6. MOM carries the same chr3q-like module on one chr9q: the SELF.p tract matches MOM chr9 h2 136,943,156-136,963,149 at 0.9988 over 20,001 bp. Three MOM reads with MAPQ 60 (270-297 kb) start 15.3-43.5 kb before the chr3q homology on SELF.p chr9q, reach the chr9 end and carry the chr9-type base at 75-78 tract sites.

Reading: DAD's real chr9q h2 end is the sequence that SELF.p carries. The drop of identity to DAD chr9q in Fig. 5B compares SELF.p with the chr1p copy in the DAD assembly. The chr9q end is an inherited haplotype that also occurs in MOM, not an exchange in this meiosis.

About `response/notes/9q3q_locus_full_report.md` section 13 (commit 8fbb742): the statement "zero chr1p alignment inside the module" does not hold, see point 2. Its read re-assignment (59 of 90 reads prefer "hap2 as assembled") has no chr1p h1 candidate, and chr1p h1 is identical to the assembled hap2 end, so a read from chr1p also prefers "hap2 as assembled".

## Patched DAD assembly, v1.1

`patch_dad_chr9.sbatch` replaces the end of `PAN011#2#chr9` with a sequence built only from DAD ONT reads. Files are in `patch_v11/`; the patched haplotype and the Fig. 5 inputs are on octopus in `$D/patch_v11/`.

Method:

1. Anchor. minimap2 asm5 of the last 400 kb of DAD chr9 h2 against the last 400 kb of SELF.p chr9. The anchor is the end of the last exact match of at least 10 kb before the nucflag MISJOIN start (133,673,490). Result: DAD 1-133,673,490 is kept; the exact match that ends there is 46,105 bp; it continues on SELF.p from 136,139,447.
2. Reads. DAD reads from the 7 DAD-on-SELF cross BAMs (`ont/SRR*.PAN011.onPAN027.v11.bam`), MAPQ >= 20, no secondary, over SELF.p chr9 136,139,447-136,209,149: 134 reads, depth min 19, median 100, max 112.
3. Consensus. `samtools consensus -m simple -c 0.6 -d 5 -H 1.0 --min-MQ 20` over that region. SELF.p is only the frame: every base comes from DAD reads. The first pass gave 69,687 bp with 93 N (positions where fewer than 60% of the reads agree).
4. Fill. `patch_read_check.sbatch` maps 2,951 DAD reads (primary records over DAD chr9q h2, chr1p h1, chr3q h1 and h2, and SELF.p chr9q) to the full DAD diploid with the patch (minimap2 map-ont, `-I 8G`). At each of the 93 N, at least 3 MAPQ >= 20 reads show a base and at least 80% of them agree; the patch takes that base. Final segment: 69,687 bp, 0 N.

Result (`patch_v11/patch.v11.tsv`, `patch_v11/patch_vs_refs.v11.tsv`):

| item | value |
|---|---|
| replaced | `PAN011#2#chr9` 133,673,491-133,744,561 (71,071 bp) |
| new chr9 h2 length | 133,743,177 bp (was 133,744,561) |
| segment md5 | `159970f65d7b7af57a67345f45b66abb` |
| patched chr9 md5 | `771cf29c513486519037008ef5d4f264` |
| segment vs SELF.p chr9q (136,139,447-136,209,149) | 4 SNV, 19 deletions (22 bp), 4 insertions (6 bp); the 4 SNV are at 136,203,675-136,205,883, in the last 5.5 kb |
| segment vs old DAD end | 286 SNV, 26 deletions, 22 insertions over 64,700 bp |
| segment vs DAD chr1p h1 | 286 SNV, 24 deletions, 21 insertions over 64,938 bp |

Read check (`patch_v11/read_check.v11.tsv`, same 2,951 reads on both assemblies, MAPQ >= 20, the 80% rule above):

| DAD reads on | wrong positions, chr9q h2 end | depth, chr9q h2 end / 100-200 kb before | wrong positions, chr1p h1 8-61,926 | depth, chr1p copy / 162-262 kb |
|---|---|---|---|---|
| original DAD v1.1 | 57 | 17 / 72 | 0 | 56 / 83 |
| patched DAD v1.1 | 0 | 77 / 72 | 0 | 67 / 83 |

The patched end is what the DAD reads support. The open test is the Fig. 5 alignment on the patched assembly: if the 9q/3q call comes from the assembly error, the PAN027pat 10:10 alignment on the patched DAD must give SELF.p chr9q its best match on DAD chr9q.

In Run B (original v1.1), `fig5_patch_compare.sh` gives an interchromosomal winner (DAD h2 chr3) in 10 of the 10 tract windows (136,168,000-136,188,000), with same-chromosome minus interchromosomal identity from -0.021347 to -0.000250.

## Reproduce

```shell
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131   # norm_v11/, inputs_v11/, ont/ (ONT BAMs, ont_map.sbatch)
sbatch -o logs/patch9.%j.log patch_dad_chr9.sbatch $D $D/patch_v11 v11
sbatch -o logs/patchchk.%j.log patch_read_check.sbatch $D $D/patch_v11 v11
sbatch -o logs/patch9.%j.log patch_dad_chr9.sbatch $D $D/patch_v11 v11 $D/patch_v11/check/wrong_positions.patched.chr9q_h2_end.v11.tsv   # fill N
sbatch -o logs/patchchk.%j.log patch_read_check.sbatch $D $D/patch_v11 v11
# Fig. 5 steps 4-7 (scripts/pedigree/fig5_stageA), PAN027pat only, on the patched target
CID=PAN027pat_vs_PAN011_joint; IN=$D/patch_v11/inputs_v11_patched; O=$D/patch_v11/out_P_v11
A=$(sbatch --parsable -p tux -w tux05 -c 48 --mem=200G -o logs/P.aln.%j.log ../fig5_stageA/fig5_align.sbatch $CID $IN $O)
sbatch -p tux -w tux05 -c 48 --mem=480G -d afterok:$A -o logs/P.sim.%j.log ../fig5_stageA/fig5_similarity.sbatch $CID $IN $O/raw_paf/$CID.frequency32.paf.gz $O
```

The patched target differs from `inputs_v11/PAN027pat_vs_PAN011_joint.target.fa` only in `PAN011#joint#h2_chr9` (same contig order, other lengths checked). The query is the same file.

## Scripts

| file | what it does |
|---|---|
| `ont_map.sbatch`, `ont_submit.sh`, `sra_verify.sbatch`, `bam_check.sbatch` | download one ONT run from SRA, map it with minimap2 to a diploid (own or cross), check the BAM |
| `assembly_error_scan.sbatch` | read-vs-assembly wrong positions and depth at the chr9q ends; assembly-vs-assembly differences |
| `igv_slices.sbatch`, `screenshot.py` | BAM slices and IGV-style views of the 9q/3q tract (`screenshots/9q3q_*.png`) |
| `synteny_chr1_chr9.sbatch`, `synteny_svbyeye.R` | DAD chr9 h2 vs chr1 h1 synteny with SVbyEye (`screenshots/synteny_*.png`) |
| `patch_dad_chr9.sbatch`, `patch_read_check.sbatch` | patched DAD chr9q end and the read check |
| `fig5_patch_compare.sh` | Fig. 5 class winners on SELF.p chr9q end, Run B (original v1.1) against the patched run |
| `diag_sites.sbatch`, `parent_sites.sbatch`, `read_calls.py`, `ont_site_check.sh`, `read_origin.sbatch`, `read_panel.sbatch`, `cross_igv.sh` | diagnostic sites (own copy vs donor copy) and per-read calls at those sites |
| `where_else.sbatch`, `tip_check.sbatch`, `synteny_extract.sbatch` | where a region occurs in a diploid; chromosome-end comparisons |
| `call_regions.py`, `run_sensitivity.py` | Fig. 5 call regions and run-calling sensitivity, v1.1 vs v1.3.1 |
| `flagger_v11/` | flagger (HiFi, ONT) and nucflag BEDs of the v1.1 release, for the four samples |
| `sites/` | diagnostic-site tables |
