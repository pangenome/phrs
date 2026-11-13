# Subtelomeric PHRs: Figure 5 (WashU pedigree)

## Paths

```shell
dir_base=/home/guarracino/git/phrs
dir_washu=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131

SWEEPGA=/home/erikg/.cargo/bin/sweepga   # sweepga 0.1.1
IMPG=/home/erikg/.cargo/bin/impg         # impg 0.4.1
SAMTOOLS=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
BGZIP=/home/guarracino/.local/bin/bgzip
FIG5_FILTERS=/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/scripts
```

## Tools

```shell
# Already on octopus01, no install
$SWEEPGA --version
    # sweepga 0.1.1
sha256sum $SWEEPGA
    # a0d7ac0c3312080d67de96d85cdcad9ce0c5a7e523897109b7f598c186ab85a6, same as the original run record (sweepga_binary.tsv, 2026-06-21)
$IMPG --version
    # impg 0.4.1
$SAMTOOLS --version | head -1
    # samtools 1.19
# Local stage B: python3 (stdlib), inkscape, Rscript with data.table and ggplot2
```

## Data

```shell
# WashU pedigree assemblies v1.1 and v1.3.1, https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/, downloaded on 2026-09-20
cd $dir_washu
bash $dir_base/scripts/pedigree/fig5_stageA/get_v11.sh    # -> assemblies_v11/, 8 haplotypes
bash $dir_base/scripts/pedigree/fig5_stageA/get_v131.sh   # -> assemblies/, 8 haplotypes

# Original Fig5 raw sweepGA PAFs (v1.1), made by erikg on 2026-06-25
RAW_ORIG=/moosefs/erikg/phrs/.wg-worktrees/agent-2762/paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/bgzf_raw_paf
```

## Inputs

Rename headers to PanSN `SAMPLE#H#chrN` (hap1 and mat are H=1, hap2 and pat are H=2), keep chr1-22, X, Y:

```shell
cd $dir_washu
bash $dir_base/scripts/pedigree/fig5_stageA/normalize.sh assemblies_v11 norm_v11
bash $dir_base/scripts/pedigree/fig5_stageA/normalize.sh assemblies norm_v131
    # 23 records per haplotype, 184 contigs per version
```

Build query (one child haplotype) and joint target (both parent haplotypes, headers `PARENT#joint#hH_chrN`):

```shell
bash $dir_base/scripts/pedigree/fig5_stageA/build_inputs.sh norm_v11 inputs_v11
bash $dir_base/scripts/pedigree/fig5_stageA/build_inputs.sh norm_v131 inputs_v131
    # 3 comparisons per version, query 23 records, target 46 records
# Check against the original Fig5 input
md5sum /moosefs/erikg/phrs/.wg-worktrees/agent-2639/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_updated_bin/inputs/PAN027pat_vs_PAN011_joint.query.fa inputs_v11/PAN027pat_vs_PAN011_joint.query.fa
    # 19f7e4bc1acaa0b58b3b12cb4dfe1292 for both
```

## Similarity

Fig5 steps 5-7 (sweepGA 10:10 filter, 2 kb depth windows, impg class winners) on the original raw PAFs, in node-local /scratch:

```shell
cd $dir_washu && mkdir -p logs
for CID in PAN027pat_vs_PAN011_joint PAN027mat_vs_PAN010_joint PAN028mat_vs_PAN027_joint; do
  sbatch -J f5A-$CID -o logs/A.$CID.%j.log $dir_base/scripts/pedigree/fig5_stageA/fig5_similarity.sbatch $CID $dir_washu/inputs_v11 $RAW_ORIG/$CID.sweepga_frequency32_many_many_j0.bgzf.paf.gz $dir_washu/out_A_orig
done
    # jobs 2861005 2861010 2861015 on highmem, 48 CPUs, 00:34:06 00:35:55 00:30:49
    # class winner lines 352941 342316 340454, 0 skipped windows, same as data/fig5_*
    # sorted content md5 equal to data/fig5_* for all 3; BED and filtered PAF equal to fig5_pre_impg_depth_filtered_similarity
```

Fig5 steps 4-7 from the v1.1 FASTAs (Run B): FastGA on workers, then steps 5-7 on highmem:

```shell
for CID in PAN027pat_vs_PAN011_joint PAN027mat_vs_PAN010_joint PAN028mat_vs_PAN027_joint; do
  B1=$(sbatch --parsable -J f5B-aln-$CID -o logs/B.aln.$CID.%j.log $dir_base/scripts/pedigree/fig5_stageA/fig5_align.sbatch $CID $dir_washu/inputs_v11 $dir_washu/out_B_v11)
  sbatch -J f5B-sim-$CID -d afterok:$B1 -o logs/B.sim.$CID.%j.log $dir_base/scripts/pedigree/fig5_stageA/fig5_similarity.sbatch $CID $dir_washu/inputs_v11 $dir_washu/out_B_v11/raw_paf/$CID.frequency32.paf.gz $dir_washu/out_B_v11
done
    # FastGA 11:57:19 05:44:10 12:21:39 on 32 CPUs (pat, mat, PAN028); steps 5-7 00:33:08 00:32:29 00:28:19
    # raw PAF lines 1453788 1469498 1591756, class winner lines 352941 342316 340454, all equal to Run A
    # sorted md5 of class winners, filtered PAF and BED equal to Run A for all 3
```

Same commands on the v1.3.1 FASTAs (Run C): `inputs_v131`, `out_C_v131`, job prefix `f5C`:

```shell
    # FastGA 11:50:54 11:11:35 14:32:36 (pat, mat, PAN028); steps 5-7 00:33:18 00:36:29 00:28:05
    # class winner lines 354277 342338 337689 vs 352941 342316 340454 for v1.1
```

## Plots

Fig5 panels C, D, E from the vendored class winners (stage B), run on the local machine:

```shell
cd $dir_base
bash submission/scripts/figures/make_fig5_panels.sh
Rscript submission/scripts/figures/make_fig5_zoom.R       # panel B
Rscript submission/scripts/figures/make_fig5_pedigree.R   # panel A, hand-drawn topology
# Panels A-E in a temp tree from one run's class winners (scp from out_*/outputs), then pixel diff against the paper
bash scripts/pedigree/fig5_stageA/fig5_panels_rebuild.sh <tmp_run_dir> $dir_base/data
    # Run A and Run B: 0 differing pixels in all 5 panels: A of 408,900, B of 1,296,000, C D E of 2,131,200
# Run C uses the v1.3.1 .fai from inputs_v131, renamed fig5_<CID>.{query,target}.fa.fai
bash scripts/pedigree/fig5_stageA/fig5_panels_rebuild.sh <tmp_runC_dir> <tmp_runC_dir>/fai
    # differing pixels: A 0, B 48,778, C 769,916, D 681,626, E 1,016,740
    # drawn runs pat 43 vs 52, mat 25 vs 25, PAN028 91 vs 104 (v1.3.1 vs v1.1)
```

## Checks

Query chromosome tip in 10 kb pieces, minimap2 asm5 against the last 3 Mb of the same chromosome in 8 haplotypes x 2 versions:

```shell
cd $dir_washu
sbatch -J tip-chr5 -o "logs/tip.PAN027_2.chr5.%j.log" $dir_base/scripts/pedigree/fig5_v131_check/tip_check.sbatch $dir_washu check/tip.PAN027_2.chr5.tsv "PAN027#2" chr5 181000000
sbatch -J tip-chr9 -o "logs/tip.PAN027_2.chr9.%j.log" $dir_base/scripts/pedigree/fig5_v131_check/tip_check.sbatch $dir_washu check/tip.PAN027_2.chr9.tsv "PAN027#2" chr9 135200000
sbatch -J tip-chr5b -o "logs/tip.PAN028_1.chr5.%j.log" $dir_base/scripts/pedigree/fig5_v131_check/tip_check.sbatch $dir_washu check/tip.PAN028_1.chr5.tsv "PAN028#1" chr5 181000000
    # PAN027#2#chr5 181.0 Mb to end vs PAN011 h2: identity 1.0000 in every 200 kb bin for v1.1 and v1.3.1
```

Diagnostic sites in a called tract: query base equal to own copy or to donor copy (minimap2 asm20 --cs), one job per tract and version:

```shell
sbatch -J d5v11 -o logs/diag.chr5a.v11.%j.log $dir_base/scripts/pedigree/fig5_v131_check/diag_sites.sbatch $dir_washu check/diag.chr5a.v11.tsv v11 "PAN027#2#chr5:182052001-182080000" "PAN011#2#chr5:181500001-181932123" "PAN011#2#chr1:1-400000"
sbatch -J d9v11 -o logs/diag.chr9t.v11.%j.log $dir_base/scripts/pedigree/fig5_v131_check/diag_sites.sbatch $dir_washu check/diag.chr9t.v11.tsv v11 "PAN027#2#chr9:136168001-136188000" "PAN011#2#chr9:133400001-133744561" "PAN011#2#chr3:205200001-205578940"
    # chr5q tract 28 kb: query like own 48, like donor 1; chr5q tract 12 kb: 16 and 0
    # chr9q tract 20 kb: query like own 9, like donor 152, neither 19
    # v1.3.1 (coordinates +19 on chr5, +7731 on chr9): same counts
```

Drawn runs under 7 run-calling rules (paper, gap1k, gap2k, gap3k, gap4k, anyhap, anyhap4k), one class-winner file per label:

```shell
python3 $dir_base/scripts/pedigree/fig5_v131_check/run_sensitivity.py $dir_base/submission/scripts/figures/fig5_ribbon_plotter.py v1.1=runB/PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz v1.3.1=runC/PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz
    # 9q/3q bp v1.1 vs v1.3.1: paper 20000 vs 10000, gap1k 20000 vs 10000, gap2k 26000 vs 10000, gap3k 26000 vs 10000, gap4k 30000 vs 22000, anyhap 20000 vs 22000, anyhap4k 30000 vs 30000
    # 5q/1p bp v1.1 vs v1.3.1: paper 40000 vs 0, anyhap4k 54000 vs 0
```

ONT UL reads of the two children, HG06807 (PAN027) and HG06808 (PAN028), 6 runs each, as .sra from the NCBI S3 bucket (EBI FASTQ gives 1 MB/s, S3 gives 30 MB/s). Each job downloads one run to /scratch, checks the NCBI md5, keeps the .sra in `ont_sra/<HG>/`, streams it with fastq-dump 2.10.6 and maps it to the child diploid v1.1 and v1.3.1 (minimap2 2.24 map-ont), keeping reads in the call regions:

```shell
# Call regions: all Fig5 drawn runs of v1.1 and v1.3.1, 1 Mb flank, partner arm end 3 Mb, both haplotypes
S=$dir_base/scripts/pedigree/fig5_v131_check/call_regions.py
for v in v11 v131; do
  python3 $S ont/beds/PAN027.$v.bed fai/$v.PAN027#1.fai fai/$v.PAN027#2.fai fai/v11.PAN027#2.fai runB/PAN027pat_vs_PAN011_joint.whole_genome_ribbon_runs.tsv fai/v11.PAN027#1.fai runB/PAN027mat_vs_PAN010_joint.whole_genome_ribbon_runs.tsv fai/v131.PAN027#2.fai runC/PAN027pat_vs_PAN011_joint.whole_genome_ribbon_runs.tsv fai/v131.PAN027#1.fai runC/PAN027mat_vs_PAN010_joint.whole_genome_ribbon_runs.tsv
  python3 $S ont/beds/PAN028.$v.bed fai/$v.PAN028#1.fai fai/$v.PAN028#2.fai fai/v11.PAN028#1.fai runB/PAN028mat_vs_PAN027_joint.whole_genome_ribbon_runs.tsv fai/v131.PAN028#1.fai runC/PAN028mat_vs_PAN027_joint.whole_genome_ribbon_runs.tsv
done
    # PAN027 BED also gets the donor side of PAN028 calls (d: spec); PAN010 and PAN011 BEDs are donor side only:
    #   python3 $S ont/beds/PAN011.$v.bed fai/$v.PAN011#1.fai fai/$v.PAN011#2.fai "d:fai/v11.PAN011#1.fai,fai/v11.PAN011#2.fai:runB/PAN027pat_vs_PAN011_joint.whole_genome_ribbon_runs.tsv" "d:fai/v131.PAN011#1.fai,fai/v131.PAN011#2.fai:runC/PAN027pat_vs_PAN011_joint.whole_genome_ribbon_runs.tsv"
    # bp per BED: PAN027 133,035,710; PAN028 117,611,126; PAN011 101,613,146; PAN010 50,029,524
# ont_runs.tsv: SRR, sample, HG, md5, bytes, S3 URL (NCBI SDL locate API, 2026-10-16)
bash $dir_base/scripts/pedigree/fig5_v131_check/ont_submit.sh $dir_washu
    # children, 876,704,093,541 bytes
bash $dir_base/scripts/pedigree/fig5_v131_check/ont_submit.sh $dir_washu ont_runs_grandparents.tsv
    # HG06803 = PAN010 (ENA sex female), HG06804 = PAN011 (male), 789,077,400,421 bytes
    # a failed job is resubmitted with the same line in a runs tsv; a kept .sra is reused
```

ONT alleles at SNV diagnostic sites of a tract, one BAM (srun on workers):

```shell
bash $dir_base/scripts/pedigree/fig5_v131_check/ont_site_check.sh check/ont/chr9t.v11.tsv ont/SRR37826069.PAN027.v11.bam check/diag.chr9t.v11.tsv "PAN027#2#chr9" 136168001 136188000
    # SRR37826069 only: chr9 tract 140 of 140 sites majority assembly base, median depth 5; chr5 28 kb 33 of 33, depth 12; chr5 12 kb 15 of 15, depth 8; same on v1.3.1
```

Cross alignments for each exchange pair (PAN011-PAN027, PAN010-PAN027, PAN027-PAN028): child reads on the parent diploid and parent reads on the child diploid, v1.1 and v1.3.1, reusing the kept .sra. Column 7 of the runs tsv is the reference sample; BAMs are named `<SRR>.<reads>.on<REF>.<v>.bam`:

```shell
cat ont_runs.tsv ont_runs_grandparents.tsv | awk -F'\t' 'BEGIN{OFS="\t"} $2=="PAN027"{print $0,"PAN011"; print $0,"PAN010"; print $0,"PAN028"} $2=="PAN011"||$2=="PAN010"||$2=="PAN028"{print $0,"PAN027"}' > ont_runs_cross.tsv
bash $dir_base/scripts/pedigree/fig5_v131_check/ont_submit.sh $dir_washu ont_runs_cross.tsv
    # 37 jobs: PAN027 on PAN011, PAN010, PAN028 (6 each); PAN011 on PAN027 (7); PAN010 on PAN027 (6); PAN028 on PAN027 (6)
```

Synteny of chromosome-end windows (500 kb), all-vs-all within each version (minimap2 asm20 -P), one panel per version:

```shell
sbatch -J synteny -o logs/synteny.5q1p.%j.log $dir_base/scripts/pedigree/fig5_v131_check/synteny_extract.sbatch $dir_washu check/synteny_5q1p 500000 "PAN011#1#chr1:p" "PAN011#2#chr1:p" "PAN011#1#chr5:q" "PAN011#2#chr5:q" "PAN027#2#chr5:q"
python3 $dir_base/scripts/pedigree/fig5_v131_check/synteny_plot.py synteny_5q1p.png 500000 "PAN011#1#chr1|p,PAN011#2#chr1|p,PAN011#1#chr5|q,PAN011#2#chr5|q,PAN027#2#chr5|q" v1.1=synteny.v11.paf v1.3.1=synteny.v131.paf --tract "PAN027#2#chr5|q:182052000-182080000"
    # job 2862388, 3708 alignments per version; v1.1 and v1.3.1 PAFs equal after renaming the PAN027 window (+19 bp)
# with MOM: 10 windows (DAD and MOM chr1p and chr5q h1/h2, SELF.m and SELF.p chr5q) -> check/synteny_5q1p_all, 16,701 alignments per version
# rows of tracks (ribbons between consecutive tracks), columns = versions; 16:9 output; rows 4-5 = MOM chr5q/chr1p chains
python3 $dir_base/scripts/pedigree/fig5_v131_check/synteny_pairs.py synteny_5q1p_pairs.png 500000 "PAN011#2#chr5|q,PAN027#2#chr5|q,PAN011#2#chr1|p;PAN011#1#chr5|q,PAN011#2#chr5|q,PAN011#2#chr1|p,PAN011#1#chr1|p;PAN011#2#chr1|p,PAN011#1#chr5|q,PAN011#1#chr1|p,PAN011#2#chr5|q" v1.1=synteny.v11.paf v1.3.1=synteny.v131.paf --note "<legend>" --tract "v1.1@PAN027#2#chr5|q:182052000-182080000" "v1.3.1@PAN027#2#chr5|q:182052019-182080019"
    # tract lifted (both versions): DAD h2 chr5q 181,839,879-181,867,879; h1 chr5q 185,005,781-185,033,557; h2 chr1p 80,597-106,596; h1 chr1p 77,044-101,785
    # DAD h1 chr1p vs h2 chr1p: 1 block 421,486 bp, query 76,938-489,288; first 76,938 bp of h1 chr1p have no h2 chr1p match >= 5 kb
    # DAD chr5q vs chr1p (>= 5 kb, >= 0.95): h1 chr5q-h2 chr1p 26 blocks 255,798 bp mean 0.9758; h2 chr5q-h1 chr1p 6 blocks 60,605 bp mean 0.9905; h2 chr5q-h2 chr1p 8 blocks 64,392 bp mean 0.9734
```

Synteny figures for the three events, written into the repo (PAFs stay on moosefs, `check/synteny_<event>_all/`):

```shell
for e in 5q1p 9q3q par1; do mkdir -p paf/synteny_${e}_all; scp "octopus01:$dir_washu/check/synteny_${e}_all/*.paf" paf/synteny_${e}_all/; done
bash $dir_base/scripts/pedigree/fig5_v131_check/make_synteny_figures.sh paf $dir_base/scripts/pedigree/fig5_v131_check/synteny
    # synteny_{5q1p,9q3q,par1}.png (16:9) + .log with lifted tract coordinates and blocks per pair
    # extraction jobs: 9q3q 2862391 (7,518 alignments per version), par1 2862392
```

Parent ONT reads at the child's diagnostic sites (sites lifted child -> parent copy by minimap2 asm20 --cs):

```shell
sbatch -J ps9v11 -o logs/psites.chr9t.v11.%j.log $dir_base/scripts/pedigree/fig5_v131_check/parent_sites.sbatch $dir_washu check/ont/dad.chr9t.v11.tsv v11 "PAN027#2#chr9:136168001-136188000" "PAN011#2#chr9:133400001-133744561" check/diag.chr9t.v11.tsv "ont/SRR378260*.PAN011.v11.bam"
    # 134 of 180 sites lifted; class donor (SELF has chr3 base): 122 sites, DAD reads majority DAD-assembly base 109, SELF base 3 (v1.3.1: 111 and 3)
```

Plot of child vs parent reads at the 9q/3q sites (per-site TSVs copied into the repo):

```shell
cd $dir_base/scripts/pedigree/fig5_v131_check
python3 plot_sites.py sites/ont_sites_9q3q.png v1.1:sites/all.chr9t.v11.tsv:sites/dad.chr9t.v11.tsv v1.3.1:sites/all.chr9t.v131.tsv:sites/dad.chr9t.v131.tsv
    # 122 sites; SELF reads with SELF base: 105 sites >= 0.9; DAD reads with SELF base: 0.0-0.65 per site, median 0.27
```

IGV-style read view at the 9q/3q sites (reads extracted on the cluster to check/ont/igv/*.sam, samtools view -F 0x904, columns 1-6 and 10):

```shell
# site files sites/sites.<event>.<self|dad>.<v>.tsv: pos, base A, base B (9q3q from dad.chr9t TSVs; 5q1p from diag.chr5a + parent_sites on DAD h2 chr5q)
bash make_read_plots.sh <dir with self/dad *.sam from check/ont/igv>
    # sites/ont_reads_9q3q.png, sites/ont_reads_5q1p.png; 5q1p: SELF 145 of 152 and DAD 92 of 96 reads mostly chr5 base, 0 mostly chr1
    # DAD: 22 reads (v1.1); MAPQ >= 50 reads (4) carry SELF's base, MAPQ 0-1 reads carry the DAD-assembly base
# same sites with a MAPQ filter, and on DAD chr3q (all lifted sites, matches included)
sbatch -J ps --export=ALL,MINQ=20 -o logs/psites.q20.chr9.%j.log $dir_base/scripts/pedigree/fig5_v131_check/parent_sites.sbatch $dir_washu check/ont/dad.q20.chr9h2.v11.tsv v11 "PAN027#2#chr9:136168001-136188000" "PAN011#2#chr9:133400001-133744561" check/diag.chr9t.v11.tsv "ont/SRR378260*.PAN011.v11.bam"
    # MAPQ >= 20: SELF base majority at 118 of 127 donor sites (median depth 2); depth on DAD h2 chr3q 110, h1 chr3q 77, h2 chr9q 11
# all 4 read-reference pairs (cross reads extracted to check/ont/igv/{dad_on_self,self_on_dad}.<v>.sam)
# output: sites/ont_reads_9q3q_4way.png (third plot in make_read_plots.sh)
    # MAPQ >= 20 reads, majority SELF base / DAD-assembly base (v1.1): DAD on SELF 21/0 (2 of 7 DAD runs), SELF on DAD 5/0, DAD own 5/1, SELF own 88/0
```

Origin of DAD reads at the 9q/3q tract (v1.1; cross reads from cross_igv.sh, per-read calls from read_calls.py with a 5th argument):

```shell
# own-assembly placement of the same reads
sbatch -o logs/read_origin.%j.log scripts_tmp/read_origin.sbatch $dir_washu check/ont/igv/dad_on_self.v11.names check/ont/igv/dad_on_self.v11.own_placement.tsv 'ont/SRR378260[57][0-9].PAN011.v11.bam'
# chr9-vs-chr3 paralog sites (SELF.p chr9q last 300 kb vs DAD h2 chr3q, DAD h1 chr3q) and vs DAD h2 chr9q
sbatch -o logs/diag.para9.v11.%j.log diag_sites.sbatch $dir_washu check/diag.para9.v11.tsv v11 "PAN027#2#chr9:135909150-136209149" "PAN011#2#chr3:205200001-205578940" "PAN011#1#chr3:202200001-202554222"
sbatch -o logs/diag.p9vs9.v11.%j.log diag_sites.sbatch $dir_washu check/diag.p9vs9.v11.tsv v11 "PAN027#2#chr9:135909150-136209149" "PAN011#2#chr9:133400001-133744561" "PAN011#2#chr3:205200001-205578940"
# each MAPQ>=20 read aligned alone to 4 candidate end copies
sbatch -o logs/panel.%j.log read_panel.sbatch
    # 3 DAD reads with MAPQ 60 on DAD h2 chr9q start in chr9q-unique sequence and carry SELF's tract (chr9 base at 5-9 of 9 paralog sites, SELF base at 29/30, 49/50, 118/118 sites)
# v1.3.1: same with dad_on_self.v131.names, 'ont/SRR378260[57][0-9].PAN011.v131.bam', check/diag.para9.v131.tsv (query PAN027#2#chr9:135916881-136216880), tract 136175732-136195731
sbatch -p lowmem -o logs/where.dad9.v131.%j.log where_else.sbatch $dir_washu check/where.dad_chr9tract.v131.paf v131 "PAN011#2#chr9:133702516-133722515" PAN011
    # v1.3.1 gives the same counts and the same 3 anchor reads; DAD h2 chr9q tract = DAD h1 chr1p 22,527-42,527 at 1.0000 in v1.3.1 too
```

Full read of every ONT BAM (quickcheck reads only header and EOF):

```shell
sbatch -o logs/bam_check.%j.log $dir_base/scripts/pedigree/fig5_v131_check/bam_check.sbatch $dir_washu check/ont/bam_check.tsv
```

Where else a region occurs in a diploid (minimap2 asm20, secondary kept):

```shell
sbatch -J where -o logs/where.dad9.v11.%j.log $dir_base/scripts/pedigree/fig5_v131_check/where_else.sbatch $dir_washu check/where.dad_chr9tract.v11.paf v11 "PAN011#2#chr9:133702042-133722041" PAN011
    # DAD h2 chr9q tract = DAD h1 chr1p 22,527-42,527 at 1.0000 over 20,000 bp; 8 more DAD ends at 0.9924-0.9986
sbatch -J where -o logs/where.self9.v11.%j.log $dir_base/scripts/pedigree/fig5_v131_check/where_else.sbatch $dir_washu check/where.self_chr9tract_in_dad.v11.paf v11 "PAN027#2#chr9:136168001-136188000" PAN011
    # SELF tract best: DAD h2 chr3q 0.9988, h1 chr3q 0.9986, h2 chr15q 0.9933; DAD h2 chr9q (assembly) 0.9923
```

Flagger and nucflag annotations (v1.1 only; none published for v1.3.1), drawn as strips on the synteny and read plots:

```shell
cd $dir_base/scripts/pedigree/fig5_v131_check
for s in PAN010 PAN011 PAN027 PAN028; do for f in flagger.$s.hifi.bed flagger.$s.ONT.bed nucflag.hifi.$s.bed; do curl -s -o flagger_v11/$f https://public.gi.ucsc.edu/~mcechova/pedigree/annotations/v1.1/$s/$f; done; done
    # PAN011.chr9.haplotype2 133,680,000-133,744,561 Dup (HiFi); PAN011.chr1.haplotype1 0-48,000 Dup (HiFi)
# make_synteny_figures.sh passes them with --annot; plot_reads.py takes --annot KIND=BED ... and an 8th panel field CONTIG
```

Plain-language explainer figures (Italian text) for the Fig5 events:

```shell
cd $dir_base/scripts/pedigree/fig5_v131_check
# SAM from cross_igv.sh (check/ont/igv/dad_on_self.v11.sam); paralog sites = diag.para9.v11 'neither' SNVs
python3 explain_events.py explain dad_on_self.v11.sam sites/sites.9q3q.para_chr9_vs_chr3.v11.tsv sites/sites.9q3q.self.v11.tsv sites/anchor_reads.9q3q.v11.txt dad.v11.sam dad.v131.sam
    # dad.<v>.sam = DAD reads on DAD's own h2 chr9q tract (check/ont/igv)
    # explain/explain_famiglia.png, explain_9q3q.png, explain_5q1p_par1.png, explain_errore_assembly.png
```

Read-vs-assembly mismatch scan (IGV-style evidence of the DAD chr9q h2 assembly error), last 330 kb of DAD h2 chr9q, first 330 kb of DAD h1 chr1p, SELF.p chr9q control:

```shell
sbatch -o logs/asmerr.%j.log assembly_error_scan.sbatch $dir_washu check/asmerr                          # DAD and SELF chr9q, MAPQ>=20, plus assembly-vs-assembly diffs
sbatch --export=ALL,PART=ends -o logs/asmerr_ends.%j.log assembly_error_scan.sbatch $dir_washu check/asmerr  # DAD chr1p h1 and chr9q h2, MAPQ>=0 and >=20
python3 plot_assembly_error.py explain/explain_errore_letture.png <local copy of check/asmerr>
    # wrong-letter columns (>=3 MAPQ>=20 reads with a base, >=80% one non-assembly base): DAD chr9q h2 106 (v1.1) / 97 (v1.3.1), all inside the 61.9 kb copied zone; DAD chr1p h1 0; SELF chr9q 0
```

BAM slices and references for IGV screenshots of the 9q/3q tract (4 read-reference pairs, v1.1 and v1.3.1; tract centred, 30 kb left flank, right flank to the chr9q end):

```shell
sbatch -o logs/igv_slices.%j.log igv_slices.sbatch $dir_washu check/igv_slices
    # {DAD,SELF}reads_on_{DAD,SELF}.<v>.bam(.bai), PAN011#2.chr9.<v>.fa.gz, PAN027#2.chr9.<v>.fa.gz (bgzip+fai), regions.{DAD,SELF}.<v>.bed; 223 MB
```

English version and IGV-style 4-pair views (no GUI; IGV desktop needs a window and JBrowse jb2export colours every ONT error):

```shell
python3 plot_assembly_error.py explain/explain_read_errors_en.png <local check/asmerr> en
python3 igv_like.py explain/igv/igv_4pairs.<v>[.en].png <local check/igv_slices> <v> [it|en]   # colours a column only where >= 3 MAPQ>=20 reads show a base and >= 80% share one non-reference base
    # disagreeing positions, v1.1 (v1.3.1): DAD reads on DAD 106 (97), SELF reads on DAD 185 (190), DAD reads on SELF 0 (0), SELF reads on SELF 0 (0)
```

Tidy-up 2026-11-05: all earlier figures (explain/, sites/*.png, synteny/) and the plot-only scripts (explain_events.py, plot_assembly_error.py, plot_reads.py, plot_sites.py, make_read_plots.sh, synteny_plot.py, synteny_pairs.py, make_synteny_figures.sh, igv_snapshots.sh, igv_like.py) moved to `scripts/pedigree/fig5_v131_check/archive/{figures,plot_scripts}/`; one-off retry run lists to `archive/run_lists/`. Commands above that call those scripts now run from `archive/plot_scripts/`.

Screenshots of the 9q/3q tract (the current figure set): DAD and SELF ONT reads on the DAD and SELF assemblies, v1.1 and v1.3.1, Fig5 tract +/- 20 kb:

```shell
cd $dir_base/scripts/pedigree/fig5_v131_check
python3 screenshot.py <local copy of check/igv_slices> screenshots 20000
    # screenshots/9q3q_{DAD,SELF}reads_on_{DAD,SELF}.{v11,v131}.png
    # disagreeing positions v1.1 (v1.3.1): DAD on DAD 97 (88), SELF on DAD 177 (182), DAD on SELF 0 (0), SELF on SELF 0 (0)
```

flagger/nucflag tracks (v1.1, from flagger_v11/) are drawn in the screenshots. Synteny of DAD chr9q h2 end vs chr1p h1 start:

```shell
sbatch -o logs/syn19.%j.log synteny_chr1_chr9.sbatch $dir_washu check/syn19 300000   # minimap2 asm20 --cs -P, last 300 kb of chr9 h2 vs first 300 kb of chr1 h1, v1.1 and v1.3.1
python3 synteny_chr1_chr9.py screenshots/synteny_DAD_chr9q_vs_chr1p.png <local copy of check/syn19>
    # 2 kb ribbons by identity; the last ~62 kb of chr9q h2 = first ~62 kb of chr1p h1 at 100%; homology continues to ~195 kb at 99-99.9%
```

Current figure set (2026-11-05), `scripts/pedigree/fig5_v131_check/screenshots/`:

```shell
sbatch -o logs/igv_slices80.%j.log igv_slices.sbatch $dir_washu check/igv_slices80 80000
python3 screenshot.py <local check/igv_slices80> screenshots 60000      # 8 PNGs, 60 kb before the Fig5 tract to the chr9q end, flagger/nucflag tracks (v1.1)
sbatch -o logs/syn19w.%j.log synteny_chr1_chr9.sbatch $dir_washu check/syn19 300000   # whole chr9 h2 vs chr1 h1 (asm20) + 300 kb ends (asm20 --cs -P); 19:37
Rscript synteny_svbyeye.R <local check/syn19> flagger_v11 screenshots/synteny_DAD_chr9_vs_chr1   # SVbyEye 0.99.0 plotMiro, p arm left; .v11.png .v131.png
    # whole chromosomes: 3,216 (v1.1) / 3,218 (v1.3.1) chains, 48 >= 10 kb
```
The matplotlib synteny script (synteny_chr1_chr9.py) and its figure moved to archive/.

## Patched DAD v1.1 chr9q end (2026-11-09)

```shell
cd $dir_washu
sbatch -J patch9 -o logs/patch9.%j.log scripts_tmp/patch_dad_chr9.sbatch $dir_washu $dir_washu/patch_v11 v11
    # job 2871258 on tux05, 00:03:17: kept PAN011#2#chr9 1-133,673,490, segment 69,687 bp (93 N), chr9 133,744,561 -> 133,743,177
sbatch -J patchchk -o logs/patchchk.%j.log scripts_tmp/patch_read_check.sbatch $dir_washu $dir_washu/patch_v11 v11
CID=PAN027pat_vs_PAN011_joint; IN=$dir_washu/patch_v11/inputs_v11_patched; O=$dir_washu/patch_v11/out_P_v11
A=$(sbatch --parsable -p tux -w tux05 -c 48 --mem=200G -J f5P-aln -o logs/P.aln.%j.log fig5_align.sbatch $CID $IN $O)
sbatch -p tux -w tux05 -c 48 --mem=480G -d afterok:$A -J f5P-sim -o logs/P.sim.%j.log fig5_similarity.sbatch $CID $IN $O/raw_paf/$CID.frequency32.paf.gz $O
```

Fill the 93 N from the read check, then recheck and compare Fig5 class winners (2026-11-09):

```shell
F=$dir_washu/patch_v11/check/wrong_positions.patched.chr9q_h2_end.v11.tsv
sbatch -J patch9b -o logs/patch9.%j.log scripts_tmp/patch_dad_chr9.sbatch $dir_washu $dir_washu/patch_v11 v11 $F
    # job 2871265, 00:03:22: filled 93 N, segment md5 159970f65d7b7af57a67345f45b66abb
sbatch -J patchchk -o logs/patchchk.%j.log scripts_tmp/patch_read_check.sbatch $dir_washu $dir_washu/patch_v11 v11
    # job 2871266: patched chr9q h2 end 0 wrong positions, depth 77 vs 72; original 57, 17 vs 72
bash scripts_tmp/fig5_patch_compare.sh $dir_washu > patch_v11/check/fig5_patch_compare.v11.tsv
    # Run B: interchrom (DAD h2 chr3) wins 10 of 10 tract windows
```

## Read-vs-assembly scan of all chromosome ends (2026-11-09)

```shell
for s in PAN011 PAN027 PAN028; do
  sbatch -J endscan-$s -o logs/endscan.$s.%j.log scripts_tmp/end_scan.sbatch $dir_washu $dir_washu/check/end_scan $s v11
done
    # jobs 2871288 2871289 2871290
```

## Whole-genome ONT BAMs (2026-11-09)

```shell
cd $dir_washu
TUX2_DEP=2871268 bash scripts_tmp/ont_submit_full.sh $dir_washu $dir_washu/ont_full > logs/ontfull.submit.txt
    # 100 jobs 2871345-2871444: own v11, cross v11, own v131, cross v131; 124 BAMs; lanes tux05 x2, octopus10, octopus11
```

## FastGA stop on SELF.p chr5 (2026-11-10)

```shell
D=$dir_washu; scp scripts/pedigree/fig5_v131_check/fastga_chr5_stop.sbatch octopus01:$D/scripts_tmp/
ssh octopus01 sbatch -D $D -J fgstop -o $D/logs/fgstop.%j.log $D/scripts_tmp/fastga_chr5_stop.sbatch $D $D/check/fastga_chr5_stop
    # job 2871645; FastGA of the sweepga 0.1.1 cache, query PAN027#2#chr5, -f32 and -f255, both versions
```

## minimap2 on 1 Mb chunks, read scans (2026-11-13)

```shell
D=$dir_washu
ssh octopus01 bash $D/scripts_tmp/fig5_mm2c_submit.sh    # aligns 2872632-2872650 + ani/cnd steps 5-7, out $D/fig5_mm2c/
scp scripts/pedigree/fig5_v131_check/fig5_runs_batch.sh octopus01:/tmp/guarracino.fig5_runs_batch.sh
ssh octopus01 bash /tmp/guarracino.fig5_runs_batch.sh    # drawn runs + comparisons, mm2c part in fig5_mm2c_runs_compare.txt
ssh octopus01 sbatch -D $D -J chr5id -o $D/logs/chr5id.%j.log $D/scripts_tmp/chr5_read_identity.sbatch $D $D/check/self_chr5    # 2872665
ssh octopus01 sbatch -D $D -J wgscan-test --array=124,136 -o "$D/logs/wgscan.%A_%a.log" $D/scripts_tmp/wg_scan.sbatch $D $D/check/wg_scan/tasks.tsv $D/check/wg_scan/bins    # 2872656
ssh octopus01 sbatch -D $D -J wgspeed -o $D/logs/wgspeed.%j.log $D/scripts_tmp/wg_speed.sbatch $D $D/check/wg_speed    # 2872695, mpileup configs
```
