#!/bin/bash
#SBATCH -p lowmem,386mem,highmem,1tbmem
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 1:00:00
# Reads vs candidate end copies (v1.1)
# input: check/ont/igv/dad_on_self.q20.fa (MAPQ>=20 DAD reads on SELF.p chr9q)
set -euo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
M=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
W=/scratch/$USER.panel.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
for r in "PAN011#2#chr9:133400001-133744561" "PAN011#2#chr3:205200001-205578940" "PAN011#1#chr3:202200001-202554222" "PAN027#2#chr9:135909150-136209149"; do
  h=${r%%#chr*}; $S faidx norm_v11/$h.fa.gz "$r"; done > $W/panel.fa
grep ">" $W/panel.fa
for t in $(grep ">" $W/panel.fa | tr -d ">"); do $S faidx $W/panel.fa "$t" > $W/t.fa
  $M -x map-ont -c --secondary=no -t 8 $W/t.fa check/ont/igv/dad_on_self.q20.fa 2>/dev/null | awk -v t="$t" 'BEGIN{OFS="\t"} {nm=0; for(i=13;i<=NF;i++) if($i~/^NM:i:/) nm=substr($i,6); print $1, t, $4-$3, $11, nm}'
done > check/ont/igv/dad_on_self.q20.panel.tsv
wc -l check/ont/igv/dad_on_self.q20.panel.tsv
