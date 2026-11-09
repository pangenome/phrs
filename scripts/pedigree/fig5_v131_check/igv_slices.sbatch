#!/bin/bash
#SBATCH -p lowmem,386mem,highmem,1tbmem
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 1:00:00
# BAM slices + refs for IGV, 9q/3q
# args: WASHU, OUT_DIR, [FLANK_BP] (default 30000, left of the Fig5 tract)
# 4 pairs: DAD/SELF reads on DAD h2 chr9q and SELF.p chr9q, v1.1 and v1.3.1
set -euo pipefail
D=$1; OUT=$2; FL=${3:-30000}; cd $D; mkdir -p $OUT
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
log () { echo "[$(date -Is)] $*"; }
for v in v11 v131; do
  if [ $v = v11 ]; then DW="PAN011#2#chr9:$((133702042 - FL))-133744561"; SW="PAN027#2#chr9:$((136168001 - FL))-136209149"; DT="133702041 133722041"; ST="136168000 136188000"; DC="133682642 133744561"
  else DW="PAN011#2#chr9:$((133702516 - FL))-133745035"; SW="PAN027#2#chr9:$((136175732 - FL))-136216880"; DT="133702515 133722515"; ST="136175731 136195731"; DC="133683116 133745035"; fi
  for pair in "DADreads_on_DAD:$DW:ont/SRR378260[57][0-9].PAN011.$v.bam" "SELFreads_on_DAD:$DW:ont/SRR378260*.PAN027.onPAN011.$v.bam" \
              "DADreads_on_SELF:$SW:ont/SRR378260*.PAN011.onPAN027.$v.bam" "SELFreads_on_SELF:$SW:ont/SRR3782606[89].PAN027.$v.bam ont/SRR3782607[0-3].PAN027.$v.bam"; do
    n=${pair%%:*}; r=$(echo "$pair" | cut -d: -f2-3); g=$(echo "$pair" | cut -d: -f4-)
    for b in $g; do $S view -b -F 0x100 $b "$r" > $OUT/tmp.$(basename $b); done
    $S merge -f -o $OUT/$n.$v.bam $OUT/tmp.*.bam; rm -f $OUT/tmp.*.bam; $S index $OUT/$n.$v.bam
    log "$n $v $r reads $($S view -c $OUT/$n.$v.bam) MAPQ>=20 $($S view -c -q 20 $OUT/$n.$v.bam) from $(echo $g | wc -w) BAM globs"
  done
  # refs: whole contigs, bgzipped, indexed
  for c in "PAN011#2#chr9:PAN011#2" "PAN027#2#chr9:PAN027#2"; do
    ctg=${c%%:*}; f=${c#*:}; $S faidx norm_$v/$f.fa.gz "$ctg" | $BGZ -c > $OUT/$f.chr9.$v.fa.gz; $S faidx $OUT/$f.chr9.$v.fa.gz
  done
  read a b <<< "$DT"; read c d <<< "$DC"; read e f <<< "$ST"
  printf "PAN011#2#chr9\t%s\t%s\twrong in the DAD assembly: copy of chr1p h1 (61.9 kb)\nPAN011#2#chr9\t%s\t%s\t9q/3q exchange called in Fig5\n" $c $d $a $b > $OUT/regions.DAD.$v.bed
  printf "PAN027#2#chr9\t%s\t%s\t9q/3q exchange called in Fig5\n" $e $f > $OUT/regions.SELF.$v.bed
done
ls -la $OUT; log done
