#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05
# NUMT test with CHM13 chrM
# args: WASHU, VERSION (v11), LOCI_TSV (sample, contig, start, end, label; 1-based), OUT_DIR
# A. CHM13 chrM (CP068254.1, mito/chm13_chrM.fa) on each of the 8 haplotypes: minimap2 2.31 -cx asm20 -N 500 -p 0.01,
#    hits >= 300 bp -> numts.<hap>.tsv (contig, start, end, chrM start, end, identity).
# B. ONT reads (own, ont_full) overlapping each locus: read length, aligned span, de:f, flag; flank 50-100 kb away as control.
# C. Those reads realigned to locus +- 100 kb plus chrM (map-ont): primary target per read, and SNV sites in the locus
#    before (locus +- 100 kb alone) and after (with chrM), bcftools call --ploidy 1, QUAL >= 20, TYPE only (no alleles).
# out: OUT_DIR/numts.*.tsv, OUT_DIR/reads.<label>.tsv, OUT_DIR/summary.tsv
set -euo pipefail
D=$1; V=$2; L=$3; O=$4; cd $D; mkdir -p $O
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2
B=/gnu/store/z0030n89k5wdgz9s3lr54rwxnh9gqgzj-bcftools-1.21/bin/bcftools; M=$D/mito/chm13_chrM.fa; T=$SLURM_CPUS_PER_TASK
W=/scratch/$USER.mito.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
for s in PAN010 PAN011 PAN027 PAN028; do for h in 1 2; do
  $MM -cx asm20 -N 500 -p 0.01 --secondary=yes -t $T norm_$V/$s#$h.fa.gz $M 2> /dev/null \
    | awk -F'\t' -v OFS='\t' '$4-$3 >= 300 {print $6, $8, $9, $3, $4, sprintf("%.4f", $10/$11)}' | sort -k1,1V -k2,2n > $O/numts.$s#$h.tsv
  log "chrM on $s#$h: hits >= 300 bp $(wc -l < $O/numts.$s#$h.tsv), bp $(awk '{s+=$3-$2} END {print s+0}' $O/numts.$s#$h.tsv), contigs with hits $(cut -f1 $O/numts.$s#$h.tsv | sort -u | wc -l)"
done; done
printf "label\tset\treads\tlen_median\tlen_max\tlen_le_17kb\tde_median\tto_chrM\tsites_before\tsites_after\n" > $O/summary.tsv
while read -r SM C A Z LAB; do H=${C%#chr*}
  ls ont_full/SRR*.$SM.$V.bam > $W/bams.txt
  for set in locus flank; do
    if [ $set = locus ]; then a=$A; z=$Z; else a=$((Z + 50000)); z=$((Z + 50000 + Z - A)); fi
    for b in $(cat $W/bams.txt); do $S view -F 0x904 $b "$C:$a-$z"; done > $W/$set.sam
    awk -F'\t' -v OFS='\t' -v s=$set '{de="NA"; for (i=12;i<=NF;i++) if ($i ~ /^de:f:/) de=substr($i,6); print s, $1, length($10), $4, $5, de}' $W/$set.sam >> $W/reads.$LAB.tsv
  done
  cut -f2 $W/reads.$LAB.tsv | sort -u > $W/names.txt
  printf "set\tread\tlength\tpos\tmapq\tde\n" | cat - $W/reads.$LAB.tsv > $O/reads.$LAB.tsv
  # C. realign locus reads: window alone vs window + chrM
  ws=$((A > 100000 ? A - 100000 : 1)); we=$((Z + 100000))
  $S faidx norm_$V/$H.fa.gz "$C:$ws-$we" | sed "1s/.*/>win/" > $W/win.fa; cat $W/win.fa $M > $W/winM.fa; $S faidx $W/win.fa; $S faidx $W/winM.fa
  awk -F'\t' '$1=="locus" {print $2}' $W/reads.$LAB.tsv | sort -u > $W/locus_names.txt
  for b in $(cat $W/bams.txt); do $S view -N $W/locus_names.txt -F 0x900 $b "$C:$A-$Z"; done | awk -F'\t' '!s[$1]++ {print "@"$1"\n"$10"\n+\n"$11}' > $W/locus.fq
  for ref in win winM; do
    $MM -ax map-ont --secondary=no -t $T $W/$ref.fa $W/locus.fq 2> /dev/null | $S sort -o $W/$ref.bam - && $S index $W/$ref.bam
    $B mpileup -X ont -Ou -f $W/$ref.fa -r "win:$((A - ws + 1))-$((Z - ws + 1))" -q 0 --max-depth 5000 $W/$ref.bam 2> /dev/null \
      | $B call -m --ploidy 1 -v -Ou 2> /dev/null | $B view -i 'QUAL>=20' -v snps -H | wc -l > $W/sites.$ref
  done
  toM=$($S view -F 0x904 $W/winM.bam chrM | cut -f1 | sort -u | wc -l)
  for set in locus flank; do
    awk -F'\t' -v s=$set '$1==s' $W/reads.$LAB.tsv | sort -u -k2,2 > $W/r.txt
    n=$(wc -l < $W/r.txt); lm=$(cut -f3 $W/r.txt | sort -n | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]+0}')
    mx=$(cut -f3 $W/r.txt | sort -n | tail -1); le=$(awk -F'\t' '$3 <= 17000' $W/r.txt | wc -l)
    dm=$(awk -F'\t' '$6!="NA" {print $6}' $W/r.txt | sort -g | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]+0}')
    if [ $set = locus ]; then tm=$toM; sb=$(cat $W/sites.win); sa=$(cat $W/sites.winM); else tm=NA; sb=NA; sa=NA; fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $LAB $set $n $lm ${mx:-0} $le $dm $tm $sb $sa | tee -a $O/summary.tsv
  done
  rm -f $W/reads.$LAB.tsv
done < $L
log done
