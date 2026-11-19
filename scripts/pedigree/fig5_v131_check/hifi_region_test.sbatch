#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH --exclude=octopus05
# HiFi vs ONT disagreement in regions
# args: WASHU, REGIONS_TSV (sample, version, region, label), OUT_TSV
# Per region: SNV sites where reads disagree with the assembly (bcftools 1.21 call --ploidy 1, QUAL >= 20, TYPE only),
# for HiFi (hifi/bam, mpileup -B -I) and ONT (ont_full, mpileup -X ont); median depth MAPQ >= 20; no alleles read.
set -euo pipefail
D=$1; R=$2; O=$3; cd $D
B=/gnu/store/z0030n89k5wdgz9s3lr54rwxnh9gqgzj-bcftools-1.21/bin/bcftools; S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
W=/scratch/$USER.hrt.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
printf "sample\tversion\tregion\tlabel\treads\tbams\tsnv\tdepth_median\n" > $O
while read -r SM V REG LAB; do H=${REG%%#chr*}
  for RD in hifi ont; do
    if [ $RD = hifi ]; then ls hifi/bam/*.$SM.$V.bam > $W/b.txt; X="-B -I"; else ls ont_full/SRR*.$SM.$V.bam > $W/b.txt; X="-X ont"; fi
    n=$($B mpileup $X -Ou -f norm_$V/$H.fa.gz -r "$REG" -q 20 --max-depth 2000 -a AD,DP -b $W/b.txt 2> /dev/null \
      | $B call -m --ploidy 1 -v -Ou 2> /dev/null | $B view -i 'QUAL>=20' -v snps -H | wc -l)
    dp=$($S depth -a -Q 20 -r "$REG" -f $W/b.txt | awk '{d=0; for (i=3;i<=NF;i++) d+=$i; print d}' | sort -n | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}')
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $SM $V "$REG" $LAB $RD $(wc -l < $W/b.txt) $n $dp | tee -a $O
  done
done < $R
