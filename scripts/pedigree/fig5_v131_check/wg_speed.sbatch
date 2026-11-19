#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 12:00:00
# Speed and recall of mpileup configs
# args: WASHU, OUT_DIR ; PAN011 v11 own ONT BAMs
# Regions: chr21 20-22 Mb (speed), chr9 h2 last 400 kb (known DAD chr9 assembly error, control).
# Configs: default (as wg_scan.sbatch), -X ont (no indels), -X ont-sup (indels-cns).
# out: OUT_DIR/speed.tsv (config, region, seconds, snv, indel), per 10 kb counts in OUT_DIR/<config>.<region>.bins.tsv
set -euo pipefail
D=$1; O=$2; cd $D; mkdir -p $O
B=/gnu/store/z0030n89k5wdgz9s3lr54rwxnh9gqgzj-bcftools-1.21/bin/bcftools
W=/scratch/$USER.wgspeed.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
ls ont_full/SRR*.PAN011.v11.bam > $W/bams.txt
L9=$(awk '$1=="PAN011#2#chr9" {print $2}' norm_v11/PAN011#2.fa.gz.fai)
printf "config\tregion\tseconds\tsnv\tindel\n" > $O/speed.tsv
for r in "PAN011#2#chr9:$((L9-400000))-$L9" "PAN011#2#chr21:20000001-22000000"; do
  for c in default ont ont-sup; do
    x=""; [ $c = default ] || x="-X $c"
    t0=$(date +%s)
    $B mpileup -Ou -f norm_v11/PAN011#2.fa.gz -r "$r" -q 20 -B --max-depth 2000 -a AD,DP $x -b $W/bams.txt 2> $W/mp.log \
      | $B call -m --ploidy 1 -v -Ou 2> $W/call.log | $B view -i 'QUAL>=20' -Ou | $B query -f '%POS\t%TYPE\n' > $W/q.tsv
    s=$(( $(date +%s) - t0 )); n=$(echo $r | tr '#:' '__')
    awk -F'\t' '{b=int(($1-1)/10000)*10000; if ($2=="SNP") s[b]++; else i[b]++; k[b]=1} END {for (b in k) print b"\t"s[b]+0"\t"i[b]+0}' $W/q.tsv | sort -n > $O/$c.$n.bins.tsv
    printf "%s\t%s\t%s\t%s\t%s\n" $c "$r" $s $(awk '$2=="SNP"' $W/q.tsv | wc -l) $(awk '$2!="SNP"' $W/q.tsv | wc -l) | tee -a $O/speed.tsv
    echo "[$(date -Is)] $c $r done; mpileup log: $(tr '\n' ' ' < $W/mp.log | cut -c1-200)"
  done
done
