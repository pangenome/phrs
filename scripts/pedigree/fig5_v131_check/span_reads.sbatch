#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 4:00:00
# Reads spanning a region, per BAM set
# args: WASHU, OUT_TSV, then triples BAM_GLOB CONTIG START-END (glob relative to WASHU, quoted)
# Counts MAPQ>=20 primary reads covering START-5 kb to END+5 kb, reads overlapping the region, and reads clipped >= 1 kb inside it.
set -euo pipefail
D=$1; O=$2; shift 2; cd $D
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
printf "bams\tn_bams\tcontig\tregion\tspanning_pm5kb\toverlapping\tclipped_ge_1kb_inside\n" > $O
while [ $# -ge 3 ]; do
  G=$1; C=$2; R=$3; shift 3; a=${R%-*}; z=${R#*-}
  bams=$(ls $G); n=$(echo "$bams" | wc -l)
  for b in $bams; do $S view -q 20 -F 0x904 $b "$C:$a-$z"; done | awk -F'\t' -v a=$a -v z=$z -v OFS='\t' '
    {e=$4; c=$6; while (match(c,/[0-9]+[MDN=X]/)) {e+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH)}
     o++; if ($4<=a-5000 && e>=z+5000) sp++
     if (($6 ~ /^[0-9]+[SH]/ && $6+0>=1000 && $4>a && $4<z) || ($6 ~ /[0-9]+[SH]$/ && e>a && e<z && match($6,/[0-9]+[SH]$/) && substr($6,RSTART)+0>=1000)) cl++}
    END {printf "%d\t%d\t%d\n", sp, o, cl}' | while IFS=$'\t' read sp o cl; do printf "%s\t%d\t%s\t%s\t%d\t%d\t%d\n" "$G" $n "$C" "$R" $sp $o $cl >> $O; done
  echo "[$(date -Is)] $G $C:$R -> $(tail -1 $O | cut -f5-7)"
done
