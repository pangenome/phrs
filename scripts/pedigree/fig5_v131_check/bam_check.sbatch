#!/bin/bash
#SBATCH -p lowmem
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 12:00:00
# Full read of every ONT BAM
# args: WASHU, OUT_TSV
# quickcheck sees only header and EOF
set -uo pipefail
D=$1; OUT=$2; cd $D
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
log () { echo "[$(date -Is)] $*"; }
log "BAMs $(ls ont/*.bam | wc -l) on $(hostname)"
echo -e "bam\tstatus\trecords\tprimary\tprimary_q20" > $OUT
for b in ont/*.bam; do
  st=OK; n=$($S view -c -@ 2 $b 2>$OUT.err) || st=READ_ERROR
  [ -s $OUT.err ] && st=READ_ERROR
  [ $st = OK ] && { p=$($S view -c -@ 2 -F 0x904 $b); q=$($S view -c -@ 2 -F 0x904 -q 20 $b); } || { n=NA; p=NA; q=NA; }
  echo -e "$(basename $b)\t$st\t$n\t$p\t$q" >> $OUT; log "$(basename $b) $st $n"
done
rm -f $OUT.err
log "done: $(grep -c READ_ERROR $OUT) READ_ERROR of $(($(wc -l < $OUT) - 1))"
