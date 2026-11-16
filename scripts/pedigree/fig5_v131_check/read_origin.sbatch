#!/bin/bash
#SBATCH -p workers
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 2:00:00
# Primary placement of named reads
# args: WASHU, NAMES (one per line), OUT_TSV, BAM_GLOB
set -euo pipefail
D=$1; N=$2; OUT=$3; GLOB=$4; cd $D
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
log () { echo "[$(date -Is)] $*"; }
log "names $(wc -l < $N), BAMs $(ls $GLOB | wc -l)"
for b in $GLOB; do $S view -@ 3 -F 0x904 -N $N $b | awk -v b=$(basename $b) 'BEGIN{OFS="\t"} {c=$6; n=0; while (match(c, /[0-9]+[MDN=X]/)) {x=substr(c, RSTART, RLENGTH); n+=substr(x, 1, length(x)-1); c=substr(c, RSTART+RLENGTH)} print $1, $3, $4, $4+n-1, $5, b}'; done > $OUT
log "placed $(cut -f1 $OUT | sort -u | wc -l) of $(wc -l < $N) names"
