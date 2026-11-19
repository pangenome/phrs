#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 32
#SBATCH --mem=240G
#SBATCH -t 2-00:00:00
# Fig5 step 4 with minimap2 instead of FastGA
# args: CID, IN_DIR, OUT_DIR, MINIMAP2
# minimap2 -c --eqx -x asm20 -N 50 -p 0.5 --secondary=yes -I 8G (whole diploid target in one index);
# every query contig >= 10 Mb must be aligned over >= 50% (union), else exit 1.
# out: OUT_DIR/raw_paf/CID.mm2.paf.gz, CID.mm2.coverage.tsv
set -euo pipefail
CID=$1; IN=$2; OUT=$3; MM=$4
W=/scratch/$USER.fig5mm2.$SLURM_JOB_ID; mkdir -p "$W" "$OUT/raw_paf"; trap 'rm -rf "$W"' EXIT
log () { echo "[$(date -Is)] $*"; }
log "host=$(hostname) cid=$CID in=$IN out=$OUT minimap2 $($MM --version)"
cp "$IN/$CID.query.fa" "$IN/$CID.query.fa.fai" "$IN/$CID.target.fa" "$W/"; cd "$W"
$MM -c --eqx -x asm20 -N 50 -p 0.5 --secondary=yes -I 8G -t $SLURM_CPUS_PER_TASK $CID.target.fa $CID.query.fa > raw.paf 2> mm2.log || { log "ERROR minimap2"; tail -5 mm2.log; exit 1; }
N=$(wc -l < raw.paf); log "raw PAF lines=$N; $(grep -E 'Real time' mm2.log | tail -1)"
cut -f1,3,4 raw.paf | sort -S 8G -k1,1 -k2,2n \
  | awk -F'\t' 'NR==FNR {len[$1]=$2; next} {if ($1!=c) {if (c!="") u[c]=t; c=$1; t=0; e=0} if ($3>e) {s=($2>e)?$2:e; t+=$3-s; e=$3}} END {u[c]=t; bad=0; for (k in len) {f=(k in u)?u[k]/len[k]:0; printf "%s\t%d\t%d\t%.3f\n", k, len[k], u[k]+0, f; if (len[k]>=10000000 && f<0.5) bad=1}; exit bad}' \
    $CID.query.fa.fai - | sort -V > cov.tsv || { cat cov.tsv; log "ERROR contig coverage below 0.5"; exit 1; }
log "coverage ok"
pigz -p 16 -c raw.paf > $CID.mm2.paf.gz 2>/dev/null || gzip -c raw.paf > $CID.mm2.paf.gz
gzip -t $CID.mm2.paf.gz && cp $CID.mm2.paf.gz "$OUT/raw_paf/" && cp cov.tsv "$OUT/raw_paf/$CID.mm2.coverage.tsv"
log "DONE $OUT/raw_paf/$CID.mm2.paf.gz lines=$N"
