#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 32
#SBATCH --mem=240G
#SBATCH -t 2-00:00:00
# Fig5 step 4, minimap2 on 1 Mb query chunks
# args: CID, IN_DIR, OUT_DIR, MINIMAP2, [CHUNK_BP] (default 1000000)
# Whole-chromosome queries let one long primary chain hide every secondary (-p is relative to the primary),
# which removed PAR1. Here each query contig is cut into CHUNK_BP pieces, aligned with
# minimap2 -c --eqx -x asm20 -N 100 -p 0.1 --secondary=yes -I 8G, and query coordinates are shifted back.
# Coverage check as fig5_align_f.sbatch. out: OUT_DIR/raw_paf/CID.mm2c.paf.gz, CID.mm2c.coverage.tsv
set -euo pipefail
CID=$1; IN=$2; OUT=$3; MM=$4; CH=${5:-1000000}
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
W=/scratch/$USER.fig5mm2c.$SLURM_JOB_ID; mkdir -p "$W" "$OUT/raw_paf"; trap 'rm -rf "$W"' EXIT
log () { echo "[$(date -Is)] $*"; }
log "host=$(hostname) cid=$CID chunk=$CH minimap2 $($MM --version), scratch $(df -h --output=avail /scratch | tail -1)"
cp "$IN/$CID.query.fa" "$IN/$CID.query.fa.fai" "$IN/$CID.target.fa" "$W/"; cd "$W"
awk -v c=$CH '{for (s=1; s<=$2; s+=c) {e=s+c-1; if (e>$2) e=$2; print $1":"s"-"e}}' $CID.query.fa.fai > regions.txt
$S faidx $CID.query.fa -r regions.txt > chunks.fa
log "query chunks: $(wc -l < regions.txt)"
$MM -c --eqx -x asm20 -N 100 -p 0.1 --secondary=yes -I 8G -t $SLURM_CPUS_PER_TASK $CID.target.fa chunks.fa 2> mm2.log \
  | awk -F'\t' -v OFS='\t' 'NR==FNR {len[$1]=$2; next} {n=split($1,a,":"); r=a[n]; name=substr($1,1,length($1)-length(r)-1); split(r,b,"-"); o=b[1]-1; $1=name; $2=len[name]; $3+=o; $4+=o; print}' $CID.query.fa.fai - > raw.paf \
  || { log "ERROR minimap2"; tail -5 mm2.log; exit 1; }
N=$(wc -l < raw.paf); log "raw PAF lines=$N; $(grep -E 'Real time' mm2.log | tail -1)"
cut -f1,3,4 raw.paf | sort -S 8G -k1,1 -k2,2n \
  | awk -F'\t' 'NR==FNR {len[$1]=$2; next} {if ($1!=c) {if (c!="") u[c]=t; c=$1; t=0; e=0} if ($3>e) {s=($2>e)?$2:e; t+=$3-s; e=$3}} END {u[c]=t; bad=0; for (k in len) {f=(k in u)?u[k]/len[k]:0; printf "%s\t%d\t%d\t%.3f\n", k, len[k], u[k]+0, f; if (len[k]>=10000000 && f<0.5) bad=1}; exit bad}' \
    $CID.query.fa.fai - | sort -V > cov.tsv || { cat cov.tsv; log "ERROR contig coverage below 0.5"; exit 1; }
log "coverage ok"
pigz -p 16 -c raw.paf > $CID.mm2c.paf.gz; gzip -t $CID.mm2c.paf.gz
cp $CID.mm2c.paf.gz "$OUT/raw_paf/" && cp cov.tsv "$OUT/raw_paf/$CID.mm2c.coverage.tsv"
log "DONE $OUT/raw_paf/$CID.mm2c.paf.gz lines=$N"
