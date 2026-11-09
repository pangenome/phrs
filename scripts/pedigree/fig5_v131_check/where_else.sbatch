#!/bin/bash
#SBATCH -p workers
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 4:00:00
# Where else a region occurs in a diploid
# args: WASHU, OUT_PAF, VERSION, REGION, SAMPLE
set -euo pipefail
WASHU=$1; OUT=$2; V=$3; REG=$4; SAMPLE=$5
SAMTOOLS=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
MM2=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
W=/scratch/$USER.where.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
h=${REG%%#chr*}
$SAMTOOLS faidx "$WASHU/norm_$V/$h.fa.gz" "$REG" > $W/q.fa
zcat "$WASHU/norm_$V/$SAMPLE#1.fa.gz" "$WASHU/norm_$V/$SAMPLE#2.fa.gz" > $W/ref.fa
log "query $REG ($V), target $SAMPLE diploid"
# all chains, secondary kept
$MM2 -cx asm20 -N 100 -p 0.5 --secondary=yes -t 16 $W/ref.fa $W/q.fa 2> $W/mm2.log > $W/hits.paf
cp $W/hits.paf $OUT
log "hits $(wc -l < $OUT)"
awk -F'\t' '$11>=2000 {de=""; for(i=13;i<=NF;i++) if($i~/^de:f:/) de=substr($i,6); tp=""; for(i=13;i<=NF;i++) if($i~/^tp:A:/) tp=substr($i,6); printf "%s:%d-%d\t%s\t%s:%d-%d\t%d\tid=%.4f\t%s\n", "q",$3,$4,$5,$6,$8,$9,$11,1-de,tp}' $OUT | sort -k4,4nr | head -30
