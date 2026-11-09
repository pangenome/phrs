#!/bin/bash
#SBATCH -p workers
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 4:00:00
# Query chrom tip vs 16 same-chrom tips
# args: WASHU dir, OUT tsv, QHAP, CHR, QSTART
set -euo pipefail
WASHU=$1; OUT=$2; QHAP=$3; CHR=$4; QSTART=$5
SAMTOOLS=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
MM2=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
W=/scratch/$USER.tip.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
TAIL=3000000; PIECE=10000
log () { echo "[$(date -Is)] $*"; }

# Queries: v1.1 and v1.3.1 tips, 10 kb pieces
for v in v11 v131; do
  f=$WASHU/norm_$v/$QHAP.fa.gz
  len=$(awk -v c="$QHAP#$CHR" '$1==c{print $2}' $f.fai)
  log "query $v $QHAP#$CHR len=$len"
  for ((s=QSTART; s<len; s+=PIECE)); do
    e=$((s+PIECE)); ((e>len)) && e=$len
    $SAMTOOLS faidx "$f" "$QHAP#$CHR:$((s+1))-$e" | sed "1s/.*/>q_${v}_$s/"
  done > $W/q_$v.fa
  log "pieces $v $(grep -c '>' $W/q_$v.fa)"
done

# Refs: last 3 Mb of CHR, 8 haplotypes x 2 versions
for v in v11 v131; do for f in $WASHU/norm_$v/*.fa.gz; do
  h=$(basename $f .fa.gz); len=$(awk -v c="$h#$CHR" '$1==c{print $2}' $f.fai)
  [ -z "$len" ] && { log "WARN no $CHR in $f"; continue; }
  $SAMTOOLS faidx "$f" "$h#$CHR:$((len-TAIL+1))-$len" | sed "1s/.*/>${h}_${v}/" > "$W/r_${h}_$v.fa"
done; done
log "refs $(ls $W/r_*.fa | wc -l)"

# Map pieces to each ref separately
echo -e "query_version\tpiece_start\tref\tref_version\tidentity\taligned_bp" > $W/out.tsv
for r in $W/r_*.fa; do
  rn=$(head -1 $r | tr -d '>'); rv=${rn##*_}; rh=${rn%_*}
  for v in v11 v131; do
    $MM2 -x asm5 -c -t 16 --secondary=no $r $W/q_$v.fa 2>/dev/null \
    | awk -v rh="$rh" -v rv="$rv" 'BEGIN{OFS="\t"} {split($1,a,"_"); k=$1; id=$10/$11; if(!(k in best)||id*$11>bid[k]*bal[k]){best[k]=1;bid[k]=id;bal[k]=$11;qv[k]=a[2];qs[k]=a[3]}} END{for(k in best) print qv[k],qs[k],rh,rv,sprintf("%.5f",bid[k]),bal[k]}' >> $W/out.tsv
  done
  log "done $rn"
done
cp $W/out.tsv $OUT
log "rows $(($(wc -l < $OUT)-1))"
