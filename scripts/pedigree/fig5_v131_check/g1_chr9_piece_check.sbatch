#!/bin/bash
#SBATCH -p tux
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 00:15:00
# G1 chr9 pieces vs targeted patch
# args: --washu DIR --run RUN_DIR (g1_locus output) --patch-fa FA.gz --orig-fa FA.gz --contig NAME --from BP --out TSV
# Aligns each candidate piece and the original tail to the targeted-patch tail (minimap2 -cx asm20 --eqx).
set -euxo pipefail
trap 'echo "exit $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --run) R=$2;; --patch-fa) PF=$2;; --orig-fa) OF=$2;; --contig) C=$2;; --from) F=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D R PF OF C F O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2
cd $D; W=/scratch/$USER/tmp/phrs/chr9cmp.$SLURM_JOB_ID; mkdir -p $W; trap 'echo "exit $? at line ${LINENO}"; rm -rf $W' EXIT
PL=$(awk -v c="$C" '$1==c {print $2}' $PF.fai); OL=$(awk -v c="$C" '$1==c {print $2}' $OF.fai)
log "contig $C: patch length $PL, original length $OL, compare from $F"
$S faidx $PF "$C:$F-$PL" > $W/patch.fa; $S faidx $OF "$C:$F-$OL" > $W/orig.fa; $S faidx $R/pieces.fa; cp $R/pieces.fa $W/pieces.fa
echo -e "query\tqlen\tqstart\tqend\ttlen\ttstart\ttend\tmatches_eq\tblock\tNM" > $O
for q in $(grep '^>' $W/pieces.fa | cut -c2- | cut -d' ' -f1) orig; do
  [ $q = orig ] && cp $W/orig.fa $W/q.fa || $S faidx $W/pieces.fa $q > $W/q.fa
  $MM -t 4 -cx asm20 --eqx $W/patch.fa $W/q.fa 2> $W/mm.log | awk -F'\t' -v OFS='\t' -v q=$q '{nm="NA"; for (i=13;i<=NF;i++) if ($i ~ /^NM:i:/) nm=substr($i,6); print q, $2, $3, $4, $7, $8, $9, $10, $11, nm}' >> $O
  log "$q: $(awk -v q=$q '$1==q' $O | wc -l) alignment records"
done
cat $O
