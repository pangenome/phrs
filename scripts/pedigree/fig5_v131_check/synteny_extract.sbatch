#!/bin/bash
#SBATCH -p workers
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 2:00:00
# Chromosome-end windows, all-vs-all within one version
# args: WASHU, OUT_DIR, WINDOW_BP, then specs HAP#CHR:p|q ...
# e.g. "PAN011#1#chr1:p" "PAN011#2#chr5:q"
set -euo pipefail
WASHU=$1; OUT=$2; WIN=$3; shift 3
SAMTOOLS=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
MM2=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
W=/scratch/$USER.synteny.$SLURM_JOB_ID; mkdir -p $W $OUT; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }

for v in v11 v131; do
  : > $W/$v.fa
  for spec in "$@"; do
    ctg=${spec%:*}; arm=${spec##*:}; hap=${ctg%%#chr*}
    f=$WASHU/norm_$v/$hap.fa.gz
    len=$(awk -v c="$ctg" '$1==c{print $2}' $f.fai)
    [ -n "$len" ] || { log "ERROR $ctg missing in $v"; exit 1; }
    if [ "$arm" = p ]; then s=1; e=$WIN; else s=$((len-WIN+1)); e=$len; fi
    # name keeps contig, arm, window start
    $SAMTOOLS faidx $f "$ctg:$s-$e" | sed "1s/.*/>$ctg|$arm|$s/" >> $W/$v.fa
    log "$v $ctg $arm $s-$e (len $len)"
  done
  # all chains, both strands, self hits dropped later
  $MM2 -cx asm20 -P --eqx -t 16 $W/$v.fa $W/$v.fa 2> $W/$v.mm2.log > $W/$v.paf
  log "$v alignments $(wc -l < $W/$v.paf)"
  cp $W/$v.paf $OUT/synteny.$v.paf
done
log "done"
