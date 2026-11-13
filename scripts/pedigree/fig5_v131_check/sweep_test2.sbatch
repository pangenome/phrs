#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 12:00:00
# Which long records sweepga ani drops
# args: WASHU, OUT_DIR
# Lists records >= 10 Mb of PAN027pat v1.1 raw f32 and f48 PAFs missing after sweepga 10:10 --scoring ani;
# keeps the chrX-query f48 subset and its filtered output as a minimal reproducer (fields 1-12 only).
set -euo pipefail
D=$1; OUT=$2; cd $D; mkdir -p $OUT
SW=/home/erikg/.cargo/bin/sweepga; c=PAN027pat_vs_PAN011_joint
W=/scratch/$USER.sweeptest2.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
long () { awk -F'\t' -v OFS='\t' '$11>=10000000 {print $1,$3,$4,$5,$6,$8,$9,$11,$10/$11}' $1 | sort; }
for F in 32 48; do
  r=regen_B_v11/raw_paf/$c.frequency32.paf.gz; [ $F = 48 ] && r=fig5_f48/out_v11/raw_paf/$c.frequency48.paf.gz
  zcat $r > $W/raw.paf
  $SW --num-mappings 10:10 --scaffold-jump 0 --scoring ani --temp-dir $W --output-file $W/ani.paf $W/raw.paf > $OUT/f$F.ani.log 2>&1
  long $W/raw.paf > $W/raw.long; long $W/ani.paf > $W/ani.long
  comm -23 $W/raw.long $W/ani.long > $OUT/f$F.long_records_dropped.tsv
  log "f$F: long records raw $(wc -l < $W/raw.long), after ani $(wc -l < $W/ani.long), dropped $(wc -l < $OUT/f$F.long_records_dropped.tsv)"
  if [ $F = 48 ]; then
    awk -F'\t' '$1=="PAN027#2#chrX"' $W/raw.paf | cut -f1-12 > $OUT/repro.f48_chrX.raw.paf
    $SW --num-mappings 10:10 --scaffold-jump 0 --scoring ani --temp-dir $W --output-file $OUT/repro.f48_chrX.ani.paf $OUT/repro.f48_chrX.raw.paf > $OUT/repro.log 2>&1
    log "repro (fields 1-12, no tags): raw $(wc -l < $OUT/repro.f48_chrX.raw.paf), ani $(wc -l < $OUT/repro.f48_chrX.ani.paf), 153 Mb record kept: $(awk -F'\t' '$11>100000000' $OUT/repro.f48_chrX.ani.paf | wc -l)"
    log "identity-1 records overlapping chrX 15.79-15.80 Mb in repro raw: $(awk -F'\t' '$3<15798000 && $4>15796000 && $10==$11' $OUT/repro.f48_chrX.raw.paf | wc -l), all: $(awk -F'\t' '$3<15798000 && $4>15796000' $OUT/repro.f48_chrX.raw.paf | wc -l)"
  fi
done
log done
