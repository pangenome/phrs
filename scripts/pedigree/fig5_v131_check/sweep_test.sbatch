#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 12:00:00
# Does sweepga 10:10 drop long records?
# args: WASHU, OUT_DIR
# Fig5 step 5 (sweepga --num-mappings 10:10 --scaffold-jump 0 --scoring ani) on PAN027pat v1.1 raw PAFs:
# full f32 (regen_B_v11) and f48; f48 subsets (query chrX; query chrX or target h1_chrX); f48 with default scoring.
# Reports records >= 10 Mb, and the chrX 153 Mb record, in each output. out: OUT_DIR/summary.tsv, logs
set -euo pipefail
D=$1; OUT=$2; cd $D; mkdir -p $OUT
SW=/home/erikg/.cargo/bin/sweepga; c=PAN027pat_vs_PAN011_joint
W=/scratch/$USER.sweeptest.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
log "host $(hostname), $($SW --version)"
zcat regen_B_v11/raw_paf/$c.frequency32.paf.gz > $W/f32.paf
zcat fig5_f48/out_v11/raw_paf/$c.frequency48.paf.gz > $W/f48.paf
awk -F'\t' '$1=="PAN027#2#chrX"' $W/f48.paf > $W/f48_qX.paf
awk -F'\t' '$1=="PAN027#2#chrX" || $6=="PAN011#joint#h1_chrX"' $W/f48.paf > $W/f48_qtX.paf
log "lines: f32 $(wc -l < $W/f32.paf), f48 $(wc -l < $W/f48.paf), f48_qX $(wc -l < $W/f48_qX.paf), f48_qtX $(wc -l < $W/f48_qtX.paf)"
log "chrX record in f48 raw (fields 1-12 and tag names, no cigar):"
awk -F'\t' '$1=="PAN027#2#chrX" && $11>100000000 {o=$1; for (i=2;i<=12;i++) o=o"\t"$i; for (i=13;i<=NF;i++) if ($i !~ /^cg:/) o=o"\t"$i; print o}' $W/f48.paf
stat () { awk -F'\t' -v t=$1 '$11>=10000000 {n++; b+=$11} $1=="PAN027#2#chrX" && $6 ~ /h1_chrX$/ && $11>100000000 {x++} END {printf "%s\t%d\t%d\t%d\t%d\n", t, NR, n, b, x}' $2; }
printf "run\tlines\trecords_ge_10Mb\tbp_ge_10Mb\tchrX_153Mb_record\n" > $OUT/summary.tsv
for i in f32.paf f48.paf f48_qX.paf f48_qtX.paf; do stat raw_${i%.paf} $W/$i >> $OUT/summary.tsv; done
filt () { # in out scoring
  $SW --num-mappings 10:10 --scaffold-jump 0 --scoring $3 --temp-dir $W --output-file $W/$2 $W/$1 > $OUT/$2.log 2>&1 || { log "ERROR sweepga $1 $3"; tail -5 $OUT/$2.log; exit 1; }
  stat ${2%.paf} $W/$2 >> $OUT/summary.tsv; log "done $2: $(tail -1 $OUT/summary.tsv)"
}
filt f32.paf f32.ani.paf ani
filt f48.paf f48.ani.paf ani
filt f48_qX.paf f48_qX.ani.paf ani
filt f48_qtX.paf f48_qtX.ani.paf ani
filt f48.paf f48.loglenani.paf log-length-ani
cat $OUT/summary.tsv
log done
