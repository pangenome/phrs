#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 6:00:00
# sweepga main vs fix: 10:10 filter only
# args: WASHU, DIR (sweepga.main, sweepga.fix, s9/s10/s11.paf)
# Synthetic cases, chrX repro, full PAN027pat v1.1 raw PAFs (f32 regen, f48) with --scoring ani.
# out: DIR/filter_summary.tsv (binary, input, scoring, kept, long records kept)
set -uo pipefail
D=$1; O=$2; cd $O; log () { echo "[$(date -Is)] $*"; }
W=/scratch/$USER.sweepfilter.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log "host $(hostname), $(./sweepga.main --version), $(./sweepga.fix --version)"
printf "binary\tinput\tscoring\tkept\tlong_kept\tlong_min_bp\tseconds\n" > filter_summary.tsv
run () { # binary input scoring min_len
  t0=$(date +%s)
  ./sweepga.$1 --num-mappings 10:10 --scaffold-jump 0 --scoring $3 --temp-dir $W --output-file $W/o.paf $2 > $W/log 2>&1 || { log "ERROR sweepga $1 $2 $3"; tail -3 $W/log; return; }
  printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\n" $1 $(basename $2) $3 $(wc -l < $W/o.paf) $(awk -F'\t' -v m=$4 '$11>=m' $W/o.paf | wc -l) $4 $(( $(date +%s) - t0 )) >> filter_summary.tsv
  log "$(tail -1 filter_summary.tsv)"
}
zcat $D/regen_B_v11/raw_paf/PAN027pat_vs_PAN011_joint.frequency32.paf.gz > $W/f32.paf
zcat $D/fig5_f48/out_v11/raw_paf/PAN027pat_vs_PAN011_joint.frequency48.paf.gz > $W/f48.paf
log "raw records >= 10 Mb: f32 $(awk -F'\t' '$11>=10000000' $W/f32.paf | wc -l), f48 $(awk -F'\t' '$11>=10000000' $W/f48.paf | wc -l)"
for v in main fix; do
  for n in 9 10 11; do run $v $O/s$n.paf ani 1000000; done
  run $v $D/check/sweep_test2/repro.f48_chrX.raw.paf ani 100000000
  for f in f32 f48; do run $v $W/$f.paf ani 10000000; done
done
log done
