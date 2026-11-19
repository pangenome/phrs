#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 12:00:00
# sweepga main vs overlap fix, octopus
# args: WASHU, OUT_DIR (holds fix_overlap.patch, s9/s10/s11.paf)
# Builds sweepga main and main + fix_overlap.patch in /scratch, runs cargo test --no-fail-fast on both,
# the synthetic cases, the chrX repro and the full PAN027pat v1.1 raw PAFs (f32, f48) with --scoring ani.
# out: OUT_DIR/tests.{main,fix}.txt, summary.tsv, sweepga.{main,fix} binaries
set -uo pipefail
D=$1; O=$2; mkdir -p $O; cd $O; log () { echo "[$(date -Is)] $*"; }
W=/scratch/$USER.sweepgafix.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
export PATH=$HOME/.cargo/bin:/gnu/store/c8zkkq4dwblbvkg2589sllnzpnk2nfav-cmake-4.1.3/bin:$PATH  # compute nodes have no cmake
log "host $(hostname), $(cargo --version), $(cmake --version | head -1)"
git clone -q https://github.com/pangenome/sweepga.git $W/main && git clone -q $W/main $W/fix && (cd $W/fix && git apply $O/fix_overlap.patch) || { log "ERROR clone or patch"; exit 1; }
log "main $(git -C $W/main log -1 --format=%h), fix: $(git -C $W/fix diff --stat | tail -1)"
for v in main fix; do
  (cd $W/$v && cargo build --release 2>&1 | grep -E "^error|Finished") ; [ -x $W/$v/target/release/sweepga ] || { log "ERROR build $v"; exit 1; }
  cp $W/$v/target/release/sweepga $O/sweepga.$v
  (cd $W/$v && cargo test --release --no-fail-fast 2>&1 | grep -E "^test .* \.\.\. |^test result") > $O/tests.$v.txt
  log "tests $v: $(grep -c ' ok$' $O/tests.$v.txt) ok, $(grep -c 'FAILED$' $O/tests.$v.txt) failed"
done
printf "binary\tinput\tscoring\tkept\tlong_kept_or_ge10Mb\n" > $O/summary.tsv
run () { # binary input scoring min_len
  $O/sweepga.$1 --num-mappings 10:10 --scaffold-jump 0 --scoring $3 --temp-dir $W --output-file $W/o.paf $2 > /dev/null 2>&1 || { log "ERROR sweepga $1 $2 $3"; return; }
  printf "%s\t%s\t%s\t%d\t%d\n" $1 $(basename $2) $3 $(wc -l < $W/o.paf) $(awk -F'\t' -v m=$4 '$11>=m' $W/o.paf | wc -l) >> $O/summary.tsv
  log "$(tail -1 $O/summary.tsv)"
}
zcat $D/regen_B_v11/raw_paf/PAN027pat_vs_PAN011_joint.frequency32.paf.gz > $W/f32.paf
zcat $D/fig5_f48/out_v11/raw_paf/PAN027pat_vs_PAN011_joint.frequency48.paf.gz > $W/f48.paf
for v in main fix; do
  for n in 9 10 11; do for sc in ani log-length-ani; do run $v $O/s$n.paf $sc 1000000; done; done
  for sc in ani log-length-ani; do run $v $D/check/sweep_test2/repro.f48_chrX.raw.paf $sc 100000000; done
  for f in f32 f48; do run $v $W/$f.paf ani 10000000; done
done
log "tests that differ (main vs fix):"; diff <(grep -v "test result" $O/tests.main.txt | sort) <(grep -v "test result" $O/tests.fix.txt | sort)
log done
