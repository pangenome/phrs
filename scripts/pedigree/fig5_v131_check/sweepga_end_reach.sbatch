#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 6
#SBATCH --mem=16G
#SBATCH -t 6:00:00
# End reach of the homolog alignment (per parent copy) in Fig5 sweepGA PAFs, one job, PAFs in parallel
# usage: sbatch sweepga_end_reach.sbatch --washu D --pafs LABEL=PAF[,LABEL=PAF...] --out TSV
# (2026-11-21 first run: the 3 original v1.1 -f32 raw PAFs, --raw-dir; now any list, e.g. the 10:10 filtered PAFs)
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --pafs) PS=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D PS O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; [ -s $O ] && { log "exists $O, skip"; exit 0; }
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/reach.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
for x in ${PS//,/ }; do
  lab=${x%%=*}; f=${x#*=}; [ -s $f ] || { echo "ERROR missing $f"; exit 1; }
  ( zcat $f | python3 scripts_tmp/sweepga_end_reach.py $lab > $W/$lab.tsv 2> $W/$lab.log ) &
done
wait
cat $W/*.log; F=$(ls $W/*.tsv | head -1); head -1 $F > $O; for f in $W/*.tsv; do tail -n +2 $f; done >> $O
log "rows $(($(wc -l < $O) - 1)) in $O"
