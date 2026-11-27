#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH --exclude=octopus05,octopus06,tux06,tux07,tux08,tux09
# Pilot for the R1-M8 PHR recall on 5 Mb flanks: does impg 0.4.1 `query` project a region only when its sequence
# is the target of an alignment, or also when it is the query? Exit code for an unknown sequence? Wall time per query?
# usage: sbatch impg_pilot.sbatch PAF_GZ IMPG_BIN [IMPG_BIN ...] (one all-vs-one PAF; each binary is tested in turn)
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
P=$1; shift
log () { echo "[$(date -Is)] $*"; }
mkdir -p /scratch/guarracino/tmp/phrs; W=$(mktemp -d /scratch/guarracino/tmp/phrs/impgpilot.XXXX)
trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
cp $P $W/; ls $W/*.paf.gz > $W/list
for IMPG in "$@"; do
log "binary $IMPG: $($IMPG --version)"; rm -f $W/idx.impg
t0=$(date +%s); $IMPG index --alignment-list $W/list -i $W/idx.impg -t 8; log "index wall $(( $(date +%s) - t0 )) s"
read T Q N NQ < <(zcat $P | awk -F'\t' 'NR == 1 {t = $6} $1 != t && q == "" {q = $1} !($1 in s) {s[$1] = 1; k++} END {print t, q, NR, k}')
log "target $T, a query $Q; PAF lines $N, distinct queries $NQ"
for r in "$T:0-5000" "$T:2000000-2005000" "$Q:0-5000" "NOSUCHSEQ:0-5000"; do
  t0=$(date +%s%N)
  rc=0; $IMPG query --alignment-list $W/list -i $W/idx.impg -r "$r" --min-identity 0.95 --min-output-length 3000 --threads 1 -v 0 > $W/q.bed 2> $W/q.err || rc=$?
  log "region $r: exit $rc, lines $(wc -l < $W/q.bed), wall $(( ($(date +%s%N) - t0) / 1000000 )) ms, stderr: $(head -c 300 $W/q.err | tr '\n' ' ')"
  head -3 $W/q.bed | cut -f1-6
done
done
