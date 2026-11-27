#!/bin/bash
# On the node of a flank5mb_phr_call.sbatch job: (1) check the --keep filter on N filtered PAFs (records whose query is in
# REMOVE must be gone, every other record kept); (2) time impg query as the caller runs it (one process per 5 kb window)
# against one process for all windows of the flank (-b BED), and print both outputs' per-window record counts.
# usage: impg_query_cost.sh WORKDIR ORIG_PAF_DIR REMOVE IMPG FLANK N
set -euo pipefail
W=$1; PD=$2; REM=$3; IMPG=$4; FL=$5; N=$6
log () { echo "[$(date -Is)] $*"; }
log "workdir $W, original PAFs $PD, remove list $REM ($(wc -l < $REM) names), impg $IMPG, flank $FL, files checked $N"
bad=0
for f in $(ls $W/paf | awk -v n=$N 'NR <= n'); do
  got=$(zcat $W/paf/$f | wc -l)
  read tot rem <<< $(zcat $PD/$f | awk -F'\t' 'NR == FNR {r[$1]; next} {n++} ($1 in r) {m++} END {print n + 0, m + 0}' $REM -)
  left=$(zcat $W/paf/$f | awk -F'\t' 'NR == FNR {r[$1]; next} ($1 in r) {m++} END {print m + 0}' $REM -)
  ok=$([ "$got" -eq $((tot - rem)) ] && [ "$left" -eq 0 ] && echo OK || echo BAD)
  [ $ok = OK ] || bad=$((bad + 1))
  echo "filter $f original $tot, with a removed query $rem, kept $got, removed queries left $left -> $ok"
done
log "filter check: $bad of $N files BAD"
O=$(dirname $W)/impgcost; mkdir -p $O; B=$O/win.bed
for i in $(seq 0 99); do printf '%s\t%d\t%d\n' "$FL" $((i * 5000)) $(((i + 1) * 5000)); done > $B
Q="--alignment-list $W/paf.list -i $W/idx.impg --min-identity 0.95 --min-output-length 3000 --threads 1 -v 0"
t0=$(date +%s.%N)
for i in $(seq 0 99); do echo "$i $($IMPG query $Q -r "$FL:$((i * 5000))-$(((i + 1) * 5000))" | wc -l)"; done > $O/perwin.txt
t1=$(date +%s.%N)
$IMPG query $Q -b $B -o bedpe > $O/batch.bedpe
t2=$(date +%s.%N)
log "100 windows, one impg process each: $(echo "$t1 - $t0" | bc) s; one impg process for all 100 (-b): $(echo "$t2 - $t1" | bc) s"
echo "per-window records (one process each), first 10: $(head -10 $O/perwin.txt | tr '\n' ' ')"
echo "batch output lines $(wc -l < $O/batch.bedpe); first 3:"; head -3 $O/batch.bedpe | cut -c1-300
log "outputs in $O"
