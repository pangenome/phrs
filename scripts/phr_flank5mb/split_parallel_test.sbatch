#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH -t 6:00:00
#SBATCH --exclude=octopus02,octopus05,octopus06,tux06,tux07,tux08,tux09
# Wall time of one all-vs-one 5 Mb target (the paper's call, wfmash -p 95 --target-prefix T on the whole flank FASTA)
# for configurations K:P:T on the same node: queries split round-robin into K -A lists, P parts at a time, T threads
# each (P x T = 16 cores; peak memory per task is about P/K of the single run). The union of the parts is compared
# record by record with REF (compare_paf.py, nothing excluded). Prints per configuration: wall, CPU (user+sys summed
# over parts), peak RSS of the largest part, records, identical to REF.
# --stage L (optional): stage FA_GZ and wfmash in the node-local dir L exactly as flank5mb_wfmash_split.sbatch does
# (same lock, md5 and marker, so later production tasks on the node reuse the copy).
# usage: sbatch split_parallel_test.sbatch --flanks FA_GZ --target NAME --ref PAF_GZ --configs "4:1:16 8:2:8" --out DIR [--stage L]
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --flanks) FA=$2;; --target) SEQ=$2;; --ref) REF=$2;; --configs) CF=$2;; --out) O=$2;; --stage) L=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in FA SEQ REF CF O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
WF=/moosefs/guarracino/pggb_wfmash023/wfmash/build/bin/wfmash; TM=/gnu/store/310psr9h9kzzd4lbssxrzjg53zhqjxrq-time-1.9/bin/time
CMP=/moosefs/guarracino/HPRCv2/PHR_III/flank5mb/compare_paf.py
if [ -n "${L:-}" ]; then
  mkdir -p $L; exec 9> $L/.lock; flock 9
  if [ ! -s $L/ok ]; then
    log "staging $FA and wfmash to $L on $(hostname)"
    ref=$(cut -d' ' -f1 $FA.md5); [ -n "$ref" ] || { echo "ERROR no md5 in $FA.md5"; exit 1; }
    cp $FA $L/tmp.fa.gz; got=$(md5sum $L/tmp.fa.gz | cut -d' ' -f1)
    [ "$got" = "$ref" ] || { echo "ERROR staged copy md5 $got, expected $ref"; exit 1; }
    for x in fai gzi; do cp $FA.$x $L/tmp.fa.gz.$x; cmp $FA.$x $L/tmp.fa.gz.$x; done
    cp $WF $L/wfmash.tmp; cmp $WF $L/wfmash.tmp
    for x in fai gzi; do mv $L/tmp.fa.gz.$x $L/$(basename $FA).$x; done; mv $L/tmp.fa.gz $L/$(basename $FA); mv $L/wfmash.tmp $L/wfmash
    echo "$ref $(basename $FA)" > $L/ok; log "staged on $(hostname): md5 $got"
  fi
  exec 9>&-
  [ "$(cut -d' ' -f1 $L/ok)" = "$(cut -d' ' -f1 $FA.md5)" ] || { echo "ERROR $L/ok does not match $FA.md5"; exit 1; }
  FA=$L/$(basename $FA); WF=$L/wfmash
fi
grep -qxF "$SEQ" <(cut -f1 $FA.fai) || { echo "ERROR $SEQ not in $FA.fai"; exit 1; }
mkdir -p $O /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/sptest.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
log "target $SEQ; configs $CF; host $(hostname); $(lscpu | grep '^Model name' | tr -s ' '); node free memory $(free -g | awk '/Mem/ {print $7}') GB"
export WF TM FA SEQ
for c in $CF; do
  IFS=: read K P T <<< "$c"; D=$W/k${K}p${P}t${T}; mkdir -p $D
  cut -f1 $FA.fai | awk -v k=$K -v w=$D '{print > (w "/q" (NR - 1) % k ".txt")}'
  t0=$(date +%s.%N)
  seq 0 $((K - 1)) | xargs -P $P -I{} bash -c 'set -euo pipefail; $TM -v -o '$D'/q{}.time $WF -p 95 --target-prefix "$SEQ" -A '$D'/q{}.txt -t '$T' $FA --quiet > '$D'/q{}.paf'
  wall=$(awk -v a=$t0 -v b=$(date +%s.%N) 'BEGIN {printf "%.0f", b - a}')
  cat $(for i in $(seq 0 $((K - 1))); do echo $D/q$i.paf; done) > $D/all.paf
  cpu=$(cat $D/q*.time | awk -F': ' '/User time|System time/ {s += $2} END {printf "%.0f", s}')
  rss=$(cat $D/q*.time | awk -F': ' '/Maximum resident/ {if ($2 > m) m = $2} END {printf "%.1f", m / 1048576}')
  cmp=$(python3 $CMP /dev/null $REF $D/all.paf | head -1 | cut -f2)
  log "config K $K P $P T $T: wall $wall s; CPU $cpu s; max part RSS $rss GB; $(wc -l < $D/all.paf) records; $cmp"
  echo -e "$SEQ\t$K\t$P\t$T\t$wall\t$cpu\t$rss\t$cmp" >> $O/results.tsv
done
log "done"
