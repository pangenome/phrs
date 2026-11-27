#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH --exclude=octopus02,octopus05,octopus06,tux06,tux07,tux08,tux09
# All-vs-one wfmash on the 5 Mb flanks with the paper's call (wfmash v0.23.0-41 -p 95 --target-prefix <flank> on the
# whole flank FASTA), with the queries split into K lists (-A, round-robin over the .fai names), P parts at a time with
# 16/P threads each.
# The target index is the same in every part and each query is mapped on its own, so the union of the parts equals the
# single run (query_split_test.sbatch: 55,925 of 55,925 records identical for NA18570#1 18p, K 4). Peak memory scales
# with the queries (about 215 GB for all 18,827 flanks, 56 GB per part with K 4), and the cluster schedules by CPU only
# (CR_CPU, --mem is not enforced), so K 4 at 16 threads lets 3 tasks share a 48-core node. One part runs on few cores
# for a long stretch, so P parts at a time fill the 16 cores: K 16, P 4 (4 threads, 13 GB each) took 366 / 380 s per
# target vs 715 / 729 s for K 4, P 1, same CPU and memory per task, identical records (split_parallel_test.sbatch,
# 2881250 / 2881251). Task i takes line i of
# TARGETS (flank names) and writes DIR/<flank>.paf.gz (parts concatenated in order) and DIR/<flank>.time (GNU time -v
# of each part, in order).
# --stage L (optional): the first task on a node copies FA_GZ, .fai, .gzi and wfmash to the node-local dir L (under a
# lock) and checks the copy against FA_GZ.md5 (md5sum of FA_GZ) and cmp for the rest; later tasks on that node read
# from L, so moosefs serves one 25 GB read per node instead of one per task.
# usage: sbatch --array=1-N%M flank5mb_wfmash_split.sbatch --flanks FA_GZ --targets TXT --parts K --parallel P --out-dir DIR [--stage L]
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --flanks) FA=$2;; --targets) TG=$2;; --parts) K=$2;; --parallel) P=$2;; --out-dir) O=$2;; --stage) L=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in FA TG K P O SLURM_ARRAY_TASK_ID; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
WF=/moosefs/guarracino/pggb_wfmash023/wfmash/build/bin/wfmash; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
TM=/gnu/store/310psr9h9kzzd4lbssxrzjg53zhqjxrq-time-1.9/bin/time
SEQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $TG); [ -n "$SEQ" ] || { echo "ERROR no line $SLURM_ARRAY_TASK_ID in $TG"; exit 1; }
LEN=$(awk -v s="$SEQ" '$1==s {print $2}' $FA.fai); [ -n "$LEN" ] || { echo "ERROR $SEQ not in $FA.fai"; exit 1; }
OUT=$O/$(echo "$SEQ" | tr ':@#-' '____').paf.gz; mkdir -p $O
[ -s $OUT ] && { log "exists $OUT, skip"; exit 0; }
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
N=${SLURM_CPUS_PER_TASK:-16}; T=$((N / P)); [ $((T * P)) -eq $N ] || { echo "ERROR $N cores not divisible by P $P"; exit 2; }
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/wfsplit.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
cut -f1 $FA.fai | awk -v k=$K -v w=$W '{print > (w "/q" (NR - 1) % k ".txt")}'
log "target $SEQ ($LEN bp), queries $(wc -l < $FA.fai) flanks in $K lists ($(wc -l $W/q*.txt | head -$K | awk '{printf "%s ", $1}')), $P at a time with $T threads, host $(hostname), $($WF --version 2>&1 | head -1), node free memory $(free -g | awk '/Mem/ {print $7}') GB"
export WF TM FA SEQ W T
seq 0 $((K - 1)) | xargs -P $P -I{} bash -c 'set -euo pipefail; $TM -v -o $W/q{}.time $WF -p 95 --target-prefix "$SEQ" -A $W/q{}.txt -t $T $FA --quiet > $W/q{}.paf'
for i in $(seq 0 $((K - 1))); do
  log "part $i: $(wc -l < $W/q$i.paf) records; $(awk -F': ' '/Elapsed|User time|System time|Maximum resident/ {printf "%s %s; ", $1, $2}' $W/q$i.time | sed 's/\t//g')"
done
cat $(for i in $(seq 0 $((K - 1))); do echo $W/q$i.paf; done) | $BGZ -l 9 -@ 4 > $W/out.paf.gz
cat $(for i in $(seq 0 $((K - 1))); do echo $W/q$i.time; done) > $W/all.time
cp $W/all.time $O/$(basename $OUT .paf.gz).time
cp $W/out.paf.gz $OUT.part
[ $(stat -c %s $W/out.paf.gz) -eq $(stat -c %s $OUT.part) ] || { echo "ERROR copy size differs for $OUT.part"; exit 1; }
mv $OUT.part $OUT
log "done: $(cat $W/q*.paf | wc -l) records in $OUT ($(stat -c %s $OUT) bytes)"
