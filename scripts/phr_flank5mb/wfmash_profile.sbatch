#!/bin/bash
#SBATCH -c 48
#SBATCH --mem=64G
#SBATCH -t 2:00:00
#SBATCH --exclude=octopus05,octopus06,octopus11,tux06,tux07,tux08,tux09
# Where does an all-vs-one wfmash run spend its time? The paper call (wfmash -p 95 --target-prefix T, all flanks as
# queries) timed with and without base-level alignment (-m: approximate mappings only), same node, same input copy.
# usage: sbatch wfmash_profile.sbatch --flanks FA_GZ --target NAME --out DIR
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --flanks) FA=$2;; --target) T=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in FA T O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
WF=/moosefs/guarracino/pggb_wfmash023/wfmash/build/bin/wfmash; N=${SLURM_CPUS_PER_TASK:-48}
mkdir -p $O /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/wfprof.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
for e in "" .fai .gzi; do cp $FA$e $W/flanks.fa.gz$e; done
log "host $(hostname), cpus $N, target $T"
TIMEFORMAT="%R s real, %U s user, %S s sys"
for mode in map full; do
  opt=$([ $mode = map ] && echo "-m" || echo "")
  { time $WF -p 95 --target-prefix "$T" -t $N $opt $W/flanks.fa.gz --quiet > $W/$mode.paf ; } 2> $W/$mode.time
  log "$mode: $(cat $W/$mode.time); $(wc -l < $W/$mode.paf) records"
done
cp $W/map.time $W/full.time $O/
gzip -c $W/full.paf > $O/full.paf.gz; gzip -c $W/map.paf > $O/map.paf.gz
log "done: $O"
