#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 48
#SBATCH -t 2-00:00:00
#SBATCH --exclude=octopus02,octopus05,octopus06,tux06,tux07,tux08,tux09
# All-vs-one wfmash, 5 Mb flanks (as flank5mb_wfmash_v2.sbatch: the paper's call, wfmash v0.23.0-41 -p 95
# --target-prefix <flank> on the whole flank set; task i takes line i of FA_GZ.fai; output DIR/<flank>.paf.gz), timed
# with GNU time -v (CPU and peak memory). The cluster schedules by CPU only (CR_CPU, no memory limits), so a whole node
# (-c 48) per task keeps two large tasks from sharing a node's RAM; one 5 Mb acrocentric task peaks near 200 GB.
# usage: sbatch --array=LIST%N flank5mb_wfmash_v3.sbatch --flanks FA_GZ --out-dir DIR
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --flanks) FA=$2;; --out-dir) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in FA O SLURM_ARRAY_TASK_ID; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
WF=/moosefs/guarracino/pggb_wfmash023/wfmash/build/bin/wfmash; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
TM=/gnu/store/310psr9h9kzzd4lbssxrzjg53zhqjxrq-time-1.9/bin/time
SEQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FA.fai | cut -f1); [ -n "$SEQ" ] || { echo "ERROR no line $SLURM_ARRAY_TASK_ID in $FA.fai"; exit 1; }
OUT=$O/$(echo "$SEQ" | tr ':@#-' '____').paf.gz; mkdir -p $O
[ -s $OUT ] && { log "exists $OUT, skip"; exit 0; }
N=${SLURM_CPUS_PER_TASK:-48}
log "target $SEQ ($(awk -v s="$SEQ" '$1==s {print $2}' $FA.fai) bp), queries $(wc -l < $FA.fai) flanks, host $(hostname), threads $N, $($WF --version 2>&1 | head -1), node free memory $(free -g | awk '/Mem/ {print $7}') GB"
mkdir -p /scratch/$USER/tmp/phrs; T=$(mktemp /scratch/$USER/tmp/phrs/wf5mb.XXXX.paf.gz); trap 'rc=$?; rm -f $T $T.time; echo "exit status $rc at line ${LINENO}"' EXIT
$TM -v -o $T.time $WF -p 95 --target-prefix "$SEQ" -t $N $FA --quiet | $BGZ -l 9 -@ 4 > $T
log "wfmash: $(awk -F': ' '/Elapsed|User time|System time|Maximum resident/ {printf "%s %s; ", $1, $2}' $T.time | sed 's/\t//g')"
cp $T.time $O/$(basename $OUT .paf.gz).time
cp $T $OUT.part && mv $OUT.part $OUT
log "done: $(zcat $OUT | wc -l) PAF lines"
