#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH --exclude=octopus05,octopus06
# All-vs-one wfmash, 5 Mb flanks
# usage: sbatch --array=1-N%20 flank5mb_wfmash.sbatch --flanks FA_GZ --out-dir DIR
# Same call as the paper (code_PHR-and-3D.md, all-vs-all subtelomeric alignment): wfmash v0.23.0-41 -p 95
# --target-prefix <flank> on the whole flank set; task i takes line i of FA_GZ.fai. Output DIR/<flank>.paf.gz.
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --flanks) FA=$2;; --out-dir) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in FA O SLURM_ARRAY_TASK_ID; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
WF=/moosefs/guarracino/pggb_wfmash023/wfmash/build/bin/wfmash; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
SEQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FA.fai | cut -f1); [ -n "$SEQ" ] || { echo "ERROR no line $SLURM_ARRAY_TASK_ID in $FA.fai"; exit 1; }
OUT=$O/$(echo "$SEQ" | tr ':@#-' '____').paf.gz; mkdir -p $O
[ -s $OUT ] && { log "exists $OUT, skip"; exit 0; }
log "target $SEQ ($(awk -v s="$SEQ" '$1==s {print $2}' $FA.fai) bp), queries $(wc -l < $FA.fai) flanks, host $(hostname), $($WF --version 2>&1 | head -1)"
mkdir -p /scratch/$USER/tmp/phrs; T=$(mktemp /scratch/$USER/tmp/phrs/wf5mb.XXXX.paf.gz); trap 'rc=$?; rm -f $T; echo "exit status $rc at line ${LINENO}"' EXIT
t0=$(date +%s)  # no /usr/bin/time on the nodes; peak memory from sacct MaxRSS
$WF -p 95 --target-prefix "$SEQ" -t ${SLURM_CPUS_PER_TASK:-16} $FA --quiet | $BGZ -l 9 -@ 4 > $T
log "wfmash wall $(( $(date +%s) - t0 )) s"
cp $T $OUT.part && mv $OUT.part $OUT
log "done: $(zcat $OUT | wc -l) PAF lines"
