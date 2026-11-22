#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH --exclude=octopus05,octopus06
# Paralog partner ends of one person (end_partners.py)
# usage: sbatch end_partners.sbatch --washu D --sample SM --version V --out TSV
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sample) SM=$2;; --version) V=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D SM V O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; [ -s $O ] && { log "exists $O, skip"; exit 0; }
mkdir -p /scratch/$USER/tmp/phrs $(dirname $O); W=$(mktemp -d /scratch/$USER/tmp/phrs/endpart.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
log "$SM $V norm_$V, host $(hostname)"
python3 scripts_tmp/end_partners.py norm_$V $SM $D/scripts_tmp/minimap2-2.31/minimap2 ${SLURM_CPUS_PER_TASK:-4} $W $O
log "rows $(($(wc -l < $O) - 1)) in $O"
