#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH --exclude=octopus05
# g1_locus over a task list
# usage: sbatch --array=1-N%20 g1_locus_array.sbatch --washu D --sample SM --version V --tasks TSV --flank BP --out DIR
#   [--asm NORM_DIR] [--mmi-dir DIR] [--chunk K (lines per task, default 1)]
# TSV lines: id, name, error (on --asm), read error (on norm_V); task i runs lines (i-1)K+1 .. iK into DIR/<id>.<name>.
# --mmi-dir: shared diploid indexes (g1_mmi.sbatch), copied once per task to node /scratch.
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
NA=""; MD=""; K=1
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sample) SM=$2;; --version) V=$2;; --tasks) TS=$2;; --flank) F=$2;; --out) O=$2;; --asm) NA=$2;; --mmi-dir) MD=$2;; --chunk) K=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D SM V TS F O SLURM_ARRAY_TASK_ID; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
A=$(( (SLURM_ARRAY_TASK_ID - 1) * K + 1 )); Z=$(( SLURM_ARRAY_TASK_ID * K )); N=$(sed -n "${A},${Z}p" $TS | wc -l)
[ $N -gt 0 ] || { echo "ERROR no lines $A-$Z in $TS"; exit 1; }
LM=""
if [ -n "$MD" ]; then
  mkdir -p /scratch/$USER/tmp/phrs; LM=$(mktemp -d /scratch/$USER/tmp/phrs/mmi.XXXX)
  for s in hifi ont; do cp $MD/$SM.$s.mmi $LM/; done
  log "indexes copied to $LM: $(du -sb $LM | cut -f1) bytes"
fi
sed -n "${A},${Z}p" $TS | while IFS=$'\t' read -r i n e r; do
  log "task $SLURM_ARRAY_TASK_ID line $i: $n error $e reads $r"
  bash $D/scripts_tmp/g1_locus.sbatch --washu $D --sample $SM --version $V --error $e --read-error $r --flank $F --out $O/$i.$n ${NA:+--asm $NA} ${LM:+--mmi-dir $LM}
done
[ -z "$LM" ] || rm -rf "$LM"
log "done $N loci"
