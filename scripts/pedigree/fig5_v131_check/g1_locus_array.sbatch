#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH --exclude=octopus05
# g1_locus over a task list
# usage: sbatch --array=1-N%20 g1_locus_array.sbatch --washu D --sample SM --version V --tasks TSV --flank BP --out DIR
#   [--asm NORM_DIR] [--mmi-dir DIR] [--chunk K (lines per task, default 1)] [--time-limit S per locus, default 7200]
#   [--rules R] [--bam-dir DIR (sims: reads aligned to --asm, see g1_locus.sbatch)] [--locus-script PATH under D (default scripts_tmp/g1_locus.sbatch)] ; a locus over the limit gets DIR/<id>.<name>/TIMEOUT and the task goes on (listed, not silent)
# TSV lines: id, name, error (on --asm), read error (on norm_V), [extra read regions R1,R2 (g1_locus --extra)];
# task i runs lines (i-1)K+1 .. iK into DIR/<id>.<name>.
# --mmi-dir: shared diploid indexes (g1_mmi.sbatch), copied once per task to node /scratch, removed on any exit
# (scancel included): leftovers filled octopus08 /scratch on 2026-11-18 (2 tasks failed).
set -euxo pipefail
LM=""; trap 'rc=$?; [ -z "$LM" ] || rm -rf "$LM"; rm -rf /scratch/$USER/tmp/phrs/locus.${SLURM_JOB_ID:-x}; echo "exit status $rc at line ${LINENO}"' EXIT
trap 'exit 143' TERM  # scancel: run the EXIT cleanup (index copy, locus dir)
NA=""; MD=""; BD=""; K=1; TL=7200; RU=""; LS=scripts_tmp/g1_locus.sbatch
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sample) SM=$2;; --version) V=$2;; --tasks) TS=$2;; --flank) F=$2;; --out) O=$2;; --asm) NA=$2;; --mmi-dir) MD=$2;; --chunk) K=$2;; --time-limit) TL=$2;; --rules) RU=$2;; --locus-script) LS=$2;; --bam-dir) BD=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D SM V TS F O SLURM_ARRAY_TASK_ID; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
A=$(( (SLURM_ARRAY_TASK_ID - 1) * K + 1 )); Z=$(( SLURM_ARRAY_TASK_ID * K )); N=$(sed -n "${A},${Z}p" $TS | wc -l)
[ $N -gt 0 ] || { echo "ERROR no lines $A-$Z in $TS"; exit 1; }
if [ -n "$MD" ]; then
  LM=/scratch/$USER/tmp/phrs/mmi.$SLURM_JOB_ID; mkdir -p $LM
  for s in hifi ont; do cp $MD/$SM.$s.mmi $LM/; done
  log "indexes copied to $LM: $(du -sb $LM | cut -f1) bytes"
fi
sed -n "${A},${Z}p" $TS | while IFS=$'\t' read -r i n e r x; do
  [ -e $O/$i.$n/TIMEOUT ] && { log "line $i: TIMEOUT marker from an earlier run, skip"; continue; }
  log "task $SLURM_ARRAY_TASK_ID line $i: $n error $e reads $r extra ${x:-none}"
  rc=0; timeout $TL bash $D/$LS --washu $D --sample $SM --version $V --error $e --read-error $r --flank $F --out $O/$i.$n ${NA:+--asm $NA} ${LM:+--mmi-dir $LM} ${RU:+--rules $RU} ${BD:+--bam-dir $BD} ${x:+--extra $x} || rc=$?
  if [ $rc -eq 124 ]; then pkill -u $USER -f "locus.$SLURM_JOB_ID" || [ $? -eq 1 ]; rm -rf /scratch/$USER/tmp/phrs/locus.$SLURM_JOB_ID; mkdir -p $O/$i.$n; echo "timeout after $TL s, $(date -Is), task $SLURM_ARRAY_TASK_ID" > $O/$i.$n/TIMEOUT; log "line $i: TIMEOUT after $TL s, written $O/$i.$n/TIMEOUT"
  elif [ $rc -ne 0 ]; then echo "ERROR line $i exit $rc"; exit $rc; fi
done
log "done $N loci"
