#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05
# g1_locus over a task list
# usage: sbatch --array=1-N g1_locus_array.sbatch --washu D --sample SM --version V --tasks TSV --flank BP --out DIR [--asm NORM_DIR]
# TSV lines: id, name, error (on --asm), read error (on norm_V); task i runs line i into DIR/<id>.<name>
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
NA=""
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sample) SM=$2;; --version) V=$2;; --tasks) TS=$2;; --flank) F=$2;; --out) O=$2;; --asm) NA=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D SM V TS F O SLURM_ARRAY_TASK_ID; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
IFS=$'\t' read -r i n e r < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $TS); [ -n "${r:-}" ] || { echo "ERROR no line $SLURM_ARRAY_TASK_ID in $TS"; exit 1; }
echo "[$(date -Is)] task $SLURM_ARRAY_TASK_ID: $i $n error $e reads $r"
bash $D/scripts_tmp/g1_locus.sbatch --washu $D --sample $SM --version $V --error $e --read-error $r --flank $F --out $O/$i.$n ${NA:+--asm $NA}
