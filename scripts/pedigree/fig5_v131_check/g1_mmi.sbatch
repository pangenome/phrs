#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 6:00:00
#SBATCH --exclude=octopus05
# Shared diploid indexes, one person
# usage: sbatch g1_mmi.sbatch --washu D --sample SM --asm NORM_DIR --out DIR
# out: DIR/SM.hifi.mmi (map-hifi), DIR/SM.ont.mmi (map-ont), minimap2 2.31 -I 8G, from NORM_DIR/SM#1 + SM#2
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sample) SM=$2;; --asm) NA=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D SM NA O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; MM=$D/scripts_tmp/minimap2-2.31/minimap2; mkdir -p $O /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/mmi.XXXX)
zcat $NA/$SM#1.fa.gz $NA/$SM#2.fa.gz > $W/dip.fa
log "$SM from $NA: $(grep -c '>' $W/dip.fa) contigs"
for s in hifi ont; do
  [ -s $O/$SM.$s.mmi ] && { log "exists $O/$SM.$s.mmi, skip"; continue; }
  $MM -t ${SLURM_CPUS_PER_TASK:-4} -I 8G -x map-$s -d $W/$SM.$s.mmi $W/dip.fa 2> $W/$s.log; tail -2 $W/$s.log
  cp $W/$SM.$s.mmi $O/$SM.$s.mmi.part && mv $O/$SM.$s.mmi.part $O/$SM.$s.mmi
  log "$s index $(stat -c %s $O/$SM.$s.mmi) bytes"
done
rm -rf "$W"
log done
