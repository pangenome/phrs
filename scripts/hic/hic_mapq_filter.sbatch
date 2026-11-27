#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH -t 1-00:00:00
#SBATCH --exclude=octopus02,octopus05,octopus06
# Hi-C contact matrices restricted by mapping quality (R1-M4 / R3-M3): from a HiC-Pro allValidPairs (columns 11 and 12
# = MAPQ of the two reads), keep the pairs with MAPQ >= T on both reads, for each T in THRESHOLDS, in one read of the
# file; then build the mcool exactly as for the paper's matrices (awk $2,$3,$5,$6 | cooler cload pairs at 5 kb on
# chrom.sizes, cooler zoomify -r 5000,10000,20000,50000,100000 --balance). Bowtie2 in this HiC-Pro setup gives
# multi-mapping reads MAPQ 1 (simulation: 37.7% of PHR reads with MAPQ 1 land on another arm, 0.38% with MAPQ 2-29),
# so T = 2 is the bowtie2 equivalent of minimap2 MAPQ >= 1. Counts of kept pairs per threshold are logged.
# usage: sbatch hic_mapq_filter.sbatch --pairs ALLVALIDPAIRS --sizes CHROM_SIZES --thresholds "1 2" --out PREFIX
#        (writes PREFIX.mapq<T>.mcool)
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --pairs) AVP=$2;; --sizes) CS=$2;; --thresholds) TH=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in AVP CS TH O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
C=/moosefs/guarracino/condatools/hicexplorer/3.7.4/bin/cooler
mkdir -p $(dirname $O) /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/mapq.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
log "host $(hostname); pairs $AVP ($(stat -c %s $AVP) bytes); thresholds $TH; $($C --version)"
awk -F'\t' -v th="$TH" -v w=$W 'BEGIN {OFS = "\t"; n = split(th, t, " ")}
  {tot++; for (i = 1; i <= n; i++) if ($11 >= t[i] && $12 >= t[i]) {k[i]++; print $2, $3, $5, $6 > (w "/mapq" t[i] ".pairs")}}
  END {printf "total %d", tot > "/dev/stderr"; for (i = 1; i <= n; i++) printf "; MAPQ >= %d on both reads %d (%.4f)", t[i], k[i], k[i] / tot > "/dev/stderr"; print "" > "/dev/stderr"}' $AVP
for t in $TH; do
  $C cload pairs -c1 1 -p1 2 -c2 3 -p2 4 $CS:5000 $W/mapq$t.pairs $W/mapq$t.5000.cool
  $C zoomify $W/mapq$t.5000.cool -r 5000,10000,20000,50000,100000 --balance -o $W/mapq$t.mcool
  cp $W/mapq$t.mcool $O.mapq$t.mcool.part && mv $O.mapq$t.mcool.part $O.mapq$t.mcool
  log "MAPQ >= $t: $O.mapq$t.mcool ($(stat -c %s $O.mapq$t.mcool) bytes)"
done
log "done"
