#!/bin/bash
#SBATCH -c 8
#SBATCH -t 1-00:00:00
# As corrected_stats.sbatch, but the expected artefact is subtracted at the resolution of each statistic (50 and 20 kb)
# instead of 5 kb: at 5 kb the expected artefact is spread in fractions of a contact over pixels that are often empty
# in the real matrix, so about half of it could not be removed (HG002 Hi-C: 142,594 of 290,017). The scale stays the
# 5 kb one of calibrate_artefact.py. Per resolution R: subtract_artefact.py --resolution R --scale-resolution 5000,
# cooler zoomify -r R --balance --field count:dtype=float64 (the same balancing as the paper's zoomify at R), then the
# paper's statistics at R into OUT/stats_R. --factor 0 removes nothing (control: must reproduce the paper's tables).
# usage: sbatch corrected_res_stats.sbatch --sample S --domains TSV --sim MCOOL --real MCOOL --fai FAI --out DIR [--factor F]
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
F=1
while [ $# -gt 0 ]; do case $1 in --sample) SM=$2;; --domains) DOM=$2;; --sim) SIM=$2;; --real) REAL=$2;; --fai) FAI=$2;; --out) O=$2;; --factor) F=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in SM DOM SIM REAL FAI O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; E=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin
[ -e $O ] && { echo "ERROR $O exists"; exit 1; }; mkdir -p $O /scratch/$USER/tmp/phrs
W=$(mktemp -d /scratch/$USER/tmp/phrs/corrres.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
log "host $(hostname); sample $SM; sim $SIM; real $REAL; out $O; factor $F"
for r in 50000 20000; do
  $E/python3 $Q/subtract_artefact.py $DOM $SIM $REAL $W/corr.$r.cool --resolution $r --scale-resolution 5000 --factor $F 2>&1 | tee $O/subtract_$r.log | grep -v '^pixels '
  $E/cooler zoomify $W/corr.$r.cool -r $r --balance --field count:dtype=float64 -o $W/corr.$r.mcool
  bash $Q/contact_stats.sbatch --sample $SM --matrix $W/corr.$r.mcool --resolution $r --fai $FAI --out $O/stats_$r
  rm -f $W/corr.$r.cool $W/corr.$r.mcool
done
log "done"
