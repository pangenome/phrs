#!/bin/bash
#SBATCH -c 8
#SBATCH -t 1-00:00:00
# The paper's contact statistics (contact_stats.sbatch: analyze_hic_communities.py, Mantel and within/between) on a real
# matrix minus the expected mapping artefact (subtract_artefact.py, upper bound), referee R1-M4 / R3-M3. In node
# /scratch: subtract at 5 kb, then cooler zoomify -r 5000,10000,20000,50000,100000 --balance as for the paper's
# matrices (float counts kept with --field count:dtype=float64), then the statistics at 50 and 20 kb into OUT/stats_<res>.
# Only the statistics go to moosefs. --self-test passes the simulation as the real matrix: no contact may be left
# between different ends' domains, and no statistics run.
# --factor 0 removes nothing: a control that the pipeline alone reproduces the paper's statistics.
# usage: sbatch corrected_stats.sbatch --sample S --domains TSV --sim MCOOL --real MCOOL --fai FAI --out DIR [--self-test 1] [--factor F]
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
ST=0; F=1
while [ $# -gt 0 ]; do case $1 in --sample) SM=$2;; --domains) DOM=$2;; --sim) SIM=$2;; --real) REAL=$2;; --fai) FAI=$2;; --out) O=$2;; --self-test) ST=$2;; --factor) F=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in SM DOM SIM REAL FAI O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; E=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin
[ -e $O ] && { echo "ERROR $O exists"; exit 1; }; mkdir -p $O /scratch/$USER/tmp/phrs
W=$(mktemp -d /scratch/$USER/tmp/phrs/corr.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
log "host $(hostname); sample $SM; sim $SIM; real $REAL; out $O; self-test $ST; factor $F; $(df -h /scratch | tail -1)"
$E/python3 $Q/subtract_artefact.py $DOM $SIM $REAL $W/corr.5000.cool --factor $F 2>&1 | tee $O/subtract.log | grep -v '^pixels '
if [ $ST = 1 ]; then
  grep -q "real contacts left between different ends' domains 0.0$" $O/subtract.log && grep -q "not removable (floor at 0) 0.0;" $O/subtract.log
  log "self-test passed"; exit 0
fi
$E/cooler zoomify $W/corr.5000.cool -r 5000,10000,20000,50000,100000 --balance --field count:dtype=float64 -o $W/corr.mcool
log "zoomify done: $(stat -c %s $W/corr.mcool) bytes"
for r in 50000 20000; do bash $Q/contact_stats.sbatch --sample $SM --matrix $W/corr.mcool --resolution $r --fai $FAI --out $O/stats_$r; done
log "done"
