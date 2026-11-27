#!/bin/bash
#SBATCH -c 4
#SBATCH -t 12:00:00
# The paper's contact statistics on a matrix rebalanced with other cooler balance filters (R1-M4 / R3-M3: with unique
# reads only, the default filters mask most PHR bins as too sparse). For each resolution R (50 and 20 kb): copy
# MCOOL::resolutions/R into a new mcool in node /scratch, rebalance it in place with `cooler balance --force BALANCE_ARGS`
# (the cooler that built the MAPQ matrices, hicexplorer 3.7.4 env), log how many bins and how many PHR windows of the
# null-simulation domains keep a weight, then run contact_stats.sbatch at R into OUT/stats_R. With the default filters
# (BALANCE_ARGS "") the statistics must equal the ones of the original matrix (control).
# usage: sbatch rebalance_stats.sbatch --sample S --mcool MCOOL --fai FAI --domains TSV --balance-args "ARGS" --out DIR
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
BA=""
while [ $# -gt 0 ]; do case $1 in --sample) SM=$2;; --mcool) MC=$2;; --fai) FAI=$2;; --domains) DOM=$2;; --balance-args) BA=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in SM MC FAI DOM O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; C=/moosefs/guarracino/condatools/hicexplorer/3.7.4/bin/cooler; PY=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3
[ -e $O ] && { echo "ERROR $O exists"; exit 1; }; mkdir -p $O /scratch/$USER/tmp/phrs
W=$(mktemp -d /scratch/$USER/tmp/phrs/rebal.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
log "host $(hostname); sample $SM; mcool $MC; balance args '$BA'; $($C --version)"
for r in 50000 20000; do
  $C cp $MC::resolutions/$r $W/m.$r.mcool::resolutions/$r
  $C balance --force $BA $W/m.$r.mcool::resolutions/$r
  $PY - $W/m.$r.mcool $MC $r $DOM <<'PY' | tee -a $O/weights.log
import sys, cooler, numpy as np, pandas as pd
new, old, r, dom = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
D = pd.read_csv(dom, sep="\t")
for nm, p in (("original", old), ("rebalanced", new)):
    c = cooler.Cooler(f"{p}::resolutions/{r}"); w = c.bins()["weight"][:].values
    k = sum(np.isfinite(c.bins().fetch(f"{d['chrom']}:{d['phr_start']}-{d['phr_end']}")["weight"].values).any() for _, d in D.iterrows())
    print(f"{r} bp {nm}: bins with a weight {int(np.isfinite(w).sum())} of {len(w)}; PHR windows with at least one weighted bin {k} of {len(D)}")
PY
  bash $Q/contact_stats.sbatch --sample $SM --matrix $W/m.$r.mcool --resolution $r --fai $FAI --out $O/stats_$r
  rm -f $W/m.$r.mcool
done
log "done"
