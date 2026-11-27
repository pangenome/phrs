#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 6:00:00
#SBATCH --exclude=octopus05,octopus06,octopus11,tux06,tux07,tux08,tux09
# Sequence-level Leiden at the paper's chosen parameters (k-NN 75, resolution 0.8; no scan), on the full RDS (control)
# and on the RDS without the dropped sequences (clean), to see what dropping them does at fixed parameters.
# usage: sbatch seq_fixed.sbatch --rds RDS --clean-rds RDS --out DIR
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --rds) RDS=$2;; --clean-rds) CR=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in RDS CR O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
P=/gnu/store/0hfxwjm2jip5i3ydmyy60w44cy466283-profile; export R_LIBS_SITE=$P/site-library; R=$P/bin/Rscript
DC=/moosefs/guarracino/HPRCv2/scripts/community/detect_communities.R
mkdir -p $O
for set in control clean; do
  in=$([ $set = control ] && echo $RDS || echo $CR)
  $R $DC --cached-rds $in --organism human --level sequence --k-nn 75 --resolution-min 0.8 --resolution-max 0.8 \
    --resolution-step 0.1 --output-prefix $O/$set.seq-k75-r0.8
done
ls -la $O
