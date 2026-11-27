#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05,octopus06,octopus11,tux06,tux07,tux08,tux09
# Communities without the dropped PHR sequences. From the cached 15,668 x 15,668 Jaccard distance RDS (pggb + odgi
# similarity of the paper) we drop the sequences in DROP and rerun, with the paper's commands, arm-level Leiden
# (detect_communities.R --level arm) and sequence-level Leiden (k-NN scan 10..125, resolution 0.1..3 by 0.1), plus
# arm-level UPGMA (upgma_arms.R, as in plot-similarity-subtelo.R). The same runs on the full RDS are the control: they
# must reproduce the March assignments with this R environment before the new ones are compared.
# usage: sbatch rerun_communities.sbatch --rds RDS --drop DROP --scripts DIR --out DIR
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --rds) RDS=$2;; --drop) DROP=$2;; --scripts) SD=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in RDS DROP SD O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
P=/gnu/store/0hfxwjm2jip5i3ydmyy60w44cy466283-profile; export R_LIBS_SITE=$P/site-library; R=$P/bin/Rscript
DC=/moosefs/guarracino/HPRCv2/scripts/community/detect_communities.R
mkdir -p $O/control $O/clean
log "host $(hostname); $($R --version 2>&1 | head -1); rds $RDS; drop $DROP ($(wc -l < $DROP) names); detect_communities.R md5 $(md5sum < $DC | cut -c1-32)"
$R $SD/subset_dist_rds.R $RDS $DROP $O/clean/dist_matrix.rds
for set in control clean; do
  in=$([ $set = control ] && echo $RDS || echo $O/clean/dist_matrix.rds)
  log "$set: arm-level Leiden"
  $R $DC --cached-rds $in --organism human --level arm --output-prefix $O/$set/arm
  log "$set: arm-level UPGMA"
  $R $SD/upgma_arms.R $O/$set/arm.dist_matrix.tsv $O/$set/arm-upgma.communities.tsv
  log "$set: sequence-level Leiden"
  $R $DC --cached-rds $in --organism human --level sequence --k-nn 50 --k-nn-scan 10,25,50,75,100,125 \
    --resolution-min 0.1 --resolution-max 3.0 --resolution-step 0.1 --output-prefix $O/$set/seq
done
ls -la $O/control $O/clean
log "done: $O"
