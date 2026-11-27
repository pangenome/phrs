#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 6:00:00
#SBATCH --exclude=octopus05,octopus06,octopus11,tux06,tux07,tux08,tux09
# arm_leiden_seeds.R on the control (full) and clean arm distance matrices, reference = the paper's arm-level Leiden.
# usage: sbatch arm_leiden_seeds.sbatch --dir OUT_DIR --ref REF.tsv --seeds N --scripts DIR
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --dir) O=$2;; --ref) REF=$2;; --seeds) N=$2;; --scripts) SD=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in O REF N SD; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
P=/gnu/store/0hfxwjm2jip5i3ydmyy60w44cy466283-profile; export R_LIBS_SITE=$P/site-library; R=$P/bin/Rscript
for set in control clean; do $R $SD/arm_leiden_seeds.R $O/$set/arm.dist_matrix.tsv $REF $N $O/$set/arm-leiden; done
