#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 4
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus02,octopus05,octopus06
# The paper's contact statistics (analyze_hic_communities.py call of PHR_III/run_exclusion_controls.sh, without
# --force-unassigned): arm Leiden communities and arm distance matrix of 2026-03-28, PHR similarity and PHR table,
# window 500 kb, 10,000 permutations, per haplotype, on one mcool; output in a new directory. Python: PY, default the
# paper env rebuilt from its exact package list (condatools/porec/1.0-rebuilt: Python 3.12.12, cooler 0.10.4, numpy
# 1.26.4, scipy 1.17.0, pandas 3.0.0), because porec/1.0 lost compiled modules (sklearn, pysam). With the hicexplorer
# 3.7.4 env the paper's HG002 Hi-C 50 kb outputs were reproduced byte for byte (2882729).
# usage: sbatch contact_stats.sbatch --sample S --matrix MCOOL --resolution R --fai FAI --out DIR
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --sample) SM=$2;; --matrix) MC=$2;; --resolution) RES=$2;; --fai) FAI=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in SM MC RES FAI O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
B=/moosefs/guarracino/HPRCv2; PY=${PY:-/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3}
SCRIPT=$B/scripts/community/analyze_hic_communities.py
SIM_TSV=$B/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.similarity.tsv.gz
PHR_TSV=$B/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv
SEQ_COMMS=$B/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv
ARM_DIST=$B/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv
[ -e $O ] && { echo "ERROR $O exists"; exit 1; }; mkdir -p $O
echo "[$(date -Is)] $(hostname); $($PY -c 'import sys, cooler, numpy, scipy, pandas, sklearn; print(sys.version.split()[0], cooler.__version__, numpy.__version__, scipy.__version__, pandas.__version__, sklearn.__version__)'); script md5 $(md5sum $SCRIPT | cut -c1-32)"
$PY $SCRIPT --organism human --sample $SM --matrix $MC --resolution $RES --fai $FAI --window-size 500000 --per-haplotype \
    --communities-tsv $SEQ_COMMS --arm-dist-matrix $ARM_DIST --phr-tsv $PHR_TSV --phr-similarity $SIM_TSV \
    --output-dir $O --n-permutations 10000
echo "[$(date -Is)] done: $(ls $O | wc -l) files"
