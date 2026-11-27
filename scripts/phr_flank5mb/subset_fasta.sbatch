#!/bin/bash
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus02,octopus05,octopus06,tux06,tux07,tux08,tux09
# Subset a bgzipped, faidx-indexed FASTA to the sequences in NAMES (in that order), on node /scratch, then bgzip -l 9,
# faidx and copy to OUT with a size check. Every name must be extracted with its full length.
# usage: sbatch subset_fasta.sbatch --fasta FA_GZ --names FILE --out FA_GZ
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --fasta) FA=$2;; --names) NM=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in FA NM O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
[ -e $O ] && { echo "ERROR $O exists"; exit 1; }
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/subfa.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
$S faidx $FA -r $NM | $BGZ -@ 8 -l 9 > $W/out.fa.gz
$S faidx $W/out.fa.gz
cut -f1 $W/out.fa.gz.fai | diff - $NM > /dev/null || { echo "ERROR names or order differ from $NM"; exit 1; }
awk 'NR==FNR {L[$1]=$2; next} L[$1] != $2 {bad++} END {exit bad > 0}' $FA.fai $W/out.fa.gz.fai || { echo "ERROR lengths differ from $FA.fai"; exit 1; }
for e in "" .fai .gzi; do cp $W/out.fa.gz$e $O$e.part; [ $(stat -c %s $W/out.fa.gz$e) = $(stat -c %s $O$e.part) ]; mv $O$e.part $O$e; done
echo "[$(date -Is)] $O: $(wc -l < $O.fai) sequences, $(awk '{s+=$2} END {print s}' $O.fai) bp, $(stat -c %s $O) bytes"
