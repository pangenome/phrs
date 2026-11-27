#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 2
#SBATCH -t 2:00:00
#SBATCH --exclude=octopus02,octopus05,octopus06
# FASTA of a list of regions (samtools region syntax, 1-based closed) from an indexed FASTA, bgzip -l 9 and faidx;
# checks that every region is present once with its expected length.
# usage: sbatch extract_loci_fasta.sbatch --fasta FA_GZ --regions TXT --out OUT_FA_GZ
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --fasta) FA=$2;; --regions) R=$2;; --out) OUT=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in FA R OUT; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
[ -e $OUT ] && { echo "ERROR $OUT exists"; exit 1; }
echo "[$(date -Is)] $(wc -l < $R) regions from $FA"
$S faidx -r $R $FA | $BGZ -l 9 > $OUT.part
$S faidx $OUT.part
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/loci.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
awk -F'[:-]' '{print $0"\t"$3-$2+1}' $R | sort > $W/expected
cut -f1,2 $OUT.part.fai | sort > $W/got
cmp $W/expected $W/got || { echo "ERROR regions or lengths differ"; exit 1; }
mv $OUT.part $OUT; mv $OUT.part.fai $OUT.fai; [ -e $OUT.part.gzi ] && mv $OUT.part.gzi $OUT.gzi
echo "[$(date -Is)] $OUT: $(wc -l < $OUT.fai) sequences, $(awk '{s+=$2} END {print s}' $OUT.fai) bp"
