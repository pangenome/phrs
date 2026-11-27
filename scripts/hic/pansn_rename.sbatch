#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH --exclude=octopus02,octopus05,octopus06
# Rename the sequences of a bgzipped FASTA to PanSN as L<n>#0#<old name> (n = order in the file, from 1), so that pggb
# treats every sequence as its own haplotype and wfmash maps every sequence against all others (pggb does not map
# sequences with the same PanSN prefix to each other). Writes OUT (bgzip -l 9, .fai, .gzi) and OUT.names.tsv (new, old);
# checks that names, order and lengths match the input .fai.
# usage: sbatch pansn_rename.sbatch --fasta FA_GZ --out OUT_FA_GZ
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --fasta) FA=$2;; --out) OUT=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in FA OUT; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
[ -e $OUT ] && { echo "ERROR $OUT exists"; exit 1; }
zcat $FA | awk '/^>/ {n++; print ">L" n "#0#" substr($1, 2); next} {print}' | $BGZ -l 9 > $OUT.part
$S faidx $OUT.part
paste <(cut -f1,2 $OUT.part.fai) <(cut -f1,2 $FA.fai) | awk -F'\t' '{split($1, a, "#"); if (a[3] != $3 || $2 != $4) bad++; n++} END {print n " sequences, " bad + 0 " name or length mismatches"; exit (bad > 0 || n == 0)}'
cut -f1 $OUT.part.fai | paste - <(cut -f1 $FA.fai) > $OUT.names.tsv
mv $OUT.part $OUT; mv $OUT.part.fai $OUT.fai; mv $OUT.part.gzi $OUT.gzi
echo "[$(date -Is)] $OUT: $(wc -l < $OUT.fai) sequences, $(awk '{s+=$2} END {print s}' $OUT.fai) bp"
