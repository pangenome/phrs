#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH --exclude=octopus05
# Fig5 KID.m query, one genome name
# usage: sbatch fig5_kidm_inputs.sbatch --washu D --version V --out INPUTS_DIR
# KID.m = hap2 for chr1 3 4 5 7 9 12 14 15 18 21 22 X, hap1 for the rest (kid_phase 2873068-9, = v1.0 trio table).
# All query contigs are named PAN028#1#chrN: with PAN028#1 and PAN028#2 names sweepga 0.1.1 sees 2 genome groups and
# aligns query vs query, dropping the target (2873257). names.tsv maps each query name to its source haplotype.
# target = inputs_V/PAN028mat_vs_PAN027_joint.target.fa (SELF joint), copied. New dir only.
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --version) V=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D V O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; C=PAN028mat_vs_PAN027_joint; S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
H2=" chr1 chr3 chr4 chr5 chr7 chr9 chr12 chr14 chr15 chr18 chr21 chr22 chrX "
[ -s $O/$C.query.fa.fai ] && { log "exists $O/$C.query.fa.fai, skip"; exit 0; }
[ -e $O ] && { echo "ERROR $O exists without query index; use a new dir"; exit 1; }
mkdir -p $O
cut -f1 norm_$V/PAN028#1.fa.gz.fai | while read -r c; do n=${c##*#}; case "$H2" in *" $n "*) h=2;; *) h=1;; esac; printf "PAN028#1#$n\tPAN028#$h#$n\n"; done > $O/names.tsv
while IFS=$'\t' read -r q s; do $S faidx norm_$V/${s%#chr*}.fa.gz "$s" | sed "1s/.*/>$q/"; done < $O/names.tsv > $O/$C.query.fa
$S faidx $O/$C.query.fa; cp inputs_$V/$C.target.fa inputs_$V/$C.target.fa.fai $O/
log "query $(wc -l < $O/$C.query.fa.fai) contigs ($(cut -f2 $O/names.tsv | grep -c '#2#') from hap2), $(awk '{s+=$2} END {print s}' $O/$C.query.fa.fai) bp; prefixes $(cut -f1 $O/$C.query.fa.fai | cut -d'#' -f1,2 | sort -u | tr '\n' ' ')"
log done
