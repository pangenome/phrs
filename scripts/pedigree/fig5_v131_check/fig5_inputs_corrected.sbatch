#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH --exclude=octopus05
# Fig5 pat inputs, corrected DAD
# usage: sbatch fig5_inputs_corrected.sbatch --washu D --version V --norm-dir CORRECTED_NORM_DIR --out INPUTS_DIR
# query = inputs_V/PAN027pat_vs_PAN011_joint.query.fa (SELF paternal, unchanged); target = corrected PAN011 h1 + h2,
# names rewritten as build_inputs.sh (PAN011#H#chrN -> PAN011#joint#hH_chrN). New dir only.
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --version) V=$2;; --norm-dir) N=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D V N O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; C=PAN027pat_vs_PAN011_joint; S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
[ -s $O/$C.target.fa.fai ] && { log "exists $O/$C.target.fa.fai, skip"; exit 0; }
[ -e $O ] && { echo "ERROR $O exists without target index; use a new dir"; exit 1; }
mkdir -p $O; cp inputs_$V/$C.query.fa inputs_$V/$C.query.fa.fai $O/
for h in 1 2; do zcat $N/PAN011#$h.fa.gz | sed "s|^>PAN011#$h#|>PAN011#joint#h${h}_|"; done > $O/$C.target.fa
$S faidx $O/$C.target.fa
log "query $(wc -l < $O/$C.query.fa.fai) contigs, target $(wc -l < $O/$C.target.fa.fai) contigs; target vs inputs_$V: $(paste <(cut -f1,2 inputs_$V/$C.target.fa.fai) <(cut -f1,2 $O/$C.target.fa.fai) | awk '$1 != $3 {n++} $2 != $4 {m++} END {print n+0" names differ, "m+0" lengths differ"}')"
log done
