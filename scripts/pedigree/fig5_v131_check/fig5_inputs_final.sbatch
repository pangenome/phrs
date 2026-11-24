#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH --exclude=octopus05,octopus06
# Fig5 inputs of the 3 comparisons from the final corrected assemblies (one folder with SM#{1,2}.fa.gz per person)
# usage: sbatch fig5_inputs_final.sbatch --washu D --mom DIR --dad DIR --self DIR --kid DIR --ref-inputs DIR --ref-kid DIR --out DIR
# Same names as the original inputs (build_inputs.sh, fig5_kidm_inputs.sbatch): query = PAN027#1 (mat) or PAN027#2 (pat)
# contigs as they are; KID.m query = hap2 for chr1 3 4 5 7 9 12 14 15 18 21 22 X, hap1 for the rest, all named
# PAN028#1#chrN (sweepga sees one query genome); target = parent h1 + h2 renamed PAN0xx#H#chrN -> PAN0xx#joint#hH_chrN.
# Check: same contig names in the same order as --ref-inputs (mat, pat) and --ref-kid (kid); contigs whose length changed
# are listed in the log. New dir only.
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --mom) M=$2;; --dad) F=$2;; --self) P=$2;; --kid) K=$2;; --ref-inputs) RI=$2;; --ref-kid) RK=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D M F P K RI RK O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
[ -e $O ] && { echo "ERROR $O exists; use a new dir"; exit 1; }
mkdir -p $O
log "inputs: MOM $M, DAD $F, SELF $P, KID $K; reference names $RI, $RK"
joint () { for h in 1 2; do zcat $1/$2#$h.fa.gz | sed "s|^>$2#$h#|>$2#joint#h${h}_|"; done; }
check () {  # $1 new fai, $2 reference fai, $3 label
  [ "$(cut -f1 $1)" = "$(cut -f1 $2)" ] || { echo "ERROR $3: contig names or order differ from $2"; exit 1; }
  log "$3: $(wc -l < $1) contigs, same names as $2; length changed: $(paste $1 $2 | awk -F'\t' '$2 != $(NF - 3) {printf "%s %+d ", $1, $2 - $(NF - 3); n++} END {print n + 0 " contigs"}')"
}
for x in "PAN027mat_vs_PAN010_joint 1 $M PAN010" "PAN027pat_vs_PAN011_joint 2 $F PAN011"; do set -- $x
  zcat $P/PAN027#$2.fa.gz > $O/$1.query.fa; joint $3 $4 > $O/$1.target.fa
  for t in query target; do $S faidx $O/$1.$t.fa; check $O/$1.$t.fa.fai $RI/$1.$t.fa.fai "$1 $t"; done
done
C=PAN028mat_vs_PAN027_joint; H2=" chr1 chr3 chr4 chr5 chr7 chr9 chr12 chr14 chr15 chr18 chr21 chr22 chrX "
cut -f1 $K/PAN028#1.fa.gz.fai | while read -r c; do n=${c##*#}; case "$H2" in *" $n "*) h=2;; *) h=1;; esac; printf "PAN028#1#$n\tPAN028#$h#$n\n"; done > $O/names.tsv
while IFS=$'\t' read -r q s; do $S faidx $K/${s%#chr*}.fa.gz "$s" | sed "1s/.*/>$q/"; done < $O/names.tsv > $O/$C.query.fa
joint $P PAN027 > $O/$C.target.fa
for t in query target; do $S faidx $O/$C.$t.fa; check $O/$C.$t.fa.fai $RK/$C.$t.fa.fai "$C $t"; done
log "KID.m query: $(cut -f2 $O/names.tsv | grep -c '#2#') contigs from hap2"
log done
