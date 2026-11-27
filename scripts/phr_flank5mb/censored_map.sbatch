#!/bin/bash
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05,octopus06,tux06,tux07,tux08,tux09
# Where does the sequence of a censored PHR end come from? Extract 5 Mb from the telomere-trimmed end of each contig,
# map it to CHM13 with wfmash -p 95 (same identity as the PHR calls), and summarise per 100 kb from the telomere the CHM13
# chromosomes it aligns to (censored_map_summary.py).
# usage: sbatch censored_map.sbatch --regions TSV (contig, start, end 1-based, arm) --fasta-dir DIR --ref FA --out PREFIX
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --regions) R=$2;; --fasta-dir) FD=$2;; --ref) REF=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in R FD REF O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; WF=/moosefs/guarracino/pggb_wfmash023/wfmash/build/bin/wfmash
log "host $(hostname); regions $(wc -l < $R); ref $REF; $($WF --version 2>&1 | head -1)"
mkdir -p /scratch/guarracino/tmp/phrs $(dirname $O); W=$(mktemp -d /scratch/guarracino/tmp/phrs/cmap.XXXX)
trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
while IFS=$'\t' read -r ctg a b arm; do
  f=$(grep -l -m1 -F "$ctg	" $FD/*.fa.gz.fai | head -1); [ -n "$f" ] || { echo "ERROR $ctg in no fai"; exit 1; }
  $S faidx ${f%.fai} "$ctg:$a-$b" | sed "1s/.*/>$ctg:$a-${b}_$arm/" >> $W/q.fa
  log "$ctg:$a-$b ($arm) from ${f%.fai}"
done < $R
$S faidx $W/q.fa; cut -f1,2 $W/q.fa.fai
t0=$(date +%s); $WF -p 95 -t ${SLURM_CPUS_PER_TASK:-32} $REF $W/q.fa > $W/q.paf; log "wfmash $(( $(date +%s) - t0 )) s, $(wc -l < $W/q.paf) records"
cp $W/q.paf $O.paf.part && mv $O.paf.part $O.paf; cp $W/q.fa.fai $O.fai
log "done: $O.paf"
