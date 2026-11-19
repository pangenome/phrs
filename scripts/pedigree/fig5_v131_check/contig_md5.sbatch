#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 0-04:00:00
#SBATCH --exclude=octopus05,octopus06
# Per-contig length and MD5, v1.1 vs v1.3.1
# args: --washu DIR --out TSV
set -euxo pipefail
trap 'echo "exit $? at line ${LINENO}"' EXIT
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
while [ $# -gt 0 ]; do case $1 in --washu) D=$2; shift 2;; --out) O=$2; shift 2;; *) echo "unknown arg $1"; exit 1;; esac; done
[ -n "${D:-}" ] && [ -n "${O:-}" ] || { echo "need --washu and --out"; exit 1; }
log() { echo "[$(date -Is)] $*"; }
[ -s "$O" ] && { log "skip: $O exists"; exit 0; }
cd $D; : > $O.tmp
for V in v11 v131; do for f in norm_$V/PAN0*#[12].fa.gz; do
  log "dict $V $f"
  $S dict $f | awk -F'\t' -v v=$V '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); sub("M5:","",$4); print v"\t"$2"\t"$3"\t"$4}' >> $O.tmp
done; done
mv $O.tmp $O
log "done: $(wc -l < $O) rows; same md5 in both: $(cut -f2,4 $O | sort | uniq -d | wc -l)"
