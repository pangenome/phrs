#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05,octopus06,tux06,tux07,tux08,tux09
# Long flanks (e.g. 5 Mb) of the telomere-anchored flank set, for the arms given (extract_long_flanks.py, which retries
# every moosefs read and checks each sequence length). The FASTA is written on node /scratch, bgzipped and indexed
# there, then copied to moosefs with retries and a size check.
# With --urls FILE (one bgzipped FASTA URL per line) the assemblies are read over HTTPS (HPRC S3) instead of --fasta-dir,
# with CURL_CA_BUNDLE set to the system CA file; every long flank must contain its 500 kb flank (--check-fasta).
# usage: sbatch flank5mb_extract.sbatch --scripts DIR --phr-dir DIR (--fasta-dir DIR | --urls FILE) --flank BP --arms "chr13_parm ..." --out FA_GZ [--jobs N]
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
J=8; FD=""; UR=""
while [ $# -gt 0 ]; do case $1 in --scripts) SD=$2;; --phr-dir) P=$2;; --fasta-dir) FD=$2;; --urls) UR=$2;; --flank) F=$2;; --arms) AR=$2;; --out) O=$2;; --jobs) J=$2;;
  *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in SD P F AR O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
[ -n "$FD$UR" ] && [ -z "$FD" -o -z "$UR" ] || { echo "ERROR give one of --fasta-dir and --urls"; exit 2; }
log () { echo "[$(date -Is)] $*"; }
retry () { local k=1; until "$@"; do [ $k -ge 6 ] && { echo "ERROR failed 6 times: $*"; return 1; }; log "retry $k/5 in $((30 * k)) s: $*"; sleep $((30 * k)); k=$((k + 1)); done; }
put () { cp "$1" "$2.part" && [ "$(stat -c %s "$1")" = "$(stat -c %s "$2.part")" ] && mv "$2.part" "$2"; }  # copy to moosefs, same size
[ -s $O.fai ] && { log "exists $O.fai, skip"; exit 0; }
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/f5mb.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
log "host $(hostname); flank $F; arms $AR; jobs $J; scripts $SD; source ${FD:-$UR ($(wc -l < $UR) URLs)}; scratch free $(df -h /scratch | tail -1)"
SRC=$([ -n "$FD" ] && echo "--fasta-dir $FD" || echo "--urls $UR --idx-dir $W/idx")
python3 $SD/extract_long_flanks.py --flank-fai $P/hprcv2.1Mb.telo_500kb_trimmed.fa.gz.fai $SRC --flank $F --arms "$AR" \
  --check-fasta $P/hprcv2.1Mb.telo_500kb_trimmed.fa.gz --samtools $S --jobs $J --out $W/out.fa --names $W/names.tsv
$BGZ -@ $J -l 9 -c $W/out.fa > $W/out.fa.gz && $S faidx $W/out.fa.gz
log "local $W/out.fa.gz: $(wc -l < $W/out.fa.gz.fai) flanks, $(awk '{s+=$2} END {print s}' $W/out.fa.gz.fai) bp, $(stat -c %s $W/out.fa.gz) bytes"
for x in "out.fa.gz $O" "out.fa.gz.fai $O.fai" "out.fa.gz.gzi $O.gzi" "names.tsv ${O%.fa.gz}.names.tsv"; do set -- $x; retry put $W/$1 $2; done
log "out $O: $(wc -l < $O.fai) flanks; length min $(cut -f2 $O.fai | sort -n | head -1), max $(cut -f2 $O.fai | sort -n | tail -1)"
