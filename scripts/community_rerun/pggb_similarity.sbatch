#!/bin/bash
#SBATCH -w octopus02
#SBATCH -p allnodes
#SBATCH -c 48
#SBATCH --mem=950G
#SBATCH -t 4-00:00:00
# PHR pangenome graph and pairwise similarity with the paper's commands (code_PHR-and-3D.md, "Pangenome graph
# construction and pairwise similarity"): pggb -p 95 -t 48 on the PHR FASTA, then odgi similarity -t 48 -P --all on the
# final graph. Tools: the paper's builds in /moosefs/guarracino/pggb_wfmash023 (pggb v0.7.4-26-g498c5d7, wfmash
# v0.23.0-41, seqwish v0.7.11-4, smoothxg v0.8.2-2 with its odgi); GFAffix v0.2.1 from the guix store, because the
# paper's GFAffix build on moosefs is no longer readable (same version, other build).
# --stage S (optional): copy the input FASTA (.fai, .gzi) and the four binaries to the node-local dir S first (same
# basename, same builds), so the run reads almost nothing from moosefs; put --out on node-local disk too.
# usage: sbatch pggb_similarity.sbatch --fasta FA_GZ --out DIR [--stage S]
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --fasta) FA=$2;; --out) O=$2;; --stage) S=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in FA O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
B=/moosefs/guarracino/pggb_wfmash023; GX=/gnu/store/4yjldswcvgyna7hzpiqdlr7q5ll50c4m-gfaffix-0.2.1/bin
TM=/gnu/store/310psr9h9kzzd4lbssxrzjg53zhqjxrq-time-1.9/bin  # pggb needs GNU time on PATH; the nodes have none now
export PATH="$B/wfmash/build/bin:$B/seqwish/bin:$B/smoothxg/bin:$GX:$B/smoothxg/deps/odgi/bin:$B/pggb:$TM:$PATH"
if [ -n "${S:-}" ]; then
  mkdir -p $S/input $S/bin
  for f in $FA $FA.fai $FA.gzi $B/wfmash/build/bin/wfmash $B/seqwish/bin/seqwish $B/smoothxg/bin/smoothxg $B/smoothxg/deps/odgi/bin/odgi; do
    d=$S/bin; case $f in $FA*) d=$S/input;; esac
    cp $f $d/ && [ $(stat -c %s $f) -eq $(stat -c %s $d/$(basename $f)) ] || { echo "ERROR copy of $f"; exit 1; }
    cmp -s $f $d/$(basename $f) || { echo "ERROR $d/$(basename $f) differs from $f"; exit 1; }
  done
  FA=$S/input/$(basename $FA); export PATH="$S/bin:$PATH"
  log "staged input and binaries in $S: $(ls -la $S/input $S/bin | grep -c '^-') files"
fi
for t in wfmash seqwish smoothxg odgi; do log "$t: $(which $t)"; done
[ "$(which time)" = $TM/time ] || { echo "ERROR GNU time not on PATH (which time: $(which time))"; exit 1; }  # pggb uses which time
log "host $(hostname); $(wfmash --version 2>&1 | head -1); seqwish $(seqwish --version 2>&1 | head -1); smoothxg $(smoothxg --version 2>&1 | head -1); $(gfaffix --version 2>&1 | head -1); odgi $(odgi version 2>&1 | head -1); $(pggb --version 2>&1 | head -1)"
log "input $FA: $(wc -l < $FA.fai) sequences, $(awk '{s+=$2} END {print s}' $FA.fai) bp"
mkdir -p $O /scratch/guarracino/tmp/phrs; T=$(mktemp -d /scratch/guarracino/tmp/phrs/pggb.XXXX); trap 'rc=$?; rm -rf $T; echo "exit status $rc at line ${LINENO}"' EXIT
t0=$(date +%s)
pggb -i $FA -o $O -p 95 -D $T -t 48
log "pggb $(( $(date +%s) - t0 )) s"
G=$(ls $O/*.smooth.final.og); [ $(echo "$G" | wc -l) -eq 1 ] || { echo "ERROR expected one final .og in $O"; exit 1; }
odgi stats -i $G -S | tail -2
t0=$(date +%s)
odgi similarity -i $G -t 48 -P --all | pigz -9 > $O/$(basename $G .og).similarity.tsv.gz.part
mv $O/$(basename $G .og).similarity.tsv.gz.part $O/$(basename $G .og).similarity.tsv.gz
log "odgi similarity $(( $(date +%s) - t0 )) s; $(ls -la $O/$(basename $G .og).similarity.tsv.gz)"
