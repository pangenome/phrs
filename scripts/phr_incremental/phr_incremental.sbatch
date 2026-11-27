#!/bin/bash
#SBATCH -c 64
#SBATCH --mem=200G
#SBATCH -t 2-00:00:00
#SBATCH --exclude=octopus05,octopus06,tux06,tux07,tux08,tux09
# Incremental PHR calls (phr_incremental.py) for a list of flanks, all on one node. moosefs reads can fail now and then
# (EIO), so we copy the flank FASTA and the wfmash binary to node /scratch with retries and a size check, work only
# there, and copy the results back to moosefs with retries at the end (also when the job fails).
# usage: sbatch phr_incremental.sbatch --flanks FA_GZ --ends FILE --scripts DIR --out DIR [--jobs 13 --threads 4 --max-rounds N --mode frontier|all]
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
J=13; TH=4; MR=0; MODE=frontier
while [ $# -gt 0 ]; do case $1 in --flanks) FA=$2;; --ends) E=$2;; --scripts) SD=$2;; --out) O=$2;; --jobs) J=$2;; --threads) TH=$2;;
  --max-rounds) MR=$2;; --mode) MODE=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in FA E SD O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
retry () { local k=1; until "$@"; do [ $k -ge 6 ] && { echo "ERROR failed 6 times: $*"; return 1; }; log "retry $k/5 in $((30 * k)) s: $*"; sleep $((30 * k)); k=$((k + 1)); done; }
get () { cp "$1" "$2" && [ "$(stat -c %s "$1")" = "$(stat -c %s "$2")" ]; }  # copy, same size
IMPG=/gnu/store/p8prk9q2g3cp4331xh0gkylyk2zv0lvn-impg-0.4.1-1.f773342/bin/impg
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
mkdir -p /scratch/guarracino/tmp/phrs; W=$(mktemp -d /scratch/guarracino/tmp/phrs/inc.XXXX); mkdir -p $W/res/logs $W/bin
save () { mkdir -p $O/logs; for f in $(find $W/res -type f); do retry get $f $O/${f#$W/res/} || echo "ERROR could not copy $f"; done; }
trap 'rc=$?; save; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
B=/moosefs/guarracino/pggb_wfmash023/wfmash/build; mkdir -p $W/lib  # wfmash finds its libraries by RUNPATH, so LD_LIBRARY_PATH wins
retry get $B/bin/wfmash $W/bin/wfmash; chmod +x $W/bin/wfmash; WF=$W/bin/wfmash
for l in libwflign.so libwfa2cpp.so.0 libwfa2.so.0; do retry get $B/lib/$l $W/lib/$l; done; export LD_LIBRARY_PATH=$W/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
[ "$(ldd $WF | grep -c "$W/lib")" = 3 ] || { echo "ERROR wfmash does not load the local libraries: $(ldd $WF | grep -E 'wflign|wfa2')"; exit 1; }
for f in phr_incremental.py split_flanks.py; do retry get $SD/$f $W/bin/$f; done
retry get $E $W/ends.txt
log "host $(hostname); flanks $FA; ends $(wc -l < $W/ends.txt); jobs $J x threads $TH; mode $MODE; $($WF --version 2>&1 | head -1); $($IMPG --version)"
t0=$(date +%s); retry get $FA $W/flanks.fa.gz; log "copied flank FASTA, $(stat -c %s $W/flanks.fa.gz) bytes, $(( $(date +%s) - t0 )) s"
t0=$(date +%s); python3 $W/bin/split_flanks.py $W/flanks.fa.gz $W/flanks; log "split $(find $W/flanks -name "*.fa" | wc -l) flanks, $(( $(date +%s) - t0 )) s"
export W WF IMPG S TH MR MODE LD_LIBRARY_PATH
san () { echo "$1" | tr ':@#-' '____'; }; export -f san
xargs -P $J -I{} bash -c 'set -eo pipefail; n=$(san "{}"); mkdir -p $W/tmp/$n; TMPDIR=$W/tmp/$n python3 $W/bin/phr_incremental.py --end "{}" --flank-dir $W/flanks --wfmash $WF --impg $IMPG --samtools $S --threads $TH --max-rounds $MR --mode $MODE --out $W/res/$n.tsv 2> $W/res/logs/$n.log' < $W/ends.txt
log "done: $(find $W/res -maxdepth 1 -name "*.tsv" | wc -l) of $(wc -l < $W/ends.txt) ends; summaries:"; grep -h "^#SUMMARY" $W/res/*.tsv
