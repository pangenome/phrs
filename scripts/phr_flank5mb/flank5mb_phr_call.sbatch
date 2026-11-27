#!/bin/bash
#SBATCH -c 48
#SBATCH --mem=340G
#SBATCH -t 3-00:00:00
#SBATCH --exclude=octopus05,octopus06,tux06,tux07,tux08,tux09
# PHR calls with the paper's caller and parameters (scripts/phr-detection/PIPELINE.md, stage 4) on a set of all-vs-one
# PAFs: copy the PAFs to node /scratch, impg index, find-multichr-regions-incremental.py with window 5 kb, step 5 kb,
# identity 0.95, output length 3 kb, >= 2 chromosomes with >= 5 records, stop after 4 failing windows.
# usage: sbatch flank5mb_phr_call.sbatch --paf-dir DIR --fai FAI --max BP --impg BIN --caller PY --out PREFIX [--sequences FILE]
#        [--keep FILE --prefix STR]: only the PAFs of the targets in FILE (named PREFIX<target, :@#- as _>.paf.gz), and in
#        them only the records whose query is in FILE (an all-vs-all among the FILE flanks)
# Outputs PREFIX.calls.tsv, PREFIX.per_window.tsv, PREFIX.caller.log
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
SEQS=""; KEEP=""; PFX=""
while [ $# -gt 0 ]; do case $1 in --paf-dir) PD=$2;; --fai) FAI=$2;; --max) MAX=$2;; --impg) IMPG=$2;; --caller) PY=$2;;
  --out) O=$2;; --sequences) SEQS=$2;; --keep) KEEP=$2;; --prefix) PFX=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in PD FAI MAX IMPG PY O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
N=${SLURM_CPUS_PER_TASK:-48}; TH=$(( N * 5 / 6 ))
log "host $(hostname), cpus $N, caller threads $TH; paf dir $PD ($(ls $PD | grep -c 'paf.gz$') PAFs); fai $FAI ($(wc -l < $FAI) sequences); max $MAX; impg $IMPG ($($IMPG --version)); caller $PY (md5 $(md5sum < $PY | cut -c1-32)); sequences ${SEQS:-all}; keep ${KEEP:-all} (prefix '$PFX')"
mkdir -p /scratch/guarracino/tmp/phrs $(dirname $O); W=$(mktemp -d /scratch/guarracino/tmp/phrs/phrcall.XXXX)
trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
t0=$(date +%s); mkdir $W/paf
if [ -z "$KEEP" ]; then find $PD -maxdepth 1 -name "*.paf.gz" -print0 | xargs -0 cp -t $W/paf/
else
  BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
  python3 -c 'import sys; [print(sys.argv[2] + "/" + sys.argv[3] + l.strip().translate(str.maketrans(":@#-", "____")) + ".paf.gz") for l in open(sys.argv[1]) if l.strip()]' $KEEP $PD "$PFX" > $W/sel.txt
  while read f; do [ -s $f ] || { echo "ERROR missing PAF $f"; exit 1; }; done < $W/sel.txt
  mkdir $W/cnt; export KEEP BGZ W
  xargs -P $N -I{} bash -c 'set -eo pipefail; b=$(basename {}); zcat {} | awk -F"\t" -v c=$W/cnt/$b "NR == FNR {k[\$1]; next} {n++} (\$1 in k) {m++; print} END {print m + 0, n + 0 > c}" $KEEP - | $BGZ > $W/paf/$b' < $W/sel.txt
  log "kept $(wc -l < $W/sel.txt) target PAFs, records with a query in $KEEP: $(cat $W/cnt/* | awk '{m += $1; n += $2} END {print m " of " n}') (counted while filtering)"
fi
log "copied $(ls $W/paf | wc -l) PAFs, $(du -sh $W/paf | cut -f1), $(( $(date +%s) - t0 )) s"
find $W/paf -name "*.paf.gz" | sort > $W/paf.list
t0=$(date +%s); $IMPG index --alignment-list $W/paf.list -i $W/idx.impg -t $N; log "impg index $(du -sh $W/idx.impg | cut -f1), $(( $(date +%s) - t0 )) s"
mkdir $W/bin; ln -s $IMPG $W/bin/impg; export PATH=$W/bin:$PATH; [ "$(command -v impg)" = $W/bin/impg ] || { echo "ERROR impg on PATH is $(command -v impg)"; exit 1; }
t0=$(date +%s)
python3 $PY --paf-list $W/paf.list -i $W/idx.impg --fai $FAI --window 5000 --step 5000 --max $MAX \
  --min-identity 0.95 --min-output-length 3000 --min-diff-chrs 2 --min-count 5 --min-consecutive-stop 4 \
  ${SEQS:+--sequences $SEQS} --per-window-output $W/per_window.tsv --threads $TH --verbose > $W/calls.tsv 2> $W/caller.log
log "caller $(( $(date +%s) - t0 )) s; $(tail -3 $W/caller.log | tr '\n' ' ')"
log "calls: $(( $(wc -l < $W/calls.tsv) - 1 )) sequences, $(awk -F'\t' 'NR > 1 && $4 != "."' $W/calls.tsv | wc -l) with a region; per-window rows $(( $(wc -l < $W/per_window.tsv) - 1 ))"
for f in calls.tsv per_window.tsv caller.log; do cp $W/$f $O.$f.part && mv $O.$f.part $O.$f; done
log "done: $O.{calls.tsv,per_window.tsv,caller.log}"
