#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 6:00:00
#SBATCH --exclude=octopus05,octopus06
# Apply accepted pieces, one person
# usage: sbatch g1_apply.sbatch --washu D --sample SM --version V --runs RUN_DIR --tasks TSV --anchor ex:5000 --vote more --out DIR
# Runs g1_apply_all.py on norm_V; checks that every task-list locus has a summary.tsv or a TIMEOUT marker first.
# out (new dir only): DIR/SM#{1,2}.fa, DIR/corrections.tsv, DIR/loci.tsv (locus status: summary or timeout)
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sample) SM=$2;; --version) V=$2;; --runs) R=$2;; --tasks) TS=$2;; --anchor) AR=$2;; --vote) VR=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D SM V R TS AR VR O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; [ -s $O/corrections.tsv ] && { log "exists $O/corrections.tsv, skip"; exit 0; }
[ -e $O ] && { echo "ERROR $O exists without corrections.tsv; use a new dir"; exit 1; }
[ -s $TS ] || { echo "ERROR empty or missing task list $TS"; exit 1; }
mkdir -p /scratch/$USER/tmp/phrs; LT=$(mktemp /scratch/$USER/tmp/phrs/loci.XXXX); miss=0
while IFS=$'\t' read -r i n e r; do
  if [ -s $R/$i.$n/summary.tsv ]; then echo -e "$i.$n\tsummary" >> $LT
  elif [ -e $R/$i.$n/TIMEOUT ]; then echo -e "$i.$n\ttimeout" >> $LT
  else echo -e "$i.$n\tmissing" >> $LT; miss=$((miss + 1)); fi
done < $TS
log "$SM: task list $TS, $(cut -f2 $LT | sort | uniq -c | tr '\n' ' ')"
[ $miss -eq 0 ] || { echo "ERROR $miss loci without summary or TIMEOUT (list $LT)"; exit 1; }
mkdir -p $O; mv $LT $O/loci.tsv
python3 scripts_tmp/g1_apply_all.py norm_$V $SM $R $AR $VR $O
log "lengths: $(for h in 1 2; do awk '/^>/ {next} {s+=length($0)} END {printf s" "}' $O/$SM#$h.fa; done) vs $(for h in 1 2; do awk '{s+=$2} END {printf s" "}' norm_$V/$SM#$h.fa.gz.fai; done)"
log done
