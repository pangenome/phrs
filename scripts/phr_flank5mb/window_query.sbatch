#!/bin/bash
#SBATCH -c 48
#SBATCH --mem=200G
#SBATCH -t 4:00:00
#SBATCH --exclude=octopus05,octopus06,tux06,tux07,tux08,tux09
# impg query of single caller windows with the caller's flags, on an index of all PAFs (bidirectional, as the caller).
# Per window: records per partner chromosome (from the _chrN_arm suffix), records from contigs whose own name holds
# 'chr' (CHM13, HG002: the March caller counted these under the contig name), and those partner names.
# usage: sbatch window_query.sbatch --paf-dir DIR --impg BIN --regions FILE (seq:start-end per line) --out TSV
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --paf-dir) PD=$2;; --impg) IMPG=$2;; --regions) R=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in PD IMPG R O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
log "host $(hostname); paf dir $PD; impg $IMPG ($($IMPG --version)); regions $(wc -l < $R)"
mkdir -p /scratch/guarracino/tmp/phrs; W=$(mktemp -d /scratch/guarracino/tmp/phrs/wq.XXXX)
trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
mkdir $W/paf; find $PD -maxdepth 1 -name "*.paf.gz" -print0 | xargs -0 cp -t $W/paf/; find $W/paf -name "*.paf.gz" | sort > $W/paf.list
log "copied $(wc -l < $W/paf.list) PAFs"
$IMPG index --alignment-list $W/paf.list -i $W/idx.impg -t ${SLURM_CPUS_PER_TASK:-48}; log "index done"
: > $O.part
while read -r reg; do
  $IMPG query --alignment-list $W/paf.list -i $W/idx.impg -r "$reg" --min-identity 0.95 --min-output-length 3000 --threads 1 -v 0 > $W/q.bed
  log "$reg: $(wc -l < $W/q.bed) records"
  awk -F'\t' -v reg="$reg" '{n = split($1, p, "_"); c = p[n-1]; k[c]++; split($1, h, "#"); split(h[3], g, ":"); if (g[1] ~ /chr/) {m[c]++; nm[c] = nm[c] " " $1}}
    END {for (c in k) print reg "\t" c "\t" k[c] "\t" (m[c] + 0) "\t" (nm[c] == "" ? "." : nm[c])}' $W/q.bed | sort -k3,3nr >> $O.part
done < $R
mv $O.part $O
log "done: $O"
