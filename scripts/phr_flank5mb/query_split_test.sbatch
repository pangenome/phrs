#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 48
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus02,octopus05,octopus06,tux06,tux07,tux08,tux09
# Does splitting the queries of one all-vs-one wfmash run (the paper's call, wfmash -p 95 --target-prefix T on the whole
# flank FASTA) into K query lists (-A) give the same PAF as the single run? The target index is the same in every run
# and each query is mapped on its own, so the records should be identical. Each part is timed with GNU time -v (CPU,
# peak memory). The union of the parts is compared with REF (the single-run PAF of the same target): identical records
# (all fields except the chain id ch:Z), only-reference and only-split counts.
# usage: sbatch query_split_test.sbatch --flanks FA_GZ --line N --parts K --ref PAF_GZ --out DIR
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --flanks) FA=$2;; --line) L=$2;; --parts) K=$2;; --ref) REF=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in FA L K REF O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
WF=/moosefs/guarracino/pggb_wfmash023/wfmash/build/bin/wfmash; TM=/gnu/store/310psr9h9kzzd4lbssxrzjg53zhqjxrq-time-1.9/bin/time
N=${SLURM_CPUS_PER_TASK:-48}; SEQ=$(sed -n "${L}p" $FA.fai | cut -f1); [ -n "$SEQ" ]
mkdir -p $O /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/qsplit.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
cut -f1 $FA.fai | awk -v k=$K -v w=$W '{print > (w "/q" (NR - 1) % k ".txt")}'
log "target $SEQ; $(wc -l < $FA.fai) queries in $K lists ($(wc -l $W/q*.txt | head -$K | awk '{printf "%s ", $1}')); threads $N; host $(hostname)"
for i in $(seq 0 $((K - 1))); do
  $TM -v -o $W/q$i.time $WF -p 95 --target-prefix "$SEQ" -A $W/q$i.txt -t $N $FA --quiet > $W/q$i.paf
  log "part $i: $(wc -l < $W/q$i.paf) records; $(awk -F': ' '/Elapsed|User time|System time|Maximum resident/ {printf "%s %s; ", $1, $2}' $W/q$i.time | sed 's/\t//g')"
done
python3 - "$REF" $W "$K" "$O/compare.txt" <<'PY'
import gzip, sys
ref, w, k, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
key = lambda l: "\t".join(x for x in l.rstrip("\n").split("\t") if not x.startswith("ch:Z:"))
a = sorted(key(l) for l in gzip.open(ref, "rt"))
b = sorted(key(l) for i in range(k) for l in open(f"{w}/q{i}.paf"))
sa, sb = set(a), set(b)
msg = (f"reference {len(a)} records, split {len(b)}; identical {len(sa & sb)}, only reference {len(sa - sb)}, only split {len(sb - sa)}; "
       f"identical multiset {a == b}")
print(msg); open(out, "w").write(msg + "\n")
PY
cp $W/q*.time $O/
log "done"
