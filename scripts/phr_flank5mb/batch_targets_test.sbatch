#!/bin/bash
#SBATCH -c 48
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05,octopus06,octopus11,tux06,tux07,tux08,tux09
# Can one wfmash run with many targets (-R list, one flank per haplotype so each target is alone in its '#' group and
# the per-group mapping filter acts per target) reproduce the all-vs-one PAFs of the paper (one run per target,
# --target-prefix)? Same wfmash binary and options (-p 95, all flanks as queries). Per target we compare the records
# (all PAF fields except the chain id ch:Z) with the original PAF and report identical, only-original and only-batch
# lines, plus the run's wall and CPU time.
# usage: sbatch batch_targets_test.sbatch --flanks FA_GZ --targets LIST --orig-dir DIR --prefix STR --out DIR
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --flanks) FA=$2;; --targets) TL=$2;; --orig-dir) OD=$2;; --prefix) PFX=$2;; --out) O=$2;;
  *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in FA TL OD PFX O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
WF=/moosefs/guarracino/pggb_wfmash023/wfmash/build/bin/wfmash; N=${SLURM_CPUS_PER_TASK:-48}
mkdir -p $O /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/btt.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
log "host $(hostname), cpus $N, $($WF --version 2>&1 | head -1); targets $(wc -l < $TL) from $(cut -d'#' -f1,2 $TL | sort -u | wc -l) haplotypes"
for e in "" .fai .gzi; do cp $FA$e $W/flanks.fa.gz$e; done
TIMEFORMAT="wfmash real %R s, user %U s, sys %S s"
{ time $WF -p 95 -R $TL -t $N $W/flanks.fa.gz > $W/batch.paf 2> $W/wfmash.err ; } 2> $W/time.txt
cat $W/time.txt; awk '/^\[wfmash\] (Parameters|Filters)/' $W/wfmash.err
log "batch PAF $(wc -l < $W/batch.paf) records"
cp $W/time.txt $O/time.txt
python3 - "$W/batch.paf" "$TL" "$OD" "$PFX" "$O/compare.tsv" <<'EOF'
import collections, gzip, os, sys
paf, tl, od, pfx, out = sys.argv[1:6]
key = lambda l: "\t".join(x for x in l.rstrip("\n").split("\t") if not x.startswith("ch:Z:"))
new = collections.defaultdict(set)
for l in open(paf):
    new[l.split("\t", 6)[5]].add(key(l))
tot = collections.Counter()
with open(out, "w") as fh:
    fh.write("target\toriginal\tbatch\tidentical\tonly_original\tonly_batch\n")
    for t in (x.strip() for x in open(tl) if x.strip()):
        f = os.path.join(od, pfx + t.translate(str.maketrans(":@#-", "____")) + ".paf.gz")
        old = {key(l) for l in gzip.open(f, "rt")}
        b = new.get(t, set()); same = len(old & b)
        r = (len(old), len(b), same, len(old - b), len(b - old)); tot.update(dict(zip(("o", "b", "s", "oo", "ob"), r)))
        fh.write(t + "\t" + "\t".join(map(str, r)) + "\n")
print(f"TOTAL original {tot['o']}, batch {tot['b']}, identical {tot['s']}, only original {tot['oo']}, only batch {tot['ob']}")
EOF
log "done: $O/compare.tsv"
