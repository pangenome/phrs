#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 2-00:00:00
#SBATCH --exclude=octopus05,octopus06
# Smallest N for the Fig5 sweepGA filter (--num-mappings N:N --scaffold-jump 0 --scoring ani) that keeps the transmitted
# parent copy at the given ends (10:10 drops it at SELF.p chr2p, chr9p, KID chr12p, chr4q on v1.1 -f32)
# usage: sbatch sweepga_nscan.sbatch --washu D --raw PAF.gz --label L --ends "CONTIG:END:COPY,..." --n-from 10 --n-to 40 --out DIR
# Same sweepga binary as Fig5 step 5 (fig5_similarity.sbatch). For N = from, from+1, ...: filter the whole raw PAF (the target
# dimension depends on every query), then coverage of the last (q) or first (p) 1 Mb by the transmitted copy of the same
# chromosome (sweepga_end_reach.py); stops at the first N where every listed end has >= 950,000 bp covered and < 10,000 bp
# uncovered at the end. N = 10 reproduces the Fig5 filter.
# out: DIR/<label>.tsv (N, seconds, filtered lines, per end: covered, uncovered at end, pass)
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --raw) RAW=$2;; --label) LAB=$2;; --ends) ENDS=$2;; --n-from) N0=$2;; --n-to) N1=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D RAW LAB ENDS N0 N1 O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; mkdir -p $O; [ -s $RAW ] || { echo "ERROR missing $RAW"; exit 1; }
[ -s $O/$LAB.tsv ] && grep -q -P "\tall_pass$" $O/$LAB.tsv && { log "exists $O/$LAB.tsv with a pass, skip"; exit 0; }
SW=/home/erikg/.cargo/bin/sweepga; T=${SLURM_CPUS_PER_TASK:-4}
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/nscan.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT; trap 'exit 143' TERM
log "$LAB: raw $RAW, ends $ENDS, N $N0..$N1, sweepga $($SW --version 2>&1 | head -1), host $(hostname)"
zcat $RAW > $W/raw.paf; log "raw lines $(wc -l < $W/raw.paf)"
printf "n\tseconds\tfiltered_lines\tends\tresult\n" > $O/$LAB.tsv
for N in $(seq $N0 $N1); do
  t0=$(date +%s)
  $SW --num-mappings $N:$N --scaffold-jump 0 --scoring ani --temp-dir $W --output-file $W/filt.paf $W/raw.paf > $W/sw.log 2>&1 || { log "ERROR sweepga N=$N"; tail -5 $W/sw.log; exit 1; }
  sec=$(( $(date +%s) - t0 )); nl=$(wc -l < $W/filt.paf)
  python3 scripts_tmp/sweepga_end_reach.py $LAB < $W/filt.paf > $W/reach.tsv 2> /dev/null
  res=$(python3 - $W/reach.tsv "$ENDS" <<'PY'
import sys
rows = {(r[1], r[3]): r for r in (l.rstrip("\n").split("\t") for l in list(open(sys.argv[1]))[1:])}
out, ok = [], True
for e in sys.argv[2].split(","):
    c, end, cp = e.rsplit(":", 2); r = rows.get((c, end))
    if r is None: sys.exit(f"ERROR end {c} {end} not in end-reach output")
    u, cov = (int(r[4]), int(r[5])) if cp == "h1" else (int(r[6]), int(r[7]))
    p = cov >= 950000 and u < 10000; ok &= p
    out.append(f"{c.split('#')[-1]}{end}:{cp}:{cov}:{u}:{'pass' if p else 'fail'}")
print(" ".join(out) + "\t" + ("all_pass" if ok else "fail"))
PY
)
  printf "$N\t$sec\t$nl\t$res\n" >> $O/$LAB.tsv
  log "$LAB N=$N ($sec s, $nl lines): $res"
  case "$res" in *all_pass) break;; esac
done
cat $O/$LAB.tsv
log done
