#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05,octopus06
# 5 Mb flanks, acrocentric p and X/Y
# usage: sbatch flank5mb_extract.sbatch --phr-dir DIR --fasta-dir DIR --flank BP --arms "chr13_parm ..." --out FA_GZ
# Same flanks as hprcv2.1Mb.telo_500kb_trimmed.fa.gz (same contigs, same telomere-trimmed start), extended inward to
# FLANK bp (capped at the contig end): p-arm a..a+FLANK-1, q-arm b-FLANK+1..b, names contig:start-end_chrN_arm (1-based).
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --phr-dir) P=$2;; --fasta-dir) FD=$2;; --flank) F=$2;; --arms) AR=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in P FD F AR O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
[ -s $O.fai ] && { log "exists $O.fai, skip"; exit 0; }
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/f5mb.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
for f in $FD/*.fa.gz; do [ -s $f.fai ] && awk -v f=$f '{print $1"\t"f"\t"$2}' $f.fai; done > $W/contigs.tsv
log "contigs indexed: $(wc -l < $W/contigs.tsv) from $(ls $FD/*.fa.gz | wc -l) fasta files"
python3 - $P/hprcv2.1Mb.telo_500kb_trimmed.fa.gz.fai $W/contigs.tsv $F "$AR" > $W/plan.tsv <<'PY'
import re, sys
fai, ct, F, arms = sys.argv[1], sys.argv[2], int(sys.argv[3]), set(sys.argv[4].split())
c = {}
for l in open(ct):
    n, f, L = l.rstrip("\n").split("\t"); c[n] = (f, int(L))
kept = capped = 0
for l in open(fai):
    name = l.split("\t")[0]
    m = re.match(r"^(.+):(\d+)-(\d+)_(chr\w+?)_(parm|qarm)$", name)
    if not m: sys.exit(f"ERROR cannot parse {name}")
    ctg, a, b, ch, arm = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4), m.group(5)
    if f"{ch}_{arm}" not in arms: continue
    if ctg not in c: sys.exit(f"ERROR contig {ctg} not in any fasta")
    f, L = c[ctg]
    if arm == "parm": na, nb = a, min(a + F - 1, L)
    else: na, nb = max(b - F + 1, 1), b
    capped += (nb - na + 1) < F; kept += 1
    print(f"{f}\t{ctg}:{na}-{nb}\t{ctg}:{na}-{nb}_{ch}_{arm}\t{name}")
print(f"flanks {kept}, shorter than {F} bp (contig end) {capped}", file=sys.stderr)
PY
log "plan: $(wc -l < $W/plan.tsv) flanks"
while IFS=$'\t' read -r f r n o; do $S faidx $f "$r" | sed "1s/.*/>$n/"; done < $W/plan.tsv > $W/out.fa
$BGZ -@ 4 -l 9 -c $W/out.fa > $O && $S faidx $O
cut -f3,4 $W/plan.tsv > ${O%.fa.gz}.names.tsv
log "out $O: $(wc -l < $O.fai) flanks, $(awk '{s+=$2} END {print s}' $O.fai) bp; length min $(cut -f2 $O.fai | sort -n | head -1), max $(cut -f2 $O.fai | sort -n | tail -1)"
