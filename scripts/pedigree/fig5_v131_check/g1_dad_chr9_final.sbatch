#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH --exclude=octopus05,octopus06
# Genome-wide G1 DAD plus targeted chr9 patch
# usage: sbatch g1_dad_chr9_final.sbatch --washu D --g1-dir DIR --patch-dir DIR --orig-dir DIR --out DIR
# DIR = folder with PAN011#{1,2}.fa.gz. Junction J = first base where patch and original PAN011#2#chr9 differ.
# In the G1 assembly the 20 kb before J must occur once, exactly; its end J' is the junction there.
# out: G1 chr9 [0, J') + patch chr9 [J, end); hap1 and other contigs as in the G1 assembly.
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --g1-dir) G=$2;; --patch-dir) P=$2;; --orig-dir) R=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D G P R O; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; [ -s $O/PAN011#2.fa.gz.fai ] && { log "exists $O/PAN011#2.fa.gz, skip"; exit 0; }
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
mkdir -p $O /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/dadfinal.XXXX); trap 'echo "exit status $? at line ${LINENO}"; rm -rf $W' EXIT
C='PAN011#2#chr9'
for x in G P R; do $S faidx ${!x}/PAN011#2.fa.gz "$C" > $W/$x.fa; done
python3 - $W <<'PY'
import sys
W = sys.argv[1]
def rd(f): return "".join(l.strip() for l in open(f) if l[0] != ">")
g, p, r = rd(f"{W}/G.fa"), rd(f"{W}/P.fa"), rd(f"{W}/R.fa")
j = next((i for i in range(min(len(p), len(r))) if p[i] != r[i]), None)
if j is None or j < 20000: sys.exit(f"ERROR no junction or junction too close to start: {j}")
a = r[j - 20000:j]; n = g.count(a)
if n != 1: sys.exit(f"ERROR the 20 kb before the junction occurs {n} times in the G1 chr9")
jg = g.find(a) + 20000
print(f"junction: original/patch {j}, G1 {jg}; original length {len(r)}, patch {len(p)}, G1 {len(g)}", file=sys.stderr)
print(f"tail replaced in G1: {len(g) - jg} bp -> {len(p) - j} bp", file=sys.stderr)
open(f"{W}/new.fa", "w").write(">PAN011#2#chr9\n" + g[:jg] + p[j:] + "\n")
PY
$S faidx $G/PAN011#2.fa.gz $(cut -f1 $G/PAN011#2.fa.gz.fai | awk -v c="$C" '$1 != c') > $W/rest.fa
cat $W/rest.fa $W/new.fa | $BGZ -@ 4 -c > $O/PAN011#2.fa.gz && $S faidx $O/PAN011#2.fa.gz
cp $G/PAN011#1.fa.gz $G/PAN011#1.fa.gz.fai $G/PAN011#1.fa.gz.gzi $O/
log "contigs $(wc -l < $O/PAN011#2.fa.gz.fai) (G1 $(wc -l < $G/PAN011#2.fa.gz.fai)); chr9 $(awk -v c="$C" '$1==c {print $2}' $O/PAN011#2.fa.gz.fai) bp"
log "non-ACGTN letters: $(zcat $O/PAN011#2.fa.gz | grep -v '>' | tr -d 'ACGTNacgtn\n' | wc -c)"
log done
