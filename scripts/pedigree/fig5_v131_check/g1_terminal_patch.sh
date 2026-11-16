#!/bin/bash
# Replace a chromosome end with an own contig
# usage: bash g1_terminal_patch.sh --in-fa IN.fa.gz --contig NAME#H#chrN --win-start BP --ctg-fa CONTIGS.fa --ctg-name NAME --out-fa OUT.fa.gz
# Aligns the contig to CONTIG[WIN_START, end) (minimap2 2.31 asm5 --eqx); the junction is the end of the first exact
# match run >= 20 kb; CONTIG[junction, end) is replaced by the contig from the matching contig position to its end.
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --in-fa) IN=$2;; --contig) C=$2;; --win-start) WS=$2;; --ctg-fa) CF=$2;; --ctg-name) CN=$2;; --out-fa) OUT=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in IN C WS CF CN OUT; do [ -n "${!v:-}" ] || { echo "ERROR missing --${v}"; exit 2; }; done
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
log () { echo "[$(date -Is)] $*"; }
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/tpatch.XXXX)
L=$(awk -v c="$C" '$1==c {print $2}' $IN.fai); [ -n "$L" ] || { echo "ERROR $C not in $IN.fai"; exit 1; }
log "in $IN, contig $C length $L, window $WS-$L, piece $CN from $CF"
$S faidx $IN "$C:$((WS + 1))-$L" > $W/t.fa; $S faidx $CF "$CN" > $W/q.fa
$MM -cx asm5 --eqx --secondary=no $W/t.fa $W/q.fa 2> $W/mm.log > $W/a.paf
read TJ QJ < <(python3 - $W/a.paf <<'PY'
import re, sys
best = None
for l in open(sys.argv[1]):
    f = l.split("\t")
    if f[4] != "+" or (best and int(f[3]) - int(f[2]) <= best[0]): continue
    best = (int(f[3]) - int(f[2]), f)
f = best[1]; q, t = int(f[2]), int(f[7])
for n, op in re.findall(r"(\d+)([=XID])", next(x[5:] for x in f[12:] if x.startswith("cg:Z:"))):
    n = int(n)
    if op == "=" and n >= 20000: print(t + n, q + n); sys.exit(0)
    if op in "=X": q += n; t += n
    elif op == "I": q += n
    else: t += n
sys.exit("ERROR no exact run >= 20 kb")
PY
)
J=$((WS + TJ)); CL=$(awk -v c="$CN" '$1==c {print $2}' $CF.fai)
log "junction $C:$J (window offset $TJ), contig position $QJ of $CL; replacing [$J, $L) = $((L - J)) bp by $((CL - QJ)) bp"
python3 $D/scripts_tmp/g1_apply_patch.py $IN $W/out.fa "$C" $J $L $CF "$CN" $QJ $CL +
$BGZ -@ 4 -c $W/out.fa > $OUT && $S faidx $OUT
log "out $OUT: $(grep -P "^$C\t" $OUT.fai | cut -f1,2 | tr '\t' ' ')"
rm -rf "$W"
