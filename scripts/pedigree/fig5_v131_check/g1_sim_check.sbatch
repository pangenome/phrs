#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 2:00:00
#SBATCH --exclude=octopus05
# Does each piece restore the original?
# usage: sbatch g1_sim_check.sbatch --washu D --sim SIM_DIR --runs RUN_DIR --sample SM --version V --hap H --out TSV
#   [--tasks TSV (default SIM_DIR/locus_tasks.tsv; first column = truth.tsv row)]
# For every piece in RUN_DIR/<event>/pieces.fa (header: id contig start end mode rule contig): patched =
# modified[start-50k, start) + piece + modified[end, end+50k) (clamped at contig ends); original = norm_V on the matching interval (shift from
# truth.tsv: x - shift before the event or inside an event of unchanged length, x - shift after it, the event start or end for x
# inside an event that changes length); minimap2 2.31 asm5 --eqx
# best record: NM and length difference.
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sim) SD=$2;; --runs) RD=$2;; --sample) SM=$2;; --version) V=$2;; --hap) H=$2;; --out) OUT=$2;; --tasks) TK=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
TK=${TK:-$SD/locus_tasks.tsv}
for v in D SD RD SM V H OUT; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
cd $D; S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/simchk.XXXX)
# terminal pieces: the telomere at a contig end (chain of TTAGGG/CCCTAA starts, each <= 50 bp from the next, the outermost
# <= 50 bp from the end) is cut from patched and original before the comparison; local assemblies carry other telomere lengths
cat > $W/trim.py <<'PY'
import re, sys
W, st, en = sys.argv[1], sys.argv[2] == "1", sys.argv[3] == "1"
def rc(x): return x[::-1].translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))
def tel(s):
    b = len(s)
    for h in reversed([m.start() for m in re.finditer(r"(?=TTAGGG|CCCTAA)", s.upper())]):
        if b - (h + 6) > 50: break
        b = h
    return len(s) - b
out = []
for f in ("p", "o"):
    h, s = open(f"{W}/{f}.fa").read().split("\n", 1); s = s.replace("\n", "")
    t0 = tel(rc(s)) if st else 0; t1 = tel(s) if en else 0
    if t0 + t1 >= len(s): sys.exit(f"ERROR {f}.fa is all telomere ({len(s)} bp)")
    open(f"{W}/{f}.fa", "w").write(f"{h}\n{s[t0:len(s) - t1]}\n"); out.append(t0 + t1)
print(*out)
PY
om () { if [ $1 -le $m0 ] || [ $((m1 - m0)) -eq $((o1 - o0)) ]; then echo $(($1 - (m0 - o0))); elif [ $1 -ge $m1 ]; then echo $(($1 - (m1 - o1))); else echo $2; fi; }  # sim -> norm_V
printf "event\tid\tcontig\tpiece_start\tpiece_end\tpatched_bp\toriginal_bp\tmatches\tnm\tlength_diff\tpatched_telo_bp\toriginal_telo_bp\n" > $OUT
while IFS=$'\t' read -r i n e r; do
  P=$RD/$i.$n/pieces.fa; [ -s $P ] || continue
  read -r o0 o1 m0 m1 <<< "$(awk -F'\t' -v i=$i 'NR == i + 1 {print $4, $5, $6, $7}' $SD/truth.tsv)"
  [ -n "${m1:-}" ] || { echo "ERROR event $i not in truth.tsv"; exit 1; }
  awk '/^>/ {h=substr($0, 2); next} {print h"\t"$0}' $P | while IFS=$'\t' read -r h seq; do
    read -r id c a z rest <<< "$h"
    # flanks up to 50 kb, clamped at the contig ends (terminal pieces start at 0 or end at the contig length)
    CL=$(awk -v c="$c" '$1 == c {print $2}' $SD/norm/$SM#$H.fa.gz.fai); L=$((a < 50000 ? a : 50000)); R=$((CL - z < 50000 ? CL - z : 50000))
    [ -n "$CL" ] && [ $L -ge 0 ] && [ $R -ge 0 ] || { echo "ERROR piece $id $c $a-$z, contig length '$CL'"; exit 1; }
    { echo ">patched"; { [ $L -eq 0 ] || $S faidx $SD/norm/$SM#$H.fa.gz "$c:$((a - L + 1))-$a" | tail -n +2; echo "$seq"
      [ $R -eq 0 ] || $S faidx $SD/norm/$SM#$H.fa.gz "$c:$((z + 1))-$((z + R))" | tail -n +2; } | tr -d '\n'; echo; } > $W/p.fa
    OS=$(om $((a - L)) $o0); OE=$(om $((z + R)) $o1)
    $S faidx norm_$V/$SM#$H.fa.gz "$c:$((OS + 1))-$OE" | sed '1s/.*/>original/' > $W/o.fa
    read -r TP TO < <(python3 $W/trim.py $W $((a == L)) $((z + R == CL))); [ -n "${TO:-}" ] || { echo "ERROR telomere trim failed on $id"; exit 1; }
    $MM -cx asm5 --eqx $W/o.fa $W/p.fa 2> /dev/null | sort -k10,10nr | head -1 > $W/a.paf
    PL=$(tail -n +2 $W/p.fa | tr -d '\n' | wc -c); OL=$(tail -n +2 $W/o.fa | tr -d '\n' | wc -c)
    MT=$(awk '{print $10 + 0}' $W/a.paf); NM=$(awk '{for (k=13; k<=NF; k++) if ($k ~ /^NM:i:/) print substr($k, 6)}' $W/a.paf)
    printf "$i.$n\t$id\t$c\t$a\t$z\t$PL\t$OL\t${MT:-0}\t${NM:-NA}\t$((PL - OL))\t$TP\t$TO\n" >> $OUT
  done
done < $TK
echo "[$(date -Is)] $(($(wc -l < $OUT) - 1)) pieces checked"
rm -rf "$W"
