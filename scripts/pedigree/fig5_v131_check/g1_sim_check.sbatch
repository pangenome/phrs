#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 2:00:00
#SBATCH --exclude=octopus05
# Does each piece restore the original?
# usage: sbatch g1_sim_check.sbatch --washu D --sim SIM_DIR --runs RUN_DIR --sample SM --version V --hap H --out TSV
# For every piece in RUN_DIR/<event>/pieces.fa (header: id contig start end mode rule contig): patched =
# modified[start-50k, start) + piece + modified[end, end+50k); original = norm_V on the matching interval (shift from
# truth.tsv); minimap2 2.31 asm5 --eqx best record: NM and length difference.
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sim) SD=$2;; --runs) RD=$2;; --sample) SM=$2;; --version) V=$2;; --hap) H=$2;; --out) OUT=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D SD RD SM V H OUT; do [ -n "${!v:-}" ] || { echo "ERROR missing $v"; exit 2; }; done
cd $D; S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/simchk.XXXX)
printf "event\tid\tcontig\tpiece_start\tpiece_end\tpatched_bp\toriginal_bp\tmatches\tnm\tlength_diff\n" > $OUT
while IFS=$'\t' read -r i n e r; do
  P=$RD/$i.$n/pieces.fa; [ -s $P ] || continue
  read -r o0 o1 m0 m1 <<< "$(awk -F'\t' -v i=$i 'NR == i + 1 {print $4, $5, $6, $7}' $SD/truth.tsv)"
  [ -n "${m1:-}" ] || { echo "ERROR event $i not in truth.tsv"; exit 1; }
  awk '/^>/ {h=substr($0, 2); next} {print h"\t"$0}' $P | while IFS=$'\t' read -r h seq; do
    read -r id c a z rest <<< "$h"
    { echo ">patched"; { $S faidx $SD/norm/$SM#$H.fa.gz "$c:$((a - 49999))-$a" | tail -n +2; echo "$seq"
      $S faidx $SD/norm/$SM#$H.fa.gz "$c:$((z + 1))-$((z + 50000))" | tail -n +2; } | tr -d '\n'; echo; } > $W/p.fa
    $S faidx norm_$V/$SM#$H.fa.gz "$c:$((a - 49999 - (m0 - o0)))-$((z - (m1 - o1) + 50000))" | sed '1s/.*/>original/' > $W/o.fa
    $MM -cx asm5 --eqx $W/o.fa $W/p.fa 2> /dev/null | sort -k10,10nr | head -1 > $W/a.paf
    PL=$(tail -n +2 $W/p.fa | tr -d '\n' | wc -c); OL=$(tail -n +2 $W/o.fa | tr -d '\n' | wc -c)
    MT=$(awk '{print $10 + 0}' $W/a.paf); NM=$(awk '{for (k=13; k<=NF; k++) if ($k ~ /^NM:i:/) print substr($k, 6)}' $W/a.paf)
    printf "$i.$n\t$id\t$c\t$a\t$z\t$PL\t$OL\t${MT:-0}\t${NM:-NA}\t$((PL - OL))\n" >> $OUT
  done
done < $SD/locus_tasks.tsv
echo "[$(date -Is)] $(($(wc -l < $OUT) - 1)) pieces checked"
rm -rf "$W"
